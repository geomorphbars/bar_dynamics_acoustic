# =============================================================================
# DATA PROCESSING 
# 
# Bedload flux over migrating bars revealed by high-resolution passive acoustic mapping
# (Loire River, France)
#
#
# This script only contains the DATA PROCESSING steps:
#   - loading raw acoustic / bathymetric files
#   - converting acoustic power -> bedload flux
#   - matching acoustic points to the nearest bathymetric point
#   - classifying points into morphological zones ( bar stoss side / wake zone /
#     low-flow channel) from digitized polygons
#   - computing the summary statistics used in the associated publication
#     (per-zone statistics, percentiles, flux contribution, proximity index,
#     width-to-depth ratio)
#
#
# >>> Before running, set DATA_DIR below and review the `campaigns` list <<<
# =============================================================================


# -----------------------------------------------------------------------------
# 1. LIBRARIES
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(sp)     # point.in.polygon(), used for morphological zone classification


# -----------------------------------------------------------------------------
# 2. CONFIGURATION - ADAPT THIS SECTION TO YOUR OWN DATA
# -----------------------------------------------------------------------------

# Folder containing the raw acoustic and bathymetric files.
# Replace with the path where you stored the data, e.g. "./data" or an
# absolute path on your own machine.
DATA_DIR <- "data/raw"

# Folder where the result tables (CSV) will be written.
OUTPUT_DIR <- "results"


# Empirical conversion constant (calibrated for this study site)
CONVERSION_P_QB    <- 7.2e-9   # Qb [g.s-1.m-1] = 7.2e-9 * P^0.75

# Maximum distance (m) allowed when matching an acoustic point to its
# nearest bathymetric point
BATHY_DISTANCE_THRESHOLD <- 5

# Fixed channel width (m) used for the width-to-depth ratio
W_FIXED <- 500

# Reference bar height (m): interpercentile range P90-P10 of bed elevations,
# averaged over the six surveys with a submerged bar. See manuscript section 2.
H_BAR <- 2.30


# -----------------------------------------------------------------------------
# 3. BASIC LOADING / UNIT CONVERSION FUNCTIONS
# -----------------------------------------------------------------------------

#' Load a raw, header-less data file (acoustic or bathymetric)
#'
#' @param file_path path to the data file
#' @param type "acoustic" (x, y, P, PdB, Fc, Fp) or "bathy" (x, y, altitude, depth)
#' @return data frame with named columns
load_raw_data <- function(file_path, type = "acoustic") {
  
  data <- read.table(file_path, header = FALSE)
  
  if (type == "acoustic") {
    colnames(data) <- c("x", "y", "P", "PdB", "Fc", "Fp")
  } else if (type == "bathy") {
    colnames(data) <- c("x", "y", "altitude", "depth")
  } else {
    stop("Unknown data type: ", type)
  }
  
  message(nrow(data), " ", type, " points loaded from ", basename(file_path))
  return(data)
}

#' Convert acoustic power to bedload flux
#' Qb = CONVERSION_P_QB * P^0.75
convert_power_to_flux <- function(data) {
  if (!"P" %in% colnames(data)) stop("Column 'P' (acoustic power) not found")
  data$P[data$P <= 0] <- NA
  data$Qb <- CONVERSION_P_QB * (data$P ^ 0.75)
  return(data)
}

#' Compute a "proximity index" and the relative submergence from altitude
#' and water depth:
#'   proximity_index       = altitude / water_surface_altitude
#'                           (1 = point near the surface, 0 = deeply submerged)
#'   relative_submergence  = water_height / H_BAR
#'                           (local flow depth in units of bar height)
compute_proximity_index <- function(data) {
  
  if (!all(c("altitude", "depth") %in% colnames(data))) {
    warning("Columns 'altitude' and/or 'depth' missing - proximity index not computed")
    return(data)
  }
  
  data$water_height   <- abs(data$depth)
  data$water_altitude  <- data$altitude + data$water_height
  data$proximity_index <- data$altitude / data$water_altitude
  
  # Relative submergence = local flow depth normalized by bar height
  # (dimensionless; the index used in the manuscript)
  data$relative_submergence <- data$water_height / H_BAR
  
  # Handle division-by-zero edge cases
  data$proximity_index[is.infinite(data$proximity_index)] <- NA
  
  return(data)
}


# -----------------------------------------------------------------------------
# 4. BUILDING THE FULL ACOUSTIC DATASET (ALL CAMPAIGNS)
# -----------------------------------------------------------------------------

#' Load, convert and stack the acoustic data of every campaign.
#'
#' @param campaigns list of campaigns, each a list with fields:
#'   name, acoustic_file, bathy_file, date, discharge
#' @return data frame with all campaigns combined (Qb already computed)
prepare_full_dataset <- function(campaigns) {
  
  full_data <- data.frame()
  
  for (camp in campaigns) {
    acoustic <- load_raw_data(camp$acoustic_file, "acoustic")
    acoustic <- convert_power_to_flux(acoustic)
    
    acoustic$campaign  <- camp$name
    acoustic$date      <- camp$date
    acoustic$discharge <- camp$discharge
    
    full_data <- rbind(full_data, acoustic)
  }
  
  return(full_data)
}

#' Match every acoustic point of every campaign to its nearest bathymetric
#' point (within `distance_threshold` meters) and compute the derived
#' proximity index on the resulting subset.
#'
#' @param full_data combined acoustic data frame (output of prepare_full_dataset(),
#'   ideally after classify_morphological_zones())
#' @param campaigns list of campaigns (see prepare_full_dataset)
#' @param distance_threshold maximum matching distance (m)
#' @return subset of full_data restricted to points with a valid bathymetric
#'   match, enriched with altitude / depth / proximity_index
associate_bathymetry <- function(full_data, campaigns,
                                 distance_threshold = BATHY_DISTANCE_THRESHOLD) {
  
  results <- data.frame()
  
  for (camp in campaigns) {
    camp_data <- full_data[full_data$campaign == camp$name, ]
    bathy <- load_raw_data(camp$bathy_file, "bathy")
    
    camp_data$altitude       <- NA
    camp_data$depth          <- NA
    camp_data$bathy_distance <- NA
    
    for (j in seq_len(nrow(camp_data))) {
      distances <- sqrt((bathy$x - camp_data$x[j])^2 + (bathy$y - camp_data$y[j])^2)
      idx_min   <- which.min(distances)
      min_dist  <- distances[idx_min]
      
      if (min_dist <= distance_threshold) {
        camp_data$altitude[j]       <- bathy$altitude[idx_min]
        camp_data$depth[j]          <- bathy$depth[idx_min]
        camp_data$bathy_distance[j] <- min_dist
      }
    }
    
    valid <- !is.na(camp_data$altitude)
    if (any(valid)) {
      camp_data_valid <- compute_proximity_index(camp_data[valid, ])
      results <- rbind(results, camp_data_valid)
    }
  }
  
  return(results)
}


# -----------------------------------------------------------------------------
# 5. MORPHOLOGICAL ZONE CLASSIFICATION
# -----------------------------------------------------------------------------

#' Classify every point into one of three morphological zones using
#' campaign-specific polygons:
#'   "Sand bar"    -> bar stoss side (banc de sable)
#'   "Wake zone"   -> lateral pool / bar wake (zone d'abris)
#'   "Outside bar" -> low-flow channel (hors banc / chenal d'etiage), default
#'
#' @param full_data combined acoustic data frame
#' @param bar_polygons  list of polygons, each list(campaign=, x=, y=)
#' @param wake_polygons list of polygons, each list(campaign=, x=, y=)
classify_morphological_zones <- function(full_data, bar_polygons, wake_polygons) {
  
  full_data$zone_morpho <- "Outside bar"
  
  # Step 1: points on the sand bar
  for (poly in bar_polygons) {
    idx <- full_data$campaign == poly$campaign
    pts <- full_data[idx, ]
    if (nrow(pts) == 0) next
    inside <- point.in.polygon(pts$x, pts$y, poly$x, poly$y) > 0
    full_data$zone_morpho[idx][inside] <- "Sand bar"
  }
  
  # Step 2: among the remaining points, those in the wake zone
  for (poly in wake_polygons) {
    idx <- full_data$campaign == poly$campaign & full_data$zone_morpho != "Sand bar"
    pts <- full_data[idx, ]
    if (nrow(pts) == 0) next
    inside <- point.in.polygon(pts$x, pts$y, poly$x, poly$y) > 0
    full_data$zone_morpho[idx][inside] <- "Wake zone"
  }
  
  return(full_data)
}


# -----------------------------------------------------------------------------
# 6. SUMMARY STATISTICS PER MORPHOLOGICAL ZONE
# -----------------------------------------------------------------------------

ZONES_3 <- c("Sand bar", "Wake zone", "Outside bar")

#' Per campaign x zone statistics on bedload flux (Qb), optionally including
#' bathymetry means when available. Equivalent to "stats_zones".
compute_zone_statistics <- function(data) {
  
  has_bathy <- all(c("altitude", "depth") %in% colnames(data))
  
  stats <- data.frame()
  
  for (camp in unique(data$campaign)) {
    for (zone in ZONES_3) {
      idx <- data$campaign == camp & data$zone_morpho == zone
      if (sum(idx) == 0) next
      
      row <- data.frame(
        campaign    = camp,
        zone_morpho = zone,
        discharge   = data$discharge[idx][1],
        n_points    = sum(idx),
        Qb_mean     = mean(data$Qb[idx], na.rm = TRUE),
        Qb_median   = median(data$Qb[idx], na.rm = TRUE),
        Qb_sd       = sd(data$Qb[idx], na.rm = TRUE),
        Qb_min      = min(data$Qb[idx], na.rm = TRUE),
        Qb_max      = max(data$Qb[idx], na.rm = TRUE)
      )
      
      if (has_bathy) {
        row$altitude_mean <- mean(data$altitude[idx], na.rm = TRUE)
        row$depth_mean    <- mean(data$depth[idx], na.rm = TRUE)
      }
      
      stats <- rbind(stats, row)
    }
  }
  
  return(stats)
}


#' Pairwise differences in mean bedload flux between the three zones,
#' for each campaign. Equivalent to "differences_flux_3_zones".
compute_zone_flux_differences <- function(zone_stats) {
  
  diffs <- data.frame()
  
  for (camp in unique(zone_stats$campaign)) {
    camp_stats <- zone_stats[zone_stats$campaign == camp, ]
    
    get_qb <- function(zone) {
      if (zone %in% camp_stats$zone_morpho) {
        camp_stats$Qb_mean[camp_stats$zone_morpho == zone]
      } else {
        NA
      }
    }
    
    qb_bar  <- get_qb("Sand bar")
    qb_wake <- get_qb("Wake zone")
    qb_out  <- get_qb("Outside bar")
    
    diffs <- rbind(diffs, data.frame(
      campaign         = camp,
      qb_bar           = qb_bar,
      qb_wake          = qb_wake,
      qb_outside       = qb_out,
      diff_bar_outside = ifelse(!is.na(qb_bar) & !is.na(qb_out),  qb_bar - qb_out,  NA),
      diff_bar_wake    = ifelse(!is.na(qb_bar) & !is.na(qb_wake), qb_bar - qb_wake, NA),
      diff_wake_outside = ifelse(!is.na(qb_wake) & !is.na(qb_out), qb_wake - qb_out, NA)
    ))
  }
  
  return(diffs)
}

#' Median Qb per campaign x zone (used e.g. as a reference value for each
#' distribution). Equivalent to "medianes_zones".
compute_zone_medians <- function(data) {
  
  medians <- data.frame()
  
  for (camp in unique(data$campaign)) {
    for (zone in ZONES_3) {
      idx <- data$campaign == camp & data$zone_morpho == zone
      if (sum(idx) == 0) next
      
      medians <- rbind(medians, data.frame(
        campaign    = camp,
        zone_morpho = zone,
        median_Qb   = median(data$Qb[idx], na.rm = TRUE)
      ))
    }
  }
  
  return(medians)
}

#' Share of total bedload flux contributed by each zone, for each campaign.
#' Equivalent to "contributions_flux_par_zone".
compute_flux_contribution <- function(data) {
  
  contribution <- data.frame()
  
  for (camp in unique(data$campaign)) {
    camp_zones <- data %>%
      filter(campaign == camp) %>%
      group_by(zone_morpho) %>%
      summarise(
        campaign    = first(campaign),
        discharge   = first(discharge),
        date        = first(date),
        n_points    = n(),
        flux_total  = sum(Qb, na.rm = TRUE),
        flux_median = median(Qb, na.rm = TRUE),
        flux_mean   = mean(Qb, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        flux_total_campaign = sum(flux_total, na.rm = TRUE),
        pct_contribution    = 100 * flux_total / flux_total_campaign
      )
    
    contribution <- rbind(contribution, camp_zones)
  }
  
  contribution <- contribution %>% arrange(discharge)
  return(contribution)
}

#' Percentile statistics (p10 to p99) of Qb per campaign x zone.
#' Equivalent to "stats_percentiles_zones".
compute_percentile_statistics <- function(data) {
  
  stats <- data.frame()
  
  for (camp in unique(data$campaign)) {
    for (zone in ZONES_3) {
      idx <- data$campaign == camp & data$zone_morpho == zone
      if (sum(idx) == 0) next
      
      vals <- data$Qb[idx]
      stats <- rbind(stats, data.frame(
        campaign  = camp,
        zone_morpho = zone,
        discharge = data$discharge[idx][1],
        n_points  = sum(idx),
        p10 = quantile(vals, 0.10, na.rm = TRUE),
        p25 = quantile(vals, 0.25, na.rm = TRUE),
        p50 = quantile(vals, 0.50, na.rm = TRUE),
        p75 = quantile(vals, 0.75, na.rm = TRUE),
        p90 = quantile(vals, 0.90, na.rm = TRUE),
        p95 = quantile(vals, 0.95, na.rm = TRUE),
        p99 = quantile(vals, 0.99, na.rm = TRUE),
        mean = mean(vals, na.rm = TRUE)
      ))
    }
  }
  
  return(stats)
}

#' Ratio of the 90th percentile (and of the mean) to the median, used to
#' characterize how "active"/skewed a flux distribution is.
compute_percentile_comparison <- function(percentile_stats) {
  
  percentile_stats %>%
    transmute(
      campaign, zone_morpho, discharge,
      median = p50,
      p90,
      mean,
      ratio_p90_median  = p90 / p50,
      ratio_mean_median = mean / p50
    )
}


# -----------------------------------------------------------------------------
# 7. PROXIMITY INDEX ANALYSIS (subset with bathymetry)
# -----------------------------------------------------------------------------

#' Median proximity index of the sand bar only, per campaign.
#' Equivalent to "proximite_mediane_banc".
compute_bar_proximity <- function(data) {
  
  result <- data.frame()
  for (camp in unique(data$campaign)) {
    idx <- data$campaign == camp & data$zone_morpho == "Sand bar"
    if (sum(idx) == 0) next
    
    result <- rbind(result, data.frame(
      campaign = camp,
      median_bar_proximity = median(data$proximity_index[idx], na.rm = TRUE)
    ))
  }
  return(result)
}

#' Median proximity index over all zones, per campaign.
#' Equivalent to "proximite_mediane_globale".
compute_global_proximity <- function(data) {
  
  result <- data.frame()
  for (camp in unique(data$campaign)) {
    idx <- data$campaign == camp
    if (sum(idx) == 0) next
    
    result <- rbind(result, data.frame(
      campaign = camp,
      median_global_proximity = median(data$proximity_index[idx], na.rm = TRUE)
    ))
  }
  return(result)
}

#' Combine the bar-only and global proximity indices into a single table,
#' alongside discharge. Equivalent to "synthese_proximite_toutes_methodes".
compute_proximity_synthesis <- function(data, bar_proximity, global_proximity) {
  
  synthesis <- data.frame()
  
  for (camp in unique(data$campaign)) {
    discharge_camp <- data$discharge[data$campaign == camp][1]
    
    prox_bar <- if (camp %in% bar_proximity$campaign) {
      bar_proximity$median_bar_proximity[bar_proximity$campaign == camp]
    } else NA
    
    prox_global <- if (camp %in% global_proximity$campaign) {
      global_proximity$median_global_proximity[global_proximity$campaign == camp]
    } else NA
    
    synthesis <- rbind(synthesis, data.frame(
      campaign = camp,
      discharge = discharge_camp,
      bar_proximity = prox_bar,
      global_proximity = prox_global
    ))
  }
  
  return(synthesis)
}


# -----------------------------------------------------------------------------
# 8. DISCHARGE-LEVEL SUMMARIES (acoustic power & width-to-depth ratio)
# -----------------------------------------------------------------------------

#' Median/mean acoustic power (P, PdB) and bedload flux per campaign,
#' across all zones combined. Equivalent to "stats_globales_Pac_debit".
compute_global_power_statistics <- function(data) {
  
  data %>%
    group_by(campaign, discharge) %>%
    summarise(
      n_points   = n(),
      P_median   = median(P, na.rm = TRUE),
      P_mean     = mean(P, na.rm = TRUE),
      PdB_median = median(PdB, na.rm = TRUE),
      PdB_mean   = mean(PdB, na.rm = TRUE),
      Qb_median  = median(Qb, na.rm = TRUE),
      Qb_mean    = mean(Qb, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(discharge)
}

#' Same as above, but broken down by morphological zone.
#' Equivalent to "stats_PdB_par_zone_debit".
compute_zone_power_statistics <- function(data) {
  
  data %>%
    group_by(campaign, discharge, zone_morpho) %>%
    summarise(
      n_points   = n(),
      PdB_median = median(PdB, na.rm = TRUE),
      P_median   = median(P, na.rm = TRUE),
      Qb_median  = median(Qb, na.rm = TRUE),
      Qb_mean    = mean(Qb, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(discharge, zone_morpho)
}

#' Width-to-depth ratio per campaign, using a fixed channel width and the
#' median depth observed across the (bathymetry) subset.
#' Equivalent to "wd_ratio_global".
compute_wd_ratio <- function(data_with_bathy, w_fixed = W_FIXED) {
  
  data_with_bathy %>%
    group_by(campaign, discharge) %>%
    summarise(
      n_points = n(),
      D_median = median(abs(depth), na.rm = TRUE),
      D_mean   = mean(abs(depth), na.rm = TRUE),
      WD_ratio = w_fixed / D_median,
      .groups = "drop"
    ) %>%
    arrange(discharge)
}

#' Difference in mean Qb between the sand bar and the low-flow channel,
#' per campaign (used to characterize at which discharge the bar becomes
#' more/less active than the channel).
compute_bar_channel_flux_gap <- function(zone_power_or_flux_stats) {
  
  zone_power_or_flux_stats %>%
    filter(zone_morpho %in% c("Sand bar", "Outside bar")) %>%
    select(campaign, discharge, zone_morpho, Qb_mean) %>%
    pivot_wider(names_from = zone_morpho, values_from = Qb_mean) %>%
    mutate(flux_gap_bar_minus_channel = `Sand bar` - `Outside bar`) %>%
    filter(!is.na(flux_gap_bar_minus_channel))
}


# =============================================================================
# 9. EXAMPLE 
#
# The block below shows how the functions above are meant to be chained.
# Replace the placeholder file paths with your own raw data files, and
# the polygon coordinates with the ones digitized for your own surveys
# (one closed polygon per campaign, in the same projected coordinate
# system as the x/y columns of the acoustic and bathymetric files).
# =============================================================================

# ---- 9.1 Campaign definitions (EDIT FILE PATHS) ---------------------------
# File paths are built from DATA_DIR (section 2) + the original file names.
# Adjust DATA_DIR to point to wherever you store these files locally.

campaigns <- list(
  list(name = "240403", date = "2024-04-03", discharge = 2520,
       acoustic_file = file.path(DATA_DIR, "2024-04-03_traces_derives.txt"),
       bathy_file    = file.path(DATA_DIR, "TXYZ_240403.txt")),
  list(name = "241020", date = "2024-10-20", discharge = 1610,
       acoustic_file = file.path(DATA_DIR, "2024-10-21_traces_derives.txt"),
       bathy_file    = file.path(DATA_DIR, "TXYZ_241021.txt")),
  list(name = "250120", date = "2025-01-20", discharge = 1270,
       acoustic_file = file.path(DATA_DIR, "2025-01-20_traces_derives.txt"),
       bathy_file    = file.path(DATA_DIR, "TXYZ_250120.txt")),
  list(name = "250131", date = "2025-01-31", discharge = 1840,
       acoustic_file = file.path(DATA_DIR, "2025-01-31_traces_derives.txt"),
       bathy_file    = file.path(DATA_DIR, "TXYZ_250131.txt")),
  list(name = "250328", date = "2025-03-28", discharge = 937,
       acoustic_file = file.path(DATA_DIR, "2025-03-28_traces_derives.txt"),
       bathy_file    = file.path(DATA_DIR, "TXYZ_250328.txt")),
  list(name = "250528", date = "2025-05-28", discharge = 343,
       acoustic_file = file.path(DATA_DIR, "2025-05-28_traces_derives.txt"),
       bathy_file    = file.path(DATA_DIR, "TXYZ_250528.txt")),
  list(name = "260217", date = "2026-02-17", discharge = 4020,
       acoustic_file = file.path(DATA_DIR, "2026-02-17_traces_derives.txt"),
       bathy_file    = file.path(DATA_DIR, "TXYZ_260217.txt"))
)

# ---- 9.2 Morphological zone polygons ---------------------------------------
# One digitized polygon per campaign (sand bar and wake zone), in the same
# projected coordinate system as the x/y columns of the acoustic and
# bathymetric files.

bar_polygons <- list(
  list(campaign = "240403",
       x = c(446074.9436, 446228.8627, 446436.4792, 446357.2781, 446270.6052, 446132.8421,
             446061.3948, 446041.6861, 446047.4944, 446074.9436),
       y = c(6707266.573, 6707230.709, 6707156.592, 6706853.237, 6706878.641, 6706952.697,
             6707056.47, 6707118.571, 6707216.673, 6707266.573)),
  list(campaign = "241020",
       x = c(446967.0767, 447010.4131, 447246.4225, 447148.9304, 446950.7185, 446916.4846,
             446902.7196, 446902.4376, 446939.6274, 446967.0767),
       y = c(6707003.566, 6707025.981, 6706939.91, 6706666.502, 6706740.144, 6706764.049,
             6706824.844, 6706874.991, 6706973.092, 6707003.566)),
  list(campaign = "250120",
       x = c(446733.9565, 446826.6068, 447073.0767, 446984.9094, 446780.1821, 446718.6313,
             446674.0825, 446672.3061, 446712.4847, 446733.9565),
       y = c(6707057.363, 6707040.925, 6706971.291, 6706682.88, 6706763.575, 6706815.216,
             6706883.124, 6706949.709, 6707037.35, 6707057.363)),
  list(campaign = "250131",
       x = c(446779.8632, 446831.329, 447086.5259, 446992.3812, 446784.7847, 446720.2452,
             446667.328, 446685.7554, 446661.796, 446679.011, 446779.8632),
       y = c(6707016.955, 6707055.57, 6706969.797, 6706687.363, 6706753.055, 6706790.111,
             6706847.439, 6706928.788, 6707023.841, 6707073.573, 6707016.955)),
  list(campaign = "250328",
       x = c(446772.2121, 446721.3441, 446736.487, 446943.6053, 446841.2117, 446646.6031,
             446533.4333, 446515.5179, 446567.1133, 446611.3592, 446747.9637, 446772.2121),
       y = c(6707043.734, 6707071.828, 6707111.343, 6707029.751, 6706740.622, 6706835.061,
             6706927.776, 6706994.779, 6707074.529, 6707051.744, 6706980.317, 6707043.734)),
  list(campaign = "250528",
       x = c(446772.2121, 446927.3467, 446868.9471, 446672.4257, 446565.9506, 446537.5148,
             446561.4675, 446772.2121),
       y = c(6707043.734, 6706984.8, 6706811.395, 6706885.75, 6706931.601, 6706999.561,
             6707111.343, 6707043.734)),
  list(campaign = "260217",
       x = c(446455.5023, 446483.8847, 446601.1684, 446635.757, 446813.3096, 446755.6007, 446693.8278,
             446544.2723, 446354.8895, 446211.8364, 446210.71, 446215.3824, 446259.0488, 446345.8224, 446455.5023),
       y = c(6707139.3987, 6707134.4192, 6707113.843, 6707079.4353, 6707020.3708, 6706800.1016, 6706791.1608,
             6706795.224, 6706827.7369, 6706862.6873, 6706884.41, 6707010.414, 6707127.4924, 6707132.499, 6707139.3987))
)

wake_polygons <- list(
  list(campaign = "240403",
       x = c(445961.3722, 446077.9323, 446056.9117, 446013.5753, 446037.485, 445999.844,
             445835.7464, 445860.8685, 445874.1486, 445961.3722),
       y = c(6707345.774, 6707293.472, 6707261.197, 6707204.412, 6707171.536, 6707103.628,
             6707199.928, 6707285.94, 6707336.222, 6707345.774)),
  list(campaign = "241020",
       x = c(446852.2377, 446782.4488, 446625.5709, 446656.3264, 446790.8489, 446992.8918,
             446934.0691, 446904.2921, 446898.6854, 446887.6241, 446852.2377),
       y = c(6706768.493, 6706875.186, 6707003.109, 6707117.6, 6707097.626, 6707032.041,
             6706966.995, 6706892.108, 6706832.93, 6706756.06, 6706768.493)),
  list(campaign = "250120",
       x = c(446680.1595, 446849.0222, 446731.4777, 446699.6318, 446672.4142, 446664.5952,
             446541.0844, 446386.8833, 446437.5224, 446680.1595),
       y = c(6707133.575, 6707046.902, 6707060.416, 6707020.46, 6706963.036, 6706924.962,
             6706957.842, 6707022.932, 6707147.932, 6707133.575)),
  list(campaign = "250131",
       x = c(446775.6431, 446680.7174, 446665.4899, 446632.2404, 446578.1614, 446515.6803,
             446450.5191, 446367.4118, 446420.9752, 446694.809, 446829.9143, 446775.6431),
       y = c(6707016.377, 6707077.212, 6707030.599, 6706989.224, 6706978.102, 6706996.696,
             6706964.227, 6707026.232, 6707164.809, 6707123.722, 6707059.997, 6707016.377)),
  list(campaign = "250328",
       x = c(446643.0994, 446732.7744, 446720.3397, 446771.5015, 446749.7125, 446564.3535,
             446510.4799, 446306.0684, 446335.6668, 446347.3954, 446547.1217, 446643.0994),
       y = c(6707130.766, 6707102.73, 6707071.89, 6707043.771, 6706986.992, 6707077.028,
             6706997.592, 6707081.81, 6707139.564, 6707189.464, 6707144.817, 6707130.766)),
  list(campaign = "250528",
       x = c(446526.4197, 446548.616, 446549.5127, 446503.1446, 446484.6574, 446425.2551,
             446434.2676, 446526.4197),
       y = c(6707082.946, 6707057.486, 6706954.853, 6706979.476, 6707024.371, 6707068.977,
             6707113.256, 6707082.946)),
  list(campaign = "260217",
       x = c(446813.3096, 446635.757, 446601.1684, 446455.5023, 446259.0488, 446215.3924, 446211.8363, 446137.0587,
             446100.2985, 446165.783, 446253.0956, 446372.1584, 446747.206, 446880.1594, 446813.3096),
       y = c(6707020.3708, 6707079.4353, 6707113.843, 6707139.3987, 6707127.4924, 6707010.414, 6706862.6873,
             6706819.6088, 6706911.1951, 6707133.4455, 6707189.0081, 6707185.0394, 6707121.5392, 6707069.9454, 6707020.3708))
)

# ---- 9.3 Run the processing --------------------------------------------------
#
# Uncomment to execute once the configuration above has been filled in.
#
# full_data        <- prepare_full_dataset(campaigns)
# full_data         <- classify_morphological_zones(full_data, bar_polygons, wake_polygons)
#
# data_with_bathy   <- associate_bathymetry(full_data, campaigns)
#
# zone_stats        <- compute_zone_statistics(full_data)
# flux_differences   <- compute_zone_flux_differences(zone_stats)
# zone_medians       <- compute_zone_medians(full_data)
# flux_contribution  <- compute_flux_contribution(full_data)
# percentile_stats   <- compute_percentile_statistics(full_data)
# percentile_compare <- compute_percentile_comparison(percentile_stats)
#
# bar_proximity      <- compute_bar_proximity(data_with_bathy)
# global_proximity   <- compute_global_proximity(data_with_bathy)
# proximity_synthesis <- compute_proximity_synthesis(data_with_bathy, bar_proximity, global_proximity)
#
# global_power_stats <- compute_global_power_statistics(full_data)
# zone_power_stats   <- compute_zone_power_statistics(full_data)
# wd_ratio           <- compute_wd_ratio(data_with_bathy)
# flux_gap           <- compute_bar_channel_flux_gap(zone_power_stats)
#

# -----------------------------------------------------------------------------
# 10. EXPORTING RESULTS TO CSV (OPTIONAL)
# -----------------------------------------------------------------------------

#' Write a named list of data frames to individual CSV files in `output_dir`.
#' @param results named list, e.g. list(stats_zones = zone_stats, ...)
export_results <- function(results, output_dir = OUTPUT_DIR) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  for (name in names(results)) {
    if (is.null(results[[name]])) next
    file_path <- file.path(output_dir, paste0(name, ".csv"))
    write.csv(results[[name]], file_path, row.names = FALSE)
    message("Exported: ", file_path)
  }
}

# ---- 10.1 Example export call ----------------------------------------------
#
# Uncomment once the pipeline (section 9.3) has been run.
#
# export_results(list(
#   donnees_morpho                     = full_data,
#   donnees_avec_rugosity              = data_with_bathy,
#   stats_zones                        = zone_stats,
#   differences_flux_3_zones           = flux_differences,
#   medianes_zones                     = zone_medians,
#   contributions_flux_par_zone        = flux_contribution,
#   stats_percentiles_zones            = percentile_stats,
#   percentile_comparison              = percentile_compare,
#   proximite_mediane_banc             = bar_proximity,
#   proximite_mediane_globale          = global_proximity,
#   synthese_proximite_toutes_methodes = proximity_synthesis,
#   stats_globales_Pac_debit           = global_power_stats,
#   stats_PdB_par_zone_debit           = zone_power_stats,
#   wd_ratio_global                    = wd_ratio,
#   flux_gap_bar_channel               = flux_gap
# ))