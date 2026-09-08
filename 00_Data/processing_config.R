# Central configuration for the ET Synchrony 2.0 data-processing workflow.
# Paths to the original project and shared data are inputs only.

cfg <- list(
  reference_project = "D:/OneDrive - UW-Madison/Research/ET Synchrony/Github repo/ET_Synchrony",
  raw_data_dir = "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/01_AMF_raw",
  hourly_input_dir = "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/02_AMF_cleaned/AMF_Hourly/AMF_sites_hourly_update",
  site_info_file = "D:/OneDrive - UW-Madison/Research/ET Synchrony/Github repo/ET_Synchrony/00_Data/ameriflux_site_info_update.csv",
  latent_heat_j_kg = 2.45e6,
  min_coverage = 0.75,
  subdaily_window_days = 5L,
  seasonal_harmonics = list(daily = 3L, weekly = 2L, monthly = 2L),
  minimum_anomaly_observations = list(daily = 30L, weekly = 20L, monthly = 18L),
  te_readiness_thresholds = c(30L, 100L, 300L, 1000L),
  scales = c("hourly", "6_hourly", "daily", "weekly", "monthly")
)
