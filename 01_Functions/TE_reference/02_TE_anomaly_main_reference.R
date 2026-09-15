# Author: Zhaozhe Chen
# Anomaly-based disturbance-scenario analysis: 2026.8.13
#
# This workflow follows the earlier project logic: calculate a centered
# five-day moving mean, subtract it from each variable, and apply TE to the
# resulting anomalies. Anomalies are calculated independently within the Pre
# and Post phases, so no moving window crosses the disturbance or transition.


# ------- Global -------
old_skip_main <- getOption("disturbance_synchrony.skip_main")
options(disturbance_synchrony.skip_main = TRUE)
source("02_Analysis/01_TE_main.R")
options(disturbance_synchrony.skip_main = old_skip_main)


# Calculate anomalies using the previous centered moving-window logic ========
calculate_TE_anomaly <- function(
    Phase_df,
    phase_name = NULL,
    scenario = NULL,
    window_size = 24 * 5) {
  if (nrow(Phase_df) < window_size) {
    stop("Each phase must contain at least window_size observations.")
  }

  Phase_df <- Phase_df[order(Phase_df$time), ]
  centered_partial_mean <- function(x, width) {
    n <- length(x)
    left_width <- floor((width - 1) / 2)
    right_width <- width - 1 - left_width
    left_idx <- pmax(1L, seq_len(n) - left_width)
    right_idx <- pmin(n, seq_len(n) + right_width)
    x_sum <- ifelse(is.na(x), 0, x)
    x_count <- as.integer(!is.na(x))
    cumulative_sum <- c(0, cumsum(x_sum))
    cumulative_count <- c(0L, cumsum(x_count))
    window_sum <- cumulative_sum[right_idx + 1L] - cumulative_sum[left_idx]
    window_count <- cumulative_count[right_idx + 1L] -
      cumulative_count[left_idx]
    ifelse(window_count > 0, window_sum / window_count, NA_real_)
  }

  source_smooth <- centered_partial_mean(
    Phase_df$source_variable, window_size
  )
  sink_smooth <- centered_partial_mean(
    Phase_df$sink_variable, window_size
  )

  Phase_df$source_variable <- Phase_df$source_variable - source_smooth
  Phase_df$sink_variable <- Phase_df$sink_variable - sink_smooth
  rownames(Phase_df) <- NULL
  Phase_df
}


# Run anomaly TE for all supplied disturbance scenarios ======================
run_anomaly_disturbance_analysis <- function() {
  anomaly_parameters <- list(
    n_bin = 11,
    max_lag = 72,
    Lag_Dependent_Crit = FALSE,
    nshuffle = 300,
    alpha = 0.05,
    ZFlagSink = FALSE,
    ZFlagSource = FALSE,
    lower_qt = 0.001,
    upper_qt = 0.999,
    seed = 111,
    parallel = FALSE,
    time_unit = "hours",
    output_figures = TRUE,
    output_data = TRUE,
    figure_dpi = 300
  )

  run_TE_main(
    scenarios = get_disturbance_scenarios(),
    Output_path = "04_Results/Anomaly",
    Report_file = "03_Reports/TE_anomaly_disturbance_report.html",
    report_title = "Transfer Entropy of Anomalies Across Disturbance Scenarios",
    phase_transform = calculate_TE_anomaly,
    phase_source_varname = "climate anomaly",
    phase_sink_varname = "log(response) anomaly",
    overlay_filename_suffix = "_Anomaly",
    analysis_note = paste(
      "For each Pre- and Post-disturbance phase separately, source and sink",
      "anomalies equal the observed value minus its centered 120-hour",
      "(five-day) moving mean. Transition observations are not analyzed.",
      "Zero-reserved binning is disabled because anomaly zero is not a",
      "structural zero."
    ),
    shared_parameters = anomaly_parameters
  )
}


# Reproduce the complete anomaly analysis after sourcing this file:
anomaly_results <- run_anomaly_disturbance_analysis()
