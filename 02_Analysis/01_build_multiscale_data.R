#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- if (length(file_arg)) normalizePath(file_arg[1], winslash = "/") else normalizePath("02_Analysis/01_build_multiscale_data.R", winslash = "/")
project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
source(file.path(project_root, "00_Data", "processing_config.R"))
source(file.path(project_root, "01_Functions", "data_processing_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
site_arg <- args[grepl("^--sites=", args)]
requested_sites <- if (length(site_arg)) strsplit(sub("^--sites=", "", site_arg[1]), ",")[[1]] else character()
worker_arg <- args[grepl("^--workers=", args)]
workers <- if (length(worker_arg)) as.integer(sub("^--workers=", "", worker_arg[1])) else 4L
if (!is.finite(workers) || workers < 1L) workers <- 1L

input_files <- list.files(cfg$hourly_input_dir, pattern = "^AMF_hourly_.*\\.csv$", full.names = TRUE)
input_sites <- sub("^AMF_hourly_(.*)\\.csv$", "\\1", basename(input_files))
if (length(requested_sites)) {
  keep <- input_sites %in% requested_sites
  missing_sites <- setdiff(requested_sites, input_sites)
  if (length(missing_sites)) warning("Requested sites not found: ", paste(missing_sites, collapse = ", "))
  input_files <- input_files[keep]
  input_sites <- input_sites[keep]
}
if (!length(input_files)) stop("No site input files selected")

site_info <- read_site_info(cfg$site_info_file)
result_root <- file.path(project_root, "04_Results")
processed_root <- file.path(result_root, "Processed_Data")
dir.create(processed_root, recursive = TRUE, showWarnings = FALSE)

process_one_site <- function(i) {
  site_id <- input_sites[i]
  started <- Sys.time()
  message(sprintf("[%d/%d] %s", i, length(input_files), site_id))
  status <- "ok"
  detail <- ""
  site_summaries <- list()
  tryCatch({
    hourly <- read_legacy_hourly_site(input_files[i], site_info, cfg$latent_heat_j_kg)
    hourly <- add_subdaily_changes(hourly, 1, cfg$subdaily_window_days)
    hourly$ET_total_mm <- hourly$ET_mm_day / 24
    hourly$expected_hours <- 1
    out_path <- file.path(processed_root, "hourly", paste0(site_id, "_hourly.csv"))
    write_csv_stable(hourly, out_path)
    site_summaries[[length(site_summaries) + 1L]] <- summarize_processed(hourly, site_id, "hourly")

    for (scale in c("6_hourly", "daily", "weekly", "monthly")) {
      scaled <- aggregate_hourly(hourly, scale, cfg$min_coverage)
      if (scale == "6_hourly") {
        scaled <- add_subdaily_changes(scaled, 6, cfg$subdaily_window_days)
      } else {
        scaled <- add_coarse_anomalies(
          scaled, scale,
          cfg$seasonal_harmonics[[scale]],
          cfg$minimum_anomaly_observations[[scale]]
        )
      }
      out_path <- file.path(processed_root, scale, paste0(site_id, "_", scale, ".csv"))
      write_csv_stable(scaled, out_path)
      site_summaries[[length(site_summaries) + 1L]] <- summarize_processed(scaled, site_id, scale)
    }
  }, error = function(e) {
    status <<- "error"
    detail <<- conditionMessage(e)
    warning(site_id, ": ", detail)
  })
  site_log <- data.frame(
    Site_ID = site_id,
    input_file = normalizePath(input_files[i], winslash = "/", mustWork = FALSE),
    status = status,
    detail = detail,
    elapsed_seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2),
    processed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
  list(summary = if (length(site_summaries)) do.call(rbind, site_summaries) else NULL,
       log = site_log)
}

workers <- min(workers, length(input_files))
if (workers > 1L) {
  message("Using ", workers, " parallel workers")
  cluster <- parallel::makeCluster(workers)
  export_names <- setdiff(ls(envir = .GlobalEnv), "cluster")
  parallel::clusterExport(cluster, export_names, envir = .GlobalEnv)
  results <- parallel::parLapplyLB(cluster, seq_along(input_files), process_one_site)
  parallel::stopCluster(cluster)
} else {
  results <- lapply(seq_along(input_files), process_one_site)
}

summaries <- lapply(results, `[[`, "summary")
summaries <- summaries[!vapply(summaries, is.null, logical(1))]
summary_df <- if (length(summaries)) do.call(rbind, summaries) else data.frame()
log_df <- do.call(rbind, lapply(results, `[[`, "log"))
write_csv_stable(summary_df, file.path(result_root, "Tables", "data_availability_summary.csv"))
write_csv_stable(log_df, file.path(result_root, "Tables", "processing_log.csv"))

manifest <- data.frame(
  parameter = c("reference_project", "raw_data_dir", "hourly_input_dir", "site_info_file",
                "latent_heat_j_kg", "min_coverage", "subdaily_window_days",
                "processed_sites", "successful_sites", "failed_sites"),
  value = c(cfg$reference_project, cfg$raw_data_dir, cfg$hourly_input_dir, cfg$site_info_file,
            cfg$latent_heat_j_kg, cfg$min_coverage, cfg$subdaily_window_days,
            length(input_files), sum(log_df$status == "ok"), sum(log_df$status != "ok")),
  stringsAsFactors = FALSE
)
write_csv_stable(manifest, file.path(result_root, "Tables", "processing_manifest.csv"))
message("Processing complete: ", sum(log_df$status == "ok"), " succeeded; ",
        sum(log_df$status != "ok"), " failed")
