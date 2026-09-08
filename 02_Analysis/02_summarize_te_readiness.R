#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- if (length(file_arg)) normalizePath(file_arg[1], winslash = "/") else normalizePath("02_Analysis/02_summarize_te_readiness.R", winslash = "/")
project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
source(file.path(project_root, "00_Data", "processing_config.R"))

args <- commandArgs(trailingOnly = TRUE)
worker_arg <- args[grepl("^--workers=", args)]
workers <- if (length(worker_arg)) as.integer(sub("^--workers=", "", worker_arg[1])) else 4L
if (!is.finite(workers) || workers < 1L) workers <- 1L

processed_root <- file.path(project_root, "04_Results", "Processed_Data")
table_dir <- file.path(project_root, "04_Results", "Tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

comparisons <- list(
  "ET + soil water potential" = c("analysis_ET_mm_day", "analysis_log10_psi_soil"),
  "ET + VPD" = c("analysis_ET_mm_day", "analysis_VPD"),
  "ET + air temperature" = c("analysis_ET_mm_day", "analysis_TA"),
  "ET + net radiation" = c("analysis_ET_mm_day", "analysis_NETRAD"),
  "ET + all core drivers" = c("analysis_ET_mm_day", "analysis_log10_psi_soil", "analysis_VPD", "analysis_TA"),
  "ET + all drivers including net radiation" = c("analysis_ET_mm_day", "analysis_log10_psi_soil", "analysis_VPD", "analysis_TA", "analysis_NETRAD")
)

longest_true_run <- function(x) {
  if (!length(x) || !any(x)) return(0L)
  runs <- rle(x)
  max(runs$lengths[runs$values])
}

summarize_file <- function(path, scale) {
  d <- read.csv(path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""), check.names = FALSE)
  absent <- setdiff(unique(unlist(comparisons)), names(d))
  if (length(absent)) stop(basename(path), " lacks ", paste(absent, collapse = ", "))
  time <- as.POSIXct(d$Time, tz = "UTC")
  site_id <- sub(paste0("_", scale, "\\.csv$"), "", basename(path))
  rows <- lapply(names(comparisons), function(label) {
    cols <- comparisons[[label]]
    good <- Reduce(`&`, lapply(d[cols], is.finite))
    where <- which(good)
    first <- if (length(where)) format(time[min(where)], tz = "UTC") else NA_character_
    last <- if (length(where)) format(time[max(where)], tz = "UTC") else NA_character_
    span_steps <- if (length(where)) max(where) - min(where) + 1L else 0L
    runs <- if (any(good)) rle(good) else list(values = logical(), lengths = integer())
    data.frame(
      Site_ID = site_id,
      scale = scale,
      variable_set = label,
      n_time_steps = nrow(d),
      n_simultaneous = sum(good),
      fraction_of_all_steps = if (nrow(d)) mean(good) else NA_real_,
      overlap_start = first,
      overlap_end = last,
      overlap_span_steps = span_steps,
      fraction_complete_within_span = if (span_steps) sum(good) / span_steps else NA_real_,
      longest_contiguous_run = longest_true_run(good),
      n_contiguous_segments = if (length(runs$values)) sum(runs$values) else 0L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

jobs <- do.call(rbind, lapply(cfg$scales, function(scale) {
  files <- list.files(file.path(processed_root, scale), pattern = "\\.csv$", full.names = TRUE)
  data.frame(path = files, scale = scale, stringsAsFactors = FALSE)
}))
if (!nrow(jobs)) stop("No processed site files found")

run_job <- function(i) summarize_file(jobs$path[i], jobs$scale[i])
workers <- min(workers, nrow(jobs))
if (workers > 1L) {
  message("Calculating exact temporal overlap with ", workers, " workers")
  cluster <- parallel::makeCluster(workers)
  parallel::clusterExport(cluster,
                          c("jobs", "comparisons", "summarize_file", "longest_true_run"),
                          envir = .GlobalEnv)
  pieces <- parallel::parLapplyLB(cluster, seq_len(nrow(jobs)), run_job)
  parallel::stopCluster(cluster)
} else {
  pieces <- lapply(seq_len(nrow(jobs)), run_job)
}
overlap <- do.call(rbind, pieces)
write.csv(overlap, file.path(table_dir, "te_input_overlap_by_site.csv"),
          row.names = FALSE, na = "NA")

groups <- split(seq_len(nrow(overlap)), interaction(overlap$scale, overlap$variable_set, drop = TRUE))
cross_site <- do.call(rbind, lapply(groups, function(i) {
  d <- overlap[i, , drop = FALSE]
  out <- data.frame(
    scale = d$scale[1],
    variable_set = d$variable_set[1],
    n_sites_total = nrow(d),
    n_sites_with_overlap = sum(d$n_simultaneous > 0),
    median_simultaneous = median(d$n_simultaneous),
    q25_simultaneous = unname(quantile(d$n_simultaneous, 0.25)),
    q75_simultaneous = unname(quantile(d$n_simultaneous, 0.75)),
    median_longest_run = median(d$longest_contiguous_run),
    stringsAsFactors = FALSE
  )
  for (threshold in cfg$te_readiness_thresholds) {
    out[[paste0("n_sites_ge_", threshold)]] <- sum(d$n_simultaneous >= threshold)
  }
  out
}))
cross_site$scale <- factor(cross_site$scale, levels = cfg$scales)
cross_site <- cross_site[order(cross_site$variable_set, cross_site$scale), ]
cross_site$scale <- as.character(cross_site$scale)
write.csv(cross_site, file.path(table_dir, "te_input_overlap_cross_site_summary.csv"),
          row.names = FALSE, na = "NA")
message("Overlap summaries written: ", nrow(overlap), " site-scale-variable-set rows")

