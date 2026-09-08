#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- if (length(file_arg)) normalizePath(file_arg[1], winslash = "/") else normalizePath("02_Analysis/02_make_processing_report.R", winslash = "/")
project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
source(file.path(project_root, "00_Data", "processing_config.R"))

result_root <- file.path(project_root, "04_Results")
table_dir <- file.path(result_root, "Tables")
figure_dir <- file.path(result_root, "Figures")
report_dir <- file.path(project_root, "03_Reports")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

availability <- read.csv(file.path(table_dir, "data_availability_summary.csv"), stringsAsFactors = FALSE)
log_df <- read.csv(file.path(table_dir, "processing_log.csv"), stringsAsFactors = FALSE)
overlap <- read.csv(file.path(table_dir, "te_input_overlap_by_site.csv"), stringsAsFactors = FALSE, na.strings = "NA")
cross_site <- read.csv(file.path(table_dir, "te_input_overlap_cross_site_summary.csv"), stringsAsFactors = FALSE)
site_info <- read.csv(cfg$site_info_file, stringsAsFactors = FALSE, check.names = FALSE)
names(site_info)[names(site_info) == "site_id"] <- "Site_ID"

availability$scale <- factor(availability$scale, levels = cfg$scales)
overlap$scale <- factor(overlap$scale, levels = cfg$scales)
cross_site$scale <- factor(cross_site$scale, levels = cfg$scales)

# ColorBrewer Set2: five scales followed by supporting colors.
set2 <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
          "#FFD92F", "#E5C494", "#B3B3B3")
scale_cols <- setNames(set2[1:5], cfg$scales)
driver_labels <- c("ET + soil water potential", "ET + VPD",
                   "ET + air temperature", "ET + net radiation")

# Select a representative site based on median all-variable data availability.
all_driver <- overlap[overlap$variable_set == "ET + all drivers including net radiation", ]
wide <- reshape(all_driver[, c("Site_ID", "scale", "n_simultaneous")],
                idvar = "Site_ID", timevar = "scale", direction = "wide")
count_cols <- paste0("n_simultaneous.", cfg$scales)
eligible <- wide[complete.cases(wide[, count_cols]) & apply(wide[, count_cols] > 0, 1, all), ]
eligible <- merge(eligible, site_info[, c("Site_ID", "latitude", "longitude")], by = "Site_ID")
eligible <- eligible[eligible$longitude >= -125 & eligible$longitude <= -66 &
                     eligible$latitude >= 24 & eligible$latitude <= 50, ]
if (!nrow(eligible)) stop("No site has all-variable overlap at every scale")
profile <- log10(as.matrix(eligible[, count_cols]) + 1)
median_profile <- apply(profile, 2, median)
eligible$distance_to_median <- sqrt(rowSums(sweep(profile, 2, median_profile)^2))
representative <- eligible$Site_ID[which.min(eligible$distance_to_median)]
representative_reason <- eligible[eligible$Site_ID == representative,
                                  c("Site_ID", "latitude", "longitude", count_cols), drop = FALSE]

read_processed <- function(site_id, scale) {
  path <- file.path(result_root, "Processed_Data", scale,
                    paste0(site_id, "_", scale, ".csv"))
  d <- read.csv(path, stringsAsFactors = FALSE,
                na.strings = c("NA", "NaN", ""), check.names = FALSE)
  d$Time <- as.POSIXct(d$Time, tz = "UTC")
  d
}

# Minimal ESRI Polygon reader for the locally available state-boundary file.
read_polygon_shp <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  readBin(con, "raw", 100)
  polygons <- list()
  repeat {
    rec_header <- readBin(con, "integer", 2, size = 4, endian = "big", signed = TRUE)
    if (length(rec_header) < 2) break
    content <- readBin(con, "raw", rec_header[2] * 2)
    rc <- rawConnection(content, "rb")
    shape_type <- readBin(rc, "integer", 1, size = 4, endian = "little", signed = TRUE)
    if (shape_type %in% c(5L, 15L, 25L)) {
      readBin(rc, "double", 4, size = 8, endian = "little")
      n_parts <- readBin(rc, "integer", 1, size = 4, endian = "little", signed = TRUE)
      n_points <- readBin(rc, "integer", 1, size = 4, endian = "little", signed = TRUE)
      starts <- readBin(rc, "integer", n_parts, size = 4, endian = "little", signed = TRUE) + 1L
      coords <- matrix(readBin(rc, "double", n_points * 2L, size = 8, endian = "little"),
                       ncol = 2, byrow = TRUE)
      ends <- c(starts[-1] - 1L, n_points)
      for (j in seq_len(n_parts)) {
        polygons[[length(polygons) + 1L]] <- coords[starts[j]:ends[j], , drop = FALSE]
      }
    }
    close(rc)
  }
  polygons
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

table_html <- function(x, digits = 2) {
  if (!nrow(x)) return("<p>No records.</p>")
  cells <- lapply(x, function(z) {
    if (is.numeric(z)) format(round(z, digits), trim = TRUE) else html_escape(as.character(z))
  })
  x[] <- cells
  paste0("<table><thead><tr>", paste0("<th>", names(x), "</th>", collapse = ""),
         "</tr></thead><tbody>",
         paste(apply(x, 1, function(r) {
           paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>")
         }), collapse = ""), "</tbody></table>")
}

display_scale <- function(x) gsub("_", "-", x)
file_slug <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))

# Figure 1: site counts at explicit simultaneous-sample thresholds.
threshold_cols <- paste0("n_sites_ge_", cfg$te_readiness_thresholds)
threshold_data <- cross_site[cross_site$variable_set %in% driver_labels, ]
png(file.path(figure_dir, "sites_meeting_sample_thresholds.png"),
    width = 1900, height = 1500, res = 170)
par(mfrow = c(2, 2), mar = c(6, 5, 4, 1), oma = c(1, 0, 3, 0))
for (driver in driver_labels) {
  d <- threshold_data[threshold_data$variable_set == driver, ]
  d <- d[match(cfg$scales, as.character(d$scale)), ]
  mat <- t(as.matrix(d[, threshold_cols, drop = FALSE]))
  colnames(mat) <- display_scale(cfg$scales)
  rownames(mat) <- paste0("at least ", cfg$te_readiness_thresholds)
  barplot(mat, beside = TRUE, col = set2[1:nrow(mat)], border = NA,
          ylim = c(0, 170), ylab = "Number of sites", las = 2, main = driver)
  abline(h = c(30, 60, 90, 120, 150), col = "grey90", lty = 3)
  legend("topright", legend = rownames(mat), fill = set2[1:nrow(mat)],
         border = NA, bty = "n", cex = 0.78)
}
mtext("Sites meeting candidate simultaneous-sample thresholds", outer = TRUE,
      font = 2, cex = 1.25)
dev.off()

# Figure 2: across-site distributions; each histogram observation is one site.
distribution_figures <- character()
for (driver in driver_labels) {
  filename <- paste0("site_sample_distribution_", file_slug(driver), ".png")
  distribution_figures[driver] <- filename
  png(file.path(figure_dir, filename), width = 2000, height = 650, res = 170)
  par(mfrow = c(1, 5), mar = c(5, 4, 4, 1), oma = c(0, 0, 3, 0))
  for (scale in cfg$scales) {
    x <- overlap$n_simultaneous[overlap$variable_set == driver & overlap$scale == scale]
    upper <- max(log10(x + 1)) + 0.2
    hist(log10(x + 1), breaks = seq(0, upper, length.out = 18),
         col = scale_cols[scale], border = "white", main = display_scale(scale),
         xlab = "log10(simultaneous observations + 1)", ylab = "Number of sites")
    mtext(paste0("median = ", median(x), "; sites > 0 = ", sum(x > 0)),
          side = 3, line = 0.25, cex = 0.68)
  }
  mtext(paste0(driver, ": distribution across 163 sites"), outer = TRUE,
        font = 2, cex = 1.15)
  dev.off()
}

# Figure 3: pairwise completeness, with dark text for visibility.
pair_only <- overlap[overlap$variable_set %in% driver_labels, ]
med_complete <- aggregate(fraction_complete_within_span ~ variable_set + scale,
                          pair_only, median, na.rm = TRUE)
mat <- xtabs(fraction_complete_within_span ~ variable_set + scale, med_complete)
mat <- mat[driver_labels, cfg$scales, drop = FALSE]
png(file.path(figure_dir, "median_pairwise_completeness.png"),
    width = 1600, height = 950, res = 170)
par(mar = c(7, 12, 4, 5))
classes <- matrix(findInterval(mat, seq(0, 1, length.out = 9), all.inside = TRUE),
                  nrow = nrow(mat))
image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(classes), axes = FALSE,
      col = set2, zlim = c(1, 8), xlab = "", ylab = "")
axis(1, seq_len(ncol(mat)), display_scale(colnames(mat)), las = 2)
axis(2, seq_len(nrow(mat)), rownames(mat), las = 2)
for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat))) {
  text(j, i, sprintf("%.2f", mat[i, j]), cex = 1.05, col = "#1F2933", font = 2)
}
title("Median simultaneous completeness within each site's overlap span")
dev.off()

# Figure 4: geographic context and representative site.
state_shp <- file.path(cfg$reference_project, "00_Data", "cb_2018_us_state_20m",
                       "cb_2018_us_state_20m.shp")
states <- read_polygon_shp(state_shp)
map_sites <- merge(unique(overlap[, "Site_ID", drop = FALSE]),
                   site_info[, c("Site_ID", "latitude", "longitude")], by = "Site_ID")
monthly_core <- overlap[overlap$scale == "monthly" &
                        overlap$variable_set == "ET + all core drivers",
                        c("Site_ID", "n_simultaneous")]
map_sites <- merge(map_sites, monthly_core, by = "Site_ID", all.x = TRUE)
map_sites$class <- cut(map_sites$n_simultaneous, breaks = c(-Inf, 0, 29, 99, Inf),
                       labels = c("0", "1-29", "30-99", "100+"))
map_cols <- setNames(set2[c(8, 7, 5, 1)], levels(map_sites$class))
rep_row <- map_sites[map_sites$Site_ID == representative, ]
png(file.path(figure_dir, "site_map_representative.png"),
    width = 1700, height = 1050, res = 170)
par(mar = c(4, 4, 4, 1))
plot(NA, xlim = c(-125, -66.5), ylim = c(24, 50.5), asp = 1.25,
     xlab = "Longitude", ylab = "Latitude",
     main = "AmeriFlux sites and the representative site")
for (poly in states) polygon(poly[, 1], poly[, 2], col = "#F4F1E8",
                             border = "#9A9A9A", lwd = 0.55)
for (cl in levels(map_sites$class)) {
  d <- map_sites[map_sites$class == cl, ]
  points(d$longitude, d$latitude, pch = 21, bg = map_cols[cl],
         col = "#37474F", cex = 1.05, lwd = 0.8)
}
points(rep_row$longitude, rep_row$latitude, pch = 8,
       col = "#D95F02", cex = 2.1, lwd = 2.2)
text(rep_row$longitude, rep_row$latitude, labels = representative,
     pos = 4, offset = 0.7, cex = 0.95, font = 2, col = "#7F2704")
legend("bottomleft", legend = paste0(names(map_cols), " monthly simultaneous steps"),
       pt.bg = map_cols, pch = 21, col = "#37474F", bty = "n", cex = 0.82)
legend("bottomright", legend = paste("Representative:", representative),
       pch = 8, col = "#D95F02", pt.cex = 1.5, bty = "n", cex = 0.9)
dev.off()

# Figures 5-14: time series and value distributions for every variable.
variable_defs <- list(
  ET = list(column = "analysis_ET_mm_day", label = "Processed ET (mm day-1)"),
  soil_water_potential = list(column = "analysis_log10_psi_soil", label = "Processed log10 soil water potential"),
  VPD = list(column = "analysis_VPD", label = "Processed VPD (hPa)"),
  air_temperature = list(column = "analysis_TA", label = "Processed air temperature (deg C)"),
  net_radiation = list(column = "analysis_NETRAD", label = "Processed net radiation (W m-2)")
)
timeseries_figures <- distribution_value_figures <- character()
for (variable in names(variable_defs)) {
  def <- variable_defs[[variable]]
  ts_file <- paste0("representative_timeseries_", variable, ".png")
  dist_file <- paste0("representative_value_distribution_", variable, ".png")
  timeseries_figures[variable] <- ts_file
  distribution_value_figures[variable] <- dist_file

  png(file.path(figure_dir, ts_file), width = 1900, height = 1300, res = 170)
  par(mfrow = c(5, 1), mar = c(2.5, 6, 2, 1), oma = c(2, 0, 3, 0))
  for (scale in cfg$scales) {
    d <- read_processed(representative, scale)
    plot(d$Time, d[[def$column]], type = "l", col = scale_cols[scale], lwd = 0.55,
         xlab = "", ylab = def$label, main = display_scale(scale))
    abline(h = 0, col = "grey60", lty = 3)
  }
  mtext(paste(def$label, "at", representative), outer = TRUE, font = 2, cex = 1.15)
  mtext("Time (UTC)", side = 1, outer = TRUE)
  dev.off()

  png(file.path(figure_dir, dist_file), width = 2000, height = 650, res = 170)
  par(mfrow = c(1, 5), mar = c(5, 4, 4, 1), oma = c(0, 0, 3, 0))
  for (scale in cfg$scales) {
    d <- read_processed(representative, scale)
    x_all <- d[[def$column]]
    x_all <- x_all[is.finite(x_all)]
    if (length(x_all)) {
      limits <- unname(quantile(x_all, c(0.005, 0.995), na.rm = TRUE))
      if (!all(is.finite(limits)) || diff(limits) == 0) {
        limits <- range(x_all) + c(-0.5, 0.5)
      }
      x <- x_all[x_all >= limits[1] & x_all <= limits[2]]
      hist(x, breaks = 40, col = scale_cols[scale], border = "white",
           main = display_scale(scale), xlab = def$label,
           ylab = "Number of time steps")
      mtext(paste0("n = ", length(x_all), "; central 99% shown"),
            side = 3, line = 0.25, cex = 0.68)
    } else {
      plot.new(); title(main = display_scale(scale)); text(0.5, 0.5, "No valid values")
    }
  }
  mtext(paste0(def$label, " distribution at ", representative),
        outer = TRUE, font = 2, cex = 1.1)
  dev.off()
}

# Tables displayed in the report.
pair_summary <- cross_site[cross_site$variable_set %in% driver_labels,
                           c("scale", "variable_set", "n_sites_with_overlap",
                             "median_simultaneous", "q25_simultaneous", "q75_simultaneous",
                             threshold_cols, "median_longest_run")]
names(pair_summary) <- c("Scale", "Variable set", "Sites with overlap", "Median n",
                         "Q25 n", "Q75 n", paste0("Sites n >= ", cfg$te_readiness_thresholds),
                         "Median longest run")
pair_summary <- pair_summary[order(match(pair_summary[["Variable set"]], driver_labels),
                                   match(as.character(pair_summary$Scale), cfg$scales)), ]

rep_table <- overlap[overlap$Site_ID == representative & overlap$variable_set %in% driver_labels,
                     c("scale", "variable_set", "n_simultaneous", "overlap_start",
                       "overlap_end", "fraction_complete_within_span",
                       "longest_contiguous_run")]
names(rep_table) <- c("Scale", "Variable set", "Simultaneous n", "First simultaneous",
                      "Last simultaneous", "Fraction within span", "Longest run")
rep_table <- rep_table[order(match(rep_table[["Variable set"]], driver_labels),
                             match(as.character(rep_table$Scale), cfg$scales)), ]

method_table <- data.frame(
  Scale = c("Hourly", "6-hourly", "Daily", "Weekly", "Monthly"),
  Bin = c("Existing UTC hourly grid", "00:00, 06:00, 12:00, 18:00 UTC",
          "UTC calendar day", "Monday-starting ISO week", "UTC calendar month"),
  Minimum.hourly.coverage = c("Source hourly value", "5 of 6 hours", "18 of 24 hours",
                              "126 of 168 hours", "75% of actual month hours"),
  Analysis.series = c(
    "One-hour first difference minus centered 5-day same-hour mean difference",
    "Difference between adjacent 6-hour means minus centered 5-day same-slot mean difference",
    "Residual from 3-harmonic cyclic day-of-year model",
    "Residual from 2-harmonic cyclic ISO-week model",
    "Residual from 2-harmonic cyclic month-of-year model"),
  stringsAsFactors = FALSE
)

figure_tags <- function(named_files, prefix) {
  paste(vapply(names(named_files), function(nm) {
    paste0("<h3>", gsub("_", " ", nm),
           "</h3><img src='../04_Results/Figures/", named_files[[nm]],
           "' alt='", prefix, " ", html_escape(nm), "'>")
  }, character(1)), collapse = "")
}

failures <- log_df[log_df$status != "ok", c("Site_ID", "detail"), drop = FALSE]
html <- paste0(
"<!doctype html><html><head><meta charset='utf-8'><title>ET Synchrony 2.0 data processing</title>
<style>body{font-family:Segoe UI,Arial,sans-serif;max-width:1180px;margin:40px auto;padding:0 24px;color:#263238;line-height:1.58}h1,h2{color:#24465f}h3{color:#3b657e;margin-top:30px}img{max-width:100%;height:auto;border:1px solid #d7d7d7;margin:10px 0 24px;background:white}table{border-collapse:collapse;width:100%;font-size:13px;margin:12px 0 28px}th,td{border:1px solid #d6dde2;padding:7px;text-align:right}th{background:#eef3f6}td:first-child,th:first-child{text-align:left}.note{background:#eef7f5;border-left:4px solid #66C2A5;padding:12px 16px;margin:14px 0}.warn{background:#fff6e6;border-left:4px solid #FC8D62;padding:12px 16px;margin:14px 0}code{background:#eef1f3;padding:2px 5px}.formula{font-family:Cambria,serif;background:#f7f8f9;padding:10px 14px}.small{font-size:13px;color:#52616b}</style></head><body>",
"<h1>ET Synchrony 2.0: multiscale data-processing and TE-readiness report</h1>",
"<p><b>Generated:</b> ", format(Sys.time(), tz = "UTC", usetz = TRUE), "</p>",
"<p>This report documents how AmeriFlux ET and hydroclimatic-driver series were transformed at five temporal scales and quantifies exact simultaneous observations for later transfer-entropy analysis. No transfer entropy is calculated here.</p>",
"<div class='note'><b>Run summary.</b> ", sum(log_df$status == "ok"), " sites completed and ", sum(log_df$status != "ok"), " failed. Site-level processed files: ", length(cfg$scales) * sum(log_df$status == "ok"), ". The example site is <b>", representative, "</b> at ", round(representative_reason$latitude, 3), " degrees N, ", round(abs(representative_reason$longitude), 3), " degrees W.</div>",
"<h2>1. Input variables and physical calculations</h2>",
"<p>The workflow starts from the original project's QC-filtered hourly site files. That earlier stage averaged available AmeriFlux channels for a variable after its fine-resolution QC filter and aggregated half-hourly records to UTC hours. No value is interpolated here. A complete UTC hourly grid is inserted so that a difference is never calculated across a missing timestamp.</p>",
"<p><b>ET.</b> Latent heat flux is converted to an equivalent daily ET rate:</p><div class='formula'>ET (mm day<sup>-1</sup>) = LE (W m<sup>-2</sup>) x 86,400 / (2.45 x 10<sup>6</sup>).</div>",
"<p>Negative calculated ET is set to zero, following the reference project. Aggregated files retain the mean ET rate and report an interval total equal to mean rate multiplied by interval duration.</p>",
"<p><b>Soil water potential.</b> Soil-water-content percentage is converted to volumetric fraction and combined with site porosity, air-entry potential, and pore-size parameter:</p><div class='formula'>psi = psi<sub>ae</sub> [(SWC/100) / porosity]<sup>-1/lambda</sup>.</div>",
"<p>The saved value is positive suction magnitude in kPa (numerically J kg<sup>-1</sup>). The future TE input is <code>log10_psi_soil</code>, because the untransformed magnitude is strongly right-skewed. VPD, air temperature, and net radiation retain hPa, degrees C, and W m<sup>-2</sup>.</p>",
"<h2>2. Exact transformation at each temporal scale</h2>", table_html(method_table),
"<p><b>Hourly.</b> For level series x at hour t, delta x<sub>t</sub> = x<sub>t</sub> - x<sub>t-1</sub>. The analysis value subtracts the mean delta x at that clock hour over offsets -2, -1, 0, +1, and +2 days. Available values are used in the local mean; output remains NA when the current difference is NA. This removes the locally repeating diurnal change while retaining departures from it.</p>",
"<p><b>6-hourly.</b> Hourly levels are averaged within four fixed UTC bins per day. A variable is valid only with at least five of six hourly values. Differences are between adjacent 6-hour bins. The same-slot adjustment compares, for example, a 06:00 bin only with 06:00 bins over the centered five-day window.</p>",
"<p><b>Daily.</b> Hourly levels are averaged when at least 18 of 24 values exist. For each variable separately, the seasonal expectation is fitted as an intercept plus three sine-cosine pairs using day of year and a 365.25-day period. Analysis value equals observed daily mean minus fitted seasonal expectation.</p>",
"<p><b>Weekly.</b> Monday-starting weeks require at least 126 of 168 hourly values. The anomaly is the residual from two sine-cosine pairs using ISO week and a 52.1775-week period.</p>",
"<p><b>Monthly.</b> Calendar months require 75% of their actual possible hours. The anomaly is the residual from two sine-cosine pairs using month and a 12-month period. At least 18 valid monthly values are required to fit an anomaly model.</p>",
"<div class='warn'><b>Important cross-scale point.</b> Subdaily products are departures in <i>change</i>; daily and coarser products are departures in <i>level</i> from a seasonal cycle. TE magnitudes should be compared as information measures, not as though processed amplitudes had identical physical meaning across scales.</div>",
"<h2>3. Simultaneous temporal coverage for TE</h2>",
"<p>A site-scale pair is counted only where ET and the named driver are both finite on the same processed timestamp. The site-level table records simultaneous n, first and last simultaneous dates, completeness within that span, number of valid segments, and longest uninterrupted run. Merely having overlapping calendar ranges is therefore not treated as sufficient.</p>",
"<p>Thresholds of 30, 100, 300, and 1,000 observations are descriptive—not automatic declarations that TE is valid. The eventual minimum depends on lag length, estimator settings, shuffle procedure, and sensitivity tests.</p>",
"<p><a href='../04_Results/Tables/te_input_overlap_by_site.csv'>Download exact site-level overlap table</a> | <a href='../04_Results/Tables/te_input_overlap_cross_site_summary.csv'>Download cross-site summary</a></p>",
"<img src='../04_Results/Figures/sites_meeting_sample_thresholds.png' alt='Sites meeting thresholds'>", table_html(pair_summary),
"<h2>4. Distribution of available observations across sites</h2>",
"<p>In these histograms, each histogram observation is one AmeriFlux site. The x-axis is the log10 simultaneous processed timestamps for the stated ET-driver pair; the y-axis is the number of sites. These are sample-size distributions, not environmental-value distributions.</p>",
figure_tags(distribution_figures, "site sample distribution"),
"<h2>5. Completeness inside temporal overlap</h2>",
"<p>Each cell is the cross-site median of simultaneous observations divided by expected time steps between the first and last simultaneous observation. Dark text is used for visibility. A high value means few internal gaps; it does not imply a long record.</p><img src='../04_Results/Figures/median_pairwise_completeness.png' alt='Pairwise completeness'>",
"<h2>6. Site map and representative-site selection</h2>",
"<p>The map colors sites by monthly simultaneous observations for ET plus core drivers (soil water potential, VPD, and air temperature). ", representative, " is marked and labeled. Among CONUS sites having ET and all four drivers—including net radiation—at every scale, it has the smallest Euclidean distance from the cross-site median log10 sample-size profile across scales. It is a median-data-availability example, not a claim that its climate or ecosystem represents all sites.</p>",
"<img src='../04_Results/Figures/site_map_representative.png' alt='Site map'>", table_html(rep_table),
"<h2>7. Processed time series at the representative site</h2>",
"<p>Each figure shows one analyzed variable at all five scales. Zero is the no-change/no-anomaly reference. Dense hourly lines are expected because the full record is shown.</p>", figure_tags(timeseries_figures, "processed time series"),
"<h2>8. Distributions of processed environmental values</h2>",
"<p>These differ from the site sample-size distributions. Here, each observation is one processed timestamp at ", representative, ". The y-axis is time steps per value bin. To stop rare extremes from flattening subdaily panels, 40 bins and the central 99% are displayed; total finite n is printed above each panel.</p>",
figure_tags(distribution_value_figures, "processed value distribution"),
"<h2>9. Limitations and next-step checks</h2><ul><li>Soil-water-potential availability is the main restriction because many sites lack usable SWC or matched hydraulic properties.</li><li>Monthly records are not suitable merely because an anomaly was fitted; simultaneous n and longest runs must be assessed by site.</li><li>Soil water potential follows the reference workflow's average SWC-channel convention. A later sensitivity analysis should test sensor depths or individual layers where metadata permit.</li><li>Before TE, source-target vectors must be subset to shared timestamps and lag-specific sample loss recalculated at every lag.</li></ul>",
"<h2>10. Processing failures</h2>", table_html(failures),
"</body></html>")

report_path <- file.path(report_dir, "Data_processing_report.html")
writeLines(html, report_path, useBytes = TRUE)
message("Report written to ", report_path, " using representative site ", representative)
