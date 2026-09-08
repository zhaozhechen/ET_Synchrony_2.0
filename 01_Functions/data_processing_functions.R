# Reusable, dependency-free functions for multiscale AmeriFlux processing.

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

safe_min_time <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_character_ else format(min(x), tz = "UTC")
}

safe_max_time <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_character_ else format(max(x), tz = "UTC")
}

parse_utc_time <- function(x) {
  x <- sub("Z$", "", as.character(x))
  out <- as.POSIXct(x, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
  missing <- is.na(out)
  if (any(missing)) {
    out[missing] <- as.POSIXct(x[missing], tz = "UTC", format = "%Y-%m-%dT%H:%M:%S")
  }
  missing <- is.na(out)
  if (any(missing)) {
    out[missing] <- as.POSIXct(paste0(x[missing], " 00:00:00"),
                              tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
  }
  out
}

read_site_info <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(x)[names(x) == "site_id"] <- "Site_ID"
  required <- c("Site_ID", "porosity", "psi.ae", "lamda")
  absent <- setdiff(required, names(x))
  if (length(absent)) stop("Site metadata lacks: ", paste(absent, collapse = ", "))
  x
}

collapse_duplicate_hours <- function(df) {
  if (!anyDuplicated(df$Time)) return(df)
  value_names <- setdiff(names(df), c("Time", "Site_ID"))
  groups <- split(seq_len(nrow(df)), df$Time)
  out <- data.frame(
    Time = parse_utc_time(names(groups)),
    stringsAsFactors = FALSE
  )
  for (nm in value_names) {
    out[[nm]] <- vapply(groups, function(i) safe_mean(df[[nm]][i]), numeric(1))
  }
  out$Site_ID <- df$Site_ID[1]
  out[, c("Time", "Site_ID", value_names), drop = FALSE]
}

read_legacy_hourly_site <- function(path, site_info, latent_heat_j_kg = 2.45e6) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                na.strings = c("NA", "NaN", ""))
  x <- x[, !grepl("^(X|X\\.[0-9]+)$", names(x)), drop = FALSE]
  if (!all(c("Time", "Site_ID") %in% names(x))) {
    stop("Missing Time or Site_ID in ", basename(path))
  }
  x$Time <- parse_utc_time(x$Time)
  x <- x[!is.na(x$Time), , drop = FALSE]
  if (!nrow(x)) stop("No parseable timestamps in ", basename(path))

  numeric_inputs <- intersect(c("SWC", "TA", "VPD", "NETRAD", "LE_F"), names(x))
  for (nm in numeric_inputs) x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
  for (nm in setdiff(c("SWC", "TA", "VPD", "NETRAD", "LE_F"), names(x))) {
    x[[nm]] <- NA_real_
  }
  x$VPD[x$VPD < 0] <- NA_real_
  x$SWC[x$SWC <= 0] <- NA_real_

  site_id <- unique(x$Site_ID)
  site_id <- site_id[!is.na(site_id) & nzchar(site_id)]
  if (length(site_id) != 1L) stop("Expected one site in ", basename(path))
  meta <- site_info[site_info$Site_ID == site_id, , drop = FALSE]
  if (nrow(meta) != 1L) stop("No unique soil metadata for ", site_id)

  x$ET_mm_day <- x$LE_F * 86400 / latent_heat_j_kg
  x$ET_mm_day[x$ET_mm_day < 0] <- 0

  porosity <- as.numeric(meta$porosity[1])
  psi_ae <- as.numeric(meta$psi.ae[1])
  lambda <- as.numeric(meta$lamda[1])
  valid_soil <- is.finite(x$SWC) & is.finite(porosity) & porosity > 0 &
    is.finite(psi_ae) & psi_ae > 0 & is.finite(lambda) & lambda > 0
  x$psi_soil_kPa <- NA_real_
  x$psi_soil_kPa[valid_soil] <- psi_ae *
    ((x$SWC[valid_soil] / 100) / porosity)^(-1 / lambda)
  x$log10_psi_soil <- ifelse(is.finite(x$psi_soil_kPa) & x$psi_soil_kPa > 0,
                             log10(x$psi_soil_kPa), NA_real_)

  keep <- c("Time", "Site_ID", "ET_mm_day", "psi_soil_kPa",
            "log10_psi_soil", "VPD", "TA", "NETRAD", "SWC", "LE_F")
  x <- x[, keep, drop = FALSE]
  x <- x[order(x$Time), , drop = FALSE]
  x <- collapse_duplicate_hours(x)

  full_time <- seq(min(x$Time), max(x$Time), by = "hour")
  full <- data.frame(Time = full_time)
  full <- merge(full, x, by = "Time", all.x = TRUE, sort = TRUE)
  full$Site_ID <- site_id
  full
}

centered_same_slot_anomaly <- function(x, slots_per_day, window_days = 5L) {
  if (window_days %% 2L != 1L) stop("subdaily_window_days must be odd")
  half <- (window_days - 1L) %/% 2L
  shifts <- (-half):half * slots_per_day
  mat <- matrix(NA_real_, nrow = length(x), ncol = length(shifts))
  n <- length(x)
  for (j in seq_along(shifts)) {
    s <- shifts[j]
    target <- seq_len(n) + s
    ok <- target >= 1L & target <= n
    mat[ok, j] <- x[target[ok]]
  }
  local_mean <- apply(mat, 1, safe_mean)
  x - local_mean
}

add_subdaily_changes <- function(df, step_hours, window_days = 5L) {
  level_vars <- c("ET_mm_day", "log10_psi_soil", "VPD", "TA", "NETRAD")
  slots_per_day <- as.integer(24 / step_hours)
  for (nm in level_vars) {
    delta <- c(NA_real_, diff(df[[nm]]))
    delta_name <- paste0("delta_", nm)
    analysis_name <- paste0("analysis_", nm)
    df[[delta_name]] <- delta
    df[[analysis_name]] <- centered_same_slot_anomaly(delta, slots_per_day,
                                                       window_days)
  }
  df
}

floor_time_seconds <- function(x, seconds) {
  as.POSIXct(floor(as.numeric(x) / seconds) * seconds,
             origin = "1970-01-01", tz = "UTC")
}

period_start <- function(time, scale) {
  if (scale == "6_hourly") return(floor_time_seconds(time, 6 * 3600))
  d <- as.Date(time, tz = "UTC")
  if (scale == "daily") return(as.POSIXct(d, tz = "UTC"))
  if (scale == "weekly") {
    monday <- d - (as.integer(format(d, "%u")) - 1L)
    return(as.POSIXct(monday, tz = "UTC"))
  }
  if (scale == "monthly") {
    first <- as.Date(format(d, "%Y-%m-01"))
    return(as.POSIXct(first, tz = "UTC"))
  }
  stop("Unknown aggregate scale: ", scale)
}

expected_hours <- function(starts, scale) {
  if (scale == "6_hourly") return(rep(6, length(starts)))
  if (scale == "daily") return(rep(24, length(starts)))
  if (scale == "weekly") return(rep(168, length(starts)))
  d <- as.Date(starts, tz = "UTC")
  next_month <- as.Date(format(seq(d[1], by = "month", length.out = 2)[2], "%Y-%m-01"))
  # Vectorized month lengths without external packages.
  nexts <- as.Date(format(seq(min(d), max(d) + 40, by = "month"), "%Y-%m-01"))
  lookup <- setNames(as.integer(diff(nexts)), format(nexts[-length(nexts)], "%Y-%m"))
  as.numeric(lookup[format(d, "%Y-%m")]) * 24
}

aggregate_hourly <- function(hourly, scale, min_coverage = 0.75) {
  hourly$Period <- period_start(hourly$Time, scale)
  idx <- split(seq_len(nrow(hourly)), hourly$Period)
  starts <- parse_utc_time(names(idx))
  expected <- expected_hours(starts, scale)
  level_vars <- c("ET_mm_day", "psi_soil_kPa", "log10_psi_soil", "VPD", "TA", "NETRAD")
  out <- data.frame(Time = starts, Site_ID = hourly$Site_ID[1],
                    expected_hours = expected, stringsAsFactors = FALSE)
  for (nm in level_vars) {
    counts <- vapply(idx, function(i) sum(is.finite(hourly[[nm]][i])), integer(1))
    values <- vapply(idx, function(i) safe_mean(hourly[[nm]][i]), numeric(1))
    coverage <- counts / expected
    values[coverage < min_coverage] <- NA_real_
    out[[nm]] <- values
    out[[paste0("coverage_", nm)]] <- coverage
  }
  duration_days <- expected / 24
  out$ET_total_mm <- out$ET_mm_day * duration_days
  out
}

cyclic_anomaly <- function(x, phase, period, harmonics, min_observations) {
  good <- is.finite(x) & is.finite(phase)
  result <- rep(NA_real_, length(x))
  fitted <- rep(NA_real_, length(x))
  needed <- max(min_observations, 2L * harmonics + 3L)
  if (sum(good) < needed || length(unique(phase[good])) < 2L * harmonics + 2L) {
    return(list(anomaly = result, seasonal = fitted, model_ok = FALSE))
  }
  design <- data.frame(y = x[good])
  for (k in seq_len(harmonics)) {
    design[[paste0("sin", k)]] <- sin(2 * pi * k * phase[good] / period)
    design[[paste0("cos", k)]] <- cos(2 * pi * k * phase[good] / period)
  }
  fit <- lm(y ~ ., data = design)
  pred <- as.numeric(fitted(fit))
  fitted[good] <- pred
  result[good] <- x[good] - pred
  list(anomaly = result, seasonal = fitted, model_ok = TRUE)
}

add_coarse_anomalies <- function(df, scale, harmonics, min_observations) {
  date <- as.Date(df$Time, tz = "UTC")
  if (scale == "daily") {
    phase <- as.integer(format(date, "%j")); period <- 365.25
  } else if (scale == "weekly") {
    phase <- as.integer(format(date, "%V")); period <- 52.1775
  } else if (scale == "monthly") {
    phase <- as.integer(format(date, "%m")); period <- 12
  } else stop("Unsupported coarse scale")
  level_vars <- c("ET_mm_day", "log10_psi_soil", "VPD", "TA", "NETRAD")
  model_status <- logical(length(level_vars))
  names(model_status) <- level_vars
  for (nm in level_vars) {
    ans <- cyclic_anomaly(df[[nm]], phase, period, harmonics, min_observations)
    df[[paste0("seasonal_", nm)]] <- ans$seasonal
    df[[paste0("analysis_", nm)]] <- ans$anomaly
    model_status[nm] <- ans$model_ok
  }
  attr(df, "anomaly_model_status") <- model_status
  df
}

summarize_processed <- function(df, site_id, scale) {
  variables <- c(ET = "ET_mm_day", soil_water_potential = "log10_psi_soil",
                 VPD = "VPD", air_temperature = "TA", net_radiation = "NETRAD")
  rows <- lapply(names(variables), function(label) {
    level <- variables[[label]]
    analysis <- paste0("analysis_", level)
    data.frame(
      Site_ID = site_id,
      scale = scale,
      variable = label,
      start = safe_min_time(df$Time),
      end = safe_max_time(df$Time),
      n_rows = nrow(df),
      n_valid_level = sum(is.finite(df[[level]])),
      n_valid_analysis = sum(is.finite(df[[analysis]])),
      fraction_valid_level = mean(is.finite(df[[level]])),
      fraction_valid_analysis = mean(is.finite(df[[analysis]])),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

write_csv_stable <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "NA")
}
