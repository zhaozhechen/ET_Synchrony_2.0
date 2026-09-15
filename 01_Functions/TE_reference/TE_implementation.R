# Author: Zhaozhe Chen
# Update Date: 2026.8.5

# Generalized codes for TE implementation

# This function is to deal with outliers before discretization of continuous data
# Assuming the upper and lower boundaries are always provided for simplification
# Input include:
# The time series need to be discritized: var
# The total number of bins: nbins
# The lower boundary for folding the first bin: lower_bd
# The upper boundary for folding the last bin: upper_bd
# Output1: the number of points in each bin after accounting for the outliers: counts
# Output2: the edges for the bins: breaks
histogram <- function(var,nbins,lower_bd,upper_bd){
  # Get # of outliers on the two sides
  lower_count <- sum(var <= lower_bd,na.rm=TRUE)
  upper_count <- sum(var >= upper_bd,na.rm=TRUE)
  # Filtered data
  filtered <- var[var > lower_bd & var < upper_bd]
  # Note: in python, np.histogram automatically find the min and max of the data, then divide the range into equal-width bins
  # To get the same behavior, manually get break points for bins
  # Get the breaking points for bins
  breaks <- seq(from = min(filtered),to = max(filtered),length.out = nbins + 1)
  # Discretize the data, left-closed bins
  h <- hist(filtered,breaks = breaks,plot=FALSE,right = FALSE)
  # Add outlier counts to the first and last bins
  h$counts[1] <- h$counts[1] + lower_count
  h$counts[length(h$counts)] <- h$counts[length(h$counts)] + upper_count
  return(list(counts = h$counts,
              breaks = breaks))
}

# This function is to discretize continuous data into bins based on bin edges
# The outliers are put into the first/last bins
# Input include:
# The time series need to be discritized: var
# The edges of the bins: binEdges
# The lower boundary for folding the first bin: lower_bd
# The upper boundary for folding the last bin: upper_bd
# Output: returns the bins each point belongs to
digitize <- function(var,binEdges,lower_bd,upper_bd){
  # Get the number of bins, this is necessary when adjusting for zero
  nbins <- length(binEdges) - 1
  # Discretize the data, left-closed
  bin_id <- cut(var,breaks = binEdges,include.lowest = TRUE,right = FALSE,labels = FALSE)
  # Assign outliers to the first and last bins
  bin_id[var <= lower_bd] <- 1
  bin_id[var >= upper_bd] <- nbins
  return(bin_id)
}

# This function calculates the upper and lower boundary based on quantile
# Input include:
# The time series need to be discritized: var
# The lower boundary for folding the first bin: lower_qt as in quantile
# The upper boundary for folding the last bin: upper_qt as in quantile
# Output: returns the upper and lower boundary as in their original scale
find_bounds <- function(var,lower_qt,upper_qt){
  lower_bd <- quantile(var,probs = lower_qt,na.rm = TRUE)
  upper_bd <- quantile(var,probs = upper_qt,na.rm = TRUE)
  return(list(lower_bd,upper_bd))
}

# This function adjusts zero values in the data when discretizing continuous data, true zero values are put in the first bin
# While the non-zero values are put into the rest n-1 bins, and shift 1 bin to the right
# Input are:
# TS data to be processed (var)
# number of total bins for descretization: nbins
# a tiny value for edge extension: ths (default 1e-6)
# lower and upper bd for the outliers
# Output is: returns the bins each point belongs to, after adjusting for zero values
ZeroAdjustment <- function(var,nbins,ths = 1e-6,lower_bd,upper_bd){
  # Initialize a vector to store output bin index
  bin_idx <- rep(NA,length(var))
  # Get index for zero and non-zero values
  zero_idx <- which(var == 0)
  nonzero_idx <- which(var !=0 & !is.na(var))
  # Get non-zero values
  nonzero_values <- var[nonzero_idx]
  # Get histogram info for nonzero values, accounting for outliers
  h <- histogram(nonzero_values,nbins = nbins - 1,lower_bd,upper_bd)
  bin_edges <- h$breaks
  # Expand both edges a little
  bin_edges[1] <- bin_edges[1] - ths
  bin_edges[length(bin_edges)] <- bin_edges[length(bin_edges)] + ths
  # Discretize nonzero values into bins
  bin_id_nonzero <- digitize(nonzero_values,bin_edges,lower_bd,upper_bd)
  # Shift all these nonzero values 1 bin to the right
  bin_idx[nonzero_idx] <- bin_id_nonzero + 1
  # All 0 are put in the first bin
  bin_idx[zero_idx] <- 1
  return(bin_idx)
}

# This function calculates joint bin counts for 3D matrix
# Input includes:
# a 3-D matrix: M (Xlagged,Yt,Yt-1)
# number of total bins for descretization: nbins
# lower and upper bd for the outliers (here they are vectors for the three columns)
# ZFlag: A logical vector of length 3, indicating which column needs zero-adjustment
# Output1: A vector of joint bin counts: N (length nbins^3)
# Output2: corr: correlation between the first two columns (Xlagged and Yt)
joint_entropy3D <- function(M,nbins,lower_bd,upper_bd,ZFlag){
  # Only keep rows that have no NA
  M <- M[complete.cases(M),]
  # Note: use 10e-4 here to be consistent with the python version
  ths <- 10e-4
  # Discretize each column with or without zero-adjustment
  bin_list <- vector("list",3)
  # Loop over each of the columns
  for(i in 1:3){
    # If zero-adjustment is needed
    if(ZFlag[i]){
      # Get bin index
      bin_list[[i]] <- ZeroAdjustment(M[,i],nbins,ths = ths,lower_bd[i],upper_bd[i])
    }else{
      # If no zero-adjustment is needed
      h <- histogram(M[,i],nbins,lower_bd[i],upper_bd[i])
      binEdges <- h$breaks
      # Expand the two edges
      binEdges[1] <- binEdges[1] - ths
      binEdges[length(binEdges)] <- binEdges[length(binEdges)] + ths
      # Get bin index
      bin_list[[i]] <- digitize(M[,i],binEdges,lower_bd[i],upper_bd[i])
    }
  }
  # Convert 3D bin_idx to 1D bin_idx
  joint_bin <- (bin_list[[1]]-1)*nbins^2 + (bin_list[[2]]-1)*nbins + bin_list[[3]]
  # Count number of obs in each bin (joint counts)
  N <- tabulate(joint_bin,nbins^3)
  # Also calculates correlation between the first and second column (Xlagged and Yt)
  corr <- cor(M[,1],M[,2])
  return(list(N=N,corr = corr))
}

# This function calculates Shannon entropy from counts
# This count could be from any dimensions, but needs to be a vector
# Input: count
# Output: Shannon entropy (Unit:bits)
cal_entropy <- function(counts){
  # Only keep positive counts
  counts <- counts[counts>0]
  # Convert to probs
  probs <- counts/sum(counts)
  # Calculate Shannon entropy
  H <- -sum(probs * log2(probs))
  return(H)
}

# This function calculates 1D entropy of X and Y, mutual information (MI(X,Y)), and TE(X->Y)
# Input include:
# a vector of joint bin counts, which is the result from joint_entropy3D
# number of total bins for descretization: nbins
# Output includes:
# Hx, Hy, MI(X,Y), and TE(X->Y)
# MI(X,Y) = H(X) + H(Y) - H(X,Y)
# TE(x->Y) = H(Yt,Yt-1) + H(Yt-1,Xt-lag) - H(Yt-1) - H(Yt,Yt-1,Xt-lag)
cal_info_metrics_3D <- function(N,nbins){
  # If N is not full length, add 0 to the end
  if(length(N) < nbins^3){
    N <- c(N,rep(0,nbins^3 - length(N)))
  }
  # Convert N to 3D array: Xlagged, Yt, and Yt-1
  N3 <- array(N,dim = c(nbins,nbins,nbins))
  
  # 1D marginal counts
  Mxt <- apply(N3, 3, sum)
  Myt <- apply(N3, 2, sum)
  Myt_1 <- apply(N3, 1, sum) 
  
  # 2D marginal counts
  Mxtyt <- apply(N3, c(2,3), sum)
  Mytyt_1 <- apply(N3, c(1,2), sum) 
  Myt_1xt <- apply(N3, c(1,3), sum)
  
  # Calculate entropies
  Hxt <- cal_entropy(Mxt)
  Hyt <- cal_entropy(Myt)
  Hyt_1 <- cal_entropy(Myt_1)
  Hxtyt <- cal_entropy(Mxtyt)
  Hytyt_1 <- cal_entropy(Mytyt_1)
  Hyt_1xt <- cal_entropy(Myt_1xt)
  Hxtytyt_1 <- cal_entropy(N)
  
  # Calculate MI and TE
  MI <- Hxt + Hyt - Hxtyt
  TE <- Hytyt_1 + Hyt_1xt - Hyt_1 - Hxtytyt_1
  return(list(Hxt = Hxt,Hyt = Hyt, MI = MI, TE = TE))
}

# This function is to shift TS input and generate a shifted matrix
# Input include:
# The original TS of the source: X
# The original TS of the sink: Y
# Lag for X: lag
# The output is a matrix of:[Xlag,Yt,Yt-1]
Lag_Data <- function(X,Y,lag){
  n <- length(X)
  # Shift the TS
  x_lag <- X[2:(n-lag)]
  yt <- Y[(lag+2):n]
  yt_1 <- Y[(lag+1):(n-1)]
  # Combine them
  M <- cbind(x_lag,yt,yt_1)
  return(M)
}

# This function is to shuffle matrix, for calculation of critical TE values
# It shuffles each column in the matrix, while preserving the locations of NA
# Input is the matrix to shuffle: M
# Output: the shuffled matrix
shuffle_matrix <- function(M){
  # Initialize a blank matrix
  M0 <- matrix(NA,nrow=nrow(M),ncol=ncol(M))
  # Loop over the columns
  for(i in 1:ncol(M)){
    col_values <- M[,i]
    non_na_idx <- which(!is.na(col_values))
    # Only shuffle non-NA values
    shuffled_values <- sample(col_values[non_na_idx])
    # Put them pack
    M0[non_na_idx,i] <- shuffled_values
  }
  return(M0)
}

# This function shuffles each column of the input matrix
# Then apply joint_entropy_3D to the shuffled matrix
# Input includes:
# a 3-D matrix: M (Xlagged,Yt,Yt-1)
# number of total bins for descretization: nbins
# lower and upper bd for the outliers (here they are vectors for the three columns)
# ZFlag: A logical vector of length 3, indicating which column needs zero-adjustment
# Output1: A vector of joint bin counts: N (length nbins^3)
# Output2: corr: correlation between the first two columns (Xlagged and Yt)
joint_3D_shuffle <- function(M,nbins,lower_bd,upper_bd,ZFlag){
  # Shuffle the matrix
  Ms <- shuffle_matrix(M)
  # Get joint counts
  results <- joint_entropy3D(Ms,nbins,lower_bd,upper_bd,ZFlag)
  return(results)
}


# ========================
# Internal helper: the original code used future_lapply; sequential execution
# is also available so the wrapper does not require a global future plan.
.TE_lapply <- function(X, FUN, parallel) {
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package 'future.apply' is required when parallel = TRUE.")
    }
    return(future.apply::future_lapply(X, FUN, future.seed = TRUE))
  }

  return(lapply(X, FUN))
}
# ============================

# This function conducts joint_3D_shuffle for nshuffle times
# Input includes:
# a 3-D matrix: M (Xlagged,Yt,Yt-1)
# number of total bins for descretization: nbins
# lower and upper bd for the outliers (here they are vectors for the three columns)
# ZFlag: A logical vector of length 3, indicating which column needs zero-adjustment
# nshuffle: number of shuffles to be conducted (bootstrapping #)
# Output1: A list of a vector of joint bin counts: N (length nbins^3)
# Output2: A list of corr: correlation between the first two columns (Xlagged and Yt)
joint3D_critical <- function(
    M, nbins, lower_bd, upper_bd, ZFlag, nshuffle, parallel = FALSE) {
  .TE_lapply(seq_len(nshuffle), function(i) {
    joint_3D_shuffle(M, nbins, lower_bd, upper_bd, ZFlag)
  }, parallel = parallel)
}

# This function is to calculate critical values for MI,TE,and Corr
# Input includes:
# a 3-D matrix: M (Xlagged,Yt,Yt-1)
# number of total bins for descretization: nbins
# lower and upper bd for the outliers (here they are vectors for the three columns)
# ZFlag: A logical vector of length 3, indicating which column needs zero-adjustment
# nshuffle: number of shuffles to be conducted (bootstrapping #)
# alpha: alpha value for statistical inference
# Output: critical values of MI, TE, and Correlation (between Xlag and Yt)
cal_critical_TE_MI_Corr <- function(
    M, nbins, lower_bd, upper_bd, ZFlag, nshuffle, alpha,
    parallel = FALSE) {
  shuffled_results <- joint3D_critical(
    M, nbins, lower_bd, upper_bd, ZFlag, nshuffle,
    parallel = parallel
  )

  info_metric_ls <- .TE_lapply(shuffled_results, function(results) {
    metrics <- cal_info_metrics_3D(results$N, nbins)
    list(
      Hxt = metrics$Hxt,
      Hyt = metrics$Hyt,
      MI = metrics$MI,
      TE = metrics$TE,
      Corr = results$corr
    )
  }, parallel = parallel)

  extract_metric <- function(field) {
    vapply(info_metric_ls, function(x) x[[field]], numeric(1))
  }
  MI_all <- extract_metric("MI")
  TE_all <- extract_metric("TE")
  Corr_all <- extract_metric("Corr")

  # Retained from the original implementation and its Python comparison
  t_stat <- qt(1 - alpha, df = 100)
  MIcrit <- mean(MI_all) + t_stat * sd(MI_all)
  TEcrit <- mean(TE_all) + t_stat * sd(TE_all)
  Corrcrit <- mean(Corr_all) + t_stat * sd(Corr_all)

  return(list(MIcrit = MIcrit, TEcrit = TEcrit, Corrcrit = Corrcrit))
}

# This is the main function for TE and MI calculation
# Input includes:
# Source: the TS of source variable
# Sink: the TS of sink variable
# nbins: number of total bins for descretization:
# nshuffle: number of shuffle for critical value calculation
# alpha: value for statistical inference
# Maxlag: allowed maximum lag
# ZFlagSink/ZFlagSource: whether the sink/source variables need zero-adjustment (TRUE or FALSE)
# Lag_Dependet_Crit: whether need lag-dependent critical values (TRUE or FALSE)
Cal_TE_MI_main <- function(
    Source, Sink, nbins, nshuffle, alpha, Maxlag,
    ZFlagSink, ZFlagSource, Lag_Dependent_Crit,
    lower_qt, upper_qt, parallel = FALSE) {
  # Get bounds of source and sink variables (only for nonzero values)
  Source_bd <- find_bounds(Source[Source != 0], lower_qt, upper_qt)
  Sink_bd <- find_bounds(Sink[Sink != 0], lower_qt, upper_qt)
  # Set bounds for joint entropy input, follow the order of Xlag,Yt,Yt-1
  lower_bd <- c(Source_bd[[1]], Sink_bd[[1]], Sink_bd[[1]])
  upper_bd <- c(Source_bd[[2]], Sink_bd[[2]], Sink_bd[[2]])
  # Set ZFlag for joint entropy input
  ZFlag <- c(ZFlagSource, ZFlagSink, ZFlagSink)
  
  # Compute critical values once, if lag-dependent critical value is not needed
  if (!Lag_Dependent_Crit) {
    M0 <- Lag_Data(Source, Sink, Maxlag)
    global_crit <- cal_critical_TE_MI_Corr(
      M0, nbins, lower_bd, upper_bd, ZFlag, nshuffle, alpha,
      parallel = parallel
    )
  }
  
  # Loop over lags to calculate information metrics
  # Initialize a list to store all results
  results <- .TE_lapply(0:Maxlag, function(lag) {
    M <- Lag_Data(Source, Sink, lag)
    lag_result <- joint_entropy3D(M, nbins, lower_bd, upper_bd, ZFlag)
    lag_metrics <- cal_info_metrics_3D(lag_result$N, nbins)

    if (Lag_Dependent_Crit) {
      metric_crit <- cal_critical_TE_MI_Corr(
        M, nbins, lower_bd, upper_bd, ZFlag, nshuffle, alpha,
        parallel = parallel
      )
    } else {
      metric_crit <- global_crit
    }

    data.frame(
      Lag = lag,
      MI = lag_metrics$MI,
      MIcrit = metric_crit$MIcrit,
      TE = lag_metrics$TE,
      TEcrit = metric_crit$TEcrit,
      Corr = lag_result$corr,
      Corrcrit = metric_crit$Corrcrit,
      Hx = lag_metrics$Hxt,
      Hy = lag_metrics$Hyt
    )
  }, parallel = parallel)

  results_df <- do.call(rbind, results)
  rownames(results_df) <- NULL

  return(results_df)
}


# Theme for TE plots
.TE_theme <- function() {
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = "black", fill = NA),
    legend.key = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 14),
    plot.title = ggplot2::element_text(size = 14),
    axis.text = ggplot2::element_text(size = 14),
    axis.title = ggplot2::element_text(size = 14)
  )
}


# Convert common character timestamps for ordering, plotting, and lag detection
.TE_parse_time <- function(Time) {
  if (!is.character(Time)) {
    return(Time)
  }

  # Normalize ISO-8601 offsets for base R's %z parser
  Time_clean <- sub(
    "([+-][0-9]{2}):([0-9]{2})$",
    "\\1\\2",
    Time
  )
  Time_clean <- sub("Z$", "+0000", Time_clean)

  Time_POSIX <- suppressWarnings(as.POSIXct(
    Time_clean,
    tz = "UTC",
    tryFormats = c(
      "%Y-%m-%d %H:%M:%OS%z",
      "%Y-%m-%d %H:%M:%OS",
      "%Y-%m-%dT%H:%M:%OS%z",
      "%Y-%m-%dT%H:%M:%OS",
      "%m/%d/%Y %H:%M:%OS",
      "%Y-%m-%d"
    )
  ))
  if (all(!is.na(Time_POSIX))) {
    return(Time_POSIX)
  }

  Time_Date <- suppressWarnings(as.Date(Time))
  if (all(!is.na(Time_Date))) {
    return(Time_Date)
  }

  return(Time)
}


# Detect the time represented by one lag (one row) and convert it to hours
.detect_lag_axis <- function(Time, time_unit = NULL) {
  unit_alias <- c(
    second = "seconds", seconds = "seconds", sec = "seconds", secs = "seconds",
    minute = "minutes", minutes = "minutes", min = "minutes", mins = "minutes",
    hour = "hours", hours = "hours", hr = "hours", hrs = "hours",
    day = "days", days = "days",
    week = "weeks", weeks = "weeks"
  )
  seconds_per_unit <- c(
    seconds = 1,
    minutes = 60,
    hours = 60 * 60,
    days = 24 * 60 * 60,
    weeks = 7 * 24 * 60 * 60
  )

  if (!is.null(time_unit)) {
    time_unit <- tolower(time_unit)
    if (!time_unit %in% names(unit_alias)) {
      stop(
        "time_unit must be one of: seconds, minutes, hours, days, or weeks."
      )
    }
    time_unit <- unname(unit_alias[time_unit])
  }

  if (inherits(Time, "POSIXt")) {
    time_diff <- diff(as.numeric(Time))
    time_diff_seconds <- time_diff
  } else if (inherits(Time, "Date")) {
    time_diff <- diff(as.numeric(Time))
    time_diff_seconds <- time_diff * seconds_per_unit["days"]
  } else if (inherits(Time, "difftime")) {
    time_diff_seconds <- diff(as.numeric(Time, units = "secs"))
  } else if (is.numeric(Time)) {
    time_diff <- diff(Time)
    time_diff <- time_diff[is.finite(time_diff) & time_diff > 0]
    if (length(time_diff) == 0) {
      stop("A positive time step could not be detected from the time column.")
    }

    step <- stats::median(time_diff)
    irregular <- any(abs(time_diff - step) > 1e-6 * max(1, abs(step)))
    if (irregular) {
      warning(
        "The time column is not evenly spaced; lag time uses the median step."
      )
    }

    if (is.null(time_unit)) {
      return(list(
        step = step,
        unit = "time units",
        x_label = "Lag (time units)"
      ))
    }

    step_hours <- step * seconds_per_unit[time_unit] /
      seconds_per_unit["hours"]
    return(list(
      step = unname(step_hours),
      unit = "hours",
      x_label = "Lag (hours)"
    ))
  } else {
    if (!is.null(time_unit)) {
      warning(
        "time_unit was supplied, but the non-numeric time spacing cannot be calculated; ",
        "lag is shown in observation steps."
      )
    }
    return(list(
      step = 1,
      unit = "observation steps",
      x_label = "Lag (observation steps)"
    ))
  }

  time_diff_seconds <- time_diff_seconds[
    is.finite(time_diff_seconds) & time_diff_seconds > 0
  ]
  if (length(time_diff_seconds) == 0) {
    stop("A positive time step could not be detected from the time column.")
  }

  step_seconds <- stats::median(time_diff_seconds)
  irregular <- any(
    abs(time_diff_seconds - step_seconds) >
      1e-6 * max(1, abs(step_seconds))
  )
  if (irregular) {
    warning(
      "The time column is not evenly spaced; lag time uses the median step."
    )
  }

  # Date-time classes have known units. Always express lag in hours so results
  # are directly comparable across input resolutions (e.g., 15 minutes = 0.25 h).
  step <- step_seconds / seconds_per_unit["hours"]
  return(list(
    step = unname(step),
    unit = "hours",
    x_label = "Lag (hours)"
  ))
}


# This function makes a time-series plot of the target variable
TS_plot <- function(varname, df, time_name, my_title, var_color) {
  plot_df <- data.frame(
    Time = df[[time_name]],
    Value = df[[varname]]
  )

  g <- ggplot2::ggplot(
    data = plot_df,
    ggplot2::aes(x = Time, y = Value)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = Time, y = 0, yend = Value),
      color = var_color
    ) +
    .TE_theme() +
    ggplot2::labs(x = "", y = my_title) +
    ggplot2::ggtitle(my_title)

  return(g)
}


# This function makes a histogram of the target variable
plot_hist <- function(
    df, var, n_bin = 11, zero_remove = FALSE, var_color, var_label = var) {
  x <- df[[var]]
  if (zero_remove) {
    plot_df <- data.frame(Value = x[x != 0 & !is.na(x)])
    plot_n_bin <- n_bin - 1
    my_title <- "0 removed"
  } else {
    plot_df <- data.frame(Value = x[!is.na(x)])
    plot_n_bin <- n_bin
    my_title <- ""
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Value)) +
    ggplot2::geom_histogram(
      bins = plot_n_bin,
      color = "black",
      fill = var_color
    ) +
    ggplot2::labs(title = my_title, x = var_label, y = "Count") +
    .TE_theme()

  return(p)
}


# Identify the strongest lag whose metric exceeds its corresponding threshold
.TE_best_metric <- function(TE_df, varname) {
  critical_name <- paste0(varname, "crit")
  metric_values <- TE_df[[varname]]
  critical_values <- TE_df[[critical_name]]
  significant_idx <- which(
    is.finite(metric_values) &
      is.finite(critical_values) &
      metric_values > critical_values
  )

  if (length(significant_idx) > 0) {
    best_idx <- significant_idx[
      which.max(metric_values[significant_idx])
    ]
    best_lag <- TE_df$Lag_time[best_idx]
    significant <- TRUE
  } else {
    best_idx <- which.max(metric_values)
    best_lag <- NA_real_
    significant <- FALSE
  }

  list(
    best_idx = best_idx,
    peak_value = metric_values[best_idx],
    critical_value = critical_values[best_idx],
    best_lag = best_lag,
    significant = significant
  )
}


# This function plots one information metric against lag time
TE_lag_plot <- function(TE_df, Type, my_color, lag_x_label, lag_unit) {
  if (Type == "TE") {
    varname <- "TE"
    y_title <- "TE (bits)"
  } else if (Type == "MI") {
    varname <- "MI"
    y_title <- "MI (bits)"
  } else if (Type == "Corr") {
    varname <- "Corr"
    y_title <- "Correlation"
  } else if (Type == "TEnorm") {
    varname <- "TEnorm"
    y_title <- "Uncertainty reduction (%)"
  } else {
    stop("Type must be TE, TEnorm, MI, or Corr.")
  }

  varname_crit <- paste0(varname, "crit")
  plot_df <- data.frame(
    Lag_time = TE_df$Lag_time,
    Metric = TE_df[[varname]],
    Metric_crit = TE_df[[varname_crit]]
  )

  g <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = Lag_time, y = Metric)
  ) +
    ggplot2::geom_line(linewidth = 0.8, color = my_color[1]) +
    ggplot2::geom_line(
      ggplot2::aes(y = Metric_crit),
      linewidth = 0.8,
      linetype = "dashed",
      color = my_color[3]
    ) +
    .TE_theme() +
    ggplot2::labs(x = lag_x_label, y = y_title) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.12))
    )

  if (Type %in% c("TE", "MI", "TEnorm", "Corr")) {
    best_metric <- .TE_best_metric(TE_df, varname)
    best_lag <- best_metric$best_lag
    metric_range <- diff(range(plot_df$Metric, na.rm = TRUE))
    label_y <- max(plot_df$Metric, na.rm = TRUE) +
      ifelse(metric_range > 0, 0.05 * metric_range, 0)
    lag_range <- diff(range(plot_df$Lag_time, na.rm = TRUE))

    if (best_metric$significant) {
      best_label <- format(best_lag, digits = 3, trim = TRUE)
      best_lag_unit <- if (isTRUE(all.equal(best_lag, 1)) &&
                           lag_unit == "hours") "hour" else lag_unit
      if (best_lag <= min(plot_df$Lag_time) + 0.25 * lag_range) {
        label_x <- best_lag + 0.04 * lag_range
        label_hjust <- 0
      } else if (best_lag >= min(plot_df$Lag_time) + 0.75 * lag_range) {
        label_x <- best_lag - 0.04 * lag_range
        label_hjust <- 1
      } else {
        label_x <- best_lag
        label_hjust <- 0.5
      }
      g <- g +
        ggplot2::geom_vline(
          xintercept = best_lag,
          color = my_color[2],
          linewidth = 0.8
        ) +
        ggplot2::annotate(
          "text",
          x = label_x,
          y = label_y,
          label = paste0("Best lag = ", best_label, " ", best_lag_unit),
          size = 5,
          hjust = label_hjust,
          color = my_color[2]
        )
    } else {
      g <- g +
        ggplot2::annotate(
          "text",
          x = min(plot_df$Lag_time, na.rm = TRUE) + 0.04 * lag_range,
          y = label_y,
          label = "Best lag = NA",
          size = 5,
          hjust = 0,
          color = my_color[2]
        )
    }
  }

  return(g)
}


# This function combines the four lag plots from the previous implementation
lag_plots_all <- function(
    TE_df, my_title, my_color, lag_x_label, lag_unit) {
  g_TE <- TE_lag_plot(TE_df, "TE", my_color, lag_x_label, lag_unit) +
    ggplot2::ggtitle(my_title) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
  g_TEnorm <- TE_lag_plot(
    TE_df, "TEnorm", my_color, lag_x_label, lag_unit
  )
  g_MI <- TE_lag_plot(TE_df, "MI", my_color, lag_x_label, lag_unit)
  g_Corr <- TE_lag_plot(TE_df, "Corr", my_color, lag_x_label, lag_unit)

  g_info <- cowplot::plot_grid(
    g_TE, g_TEnorm, g_MI, g_Corr,
    nrow = 1,
    align = "hv"
  )

  return(g_info)
}


# This function writes both PDF and PNG versions of a figure
print_g <- function(g, title, w, h, Output_path, dpi = 600) {
  dir.create(Output_path, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(Output_path, paste0(title, ".pdf"))
  png_file <- file.path(Output_path, paste0(title, ".png"))

  save_with_retry <- function(filename, ..., attempts = 5) {
    last_error <- NULL
    for (attempt in seq_len(attempts)) {
      saved <- tryCatch(
        {
          ggplot2::ggsave(filename, ...)
          TRUE
        },
        error = function(error) {
          last_error <<- error
          FALSE
        }
      )
      if (saved) {
        return(invisible(TRUE))
      }
      if (attempt < attempts) {
        Sys.sleep(1)
      }
    }
    stop(last_error)
  }

  save_with_retry(
    pdf_file, plot = g, width = w, height = h, units = "in"
  )
  save_with_retry(
    png_file, plot = g, width = w, height = h, units = "in", dpi = dpi
  )

  return(c(pdf = pdf_file, png = png_file))
}


# General wrapper for a three-column time series data frame
#
# Site_df must contain time, source_variable, and sink_variable columns.
# The input df is sorted by time, but otherwise values are not detrended, interpolated,
# aggregated, or transformed by this function.
run_TE <- function(
    Site_df,
    time_name = "time",
    source_name = "source_variable",
    sink_name = "sink_variable",
    n_bin = 11, # Number of bins for TE discritization of continuous data (e.g., SM)
    max_lag = 72, # Maximum lag to consider (This should be adjusted according to the processes and the temporal resolution of data)
    Lag_Dependent_Crit = FALSE, # Determine if critical TE is lag-dependent
    nshuffle = 300, # Number of shuffles (bootstrap) for critical TE for statistical inference
    alpha = 0.05, # Confidence level for critical TE
    ZFlagSink = TRUE, # Determines if zero should be adjusted for the Sink and Source variables
    ZFlagSource = TRUE,
    lower_qt = 0.001, # These are folding parameters to deal with extreme values (outliers) in the time series
    upper_qt = 1 - lower_qt, # i.e., extreme values will be binned into the first or last bin
    seed = 111,
    parallel = FALSE,
    Site_ID = paste0(source_name, "_to_", sink_name),
    source_varname = source_name,
    sink_varname = sink_name,
    time_unit = NULL,
    Output_path = "04_Results",
    output_figures = TRUE,
    output_data = TRUE,
    figure_width = 16,
    figure_height = 12,
    figure_dpi = 600) {
  if (!is.data.frame(Site_df)) {
    stop("Site_df must be a data frame.")
  }

  required_names <- c(time_name, source_name, sink_name)
  if (length(unique(required_names)) != 3 ||
      !all(required_names %in% names(Site_df))) {
    stop(
      "Site_df must contain three distinct columns named by time_name, ",
      "source_name, and sink_name."
    )
  }

  if (!is.numeric(Site_df[[source_name]]) ||
      !is.numeric(Site_df[[sink_name]])) {
    stop("The source and sink columns must both be numeric.")
  }
  if (any(!is.finite(Site_df[[source_name]]), na.rm = TRUE) ||
      any(!is.finite(Site_df[[sink_name]]), na.rm = TRUE)) {
    stop("The source and sink columns cannot contain Inf or -Inf.")
  }
  if (anyNA(Site_df[[time_name]])) {
    stop("The time column cannot contain missing values.")
  }
  if (anyDuplicated(Site_df[[time_name]])) {
    stop("The time column must contain unique values.")
  }
  if (length(max_lag) != 1 || max_lag < 0 || max_lag != as.integer(max_lag)) {
    stop("max_lag must be one non-negative integer.")
  }
  if (nrow(Site_df) <= max_lag + 2) {
    stop("Site_df needs more than max_lag + 2 rows.")
  }
  if (length(n_bin) != 1 || n_bin < 3 || n_bin != as.integer(n_bin)) {
    stop("n_bin must be one integer greater than or equal to 3.")
  }
  if (length(nshuffle) != 1 || nshuffle < 2 ||
      nshuffle != as.integer(nshuffle)) {
    stop("nshuffle must be one integer greater than or equal to 2.")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1.")
  }
  if (lower_qt < 0 || upper_qt > 1 || lower_qt >= upper_qt) {
    stop("lower_qt and upper_qt must satisfy 0 <= lower_qt < upper_qt <= 1.")
  }
  if (!is.character(Site_ID) || length(Site_ID) != 1 || !nzchar(Site_ID)) {
    stop("Site_ID must be one non-empty character string.")
  }
  if (!is.character(Output_path) || length(Output_path) != 1 ||
      !nzchar(Output_path)) {
    stop("Output_path must be one non-empty character string.")
  }
  if (!is.logical(output_figures) || length(output_figures) != 1 ||
      is.na(output_figures) ||
      !is.logical(output_data) || length(output_data) != 1 ||
      is.na(output_data)) {
    stop("output_figures and output_data must each be TRUE or FALSE.")
  }
  if (output_figures &&
      (!requireNamespace("ggplot2", quietly = TRUE) ||
       !requireNamespace("cowplot", quietly = TRUE))) {
    stop(
      "Packages 'ggplot2' and 'cowplot' are required when output_figures = TRUE."
    )
  }

  Site_df <- Site_df[, required_names, drop = FALSE]
  Site_df[[time_name]] <- .TE_parse_time(Site_df[[time_name]])
  if (anyDuplicated(Site_df[[time_name]])) {
    stop("The parsed time column must contain unique values.")
  }
  Site_df <- Site_df[order(Site_df[[time_name]]), , drop = FALSE]
  Source <- Site_df[[source_name]]
  Sink <- Site_df[[sink_name]]
  lag_info <- .detect_lag_axis(Site_df[[time_name]], time_unit = time_unit)

  if (sum(!is.na(Source)) <= max_lag + 2 ||
      sum(!is.na(Sink)) <= max_lag + 2) {
    stop("The source and sink need more than max_lag + 2 non-missing values.")
  }
  if (length(unique(Source[!is.na(Source) & Source != 0])) < 3 ||
      length(unique(Sink[!is.na(Sink) & Sink != 0])) < 3) {
    stop("The source and sink each need at least three distinct non-zero values.")
  }

  set.seed(seed)
  TE_df <- Cal_TE_MI_main(
    Source = Source,
    Sink = Sink,
    nbins = n_bin,
    Maxlag = max_lag,
    alpha = alpha,
    nshuffle = nshuffle,
    upper_qt = upper_qt,
    lower_qt = lower_qt,
    ZFlagSource = ZFlagSource,
    ZFlagSink = ZFlagSink,
    Lag_Dependent_Crit = Lag_Dependent_Crit,
    parallel = parallel
  )

  TE_df$TEnorm <- TE_df$TE / TE_df$Hy * 100
  TE_df$TEnormcrit <- TE_df$TEcrit / TE_df$Hy * 100
  TE_df$TE_significant <- TE_df$TE > TE_df$TEcrit
  TE_df$MI_significant <- TE_df$MI > TE_df$MIcrit
  TE_df$Lag_time <- TE_df$Lag * lag_info$step
  TE_df$Lag_unit <- lag_info$unit

  safe_Site_ID <- gsub("[^A-Za-z0-9._-]+", "_", Site_ID)
  Table_path <- file.path(Output_path, "Tables")
  Figure_path <- file.path(Output_path, "Figures")
  output_files <- character(0)

  if (output_data) {
    dir.create(Table_path, recursive = TRUE, showWarnings = FALSE)
    csv_file <- file.path(Table_path, paste0("TE_df_", safe_Site_ID, ".csv"))
    utils::write.csv(TE_df, csv_file, row.names = FALSE)
    output_files <- c(output_files, csv = csv_file)
  }

  if (output_figures) {
    # Same Set2 colors used in the previous implementation
    my_color <- c("#66C2A5", "#FC8D62", "#8DA0CB")

    g_source_TS <- TS_plot(
      source_name, Site_df, time_name, source_varname, my_color[3]
    )
    g_sink_TS <- TS_plot(
      sink_name, Site_df, time_name, sink_varname, my_color[2]
    )
    g_source_hist <- plot_hist(
      Site_df, source_name, n_bin, FALSE, my_color[3], source_varname
    )
    g_source_hist_no0 <- plot_hist(
      Site_df, source_name, n_bin, TRUE, my_color[3], source_varname
    )
    g_sink_hist <- plot_hist(
      Site_df, sink_name, n_bin, FALSE, my_color[2], sink_varname
    )
    g_sink_hist_no0 <- plot_hist(
      Site_df, sink_name, n_bin, TRUE, my_color[2], sink_varname
    )

    g_data <- cowplot::plot_grid(
      g_source_TS, g_source_hist, g_source_hist_no0,
      g_sink_TS, g_sink_hist, g_sink_hist_no0,
      nrow = 2,
      align = "hv",
      axis = "btlr",
      rel_widths = c(1, 0.5, 0.5)
    )
    g_TE <- lag_plots_all(
      TE_df,
      paste(source_varname, "->", sink_varname),
      my_color,
      lag_info$x_label,
      lag_info$unit
    )
    g_all <- cowplot::plot_grid(
      g_data, g_TE,
      nrow = 2,
      align = "hv",
      axis = "tblr",
      rel_heights = c(2, 1)
    )

    figure_files <- print_g(
      g_all,
      paste0("TE_lag_", safe_Site_ID),
      figure_width,
      figure_height,
      Figure_path,
      dpi = figure_dpi
    )
    output_files <- c(output_files, figure_files)
  }

  attr(TE_df, "direction") <- paste0(source_name, " -> ", sink_name)
  attr(TE_df, "time_name") <- time_name
  attr(TE_df, "lag_step") <- lag_info$step
  attr(TE_df, "lag_unit") <- lag_info$unit
  attr(TE_df, "output_files") <- output_files
  attr(TE_df, "parameters") <- list(
    n_bin = n_bin,
    max_lag = max_lag,
    Lag_Dependent_Crit = Lag_Dependent_Crit,
    nshuffle = nshuffle,
    alpha = alpha,
    ZFlagSink = ZFlagSink,
    ZFlagSource = ZFlagSource,
    lower_qt = lower_qt,
    upper_qt = upper_qt,
    seed = seed,
    input_time_unit = time_unit,
    lag_unit = lag_info$unit
  )

  message("TE analysis done. Outputs: ", normalizePath(Output_path))
  return(TE_df)
}
