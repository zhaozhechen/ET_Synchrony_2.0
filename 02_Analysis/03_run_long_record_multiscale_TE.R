# Long-record AmeriFlux TE, using the unchanged generalized reference estimator.
# Run from any directory with Rscript. Original processed data are read-only.
arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])[1]
root <- normalizePath(file.path(dirname(arg), '..'), winslash = '/')
source(file.path(root, '01_Functions/data_processing_functions.R'))
source(file.path(root, '01_Functions/TE_reference/TE_implementation.R'))
out <- file.path(root, '04_Results/Long_record_TE')
for (sub in c('Tables', 'Figures', 'Inputs')) dir.create(file.path(out, sub), recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) write.csv(x, file.path(out, 'Tables', name), row.names = FALSE, na = 'NA')
scales <- c('hourly', '6_hourly', 'daily', 'weekly', 'biweekly', 'monthly')
labels <- c('Hourly', '6-hourly', 'Daily', 'Weekly', 'Biweekly', 'Monthly')
drivers <- c(VPD = 'analysis_VPD', Soil = 'analysis_log10_psi_soil')
cols <- c('analysis_ET_mm_day', unname(drivers))
read_scale <- function(site, scale) read.csv(file.path(root, '04_Results/Processed_Data', scale, paste0(site, '_', scale, '.csv')))

# Pair-count upper bounds let us select the site with the largest actual
# three-variable overlap without reading all 2.6 GB of hourly files.
av <- read.csv(file.path(root, '04_Results/Tables/te_input_overlap_by_site.csv'))
a <- subset(av, scale == 'hourly' & variable_set == 'ET + VPD', c(Site_ID, n_simultaneous))
b <- subset(av, scale == 'hourly' & variable_set == 'ET + soil water potential', c(Site_ID, n_simultaneous))
rank <- merge(a, b, by = 'Site_ID', suffixes = c('_VPD', '_soil'))
rank$upper_bound <- pmin(rank$n_simultaneous_VPD, rank$n_simultaneous_soil)
rank <- rank[order(-rank$upper_bound, rank$Site_ID), ]
rank$actual_joint_hourly <- NA_integer_
best <- -1
for (i in seq_len(nrow(rank))) {
  if (rank$upper_bound[i] < best) break
  z <- read_scale(rank$Site_ID[i], 'hourly')
  n <- sum(complete.cases(z[, cols]))
  rank$actual_joint_hourly[i] <- n
  if (n > best) { best <- n; site <- rank$Site_ID[i]; hourly <- z }
}
save_table(rank, 'site_selection_ranking.csv')
message('Selected ', site, ': ', best, ' joint processed hourly observations')
hourly$Time <- parse_utc_time(hourly$Time)
stopifnot(all(diff(as.numeric(hourly$Time)) == 3600))

# Biweekly is a new scale only: fixed 14-day Monday bins anchored 1970-01-05.
# Mean available hourly levels, >=75% of 336 hours per variable; seasonal
# harmonic residuals (two annual harmonics) evaluated at bin midpoint DOY.
anchor <- as.POSIXct('1970-01-05', tz = 'UTC')
period <- as.POSIXct(as.numeric(anchor) + floor((as.numeric(hourly$Time)-as.numeric(anchor))/(14*86400))*(14*86400), origin='1970-01-01', tz='UTC')
groups <- split(seq_len(nrow(hourly)), period)
bi <- data.frame(Time = parse_utc_time(names(groups)), expected_hours = 336)
for (nm in c('ET_mm_day', 'log10_psi_soil', 'VPD')) {
  count <- vapply(groups, function(i) sum(is.finite(hourly[[nm]][i])), integer(1))
  bi[[nm]] <- vapply(groups, function(i) safe_mean(hourly[[nm]][i]), numeric(1))
  bi[[paste0('coverage_', nm)]] <- count/336
  bi[[nm]][count/336 < .75] <- NA_real_
  phase <- as.integer(format(bi$Time + 7*86400, '%j'))
  fit <- cyclic_anomaly(bi[[nm]], phase, 365.25, harmonics=2, min_observations=20)
  bi[[paste0('seasonal_', nm)]] <- fit$seasonal
  bi[[paste0('analysis_', nm)]] <- fit$anomaly
}
bi$ET_total_mm <- bi$ET_mm_day * 14
write.csv(bi, file.path(out,'Inputs',paste0(site,'_biweekly_processed.csv')), row.names=FALSE)

# Exact anomaly parameters from reference 02_TE_anomaly_main.R; max_lag is
# the requested override. Use base-R graphics instead of wrapper diagnostics.
params <- list(n_bin=11, max_lag=48, Lag_Dependent_Crit=FALSE, nshuffle=300,
               alpha=.05, ZFlagSink=FALSE, ZFlagSource=FALSE,
               lower_qt=.001, upper_qt=.999, seed=111, parallel=FALSE,
               output_figures=FALSE, output_data=FALSE)
saveRDS(params, file.path(out,'Tables','TE_parameters.rds'))
save_table(data.frame(parameter=names(params), value=vapply(params, as.character, character(1))), 'TE_parameters.csv')
all_te <- list(); peaks <- list(); counts <- list()
for (s in seq_along(scales)) {
  scale <- scales[s]
  d <- if(scale == 'hourly') hourly else if(scale == 'biweekly') bi else read_scale(site, scale)
  d$Time <- parse_utc_time(d$Time)
  stopifnot(!anyNA(d$Time), !anyDuplicated(d$Time))
  if(scale == 'monthly') {
    index <- as.integer(format(d$Time,'%Y'))*12 + as.integer(format(d$Time,'%m'))
    stopifnot(all(diff(index)==1))
  } else stopifnot(all(diff(as.numeric(d$Time)) == c(hourly=3600,'6_hourly'=21600,daily=86400,weekly=604800,biweekly=1209600)[scale]))
  good <- complete.cases(d[,cols])
  # Match both driver comparisons exactly: mask ALL three variables whenever
  # any is missing; retain internal NA rows, trim only leading/trailing gaps.
  span <- range(which(good)); d <- d[seq(span[1], span[2]), ]
  good <- complete.cases(d[,cols]); d[!good, cols] <- NA_real_
  input <- data.frame(time=d$Time, d[,cols], check.names=FALSE)
  write.csv(input, file.path(out,'Inputs',paste0(site,'_',scale,'_aligned.csv')), row.names=FALSE)
  counts[[s]] <- data.frame(Site_ID=site,scale=scale,start=as.character(min(d$Time)),end=as.character(max(d$Time)),
                            n_grid=nrow(d),n_joint=sum(good),n_missing=sum(!good))
  for (driver in names(drivers)) {
    message(format(Sys.time()), ': ', scale, ' ', driver, ' -> ET, N=',sum(good))
    x <- d[[drivers[[driver]]]]; y <- d$analysis_ET_mm_day
    # Calendar months are row-based, not fixed-duration days. Explicit
    # calendar lag labels are attached below without median-month inference.
    # Reference wrapper's !is.finite check incorrectly rejects NA. Call its
    # unchanged core directly; validate infinities without dropping NA rows.
    stopifnot(!any(is.infinite(x)), !any(is.infinite(y)),
              sum(is.finite(x))>50, sum(is.finite(y))>50)
    set.seed(params$seed)
    te <- Cal_TE_MI_main(Source=x, Sink=y, nbins=params$n_bin,
       nshuffle=params$nshuffle, alpha=params$alpha, Maxlag=params$max_lag,
       ZFlagSink=params$ZFlagSink, ZFlagSource=params$ZFlagSource,
       Lag_Dependent_Crit=params$Lag_Dependent_Crit, lower_qt=params$lower_qt,
       upper_qt=params$upper_qt, parallel=params$parallel)
    te$TEnorm <- te$TE/te$Hy*100
    te$TEnormcrit <- te$TEcrit/te$Hy*100
    te$TE_significant <- te$TE>te$TEcrit
    te$MI_significant <- te$MI>te$MIcrit
    te$Site_ID <- site; te$scale <- scale; te$driver <- driver
    te$Lag_time <- te$Lag * c(hourly=1,'6_hourly'=6,daily=1,weekly=1,biweekly=2,monthly=1)[scale]
    te$Lag_unit <- c(hourly='hours','6_hourly'='hours',daily='days',weekly='weeks',biweekly='weeks',monthly='months')[scale]
    te$record_start <- as.character(min(d$Time)); te$record_end <- as.character(max(d$Time))
    te$n_grid <- nrow(d); te$n_joint <- sum(good)
    te$n_valid_triplets <- vapply(te$Lag, function(l) sum(complete.cases(Lag_Data(x,y,l))), integer(1))
    te$target_start <- te$target_end <- NA_character_
    for(j in seq_len(nrow(te))) {
      l <- te$Lag[j]; valid <- complete.cases(Lag_Data(x,y,l))
      tt <- d$Time[(l+2):nrow(d)][valid]
      te$target_start[j] <- as.character(min(tt)); te$target_end[j] <- as.character(max(tt))
    }
    stopifnot(nrow(te)==49, identical(te$Lag,0:48), length(unique(te$TEcrit))==1,
              all(is.finite(te$TE)), all(te$n_valid_triplets>0))
    save_table(te, paste0('TE_df_',site,'_',scale,'_',driver,'_to_ET.csv'))
    all_te[[paste(scale,driver)]] <- te
    # Match reference significant-peak convention, separately for raw and
    # normalized TE. If none pass, report peak amplitude but lag is NA.
    p <- .TE_best_metric(te,'TEnorm'); r <- .TE_best_metric(te,'TE')
    peaks[[length(peaks)+1]] <- data.frame(Site_ID=site,scale=scale,driver=driver,
      peak_TEnorm=p$peak_value,critical_TEnorm=p$critical_value,significant=p$significant,
      peak_lag_steps=if(p$significant) te$Lag[p$best_idx] else NA,
      peak_lag_time=p$best_lag,lag_unit=te$Lag_unit[1],
      peak_TE_bits=r$peak_value,peak_raw_TE_lag_steps=if(r$significant) te$Lag[r$best_idx] else NA,
      unfiltered_max_TEnorm=max(te$TEnorm),unfiltered_max_lag_steps=te$Lag[which.max(te$TEnorm)],
      n_joint=sum(good),min_triplets=min(te$n_valid_triplets),max_triplets=max(te$n_valid_triplets))
  }
  save_table(do.call(rbind,counts),'sample_availability.csv')
  save_table(do.call(rbind,peaks),'peak_TE_by_scale.csv')
}
save_table(do.call(rbind,all_te),'TE_all_scales.csv')
saveRDS(all_te,file.path(out,'Tables','TE_all_scales.rds'))
capture.output(sessionInfo(),file=file.path(out,'Tables','sessionInfo.txt'))
message('All 12 analyses complete.')
source(file.path(root,'02_Analysis/04_plot_long_record_TE.R'))
