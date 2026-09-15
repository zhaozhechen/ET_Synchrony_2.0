# Structural, sample-alignment, estimator-equivalence and provenance checks.
arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])[1]
root <- normalizePath(file.path(dirname(arg),'..'),winslash='/')
out <- file.path(root,'04_Results/Long_record_TE')
source(file.path(root,'01_Functions/TE_reference/TE_implementation.R'))
reference <- 'D:/OneDrive - UW-Madison/Research/Disturbance synchrony/Disturbance_Synchrony'
paths <- c(file.path(reference,'01_Functions/TE_implementation.R'),
           file.path(root,'01_Functions/TE_reference/TE_implementation.R'),
           file.path(reference,'02_Analysis/02_TE_anomaly_main.R'),
           file.path(root,'01_Functions/TE_reference/02_TE_anomaly_main_reference.R'))
hash <- tools::md5sum(paths)
stopifnot(hash[1]==hash[2],hash[3]==hash[4])
write.csv(data.frame(path=names(hash),md5=unname(hash)),file.path(out,'Tables/reference_checksums.csv'),row.names=FALSE)
set.seed(91)
d <- data.frame(time=1:700,source_variable=rnorm(700),sink_variable=rnorm(700))
params <- readRDS(file.path(out,'Tables/TE_parameters.rds'))
r <- do.call(run_TE,c(list(Site_df=d,Output_path=out),params))
set.seed(params$seed)
z <- Cal_TE_MI_main(d$source_variable,d$sink_variable,nbins=11,nshuffle=300,
                   alpha=.05,Maxlag=48,ZFlagSink=FALSE,ZFlagSource=FALSE,
                   Lag_Dependent_Crit=FALSE,lower_qt=.001,upper_qt=.999)
stopifnot(isTRUE(all.equal(r[,names(z)],z,check.attributes=FALSE)))
t <- read.csv(file.path(out,'Tables/TE_all_scales.csv'))
stopifnot(nrow(t)==588,all(is.finite(t$TE)),all(t$n_valid_triplets>0),
          all(abs(t$TEnorm-100*t$TE/t$Hy)<1e-10))
for(scale in unique(t$scale)) {
  a <- t[t$scale==scale & t$driver=='VPD', ]
  b <- t[t$scale==scale & t$driver=='Soil', ]
  stopifnot(identical(a$Lag,0:48),identical(b$Lag,0:48),
            identical(a$n_valid_triplets,b$n_valid_triplets),
            identical(a$target_start,b$target_start),identical(a$target_end,b$target_end),
            length(unique(a$TEcrit))==1,length(unique(b$TEcrit))==1)
}
writeLines(c('PASS: exact MD5 equality of both reference snapshots.',
             'PASS: unchanged core reproduces reference wrapper metrics on complete-data fixture (11 bins, 300 shuffles, 49 lags).',
             'PASS: 12 analyses x 49 lags = 588 finite TE rows.',
             'PASS: both drivers use identical valid triplet counts and target timestamp spans at every lag.',
             'PASS: exactly one critical TE in bits per analysis.',
             'PASS: normalized TE exactly matches 100*TE/Hy.'),file.path(out,'Tables/validation.txt'))
message('All validation checks passed.')
