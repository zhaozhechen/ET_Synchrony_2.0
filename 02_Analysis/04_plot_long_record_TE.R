# Can also run independently after 03_run_long_record_multiscale_TE.R.
arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])[1]
root <- normalizePath(file.path(dirname(arg),'..'),winslash='/')
out <- file.path(root,'04_Results/Long_record_TE')
p <- read.csv(file.path(out,'Tables/peak_TE_by_scale.csv'))
t <- read.csv(file.path(out,'Tables/TE_all_scales.csv'))
a <- read.csv(file.path(out,'Tables/sample_availability.csv'))
scales <- c('hourly','6_hourly','daily','weekly','biweekly','monthly')
labels <- c('Hourly','6-hourly','Daily','Weekly','Biweekly','Monthly')
colors <- c(VPD='#66C2A5',Soil='#FC8D62')
driver_labels <- c('VPD -> ET','Soil water potential -> ET')
draw_summary <- function() {
  par(mfrow=c(2,1),mar=c(4.3,6,2.8,1.2),oma=c(2,0,4,0),las=1,cex.axis=1.15,cex.lab=1.2)
  for(metric in c('peak_TEnorm','peak_lag_steps')) {
    ylim <- if(metric=='peak_lag_steps') c(-8,48) else c(0,max(pretty(p$peak_TEnorm))*1.1)
    plot(NA,xlim=c(.7,6.3),ylim=ylim,xaxt='n',yaxt='n',xlab='',ylab=if(metric=='peak_TEnorm') 'Peak normalized TE (%)' else 'Lag at normalized TE peak (steps)',
         main=if(metric=='peak_TEnorm') 'A. Relative uncertainty reduction' else 'B. Peak lag at each sampling scale')
    ticks <- if(metric=='peak_lag_steps') seq(0,40,10) else axTicks(2)
    abline(h=ticks,col='#E8E8E8'); axis(2,at=ticks); axis(1,1:6,labels)
    for(driver in names(colors)) {
      d <- p[p$driver==driver,]; d <- d[match(scales,d$scale),]
      xx <- 1:6 + if(driver=='VPD') -.035 else .035
      lines(xx,d[[metric]],col=colors[driver],lwd=3)
      points(xx,d[[metric]],col=colors[driver],pch=ifelse(d$significant,16,1),cex=1.35,lwd=2)
      if(metric=='peak_lag_steps') {
        txt <- ifelse(d$significant,paste(d$peak_lag_time,d$lag_unit),'')
        text(xx,d[[metric]],txt,pos=if(driver=='VPD')3 else 1,cex=.85,col=colors[driver],offset=.7)
      }
    }
    if(metric=='peak_TEnorm') legend('topleft',driver_labels,col=colors,lwd=3,pch=16,bty='n',cex=1.05)
  }
  mtext(paste(unique(p$Site_ID),' | Vaira Ranch, California | Matched three-variable coverage'),outer=TRUE,side=3,line=2,cex=1.3,font=2)
  mtext('11 bins; 300 shuffles; lags 0-48; one critical TE per direction and scale',outer=TRUE,side=3,line=.7,cex=1.05)
  mtext('Filled: exceeds critical TE. Open: no significant lag (peak lag omitted). Monthly estimates have limited samples.',outer=TRUE,side=1,line=.5,cex=.9)
}
draw_lags <- function() {
  par(mfrow=c(2,3),mar=c(4.2,4.8,3,1),oma=c(3,0,4,0),las=1,cex.axis=1.05,cex.lab=1.05)
  for(s in seq_along(scales)) {
    d <- t[t$scale==scales[s],]
    plot(NA,xlim=c(0,48),ylim=c(0,max(pretty(c(d$TEnorm,d$TEnormcrit)))),xlab='Lag (steps)',ylab='Normalized TE (%)',
         main=paste0(labels[s],' | N = ',a$n_joint[match(scales[s],a$scale)]))
    abline(h=axTicks(2),col='#EEEEEE')
    for(driver in names(colors)) {
      q <- d[d$driver==driver,]
      lines(q$Lag,q$TEnorm,col=colors[driver],lwd=2.4)
      lines(q$Lag,q$TEnormcrit,col=colors[driver],lty=2,lwd=1.4)
    }
  }
  mtext(paste(unique(p$Site_ID),' | Transfer entropy across all tested lags'),outer=TRUE,side=3,line=2,cex=1.4,font=2)
  mtext('Green: VPD -> ET     Orange: soil water potential -> ET     Solid: TE     Dashed: critical TE',outer=TRUE,side=3,line=.4,cex=1.1)
  mtext('Critical TE is constant in bits; its normalized value varies with sink entropy at each lag.',outer=TRUE,side=1,line=1,cex=1)
}
for(name in c('peak_TE_comparison','TE_lag_profiles')) {
  draw <- if(name=='peak_TE_comparison') draw_summary else draw_lags
  png(file.path(out,'Figures',paste0(name,'.png')),width=2400,height=1800,res=180); draw(); dev.off()
  pdf(file.path(out,'Figures',paste0(name,'.pdf')),width=13.33,height=10); draw(); dev.off()
}
html_table <- function(df) {
  for(nm in names(df)) if(is.numeric(df[[nm]]) && !is.integer(df[[nm]])) df[[nm]] <- signif(df[[nm]],4)
  paste0('<table><tr>',paste0('<th>',names(df),'</th>',collapse=''),'</tr>',
  paste(apply(df,1,function(row)paste0('<tr>',paste0('<td>',row,'</td>',collapse=''),'</tr>')),collapse=''),'</table>')
}
report <- c('<!doctype html><meta charset="utf-8"><title>Long-record multiscale TE</title>',
  '<style>body{font:17px Arial;max-width:1200px;margin:40px auto;line-height:1.5;color:#222}table{border-collapse:collapse;font-size:14px}td,th{border:1px solid #bbb;padding:7px}img{width:100%}code{background:#eee}</style>',
  '<h1>US-Var: long-record multiscale transfer entropy</h1>',
  '<p>Vaira Ranch, California (38.4133 N, 120.9508 W). Selected by the largest simultaneous hourly ET/VPD/soil-potential processed sample count among available sites (180,225), not simply the longest nominal site metadata range. Pair-count upper bounds and actual checks are in Tables/site_selection_ranking.csv. The selected record spans approximately 21 years.</p>',
  '<h2>Methods and exact algorithm provenance</h2>',
  '<p>The unchanged TE_implementation.R was copied from the user-specified Disturbance_Synchrony/01_Functions folder. Numerical parameters follow that project\'s 02_Analysis/02_TE_anomaly_main.R: 11 bins, 300 shuffles, alpha=0.05, lower/upper folding quantiles=0.001/0.999, seed=111, both zero-bin flags FALSE. The requested maximum lag is 48 steps and Lag_Dependent_Crit=FALSE. Sequential execution preserves the reference randomization behavior.</p>',
  '<p>TE = H(Yt,Yt-1) + H(Yt-1,Xt-lag) - H(Yt-1) - H(Yt,Yt-1,Xt-lag), in bits. Normalized TE = 100 TE/H(Yt). The original routine excludes the first source row and uses target rows lag+2 through N. It derives quantile bounds from nonzero data, folds tails into equal-width edge bins, and removes missing triplets only after lagging. No estimator corrections were substituted.</p>',
  '<p>One critical TE is computed from the lag-48 matrix: each column is independently shuffled with missing positions retained, 300 times. Threshold = mean(shuffled TE) + qt(0.95, df=100) times SD(shuffled TE), exactly as in the reference. Normalized thresholds can vary with lag because H(Yt) varies. Significance is TE &gt; TEcrit. These are pointwise tests, not a multiple-lag familywise test, and shuffling does not preserve serial correlation.</p>',
  '<p><strong>Wrapper compatibility:</strong> run_TE() in the reference rejects NA values because !is.finite(NA) is TRUE. To retain gaps, this workflow calls its unchanged Cal_TE_MI_main() core directly with the same seed and parameters, then reproduces the wrapper\'s normalization/significance columns. The external reference and the local snapshot are not patched.</p>',
  '<h2>Inputs across scales</h2>',
  '<p>Existing processed inputs are unchanged. ET is the nonnegative latent-heat-derived rate (LE_F times 86400/2.45e6, mm/day). Soil input is log10 positive suction magnitude, estimated from SWC and site soil hydraulic parameters, not a direct water-potential observation. Higher values indicate drier soil.</p>',
  '<ul><li>Hourly: adjacent-hour difference, then subtract the centered five-day same-hour mean of those differences.</li><li>6-hourly: mean hourly levels in fixed six-hour bins (75% coverage), adjacent-bin differences, then the centered five-day same-six-hour-slot mean difference is subtracted.</li><li>Daily: mean hourly levels (75% coverage), minus a fitted annual seasonal cycle with three harmonics using DOY/365.25.</li><li>Weekly: Monday-start seven-day means (75% coverage), minus two annual harmonics using ISO week/52.1775.</li><li>Biweekly: new nonoverlapping 14-day Monday bins anchored 1970-01-05; mean hourly levels with at least 252/336 valid hours; minus two annual harmonics using midpoint DOY/365.25. Saved separately; previous processing code is untouched.</li><li>Monthly: calendar-month mean rates (75% coverage), minus two seasonal harmonics using month/12. These are anomalies of ET rate, not monthly ET totals.</li></ul>',
  '<p>At each scale, all three analysis variables share a common finite-observation mask and time window. Leading/trailing jointly missing bins are trimmed; internal NA bins remain. Both directions therefore use identical valid triplet timestamps at every lag. Seasonal fits are those from the original full records (biweekly fits use its full record), not refitted after masking. This is a full-record lag analysis, not a rolling-window analysis.</p>',
  '<h2>Availability</h2>',html_table(a),
  '<p>Each direction has 49 rows (lags 0-48). Lag_time/Lag_unit give physical delay; target_start/target_end give the usable target-time span at each lag, and n_valid_triplets gives its actual sample count. Full regular aligned inputs are saved in Inputs. Monthly delays use calendar months, never an assumed fixed month length. At maximum lag: 48 hours, 288 hours, 48 days, 48 weeks, 96 weeks, and 48 months respectively.</p>',
  '<h2>Peak comparison</h2><img src="../04_Results/Long_record_TE/Figures/peak_TE_comparison.png">',
  '<p>Primary peak is the maximum significant normalized TE, matching the reference peak-selection convention. If no lag is significant, the maximum amplitude is shown open and its peak lag is NA. Raw-TE peaks and unfiltered maxima are also saved in peak_TE_by_scale.csv. Lag zero is included, as in the reference, and indicates contemporaneous association rather than delayed predictive transfer.</p>',
  html_table(p[,c('scale','driver','peak_TEnorm','significant','peak_lag_steps','peak_lag_time','lag_unit','n_joint','min_triplets')]),
  '<h2>Full lag profiles</h2><img src="../04_Results/Long_record_TE/Figures/TE_lag_profiles.png">',
  '<h2>Findings and interpretation limits</h2><p>VPD has the larger peak normalized TE at every scale. Both drivers exceed critical TE at hourly, six-hourly and daily scales; only VPD does at weekly scale (32-week peak). Neither driver exceeds critical TE at biweekly or monthly scales. Most fine-scale peaks are lag zero; the hourly soil peak is at two hours. These results do not establish increasing significant soil-potential control at coarse scales.</p>',
  '<p>Monthly and biweekly estimates have far fewer samples than hourly estimates. With 11 bins, the joint histogram has 1,331 possible cells, so sparse-sample bias can be substantial. Passing the computational minimum is not proof of reliable TE estimation. Changing sample size, autocorrelation and preprocessing across scales can change TE without demonstrating a change in causal importance. Centered subdaily preprocessing also uses neighboring future days; interpret results as processed-series synchrony, not real-time causal forecasts. A single site is a pilot, not a cross-site hypothesis test.</p>',
  '<p>Reproduce with Rscript 02_Analysis/03_run_long_record_multiscale_TE.R. Replot with Rscript 02_Analysis/04_plot_long_record_TE.R. Exact reference snapshots and file checksums document provenance. No files were committed or pushed.</p>')
dir.create(file.path(root,'03_Reports'),recursive=TRUE,showWarnings=FALSE)
writeLines(report,file.path(root,'03_Reports','Long_record_TE_report.html'))
