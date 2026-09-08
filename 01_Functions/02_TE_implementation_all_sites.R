# Author: Zhaozhe Chen
# Update date: 2025.8.17

# This code is to run TE at hourly scale for all AMF sites, across different time periods
# Including full TS, GS during full TS, Non-GS during full TS, full TS across years, GS across years, and Non-GS across years
# Calculate TE between variable pairs of:
# Soil water potential, ET, and T
# For both directions
# Note: The variables are pre-processed to get diurnal anomaly of change in time series (delta_TS)

# Output includes:
# One figure of all target variables for each site, to Output/Var_plots
# One figure of Lag plots for each of full TS, GS during full TS, and non-GS during full TS, to Output/Lag_plots
# One TE_df_ls for each of above, to Output/TE_df
# One figure of normalized TE vs lag for each of above across years, to Output/Lag_plots
# One TE_df_ls for each of above across years, to Output/TE_df

# -------- Global ------------
library(here)

# Server parameters ================
# Determines which site to process
#args <- commandArgs(trailingOnly = TRUE)
#arrayid <- as.integer(args[1]) + 1
#out_base <- args[2]

# The root directory for output
#Output_path <- out_base
#Output_rdata_path <- file.path(out_base,"TE_df")

# Create output folders
#dir.create(Output_path, recursive = TRUE, showWarnings = FALSE)
#dir.create(file.path(out_base,"Var_plots"), recursive = TRUE, showWarnings = FALSE)
#dir.create(file.path(out_base,"Lag_plots"), recursive = TRUE, showWarnings = FALSE)
#dir.create(Output_rdata_path, recursive = TRUE, showWarnings = FALSE)

# Data path ===================
# Source functions for data processing, TE implementation, and visualization
source(here("01_Data_processing","AMF_processing_functions.R"))
source(here("02_TE_implementation","TE_core_codes.R"))
source(here("05_Visualization","Plotting_functions.R"))
# Input path for hourly AMF data
AMF_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/02_AMF_cleaned/AMF_Hourly/AMF_sites_hourly_update/"
# Input path to AMF site info, which also includes soil info
site_info <- read.csv(here("00_data","ameriflux_site_info_update_GS.csv"))
# Output path for the results
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites/"

# Parameters for TE implementation ===============
n_bin <- 11 # Number of bins for TE discritization of continuous data (e.g., SM)
max_lag <- 72 # Maximum lag to consider (This should be adjusted according to the processes and the temporal resolution of data)
Lag_Dependent_Crit <- FALSE # Determine if critical TE is lag-dependent
nshuffle <- 300 # Number of shuffles (bootstrap) for critical TE for statistical inference
alpha <- 0.05 # Confidence level for critical TE
# Set parallel session
plan(multisession,workers = availableCores()-1)
# Ensure reproducibility
set.seed(111)
# Determines if zero should be adjusted for the Sink and Source variables
ZFlagSink <- FALSE
ZFlagSource <- FALSE
# Threshold for calculating TE, only calculate TE if number of observations greater than obs_th
obs_th <- 1000

# These are folding parameters to deal with extreme values (outliers) in the time series
# i.e., extreme values will be binned into the first or last bin
lower_qt <- 0.001
upper_qt <- 1-lower_qt

# 2 colors for growing season and non-growing season
season_color <- brewer.pal(3,"Set2")[1:2]
# 3 colors for Lag plots
my_color <- brewer.pal(3,"Set2")

# All variable pairs to consider
var_ls <- c("delta_ET_anomaly","delta_log10_psi_soil_anomaly","delta_VPD_anomaly","delta_TA_anomaly")
# Their names for simplicity
vartitle_ls <- c("ET","psi","VPD","TA")

# -------- Main ----------
# US-CGG (20), US-CS2 (24), US-Hn3 (41), US-KS1 (50), US-RGF (83)
arrayid <- 20

Site_ID <- site_info$site_id[arrayid]

# Step 1. Data processing and Calculation of secondary variables ---------------- 
# Read in hourly data for this site
AMF_df <- read.csv(paste0(AMF_path,"AMF_hourly_",Site_ID,".csv"))
# Standardize the time of the df
AMF_df <- Standardize_time(AMF_df) %>%
  select(-c(X.1,X))
# Complete time to avoid mismatch in time series alignment
AMF_df <- AMF_df %>%
  complete(Time = seq(min(Time),max(Time),by="1 hour"))

# Convert psi_soil to log(psi_soil) to reduce skewness
AMF_df$log10_psi_soil <- log10(AMF_df$psi_soil)
# Calculate change in delta_log10_psi_soil
delta_log10_psi_soil <- delta_TS(AMF_df,"log10_psi_soil")
# Calculate delta_ET
delta_ET <- delta_TS(AMF_df,"ET")
# Calculate delta_VPD
delta_VPD <- delta_TS(AMF_df,"VPD")
# Calculate delta_TA
delta_TA <- delta_TS(AMF_df,"TA")

# Put them in a df
df <- data.frame(
  Time = AMF_df$Time[2:nrow(AMF_df)],
  delta_log10_psi_soil,
  delta_VPD,
  delta_ET,
  delta_TA)

# Get diurnal anomaly
df <- Cal_diurnal_anomaly(df,"delta_log10_psi_soil",5)
df <- Cal_diurnal_anomaly(df,"delta_VPD",5)
df <- Cal_diurnal_anomaly(df,"delta_ET",5)
df <- Cal_diurnal_anomaly(df,"delta_TA",5)

# Add growing season (GS), the time between SOS and EOS
# Get SOS for this site
site_SOS <- format(as.Date(site_info$SOS[site_info$site_id == Site_ID]),"%m-%d")
site_EOS <- format(as.Date(site_info$EOS[site_info$site_id == Site_ID]),"%m-%d")

df <- df %>%
  mutate(
    # Get month and day
    md = format(Time,"%m-%d"),
    GS = if_else((site_SOS <= md & md <= site_EOS),"GS","Non-GS")
  )

# Step 2. Make TS and histogram plots of input variables ----------------------
# Get TS plots and and distribution for ET, psi, VPD, TA
# Note: all the first four variables are diurnal anomaly of delta_TS
g_ET <- var_plot_TS_Hist("delta_ET_anomaly",y_title = bquote(Delta~ET~"(mmday-1)"),
                         df,my_color = season_color,ZFlag = FALSE,nbins=n_bin)

g_psi <- var_plot_TS_Hist("delta_log10_psi_soil_anomaly",y_title = bquote(Delta~psi~"(Jkg-1 log scale)"),
                          df,my_color = season_color,ZFlag = FALSE,nbins=n_bin)

g_VPD <- var_plot_TS_Hist("delta_VPD_anomaly",y_title = bquote(Delta~VPD~"(hPa)"),
                          df,my_color = season_color,ZFlag = FALSE,nbins=n_bin)

g_TA <- var_plot_TS_Hist("delta_TA_anomaly",y_title = bquote(Delta~T~"("~degree~C~")"),
                         df,my_color = season_color,ZFlag = FALSE,nbins=n_bin)

# Combine all plots
g_var_all <- plot_grid(plotlist = c(g_ET,g_psi,g_VPD,g_TA),
                       ncol=4,align="hv",
                       axis="tblr",rel_widths = c(1.5,1,1,1))

# Add Site ID to the top of the plot
g_title <- ggdraw() + draw_label(Site_ID,fontface = "bold",size=20)
g_var_all <- plot_grid(g_title,g_var_all,
                       ncol=1,
                       rel_heights = c(0.1,2))
# Output this figure, one for each site
print_g(g_var_all,paste0("/Var_plots/Var_plots_",Site_ID),18,14)

# Step 3. Run hourly TE between each pair of variables -----------------------------
# Get all combinations of variable pairs, order matters
var_comb <- expand.grid(from = var_ls,
                        to = var_ls) %>%
  filter(from != to)

# Timing the TE calculation
start_time <- Sys.time()
# For all data in full TS =============
TE_results <- TE_all_var_pairs(var_comb,df,
                               Maxlag = max_lag,nbins = n_bin,alpha = alpha,nshuffle = nshuffle,
                               ZFlagSource = ZFlagSource,ZFlagSink = ZFlagSink,Lag_Dependent_Crit = Lag_Dependent_Crit)
TE_df_ls <- TE_results[[1]]
TE_g <- TE_results[[2]]
# Output these
saveRDS(TE_df_ls,paste0(Output_path,"TE_df/TE_df_ls_full_TS_",Site_ID,".rds"))
print_g(TE_g,paste0("/Lag_plots/Lag_plots_full_TS_",Site_ID),
        18,40)
message("Complete full TS")

# For GS in full TS =============
# Process df and assign NA to non-GS values
# Note: should not remove them, to keep temporal dependence
df_GS <- df %>%
  mutate(across(
    starts_with("delta_"),
    ~if_else(GS=="GS", ., NA)
  ))

TE_results_GS <- TE_all_var_pairs(var_comb,df_GS,
                                  Maxlag = max_lag,nbins = n_bin,alpha = alpha,nshuffle = nshuffle,
                                  ZFlagSource = ZFlagSource,ZFlagSink = ZFlagSink,Lag_Dependent_Crit = Lag_Dependent_Crit)
TE_df_ls_GS <- TE_results_GS[[1]]
TE_g_GS <- TE_results_GS[[2]]
# Output these
saveRDS(TE_df_ls_GS,paste0(Output_path,"TE_df/TE_df_ls_GS_",Site_ID,".rds"))
print_g(TE_g_GS,paste0("/Lag_plots/Lag_plots_GS_",Site_ID),
        18,40)
message("Complete full TS GS")

# For Non-GS in full TS =============
# Process df and assign NA to GS values
# Note: should not remove them, to keep temporal dependence
df_NGS <- df %>%
  mutate(across(
    starts_with("delta_"),
    ~if_else(GS=="Non-GS", ., NA)
  ))

TE_results_NGS <- TE_all_var_pairs(var_comb,df_NGS,
                                   Maxlag = max_lag,nbins = n_bin,alpha = alpha,nshuffle = nshuffle,
                                   ZFlagSource = ZFlagSource,ZFlagSink = ZFlagSink,Lag_Dependent_Crit = Lag_Dependent_Crit)
TE_df_ls_NGS <- TE_results_NGS[[1]]
TE_g_NGS <- TE_results_NGS[[2]]
# Output these
saveRDS(TE_df_ls_NGS,paste0(Output_path,"TE_df/TE_df_ls_NGS_",Site_ID,".rds"))
print_g(TE_g_NGS,paste0("/Lag_plots/Lag_plots_NGS_",Site_ID),
        18,40)
message("Complete full TS Non-GS")

# For all data across years ===========
TE_results_years <- TE_all_var_pairs_year(var_comb,df,
                                          Maxlag = max_lag,nbins = n_bin,alpha = alpha,nshuffle = nshuffle,
                                          ZFlagSource = ZFlagSource,ZFlagSink = ZFlagSink,Lag_Dependent_Crit = Lag_Dependent_Crit)
TE_df_ls_years <- TE_results_years[[1]]
TE_g_ls_years <- TE_results_years[[2]]
# Output these
saveRDS(TE_df_ls_years,paste0(Output_path,"TE_df/TE_df_ls_full_TS_years_",Site_ID,".rds"))
# Combine all plots
TE_g_years <- plot_grid(plotlist = TE_g_ls_years,ncol=6,
                        align="hv")
# Calculate height of the figure
g_h <- 4*length(TE_g_ls_years)/6
pdf(paste0(Output_path,"/Lag_plots/Lag_plots_full_TS_years_",Site_ID,".pdf"),
    width=18,height=g_h)
print(TE_g_years)
dev.off()
message("Complete full TS across years")

# For GS across years =================
TE_results_GS_years <- TE_all_var_pairs_year(var_comb,df_GS,
                                             Maxlag = max_lag,nbins = n_bin,alpha = alpha,nshuffle = nshuffle,
                                             ZFlagSource = ZFlagSource,ZFlagSink = ZFlagSink,Lag_Dependent_Crit = Lag_Dependent_Crit)
TE_df_ls_GS_years <- TE_results_GS_years[[1]]
TE_g_ls_GS_years <- TE_results_GS_years[[2]]
# Output these
saveRDS(TE_df_ls_GS_years,paste0(Output_path,"TE_df/TE_df_ls_GS_years_",Site_ID,".rds"))
# Combine all plots
TE_g_GS_years <- plot_grid(plotlist = TE_g_ls_GS_years,ncol=6,
                           align="hv")
# Calculate height of the figure
g_h <- 4*length(TE_g_ls_GS_years)/6
pdf(paste0(Output_path,"/Lag_plots/Lag_plots_GS_years_",Site_ID,".pdf"),
    width=18,height=g_h)
print(TE_g_GS_years)
dev.off()
message("Complete GS across years")

# For Non-GS across years ==============
TE_results_NGS_years <- TE_all_var_pairs_year(var_comb,df_NGS,
                                              Maxlag = max_lag,nbins = n_bin,alpha = alpha,nshuffle = nshuffle,
                                              ZFlagSource = ZFlagSource,ZFlagSink = ZFlagSink,Lag_Dependent_Crit = Lag_Dependent_Crit)
TE_df_ls_NGS_years <- TE_results_NGS_years[[1]]
TE_g_ls_NGS_years <- TE_results_NGS_years[[2]]
# Output these
saveRDS(TE_df_ls_NGS_years,paste0(Output_path,"TE_df/TE_df_ls_NGS_years_",Site_ID,".rds"))
# Combine all plots
TE_g_NGS_years <- plot_grid(plotlist = TE_g_ls_NGS_years,ncol=6,
                            align="hv")
# Calculate height of the figure
g_h <- 4*length(TE_g_ls_NGS_years)/6
pdf(paste0(Output_path,"/Lag_plots/Lag_plots_NGS_years_",Site_ID,".pdf"),
    width=18,height=g_h)
print(TE_g_NGS_years)
dev.off()
message("Complete Non-GS across years")

end_time <- Sys.time()
run_time <- as.character(end_time - start_time)
run_time <- paste(Site_ID,run_time)
message(run_time)
writeLines(run_time,paste0(Output_path,"Run_time_",Site_ID,".txt"))
message(paste("All done!!!", Site_ID))
  







