# Author: Zhaozhe Chen
# Update date: 2025.7.20

# This code is to run TE at hourly scale
# Test at Site US-Ne1

# -------- Global ------------
library(here)
# Source functions for data processing, TE implementation, and visualization
source(here("01_Data_processing","AMF_processing_functions.R"))
source(here("02_TE_implementation","TE_core_codes.R"))
source(here("05_Visualization","Plotting_functions.R"))

# Input path for hourly AMF data
# Note: this is a test dataset at US-Ne1, need to update
AMF_df <- read.csv(here("00_Data","Data_processed","AMF_hourly_US-Ne1.csv"))
# Output path for the results
Output_path <- here("02_TE_implementation","Results","Hourly_TE_test_US-Ne1")
# Test at US-Ne1
Site_ID <- "US-Ne1"

# Parameters for TE implementation
n_bin <- 11 # Number of bins for TE discritization of continuous data (e.g., SM)
max_lag <- 72 # Maximum lag to consider (This should be adjusted according to the processes and the temporal resolution of data)
Lag_Dependent_Crit <- TRUE # Determine if critical TE is lag-dependent
nshuffle <- 300 # Number of shuffles (bootstrap) for critical TE for statistical inference
alpha <- 0.05 # Confidence level for critical TE
# Set parallel session
plan(multisession,workers = availableCores()-1)
# Ensure reproducibility
set.seed(111)
# Determines if zero should be adjusted for the Sink and Source variables
ZFlagSink <- FALSE
ZFlagSource <- FALSE

# Whether to plot TS and histogram of the source and sink
TS_Hist_plot <- TRUE

# These are folding parameters to deal with extreme values (outliers) in the time series
# i.e., extreme values will be binned into the first or last bin
#lower_qt <- 0.001
#upper_qt <- 1-lower_qt

# 3 colors for var, var_diurnal_mean, and var_diurnal_anomaly
my_color <- brewer.pal(3,"Set2")

# ----------- Main -----------
# Step 1. Data processing and Calculation of secondary variables ---------------- 
# Standardize the time of the df
AMF_df <- Standardize_time(AMF_df)
# Convert psi_soil to log(psi_soil) to reduce skewness
AMF_df$log10_psi_soil <- log10(AMF_df$psi_soil)
# Calculate change in delta_log10_psi_soil
delta_log10_psi_soil <- delta_TS(AMF_df,"log10_psi_soil")
# Calculate delta_ET
delta_ET <- delta_TS(AMF_df,"ET")
# Calculate delta_VPD
delta_VPD <- delta_TS(AMF_df,"VPD")

# Put them in a df
df <- data.frame(
  Time = AMF_df$Time[2:nrow(AMF_df)],
  delta_log10_psi_soil,
  delta_VPD,
  delta_ET)

# Get diurnal anomaly
df<- Cal_diurnal_anomaly(df,"delta_log10_psi_soil",5)
df<- Cal_diurnal_anomaly(df,"delta_VPD",5)
df<- Cal_diurnal_anomaly(df,"delta_ET",5)

# Step 2. Make TS and histogram plots of input variables ----------------------
# Get TS plots and and distribution for psi, VPD and ET
# Output distributions of target variables with multiple quantile for outlier handling
if(TS_Hist_plot == TRUE){
  g_psi <- var_plots_all("delta_log10_psi_soil",bquote(Delta~psi),df,my_color,ZFlag = FALSE,nbins=n_bin)
  g_VPD <- var_plots_all("delta_VPD",bquote(Delta~VPD),df,my_color,ZFlag = FALSE,nbins=n_bin)
  g_ET <- var_plots_all("delta_ET",bquote(Delta~ET),df,my_color,ZFlag = FALSE,nbins=n_bin)
  print_g(g_psi,"TS_Hist_delta_psi",18,18)
  print_g(g_VPD,"TS_Hist_delta_VPD",18,18)
  print_g(g_ET,"TS_Hist_delta_ET",18,18)  
}

# Step 3. Run hourly TE at multiple quantiles --------------------
# Test multiple lower_qt
lower_qt_ls <- c(0.001,0.005,0.01,0.05,0.1)
# TE from hourly delta_log10_psi_soil_anomaly -> delta_ET_anomaly
psi_TE_ls <- TE_quantile(Source = df$delta_log10_psi_soil_anomaly,
                         Sink = df$delta_ET_anomaly,
                         nbins = n_bin,
                         Maxlag = max_lag,
                         alpha = alpha,
                         nshuffle = nshuffle,
                         ZFlagSource = ZFlagSource,
                         ZFlagSink = ZFlagSink,
                         Lag_Dependent_Crit = Lag_Dependent_Crit,
                         lower_qt_ls = lower_qt_ls,
                         title = "psi")

# TE from hourly delta_VPD_anomaly -> delta_ET_anomaly
VPD_TE_ls <- TE_quantile(Source = df$delta_VPD_anomaly,
                         Sink = df$delta_ET_anomaly,
                         nbins = n_bin,
                         Maxlag = max_lag,
                         alpha = alpha,
                         nshuffle = nshuffle,
                         ZFlagSource = ZFlagSource,
                         ZFlagSink = ZFlagSink,
                         Lag_Dependent_Crit = Lag_Dependent_Crit,
                         lower_qt_ls = lower_qt_ls,
                         title = "VPD")

# Step 4. Make plots of information metrics vs lag for multiple quantiles ---------------
# Loop over the tested quantiles
# Initialize a list to store the plots
g_psi_all <- list()
g_VPD_all <- list()

for(i in 1:length(lower_qt_ls)){
  # Extract parameters and results for this tested quantile
  lower_qt <- lower_qt_ls[i]
  psi_TE_df <- psi_TE_ls[[i]]
  VPD_TE_df <- VPD_TE_ls[[i]]
  # Get quantile in the title
  qt <- lower_qt * 100
  # Get title for psi and VPD
  psi_title <- bquote(Delta~psi~"\u2192"~Delta~ET~"(q ="~.(qt)*"%)")
  VPD_title <- bquote(Delta~VPD~"\u2192"~Delta~ET~"(q ="~.(qt)*"%)")
  g_psi <- lag_plots_all(psi_TE_df,psi_title)
  g_VPD <- lag_plots_all(VPD_TE_df,VPD_title)
  g_psi_all[[i]] <- g_psi
  g_VPD_all[[i]] <- g_VPD
}

# Output these plots
g_psi <- plot_grid(plotlist = g_psi_all,ncol=1)
print_g(g_psi,"TE_lag_psi",18,18)
g_VPD <- plot_grid(plotlist = g_VPD_all,ncol=1)
print_g(g_VPD,"TE_lag_VPD",18,18)

# Step 5. Test lag from the opposite direction --------------------------------
# Test the other way, only do qt at 0.5%
lower_qt_ls <- 0.005
# TE from hourly delta_ET_anomaly -> delta_log10_psi_soil_anomaly -----
psi_TE_ls <- TE_quantile(Source = df$delta_ET_anomaly,
                         Sink = df$delta_log10_psi_soil_anomaly,
                         nbins = n_bin,
                         Maxlag = max_lag,
                         alpha = alpha,
                         nshuffle = nshuffle,
                         ZFlagSource = ZFlagSource,
                         ZFlagSink = ZFlagSink,
                         Lag_Dependent_Crit = Lag_Dependent_Crit,
                         lower_qt_ls = lower_qt_ls,
                         title = "ET_to_psi")
# Make lag plots
# Get quantile in the title
qt <- lower_qt_ls * 100
# Get title for psi and VPD
psi_title <- bquote(Delta~ET~"\u2192"~Delta~psi~"(q ="~.(qt)*"%)")
g_ET_psi <- lag_plots_all(psi_TE_ls[[1]],psi_title)
print_g(g_ET_psi,"TE_lag_ET_to_psi",18,3.6)

# TE from hourly delta_ET_anomaly -> delta_VPD ----
VPD_TE_ls <- TE_quantile(Source = df$delta_ET_anomaly,
                         Sink = df$delta_VPD_anomaly,
                         nbins = n_bin,
                         Maxlag = max_lag,
                         alpha = alpha,
                         nshuffle = nshuffle,
                         ZFlagSource = ZFlagSource,
                         ZFlagSink = ZFlagSink,
                         Lag_Dependent_Crit = Lag_Dependent_Crit,
                         lower_qt_ls = lower_qt_ls,
                         title = "ET_to_VPD")
# Make lag plots
# Get quantile in the title
qt <- lower_qt_ls * 100
# Get title for psi and VPD
VPD_title <- bquote(Delta~ET~"\u2192"~Delta~VPD~"(q ="~.(qt)*"%)")
g_ET_VPD <- lag_plots_all(VPD_TE_ls[[1]],VPD_title)
print_g(g_ET_VPD,"TE_lag_ET_to_VPD",18,3.6)

# Step 6. Get TE between VPD and psi ---------------------
lower_qt_ls <- 0.005
# TE from hourly delta_VPD_anomaly -> delta_psi_anomaly
VPD_psi_TE_ls <- TE_quantile(Source = df$delta_VPD_anomaly,
                             Sink = df$delta_log10_psi_soil_anomaly,
                             nbins = n_bin,
                             Maxlag = max_lag,
                             alpha = alpha,
                             nshuffle = nshuffle,
                             ZFlagSource = ZFlagSource,
                             ZFlagSink = ZFlagSink,
                             Lag_Dependent_Crit = Lag_Dependent_Crit,
                             lower_qt_ls = lower_qt_ls,
                             title = "VPD_to_psi")
# Make lag plots
# Get quantile in the title
qt <- lower_qt_ls * 100
# Get title for psi and VPD
VPD_psi_title <- bquote(Delta~VPD~"\u2192"~Delta~psi~"(q ="~.(qt)*"%)")
g_VPD_psi <- lag_plots_all(VPD_psi_TE_ls[[1]],VPD_psi_title)
print_g(g_VPD_psi,"TE_lag_VPD_to_psi",18,3.6)

# TE from hourly delta_psi_anomaly -> delta_VPD_anomaly
psi_VPD_TE_ls <- TE_quantile(Source = df$delta_log10_psi_soil_anomaly,
                             Sink = df$delta_VPD_anomaly,
                             nbins = n_bin,
                             Maxlag = max_lag,
                             alpha = alpha,
                             nshuffle = nshuffle,
                             ZFlagSource = ZFlagSource,
                             ZFlagSink = ZFlagSink,
                             Lag_Dependent_Crit = Lag_Dependent_Crit,
                             lower_qt_ls = lower_qt_ls,
                             title = "psi_to_VPD")
# Make lag plots
# Get quantile in the title
qt <- lower_qt_ls * 100
# Get title for psi and VPD
psi_VPD_title <- bquote(Delta~psi~"\u2192"~Delta~VPD~"(q ="~.(qt)*"%)")
g_psi_VPD <- lag_plots_all(psi_VPD_TE_ls[[1]],psi_VPD_title)
print_g(g_psi_VPD,"TE_lag_psi_to_VPD",18,3.6)

# Step 7. Get peak TE between each two processes --------------------------
# Only keep qt = 0.5%

psi_ET <- readRDS(paste0(Output_path,"/TE_df_ls_psi.rds"))[[2]]
print(peak_lag(psi_ET))
VPD_ET <- readRDS(paste0(Output_path,"/TE_df_ls_VPD.rds"))[[2]]
print(peak_lag(VPD_ET))
ET_psi <- readRDS(paste0(Output_path,"/TE_df_ls_ET_to_psi.rds"))[[1]]
print(peak_lag(ET_psi))
ET_VPD <- readRDS(paste0(Output_path,"/TE_df_ls_ET_to_VPD.rds"))[[1]]
print(peak_lag(ET_VPD))
psi_VPD <- readRDS(paste0(Output_path,"/TE_df_ls_psi_to_VPD.rds"))[[1]]
print(peak_lag(psi_VPD))
VPD_psi <- readRDS(paste0(Output_path,"/TE_df_ls_VPD_to_psi.rds"))[[1]]
print(peak_lag(VPD_psi))




