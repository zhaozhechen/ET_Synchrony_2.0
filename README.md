# ET Synchrony 2.0

This repository prepares AmeriFlux time series for a later transfer-entropy
analysis of how synchrony between evapotranspiration (ET) and hydroclimatic
drivers changes among sites and temporal scales.

The current phase is **data processing only**. It produces aligned time series
at hourly, 6-hourly, daily, weekly, and monthly scales for:

- ET derived from latent heat flux (`LE_F`)
- soil-water-potential magnitude derived from soil water content and site soil
  hydraulic properties
- vapor-pressure deficit (VPD)
- air temperature
- net radiation, where available

## Processing choices

- The source is the QC-filtered hourly output from the original
  `ET_Synchrony` workflow. The original repository and raw data are read-only.
- Hourly and 6-hourly analysis series are first differences followed by removal
  of the local, same-time-of-day mean change over a centered 5-day window.
  Both the unadjusted differences and the adjusted analysis series are saved.
- Daily, weekly, and monthly analysis series are anomalies from a cyclic
  seasonal model using day of year, ISO week, or month of year, respectively.
- Aggregates require at least 75% valid hourly observations for each variable.
  Incomplete values are retained as `NA`, not silently filled.
- ET is stored as a mean rate (`mm day-1`) at every scale and as an estimated
  interval total (`mm`) for aggregated scales. The mean rate is the canonical
  cross-scale ET variable.
- Soil-water potential is stored as a positive suction magnitude in kPa
  (numerically equivalent to J kg-1) and as `log10` magnitude, matching the
  original project convention.

## Run

From the repository root in PowerShell:

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" "02_Analysis\01_build_multiscale_data.R"
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" "02_Analysis\02_summarize_te_readiness.R"
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" "02_Analysis\02_make_processing_report.R"
```

For a quick test on selected sites:

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" "02_Analysis\01_build_multiscale_data.R" --sites=US-Ne1,US-Ha1
```

Outputs are written beneath `04_Results/Processed_Data`, `04_Results/Tables`,
and `04_Results/Figures`. The report is
`03_Reports/Data_processing_report.html`.

## Repository map

- `00_Data/processing_config.R`: paths and scientific processing parameters
- `01_Functions/data_processing_functions.R`: reusable processing functions
- `02_Analysis/01_build_multiscale_data.R`: multiscale data builder
- `02_Analysis/02_summarize_te_readiness.R`: exact simultaneous-coverage summaries for future TE inputs
- `02_Analysis/02_make_processing_report.R`: figures and HTML report
- `03_Reports/`: generated report
- `04_Results/`: generated data, tables, logs, and figures
