# Long-record multiscale TE pilot

Results and detailed methods: [Long-record TE report](Long_record_TE_report.html).

Run in order with Rscript:

1. `02_Analysis/03_run_long_record_multiscale_TE.R`: select the maximum jointly usable hourly-record site, prepare matched inputs, run both directions at six scales, save tables and figures.
2. `02_Analysis/05_validate_long_record_TE.R`: reference checksums, core/wrapper equivalence, all-lag sample alignment and normalization tests.

`02_Analysis/04_plot_long_record_TE.R` independently redraws the two PNG/PDF comparisons and HTML report from completed tables.

The previous data-processing scripts, processed time series, figures and reports are unchanged. Today's unrelated ET/USGS plotting work was relocated to `D:/OneDrive - UW-Madison/Research/Data Center` before this pilot was added. No Git operations publish either workflow.

The TE source is a byte-identical snapshot of the generalized reference; the anomaly parameter script is retained for provenance only and must not be executed here. Its source/zero-bin flags are FALSE for anomaly data. We bypass only a wrapper NA-validation bug, not the estimator. Missing intervals remain on a regular grid.

The Inputs directory contains aligned full-grid timestamps and three analysis variables for reproducing individual lagged samples. TE tables contain full-record metrics by lag, not evolving-window estimates. Monthly estimates are computationally available but have sparse histogram support; do not treat a rising uncorrected coarse-scale peak as confirmation of the hypothesis.
