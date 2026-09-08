## run_all.R -- master pipeline
## Reproduces Table 1, Table 2, Figure 1, and Figure 2 of the manuscript.
## Requires OPENGWAS_JWT to be set (see 00_setup.R). Instrument extraction
## against the live OpenGWAS API (steps 01 and 04) takes several hours for
## the full 211-taxon panel; results are cached incrementally under data/.

steps <- c(
  "01_instruments_and_mr.R",
  "02_sensitivity_diagnostics.R",
  "03_fdr_correction.R",
  "04_reverse_mr_lead_signal.R",
  "05_mvmr_exploratory.R"
)
for (s in steps) {
  cat("\n############", s, "############\n")
  source(s)
}
cat("\nPipeline complete.\n")
