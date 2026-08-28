# Gut Microbiota–Alzheimer's Disease MR Screen

Analysis code for:

> Chu K-Y. A Correction-Aware Mendelian Randomization Screen of 211 Gut Microbial Taxa for Alzheimer's Disease Risk. *Human Genetics and Genomics Advances* (submitted).

## What this reproduces

A systematic, pre-specified two-sample Mendelian randomization (MR) screen of
all 211 gut microbial taxa profiled by the MiBioGen Consortium against
Alzheimer's disease (AD) risk, using the Bellenguez et al. (2022) EADB GWAS.
For each taxon: instrument extraction (P < 1×10⁻⁵, LD-clumped r² < 0.001 /
10,000 kb), harmonization against the AD outcome, and causal-effect
estimation via inverse-variance-weighted (IVW) regression, MR-Egger, and
weighted median. Pleiotropy (MR-Egger intercept), heterogeneity (Cochran's
Q), outlier bias (MR-PRESSO), and causal direction (Steiger filtering) are
assessed for every nominally significant (P < 0.05) taxon. Benjamini-Hochberg
FDR correction is applied across independent signal clusters, collapsing
taxonomically-nested duplicate signals (e.g. a genus that is the sole member
of its family) to a single test — see `03_fdr_correction.R`. Reverse MR
(AD as exposure) is run for the single taxon with the smallest forward P-value.

## Data sources

All data are publicly available, summary-level GWAS statistics accessed via
the [MRC IEU OpenGWAS platform](https://gwas.mrcieu.ac.uk):

- **Exposure** — MiBioGen Consortium gut microbiota GWAS (Kurilshikov et al.
  2021, *Nat Genet*), OpenGWAS accessions `ebi-a-GCST90016908` through
  `ebi-a-GCST90017118` (211 taxa).
- **Outcome** — Bellenguez et al. 2022 EADB Alzheimer's disease GWAS, OpenGWAS
  accession `ebi-a-GCST90027158` (N = 487,511).

No individual-level data are used or distributed by this code.

## Requirements

- R ≥ 4.5
- Packages: `TwoSampleMR` (0.7.4), `ieugwasr` (1.1.0), `MRPRESSO` (1.0),
  `dplyr`, `data.table`
- An OpenGWAS API token (free, https://api.opengwas.io), set as the
  environment variable `OPENGWAS_JWT` (e.g. in `~/.Renviron`). **Do not**
  hard-code the token in any script — `00_setup.R` reads it from the
  environment and will stop with an error if it is unset.

## Running

```r
source("run_all.R")
```

Runs `01_instruments_and_mr.R` → `02_sensitivity_diagnostics.R` →
`03_fdr_correction.R` → `04_reverse_mr_lead_signal.R` in sequence. Step 1 and
step 4 query the live OpenGWAS API and are the slow steps (several hours for
the full 211-taxon panel on a typical connection); both cache per-taxon
results incrementally under `data/`, so an interrupted run can simply be
re-sourced and will skip already-completed taxa.

## Outputs

- `results/ivw_by_taxon_with_fdr.csv` — per-taxon IVW estimate and FDR
  q-value; underlies manuscript Table 1 and Figure 1.
- `results/sensitivity_summary.csv` — Egger intercept, Cochran's Q, Steiger
  P, and MR-PRESSO global P for every nominally significant taxon; the lead
  signal's row underlies manuscript Table 2.
- `results/reverse_mr_lead_signal.csv` — reverse-MR estimates for the lead
  signal, underlying the last two rows of Table 2.

## License

MIT (code only; the GWAS summary statistics used are governed by their
respective original data-use agreements — see the accessions above).
