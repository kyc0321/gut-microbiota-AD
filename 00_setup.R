## 00_setup.R -- environment and OpenGWAS authentication
##
## This pipeline reproduces the taxon-wide Mendelian randomization screen
## reported in "A Correction-Aware Mendelian Randomization Screen of 211 Gut
## Microbial Taxa for Alzheimer's Disease Risk" (Human Genetics and Genomics
## Advances).
##
## Requires an OpenGWAS API token (https://api.opengwas.io) set as the
## environment variable OPENGWAS_JWT -- e.g. in your ~/.Renviron file:
##   OPENGWAS_JWT=your_token_here
## Do NOT hard-code the token in any script committed to this repository.

suppressPackageStartupMessages({
  library(TwoSampleMR)   # v0.7.4
  library(ieugwasr)      # v1.1.0
  library(MRPRESSO)      # v1.0
  library(dplyr)
  library(data.table)
})

stopifnot("Set OPENGWAS_JWT before running this pipeline (see header of this file)." =
            nzchar(Sys.getenv("OPENGWAS_JWT")))

for (d in c("data", "results", "figures")) dir.create(d, showWarnings = FALSE)

# --- Data sources (OpenGWAS accession IDs) ---
MIBIOGEN_ID_FIRST <- "ebi-a-GCST90016908"
MIBIOGEN_ID_LAST  <- "ebi-a-GCST90017118"   # 211 taxa, MiBioGen Consortium (Kurilshikov et al. 2021)
AD_OUTCOME_ID     <- "ebi-a-GCST90027158"   # Bellenguez et al. 2022 EADB GWAS, N = 487,511
AD_REVERSE_INSTR_P <- 5e-8                  # genome-wide-significant threshold for reverse-MR instruments

cat("[OK] 00_setup.R complete.\n")
