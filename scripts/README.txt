Code S1. Analysis scripts for the periodontitis–gut microbiota–Alzheimer's disease MR study.

Language: R 4.5.3. Packages: TwoSampleMR 0.7.4, ieugwasr 1.1.0 (MIT licence), MRPRESSO 1.0
(GPL-3), mr.raps 0.4.3, MRlap 0.0.3.3, MRcML, MVMR, cause 1.2.0, coloc, jsonlite, data.table, dplyr, ggplot2, patchwork.

Before running
  1. Register at https://api.opengwas.io and set the environment variable OPENGWAS_JWT
     (for example in ~/.Renviron). No access token is stored in these scripts.
  2. Download the FinnGen periodontitis summary statistics (Salminen et al., 2025) from
     https://storage.googleapis.com/fg-publication-green-public/F_2023_026_20250625/summary_statistics_periodontitis.zip
     and place periodontitis_diagnosis_results.txt.gz at
     ../FinGen Download/salminen_2025_sumstats/ relative to the project root.
  3. For scripts 22a/22b: download the harmonised Bellenguez et al. (2022) summary statistics
     (GWAS Catalog GCST90027158, file 35379992-GCST90027158-MONDO_0004975.h.tsv.gz) to
     ext_data/bellenguez/, and the European LD scores (eur_w_ld_chr, with w_hm3.snplist)
     distributed with LD score regression to ext_data/eur_w_ld_chr/.
  4. For scripts 23-29: PLINK 1.9 at ext_data/plink/plink; the 1000 Genomes EUR reference (EUR.bed/bim/fam,
     http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz) at ext_data/ld_ref/; risk-factor GWAS from the GWAS Catalog
     (GCST007474, GCST006867, GCST006900, GCST90239658, GCST90029013, GCST90014023) at ext_data/mvmr_rf/
     (file names in script 23a); GLIDE EUR_perio_excl_HCHSSOL.txt (doi:10.5523/bris.2j2rqgzedxlq02oqbb4vmycnc2)
     at ext_data/glide/.
  5. Run all scripts from the project root (folders data/, results/, scripts/).

Run order
  00_setup.R                        dataset discovery on OpenGWAS (writes data/datasets_index.rds)
  _helpers.R                        token handling, sourced by other scripts
  09_new_instruments_finngen_relaxed.R   periodontitis instruments (P<1e-5, r2<0.001, 10 Mb) -> Table S2
  10_primary_total_effect_relaxed.R      total effect on Alzheimer's disease, heterogeneity, Egger,
                                         Steiger, leave-one-out, Kunkle sensitivity -> Table 2
  11_step1_relaxed.R                     Step 1: periodontitis -> 211 MiBioGen taxa -> Table 1, Table S4
  14_taxa_independence_meff_relaxed.R    genetic-profile correlation and M_eff -> Tables S5, S6
  11_step1_relaxed.R (second pass)       applies the M_eff threshold from script 14
  12_reverse_mr_relaxed.R                reverse MR (Alzheimer's disease -> periodontitis) -> Table 3
  15_mrpresso_primary.R                  MR-PRESSO global and outlier tests
  16_fig2_volcano_ascii_pthr.R           Figure 2
  17_fig3_total_effect_composite.R       Figure 3
  18_power_corrected_neff.R              analytic power and minimum detectable effects (Methods S1)
  19_table_s7_snp_associations.R         harmonized SNP-level associations -> Table S7
  20_pleiotropy_lookup.R                 instrument lookup in OpenGWAS PheWAS and GWAS Catalog v2 -> Table S8
  21_raps_pleio_exclusion_scaling.R      MR-RAPS; re-analysis excluding pleiotropic SNPs; per-doubling ORs -> Table S9
  22a_prep_mrlap_exposure.R              formats FinnGen summary statistics for MRlap
  22b_mrlap.R                            MRlap correction for winner's curse, weak instruments and overlap -> Table S9
  23a_mvmr_data_local.R                  MVMR data: risk-factor GWAS (GWAS Catalog files), local clumping, harmonisation
  23b_mvmr_analysis.R                    MV-IVW, MV-Egger, MV-median, conditional F -> Table S13
  23c_mvmr_qhet.R                        weak-instrument-robust MVMR (Q-statistic minimisation) -> Table S13
  24_mrcml_steiger_tost.R                MR-cML-BIC/DP, Steiger filtering, equivalence tests -> Table S10
  25_coloc_gw_loci.R                     colocalisation at the three genome-wide significant loci -> Table S11
  26_lcv.R                               latent causal variable model (needs RunLCV.R/MomentFunctions.R from github.com/lukejoconnor/LCV)
  27_cause.R                             CAUSE (with loo compatibility patch) -> Table S10
  28_glide_replication.R                 replication with GLIDE periodontitis; FinnGen instruments in GLIDE -> Tables S10, S12
  29_conmix_lasso.R                      MR-ConMix and MR-Lasso -> Table S10
  30_positive_controls.R                 positive controls: smoking initiation and T2D -> periodontitis (FinnGen, GLIDE) -> Table S14

The 3-SNP genome-wide significant sensitivity instrument is the P<5e-8 subset of the same
FinnGen summary statistics, clumped with the same settings.
