# Periodontitis, gut microbiota and Alzheimer's disease: Mendelian randomization pipeline

Analysis code and result tables for a two-sample Mendelian randomization (MR) study of
periodontitis, gut microbiota and Alzheimer's disease (AD), by Kuan-Yu Chu (manuscript in preparation).

Version 2.0.0 replaces the earlier pipeline in this repository (v1.0.0, v1.1.0, which accompanied
an earlier, differently designed version of the manuscript). Those versions remain available in the
git history and in their Zenodo records.

## What the pipeline does

1. **Exposure instruments**: periodontitis (FinnGen R11 diagnosis-based GWAS, Salminen et al. 2025),
   P < 1×10⁻⁵, LD-clumped (r² < 0.001, 10 Mb); genome-wide significant subset as sensitivity.
2. **Step 1**: periodontitis → 211 MiBioGen gut taxa; genetic-profile correlation between taxa and
   effective number of independent tests (M_eff; Li–Ji, Galwey).
3. **Total effect**: periodontitis → AD (Bellenguez et al. 2022; Kunkle et al. 2019 as sensitivity outcome);
   IVW, MR-Egger, weighted median/mode, MR-PRESSO, Steiger, leave-one-out; reverse MR.
4. **Robustness**: instrument PheWAS/GWAS Catalog lookup and pleiotropy-based exclusions; MR-RAPS;
   MRlap (winner's curse, weak instruments, sample overlap); MR-cML-BIC/DP; MR-ConMix; MR-Lasso;
   Steiger filtering; equivalence tests; multivariable MR with six AD risk factors (MV-IVW/Egger/median,
   conditional F, Q-statistic minimisation); CAUSE; colocalisation at genome-wide significant loci;
   latent causal variable model; replication with GLIDE periodontitis; positive controls
   (smoking initiation and type 2 diabetes → periodontitis).

## Repository layout

```
scripts/   R scripts, numbered in run order (see scripts/README.txt for inputs/outputs of each)
supplementary_tables/   Supplementary Tables S1–S14 (CSV) produced by the scripts
```

## Requirements

- R ≥ 4.5 with TwoSampleMR 0.7.4, ieugwasr 1.1.0, MRPRESSO 1.0, MendelianRandomization, mr.raps 0.4.3,
  MRlap 0.0.3.3, MRcML, MVMR, cause 1.2.0 (script 27 applies a small compatibility patch for loo ≥ 2.8),
  coloc, data.table, dplyr, ggplot2, jsonlite.
- PLINK 1.9 and the 1000 Genomes European reference panel (`EUR.bed/bim/fam`).
- LD score regression European LD scores (`eur_w_ld_chr/`, with `w_hm3.snplist`).
- An OpenGWAS access token in the environment variable `OPENGWAS_JWT` (free registration at
  https://api.opengwas.io). **No token is stored in this repository.**

## Input data (all public; not redistributed here)

| Data | Source |
|---|---|
| Periodontitis, FinnGen (Salminen et al. 2025) | https://storage.googleapis.com/fg-publication-green-public/F_2023_026_20250625/summary_statistics_periodontitis.zip |
| AD, Bellenguez et al. 2022 | GWAS Catalog GCST90027158 (harmonised file) / OpenGWAS ebi-a-GCST90027158 |
| AD, Kunkle et al. 2019 | OpenGWAS ieu-b-2 |
| Gut microbiota, MiBioGen (Kurilshikov et al. 2021) | OpenGWAS ebi-a-GCST90016908 … ebi-a-GCST90017118 |
| Periodontitis, GLIDE (Shungin et al. 2019) | doi:10.5523/bris.2j2rqgzedxlq02oqbb4vmycnc2 (`EUR_perio_excl_HCHSSOL.txt`) |
| Smoking initiation (Liu et al. 2019, GSCAN excl. 23andMe) | GWAS Catalog GCST007474 |
| Type 2 diabetes (Xue et al. 2018) | GWAS Catalog GCST006867 |
| Body mass index (Yengo et al. 2018) | GWAS Catalog GCST006900 |
| LDL cholesterol (Graham et al. 2021, European) | GWAS Catalog GCST90239658 |
| Years of education (UK Biobank; Loh et al. 2018) | GWAS Catalog GCST90029013 |
| Type 1 diabetes (Chiou et al. 2021) | GWAS Catalog GCST90014023 |
| 1000 Genomes EUR reference | http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz |

Expected local paths are given in `scripts/README.txt`. Run all scripts from the project root
(the folder containing `data/`, `results/`, `scripts/`, `ext_data/`).

## Licence

MIT (see `LICENSE`). Third-party summary statistics remain under their original terms.
