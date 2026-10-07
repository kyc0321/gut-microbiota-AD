## 22a: format FinnGen periodontitis (Salminen 2025) for MRlap — HapMap3 SNPs + instrument SNPs
suppressPackageStartupMessages(library(data.table))
hm3 <- fread("ext_data/eur_w_ld_chr/w_hm3.snplist")
iv  <- readRDS("data/perio_instruments_finngen2025_relaxed.rds")
dt  <- fread("../FinGen Download/salminen_2025_sumstats/periodontitis_diagnosis_results.txt.gz",
             select = c("chrom","pos","ref","alt","rsids","beta","sebeta","pval"))
dt[, rsid := sub(",.*", "", rsids)]
dt <- dt[rsid %in% c(hm3$SNP, iv$SNP)][!duplicated(rsid)]
exp <- dt[, .(rsid, chr = chrom, pos, alt, ref, beta, se = sebeta)]
saveRDS(exp, "ext_data/mrlap_exposure_perio.rds")
cat("Exposure rows:", nrow(exp), " instruments present:", sum(iv$SNP %in% exp$rsid), "\n")
