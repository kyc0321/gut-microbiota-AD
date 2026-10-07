## 19_table_s7_snp_associations.R -- Table S7: harmonised SNP-exposure and
## SNP-outcome associations underlying the primary total-effect analysis.
suppressPackageStartupMessages(library(data.table))
d <- as.data.table(readRDS("data/PRIMARY_total_perio_ad_relaxed.rds")$dat)
out <- d[, .(SNP, chr = chr.exposure, pos = pos.exposure,
             effect_allele = effect_allele.exposure, other_allele = other_allele.exposure,
             EAF_exposure = eaf.exposure,
             beta_periodontitis = beta.exposure, se_periodontitis = se.exposure, p_periodontitis = pval.exposure,
             beta_AD = beta.outcome, se_AD = se.outcome, p_AD = pval.outcome,
             palindromic, used_in_analysis = mr_keep)]
setorder(out, chr, pos)
fwrite(out, "results/Table_S7_harmonised_SNP_associations.csv")
cat("rows:", nrow(out), " used:", sum(out$used_in_analysis), "\n")
