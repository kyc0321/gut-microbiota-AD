## 28_glide_replication.R — replication with an independent exposure GWAS:
## GLIDE clinical periodontitis, European ancestry excluding HCHS/SOL (Shungin et al. 2019;
## 12,289 cases / 22,326 controls). Markers (hg19 chr:pos) mapped to rsIDs via the 1000 Genomes
## EUR reference; palindromic SNPs dropped (no allele frequencies released); variants analysed in
## <80% of the maximum N dropped. Instruments P<1e-5, clumped locally (r2<0.001, 10 Mb).
suppressPackageStartupMessages({ library(TwoSampleMR); library(mr.raps); library(data.table) })
g <- fread("ext_data/glide/EUR_perio_excl_HCHSSOL.txt")
setnames(g, c("MarkerName","Allele1","Allele2","Effect","StdErr","P-value"), c("id","a1","a2","beta","se","p"))
g <- g[N >= 0.8 * max(N)]
g[, `:=`(a1 = toupper(a1), a2 = toupper(a2))]
g <- g[!((a1 %in% c("A","T") & a2 %in% c("A","T")) | (a1 %in% c("C","G") & a2 %in% c("C","G")))]
bim <- fread("ext_data/ld_ref/EUR.bim", col.names = c("chr","rsid","cm","pos","b1","b2"))
bim[, id := paste(chr, pos, sep = ":")]
g <- merge(g, bim[, .(id, rsid, b1, b2)], by = "id")
g <- g[(a1 == b1 & a2 == b2) | (a1 == b2 & a2 == b1)][!duplicated(rsid)]
cat("GLIDE variants mapped:", nrow(g), "; max N:", max(g$N), "\n")
cand <- g[p < 1e-5]
fwrite(cand[, .(SNP = rsid, P = p)], "ext_data/glide_clump_input.txt", sep = "\t")
system2("ext_data/plink/plink", c("--bfile", "ext_data/ld_ref/EUR", "--clump", "ext_data/glide_clump_input.txt",
        "--clump-p1", "1e-5", "--clump-r2", "0.001", "--clump-kb", "10000", "--out", "ext_data/glide_clump"), stdout = FALSE)
ivs <- fread("ext_data/glide_clump.clumped")$SNP
iv <- cand[rsid %in% ivs]
cat("GLIDE instruments:", nrow(iv), " mean F:", round(mean((iv$beta/iv$se)^2), 1), "\n")
exp_dat <- data.frame(SNP = iv$rsid, beta.exposure = iv$beta, se.exposure = iv$se, effect_allele.exposure = iv$a1,
                      other_allele.exposure = iv$a2, pval.exposure = iv$p, exposure = "Periodontitis (GLIDE)",
                      id.exposure = "GLIDE", samplesize.exposure = iv$N, eaf.exposure = NA_real_, stringsAsFactors = FALSE)
ad <- fread("ext_data/bellenguez/35379992-GCST90027158-MONDO_0004975.h.tsv.gz",
            select = c("hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","p_value","hm_effect_allele_frequency"))
ad <- ad[hm_rsid %in% iv$rsid & !is.na(hm_beta)][!duplicated(hm_rsid)]
out_dat <- data.frame(SNP = ad$hm_rsid, beta.outcome = ad$hm_beta, se.outcome = ad$standard_error,
                      effect_allele.outcome = ad$hm_effect_allele, other_allele.outcome = ad$hm_other_allele,
                      eaf.outcome = ad$hm_effect_allele_frequency, pval.outcome = ad$p_value,
                      outcome = "AD (Bellenguez 2022)", id.outcome = "GCST90027158", stringsAsFactors = FALSE)
dat <- harmonise_data(exp_dat, out_dat, action = 1); dat <- subset(dat, mr_keep)   # palindromes already removed
res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
rp <- mr.raps.overdispersed.robust(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, loss.function = "huber")
res <- rbind(as.data.table(res)[, .(method, nsnp, b, se, pval)],
             data.table(method = "MR-RAPS", nsnp = nrow(dat), b = rp$beta.hat, se = rp$beta.se, pval = 2*pnorm(-abs(rp$beta.hat/rp$beta.se))))
res[, `:=`(OR = exp(b), LCI = exp(b - 1.96*se), UCI = exp(b + 1.96*se), OR_doubling = exp(b*log(2)))]
het <- mr_heterogeneity(dat, method_list = "mr_ivw"); ple <- mr_pleiotropy_test(dat)
fwrite(res, "results/R_AHG/glide_replication.csv"); fwrite(dat, "results/R_AHG/glide_harmonised.csv")
print(res); print(het[, c("Q","Q_df","Q_pval")]); print(ple[, c("egger_intercept","pval")])
## overlap check: FinnGen instruments vs GLIDE effects (concordance of direction)
fi <- readRDS("data/perio_instruments_finngen2025_relaxed.rds")
cc <- merge(as.data.table(fi)[, .(rsid = SNP, ea = effect_allele.exposure, oa = other_allele.exposure, bf = beta.exposure)],
            g[, .(rsid, a1, a2, bg = beta, pg = p)], by = "rsid")
cc[, bg2 := fifelse(a1 == ea, bg, fifelse(a1 == oa, -bg, NA_real_))]
cat("FinnGen instruments found in GLIDE:", nrow(cc), "; same direction:", sum(sign(cc$bf) == sign(cc$bg2), na.rm = TRUE),
    "; binomial P:", signif(binom.test(sum(sign(cc$bf) == sign(cc$bg2), na.rm = TRUE), sum(!is.na(cc$bg2)))$p.value, 2), "\n")
fwrite(cc, "results/R_AHG/finngen_ivs_in_glide.csv")
