## 25_coloc_gw_loci.R — colocalisation (coloc.abf) of periodontitis and AD at the three
## genome-wide significant periodontitis loci (FST, HLA, MFHAS1), +/-500 kb (GRCh38).
## H3 = distinct causal variants (LD-confounding), H4 = shared causal variant.
suppressPackageStartupMessages({ library(coloc); library(data.table) })
loci <- fread("data/perio_instruments_finngen2025.csv")[, .(SNP, chr = chr.exposure, pos = pos.exposure)]
loci[, gene := c("FST", "HLA region", "MFHAS1")[match(chr, c(5, 6, 8))]]
W <- 5e5
pe <- fread("../FinGen Download/salminen_2025_sumstats/periodontitis_diagnosis_results.txt.gz",
            select = c("chrom","pos","ref","alt","rsids","beta","sebeta","af_alt"))
ad <- fread("ext_data/bellenguez/35379992-GCST90027158-MONDO_0004975.h.tsv.gz",
            select = c("hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","standard_error"))
setnames(ad, c("chrom","pos","oa","ea","b_ad","se_ad")); ad[, chrom := as.integer(chrom)]
res <- list()
for (i in seq_len(nrow(loci))) {
  L <- loci[i]
  a <- pe[chrom == L$chr & abs(pos - L$pos) <= W]
  b <- ad[chrom == L$chr & abs(pos - L$pos) <= W & !is.na(b_ad)]
  m <- merge(a, b, by = c("chrom","pos"))
  m[, b_ad2 := fifelse(ea == alt & oa == ref, b_ad, fifelse(ea == ref & oa == alt, -b_ad, NA_real_))]
  m <- m[!is.na(b_ad2) & af_alt > 0.01 & af_alt < 0.99]
  m[, snp := paste(chrom, pos, ref, alt, sep = ":")]
  m <- m[!duplicated(snp)]
  d1 <- list(beta = m$beta, varbeta = m$sebeta^2, snp = m$snp, MAF = pmin(m$af_alt, 1 - m$af_alt), type = "cc", s = 38157/250000)
  d2 <- list(beta = m$b_ad2, varbeta = m$se_ad^2, snp = m$snp, MAF = pmin(m$af_alt, 1 - m$af_alt), type = "cc", s = (39106 + 46828)/487511)
  for (p12 in c(1e-5, 5e-6)) {
    cc <- suppressMessages(coloc.abf(d1, d2, p12 = p12))
    pp <- cc$summary
    res[[length(res) + 1]] <- data.table(locus = L$gene, lead = L$SNP, chr = L$chr, nsnps = pp["nsnps"], p12,
      PP.H0 = pp["PP.H0.abf"], PP.H1 = pp["PP.H1.abf"], PP.H2 = pp["PP.H2.abf"], PP.H3 = pp["PP.H3.abf"], PP.H4 = pp["PP.H4.abf"],
      min_p_AD_region = min(2 * pnorm(-abs(m$b_ad2 / m$se_ad))))
  }
}
res <- rbindlist(res); fwrite(res, "results/R_AHG/coloc_gw_loci.csv")
print(res[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
