## 20_pleiotropy_lookup.R  (post-AHG editor comment 2)
## Look up every primary instrument SNP (P<1e-5, n=50) in OpenGWAS PheWAS
## (P<5e-8) and the NHGRI-EBI GWAS Catalog; flag associations with
## established AD risk factors / AD itself.
suppressPackageStartupMessages({ library(ieugwasr); library(data.table); library(dplyr); library(jsonlite) })
source("scripts/_helpers.R")
iv <- readRDS("data/perio_instruments_finngen2025_relaxed.rds")
snps <- iv$SNP

## 1. OpenGWAS PheWAS
pw <- rbindlist(lapply(split(snps, ceiling(seq_along(snps)/10)), function(s)
  tryCatch(as.data.table(phewas(s, pval = 5e-8)), error = function(e) {message(e$message); NULL})), fill = TRUE)
fwrite(pw, "results/R_AHG/phewas_opengwas_raw.csv")

## 2. GWAS Catalog (REST API)
gc_one <- function(rs) {
  out <- list(); pg <- 0
  repeat {
    u <- sprintf("https://www.ebi.ac.uk/gwas/rest/api/v2/associations?rs_id=%s&size=500&page=%d", rs, pg)
    j <- tryCatch(fromJSON(u, simplifyVector = FALSE), error = function(e) NULL)
    a <- j$`_embedded`$associations
    if (is.null(a)) break
    out[[length(out) + 1]] <- rbindlist(lapply(a, function(x) data.table(
      SNP = rs, pvalue = as.numeric(x$p_value),
      trait = paste(c(vapply(x$efo_traits, function(t) t$efo_trait, ""), unlist(x$reported_trait)), collapse = "; "))))
    pg <- pg + 1
    if (pg >= j$page$totalPages) break
  }
  rbindlist(out)
}
gc <- rbindlist(lapply(snps, function(s) { Sys.sleep(0.2); gc_one(s) }), fill = TRUE)
gc <- if (nrow(gc)) gc[pvalue < 5e-8] else data.table(SNP=character(), pvalue=numeric(), trait=character())
fwrite(gc, "results/R_AHG/gwascatalog_raw.csv")

## 3. Flag AD-relevant traits (pre-specified keyword list)
kw <- paste(c("alzheim","dementia","cognit","educat","intellig","apolipo","cholesterol","ldl","hdl","triglycer",
              "lipid","body mass","bmi","obes","waist","diabet","glucose","hba1c","insulin","blood pressure","hypertens",
              "smok","cigarette","alcohol","stroke","coronary","myocard","cardiovasc","sleep","depress",
              "c-reactive","crp","hearing","physical activ"), collapse = "|")
op <- if (nrow(pw)) pw[, .(SNP = rsid, trait, p, source = "OpenGWAS")] else data.table()
gcx <- if (nrow(gc)) gc[, .(SNP, trait, p = pvalue, source = "GWAS Catalog")] else data.table()
all <- rbind(op, gcx, fill = TRUE)
all <- all[!grepl("periodont|gingiv|tooth|dental|mouth|oral", trait, ignore.case = TRUE)]
all[, ad_relevant := grepl(kw, trait, ignore.case = TRUE)]
fwrite(all, "results/R_AHG/Table_SX_instrument_pleiotropy_lookup.csv")

summ <- all[, .(n_assoc = .N, n_ad_relevant = sum(ad_relevant),
                ad_traits = paste(unique(trait[ad_relevant]), collapse = " | "),
                other_traits = paste(head(unique(trait[!ad_relevant]), 6), collapse = " | ")), by = SNP]
summ <- merge(data.table(SNP = snps, chr = iv$chr.exposure, pos = iv$pos.exposure), summ, by = "SNP", all.x = TRUE)
summ[is.na(n_assoc), `:=`(n_assoc = 0L, n_ad_relevant = 0L)]
fwrite(summ, "results/R_AHG/instrument_pleiotropy_summary.csv")
cat("SNPs with any other GW-significant association:", sum(summ$n_assoc > 0), "/", nrow(summ), "\n")
cat("SNPs with AD-relevant association:", sum(summ$n_ad_relevant > 0), "\n")
print(summ[n_assoc > 0, .(SNP, chr, n_assoc, n_ad_relevant, ad_traits = substr(ad_traits, 1, 90))])
