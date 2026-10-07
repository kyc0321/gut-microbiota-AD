## 14_taxa_independence_meff_relaxed.R
## Redo of 06_taxa_independence_meff.R using the NEW 50-SNP relaxed periodontitis
## instrument (data/step1_chunks_relaxed/, from 11_step1_relaxed.R) instead of the
## original 15-SNP instrument, for consistency with the rest of the v2 analysis and
## improved precision of the 211x211 correlation-matrix estimate (more SNPs = more
## data points per pairwise correlation).
suppressPackageStartupMessages({ library(data.table); library(dplyr) })

chunk_dir <- "data/step1_chunks_relaxed"
files <- list.files(chunk_dir, pattern = "\\.rds$", full.names = TRUE)
cat("Loading", length(files), "chunk files...\n")
dat <- rbindlist(lapply(files, readRDS), fill = TRUE)
cat("Rows (SNP x taxon pairs):", nrow(dat), "\n")
cat("Unique SNPs:", length(unique(dat$SNP)), "\n")
cat("Unique taxa (id.outcome):", length(unique(dat$id.outcome)), "\n")

dat[, z := beta.outcome / se.outcome]
wide <- dcast(dat, SNP ~ id.outcome, value.var = "z", fun.aggregate = mean)
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- wide$SNP
cat("Matrix dim (SNP x taxon):", paste(dim(mat), collapse=" x "), "\n")

miss_frac <- colMeans(is.na(mat))
keep <- miss_frac <= 0.2
cat("Taxa dropped for >20% missing SNPs:", sum(!keep), "of", ncol(mat), "\n")
mat <- mat[, keep]

R <- suppressWarnings(cor(mat, use = "pairwise.complete.obs"))
M <- ncol(R)
cat("Final taxon x taxon correlation matrix:", M, "x", M, "\n")

eig <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
eig[eig < 0] <- 0
var_lambda <- sum((eig - 1)^2) / (M - 1)
Meff_Nyholt <- 1 + (M - 1) * (1 - var_lambda / M)
Meff_LiJi <- sum(ifelse(eig >= 1, 1, eig - floor(eig)))
Meff_Galwey <- (sum(sqrt(eig)))^2 / sum(eig)

cat("\n=== Effective number of independent tests (raw M =", M, ") ===\n")
cat(sprintf("Nyholt (2004):  Meff = %.1f  -> Bonferroni alpha = %.2e\n", Meff_Nyholt, 0.05/Meff_Nyholt))
cat(sprintf("Li & Ji (2005): Meff = %.1f  -> Bonferroni alpha = %.2e\n", Meff_LiJi, 0.05/Meff_LiJi))
cat(sprintf("Galwey (2009):  Meff = %.1f  -> Bonferroni alpha = %.2e\n", Meff_Galwey, 0.05/Meff_Galwey))
cat(sprintf("(For comparison, naive Bonferroni on raw M: alpha = %.2e)\n", 0.05/M))

R_upper <- R
R_upper[lower.tri(R_upper, diag = TRUE)] <- NA
dup_idx <- which(abs(R_upper) > 0.95, arr.ind = TRUE)
dup_pairs <- data.frame(
  taxon_A = rownames(R)[dup_idx[,1]], taxon_B = colnames(R)[dup_idx[,2]], r = R[dup_idx]
)
id2name <- dat[!duplicated(id.outcome), .(id.outcome, originalname.outcome)]
dup_pairs <- dup_pairs %>%
  left_join(id2name, by = c("taxon_A" = "id.outcome")) %>% rename(name_A = originalname.outcome) %>%
  left_join(id2name, by = c("taxon_B" = "id.outcome")) %>% rename(name_B = originalname.outcome) %>%
  arrange(desc(abs(r)))

cat("\n=== Near-duplicate taxa pairs (|r| > 0.95), n =", nrow(dup_pairs), "===\n")
print(head(dup_pairs, 15), row.names = FALSE)

saveRDS(R, "data/taxa_correlation_matrix_relaxed.rds")
fwrite(dup_pairs, "results/taxa_near_duplicate_pairs_relaxed.csv")
summary_tbl <- data.frame(
  method = c("Naive (raw M)", "Nyholt (2004)", "Li & Ji (2005)", "Galwey (2009)"),
  Meff = c(M, round(Meff_Nyholt,1), round(Meff_LiJi,1), round(Meff_Galwey,1)),
  bonferroni_alpha = c(0.05/M, 0.05/Meff_Nyholt, 0.05/Meff_LiJi, 0.05/Meff_Galwey)
)
fwrite(summary_tbl, "results/taxa_meff_summary_relaxed.csv")
cat("\n[OK] Saved *_relaxed outputs.\n")
