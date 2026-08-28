## 03_fdr_correction.R -- Benjamini-Hochberg FDR correction across independent
## signal clusters, collapsing taxonomically-nested duplicate signals.
##
## The MiBioGen taxa are taxonomically nested (phylum/class/order/family/genus
## derived from the same underlying compositional abundance profiles), so a
## lineage instrumented by an identical SNP set at two adjacent taxonomic
## levels (e.g. a genus that is the sole member of its family) yields
## numerically identical MR estimates reported under two different taxon
## names. Treating these as two independent tests would inflate the
## multiple-testing correction; this script collapses them to one signal
## cluster per unique (beta, se, p) triple before applying BH-FDR, and reports
## the resulting q-value at every nested taxonomic level within a cluster.

source("00_setup.R")

ivw_tbl <- fread("results/ivw_by_taxon.csv")
n_tested <- nrow(ivw_tbl)

# cluster key: numerically identical estimates (rounded) = same underlying signal
ivw_tbl[, cluster_key := paste(round(beta, 8), round(beta_se, 8), round(beta_p, 10))]
cluster_sizes <- ivw_tbl[, .N, by = cluster_key]
n_clusters <- nrow(cluster_sizes)
cat(sprintf("Tested taxa: %d  |  independent signal clusters: %d  |  duplicate rows collapsed: %d\n",
            n_tested, n_clusters, n_tested - n_clusters))

# one representative P per cluster (all members share the same P by construction)
cluster_p <- ivw_tbl[, .(beta_p = beta_p[1]), by = cluster_key]
cluster_p <- cluster_p[order(beta_p)]
cluster_p[, rank := .I]
cluster_p[, raw_q := beta_p * n_clusters / rank]
# standard BH step-up: running minimum from largest P to smallest
cluster_p <- cluster_p[order(-beta_p)]
cluster_p[, q := cummin(raw_q)]
cluster_p <- cluster_p[order(beta_p)]

ivw_tbl <- merge(ivw_tbl, cluster_p[, .(cluster_key, q)], by = "cluster_key")
ivw_tbl <- ivw_tbl[order(beta_p)]
fwrite(ivw_tbl[, .(taxon, taxon_id, beta, beta_se, beta_p, beta_nsnp, cluster_key, q)],
       "results/ivw_by_taxon_with_fdr.csv")

nominal <- ivw_tbl[beta_p < 0.05]
cat(sprintf("\nNominal taxa (P<0.05): %d raw rows, %d independent clusters\n",
            nrow(nominal), uniqueN(nominal$cluster_key)))
cat("Smallest q:", round(min(ivw_tbl$q), 3), "\n")
cat("[OK] 03_fdr_correction.R complete.\n")
