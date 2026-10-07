## 17_fig3_total_effect_composite.R -- Clean 4-panel Figure 3 (total effect of
## periodontitis on AD) for the Annals of Human Genetics submission: short,
## untruncated axis labels; no stray as.factor legends; panels a-d via patchwork.
suppressPackageStartupMessages({ library(TwoSampleMR); library(ggplot2); library(patchwork) })

obj <- readRDS("data/PRIMARY_total_perio_ad_relaxed.rds")
dat <- obj$dat; res <- obj$res; loo <- obj$loo
dat$exposure <- "Periodontitis"; dat$outcome <- "Alzheimer's disease"
res$exposure <- "Periodontitis"; res$outcome <- "Alzheimer's disease"
loo$exposure <- "Periodontitis"; loo$outcome <- "Alzheimer's disease"
res_single <- mr_singlesnp(dat)

th <- theme_bw(base_size = 9) + theme(plot.title = element_blank())

p_a <- mr_scatter_plot(res, dat)[[1]] + th +
  labs(x = "SNP effect on periodontitis", y = "SNP effect on Alzheimer's disease") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(colour = guide_legend(ncol = 2))
p_b <- mr_forest_plot(res_single)[[1]] + th +
  labs(x = "MR effect size (log odds ratio)") + theme(legend.position = "none",
  axis.text.y = element_text(size = 5))
p_c <- mr_leaveoneout_plot(loo)[[1]] + th +
  labs(x = "IVW estimate with SNP omitted (log odds ratio)") + theme(legend.position = "none",
  axis.text.y = element_text(size = 5))
p_d <- mr_funnel_plot(res_single)[[1]] + th +
  theme(legend.position = "bottom", legend.title = element_blank())

fig <- (p_a | p_b) / (p_c | p_d) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 12))
dir.create("figures_ahg", showWarnings = FALSE)
ggsave("figures_ahg/Figure3_TotalEffect.pdf", fig, width = 7.5, height = 10)
ggsave("figures_ahg/Figure3_TotalEffect.png", fig, width = 7.5, height = 10, dpi = 300)
cat("[OK] Figure 3 composite saved\n")
