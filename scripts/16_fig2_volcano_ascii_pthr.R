## 16_fig2_volcano_ascii_pthr.R
## Regenerate the Step-1 volcano plot (current manuscript Figure 2) with the
## corrected threshold notation: "P<X" instead of "a~X" (the earlier ASCII
## rendering still used "a" -- alpha -- as a prefix for the significance
## threshold, clashing with alpha's other use in this paper as the
## periodontitis-on-taxon effect-size symbol; main text/legend already fixed
## to "P<..." wording, this brings the baked-in image annotation into line).
## Uses default pdf() device + pure-ASCII text (no Unicode/Greek), matching
## this machine's established workaround for cairo_pdf failures.
suppressPackageStartupMessages({ library(ggplot2) })

step1 <- read.csv("results/EXPLORATORY_step1_alpha_IVW_relaxed.csv")
step1$neglog10p <- -log10(step1$alpha_p)
step1$highlight <- ifelse(grepl("Mollicutes|Tenericutes", step1$taxon),
                           "Prior candidate\n(Mollicutes/Tenericutes)", "Other taxa")

alpha_naive <- 0.05/211
alpha_meff  <- 0.05/33  # M_eff = 32-33 (Li&Ji/Galwey); matches P<1.5e-3 reported in text

p_volcano <- ggplot(step1, aes(x = alpha, y = neglog10p)) +
  geom_point(aes(color = highlight, size = highlight), alpha = 0.75) +
  scale_color_manual(values = c("Other taxa" = "grey40", "Prior candidate\n(Mollicutes/Tenericutes)" = "firebrick")) +
  scale_size_manual(values = c("Other taxa" = 1.6, "Prior candidate\n(Mollicutes/Tenericutes)" = 3.2)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = -log10(alpha_meff), linetype = "dashed", color = "steelblue") +
  geom_hline(yintercept = -log10(alpha_naive), linetype = "solid", color = "black") +
  annotate("text", x = min(step1$alpha, na.rm=TRUE), y = -log10(0.05)+0.15,
           label = "nominal P=0.05", hjust=0, size=3, color="grey50") +
  annotate("text", x = min(step1$alpha, na.rm=TRUE), y = -log10(alpha_meff)+0.15,
           label = "M_eff-corrected (P<1.5x10^-3)", hjust=0, size=3, color="steelblue") +
  annotate("text", x = min(step1$alpha, na.rm=TRUE), y = -log10(alpha_naive)+0.15,
           label = "naive Bonferroni (P<2.4x10^-4)", hjust=0, size=3, color="black") +
  labs(title = "Step 1: Periodontitis (Salminen 2025) -> 211 gut microbiota taxa",
       subtitle = "No taxon survives correction; prior Mollicutes/Tenericutes candidate not replicated",
       x = "alpha (IVW effect of periodontitis on taxon abundance)",
       y = "-log10(P)", color = NULL, size = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")

ggsave("figures_v2_ascii/Figure2_Step1_volcano.pdf", p_volcano, width = 7.5, height = 6.5, device = pdf)
cat("[OK] Figure 2 volcano (ASCII, P<threshold) saved.\n")
