# ==============================================================================
# 03b: Re-draw Figure 3A-B (BAS vs Random vs Oracle) from saved results
#
# Plot-only companion to 03_BAS_simulation_Fig3AB.R.  Reads
# simulation_results_Fig3AB_v2.rds (100 simulations, produced by 03) and
# reproduces Steps 3-5 of that script verbatim, so the figure can be
# regenerated without re-running the ~30 min simulation.
#
# Output:
#   BAS_GlobalFungi_Plot_v2_withOracle.png   (Fig. 3A-B, main text)
#   BAS_GlobalFungi_Plot_v2_OracleLine.pdf   (panel A alone; Fig. S2)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

results_df <- readRDS("simulation_results_Fig3AB_v2.rds")

# ==============================================================================
# Step 3: Figure A — line plot
# ==============================================================================

target_box    <- c(1200, 1500, 2000, 2500, 3000)
target_levels <- c("Oracle", "BAS", "Random")

summ_line <- results_df %>%
  group_by(n_samples, Method) %>%
  summarise(Mean = mean(Richness), SD = sd(Richness),
            ymin = Mean - SD, ymax = Mean + SD, .groups = "drop") %>%
  mutate(Method = factor(Method, levels = target_levels))

df_box <- results_df %>%
  filter(n_samples %in% target_box) %>%
  mutate(
    n_label = factor(paste("n =", n_samples), levels = paste("n =", target_box)),
    Method  = factor(Method, levels = target_levels)
  )

cols      <- c("Oracle" = "#009E73", "BAS" = "#0072B2", "Random" = "#D55E00")
fills     <- c("Oracle" = "#009E73", "BAS" = "#0072B2", "Random" = "#D55E00")
linetypes <- c("Oracle" = "dashed",  "BAS" = "solid",   "Random" = "solid")

p_line <- ggplot(summ_line, aes(x = n_samples, y = Mean, color = Method, fill = Method, linetype = Method)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, linetype = 0) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = fills) +
  scale_linetype_manual(values = linetypes) +
  labs(tag = "A", x = "Number of Samples", y = "Cumulative Richness",
       subtitle = "Dashed line = theoretical maximum (Oracle)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", legend.title = element_blank())

# ==============================================================================
# Step 4: Figure B — box plot with significance brackets
# ==============================================================================

get_sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

anno_df <- data.frame()
for (n_lab in levels(df_box$n_label)) {
  dat_n  <- df_box %>% filter(n_label == n_lab)
  y_span <- diff(range(dat_n$Richness))
  gap    <- y_span * 0.1
  y_base <- max(dat_n$Richness)

  # BAS (x=2) vs Random (x=3): 下段ブラケット
  paired_br <- inner_join(
    dat_n %>% filter(Method == "BAS")    %>% dplyr::select(SimID, BAS = Richness),
    dat_n %>% filter(Method == "Random") %>% dplyr::select(SimID, Random = Richness),
    by = "SimID"
  )
  if (nrow(paired_br) > 1) {
    p_br <- t.test(paired_br$BAS, paired_br$Random, paired = TRUE, alt = "greater")$p.value
    y_br <- y_base + gap
    anno_df <- rbind(anno_df,
      data.frame(n_label = n_lab, x = 2, xend = 3, y = y_br, label = get_sig(p_br))
    )
  }

  # Oracle (x=1) vs BAS (x=2): 上段ブラケット
  paired_ob <- inner_join(
    dat_n %>% filter(Method == "Oracle") %>% dplyr::select(SimID, Oracle = Richness),
    dat_n %>% filter(Method == "BAS")    %>% dplyr::select(SimID, BAS = Richness),
    by = "SimID"
  )
  if (nrow(paired_ob) > 1) {
    p_ob <- t.test(paired_ob$Oracle, paired_ob$BAS, paired = TRUE, alt = "greater")$p.value
    y_ob <- y_base + gap * 2.2
    anno_df <- rbind(anno_df,
      data.frame(n_label = n_lab, x = 1, xend = 2, y = y_ob, label = get_sig(p_ob))
    )
  }
}

p_box <- ggplot(df_box, aes(x = Method, y = Richness)) +
  geom_point(aes(color = Method),
             position = position_jitterdodge(0.2, 0, 0.6), size = 0.5, alpha = 0.3) +
  geom_boxplot(aes(color = Method, fill = Method), alpha = 0.6, width = 0.6, outlier.shape = NA) +
  facet_wrap(~n_label, scales = "free_y", nrow = 1) +
  scale_color_manual(values = cols, drop = FALSE) +
  scale_fill_manual(values = fills, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(tag = "B", x = NULL, y = "Cumulative Richness") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_blank(),
    strip.text       = element_text(face = "bold"),
    panel.border     = element_rect(color = "grey", fill = NA),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    plot.margin      = margin(t = 20)
  )

if (nrow(anno_df) > 0) {
  p_box <- p_box +
    geom_segment(data = anno_df, aes(x = x,    xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.5) +
    geom_segment(data = anno_df, aes(x = x,    xend = x,    y = y, yend = y - (y * 0.002)),
                 inherit.aes = FALSE, linewidth = 0.5) +
    geom_segment(data = anno_df, aes(x = xend, xend = xend, y = y, yend = y - (y * 0.002)),
                 inherit.aes = FALSE, linewidth = 0.5) +
    geom_text(data = anno_df, aes(x = (x + xend) / 2, y = y, label = label),
              inherit.aes = FALSE, vjust = -0.2, size = 3.5)
}

# ==============================================================================
# Step 5: Combine and save
# ==============================================================================

p_comb <- p_line / p_box + plot_layout(heights = c(1, 1.2))
ggsave("BAS_GlobalFungi_Plot_v2_withOracle.png", p_comb, width = 10, height = 14, dpi = 300)
ggsave("BAS_GlobalFungi_Plot_v2_OracleLine.pdf", p_line, width = 8, height = 6)
message("Saved: BAS_GlobalFungi_Plot_v2_withOracle.png, BAS_GlobalFungi_Plot_v2_OracleLine.pdf")
