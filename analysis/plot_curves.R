# ==============================================================================
# plot_curves.R — shared figure code for Fig. 2 and Fig. 3
#   Panel A: mean cumulative richness ± 1 SD, Oracle as dashed line
#   Panel B: BAS vs Random across simulation runs at selected n, one-sided paired t-test
# ==============================================================================
suppressMessages({library(dplyr); library(ggplot2); library(patchwork)})
bas_cols <- c("Oracle" = "#009E73", "BAS" = "#0072B2", "Random" = "#D55E00")
bas_lts  <- c("Oracle" = "dashed", "BAS" = "solid", "Random" = "solid")
sig_star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

plot_bas_figure <- function(df, target_box, xlim, file, width = 10, height = 12, dpi = 300) {
  lv <- c("Oracle", "BAS", "Random"); df <- df %>% filter(Method %in% lv)
  summ <- df %>% filter(n_samples >= xlim[1], n_samples <= xlim[2]) %>% group_by(n_samples, Method) %>%
    summarise(Mean = mean(Richness), SD = sd(Richness), .groups = "drop") %>% mutate(Method = factor(Method, lv))
  pA <- ggplot(summ, aes(n_samples, Mean, color = Method, fill = Method, linetype = Method)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.15, linetype = 0) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = bas_cols) + scale_fill_manual(values = bas_cols) + scale_linetype_manual(values = bas_lts) +
    labs(tag = "A", x = "Number of Samples", y = "Cumulative Richness") +
    theme_minimal(base_size = 14) + theme(legend.position = "bottom", legend.title = element_blank())
  box <- df %>% filter(n_samples %in% target_box, Method %in% c("BAS", "Random")) %>%
    mutate(n_label = factor(paste("n =", n_samples), paste("n =", target_box)), Method = factor(Method, c("BAS", "Random")))
  anno <- box %>% group_by(n_label) %>% group_modify(function(d, k) {
    w <- inner_join(filter(d, Method == "BAS") %>% select(SimID, BAS = Richness), filter(d, Method == "Random") %>% select(SimID, Random = Richness), by = "SimID")
    p <- t.test(w$BAS, w$Random, paired = TRUE, alternative = "greater")$p.value
    data.frame(x = 1, xend = 2, y = max(d$Richness) + diff(range(d$Richness)) * 0.1, label = sig_star(p)) }) %>% ungroup()
  pB <- ggplot(box, aes(Method, Richness)) +
    geom_point(aes(color = Method), position = position_jitterdodge(0.2, 0, 0.6), size = 0.5, alpha = 0.3) +
    geom_boxplot(aes(color = Method, fill = Method), alpha = 0.6, width = 0.6, outlier.shape = NA) +
    facet_wrap(~n_label, scales = "free_y", nrow = 1) +
    scale_color_manual(values = bas_cols) + scale_fill_manual(values = bas_cols) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    geom_segment(data = anno, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.5) +
    geom_text(data = anno, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, vjust = -0.2, size = 3.5) +
    labs(tag = "B", x = NULL, y = "Cumulative Richness") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom", legend.title = element_blank(), strip.text = element_text(face = "bold"),
          panel.border = element_rect(color = "grey", fill = NA), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(t = 20))
  p <- pA / pB + plot_layout(heights = c(1, 1.2))
  ggsave(paste0(file, ".png"), p, width = width, height = height, dpi = dpi)
  ggsave(paste0(file, ".pdf"), p, width = width, height = height)
  invisible(p)
}

summarise_bas <- function(df, ns) {
  wide <- df %>% filter(n_samples %in% ns, Method %in% c("Oracle", "BAS", "Random")) %>% tidyr::pivot_wider(names_from = Method, values_from = Richness)
  wide %>% group_by(n_samples) %>% summarise(
    BAS_mean = mean(BAS), BAS_sd = sd(BAS), Random_mean = mean(Random), Random_sd = sd(Random), Oracle_mean = mean(Oracle),
    gain = mean(BAS - Random), pct_of_oracle = 100 * mean(BAS) / mean(Oracle),
    gap_closed_pct = 100 * mean((BAS - Random) / (Oracle - Random)),
    p_one_sided = stats::t.test(BAS, Random, paired = TRUE, alternative = "greater")$p.value, .groups = "drop")
}
