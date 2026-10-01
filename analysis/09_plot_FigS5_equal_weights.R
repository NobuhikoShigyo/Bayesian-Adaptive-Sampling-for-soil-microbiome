# ==============================================================================
# 09: Figure S5 — sensitivity to the acquisition weights
#   Tuned weights (main analysis) vs equal weights (1/3, 1/3, 1/3) vs Random,
#   on the evaluation pilot sets of 01 (local) and 03 (global).
# Input : simulation_results_Fig2AB.rds, simulation_results_Fig3AB.rds (both contain the "BAS (equal weights)" arm)
# Output: FigS5_EqualWeights.png / .pdf
# ==============================================================================
suppressMessages({library(dplyr); library(ggplot2); library(patchwork)})
cols <- c("BAS" = "#0072B2", "BAS (equal weights)" = "#56B4E9", "Random" = "#D55E00")
lts  <- c("BAS" = "solid", "BAS (equal weights)" = "dashed", "Random" = "solid")
panel <- function(df, xlim, tag, sub) {
  s <- df %>% filter(Method %in% names(cols), n_samples >= xlim[1], n_samples <= xlim[2]) %>%
    group_by(n_samples, Method) %>% summarise(Mean = mean(Richness), SD = sd(Richness), .groups = "drop") %>%
    mutate(Method = factor(Method, names(cols)))
  ggplot(s, aes(n_samples, Mean, color = Method, fill = Method, linetype = Method)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.15, linetype = 0) + geom_line(linewidth = 1.1) +
    scale_color_manual(values = cols) + scale_fill_manual(values = cols) + scale_linetype_manual(values = lts) +
    labs(tag = tag, x = "Number of Samples", y = "Cumulative Richness", subtitle = sub) +
    theme_minimal(base_size = 14) + theme(legend.position = "bottom", legend.title = element_blank())
}
L <- readRDS("simulation_results_Fig2AB.rds"); G <- readRDS("simulation_results_Fig3AB.rds")
p <- (panel(L, c(5, 30), "A", "Local (tuned weights 0.4 / 0.2 / 0.4)") | panel(G, c(1000, 3000), "B", "Global (tuned weights 0.2 / 0.6 / 0.2)")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("FigS5_EqualWeights.png", p, width = 12, height = 5.5, dpi = 300); ggsave("FigS5_EqualWeights.pdf", p, width = 12, height = 5.5)

summ <- function(df, ns) df %>% filter(n_samples %in% ns, Method %in% names(cols)) %>% group_by(n_samples, Method) %>%
  summarise(m = mean(Richness), sd = sd(Richness), .groups = "drop") %>% tidyr::pivot_wider(names_from = Method, values_from = c(m, sd))
cat("\n--- local ---\n");  print(as.data.frame(summ(L, c(10, 20, 30))), digits = 5)
cat("\n--- global ---\n"); print(as.data.frame(summ(G, c(1200, 1500, 2000, 2500, 3000))), digits = 6)
pt <- function(df, n) { w <- df %>% filter(n_samples == n, Method %in% c("BAS", "BAS (equal weights)")) %>% tidyr::pivot_wider(names_from = Method, values_from = Richness)
  t.test(w$BAS, w[["BAS (equal weights)"]], paired = TRUE)$p.value }
cat("\nP (tuned vs equal, two-sided paired): local n=20", signif(pt(L, 20), 3), " n=30", signif(pt(L, 30), 3), " global n=1500", signif(pt(G, 1500), 3), " n=3000", signif(pt(G, 3000), 3), "\n")
