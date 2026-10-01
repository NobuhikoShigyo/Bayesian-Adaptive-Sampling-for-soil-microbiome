# ==============================================================================
# 06: Ternary plots of the weight grid searches (Fig. S1 local, Fig. S2 global)
#
# Input : gridsearch_results_local.rds (05), gridsearch_results_global.rds (07)
#         columns Wr, Wu, Wunc, sim, n_samples, Richness
# Output: FigS1_WeightGridSearch_local.png, FigS2_WeightGridSearch_global.png
# NOTE: ggtern is loaded after all ggplot2 code to avoid namespace conflicts.
# ==============================================================================
suppressMessages({library(dplyr); library(grid); library(gridExtra)})
if (!requireNamespace("ggtern", quietly = TRUE)) install.packages("ggtern")
library(ggtern)

plot_grid_ternary <- function(rds, eval_ns, used_w, title, subtitle, file, ncol = 3) {
  res <- readRDS(rds)
  summ <- res %>% filter(n_samples %in% eval_ns) %>% group_by(eval_n = n_samples, Wr, Wu, Wunc) %>%
    summarise(mean_richness = mean(Richness), .groups = "drop") %>%
    mutate(eval_n = factor(eval_n, levels = eval_ns, labels = paste0("n = ", format(eval_ns, big.mark = ","))))
  best <- summ %>% group_by(eval_n) %>% slice_max(mean_richness, n = 1, with_ties = FALSE) %>% ungroup()
  used <- summ %>% filter(abs(Wr - used_w[1]) < 1e-6, abs(Wu - used_w[2]) < 1e-6, abs(Wunc - used_w[3]) < 1e-6)
  cat("\n", title, "\n"); for (lv in levels(summ$eval_n)) { cat("--", lv, "top 3 --\n")
    print(as.data.frame(summ %>% filter(eval_n == lv) %>% slice_max(mean_richness, n = 3, with_ties = FALSE) %>% select(Wr, Wu, Wunc, mean_richness)), digits = 5) }
  cat("weights used in the main analysis:\n"); print(as.data.frame(used %>% select(eval_n, Wr, Wu, Wunc, mean_richness)), digits = 5)
  p_facet <- ggtern(summ, aes(x = Wr, y = Wu, z = Wunc, color = mean_richness)) +
    geom_point(size = 2.8, alpha = 0.88) +
    geom_point(data = best, color = "red", fill = "red", shape = 23, size = 5, show.legend = FALSE) +
    geom_point(data = used, color = "black", fill = "white", shape = 21, size = 4, stroke = 1, show.legend = FALSE) +
    scale_color_viridis_c(option = "plasma", name = "Mean\nrichness", limits = range(summ$mean_richness)) +
    facet_wrap(~eval_n, ncol = ncol) +
    labs(title = title, subtitle = subtitle, x = expression(W[rich]), y = expression(W[uniq]), z = expression(W[unc])) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 9, color = "grey40"),
          legend.position = "right", strip.text = element_text(face = "bold", size = 11))
  traj <- best %>% arrange(eval_n)
  p_traj <- ggtern(traj, aes(x = Wr, y = Wu, z = Wunc)) +
    geom_path(aes(group = 1), arrow = arrow(length = unit(0.3, "cm"), ends = "last", type = "closed"), color = "grey30", linewidth = 0.9) +
    geom_point(aes(fill = eval_n), shape = 21, size = 6, color = "black") +
    scale_fill_brewer(palette = "RdYlBu", name = "Eval n") +
    labs(title = "Trajectory of optimal weights", x = expression(W[rich]), y = expression(W[uniq]), z = expression(W[unc])) +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 13))
  png(file, width = 14, height = 12, units = "in", res = 300)
  grid.arrange(p_facet, p_traj, ncol = 1, heights = c(1.6, 1))
  dev.off()
  message("Saved: ", file)
}

W_LOCAL  <- as.numeric(strsplit(Sys.getenv("BAS_WEIGHTS_LOCAL",  "0.4,0.3,0.3"), ",")[[1]])
W_GLOBAL <- as.numeric(strsplit(Sys.getenv("BAS_WEIGHTS_GLOBAL", "0.2,0.4,0.4"), ",")[[1]])

plot_grid_ternary("gridsearch_results_local.rds", c(10, 15, 20, 25, 30), W_LOCAL,
  "BAS weight grid search — local catchment (Case 1)",
  "◆ red = best combination  |  ○ white = weights used in the main analysis\n66 combinations x 50 simulations; regression kriging; uniqueness = novelty",
  "FigS1_WeightGridSearch_local.png")
plot_grid_ternary("gridsearch_results_global.rds", c(1200, 1500, 2000, 2500, 3000), W_GLOBAL,
  "BAS weight grid search — GlobalFungi (Case 2)",
  "◆ red = best combination  |  ○ white = weights used in the main analysis\n66 combinations x 10 simulations; QRF (100 trees); uniqueness = novelty",
  "FigS2_WeightGridSearch_global.png")
