# ==============================================================================
# BAS Weight Grid Search — Multi-N Post-Processing
#
# Input:  gridsearch_results.rds  (produced by BAS_WeightGridSearch.R)
#         Columns: Wr, Wu, Wunc, sim, n_samples, Richness
#
# Output: BAS_WeightOptimization_MultiN_v2.png  (width=14, height=12)
#         Console: top-3 optimal weights per eval_n
#
# NOTE: ggtern は ggplot2 系パッケージのロード後に読み込む（名前空間競合を回避）
# ==============================================================================

library(dplyr)
library(grid)        # unit()
library(gridExtra)   # grid.arrange

# ==============================================================================
# Step 1: Load raw grid-search results
# ==============================================================================

# setwd("<path>/BAS_publication/analysis")   # run from the analysis/ directory

results_grid <- readRDS("gridsearch_results.rds")
message("Loaded gridsearch_results.rds: ", nrow(results_grid), " rows")
message("n_samples values present: ",
        paste(sort(unique(results_grid$n_samples)), collapse = ", "))

# ==============================================================================
# Step 2: Aggregate mean_richness per weight combo for each eval_n
# ==============================================================================

eval_ns <- c(10, 15, 20, 25, 30, 35, 40)

summary_by_n <- results_grid %>%
  dplyr::filter(n_samples %in% eval_ns) %>%
  dplyr::group_by(eval_n = n_samples, Wr, Wu, Wunc) %>%
  dplyr::summarise(
    mean_richness = mean(Richness),
    sd_richness   = sd(Richness),
    .groups       = "drop"
  ) %>%
  dplyr::mutate(
    eval_n = factor(eval_n,
                    levels = eval_ns,
                    labels = paste0("n = ", eval_ns))
  )

message("summary_by_n rows: ", nrow(summary_by_n),
        "  (", length(levels(summary_by_n$eval_n)),
        " eval_n \u00d7 66 weight combos)")

# ==============================================================================
# Step 3: Best row per eval_n  +  equal-weight reference  +  trajectory
# ==============================================================================

best_by_n <- summary_by_n %>%
  dplyr::group_by(eval_n) %>%
  dplyr::slice_max(mean_richness, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# Equal-weight grid approximation: round(1/3, 1) == 0.3 → (0.3, 0.3, 0.4)
current_by_n <- summary_by_n %>%
  dplyr::filter(
    round(Wr, 1) == round(1 / 3, 1),
    round(Wu, 1) == round(1 / 3, 1)
  )

# Trajectory: optimal weights ordered n=10→20→30→40
trajectory <- best_by_n %>%
  dplyr::arrange(eval_n) %>%          # factor levels already ordered
  dplyr::select(eval_n, Wr, Wu, Wunc, mean_richness)

# ==============================================================================
# Step 7: Console output — top 3 per eval_n
# ==============================================================================

cat("\n========================================================\n")
cat("  Top-3 optimal weight combinations per eval_n\n")
cat("========================================================\n")

for (en in levels(summary_by_n$eval_n)) {
  cat("\n--- eval_n =", en, "---\n")
  top3 <- summary_by_n %>%
    dplyr::filter(eval_n == en) %>%
    dplyr::slice_max(mean_richness, n = 3, with_ties = FALSE) %>%
    dplyr::select(Wr, Wu, Wunc, mean_richness, sd_richness)
  print(top3, digits = 4)
}

cat("\n--- Optimal weight trajectory (n = 10 → 15 → 20 → 25 → 30 → 35 → 40) ---\n")
print(trajectory %>% dplyr::select(eval_n, Wr, Wu, Wunc, mean_richness),
      digits = 4)

# ==============================================================================
# Step 4 & 5: ggtern plots
# Load ggtern HERE — after all ggplot2-based code — to avoid conflicts
# ==============================================================================

if (!requireNamespace("ggtern", quietly = TRUE)) install.packages("ggtern")
library(ggtern)

rich_range <- range(summary_by_n$mean_richness)

# ------------------------------------------------------------------------------
# Plot A: 2×2 faceted ternary plot (one panel per eval_n)
# Overlay の geom_point では inherit.aes = FALSE + x/y/z を明示
# ------------------------------------------------------------------------------

p_facet <- ggtern(
  summary_by_n,
  aes(x = Wr, y = Wu, z = Wunc, color = mean_richness)
) +
  # 全66点を plasma カラースケールで表示
  geom_point(size = 2.8, alpha = 0.88) +
  # 最適点：赤ダイヤモンド
  # inherit.aes=TRUE（デフォルト）でterenary coord を継承し、colorは定数上書き
  geom_point(
    data        = best_by_n,
    color       = "red",
    fill        = "red",
    shape       = 23,
    size        = 5,
    show.legend = FALSE
  ) +
  # 均等重み近似点 (0.3, 0.3, 0.4)：白丸
  geom_point(
    data        = current_by_n,
    color       = "white",
    fill        = "white",
    shape       = 21,
    size        = 3.5,
    show.legend = FALSE
  ) +
  scale_color_viridis_c(
    option = "plasma",
    name   = "Mean\nrichness",
    limits = rich_range
  ) +
  facet_wrap(~ eval_n, ncol = 2) +   # 2列×4行、最後の1セルは自動で空白
  labs(
    title    = "BAS weight optimization: cumulative richness at multiple n",
    subtitle = paste0(
      "\u25c6 red = optimal  |",
      "  \u25cb white = equal-weight approx. (0.3, 0.3, 0.4)",
      "\n100 simulations per weight combo; local dataset (N = 53)"
    ),
    x = expression(W[rich]),
    y = expression(W[uniq]),
    z = expression(W[unc])
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "right",
    strip.text    = element_text(face = "bold", size = 11)
  )

# ------------------------------------------------------------------------------
# Plot B: 最適重みの軌跡 (n=10→20→30→40) を矢印でつないだ三角図
# geom_path に arrow= を渡して方向を明示
# ------------------------------------------------------------------------------

p_traj <- ggtern(
  trajectory,
  aes(x = Wr, y = Wu, z = Wunc)
) +
  # 軌跡を矢印つき折れ線で描画
  geom_path(
    aes(group = 1),
    arrow     = arrow(length = unit(0.3, "cm"),
                      ends   = "last",
                      type   = "closed"),
    color     = "grey30",
    linewidth = 0.9,
    linejoin  = "round"
  ) +
  # 各 n での最適点を eval_n でカラー分け
  geom_point(
    aes(fill = eval_n),
    shape = 21,
    size  = 6,
    color = "black"
  ) +
  # ラベル：ggtern では geom_text が内部で PositionNudge を生成するため除去
  # → 各点の n は fill 凡例で判別する
  scale_fill_brewer(palette = "RdYlBu", name = "Eval n") +
  labs(
    title    = "Trajectory of optimal weights  (n = 10 \u2192 15 \u2192 20 \u2192 25 \u2192 30 \u2192 35 \u2192 40)",
    subtitle = "Arrows show how the optimal weight combination shifts as n increases",
    x = expression(W[rich]),
    y = expression(W[uniq]),
    z = expression(W[unc])
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "right"
  )

# ==============================================================================
# Step 6: PNG 保存 (width=14, height=12, dpi=300)
# grid.arrange で上段(facet 2×2) + 下段(軌跡) を縦積み
# ==============================================================================

out_file <- "BAS_WeightOptimization_MultiN_v2.png"

png(out_file, width = 14, height = 12, units = "in", res = 300)
grid.arrange(
  p_facet,
  p_traj,
  nrow    = 2,
  heights = c(2.5, 1.5)
)
dev.off()

message("Saved: ", out_file)
message("\n=== All done ===")
