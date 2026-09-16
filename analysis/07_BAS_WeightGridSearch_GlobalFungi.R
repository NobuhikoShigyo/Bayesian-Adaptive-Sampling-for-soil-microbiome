# ==============================================================================
# BAS Weight Grid Search — GlobalFungi Dataset (v2: Fig3準拠)
#
# v1 からの主な変更点:
#   1. LCBD計算: comm_pa + 毎バッチdecostand → 事前計算済みcomm_hel を使用
#      (Fig3と同一: 全49245種で距離計算, 種フィルタなし)
#   2. ranger に quantreg = TRUE を追加 (Fig3と同一)
#   3. N_TREES: 100 → 200 (Fig3と同一)
#   4. ワーカーが comm_hel + comm_pa の両方を独立ロード
#
# メモリ注意:
#   comm_hel: ~2GB/ワーカー, comm_pa(logical): ~250MB/ワーカー
#   4コア合計 ≈ 9GB (ワーカー分) + 親プロセス → 32GB RAM 環境では num_cores=3 推奨
#
# Input:
#   master_data_GlobalFungi.rds
#   comm_pa_GlobalFungi.rds   (5000 × 49245, logical)
#   comm_hel_GlobalFungi.rds  (5000 × 49245, Hellinger変換済み double)
#
# Output:
#   gridsearch_results_GlobalFungi_light.rds
#   BAS_WeightOptimization_GlobalFungi_light.png
# ==============================================================================

library(dplyr)
library(vegan)
library(ranger)
library(foreach)
library(doParallel)
library(grid)
library(gridExtra)

# setwd("<path>/BAS_publication/analysis")   # run from the analysis/ directory

# ==============================================================================
# ★ 調整可能パラメータ
# ==============================================================================

N_SIMULATIONS  <- 20      # シミュレーション回数  (精度↑→100, 速度優先→20)
N_TREES        <- 50      # ranger 木の本数       (light版: 50)
N_INITIAL      <- 1000    # 初期サンプル数
N_TOTAL        <- 3000    # 合計サンプル数
BATCH_SIZE     <- 100     # バッチサイズ

eval_ns <- c(1200, 1500, 2000, 2500, 3000)

out_rds <- "gridsearch_results_GlobalFungi_light.rds"
out_png <- "BAS_WeightOptimization_GlobalFungi_light.png"

# ==============================================================================
# 内部定数
# ==============================================================================

N_BATCHES <- (N_TOTAL - N_INITIAL) / BATCH_SIZE   # = 20
DATA_PATH <- getwd()

# ==============================================================================
# Step 1: master_data 読み込み・検証
# ==============================================================================

master_data <- readRDS(file.path(DATA_PATH, "master_data_GlobalFungi.rds"))

req_cols <- c("pH", "MAT", "MAP", "SOC", "X", "Y", "TrueRichness",
              "latitude", "sample_ID")
missing_cols <- setdiff(req_cols, names(master_data))
if (length(missing_cols) > 0)
  stop("master_data に必要な列がありません: ", paste(missing_cols, collapse = ", "))

na_counts <- sapply(master_data[, req_cols], function(x) sum(is.na(x)))
if (any(na_counts > 0)) {
  message("NA を含む行を除外します:")
  print(na_counts[na_counts > 0])
  master_data <- master_data %>%
    filter(complete.cases(dplyr::select(., all_of(req_cols))))
}

message("Dataset ready: ", nrow(master_data), " sites")

# ==============================================================================
# Step 2: Weight grid (66 combinations)
# ==============================================================================

weight_grid <- expand.grid(
  Wr   = seq(0, 1, by = 0.1),
  Wu   = seq(0, 1, by = 0.1),
  Wunc = seq(0, 1, by = 0.1)
) %>%
  dplyr::filter(round(Wr + Wu + Wunc, 10) == 1.0)

stopifnot(nrow(weight_grid) == 66)

# ==============================================================================
# Step 3: コア数・見積もり
# ==============================================================================

# メモリ: comm_hel(double ~2GB) + comm_pa(logical ~250MB) = ~2.25GB/ワーカー
# 32GB RAM 環境では 3コア推奨（余裕を持って）
num_cores <- min(max(1L, parallel::detectCores() - 1L), 2L)

cat("\n=== 実行設定 (v2: Fig3準拠) ===\n")
cat(sprintf("  使用コア数 : %d\n", num_cores))
cat(sprintf("  N_SIM=%d  N_TREES=%d  (Fig3準拠)\n", N_SIMULATIONS, N_TREES))
cat(sprintf("  LCBD: comm_hel (事前計算済み, 全種)\n"))
cat(sprintf("  quantreg: TRUE  (Fig3準拠)\n"))
cat(sprintf("  推定ワーカーRAM: %.1f GB/コア × %d = %.1f GB\n",
            2.25, num_cores, 2.25 * num_cores))
cat("================================\n\n")

# ==============================================================================
# Step 4: Grid search simulation
# ==============================================================================

if (file.exists(out_rds)) {
  message(out_rds, " が存在するためスキップ。読み込み中...")
  results_grid <- readRDS(out_rds)
  message("読み込み完了: ", nrow(results_grid), " 行")

} else {

  message("クラスター起動: ", num_cores, " コア")
  cl <- makeCluster(num_cores, type = "PSOCK", outfile = "")

  # 各ワーカーが独立して comm_pa + comm_hel をディスクから読み込む
  clusterExport(cl, varlist = "DATA_PATH", envir = environment())
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(vegan)
      library(ranger)
    })
    # comm_pa: 存在・不在（logical, ~250MB）← 種蓄積の計算に使用
    comm_pa_w  <- readRDS(file.path(DATA_PATH, "comm_pa_GlobalFungi.rds"))
    comm_pa_w  <- comm_pa_w > 0L   # double → logical

    # comm_hel: Hellinger変換済み（double, ~2GB）← LCBD計算に使用（Fig3準拠）
    comm_hel_w <- readRDS(file.path(DATA_PATH, "comm_hel_GlobalFungi.rds"))

    master_data_w <- readRDS(file.path(DATA_PATH, "master_data_GlobalFungi.rds"))
    invisible(NULL)
  })
  message("ワーカー初期化完了")

  registerDoParallel(cl)
  message("グリッドサーチ開始 ...")
  start_time <- Sys.time()

  results_grid <- foreach(
    w         = seq_len(nrow(weight_grid)),
    .combine  = rbind,
    .packages = c("vegan", "dplyr", "ranger"),
    .export   = c("weight_grid",
                  "N_INITIAL", "N_TOTAL", "BATCH_SIZE",
                  "N_BATCHES", "N_SIMULATIONS", "N_TREES",
                  "start_time")
  ) %dopar% {

    Wr   <- weight_grid$Wr[w]
    Wu   <- weight_grid$Wu[w]
    Wunc <- weight_grid$Wunc[w]

    norm01 <- function(x) {
      x[is.na(x)] <- 0
      rng <- range(x, na.rm = TRUE)
      if (diff(rng) == 0) return(rep(0, length(x)))
      (x - rng[1]) / diff(rng)
    }

    sim_rows <- vector("list", N_SIMULATIONS)

    for (sim in seq_len(N_SIMULATIONS)) {
      if (sim %% 5 == 0)
        message(sprintf("  [Weight %d/%d] sim %d/%d ...",
                        w, nrow(weight_grid), sim, N_SIMULATIONS))
      set.seed(123 + sim)

      prob_w   <- ifelse(master_data_w$latitude > 20, 1.0, 0.05)
      curr_sam <- sample(seq_len(nrow(master_data_w)), N_INITIAL, prob = prob_w)
      curr_uns <- setdiff(seq_len(nrow(master_data_w)), curr_sam)

      # 初期豊富度（comm_pa で計算）
      init_spc <- colSums(comm_pa_w[curr_sam, , drop = FALSE])
      init_r   <- sum(init_spc > 0)

      curve_rows <- vector("list", N_BATCHES + 1L)
      curve_rows[[1L]] <- data.frame(
        Wr = Wr, Wu = Wu, Wunc = Wunc, sim = sim,
        n_samples = N_INITIAL, Richness = init_r
      )

      curr_spc <- init_spc

      for (b in seq_len(N_BATCHES)) {
        if (length(curr_uns) < BATCH_SIZE) break
        if (b == 1 && sim == 1)
          message(sprintf("  [Weight %d] 最初のバッチ完了", w))

        # ---- LCBD計算: Fig3準拠 ----
        # comm_hel（事前Hellinger変換済み）をそのまま使用
        # 種フィルタなし（全49245種で距離計算）← v1 との最大の違い
        # centered行列を作らずHuygens分解で計算（メモリ削減）
        cm       <- comm_hel_w[curr_sam, , drop = FALSE]
        ss_rows  <- rowSums(cm^2)
        col_mean <- colMeans(cm)
        ss_total <- sum(ss_rows) - nrow(cm) * sum(col_mean^2)
        curr_lcbd <- (ss_rows - 2 * as.vector(cm %*% col_mean) + sum(col_mean^2)) /
                     max(ss_total, 1e-12)
        rm(cm, ss_rows, col_mean)

        train_d      <- master_data_w[curr_sam, ]
        train_d$LCBD <- curr_lcbd
        cand_d       <- master_data_w[curr_uns, ]

        # ---- ranger: quantreg = FALSE (light版) ----
        mod_r <- ranger::ranger(
          TrueRichness ~ pH + MAT + MAP + SOC + X + Y,
          data        = train_d,
          num.trees   = N_TREES,
          quantreg    = FALSE,   # light版
          keep.inbag  = TRUE,    # Jackknife SE に必要
          num.threads = 1,
          save.memory = TRUE
        )
        mod_l <- ranger::ranger(
          LCBD ~ pH + MAT + MAP + SOC + X + Y,
          data        = train_d,
          num.trees   = N_TREES,
          quantreg    = FALSE,   # light版
          keep.inbag  = TRUE,
          num.threads = 1,
          save.memory = TRUE
        )

        pred_r_val <- predict(mod_r, data = cand_d)$predictions
        pred_l_val <- predict(mod_l, data = cand_d)$predictions
        pred_r_se  <- predict(mod_r, data = cand_d, type = "se")$se

        sc       <- Wr   * norm01(pred_r_val) +
                    Wu   * norm01(pred_l_val) +
                    Wunc * norm01(pred_r_se)
        best_idx <- order(sc, decreasing = TRUE)[seq_len(BATCH_SIZE)]

        best_global <- curr_uns[best_idx]
        curr_sam    <- c(curr_sam, best_global)
        curr_uns    <- curr_uns[-best_idx]

        new_spc  <- colSums(comm_pa_w[best_global, , drop = FALSE])
        curr_spc <- curr_spc + new_spc
        cur_r    <- sum(curr_spc > 0)

        curve_rows[[b + 1L]] <- data.frame(
          Wr = Wr, Wu = Wu, Wunc = Wunc, sim = sim,
          n_samples = N_INITIAL + b * BATCH_SIZE,
          Richness  = cur_r
        )


      } # end batch loop

      gc()
      sim_rows[[sim]] <- do.call(rbind, curve_rows)
    } # end sim loop

    elapsed   <- difftime(Sys.time(), start_time, units = "mins")
    rate      <- as.numeric(elapsed) / w
    remaining <- rate * (nrow(weight_grid) - w)
    message(sprintf(
      "[%s] Weight %d/%d done (Wr=%.1f Wu=%.1f Wunc=%.1f) | 経過: %.1f分 | 残り推定: %.1f分 (%.1f時間) | 完了予定: %s",
      format(Sys.time(), "%H:%M:%S"),
      w, nrow(weight_grid),
      Wr, Wu, Wunc,
      as.numeric(elapsed),
      remaining,
      remaining / 60,
      format(Sys.time() + as.difftime(remaining, units = "mins"), "%m/%d %H:%M")
    ))

    do.call(rbind, sim_rows)
  } # end foreach

  stopCluster(cl)
  message("グリッドサーチ完了: ", nrow(results_grid), " 行")

  saveRDS(results_grid, out_rds)
  message("保存: ", out_rds)
}

# ==============================================================================
# Step 5: 集計
# ==============================================================================

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

message("summary_by_n: ", nrow(summary_by_n), " 行")

# ==============================================================================
# Step 6: 最適重み
# ==============================================================================

best_by_n <- summary_by_n %>%
  dplyr::group_by(eval_n) %>%
  dplyr::slice_max(mean_richness, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

current_by_n <- summary_by_n %>%
  dplyr::filter(
    round(Wr, 1) == round(1 / 3, 1),
    round(Wu, 1) == round(1 / 3, 1)
  )

trajectory <- best_by_n %>%
  dplyr::arrange(eval_n) %>%
  dplyr::select(eval_n, Wr, Wu, Wunc, mean_richness)

# ==============================================================================
# Step 7: コンソール出力
# ==============================================================================

cat("\n========================================================\n")
cat("  Top-3 optimal weight combinations per eval_n  (v2: Fig3準拠)\n")
cat("========================================================\n")

for (en in levels(summary_by_n$eval_n)) {
  cat("\n--- eval_n =", en, "---\n")
  top3 <- summary_by_n %>%
    dplyr::filter(eval_n == en) %>%
    dplyr::slice_max(mean_richness, n = 3, with_ties = FALSE) %>%
    dplyr::select(Wr, Wu, Wunc, mean_richness, sd_richness)
  print(top3, digits = 4)
}

cat("\n--- Optimal weight trajectory ---\n")
print(trajectory %>% dplyr::select(eval_n, Wr, Wu, Wunc, mean_richness), digits = 4)

# ==============================================================================
# Step 8: ggtern プロット
# ==============================================================================

if (!requireNamespace("ggtern", quietly = TRUE)) install.packages("ggtern")
library(ggtern)

rich_range <- range(summary_by_n$mean_richness)

p_facet <- ggtern(
  summary_by_n,
  aes(x = Wr, y = Wu, z = Wunc, color = mean_richness)
) +
  geom_point(size = 2.8, alpha = 0.88) +
  geom_point(
    data        = best_by_n,
    color       = "red", fill = "red",
    shape       = 23, size = 5, show.legend = FALSE
  ) +
  geom_point(
    data        = current_by_n,
    color       = "white", fill = "white",
    shape       = 21, size = 3.5, show.legend = FALSE
  ) +
  scale_color_viridis_c(option = "plasma", name = "Mean\nrichness", limits = rich_range) +
  facet_wrap(~ eval_n, ncol = 3) +
  labs(
    title    = "BAS weight optimization — GlobalFungi (v2: Fig3準拠)",
    subtitle = paste0(
      "\u25c6 red = optimal  |  \u25cb white = equal-weight approx. (0.3, 0.3, 0.4)",
      "\n", N_SIMULATIONS, " simulations; ranger RF (", N_TREES,
      " trees, quantreg=TRUE, jackknife SE)"
    ),
    x = expression(W[rich]),
    y = expression(W[uniq]),
    z = expression(W[unc])
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 9, color = "grey40"),
    legend.position = "right",
    strip.text      = element_text(face = "bold", size = 11)
  )

p_traj <- ggtern(trajectory, aes(x = Wr, y = Wu, z = Wunc)) +
  geom_path(
    aes(group = 1),
    arrow     = arrow(length = unit(0.3, "cm"), ends = "last", type = "closed"),
    color     = "grey30", linewidth = 0.9, linejoin = "round"
  ) +
  geom_point(
    aes(fill = eval_n), shape = 21, size = 6, color = "black"
  ) +
  scale_fill_brewer(palette = "RdYlBu", name = "Eval n") +
  labs(
    title    = paste0("Trajectory of optimal weights  (",
                      paste(paste0("n=", eval_ns), collapse = " \u2192 "), ")"),
    subtitle = paste0("GlobalFungi (v2: Fig3準拠)  |  ",
                      N_SIMULATIONS, " sims \u00d7 66 weight combos"),
    x = expression(W[rich]),
    y = expression(W[uniq]),
    z = expression(W[unc])
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 9, color = "grey40"),
    legend.position = "right"
  )

# ==============================================================================
# Step 9: PNG 保存
# ==============================================================================

png(out_png, width = 14, height = 12, units = "in", res = 300)
grid.arrange(p_facet, p_traj, nrow = 2, heights = c(2.5, 1.5))
dev.off()

message("保存: ", out_png)
message("\n=== All done (v2) ===")
