# ==============================================================================
# BAS (Balanced) vs Random vs Oracle: Simulation for Figure 3A-B  (v2: Windows版)
#
# v1 からの変更点:
#   1. comm_hel (2GB/ワーカー) → K カーネル行列 (200MB, 親で1回) に置換
#   2. comm_pa (dense) → comm_pa_sp (sparse) に置換
#   3. 並列方式: Windows PSOCK 固定 (fork 非対応)
#   4. predict(mod_r) を1回に統合 (2回 → 1回)
#   5. Oracle 戦略 (greedy 貪欲法) を追加
#
# 前提: BAS_Precompute.R を先に実行して以下を生成済みであること
#   K_GlobalFungi.rds          (5000×5000, ~200MB)
#   comm_pa_sp_GlobalFungi.rds (sparse, ~数MB〜数十MB)
#
# メモリ見積もり (Windows, 16コア / 25GB RAM):
#   PSOCK方式: 各ワーカーがデータのコピーを保持
#   K_w    : 200MB × 14ワーカー = ~2.8GB
#   comm_pa: ~50MB × 14ワーカー = ~0.7GB
#   ranger : ~30MB × 2モデル × 14ワーカー = ~0.8GB
#   親プロセス + OS         = ~1.0GB
#   合計ピーク              : ~6GB / 25GB  ✓
#
# 所要時間見積もり (Windows 16コア, 14並列):
#   ranger(6特徴, 200木, quantreg=TRUE): ~1-3秒/モデル
#   BAS 1sim = 2モデル × 20バッチ = ~40-120秒
#   Oracle 1sim = sparse行列演算 × 20バッチ = ~数秒 (無視できる)
#   100sim / 14コア = 8ラウンド
#   PSCOKデータ転送 (K=200MB × 14): ~2-3分
#   合計: 約 15-35 分
#
# 科学的設定 (Fig3準拠, 変更不可):
#   N_TREES=200, quantreg=TRUE, N_SIMULATIONS=100
# ==============================================================================

# ==============================================================================
# Step 1: Build dataset from GlobalFungi — apply WorldClim, subsample 5000 sites
# (v1 と同一。master_data / comm_pa / comm_hel の .rds が既にある場合はスキップ)
# ==============================================================================

if (!all(file.exists(c("master_data_GlobalFungi.rds",
                       "comm_pa_GlobalFungi.rds",
                       "comm_hel_GlobalFungi.rds")))) {

  library(dplyr)
  library(data.table)
  library(geodata)
  library(terra)
  library(vegan)
  library(sf)

  coords <- fread("GlobalFungi_5_sample_metadata.txt", quote = "", nThread = 4)

  soil_data <- coords %>%
    mutate(
      latitude  = as.numeric(latitude),
      longitude = as.numeric(longitude),
      pH        = as.numeric(pH)
    ) %>%
    filter(!is.na(latitude) & !is.na(longitude) & !is.na(pH)) %>%
    filter(sample_type %in% c("soil", "topsoil"))

  clim_data   <- worldclim_global(var = "bio", res = 10, path = tempdir())
  coords_ll   <- soil_data %>% dplyr::select(longitude, latitude)
  clim_vals   <- terra::extract(clim_data, coords_ll)
  soil_data$MAT_wc <- clim_vals$wc2.1_10m_bio_1
  soil_data$MAP_wc <- clim_vals$wc2.1_10m_bio_12
  soil_data_valid  <- soil_data %>% filter(!is.na(MAT_wc) & !is.na(MAP_wc))

  set.seed(123)
  soil_data_5000 <- if (nrow(soil_data_valid) > 5000) {
    soil_data_valid %>% sample_n(5000)
  } else {
    soil_data_valid
  }
  target_ids <- as.character(soil_data_5000$sample_ID)

  message("Loading community data (~12 GB)...")
  comm    <- fread("GlobalFungi_5_SH_abundance_ITS1_ITS2.txt", nThread = 4)
  comm_filtered <- comm[sample_ID %in% target_ids]
  rm(comm, coords, soil_data, soil_data_valid, clim_data); gc()

  comm_mat <- as.matrix(comm_filtered[, -1, with = FALSE])
  rownames(comm_mat) <- comm_filtered$sample_ID
  rm(comm_filtered); gc()

  common_ids <- intersect(rownames(comm_mat), soil_data_5000$sample_ID)
  site_data  <- soil_data_5000 %>%
    filter(sample_ID %in% common_ids) %>%
    arrange(match(sample_ID, common_ids))
  comm_data  <- comm_mat[common_ids, ]
  comm_data  <- comm_data[, colSums(comm_data) > 0]
  comm_pa    <- (comm_data > 0) * 1L

  richness_df <- data.frame(SiteID = rownames(comm_pa), TrueRichness = rowSums(comm_pa))
  comm_hel    <- vegan::decostand(comm_data, "hellinger")

  sites_sf    <- st_as_sf(site_data, coords = c("longitude", "latitude"), crs = 4326)
  sites_proj  <- st_transform(sites_sf, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  coords_moll <- st_coordinates(sites_proj)

  master_data <- site_data %>%
    mutate(
      X         = coords_moll[, 1],
      Y         = coords_moll[, 2],
      pH        = as.numeric(pH),
      MAT       = as.numeric(MAT_wc),
      MAP       = as.numeric(MAP_wc),
      SOC       = as.numeric(SOC),
      latitude  = as.numeric(latitude),
      longitude = as.numeric(longitude)
    ) %>%
    left_join(richness_df, by = c("sample_ID" = "SiteID")) %>%
    as.data.frame()

  saveRDS(master_data, "master_data_GlobalFungi.rds")   # NB: the copy in the repository is trimmed to the 12 columns used downstream
  saveRDS(comm_pa,     "comm_pa_GlobalFungi.rds")
  saveRDS(comm_hel,    "comm_hel_GlobalFungi.rds")
  message("Saved: master_data / comm_pa / comm_hel")

} else {
  message("既存の .rds ファイルを使用します (Step 1 スキップ)")
}

# ==============================================================================
# Step 2: BAS Simulation (v2: RAM最適化)
# ==============================================================================

library(dplyr)
library(vegan)
library(ggplot2)
library(patchwork)
library(Matrix)
library(ranger)
library(foreach)

# setwd("<path>/BAS_publication/analysis")   # run from the analysis/ directory

# --- スレッド制御（ranger は num.threads=1 で個別制御）---
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}

# --- シミュレーション設定 (Fig3準拠) ---
N_SIMULATIONS <- 100
N_INITIAL     <- 1000
N_TOTAL       <- 3000
BATCH_SIZE    <- 100
N_BATCHES     <- floor((N_TOTAL - N_INITIAL) / BATCH_SIZE)
N_TREES       <- 200      # Fig3準拠 (変更不可)

strategies <- list(
  BAS_Balanced = list(Wr = 0.25, Wu = 0.25, Wunc = 0.5, Name = "BAS",    Type = "Env"),
  Random       = list(                                   Name = "Random", Type = "Random"),
  Oracle       = list(                                   Name = "Oracle", Type = "Oracle")
)

norm01 <- function(x) {
  if (length(x[!is.na(x)]) == 0 || max(x, na.rm = TRUE) == min(x, na.rm = TRUE))
    return(rep(0, length(x)))
  out <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  out[is.na(out)] <- 0
  out
}

# LCBD カーネル計算（comm_hel 不要, K のみ使用）
# v1: rowSums(centered^2) / ss_total
# v2: 数学的等価 — K[curr_sam, curr_sam] から計算
calc_lcbd_kernel <- function(K, idx) {
  K_sub  <- K[idx, idx, drop = FALSE]
  n      <- length(idx)
  rs     <- diag(K_sub)
  rs_K   <- rowSums(K_sub)
  tot_K  <- sum(K_sub)
  dev_sq <- rs - 2 * rs_K / n + tot_K / n^2   # ||x_i - centroid||^2
  ss_total <- sum(dev_sq)
  if (ss_total > 0) dev_sq / ss_total else rep(0, n)
}

# ==============================================================================
# 親プロセスでデータを一度だけロード（fork で全ワーカーが共有）
# ==============================================================================

DATA_PATH <- getwd()

# K と sparse comm_pa の存在確認
k_path  <- file.path(DATA_PATH, "K_GlobalFungi.rds")
sp_path <- file.path(DATA_PATH, "comm_pa_sp_GlobalFungi.rds")

if (!file.exists(k_path) || !file.exists(sp_path)) {
  stop(
    "事前計算ファイルが見つかりません。先に BAS_Precompute.R を実行してください。\n",
    "  必要: ", k_path, "\n",
    "  必要: ", sp_path
  )
}

message("K 行列読み込み中 (~200MB)...")
K_w <- readRDS(k_path)

message("sparse comm_pa 読み込み中...")
comm_pa_w <- readRDS(sp_path)

message("master_data 読み込み中...")
master_data <- readRDS(file.path(DATA_PATH, "master_data_GlobalFungi.rds"))

message(sprintf("  K: %d×%d (%.0f MB)", nrow(K_w), ncol(K_w), object.size(K_w) / 1e6))
message(sprintf("  comm_pa_sp: %d×%d (%.0f MB)", nrow(comm_pa_w), ncol(comm_pa_w), object.size(comm_pa_w) / 1e6))
message(sprintf("  master_data: %d 行", nrow(master_data)))

# ==============================================================================
# 並列バックエンド (Windows: PSOCK固定, 16コア中14コア使用)
# ==============================================================================

# 100 sims / 14コア = 8 ラウンド
num_cores <- 14L

library(doParallel)
cat("\n=== 実行設定 (Fig3AB v2: Windows 16コア / 25GB RAM) ===\n")
cat(sprintf("  使用コア数    : %d / 16\n", num_cores))
cat(sprintf("  並列方式      : PSOCK (Windows固定)\n"))
cat(sprintf("  N_SIM=%d  N_TREES=%d  quantreg=TRUE\n", N_SIMULATIONS, N_TREES))
cat(sprintf("  LCBD          : K カーネル行列 (数学的等価)\n"))
cat(sprintf("  メモリ推定    : K=%.0fMB × %dワーカー = %.0fMB\n",
            object.size(K_w) / 1e6, num_cores,
            object.size(K_w) / 1e6 * num_cores))
cat(sprintf("  ピーク推定    : ~6GB / 25GB\n"))
cat(sprintf("  所要時間目安  : 15-35分\n"))
cat("=======================================================\n\n")

cl <- makeCluster(num_cores, type = "PSOCK", outfile = "")
clusterExport(cl, varlist = c("K_w", "comm_pa_w", "master_data",
                              "strategies", "norm01", "calc_lcbd_kernel",
                              "N_INITIAL", "N_TOTAL", "BATCH_SIZE",
                              "N_BATCHES", "N_TREES"),
              envir = environment())
clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(ranger); library(Matrix) })
  invisible(NULL)
})
registerDoParallel(cl)
message("PSOCK ワーカー初期化完了 (", num_cores, " workers)")

message("Starting simulation: N_INITIAL=", N_INITIAL, " N_TOTAL=", N_TOTAL,
        " BATCH=", BATCH_SIZE, " (Fig3準拠)")

start_time <- Sys.time()

# ==============================================================================
# メインループ（sim 単位で並列）
# ==============================================================================

results_df <- foreach(
  sim            = 1:N_SIMULATIONS,
  .combine       = rbind,
  .errorhandling = "remove",
  .packages      = c("ranger", "Matrix")
) %dopar% {

  set.seed(123 + sim)

  prob_weights <- ifelse(master_data$latitude > 20, 1.0, 0.05)
  initial_idx  <- sample(seq_len(nrow(master_data)), N_INITIAL, prob = prob_weights)
  remain_idx   <- setdiff(seq_len(nrow(master_data)), initial_idx)
  random_order <- sample(remain_idx)

  # sparse colSums で初期豊富度
  init_species_count <- Matrix::colSums(comm_pa_w[initial_idx, , drop = FALSE])
  init_r             <- sum(init_species_count > 0)

  sim_res <- list()

  for (strat_key in names(strategies)) {
    conf <- strategies[[strat_key]]

    curr_sam           <- initial_idx
    curr_uns           <- remain_idx
    curr_species_count <- init_species_count

    n_vec    <- integer(N_BATCHES + 1)
    r_vec    <- numeric(N_BATCHES + 1)
    n_vec[1] <- N_INITIAL
    r_vec[1] <- init_r

    for (b in seq_len(N_BATCHES)) {

      next_nodes <- NULL

      if (conf$Type == "Random") {
        idx_start <- (b - 1) * BATCH_SIZE + 1
        idx_end   <- min(b * BATCH_SIZE, length(random_order))
        if (length(random_order) >= idx_start)
          next_nodes <- random_order[idx_start:idx_end]

      } else if (conf$Type == "Oracle") {
        if (length(curr_uns) >= 1) {
          already_found <- Matrix::colSums(comm_pa_w[curr_sam, , drop = FALSE]) > 0
          new_sp        <- Matrix::rowSums(comm_pa_w[curr_uns, !already_found, drop = FALSE])
          n_pick        <- min(BATCH_SIZE, length(curr_uns))
          best_idx      <- order(new_sp, decreasing = TRUE)[seq_len(n_pick)]
          next_nodes    <- curr_uns[best_idx]
          curr_uns      <- curr_uns[-best_idx]
        }

      } else if (length(curr_uns) >= BATCH_SIZE) {

        # ---- LCBD: K カーネル行列で計算（comm_hel 不要）----
        curr_lcbd <- calc_lcbd_kernel(K_w, curr_sam)

        train_d      <- master_data[curr_sam, ]
        train_d$LCBD <- curr_lcbd
        cand_d       <- master_data[curr_uns, ]

        # ---- ranger (Fig3準拠: quantreg=TRUE, N_TREES=200) ----
        mod_r <- ranger::ranger(
          TrueRichness ~ pH + MAT + MAP + SOC + X + Y,
          data        = train_d,
          num.trees   = N_TREES,
          quantreg    = TRUE,
          keep.inbag  = TRUE,
          num.threads = 1,
          save.memory = TRUE
        )
        mod_l <- ranger::ranger(
          LCBD ~ pH + MAT + MAP + SOC + X + Y,
          data        = train_d,
          num.trees   = N_TREES,
          quantreg    = TRUE,
          keep.inbag  = TRUE,
          num.threads = 1,
          save.memory = TRUE
        )

        # ---- predict を1回に統合（v1は2回呼び出し）----
        pred_r     <- predict(mod_r, data = cand_d, type = "se")
        pred_r_val <- pred_r$predictions
        pred_r_se  <- pred_r$se
        pred_l_val <- predict(mod_l, data = cand_d)$predictions

        sc       <- conf$Wr * norm01(pred_r_val) +
                    conf$Wu * norm01(pred_l_val) +
                    conf$Wunc * norm01(pred_r_se)
        best_idx   <- order(sc, decreasing = TRUE)[seq_len(BATCH_SIZE)]
        next_nodes <- curr_uns[best_idx]
        curr_uns   <- curr_uns[-best_idx]
      }

      if (!is.null(next_nodes)) {
        curr_sam           <- c(curr_sam, next_nodes)
        new_counts         <- Matrix::colSums(comm_pa_w[next_nodes, , drop = FALSE])
        curr_species_count <- curr_species_count + new_counts
        n_vec[b + 1]       <- N_INITIAL + b * BATCH_SIZE
        r_vec[b + 1]       <- sum(curr_species_count > 0)
      }
    }

    sim_res[[length(sim_res) + 1]] <- data.frame(
      n_samples = n_vec, Richness = r_vec, SimID = sim, Method = conf$Name
    )
  }

  if (sim %% 10 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    message(sprintf("[%s] sim %d/%d 完了 | 経過: %.1f分",
                    format(Sys.time(), "%H:%M:%S"), sim, N_SIMULATIONS, elapsed))
  }

  gc()
  do.call(rbind, sim_res)
}

if (!is.null(cl)) stopCluster(cl)
message("All simulations complete.")

saveRDS(results_df, "simulation_results_Fig3AB_v2.rds")
message("保存: simulation_results_Fig3AB_v2.rds")

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
message("Saved: BAS_GlobalFungi_Plot_v2_withOracle.png")
message("\n=== All done (Fig3AB v2) ===")
