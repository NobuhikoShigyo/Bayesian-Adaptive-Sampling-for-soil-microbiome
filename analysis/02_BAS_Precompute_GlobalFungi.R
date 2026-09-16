# ==============================================================================
# BAS 事前計算スクリプト（一回だけ実行）
#
# 生成ファイル:
#   K_GlobalFungi.rds       … tcrossprod(comm_hel) の 5000×5000 カーネル行列
#                              (~200MB) ← comm_hel の代わりに使用
#   comm_pa_sp_GlobalFungi.rds … sparse logical 行列 (~数MB〜数十MB)
#
# 実行後、BAS_WeightGridSearch_GlobalFungi_v3.R でこれらを使用。
# ==============================================================================

library(Matrix)

DATA_PATH <- getwd()   # run from the analysis/ directory

# ------------------------------------------------------------------------------
# 1. カーネル行列 K = comm_hel %*% t(comm_hel)  (5000×5000)
# ------------------------------------------------------------------------------

message("comm_hel 読み込み中 (~2GB)...")
comm_hel <- readRDS(file.path(DATA_PATH, "comm_hel_GlobalFungi.rds"))
message("行列サイズ: ", nrow(comm_hel), " × ", ncol(comm_hel))

message("K = tcrossprod(comm_hel) 計算中 ...")
t0 <- proc.time()
K <- tcrossprod(comm_hel)   # 5000×5000, BLAS を使用するため通常 <2分
elapsed <- proc.time() - t0
message(sprintf("  完了: %.1f秒", elapsed["elapsed"]))
message(sprintf("  K サイズ: %.1f MB", object.size(K) / 1e6))

k_path <- file.path(DATA_PATH, "K_GlobalFungi.rds")
saveRDS(K, k_path, compress = FALSE)   # compress=FALSE で高速保存
message("保存: ", k_path)

rm(comm_hel, K); gc()

# ------------------------------------------------------------------------------
# 2. sparse comm_pa
# ------------------------------------------------------------------------------

message("comm_pa 読み込み中 (~250MB)...")
comm_pa <- readRDS(file.path(DATA_PATH, "comm_pa_GlobalFungi.rds"))
comm_pa <- comm_pa > 0L

fill_rate <- mean(comm_pa)
message(sprintf("  fill rate: %.4f%% (%s 非ゼロ要素)",
                fill_rate * 100,
                format(sum(comm_pa), big.mark = ",")))

message("sparse 変換中...")
comm_pa_sp <- Matrix::Matrix(comm_pa, sparse = TRUE)
message(sprintf("  sparse サイズ: %.1f MB (dense: %.1f MB)",
                object.size(comm_pa_sp) / 1e6,
                object.size(comm_pa) / 1e6))

sp_path <- file.path(DATA_PATH, "comm_pa_sp_GlobalFungi.rds")
saveRDS(comm_pa_sp, sp_path, compress = FALSE)
message("保存: ", sp_path)

rm(comm_pa, comm_pa_sp); gc()

message("\n=== 事前計算完了 ===")
message("次に BAS_WeightGridSearch_GlobalFungi_v3.R を実行してください。")
