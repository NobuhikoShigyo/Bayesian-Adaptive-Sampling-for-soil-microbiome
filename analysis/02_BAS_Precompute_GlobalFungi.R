# ==============================================================================
# 02: Precompute GlobalFungi objects (run once)
#
# Output files:
#   K_GlobalFungi.rds          5000 x 5000 kernel matrix tcrossprod(comm_hel) (~200 MB);
#                              not needed by the v3 scripts, kept for the LCBD comparison
#   comm_pa_sp_GlobalFungi.rds sparse logical presence/absence matrix (a few MB)
#
# comm_pa_sp_GlobalFungi.rds is read by 03, 04 and 07.
# ==============================================================================

library(Matrix)

DATA_PATH <- getwd()   # run from the analysis/ directory

# ------------------------------------------------------------------------------
# 1. Kernel matrix K = comm_hel %*% t(comm_hel)  (5000 x 5000)
# ------------------------------------------------------------------------------

message("Reading comm_hel (~2 GB) ...")
comm_hel <- readRDS(file.path(DATA_PATH, "comm_hel_GlobalFungi.rds"))
message("Matrix size: ", nrow(comm_hel), " x ", ncol(comm_hel))

message("Computing K = tcrossprod(comm_hel) ...")
t0 <- proc.time()
K <- tcrossprod(comm_hel)   # 5000 x 5000; usually < 2 min with BLAS
elapsed <- proc.time() - t0
message(sprintf("  done in %.1f s", elapsed["elapsed"]))
message(sprintf("  K size: %.1f MB", object.size(K) / 1e6))

k_path <- file.path(DATA_PATH, "K_GlobalFungi.rds")
saveRDS(K, k_path, compress = FALSE)   # uncompressed for speed
message("Saved: ", k_path)

rm(comm_hel, K); gc()

# ------------------------------------------------------------------------------
# 2. Sparse presence/absence matrix
# ------------------------------------------------------------------------------

message("Reading comm_pa (~250 MB) ...")
comm_pa <- readRDS(file.path(DATA_PATH, "comm_pa_GlobalFungi.rds"))
comm_pa <- comm_pa > 0L

fill_rate <- mean(comm_pa)
message(sprintf("  fill rate: %.4f%% (%s non-zero cells)",
                fill_rate * 100,
                format(sum(comm_pa), big.mark = ",")))

message("Converting to sparse ...")
comm_pa_sp <- Matrix::Matrix(comm_pa, sparse = TRUE)
message(sprintf("  sparse size: %.1f MB (dense: %.1f MB)",
                object.size(comm_pa_sp) / 1e6,
                object.size(comm_pa) / 1e6))

sp_path <- file.path(DATA_PATH, "comm_pa_sp_GlobalFungi.rds")
saveRDS(comm_pa_sp, sp_path, compress = FALSE)
message("Saved: ", sp_path)

rm(comm_pa, comm_pa_sp); gc()

message("\n=== Precomputation finished ===")
message("Next: 03_BAS_simulation_Fig3AB.R (and 07 for the global weight grid search).")
