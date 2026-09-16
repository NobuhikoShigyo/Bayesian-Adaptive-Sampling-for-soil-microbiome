# ==============================================================================
# BAS Smoke Test — RobustGaSP 6D covariate-aware GP replacement
# N_SIMULATIONS = 10  (smoke test only; production = 100)
# ==============================================================================
# This script replaces gstat local Universal Kriging with RobustGaSP::rgasp
# using a 6D standardized kernel input space: (X, Y, pH, MAT, MAP, SOC).
#
# Key design decisions (see Q&A below):
#
# Q1 — Trend in rgasp():
#   trend = matrix(1, n, 1)   (intercept-only constant mean)
#   Rationale: the 'trend' argument in rgasp() is the parametric GLS mean
#   function matrix (stored as model@X, separate from model@input which is
#   the kernel input).  If covariates are passed as BOTH kernel inputs and
#   trend regressors, the marginal likelihood sees them twice: once through
#   the covariance structure (range/scale hyperparameters fitted by MLE) and
#   once as GLS regression coefficients.  Identifiability breaks — the
#   optimizer can explain any systematic covariate effect either as a mean
#   term or as a long-range correlation.  Passing only an intercept forces all
#   systematic environmental structure into the covariance, which is what we
#   want for covariate-aware posterior variance.
#
# Q2 — predict.rgasp return slot for uncertainty:
#   pred$sd  is the posterior predictive STANDARD DEVIATION (not variance).
#   Source confirmed: in predict.rgasp source, output.list$sd = sqrt(pred_list[[4]])
#   where pred_list[[4]] is the variance from the C++ backend.
#   The existing bas_select_batch() at line 323 uses sqrt(gp_r$var) on gstat's
#   var1.var (a variance slot).  With RobustGaSP the replacement is:
#     unc <- pred$sd        # already SD — no sqrt() needed
#   The new gp_predict_rgasp() returns list(mean = ..., var = pred$sd^2, sd = pred$sd)
#   so the caller line `unc <- sqrt(gp_r$var)` remains correct (sqrt of variance = SD).
#
# Q3 — Both gp_r and gp_l use 6D GP:
#   Yes.  Both TrueRichness (gp_r) and LCBD (gp_l) are fitted with the same
#   gp_predict_rgasp() function.  The training matrices are assembled
#   identically for both: rows = sampled sites, columns = standardized
#   (X, Y, pH, MAT, MAP, SOC).  The only difference is the response vector
#   (TrueRichness vs curr_lcbd).  LCBD is a scalar per sampled site computed
#   from comm_hel BEFORE calling gp_predict_rgasp(), so by the time the GP
#   sees it, it is just another numeric response — no asymmetry in the
#   function path.
#
# Q4 — Memory budget (post-2026-05-13 crash):
#   Original worker bundle = 3.95 GB uncompressed × 4 PSOCK workers ≈ 16 GB
#   resident across workers, plus another ~4 GB on master.  PC crashed mid-run.
#   Mitigations applied here:
#     (a) comm_pa and comm_hel converted to Matrix::dgCMatrix sparse storage
#         BEFORE being saved to the worker bundle (decostand and rowsq are
#         still computed on dense, then dense is freed via gc()).
#     (b) saveRDS(... compress = TRUE) for the worker bundle on disk.
#     (c) Per-worker and per-5-batch memory snapshots logged via gc()[,2]
#         summed → GB, so a future crash can be diagnosed quantitatively.
# ==============================================================================

# --- Libraries ---
if (!require(RobustGaSP,    quietly = TRUE)) install.packages("RobustGaSP")
if (!require(RhpcBLASctl,   quietly = TRUE)) install.packages("RhpcBLASctl")

library(RobustGaSP)
library(RhpcBLASctl)
library(vegan)
library(dplyr)
library(ggplot2)
library(adespatial)
library(sp)
library(sf)
library(tidyr)
library(foreach)
library(doParallel)
library(data.table)
library(Matrix)   # sparse matrix storage for comm_pa / comm_hel (Q4)

# --- Memory snapshot helper (returns total R-heap usage in GB) ---
mem_gb <- function() round(sum(gc(verbose = FALSE)[, 2]) / 1024, 2)

# set.seed(123)  # Reproducibility handled by clusterSetRNGStream below

setDTthreads(1)
blas_set_num_threads(1)
omp_set_num_threads(1)

# ==============================================================================
# Step 0: Load the Case 2 objects written by 03_BAS_simulation_Fig3AB.R Step 1
# ==============================================================================
# The public repository ships the transformed objects (presence/absence and
# Hellinger matrices) and a column-trimmed site table instead of the raw
# GlobalFungi read-count subset this script originally read.  They are exactly
# what this script used to derive from the raw table (same row/column order;
# Hellinger values agree to 1e-16), so the results are unchanged.
# Note: 03 stores X/Y in metres; this script works in km.
DATA_PATH <- ".."
for (f in c("master_data_GlobalFungi.rds", "comm_pa_GlobalFungi.rds",
            "comm_hel_GlobalFungi.rds")) {
  if (!file.exists(file.path(DATA_PATH, f)))
    stop(f, " not found in ", normalizePath(DATA_PATH),
         " -- run 03_BAS_simulation_Fig3AB.R Step 1 first.")
}

master_data <- readRDS(file.path(DATA_PATH, "master_data_GlobalFungi.rds"))
comm_pa     <- readRDS(file.path(DATA_PATH, "comm_pa_GlobalFungi.rds"))
comm_hel    <- readRDS(file.path(DATA_PATH, "comm_hel_GlobalFungi.rds"))
stopifnot(identical(master_data$sample_ID, rownames(comm_pa)),
          identical(rownames(comm_pa),     rownames(comm_hel)))
storage.mode(comm_pa) <- "double"   # as produced by (comm_mat > 0) * 1

master_data <- master_data %>%
  mutate(X = X / 1000, Y = Y / 1000) %>%   # metres -> km (kernel input scale)
  as.data.frame()

drift_vars <- c("pH", "MAT", "MAP", "SOC", "X", "Y", "TrueRichness")
ok_rows    <- stats::complete.cases(master_data[, drift_vars])
if (!all(ok_rows)) {
  master_data <- master_data[ok_rows, , drop = FALSE]
  comm_pa     <- comm_pa[ok_rows, , drop = FALSE]
  comm_hel    <- comm_hel[ok_rows, , drop = FALSE]
  message("Dropped ", sum(!ok_rows), " rows with NA covariates -- N = ",
          nrow(master_data))
}

comm_hel_rowsq <- rowSums(comm_hel^2)

# --- Sparsify comm_pa and comm_hel before sending to workers (Q4) ---
# decostand and rowsq above are computed on the dense forms (faster), and we
# convert in place + gc() to release the dense memory before pickling.
message("Master mem before sparsify: ", mem_gb(), " GB")
comm_pa  <- Matrix::Matrix(comm_pa,  sparse = TRUE)
comm_hel <- Matrix::Matrix(comm_hel, sparse = TRUE)
gc(FALSE)
message("Master mem after  sparsify: ", mem_gb(), " GB",
        "  | comm_pa nnz=", length(comm_pa@x),
        "  comm_hel nnz=", length(comm_hel@x))

# --- Dump large matrices to disk so workers load via RDS (avoids socket serialization) ---
.worker_bundle_path <- "_worker_data.rds"
saveRDS(
  list(master_data     = master_data,
       comm_pa         = comm_pa,
       comm_hel        = comm_hel,
       comm_hel_rowsq  = comm_hel_rowsq),
  file = .worker_bundle_path,
  compress = TRUE   # Q4: was FALSE — disk I/O slightly slower, ~10× smaller
)
message("Worker bundle saved: ",
        round(file.size(.worker_bundle_path) / 1e6, 1), " MB",
        " (compressed)")

# ==============================================================================
# Simulation settings — SMOKE TEST values
# ==============================================================================
N_SIMULATIONS <- 10          # SMOKE TEST: 10 (production = 100)
N_INITIAL     <- 100         # matches production
N_TOTAL       <- 3000
BATCH_SIZE    <- 100
N_BATCHES     <- floor((N_TOTAL - N_INITIAL) / BATCH_SIZE)   # = 29

# RobustGaSP scalability knob:
# Training is capped at RGASP_FIT_N rows per batch.  At N=1000, Cholesky of
# the 1000×1000 correlation matrix takes ~0.1 s; at n_sampled growing to
# 3000 it would be ~2 s.  We cap at 1000 to keep each batch under ~2 min.
RGASP_FIT_N   <- 1000

# Weights: same as the v2 GP script
W_R   <- 0.35
W_U   <- 0.35
W_UNC <- 0.30

# ==============================================================================
# Helper functions
# ==============================================================================

norm01 <- function(x) {
  finite_x <- x[is.finite(x)]
  if (length(finite_x) == 0 || max(finite_x) == min(finite_x))
    return(rep(0, length(x)))
  out <- (x - min(finite_x)) / (max(finite_x) - min(finite_x))
  out[!is.finite(out)] <- 0
  out
}

compute_lcbd_fast <- function(sampled_idx) {
  sub      <- comm_hel[sampled_idx, , drop = FALSE]
  centroid <- colMeans(sub)
  xmu      <- as.vector(sub %*% centroid)
  rm(sub); invisible(gc(FALSE))
  per_row_ss <- comm_hel_rowsq[sampled_idx] - 2 * xmu + sum(centroid^2)
  per_row_ss[per_row_ss < 0] <- 0
  ss_total <- sum(per_row_ss)
  if (ss_total > 0) per_row_ss / ss_total else rep(0, length(sampled_idx))
}

# ------------------------------------------------------------------------------
# gp_predict_rgasp():
#   Fits RobustGaSP with 6D standardized kernel input (X, Y, pH, MAT, MAP, SOC)
#   and intercept-only trend (constant mean function).
#
#   Args:
#     train_X6  : n_train × 6 matrix of raw covariate values (training)
#     train_y   : length-n_train response vector
#     cand_X6   : n_cand × 6 matrix of raw covariate values (candidates)
#     fit_n     : max training rows (subsample if n_train > fit_n)
#
#   Returns: list(mean, var, sd) at candidate sites, or NULL on failure.
#     'var'  = pred$sd^2  (posterior predictive variance)
#     'sd'   = pred$sd    (posterior predictive standard deviation)
#   The caller uses sqrt(gp_r$var) = gp_r$sd as the uncertainty signal,
#   consistent with the existing bas_select_batch() code structure.
# ------------------------------------------------------------------------------
gp_predict_rgasp <- function(train_X6, train_y, cand_X6, fit_n = RGASP_FIT_N) {

  n_tr <- nrow(train_X6)

  # --- Subsample training rows if needed ---
  if (n_tr > fit_n) {
    idx_fit   <- sample.int(n_tr, fit_n)
    train_X6  <- train_X6[idx_fit, , drop = FALSE]
    train_y   <- train_y[idx_fit]
    n_tr      <- fit_n
  }

  # --- Standardize kernel input using training-set moments ---
  # Centering and scaling in 6D keeps all dimensions numerically comparable.
  # We apply the SAME transformation to candidates.
  col_means <- colMeans(train_X6)
  col_sds   <- apply(train_X6, 2, sd)
  col_sds[col_sds == 0] <- 1   # guard against degenerate columns

  X_sc      <- sweep(sweep(train_X6, 2, col_means, "-"), 2, col_sds, "/")
  X_cand_sc <- sweep(sweep(cand_X6,  2, col_means, "-"), 2, col_sds, "/")

  # --- Trend matrix: intercept only (see Q1 note at top of file) ---
  trend_train <- matrix(1, n_tr,          1)
  trend_cand  <- matrix(1, nrow(X_cand_sc), 1)

  # --- Fit RobustGaSP ---
  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      rgasp(
        design      = X_sc,
        response    = train_y,
        trend       = trend_train,
        zero.mean   = "No",
        nugget.est  = TRUE,     # estimate noise floor; important at global scale
        kernel_type = "matern_5_2",
        isotropic   = FALSE     # separate range per dimension → env vs space
      )
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  # --- Predict at candidates ---
  pred <- tryCatch(
    suppressMessages(suppressWarnings(
      predict(fit, testing_input = X_cand_sc, testing_trend = trend_cand)
    )),
    error = function(e) NULL
  )
  if (is.null(pred)) return(NULL)

  # pred$sd is the posterior predictive standard deviation (confirmed from
  # predict.rgasp source: output.list$sd = sqrt(pred_list[[4]]) where
  # pred_list[[4]] is the C++ variance output).
  sd_vec <- pred$sd
  sd_vec[is.na(sd_vec) | sd_vec < 0] <- 0

  list(
    mean = pred$mean,
    var  = sd_vec^2,   # variance  — sqrt(gp_r$var) in caller = sd_vec
    sd   = sd_vec      # SD        — convenience slot
  )
}

# ------------------------------------------------------------------------------
# bas_select_batch(): identical logic to production script but calls
# gp_predict_rgasp() instead of gp_predict_local().
# Per-batch unc_cv is computed from the raw SD vector (= gp_r$sd).
# ------------------------------------------------------------------------------
bas_select_batch <- function(curr_sam, curr_uns, conf) {

  curr_lcbd <- compute_lcbd_fast(curr_sam)

  # --- Assemble 6D covariate matrix for training and candidate sets ---
  COV_COLS <- c("X", "Y", "pH", "MAT", "MAP", "SOC")

  train_X6  <- as.matrix(master_data[curr_sam, COV_COLS])
  cand_X6   <- as.matrix(master_data[curr_uns,  COV_COLS])

  train_r   <- master_data$TrueRichness[curr_sam]
  train_l   <- curr_lcbd

  # --- Fit GP for richness (gp_r) and LCBD (gp_l) ---
  # Both use the identical gp_predict_rgasp() path (see Q3 note).
  gp_r <- gp_predict_rgasp(train_X6, train_r, cand_X6)
  gp_l <- gp_predict_rgasp(train_X6, train_l, cand_X6)

  rm(train_X6, cand_X6); invisible(gc(FALSE))

  if (is.null(gp_r) || is.null(gp_l)) return(NULL)

  # unc = sqrt(variance) = gp_r$sd  (posterior predictive SD for richness)
  # Using sqrt(gp_r$var) preserves the existing caller convention and is
  # numerically identical to gp_r$sd.
  unc <- sqrt(gp_r$var)

  sc <- conf$Wr   * norm01(gp_r$mean) +
        conf$Wu   * norm01(gp_l$mean) +
        conf$Wunc * norm01(unc)

  cv_unc <- if (mean(unc) > 0) sd(unc) / mean(unc) else 0
  attr(sc, "unc_cv") <- cv_unc

  ord <- order(sc, decreasing = TRUE)[1:BATCH_SIZE]
  attr(ord, "unc_cv") <- cv_unc
  ord
}

# ==============================================================================
# Simulation loop — BAS arm only (no QRF comparison in smoke test)
# ==============================================================================

strategies <- list(
  BAS_Balanced = list(Wr = W_R, Wu = W_U, Wunc = W_UNC, Name = "BAS", Type = "Env")
)

n_workers <- 4L
cl <- parallel::makeCluster(n_workers, type = "PSOCK")
doParallel::registerDoParallel(cl)

RNGkind("L'Ecuyer-CMRG")
parallel::clusterSetRNGStream(cl, iseed = 123)

parallel::clusterExport(cl, varlist = ".worker_bundle_path", envir = environment())

parallel::clusterEvalQ(cl, {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  data.table::setDTthreads(1)
  library(Matrix)   # needed to deserialize dgCMatrix S4 objects
  .worker_bundle <- readRDS(.worker_bundle_path)
  master_data    <- .worker_bundle$master_data
  comm_pa        <- .worker_bundle$comm_pa
  comm_hel       <- .worker_bundle$comm_hel
  comm_hel_rowsq <- .worker_bundle$comm_hel_rowsq
  rm(.worker_bundle); gc(FALSE)
  .worker_mem_after_load <- round(sum(gc(verbose = FALSE)[, 2]) / 1024, 2)
  message("  Worker init: bundle loaded, mem = ",
          .worker_mem_after_load, " GB  pid = ", Sys.getpid())
  TRUE
})

message("Starting SMOKE TEST: N_SIM=", N_SIMULATIONS,
        " N_INITIAL=", N_INITIAL, " N_TOTAL=", N_TOTAL,
        " BATCH=", BATCH_SIZE, " N_BATCHES=", N_BATCHES,
        " RGASP_FIT_N=", RGASP_FIT_N)

t0 <- Sys.time()

parallel::clusterExport(cl, varlist = c(
  "gp_predict_rgasp", "bas_select_batch", "norm01", "compute_lcbd_fast",
  "strategies", "N_INITIAL", "N_TOTAL", "BATCH_SIZE", "N_BATCHES", "RGASP_FIT_N",
  "W_R", "W_U", "W_UNC"
), envir = environment())

sim_results <- foreach::foreach(
  sim = seq_len(N_SIMULATIONS),
  .packages = c("RobustGaSP", "RhpcBLASctl", "vegan", "dplyr", "adespatial", "data.table", "Matrix"),
  .noexport = c("master_data", "comm_pa", "comm_hel", "comm_hel_rowsq")
) %dopar% {

  set.seed(123 + sim)   # exactly as production script

  prob_weights <- ifelse(master_data$latitude > 20, 1.0, 0.05)
  initial_idx  <- sample.int(nrow(master_data), N_INITIAL, prob = prob_weights)
  remain_idx   <- (1:nrow(master_data))[-initial_idx]

  init_species_count <- colSums(comm_pa[initial_idx, , drop = FALSE])
  init_r             <- sum(init_species_count > 0)

  sim_rows        <- list()
  unc_cv_rows_sim <- list()

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

      if (length(curr_uns) < BATCH_SIZE) break

      best_idx <- tryCatch(
        bas_select_batch(curr_sam, curr_uns, conf),
        error = function(e) NULL
      )

      cv_val <- NA_real_
      if (is.null(best_idx)) {
        best_idx <- sample.int(length(curr_uns), BATCH_SIZE)
        message(sprintf("  sim %d  batch %d  FALLBACK to random (GP failed)", sim, b))
      } else {
        cv_val <- attr(best_idx, "unc_cv")
        # Log EVERY batch (not every 5th — smoke-test diagnostic requirement)
        message(sprintf("  sim %d  batch %2d  unc_cv = %.4f", sim, b, cv_val))
      }

      # Record unc_cv for this sim × batch
      unc_cv_rows_sim[[length(unc_cv_rows_sim) + 1]] <- data.frame(
        SimID  = sim,
        Batch  = b,
        unc_cv = cv_val
      )

      next_nodes         <- curr_uns[best_idx]
      curr_uns           <- curr_uns[-best_idx]
      curr_sam           <- c(curr_sam, next_nodes)
      new_counts         <- colSums(comm_pa[next_nodes, , drop = FALSE])
      curr_species_count <- curr_species_count + new_counts
      n_vec[b + 1]       <- N_INITIAL + b * BATCH_SIZE
      r_vec[b + 1]       <- sum(curr_species_count > 0)

      if (b %% 5 == 0) {
        g <- gc(verbose = FALSE)
        message(sprintf("  sim %d  batch %2d  MEM = %.2f GB",
                        sim, b, sum(g[, 2]) / 1024))
      }
    }

    sim_rows[[length(sim_rows) + 1]] <- data.frame(
      n_samples = n_vec,
      Richness  = r_vec,
      SimID     = sim,
      Method    = conf$Name
    )
  }

  invisible(gc(FALSE))

  list(
    result_row   = do.call(rbind, sim_rows),
    unc_cv_rows  = do.call(rbind, unc_cv_rows_sim)
  )
}

parallel::stopCluster(cl)
unlink(.worker_bundle_path)
message("Worker bundle cleaned up.")

total_time <- round(difftime(Sys.time(), t0, units = "mins"), 1)
message("All simulations complete. Total time: ", total_time, " min")

results_list  <- lapply(sim_results, `[[`, "result_row")
unc_cv_records <- do.call(rbind, lapply(sim_results, `[[`, "unc_cv_rows"))

# ==============================================================================
# Assemble outputs
# ==============================================================================

results_df  <- do.call(rbind, results_list)
unc_cv_df   <- unc_cv_records   # already a data.frame from do.call(rbind, ...) above

# --- Save RDS outputs ---
OUT_DIR <- "../output"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
saveRDS(results_df, file.path(OUT_DIR, "smoke_test_results.rds"))
saveRDS(unc_cv_df,  file.path(OUT_DIR, "smoke_test_unc_cv.rds"))
message("Saved: ", OUT_DIR, "/smoke_test_results.rds  smoke_test_unc_cv.rds")

# ==============================================================================
# Summary printout
# ==============================================================================

# 1. Mean cumulative richness at n=3000 for BAS arm
rich_3000 <- results_df %>%
  filter(Method == "BAS", n_samples == N_TOTAL) %>%
  summarise(
    Mean_Richness = mean(Richness),
    SD_Richness   = sd(Richness),
    N_sims        = n()
  )
cat("\n=== 1. Mean cumulative richness at n=3000 (BAS) ===\n")
print(rich_3000)

# 2. Mean unc_cv per batch (all 20 batches)
cv_by_batch <- unc_cv_df %>%
  filter(!is.na(unc_cv)) %>%
  group_by(Batch) %>%
  summarise(
    Mean_unc_cv = mean(unc_cv),
    Min_unc_cv  = min(unc_cv),
    Max_unc_cv  = max(unc_cv),
    .groups     = "drop"
  ) %>%
  arrange(Batch)

cat("\n=== 2. Mean unc_cv per batch (all batches) ===\n")
print(as.data.frame(cv_by_batch), digits = 4)

# 3. Batches where mean unc_cv < 0.10 (degeneration criterion)
degenerate_batches <- cv_by_batch %>% filter(Mean_unc_cv < 0.10)
cat("\n=== 3. Batches with mean unc_cv < 0.10 (degeneration criterion) ===\n")
if (nrow(degenerate_batches) == 0) {
  cat("None — uncertainty signal remains discriminative through all batches.\n")
} else {
  cat("WARNING: the following batches show degeneration:\n")
  print(as.data.frame(degenerate_batches), digits = 4)
}

cat("\n=== Wall time: ", as.character(total_time), " min ===\n")
cat("=== Output files: smoke_test_results.rds  smoke_test_unc_cv.rds ===\n")
