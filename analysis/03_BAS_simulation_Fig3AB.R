# ==============================================================================
# 03: Case 2 — GlobalFungi (5,000 sites), BAS vs Random vs Oracle
#
# v3 (2026-10): engine = Quantile Regression Forest (ranger, 200 trees, jackknife SE)
#               uniqueness term = novelty (unique-taxa count), see bas_core.R
# Requires: master_data_GlobalFungi.rds, comm_pa_sp_GlobalFungi.rds (from 02)
# Output:   simulation_results_Fig3AB.rds, Fig3AB.png / .pdf
# ==============================================================================
suppressMessages({library(dplyr); library(foreach); library(doParallel); library(Matrix); library(ggplot2); library(patchwork)})
source("bas_core.R")

master_data <- readRDS("master_data_GlobalFungi.rds")
comm_pa     <- readRDS("comm_pa_sp_GlobalFungi.rds")      # sparse presence/absence (02)

N_SIMULATIONS <- as.integer(Sys.getenv("N_SIM", "100")); N_TREES <- as.integer(Sys.getenv("N_TREES", "200"))
N_INITIAL <- 1000; N_TOTAL <- 3000; BATCH_SIZE <- 100; N_BATCHES <- (N_TOTAL - N_INITIAL) / BATCH_SIZE
SIM_OFFSET <- as.integer(Sys.getenv("SIM_OFFSET", "10"))   # evaluation runs use pilot sets 11-110; the grid search (07) used 1-10
W_EQUAL    <- c(1/3, 1/3, 1/3)                               # sensitivity arm with equal weights (Fig. S5)
W_BAS     <- as.numeric(strsplit(Sys.getenv("BAS_WEIGHTS", "0.2,0.6,0.2"), ",")[[1]])   # (W_rich, W_uniq, W_unc) from grid search (07)
N_WORKERS <- as.integer(Sys.getenv("N_WORKERS", "10"))
message(sprintf("N_SIM=%d trees=%d weights: rich=%.2f uniq=%.2f unc=%.2f", N_SIMULATIONS, N_TREES, W_BAS[1], W_BAS[2], W_BAS[3]))

# NOTE (Windows/OpenBLAS): set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 before starting R.
cl <- makeCluster(N_WORKERS, type = "PSOCK", outfile = ""); registerDoParallel(cl)
clusterEvalQ(cl, { source("bas_core.R"); library(Matrix) })
clusterExport(cl, c("master_data", "comm_pa", "N_INITIAL", "BATCH_SIZE", "N_BATCHES", "W_BAS", "W_EQUAL", "N_TREES", "SIM_OFFSET"))
t0 <- Sys.time()
results_df <- foreach(sim = SIM_OFFSET + (1:N_SIMULATIONS), .combine = rbind, .packages = c("dplyr", "Matrix")) %dopar% {
  pl <- draw_pilot_global(sim, master_data, N_INITIAL)
  out <- rbind(campaign_curve("BAS",    "BAS",    sim, master_data, comm_pa, pl$init, pl$rnd_order, BATCH_SIZE, N_BATCHES, engine = "qrf", w = W_BAS, n_trees = N_TREES),
               campaign_curve("BAS",    "BAS (equal weights)", sim, master_data, comm_pa, pl$init, pl$rnd_order, BATCH_SIZE, N_BATCHES, engine = "qrf", w = W_EQUAL, n_trees = N_TREES),
               campaign_curve("Random", "Random", sim, master_data, comm_pa, pl$init, pl$rnd_order, BATCH_SIZE, N_BATCHES),
               campaign_curve("Oracle", "Oracle", sim, master_data, comm_pa, pl$init, pl$rnd_order, BATCH_SIZE, N_BATCHES))
  cat(sprintf("sim %d done (%s)\n", sim, format(Sys.time()))); out
}
stopCluster(cl); message("elapsed: ", format(Sys.time() - t0))
saveRDS(results_df, "simulation_results_Fig3AB.rds")

source("plot_curves.R")
plot_bas_figure(results_df, target_box = c(1200, 1500, 2000, 2500, 3000), xlim = c(N_INITIAL, N_TOTAL), file = "Fig3AB", height = 14)
print(summarise_bas(results_df, c(1200, 1500, 2000, 2500, 3000)))
