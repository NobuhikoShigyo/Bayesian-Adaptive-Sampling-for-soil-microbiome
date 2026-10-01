# ==============================================================================
# 07: Weight grid search — Case 2 (GlobalFungi), 66 points on the simplex (step 0.1)
#
# v3 (2026-10): engine = QRF (100 trees for run time), uniqueness = novelty (bas_core.R).
# Pilot sets identical to 03 (seed 123 + sim). 10 simulations per combination.
# Output: gridsearch_results_global.rds.  Plot with 06_plot_WeightGridSearch.R.
# NOTE (Windows/OpenBLAS): set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 before starting R.
# ==============================================================================
suppressMessages({library(dplyr); library(foreach); library(doParallel); library(Matrix)})
source("bas_core.R")
master_data <- readRDS("master_data_GlobalFungi.rds"); comm_pa <- readRDS("comm_pa_sp_GlobalFungi.rds")
N_SIMULATIONS <- as.integer(Sys.getenv("N_SIM", "10")); N_TREES <- as.integer(Sys.getenv("N_TREES", "100"))
N_INITIAL <- 1000; BATCH_SIZE <- 100; N_BATCHES <- 20
N_WORKERS <- as.integer(Sys.getenv("N_WORKERS", "8"))
grid <- weight_grid(0.1); message("weight combinations: ", nrow(grid), " x ", N_SIMULATIONS, " simulations, ", N_TREES, " trees")

cl <- makeCluster(N_WORKERS, type = "PSOCK", outfile = ""); registerDoParallel(cl)
clusterEvalQ(cl, { source("bas_core.R"); library(Matrix) })
clusterExport(cl, c("master_data", "comm_pa", "N_SIMULATIONS", "N_TREES", "N_INITIAL", "BATCH_SIZE", "N_BATCHES", "grid"))
t0 <- Sys.time()
results_grid <- foreach(w = seq_len(nrow(grid)), .combine = rbind, .packages = c("dplyr", "Matrix")) %dopar% {
  wv <- c(grid$Wr[w], grid$Wu[w], grid$Wunc[w]); out <- vector("list", N_SIMULATIONS)
  for (sim in seq_len(N_SIMULATIONS)) {
    pl <- draw_pilot_global(sim, master_data, N_INITIAL)
    cc <- campaign_curve("BAS", "BAS", sim, master_data, comm_pa, pl$init, pl$rnd_order, BATCH_SIZE, N_BATCHES, engine = "qrf", w = wv, n_trees = N_TREES)
    out[[sim]] <- data.frame(Wr = wv[1], Wu = wv[2], Wunc = wv[3], sim = sim, n_samples = cc$n_samples, Richness = cc$Richness)
  }
  cat(sprintf("combo %d/%d done (%s)\n", w, nrow(grid), paste(wv, collapse = "/"))); do.call(rbind, out)
}
stopCluster(cl); message("elapsed: ", format(Sys.time() - t0))
saveRDS(results_grid, "gridsearch_results_global.rds")
print(results_grid %>% filter(n_samples == 3000) %>% group_by(Wr, Wu, Wunc) %>% summarise(mean_richness = mean(Richness), .groups = "drop") %>% arrange(desc(mean_richness)) %>% head(10))
