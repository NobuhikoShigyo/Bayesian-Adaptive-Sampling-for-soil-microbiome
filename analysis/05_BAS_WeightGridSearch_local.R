# ==============================================================================
# 05: Weight grid search — Case 1 (local), 66 points on the simplex (step 0.1)
#
# v3 (2026-10): engine = regression kriging, uniqueness = novelty (bas_core.R).
# Pilot sets identical to 01 (seed 123 + sim). Output: gridsearch_results_local.rds
# Plot with 06_plot_WeightGridSearch.R.
# ==============================================================================
suppressMessages({library(dplyr); library(foreach); library(doParallel); library(Matrix)})
source("bas_core.R")
master_data <- readRDS("master_data_local.rds"); comm_pa <- readRDS("comm_data_pa_local.rds")
N_SIMULATIONS <- as.integer(Sys.getenv("N_SIM", "50")); N_INITIAL <- 5; BATCH_SIZE <- 5; N_BATCHES <- 5
N_WORKERS <- as.integer(Sys.getenv("N_WORKERS", "6"))
grid <- weight_grid(0.1); message("weight combinations: ", nrow(grid), " x ", N_SIMULATIONS, " simulations")

cl <- makeCluster(N_WORKERS, type = "PSOCK"); registerDoParallel(cl)
clusterEvalQ(cl, { source("bas_core.R"); library(Matrix) })
clusterExport(cl, c("master_data", "comm_pa", "N_SIMULATIONS", "N_INITIAL", "BATCH_SIZE", "N_BATCHES", "grid"))
t0 <- Sys.time()
results_grid <- foreach(w = seq_len(nrow(grid)), .combine = rbind, .packages = c("dplyr", "Matrix")) %dopar% {
  wv <- c(grid$Wr[w], grid$Wu[w], grid$Wunc[w]); out <- vector("list", N_SIMULATIONS)
  for (sim in seq_len(N_SIMULATIONS)) {
    pl <- draw_pilot_local(sim, nrow(master_data), N_INITIAL)
    cc <- campaign_curve("BAS", "BAS", sim, master_data, comm_pa, pl$init, pl$rnd_order, BATCH_SIZE, N_BATCHES, engine = "rk", w = wv)
    out[[sim]] <- data.frame(Wr = wv[1], Wu = wv[2], Wunc = wv[3], sim = sim, n_samples = cc$n_samples, Richness = cc$Richness)
  }
  do.call(rbind, out)
}
stopCluster(cl); message("elapsed: ", format(Sys.time() - t0))
saveRDS(results_grid, "gridsearch_results_local.rds")
print(results_grid %>% filter(n_samples == 30) %>% group_by(Wr, Wu, Wunc) %>% summarise(mean_richness = mean(Richness), .groups = "drop") %>% arrange(desc(mean_richness)) %>% head(10))
