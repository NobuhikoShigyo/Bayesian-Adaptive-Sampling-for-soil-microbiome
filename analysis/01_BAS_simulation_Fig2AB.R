# ==============================================================================
# 01: Case 1 — local catchment (Shigyo et al. 2022), BAS vs Random vs Oracle
#
# v3 (2026-10): engine = regression kriging (automap variogram, gstat::krige)
#               uniqueness term = novelty (unique-taxa count), see bas_core.R
# Output: simulation_results_Fig2AB.rds, Fig2AB.png / .pdf (A: curves incl. Oracle, B: BAS vs Random boxes)
# ==============================================================================
suppressMessages({library(vegan); library(dplyr); library(foreach); library(doParallel); library(Matrix); library(ggplot2); library(patchwork)})
source("bas_core.R")

# ---- Data (raw CSVs -> master_data / comm_pa; cached as rds) -----------------
if (!file.exists("master_data_local.rds")) {
  asv_raw   <- read.csv("comm_data_Shigyo_et_al_2022.csv", row.names = 1, header = TRUE, check.names = FALSE)
  site_data <- read.csv("site_data_Shigyo_et_al_2022.csv", header = TRUE)
  common    <- intersect(rownames(asv_raw), site_data$SiteID)
  comm_data <- as.matrix(asv_raw[common, ]); comm_data <- comm_data[, colSums(comm_data) > 0]
  comm_pa   <- (comm_data > 0) * 1
  master_data <- site_data %>% filter(SiteID %in% common) %>% arrange(match(SiteID, common)) %>%
    left_join(data.frame(SiteID = rownames(comm_pa), TrueRichness = rowSums(comm_pa)), by = "SiteID") %>%
    rename(X = lon, Y = lat)
  saveRDS(master_data, "master_data_local.rds"); saveRDS(comm_data, "comm_data_local.rds"); saveRDS(comm_pa, "comm_data_pa_local.rds")
}
master_data <- readRDS("master_data_local.rds"); comm_pa <- readRDS("comm_data_pa_local.rds")

# ---- Settings ----------------------------------------------------------------
N_SIMULATIONS <- 100; N_INITIAL <- 5; BATCH_SIZE <- 5; N_TOTAL <- 30
SIM_OFFSET    <- as.integer(Sys.getenv("SIM_OFFSET", "50"))   # evaluation runs use pilot sets 51-150; the grid search (05) used 1-50
W_EQUAL       <- c(1/3, 1/3, 1/3)                                # sensitivity arm with equal weights (Fig. S5)
N_BATCHES     <- (N_TOTAL - N_INITIAL) / BATCH_SIZE
W_BAS         <- as.numeric(strsplit(Sys.getenv("BAS_WEIGHTS", "0.4,0.2,0.4"), ",")[[1]])   # (W_rich, W_uniq, W_unc) from grid search (05)
N_WORKERS     <- as.integer(Sys.getenv("N_WORKERS", "6"))
message(sprintf("weights: rich=%.2f uniq=%.2f unc=%.2f", W_BAS[1], W_BAS[2], W_BAS[3]))

# ---- Simulation ---------------------------------------------------------------
cl <- makeCluster(N_WORKERS, type = "PSOCK"); registerDoParallel(cl)
clusterEvalQ(cl, { source("bas_core.R"); library(Matrix) })
clusterExport(cl, c("master_data", "comm_pa", "N_INITIAL", "BATCH_SIZE", "N_BATCHES", "W_BAS", "W_EQUAL", "SIM_OFFSET"))
results_df <- foreach(sim = SIM_OFFSET + (1:N_SIMULATIONS), .combine = rbind, .packages = c("dplyr", "Matrix")) %dopar% {
  pl <- draw_pilot_local(sim, nrow(master_data), N_INITIAL)
  rbind(campaign_curve("BAS",    "BAS",    sim, master_data, comm_pa, pl$init, pl$rnd_order, BATCH_SIZE, N_BATCHES, engine = "rk", w = W_BAS),
        campaign_curve("BAS",    "BAS (equal weights)", sim, master_data, comm_pa, pl$init, pl$rnd_order, BATCH_SIZE, N_BATCHES, engine = "rk", w = W_EQUAL),
        campaign_curve("Random", "Random", sim, master_data, comm_pa, pl$init, pl$rnd_order, BATCH_SIZE, N_BATCHES),
        campaign_curve("Oracle", "Oracle", sim, master_data, comm_pa, pl$init, pl$rnd_order, BATCH_SIZE, N_BATCHES))
}
stopCluster(cl)
saveRDS(results_df, "simulation_results_Fig2AB.rds")

# ---- Figure ------------------------------------------------------------------
source("plot_curves.R")
plot_bas_figure(results_df, target_box = c(10, 15, 20, 25, 30), xlim = c(N_INITIAL, N_TOTAL), file = "Fig2AB")
print(summarise_bas(results_df, c(10, 15, 20, 25, 30)))
