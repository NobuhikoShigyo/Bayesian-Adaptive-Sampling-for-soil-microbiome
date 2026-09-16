# ==============================================================================
# BAS Acquisition Function Weight Grid Search — Ternary Plot Visualization
#
# Grid search: W_rich + W_uniq + W_unc = 1.0 (step 0.1) → 66 combinations
# Each combination: 100 simulations on local dataset (N=53, Shigyo et al. 2022)
# Evaluation metric: cumulative richness at n = 30
#
# Output:
#   BAS_WeightOptimization_TernaryPlot.png
#   gridsearch_results.rds
#   gridsearch_summary.csv
# ==============================================================================

# NOTE: ggtern conflicts with ggplot2 if loaded at the top.
#       Load order: ggplot2-based packages first, ggtern only at plot section.

library(vegan)
library(dplyr)
library(gstat)
library(adespatial)
library(sp)
library(foreach)
library(doParallel)

set.seed(123)

# ==============================================================================
# Step 1: Load and prepare empirical data (identical to BAS_simulation_Fig2AB.R)
# ==============================================================================

asv_table_raw <- read.csv("comm_data_Shigyo_et_al_2022.csv",
                           row.names = 1, header = TRUE, check.names = FALSE)
site_data_raw <- read.csv("site_data_Shigyo_et_al_2022.csv", header = TRUE)

common_sites <- intersect(rownames(asv_table_raw), site_data_raw$SiteID)
if (length(common_sites) == 0) stop("Error: no matching SiteIDs.")

asv_table <- asv_table_raw[common_sites, ]
site_data_raw <- site_data_raw %>%
  filter(SiteID %in% common_sites) %>%
  arrange(match(SiteID, common_sites))
asv_table <- asv_table[, colSums(asv_table) > 0]

comm_data    <- as.matrix(asv_table)
comm_data_pa <- (comm_data > 0) * 1L
all_richness <- data.frame(
  SiteID       = rownames(comm_data_pa),
  TrueRichness = rowSums(comm_data_pa)
)

master_data <- site_data_raw %>%
  left_join(all_richness, by = "SiteID") %>%
  rename(X = lon, Y = lat)

message("Dataset ready: ", nrow(master_data), " sites")

# ==============================================================================
# Step 2: Simulation settings & weight grid
# ==============================================================================

N_SIMULATIONS     <- 100
N_INITIAL_SAMPLES <- 5
N_TOTAL_SAMPLES   <- 50   # unchanged from Fig2AB
BATCH_SIZE        <- 5    # unchanged from Fig2AB
N_BATCHES         <- floor((N_TOTAL_SAMPLES - N_INITIAL_SAMPLES) / BATCH_SIZE)
EVAL_N            <- 30   # richness evaluation point

# Generate all (Wr, Wu, Wunc) triples on 0.1 grid that sum to 1.0 → 66 rows
weight_grid <- expand.grid(
  Wr   = seq(0, 1, by = 0.1),
  Wu   = seq(0, 1, by = 0.1),
  Wunc = seq(0, 1, by = 0.1)
) %>%
  dplyr::filter(round(Wr + Wu + Wunc, 10) == 1.0)

stopifnot(nrow(weight_grid) == 66)
message("Weight combinations: ", nrow(weight_grid),
        "  (", N_SIMULATIONS, " sims each = ",
        nrow(weight_grid) * N_SIMULATIONS, " total runs)")

# ==============================================================================
# Step 3: Grid search simulation (outer loop: weight combinations, parallelised)
# ==============================================================================

num_cores <- max(1L, parallel::detectCores() - 1L)
cl <- makeCluster(num_cores, type = "PSOCK", outfile = "")
registerDoParallel(cl)

message("Starting grid search on ", num_cores, " cores ...")

results_grid <- foreach(
  w         = seq_len(nrow(weight_grid)),
  .combine  = rbind,
  .packages = c("vegan", "dplyr", "gstat", "adespatial", "sp"),
  .export   = c("comm_data", "comm_data_pa", "master_data",
                "weight_grid",
                "N_INITIAL_SAMPLES", "N_TOTAL_SAMPLES",
                "BATCH_SIZE", "N_BATCHES", "N_SIMULATIONS", "EVAL_N")
) %dopar% {

  Wr   <- weight_grid$Wr[w]
  Wu   <- weight_grid$Wu[w]
  Wunc <- weight_grid$Wunc[w]

  # Normalise a numeric vector to [0, 1]; returns zeros for constant/empty input
  normalize <- function(x) {
    if (length(x) == 0 || all(is.na(x))) return(rep(0, length(x)))
    x[is.na(x)] <- min(x, na.rm = TRUE)
    if (max(x) == min(x)) return(rep(0, length(x)))
    (x - min(x)) / (max(x) - min(x))
  }

  # ---- Inner loop: N_SIMULATIONS repetitions for this weight combination ----
  sim_rows <- vector("list", N_SIMULATIONS)

  for (sim in seq_len(N_SIMULATIONS)) {
    set.seed(123 + sim)

    initial_indices   <- sample(seq_len(nrow(master_data)), N_INITIAL_SAMPLES)
    remaining_pool    <- setdiff(seq_len(nrow(master_data)), initial_indices)
    initial_richness  <- sum(colSums(comm_data_pa[initial_indices, , drop = FALSE]) > 0)

    current_sampled   <- initial_indices
    current_unsampled <- remaining_pool

    curve_rows <- vector("list", N_BATCHES + 1L)
    curve_rows[[1L]] <- data.frame(
      Wr = Wr, Wu = Wu, Wunc = Wunc, sim = sim,
      n_samples = N_INITIAL_SAMPLES, Richness = initial_richness
    )

    for (b in seq_len(N_BATCHES)) {
      if (length(current_unsampled) < BATCH_SIZE) break

      # Compute LCBD on currently sampled sites
      curr_comm <- comm_data[current_sampled, , drop = FALSE]
      curr_lcbd <- tryCatch({
        if (nrow(curr_comm) < 3) rep(0, nrow(curr_comm))
        else adespatial::beta.div(
          vegan::decostand(curr_comm, "hellinger"),
          method = "hellinger", nperm = 0)$LCBD
      }, error = function(e) rep(0, length(current_sampled)))

      train_d       <- master_data[current_sampled, ]
      train_d$LCBD  <- curr_lcbd
      coordinates(train_d) <- ~X + Y

      cand_d <- master_data[current_unsampled, ]
      coordinates(cand_d) <- ~X + Y

      g_r <- try(gstat(formula = TrueRichness ~ pH + WC,
                       locations = train_d, nmax = 10, set = list(idp = .5)),
                 silent = TRUE)
      g_l <- try(gstat(formula = LCBD ~ pH + WC,
                       locations = train_d, nmax = 10, set = list(idp = .5)),
                 silent = TRUE)

      best_idx <- NULL
      if (!inherits(g_r, "try-error") && !inherits(g_l, "try-error")) {
        pred_r <- try(
          suppressMessages(predict(g_r, newdata = cand_d, debug.level = -1)),
          silent = TRUE)
        pred_l <- try(
          suppressMessages(predict(g_l, newdata = cand_d, debug.level = -1)),
          silent = TRUE)

        if (!inherits(pred_r, "try-error") && !inherits(pred_l, "try-error")) {
          var_r <- if ("var1.var" %in% names(pred_r)) pred_r$var1.var
                   else rep(0, length(pred_r$var1.pred))
          var_r[var_r < 0] <- 0

          sc <- Wr   * normalize(pred_r$var1.pred) +
                Wu   * normalize(pred_l$var1.pred) +
                Wunc * normalize(sqrt(var_r))
          best_idx <- order(sc, decreasing = TRUE)[seq_len(BATCH_SIZE)]
        }
      }
      if (is.null(best_idx))
        best_idx <- sample(seq_along(current_unsampled), BATCH_SIZE)

      best_global       <- current_unsampled[best_idx]
      current_sampled   <- c(current_sampled, best_global)
      current_unsampled <- current_unsampled[-best_idx]

      cur_r <- sum(colSums(comm_data_pa[current_sampled, , drop = FALSE]) > 0)
      curve_rows[[b + 1L]] <- data.frame(
        Wr = Wr, Wu = Wu, Wunc = Wunc, sim = sim,
        n_samples = N_INITIAL_SAMPLES + b * BATCH_SIZE,
        Richness  = cur_r
      )
    }

    sim_rows[[sim]] <- do.call(rbind, curve_rows)
  } # end sim loop

  message("Weight ", w, "/", nrow(weight_grid),
          "  Wr=", Wr, "  Wu=", Wu, "  Wunc=", Wunc, "  done")

  do.call(rbind, sim_rows)
} # end foreach

stopCluster(cl)
message("Grid search complete. Total rows: ", nrow(results_grid))

# --- Save raw results (optional but recommended for large runs) ---
saveRDS(results_grid, "gridsearch_results.rds")
message("Saved: gridsearch_results.rds")

# ==============================================================================
# Step 4: Aggregate — mean cumulative richness at EVAL_N per weight combination
# ==============================================================================

gridsearch_summary <- results_grid %>%
  dplyr::filter(n_samples == EVAL_N) %>%
  dplyr::group_by(Wr, Wu, Wunc) %>%
  dplyr::summarise(
    mean_richness = mean(Richness),
    sd_richness   = sd(Richness),
    .groups       = "drop"
  )

write.csv(gridsearch_summary, "gridsearch_summary.csv", row.names = FALSE)
message("Saved: gridsearch_summary.csv")

# ==============================================================================
# Step 5: Ternary plot
# Load ggtern HERE (after all ggplot2-based code) to avoid namespace conflicts
# ==============================================================================

if (!requireNamespace("ggtern", quietly = TRUE)) install.packages("ggtern")
library(ggtern)

best_row    <- gridsearch_summary %>% dplyr::slice_max(mean_richness, n = 1)

# Nearest grid approximation to equal weights (1/3, 1/3, 1/3):
# round(1/3, 1) == 0.3 → closest grid point is (0.3, 0.3, 0.4)
current_row <- gridsearch_summary %>%
  dplyr::filter(
    round(Wr,   1) == round(1 / 3, 1),
    round(Wu,   1) == round(1 / 3, 1)
  )

rich_range <- range(gridsearch_summary$mean_richness)

p_tern <- ggtern(gridsearch_summary,
                 aes(x = Wr, y = Wu, z = Wunc, color = mean_richness)) +
  # All weight combinations
  geom_point(size = 3.5, alpha = 0.88) +
  # Best weight combination (red diamond)
  geom_point(
    data        = best_row,
    aes(x = Wr, y = Wu, z = Wunc),
    color       = "red", fill = "red",
    shape       = 23, size = 5,
    inherit.aes = FALSE
  ) +
  # Nearest-grid equal-weight reference (white circle)
  geom_point(
    data        = current_row,
    aes(x = Wr, y = Wu, z = Wunc),
    color       = "white", fill = "white",
    shape       = 21, size = 4,
    inherit.aes = FALSE
  ) +
  scale_color_viridis_c(
    option = "plasma",
    name   = paste0("Mean richness\nat n = ", EVAL_N),
    limits = rich_range
  ) +
  labs(
    title    = "Weight optimization: BAS acquisition function",
    subtitle = paste0(
      "Local dataset (N = 53); n = ", EVAL_N,
      "; 100 simulations per weight combination\n",
      "\u25c6 red = optimal  |  \u25cb white = equal-weight approx. (0.3, 0.3, 0.4)"
    ),
    x = expression(W[rich]),
    y = expression(W[uniq]),
    z = expression(W[unc])
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "right"
  )

ggsave("BAS_WeightOptimization_TernaryPlot.png",
       plot = p_tern, width = 8, height = 7, dpi = 300)
message("Saved: BAS_WeightOptimization_TernaryPlot.png")

# ==============================================================================
# Step 6: Console summary
# ==============================================================================

cat("\n=== Optimal weight combinations (top 5) ===\n")
print(gridsearch_summary %>% dplyr::slice_max(mean_richness, n = 5))

cat("\n=== Current equal-weight approximation (0.3, 0.3, 0.4) ===\n")
print(current_row)

cat("\n=== Improvement over equal-weight approx ===\n")
if (nrow(current_row) > 0 && nrow(best_row) > 0) {
  delta <- best_row$mean_richness[1] - current_row$mean_richness[1]
  cat(sprintf("Best: %.2f  |  Equal-weight approx: %.2f  |  Delta: %+.2f\n",
              best_row$mean_richness[1],
              current_row$mean_richness[1],
              delta))
}

message("\n=== All done ===")
