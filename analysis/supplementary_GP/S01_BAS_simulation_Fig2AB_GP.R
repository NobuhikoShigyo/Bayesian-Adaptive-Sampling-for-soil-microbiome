# ==============================================================================
# BAS (+Env) vs Random: Simulation for Figure 2A-B  — Gaussian Process version
# (Memory-safe / Windows-safe revision)
# ==============================================================================
# Predictor: Universal Kriging via gstat + automap  (kriging = GP regression
# with a parametric covariance; auto-fit Matern/Sph/Exp by REML).
# Acquisition uses GP posterior mean (richness / LCBD) and variance (W_unc).
#
# Windows / 32 GB RAM notes:
#   - Capped parallel workers (PSOCK copies all globals per worker).
#   - USE_PARALLEL = FALSE falls back to sequential — safest on any machine.
#   - Variogram fit failures → nugget-only GP fallback → random fallback.
#
# Output: BAS_Final_Plot_GP.png
# ==============================================================================

# --- Libraries ---
library(vegan)
library(dplyr)
library(gstat)
library(automap)
library(ggplot2)
library(adespatial)
library(sp)
library(tidyr)
library(foreach)
library(doParallel)
library(patchwork)

set.seed(123)

# ==============================================================================
# Step 1: Load and prepare empirical data  (identical to IDW version)
# ==============================================================================

asv_table_raw <- read.csv("comm_data_Shigyo_et_al_2022.csv",
                          row.names = 1, header = TRUE, check.names = FALSE)
site_data     <- read.csv("site_data_Shigyo_et_al_2022.csv", header = TRUE)

common_sites <- intersect(rownames(asv_table_raw), site_data$SiteID)
if (length(common_sites) == 0) stop("Error: no matching SiteIDs.")

asv_table <- asv_table_raw[common_sites, ]
site_data <- site_data %>%
  filter(SiteID %in% common_sites) %>%
  arrange(match(SiteID, common_sites))
asv_table <- asv_table[, colSums(asv_table) > 0]

comm_data    <- as.matrix(asv_table)
comm_data_pa <- (comm_data > 0) * 1
all_richness <- data.frame(SiteID = rownames(comm_data_pa),
                           TrueRichness = rowSums(comm_data_pa))

master_data <- site_data %>%
  left_join(all_richness, by = "SiteID") %>%
  rename(X = lon, Y = lat)

# Drop rows with NA in any quantity used downstream (kriging fails on NAs)
ok_rows <- stats::complete.cases(master_data[, c("X", "Y", "pH", "WC",
                                                 "TrueRichness")])
if (!all(ok_rows)) {
  master_data  <- master_data[ok_rows, , drop = FALSE]
  comm_data    <- comm_data[master_data$SiteID, , drop = FALSE]
  comm_data_pa <- comm_data_pa[master_data$SiteID, , drop = FALSE]
  message("Dropped ", sum(!ok_rows), " rows with NA covariates.")
}

N_SITES_TOTAL <- nrow(master_data)
message("Dataset ready: ", N_SITES_TOTAL, " sites")

saveRDS(master_data,  "master_data_local.rds")
saveRDS(comm_data,    "comm_data_local.rds")
saveRDS(comm_data_pa, "comm_data_pa_local.rds")
message("Saved: master_data / comm_data / comm_data_pa for Fig4.")

rm(asv_table, asv_table_raw, site_data, all_richness); invisible(gc(FALSE))

# ==============================================================================
# Step 2: Simulation settings
# ==============================================================================

N_SIMULATIONS     <- 100
N_INITIAL_SAMPLES <- 5
N_TOTAL_SAMPLES   <- 50
BATCH_SIZE        <- 5
N_BATCHES         <- floor((N_TOTAL_SAMPLES - N_INITIAL_SAMPLES) / BATCH_SIZE)

# --- Parallel settings (Windows-friendly) ---
# Each PSOCK worker copies globals; keep workers modest.
USE_PARALLEL <- TRUE
MAX_CORES    <- 4
num_cores    <- if (USE_PARALLEL) {
  min(MAX_CORES, max(1, parallel::detectCores(logical = FALSE) - 1))
} else 1
message("Parallel workers: ", num_cores)

strategies <- list(
  BAS_Plus  = list(Wr = 1/3, Wu = 1/3, Wunc = 1/3, Name = "BAS",    Type = "Env"),
  Random    = list(                                  Name = "Random", Type = "Random")
)

# ==============================================================================
# Step 3: GP predictor (universal kriging with auto-fit variogram)
# ==============================================================================

normalize01 <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(rep(0, length(x)))
  x[is.na(x)] <- min(x, na.rm = TRUE)
  if (max(x) == min(x)) return(rep(0, length(x)))
  (x - min(x)) / (max(x) - min(x))
}

# Fit GP (universal kriging) on train_sp and predict at cand_sp.
# Returns list(mean, var) or NULL if fitting/prediction fail irrecoverably.
gp_predict <- function(formula, train_sp, cand_sp) {
  vg_model <- tryCatch(
    suppressWarnings(
      automap::autofitVariogram(formula, train_sp,
                                model = c("Sph", "Exp", "Mat"))$var_model
    ),
    error = function(e) NULL
  )
  if (is.null(vg_model)) {
    resp <- all.vars(formula)[1]
    v0   <- stats::var(train_sp@data[[resp]], na.rm = TRUE)
    vg_model <- gstat::vgm(psill  = if (is.finite(v0) && v0 > 0) v0 else 1,
                           model  = "Nug", nugget = 0)
  }

  pred <- tryCatch(
    suppressWarnings(suppressMessages(
      gstat::krige(formula, train_sp, cand_sp, model = vg_model,
                   debug.level = -1)
    )),
    error = function(e) NULL
  )
  if (is.null(pred)) return(NULL)

  v <- if ("var1.var" %in% names(pred)) pred$var1.var else rep(0, length(pred$var1.pred))
  v[is.na(v) | v < 0] <- 0
  list(mean = pred$var1.pred, var = v)
}

# ==============================================================================
# Step 4: Main simulation loop
# ==============================================================================

if (num_cores > 1) {
  cl <- makeCluster(num_cores, type = "PSOCK")
  registerDoParallel(cl)
  `%loop%` <- `%dopar%`
} else {
  registerDoSEQ()
  `%loop%` <- `%do%`
}

results_df <- foreach(
  sim = 1:N_SIMULATIONS,
  .combine       = rbind,
  .errorhandling = "remove",
  .packages      = c("vegan", "dplyr", "gstat", "automap", "adespatial", "sp")
) %loop% {

  set.seed(123 + sim)

  initial_indices     <- sample.int(nrow(master_data), N_INITIAL_SAMPLES)
  remaining_pool_base <- (1:nrow(master_data))[-initial_indices]
  initial_richness    <- sum(colSums(comm_data_pa[initial_indices, , drop = FALSE]) > 0)
  random_pool_order   <- sample(remaining_pool_base)

  sim_results <- list()

  for (strat_key in names(strategies)) {
    conf <- strategies[[strat_key]]

    current_sampled   <- initial_indices
    current_unsampled <- remaining_pool_base

    curve_data <- data.frame(
      n_samples = N_INITIAL_SAMPLES,
      Richness  = initial_richness,
      SimID     = sim,
      Method    = conf$Name
    )

    for (b in 1:N_BATCHES) {

      if (conf$Type == "Random") {
        idx_start <- (b - 1) * BATCH_SIZE + 1
        idx_end   <- b * BATCH_SIZE
        if (length(random_pool_order) >= BATCH_SIZE && idx_end <= length(random_pool_order)) {
          next_idx <- random_pool_order[idx_start:idx_end]
          current_sampled <- c(current_sampled, next_idx)
        }

      } else {
        if (length(current_unsampled) >= BATCH_SIZE) {

          best_idx <- tryCatch({
            curr_comm <- comm_data[current_sampled, , drop = FALSE]
            curr_lcbd <- tryCatch({
              if (nrow(curr_comm) < 3) rep(0, nrow(curr_comm))
              else adespatial::beta.div(vegan::decostand(curr_comm, "hellinger"),
                                        method = "hellinger", nperm = 0)$LCBD
            }, error = function(e) rep(0, length(current_sampled)))

            train_d <- master_data[current_sampled, ] %>% mutate(LCBD = curr_lcbd)
            coordinates(train_d) <- ~X + Y
            cand_d <- master_data[current_unsampled, ]
            coordinates(cand_d) <- ~X + Y

            if (conf$Type == "Env") {
              f_r <- formula(TrueRichness ~ pH + WC)
              f_l <- formula(LCBD          ~ pH + WC)
            } else {
              f_r <- formula(TrueRichness ~ 1)
              f_l <- formula(LCBD          ~ 1)
            }

            gp_r <- gp_predict(f_r, train_d, cand_d)
            gp_l <- gp_predict(f_l, train_d, cand_d)

            if (!is.null(gp_r) && !is.null(gp_l)) {
              sc <- conf$Wr   * normalize01(gp_r$mean) +
                    conf$Wu   * normalize01(gp_l$mean) +
                    conf$Wunc * normalize01(sqrt(gp_r$var))
              order(sc, decreasing = TRUE)[1:BATCH_SIZE]
            } else {
              sample.int(length(current_unsampled), BATCH_SIZE)
            }
          }, error = function(e) sample.int(length(current_unsampled), BATCH_SIZE))

          best_global       <- current_unsampled[best_idx]
          current_sampled   <- c(current_sampled, best_global)
          current_unsampled <- current_unsampled[-best_idx]
        }
      }

      cur_r <- sum(colSums(comm_data_pa[current_sampled, , drop = FALSE]) > 0)
      curve_data <- rbind(curve_data, data.frame(
        n_samples = N_INITIAL_SAMPLES + b * BATCH_SIZE,
        Richness  = cur_r,
        SimID     = sim,
        Method    = conf$Name
      ))
    }
    sim_results[[length(sim_results) + 1]] <- curve_data
  }

  invisible(gc(FALSE))
  do.call(rbind, sim_results)
}

if (num_cores > 1) stopCluster(cl)
message("Simulation complete.")

# ==============================================================================
# Step 5: Figure A — line plot
# ==============================================================================

target_samples_line <- 5:30
target_samples_box  <- c(10, 15, 20, 25, 30)
target_methods      <- c("BAS", "Random")

df_line_summary <- results_df %>%
  filter(n_samples %in% target_samples_line, Method %in% target_methods) %>%
  group_by(n_samples, Method) %>%
  summarise(
    Mean = mean(Richness),
    SD   = sd(Richness),
    ymin = Mean - SD,
    ymax = Mean + SD,
    .groups = "drop"
  )

df_box <- results_df %>%
  filter(n_samples %in% target_samples_box, Method %in% target_methods)

df_line_summary$Method    <- factor(df_line_summary$Method, levels = target_methods)
df_box$Method             <- factor(df_box$Method,          levels = target_methods)
df_box$n_samples_label    <- factor(paste("n =", df_box$n_samples),
                                    levels = paste("n =", target_samples_box))

cols      <- c("BAS" = "#0072B2", "Random" = "#D55E00")
fills     <- c("BAS" = "#0072B2", "Random" = "#D55E00")
linetypes <- c("BAS" = "solid", "Random" = "solid")

p_line <- ggplot(df_line_summary,
                 aes(x = n_samples, y = Mean, color = Method, fill = Method, linetype = Method)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, linetype = 0) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_linetype_manual(values = linetypes) +
  coord_cartesian(xlim = c(5, 30)) +
  labs(tag = "A", x = "Number of Samples", y = "Cumulative Richness") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.margin = margin(b = 20))

# ==============================================================================
# Step 6: Figure B — box plot with significance brackets
# ==============================================================================

get_sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

anno_df <- data.frame()
for (n_lab in levels(df_box$n_samples_label)) {
  dat_n  <- df_box %>% filter(n_samples_label == n_lab)
  y_span <- diff(range(dat_n$Richness))
  gap    <- y_span * 0.1

  paired_dat <- inner_join(
    dat_n %>% filter(Method == "BAS")    %>% select(SimID, BAS = Richness),
    dat_n %>% filter(Method == "Random") %>% select(SimID, Random = Richness),
    by = "SimID"
  )
  if (nrow(paired_dat) > 1) {
    p1 <- t.test(paired_dat$BAS, paired_dat$Random, paired = TRUE, alt = "greater")$p.value
    y1 <- max(dat_n$Richness) + gap

    anno_df <- rbind(anno_df,
      data.frame(n_samples_label = n_lab, x = 1, xend = 2, y = y1, label = get_sig(p1))
    )
  }
}

p_box <- ggplot(df_box, aes(x = Method, y = Richness)) +
  geom_point(aes(color = Method),
             position = position_jitterdodge(0.2, 0, 0.6), size = 0.5, alpha = 0.3) +
  geom_boxplot(aes(color = Method, fill = Method), alpha = 0.6, width = 0.6, outlier.shape = NA) +
  facet_wrap(~n_samples_label, scales = "free_y", nrow = 1) +
  geom_segment(data = anno_df, aes(x = x, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, size = 0.5) +
  geom_segment(data = anno_df, aes(x = x,    xend = x,    y = y, yend = y - (y * 0.002)),
               inherit.aes = FALSE, size = 0.5) +
  geom_segment(data = anno_df, aes(x = xend, xend = xend, y = y, yend = y - (y * 0.002)),
               inherit.aes = FALSE, size = 0.5) +
  geom_text(data = anno_df, aes(x = (x + xend) / 2, y = y, label = label),
            inherit.aes = FALSE, vjust = -0.2, size = 3.5) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = fills) +
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

# ==============================================================================
# Step 7: Combine and save
# ==============================================================================

p_A <- p_line + theme(legend.position = "bottom")
p_B <- p_box  + theme(legend.position = "bottom")

p_combined <- p_A / p_B + plot_layout(heights = c(1, 1.2))

ggsave("BAS_Final_Plot_GP.png", plot = p_combined, width = 10, height = 12, dpi = 300)
message("Saved: BAS_Final_Plot_GP.png")
