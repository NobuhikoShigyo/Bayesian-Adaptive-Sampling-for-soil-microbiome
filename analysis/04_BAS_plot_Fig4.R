# ==============================================================================
# Figure 4: Symmetric rare biosphere comparison
#   "Unique to BAS" vs "Unique to Random"
#
# Run AFTER BAS_simulation_Fig2AB.R and BAS_simulation_Fig3AB.R.
# Required session objects from Fig3: master_data, comm_pa, comm_hel, norm01
# Local dataset (comm_data / master_data overwritten by Fig3) is reloaded
# from the original CSVs at the top of this script.
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

# ==============================================================================
# Setup: reload local dataset
# (master_data and comm_data are overwritten when Fig3 runs)
# ==============================================================================

if (file.exists("master_data_local.rds")) {
  master_data_L <- readRDS("master_data_local.rds")
  comm_data_L   <- readRDS("comm_data_local.rds")
  comm_pa_L     <- readRDS("comm_data_pa_local.rds")
} else {
  site_data_local <- read.csv("site_data_Shigyo_et_al_2022.csv", header = TRUE)
  asv_raw_L <- as.matrix(read.csv("comm_data_Shigyo_et_al_2022.csv",
                                   row.names = 1, header = TRUE, check.names = FALSE))
  common_L     <- intersect(rownames(asv_raw_L), site_data_local$SiteID)
  comm_data_L  <- asv_raw_L[common_L, , drop = FALSE]
  comm_data_L  <- comm_data_L[, colSums(comm_data_L) > 0, drop = FALSE]
  comm_pa_L    <- (comm_data_L > 0) * 1
  richness_L   <- data.frame(SiteID = rownames(comm_pa_L),
                              TrueRichness = rowSums(comm_pa_L))
  master_data_L <- site_data_local %>%
    filter(SiteID %in% common_L) %>%
    arrange(match(SiteID, common_L)) %>%
    left_join(richness_L, by = "SiteID") %>%
    rename(X = lon, Y = lat)
}

message("Local dataset ready: ", nrow(master_data_L), " sites, ",
        ncol(comm_data_L), " OTUs")

# ==============================================================================
# Helpers
# ==============================================================================

# Mean relative abundance per species across all sites
calc_mean_rel_abund <- function(mat) {
  row_sums <- rowSums(mat)
  valid    <- row_sums > 0
  rel      <- mat[valid, , drop = FALSE] / row_sums[valid]
  colMeans(rel, na.rm = TRUE)
}

# Symmetric group comparison data frame
make_sym_df <- function(idx_bas, idx_rnd, mean_abund) {
  idx_unique_bas <- setdiff(idx_bas, idx_rnd)
  idx_unique_rnd <- setdiff(idx_rnd, idx_bas)
  data.frame(
    Group            = c(rep("Unique to BAS",    length(idx_unique_bas)),
                         rep("Unique to Random", length(idx_unique_rnd))),
    MeanRelAbundance = c(mean_abund[idx_unique_bas],
                         mean_abund[idx_unique_rnd])
  )
}

# Significance stars from p-value
get_sig <- function(p) {
  ifelse(p < 0.0001, "****",
  ifelse(p < 0.001,  "***",
  ifelse(p < 0.01,   "**",
  ifelse(p < 0.05,   "*", "ns"))))
}

# Normalize to [0, 1]; NAs replaced by column minimum before scaling
norm01_L <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(rep(0, length(x)))
  x[is.na(x)] <- min(x, na.rm = TRUE)
  if (max(x) == min(x)) return(rep(0, length(x)))
  (x - min(x)) / (max(x) - min(x))
}

# Shared plot theme and aesthetics
col_vals <- c("Unique to BAS" = "#0072B2", "Unique to Random" = "#D55E00")
y_expand <- expansion(mult = c(0.05, 0.25))

p_theme <- theme_minimal(base_size = 14) +
  theme(legend.position  = "bottom",
        legend.title      = element_blank(),
        plot.margin       = margin(t = 20, r = 10, b = 10, l = 10))

# ==============================================================================
# Panel A  —  Local scale
# ==============================================================================

message("Panel A: running one representative simulation on local dataset...")

TARGET_N_LOCAL <- 30
N_INIT_LOCAL   <- 5
BATCH_LOCAL    <- 5
set.seed(42)

init_idx_L  <- sample(seq_len(nrow(master_data_L)), N_INIT_LOCAL)
remain_L    <- setdiff(seq_len(nrow(master_data_L)), init_idx_L)
rnd_order_L <- sample(remain_L)

# BAS (+Env) trajectory
curr_bas_L <- init_idx_L
curr_uns_L <- remain_L

while (length(curr_bas_L) < TARGET_N_LOCAL) {

  curr_comm <- comm_data_L[curr_bas_L, , drop = FALSE]
  curr_lcbd <- tryCatch(
    adespatial::beta.div(vegan::decostand(curr_comm, "hellinger"),
                         method = "hellinger", nperm = 0)$LCBD,
    error = function(e) rep(0, nrow(curr_comm))
  )

  train_d <- master_data_L[curr_bas_L, ] %>% mutate(LCBD = curr_lcbd)
  cand_d  <- master_data_L[curr_uns_L, ]
  sp::coordinates(train_d) <- ~X + Y
  sp::coordinates(cand_d)  <- ~X + Y

  g_r <- try(gstat::gstat(formula = TrueRichness ~ pH + WC,
                           locations = train_d, nmax = 10, set = list(idp = 0.5)),
             silent = TRUE)
  g_l <- try(gstat::gstat(formula = LCBD ~ pH + WC,
                           locations = train_d, nmax = 10, set = list(idp = 0.5)),
             silent = TRUE)

  n_pick     <- min(BATCH_LOCAL, TARGET_N_LOCAL - length(curr_bas_L), length(curr_uns_L))
  best_local <- NULL

  if (!inherits(g_r, "try-error") && !inherits(g_l, "try-error")) {
    p_r <- try(suppressMessages(predict(g_r, newdata = cand_d, debug.level = -1)),
               silent = TRUE)
    p_l <- try(suppressMessages(predict(g_l, newdata = cand_d, debug.level = -1)),
               silent = TRUE)
    if (!inherits(p_r, "try-error") && !inherits(p_l, "try-error")) {
      var_r      <- pmax(p_r$var1.var, 0)
      sc         <- (1/3) * norm01_L(p_r$var1.pred) +
                    (1/3) * norm01_L(p_l$var1.pred) +
                    (1/3) * norm01_L(sqrt(var_r))
      best_local <- order(sc, decreasing = TRUE)[seq_len(n_pick)]
    }
  }
  if (is.null(best_local)) best_local <- sample(length(curr_uns_L), n_pick)

  curr_bas_L <- c(curr_bas_L, curr_uns_L[best_local])
  curr_uns_L <- curr_uns_L[-best_local]
}

# Random trajectory
n_needed_L <- TARGET_N_LOCAL - N_INIT_LOCAL
curr_rnd_L <- c(init_idx_L, rnd_order_L[seq_len(n_needed_L)])

# Detected species (column indices)
idx_bas_L <- which(colSums(comm_pa_L[curr_bas_L, , drop = FALSE]) > 0)
idx_rnd_L <- which(colSums(comm_pa_L[curr_rnd_L, , drop = FALSE]) > 0)

mra_L    <- calc_mean_rel_abund(comm_data_L)
df_local <- make_sym_df(idx_bas_L, idx_rnd_L, mra_L) %>%
  filter(MeanRelAbundance > 0)

cat(sprintf("[Local]  Unique to BAS: %d  |  Unique to Random: %d\n",
            sum(df_local$Group == "Unique to BAS"),
            sum(df_local$Group == "Unique to Random")))

# Panel A plot
p_val_A <- wilcox.test(
  df_local$MeanRelAbundance[df_local$Group == "Unique to BAS"],
  df_local$MeanRelAbundance[df_local$Group == "Unique to Random"],
  alternative = "less"
)$p.value

p4_A <- ggplot(df_local, aes(x = MeanRelAbundance, color = Group)) +
  stat_ecdf(linewidth = 1.2) +
  scale_x_log10(labels = scales::label_scientific()) +
  scale_color_manual(values = col_vals) +
  labs(x = "Mean Relative Abundance", y = "Cumulative Proportion") +
  annotate("text", x = Inf, y = 0.02, hjust = 1.1,
           label = paste0("Wilcoxon ", get_sig(p_val_A)), size = 4) +
  p_theme +
  theme(legend.position = c(0.72, 0.18))

message("Panel A complete.")

# ==============================================================================
# Panel B  —  Global scale
# (uses master_data, comm_pa, comm_hel, norm01 from Fig3 session)
# ==============================================================================

message("Panel B: running one representative simulation on GlobalFungi dataset...")

TARGET_N_GLOBAL <- 2000
N_INIT_GLOBAL   <- 1000
BATCH_GLOBAL    <- 100

# Load GlobalFungi session objects from RDS if not available from Fig3 session
if (!exists("master_data") || nrow(master_data) < N_INIT_GLOBAL)
  master_data <- readRDS("master_data_GlobalFungi.rds")
if (!exists("comm_pa"))
  comm_pa  <- readRDS("comm_pa_GlobalFungi.rds")
if (!exists("comm_hel"))
  comm_hel <- readRDS("comm_hel_GlobalFungi.rds")
if (!exists("norm01"))
  norm01 <- function(x) {
    if (length(x[!is.na(x)]) == 0 || max(x, na.rm = TRUE) == min(x, na.rm = TRUE))
      return(rep(0, length(x)))
    out <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    out[is.na(out)] <- 0
    out
  }

set.seed(42)

prob_w      <- ifelse(master_data$latitude > 20, 1.0, 0.05)
init_idx_G  <- sample(seq_len(nrow(master_data)), N_INIT_GLOBAL, prob = prob_w)
remain_G    <- setdiff(seq_len(nrow(master_data)), init_idx_G)
rnd_order_G <- sample(remain_G)

# BAS trajectory
curr_bas_G <- init_idx_G
curr_uns_G <- remain_G

while (length(curr_bas_G) < TARGET_N_GLOBAL &&
       length(curr_uns_G) >= BATCH_GLOBAL) {

  curr_hel  <- comm_hel[curr_bas_G, , drop = FALSE]
  centroid  <- colMeans(curr_hel)
  centered  <- sweep(curr_hel, 2, centroid, "-")
  ss_total  <- sum(rowSums(centered^2))
  curr_lcbd <- if (ss_total > 0) rowSums(centered^2) / ss_total else
               rep(0, length(curr_bas_G))

  train_d      <- master_data[curr_bas_G, ]
  train_d$LCBD <- curr_lcbd
  cand_d       <- master_data[curr_uns_G, ]

  mod_r <- ranger::ranger(TrueRichness ~ pH + MAT + MAP + SOC + X + Y,
                          data = train_d, num.trees = 200,
                          quantreg = TRUE, keep.inbag = TRUE,
                          num.threads = 1, save.memory = TRUE)
  mod_l <- ranger::ranger(LCBD ~ pH + MAT + MAP + SOC + X + Y,
                          data = train_d, num.trees = 200,
                          quantreg = TRUE, keep.inbag = TRUE,
                          num.threads = 1, save.memory = TRUE)

  pred_r_val <- predict(mod_r, data = cand_d)$predictions
  pred_l_val <- predict(mod_l, data = cand_d)$predictions
  pred_r_se  <- predict(mod_r, data = cand_d, type = "se")$se

  sc       <- norm01(pred_r_val) * 0.25 +
              norm01(pred_l_val) * 0.25 +
              norm01(pred_r_se)  * 0.50
  best_idx <- order(sc, decreasing = TRUE)[seq_len(BATCH_GLOBAL)]

  curr_bas_G <- c(curr_bas_G, curr_uns_G[best_idx])
  curr_uns_G <- curr_uns_G[-best_idx]

  cat(sprintf("  BAS samples: %d / %d\n", length(curr_bas_G), TARGET_N_GLOBAL))
  gc()
}

# Random trajectory
n_needed_G <- TARGET_N_GLOBAL - N_INIT_GLOBAL
curr_rnd_G <- c(init_idx_G, rnd_order_G[seq_len(n_needed_G)])

# Detected species (column indices)
idx_bas_G <- which(colSums(comm_pa[curr_bas_G, , drop = FALSE]) > 0)
idx_rnd_G <- which(colSums(comm_pa[curr_rnd_G, , drop = FALSE]) > 0)

mra_G     <- calc_mean_rel_abund(comm_data)
df_global <- make_sym_df(idx_bas_G, idx_rnd_G, mra_G) %>%
  filter(MeanRelAbundance > 0)

cat(sprintf("[Global] Unique to BAS: %d  |  Unique to Random: %d\n",
            sum(df_global$Group == "Unique to BAS"),
            sum(df_global$Group == "Unique to Random")))

# Panel B plot
p_val_B <- wilcox.test(
  df_global$MeanRelAbundance[df_global$Group == "Unique to BAS"],
  df_global$MeanRelAbundance[df_global$Group == "Unique to Random"],
  alternative = "less"
)$p.value

p4_B <- ggplot(df_global, aes(x = MeanRelAbundance, color = Group)) +
  stat_ecdf(linewidth = 1.2) +
  scale_x_log10(labels = scales::label_scientific()) +
  scale_color_manual(values = col_vals) +
  labs(x = "Mean Relative Abundance", y = "Cumulative Proportion") +
  annotate("text", x = Inf, y = 0.02, hjust = 1.1,
           label = paste0("Wilcoxon ", get_sig(p_val_B)), size = 4) +
  p_theme +
  theme(legend.position = c(0.72, 0.18))

message("Panel B complete.")

# ==============================================================================
# Combine and save Figure 4
# ==============================================================================

fig4 <- p4_A + p4_B +
  plot_annotation(tag_levels = "A")

print(fig4)

ggsave("Figure4_RareBiosphere1.png", fig4, width = 10, height = 6, dpi = 300)
ggsave("Figure4_RareBiosphere1.pdf", fig4, width = 10, height = 6)
message("Saved: Figure4_RareBiosphere.png / .pdf")

# ==============================================================================
# Summary statistics
# ==============================================================================

cat("\n--- Local scale (n = 30) ---\n")
print(df_local  %>% group_by(Group) %>%
      summarise(Mean   = mean(MeanRelAbundance),
                Median = median(MeanRelAbundance),
                N      = n(), .groups = "drop"))

cat("\n--- Global scale (n = 2000) ---\n")
print(df_global %>% group_by(Group) %>%
      summarise(Mean   = mean(MeanRelAbundance),
                Median = median(MeanRelAbundance),
                N      = n(), .groups = "drop"))
