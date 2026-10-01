# ==============================================================================
# 04: Figure 4 — rare-biosphere comparison, "unique to BAS" vs "unique to Random"
#
# One representative campaign per scale (seed 42), same engines/weights as 01 and 03.
# Taxa detected by one strategy but not the other are compared by their mean
# relative abundance across all sites (ECDF; one-sided Wilcoxon, BAS rarer).
# v3 (2026-10): uses bas_core.R (novelty uniqueness; RK engine locally, QRF globally)
# Output: Fig4AB.png / .pdf and summary statistics on the console
# ==============================================================================
suppressMessages({library(dplyr); library(ggplot2); library(patchwork); library(Matrix)})
source("bas_core.R")

W_LOCAL  <- as.numeric(strsplit(Sys.getenv("BAS_WEIGHTS_LOCAL",  "0.4,0.3,0.3"), ",")[[1]])
W_GLOBAL <- as.numeric(strsplit(Sys.getenv("BAS_WEIGHTS_GLOBAL", "0.2,0.4,0.4"), ",")[[1]])
N_TREES  <- as.integer(Sys.getenv("N_TREES", "200"))

mean_rel_abund <- function(rel) colMeans(rel)          # rel = relative-abundance matrix (rows sum to 1)
sym_df <- function(idx_bas, idx_rnd, mra) {
  ub <- setdiff(idx_bas, idx_rnd); ur <- setdiff(idx_rnd, idx_bas)
  data.frame(Group = c(rep("Unique to BAS", length(ub)), rep("Unique to Random", length(ur))),
             MeanRelAbundance = c(mra[ub], mra[ur])) %>% filter(MeanRelAbundance > 0)
}
sig_star4 <- function(p) ifelse(p < 1e-4, "****", ifelse(p < 1e-3, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
cols <- c("Unique to BAS" = "#0072B2", "Unique to Random" = "#D55E00")
panel <- function(df, tag) {
  p <- wilcox.test(df$MeanRelAbundance[df$Group == "Unique to BAS"], df$MeanRelAbundance[df$Group == "Unique to Random"], alternative = "less")$p.value
  g <- ggplot(df, aes(MeanRelAbundance, color = Group)) + stat_ecdf(linewidth = 1.2) +
    scale_x_log10(labels = scales::label_scientific()) + scale_color_manual(values = cols) +
    labs(tag = tag, x = "Mean Relative Abundance", y = "Cumulative Proportion") +
    annotate("text", x = Inf, y = 0.02, hjust = 1.1, label = paste0("Wilcoxon ", sig_star4(p)), size = 4) +
    theme_minimal(base_size = 14) + theme(legend.position = c(0.72, 0.18), legend.title = element_blank())
  list(plot = g, p = p)
}
stats_of <- function(df) df %>% group_by(Group) %>% summarise(N = n(), Mean = mean(MeanRelAbundance), Median = median(MeanRelAbundance), .groups = "drop")

# ---- Local -------------------------------------------------------------------
md_L <- readRDS("master_data_local.rds"); comm_L <- readRDS("comm_data_local.rds"); pa_L <- readRDS("comm_data_pa_local.rds")
set.seed(42); init_L <- sample(seq_len(nrow(md_L)), 5); rnd_L <- sample(setdiff(seq_len(nrow(md_L)), init_L))
S_bas_L <- run_campaign("BAS",    md_L, pa_L, init_L, rnd_L, batch = 5, n_batches = 5, engine = "rk", w = W_LOCAL)
S_rnd_L <- run_campaign("Random", md_L, pa_L, init_L, rnd_L, batch = 5, n_batches = 5)
mra_L   <- mean_rel_abund(comm_L / rowSums(comm_L))
df_L    <- sym_df(which(colSums(pa_L[S_bas_L, ]) > 0), which(colSums(pa_L[S_rnd_L, ]) > 0), mra_L)
A <- panel(df_L, "A")

# ---- Global ------------------------------------------------------------------
md_G <- readRDS("master_data_GlobalFungi.rds"); pa_G <- readRDS("comm_pa_sp_GlobalFungi.rds"); hel_G <- readRDS("comm_hel_GlobalFungi.rds")
set.seed(42); prob_w <- ifelse(md_G$latitude > 20, 1.0, 0.05)
init_G <- sample(seq_len(nrow(md_G)), 1000, prob = prob_w); rnd_G <- sample(setdiff(seq_len(nrow(md_G)), init_G))
S_bas_G <- run_campaign("BAS",    md_G, pa_G, init_G, rnd_G, batch = 100, n_batches = 10, engine = "qrf", w = W_GLOBAL, n_trees = N_TREES)
S_rnd_G <- run_campaign("Random", md_G, pa_G, init_G, rnd_G, batch = 100, n_batches = 10)
mra_G   <- mean_rel_abund(hel_G^2)                      # Hellinger^2 = relative abundance
df_G    <- sym_df(which(Matrix::colSums(pa_G[S_bas_G, ]) > 0), which(Matrix::colSums(pa_G[S_rnd_G, ]) > 0), mra_G)
B <- panel(df_G, "B")

fig4 <- A$plot + B$plot
ggsave("Fig4AB.png", fig4, width = 10, height = 6, dpi = 300); ggsave("Fig4AB.pdf", fig4, width = 10, height = 6)
cat("\n--- Local (n = 30) ---\n"); print(stats_of(df_L)); cat("P =", signif(A$p, 3), "\n")
cat("\n--- Global (n = 2000) ---\n"); print(stats_of(df_G)); cat("P =", signif(B$p, 3), "\n")
