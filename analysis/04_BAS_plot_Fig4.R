# ==============================================================================
# 04: Figure 4 — what kind of diversity does BAS find?
#   (A) local  : cumulative number of range-restricted ASVs (present at <= RARE_LOCAL of 53 sites)
#   (B) global : cumulative number of range-restricted SHs  (present at <= RARE_GLOBAL of 5,000 sites)
#   (C) global : share of sampled sites south of 20°N (the region under-represented in the pilot)
#   BAS vs Random vs Oracle, mean ± SD across 100 simulation runs; same pilot sets as 01 / 03.
# v4 (2026-10): replaces the single-campaign abundance ECDF.
# Output: fig4_trajectories.rds, Fig4ABC.png / .pdf
# ==============================================================================
suppressMessages({library(dplyr); library(foreach); library(doParallel); library(Matrix); library(ggplot2); library(patchwork)})
source("bas_core.R")

W_LOCAL  <- as.numeric(strsplit(Sys.getenv("BAS_WEIGHTS_LOCAL",  "0.4,0.2,0.4"), ",")[[1]])
W_GLOBAL <- as.numeric(strsplit(Sys.getenv("BAS_WEIGHTS_GLOBAL", "0.2,0.6,0.2"), ",")[[1]])
N_TREES  <- as.integer(Sys.getenv("N_TREES", "200")); N_SIM <- as.integer(Sys.getenv("N_SIM", "100"))
OFF_L <- as.integer(Sys.getenv("SIM_OFFSET_LOCAL", "50")); OFF_G <- as.integer(Sys.getenv("SIM_OFFSET_GLOBAL", "10"))   # same evaluation pilot sets as 01 / 03
RARE_LOCAL <- 2; RARE_GLOBAL <- 5
cols <- c("Oracle" = "#009E73", "BAS" = "#0072B2", "Random" = "#D55E00"); lts <- c("Oracle" = "dashed", "BAS" = "solid", "Random" = "solid")

if (!file.exists("fig4_trajectories.rds")) {
  # ---- local ---------------------------------------------------------------
  mdL <- readRDS("master_data_local.rds"); paL <- readRDS("comm_data_pa_local.rds"); rareL <- colSums(paL) <= RARE_LOCAL
  cl <- makeCluster(as.integer(Sys.getenv("N_WORKERS", "10")), type = "PSOCK", outfile = ""); registerDoParallel(cl)
  clusterEvalQ(cl, { source("bas_core.R"); library(Matrix) }); clusterExport(cl, c("mdL", "paL", "rareL", "W_LOCAL", "OFF_L"))
  resL <- foreach(sim = OFF_L + (1:N_SIM), .combine = rbind, .packages = "Matrix") %dopar% {
    pl <- draw_pilot_local(sim, nrow(mdL), 5); out <- list()
    for (st in c("BAS", "Random", "Oracle")) {
      rec <- function(S, b) out[[length(out) + 1]] <<- data.frame(SimID = sim, Method = st, n_samples = length(S),
                                                                   Rare = sum(colSums(paL[S, rareL, drop = FALSE]) > 0))
      run_campaign(st, mdL, paL, pl$init, pl$rnd_order, 5, 5, engine = "rk", w = W_LOCAL, record = rec) }
    do.call(rbind, out) }
  # ---- global --------------------------------------------------------------
  mdG <- readRDS("master_data_GlobalFungi.rds"); paG <- readRDS("comm_pa_sp_GlobalFungi.rds")
  rareG <- which(Matrix::colSums(paG) <= RARE_GLOBAL); south <- mdG$latitude < 20
  clusterExport(cl, c("mdG", "paG", "rareG", "south", "W_GLOBAL", "N_TREES", "OFF_G"))
  resG <- foreach(sim = OFF_G + (1:N_SIM), .combine = rbind, .packages = "Matrix") %dopar% {
    pl <- draw_pilot_global(sim, mdG, 1000); out <- list()
    for (st in c("BAS", "Random", "Oracle")) {
      rec <- function(S, b) out[[length(out) + 1]] <<- data.frame(SimID = sim, Method = st, n_samples = length(S),
                                                                   Rare = sum(Matrix::colSums(paG[S, rareG, drop = FALSE]) > 0), South = mean(south[S]))
      run_campaign(st, mdG, paG, pl$init, pl$rnd_order, 100, 20, engine = "qrf", w = W_GLOBAL, n_trees = N_TREES, record = rec) }
    cat(sprintf("sim %d done\n", sim)); do.call(rbind, out) }
  stopCluster(cl)
  saveRDS(list(local = resL, global = resG, rare_local = RARE_LOCAL, rare_global = RARE_GLOBAL,
               n_rare_local = sum(rareL), n_rare_global = length(rareG), south_pool = mean(south)), "fig4_trajectories.rds")
}
tr <- readRDS("fig4_trajectories.rds")

curve_panel <- function(df, yvar, ylab, tag, hline = NULL, pct = FALSE) {
  s <- df %>% group_by(n_samples, Method) %>% summarise(Mean = mean(.data[[yvar]]), SD = sd(.data[[yvar]]), .groups = "drop") %>%
    mutate(Method = factor(Method, c("Oracle", "BAS", "Random")))
  p <- ggplot(s, aes(n_samples, Mean, color = Method, fill = Method, linetype = Method)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.15, linetype = 0) + geom_line(linewidth = 1.2) +
    scale_color_manual(values = cols) + scale_fill_manual(values = cols) + scale_linetype_manual(values = lts) +
    labs(tag = tag, x = "Number of Samples", y = ylab) +
    theme_minimal(base_size = 14) + theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(10, 15, 5, 5))
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, color = "grey40", linetype = "dotted") +
    annotate("text", x = Inf, y = hline, label = "share among all candidates", hjust = 1.05, vjust = -0.4, size = 3.5, color = "grey40")
  if (pct) p <- p + scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  p
}
pA <- curve_panel(tr$local,  "Rare",  sprintf("Range-restricted ASVs found
(present at ≤ %d of 53 sites)", tr$rare_local), "A")
pB <- curve_panel(tr$global, "Rare",  sprintf("Range-restricted SHs found
(present at ≤ %d of 5,000 sites)", tr$rare_global), "B")
pC <- curve_panel(tr$global, "South", "Sampled sites south of 20°N
(share of sites sampled so far)", "C", hline = tr$south_pool, pct = TRUE)
fig4 <- (pA | pB | pC) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("Fig4ABC.png", fig4, width = 15, height = 6, dpi = 300); ggsave("Fig4ABC.pdf", fig4, width = 15, height = 6)

# ---- numbers for the text ----------------------------------------------------
summ <- function(df, yvar, ns) df %>% filter(n_samples %in% ns) %>% group_by(n_samples, Method) %>%
  summarise(m = mean(.data[[yvar]]), sd = sd(.data[[yvar]]), .groups = "drop") %>% tidyr::pivot_wider(names_from = Method, values_from = c(m, sd))
cat("\n--- local: range-restricted ASVs ---\n");  print(as.data.frame(summ(tr$local, "Rare", c(10, 20, 30))), digits = 5)
cat("\n--- global: range-restricted SHs ---\n");  print(as.data.frame(summ(tr$global, "Rare", c(1200, 1500, 2000, 2500, 3000))), digits = 5)
cat("\n--- global: share south of 20N (pool =", round(tr$south_pool, 3), ") ---\n"); print(as.data.frame(summ(tr$global, "South", c(1000, 1200, 1500, 2000, 3000))), digits = 3)
pt <- function(df, yvar, n) { w <- df %>% filter(n_samples == n) %>% select(SimID, Method, all_of(yvar)) %>% tidyr::pivot_wider(names_from = Method, values_from = all_of(yvar))
  t.test(w$BAS, w$Random, paired = TRUE, alternative = "greater")$p.value }
cat("\nP (BAS > Random, paired one-sided): local rare n=20:", signif(pt(tr$local, "Rare", 20), 3), " global rare n=1500:", signif(pt(tr$global, "Rare", 1500), 3), " n=3000:", signif(pt(tr$global, "Rare", 3000), 3), " south n=1500:", signif(pt(tr$global, "South", 1500), 3), "\n")
