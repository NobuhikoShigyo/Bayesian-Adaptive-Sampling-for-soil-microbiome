# ==============================================================================
# bas_core.R — shared engine for the BAS simulations (v3, 2026-10)
#
# Acquisition:  a(x) = W_rich * mu_rich(x) + W_uniq * mu_uniq(x) + W_unc * sigma(x)
#   mu_rich : predicted richness
#   mu_uniq : predicted "novelty" — number of taxa present at a site and absent
#             from every other site sampled so far (presence/absence counterpart
#             of LCBD, in richness units; recomputed each batch from observed data)
#   sigma   : prediction uncertainty of the richness model
# Engines:  "rk"  = regression kriging (automap-fitted variogram, gstat::krige)   — Case 1
#           "qrf" = Quantile Regression Forest (ranger, jackknife SE)             — Case 2
# ==============================================================================

norm01 <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(rep(0, length(x)))
  x[is.na(x)] <- min(x, na.rm = TRUE)
  if (max(x) == min(x)) return(rep(0, length(x)))
  (x - min(x)) / (max(x) - min(x))
}

# novelty (unique-taxa count) for sampled sites S; pa = dense 0/1 matrix or sparse Matrix
novelty_of <- function(pa, S) {
  sub <- pa[S, , drop = FALSE]
  singles <- Matrix::colSums(sub) == 1
  as.numeric(Matrix::rowSums(sub[, singles, drop = FALSE]))
}
richness_of <- function(pa, S) sum(Matrix::colSums(pa[S, , drop = FALSE]) > 0)

# ---- Engine: regression kriging with automap variogram -----------------------
# master_data must contain X, Y and the covariates named in `covs`
# NOTE: automap::autofitVariogram can crash the R session with very few training points
# (observed with 5 points); below MIN_AUTOFIT sites a nugget-only model (the covariate regression) is used.
MIN_AUTOFIT <- 10
SMALL_N_RULE <- "fixed_exp"   # "fixed_exp": exponential variogram with psill = var(y), range = median pairwise distance; "nugget": covariate regression only
engine_rk <- function(master_data, S, U, y_rich, y_uniq, covs = c("pH", "WC")) {
  tr <- master_data[S, ]; tr$RICH <- y_rich; tr$UNIQ <- y_uniq
  cd <- master_data[U, ]
  sp::coordinates(tr) <- ~X + Y; sp::coordinates(cd) <- ~X + Y
  rhs <- paste(covs, collapse = " + ")
  one <- function(resp) {
    fm <- stats::as.formula(paste(resp, "~", rhs))
    vg <- if (length(S) >= MIN_AUTOFIT) tryCatch(suppressWarnings(automap::autofitVariogram(fm, tr, model = c("Sph", "Exp", "Mat"))$var_model),
                                                  error = function(e) NULL) else NULL
    if (is.null(vg)) { v0 <- stats::var(tr@data[[resp]], na.rm = TRUE); v0 <- if (is.finite(v0) && v0 > 0) v0 else 1
      if (length(S) < MIN_AUTOFIT && SMALL_N_RULE == "fixed_exp") {
        d <- as.matrix(stats::dist(sp::coordinates(tr))); rng <- stats::median(d[upper.tri(d)]); if (!is.finite(rng) || rng <= 0) rng <- 1
        vg <- gstat::vgm(psill = v0, model = "Exp", range = rng, nugget = 0)
      } else vg <- gstat::vgm(psill = v0, model = "Nug", nugget = 0) }
    p <- tryCatch(suppressWarnings(suppressMessages(gstat::krige(fm, tr, cd, model = vg, debug.level = -1))),
                  error = function(e) NULL)
    if (is.null(p)) return(NULL)
    v <- if ("var1.var" %in% names(p)) p$var1.var else rep(0, length(p$var1.pred)); v[is.na(v) | v < 0] <- 0
    list(mean = p$var1.pred, sd = sqrt(v))
  }
  r <- one("RICH"); u <- one("UNIQ")
  if (is.null(r) || is.null(u)) return(NULL)
  list(rich = r$mean, uniq = u$mean, unc = r$sd)
}

# ---- Engine: Quantile Regression Forest --------------------------------------
engine_qrf <- function(master_data, S, U, y_rich, y_uniq, covs = c("pH", "MAT", "MAP", "SOC", "X", "Y"), n_trees = 200) {
  tr <- master_data[S, covs, drop = FALSE]; cd <- master_data[U, covs, drop = FALSE]
  one <- function(y, se) {
    d <- tr; d$y <- y
    m <- tryCatch(ranger::ranger(y ~ ., data = d, num.trees = n_trees, quantreg = TRUE, keep.inbag = TRUE,
                                 num.threads = 1, save.memory = TRUE), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    if (se) { p <- predict(m, data = cd, type = "se"); list(mean = p$predictions, sd = p$se) }
    else list(mean = predict(m, data = cd)$predictions, sd = NULL)
  }
  r <- one(y_rich, TRUE); u <- one(y_uniq, FALSE)
  if (is.null(r) || is.null(u)) return(NULL)
  list(rich = r$mean, uniq = u$mean, unc = r$sd)
}

acquisition <- function(p, w) w[1] * norm01(p$rich) + w[2] * norm01(p$uniq) + w[3] * norm01(p$unc)

# ---- One campaign ------------------------------------------------------------
# strategy: "BAS" (uses engine + weights), "Random" (rnd_order), "Oracle" (greedy, full knowledge)
# returns integer vector of the final sampled indices in order of acquisition
run_campaign <- function(strategy, master_data, pa, init, rnd_order, batch, n_batches,
                         engine = c("rk", "qrf"), w = c(1/3, 1/3, 1/3), covs = NULL, n_trees = 200,
                         record = NULL) {
  engine <- match.arg(engine)
  S <- init; U <- setdiff(seq_len(nrow(master_data)), init)
  if (!is.null(record)) record(S, 0)
  for (b in seq_len(n_batches)) {
    if (length(U) < 1) break
    k <- min(batch, length(U))
    if (strategy == "Random") {
      nxt <- rnd_order[((b - 1) * batch + 1):((b - 1) * batch + k)]
      S <- c(S, nxt); U <- setdiff(U, nxt)
    } else if (strategy == "Oracle") {
      found <- Matrix::colSums(pa[S, , drop = FALSE]) > 0
      new_sp <- Matrix::rowSums(pa[U, !found, drop = FALSE])
      j <- order(new_sp, decreasing = TRUE)[seq_len(k)]
      S <- c(S, U[j]); U <- U[-j]
    } else {
      y_rich <- master_data$TrueRichness[S]; y_uniq <- novelty_of(pa, S)
      p <- if (engine == "rk") engine_rk(master_data, S, U, y_rich, y_uniq, covs = if (is.null(covs)) c("pH", "WC") else covs)
           else engine_qrf(master_data, S, U, y_rich, y_uniq, covs = if (is.null(covs)) c("pH", "MAT", "MAP", "SOC", "X", "Y") else covs, n_trees = n_trees)
      j <- if (is.null(p)) sample(length(U), k) else order(acquisition(p, w), decreasing = TRUE)[seq_len(k)]
      S <- c(S, U[j]); U <- U[-j]
    }
    if (!is.null(record)) record(S, b)
  }
  S
}

# cumulative-richness curve for a campaign (data.frame rows per batch)
campaign_curve <- function(strategy, label, sim, master_data, pa, init, rnd_order, batch, n_batches, ...) {
  rows <- list()
  rec <- function(S, b) rows[[length(rows) + 1]] <<- data.frame(SimID = sim, Method = label, n_samples = length(S),
                                                                 Richness = richness_of(pa, S))
  run_campaign(strategy, master_data, pa, init, rnd_order, batch, n_batches, record = rec, ...)
  do.call(rbind, rows)
}

# pilot draw identical to 01 (local) and 03 (global)
draw_pilot_local  <- function(sim, n_sites, n_init) { set.seed(123 + sim); init <- sample(seq_len(n_sites), n_init)
  list(init = init, rnd_order = sample(setdiff(seq_len(n_sites), init))) }
draw_pilot_global <- function(sim, master_data, n_init) { set.seed(123 + sim)
  prob_w <- ifelse(master_data$latitude > 20, 1.0, 0.05)
  init <- sample(seq_len(nrow(master_data)), n_init, prob = prob_w)
  list(init = init, rnd_order = sample(setdiff(seq_len(nrow(master_data)), init))) }

weight_grid <- function(step = 0.1) { g <- expand.grid(Wr = seq(0, 1, step), Wu = seq(0, 1, step), Wunc = seq(0, 1, step))
  g[abs(g$Wr + g$Wu + g$Wunc - 1) < 1e-9, ] }
