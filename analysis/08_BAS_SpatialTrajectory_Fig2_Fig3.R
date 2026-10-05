# ==============================================================================
# BAS Spatial Trajectory Visualization
# Fig2: Local dataset (N=53), snapshots at n = 5, 10, 20, 30
# Fig3: GlobalFungi (N=5000), snapshots at n = 1000, 1500, 2000, 3000
#
# Output files:
#   BAS_SpatialTrajectory_Fig2.png
#   BAS_SpatialTrajectory_Fig2.gif  (if gganimate + gifski installed)
#   BAS_SpatialTrajectory_Fig3.png  (if GlobalFungi RDS files exist)
#   BAS_SpatialTrajectory_Fig3.gif  (if above + gganimate + gifski)
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))   # RStudio
# setwd("/path/to/BAS")                                       # script/terminal

library(vegan)
library(dplyr)
library(gstat)
library(ggplot2)
library(adespatial)
library(sp)
library(sf)
library(patchwork)
library(Matrix)
source("bas_core.R")     # shared engine (v3: regression kriging / QRF, novelty uniqueness)

set.seed(42)

# ==============================================================================
# Shared helpers
# ==============================================================================

# Build a data frame that tags every site with its selection step.
# traj_list: named list with elements $n, $sampled (cumulative), $new (added this step)
build_panel_df <- function(master_df, traj_list, step_labels,
                           lon_col = "X", lat_col = "Y") {
  panels <- lapply(seq_along(traj_list), function(s) {
    master_df %>%
      mutate(
        site_idx = seq_len(n()),
        Status   = case_when(
          site_idx %in% traj_list[[s]]$new     ~ "New this step",
          site_idx %in% traj_list[[s]]$sampled ~ "Previously selected",
          TRUE                                  ~ "Unsampled"
        ),
        Panel = factor(step_labels[s], levels = step_labels)
      )
  })
  bind_rows(panels)
}

# Acquire a world border layer (tries rnaturalearth first, falls back to maps)
get_world_layer <- function() {
  if (requireNamespace("rnaturalearth", quietly = TRUE) &&
      requireNamespace("rnaturalearthdata", quietly = TRUE)) {
    rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  } else if (requireNamespace("maps", quietly = TRUE)) {
    maps::map("world", plot = FALSE, fill = TRUE) %>%
      sf::st_as_sf()
  } else {
    NULL
  }
}

# ==============================================================================
# SECTION A — Fig2: Local dataset (N = 53, Shigyo et al. 2022)
# ==============================================================================
message("\n=== SECTION A: Fig2 Local Dataset ===")

# --- Load data ----------------------------------------------------------------
asv_raw  <- read.csv("comm_data_Shigyo_et_al_2022.csv",
                     row.names = 1, header = TRUE, check.names = FALSE)
site_raw <- read.csv("site_data_Shigyo_et_al_2022.csv", header = TRUE)

common_sites <- intersect(rownames(asv_raw), site_raw$SiteID)
asv_mat      <- asv_raw[common_sites, ]
site_f2      <- site_raw %>%
  filter(SiteID %in% common_sites) %>%
  arrange(match(SiteID, common_sites))
asv_mat      <- asv_mat[, colSums(asv_mat) > 0]

comm_pa_f2   <- (as.matrix(asv_mat) > 0) * 1L
richness_f2  <- data.frame(SiteID      = rownames(comm_pa_f2),
                            TrueRichness = rowSums(comm_pa_f2))

master_f2 <- site_f2 %>%
  left_join(richness_f2, by = "SiteID") %>%
  rename(X = lon, Y = lat)          # X = longitude, Y = latitude

N_F2 <- nrow(master_f2)
message("N = ", N_F2, " sites loaded")

# --- Load shapefiles for Fig2 map background ----------------------------------
# contour1.shp  : elevation contours  (JGD2011 Zone IX, EPSG:6677)
# boundary_polygon.shp: outer boundary polygon (UTM Zone 54N,  EPSG:32654)
shp_boundary <- sf::st_read("boundary_polygon.shp", quiet = TRUE) |>
  sf::st_set_crs(32654) |> sf::st_transform(4326)
shp_contour  <- sf::st_read("contour1.shp",    quiet = TRUE) |>
  sf::st_set_crs(6677)  |> sf::st_transform(4326) |>
  sf::st_intersection(shp_boundary)   # clip to study-area boundary

# --- BAS trajectory for Fig2 (same hyper-parameters as 01) ----------------------
BATCH_F2       <- 5
N_INIT_F2      <- 5
SNAPSHOTS_F2   <- seq(N_INIT_F2, 30, by = BATCH_F2)  # every batch: 5,10,15,20,25,30
PANEL_STEPS_F2 <- c(5, 10, 20, 30)                    # milestone steps for static plot
W_F2 <- as.numeric(strsplit(Sys.getenv("BAS_WEIGHTS_LOCAL", "0.4,0.2,0.4"), ",")[[1]])   # same weights as 01

set.seed(42)
init_idx_f2 <- sample(seq_len(N_F2), N_INIT_F2)
rnd_f2      <- sample(setdiff(seq_len(N_F2), init_idx_f2))

# run one BAS campaign with the shared engine (bas_core.R) and record every batch
traj_f2 <- list(); prev_f2 <- integer(0)
rec_f2  <- function(S, b) { traj_f2[[b + 1L]] <<- list(n = length(S), sampled = S, new = setdiff(S, prev_f2)); prev_f2 <<- S }
invisible(run_campaign("BAS", master_f2, comm_pa_f2, init_idx_f2, rnd_f2, batch = BATCH_F2,
                       n_batches = (max(SNAPSHOTS_F2) - N_INIT_F2) / BATCH_F2, engine = "rk", w = W_F2, record = rec_f2))

panel_idx_f2        <- which(SNAPSHOTS_F2 %in% PANEL_STEPS_F2)
step_labels_f2      <- paste0("n = ", SNAPSHOTS_F2[panel_idx_f2])  # static 4-panel
step_labels_anim_f2 <- paste0("n = ", SNAPSHOTS_F2)                # GIF all steps
df_panels_f2        <- build_panel_df(master_f2, traj_f2[panel_idx_f2], step_labels_f2)

# --- 4-panel map with shapefile background ------------------------------------
p_fig2 <- ggplot() +
  # Outer boundary polygon
  geom_sf(data = shp_boundary, fill = "grey96", color = "black",
          linewidth = 0.7, inherit.aes = FALSE) +
  # Elevation contours
  geom_sf(data = shp_contour, color = "grey55", linewidth = 0.25,
          inherit.aes = FALSE) +
  # Unsampled sites — solid black
  geom_point(
    data   = df_panels_f2 %>% filter(Status == "Unsampled"),
    aes(x = X, y = Y),
    fill = "black", shape = 21, size = 1.8, color = "black", stroke = 0.2, alpha = 0.6
  ) +
  # Previously selected
  geom_point(
    data   = df_panels_f2 %>% filter(Status == "Previously selected"),
    aes(x = X, y = Y),
    fill = "#56B4E9", shape = 21, size = 3.0, color = "white", stroke = 0.6
  ) +
  # Newly selected this step — emphasised with diamond
  geom_point(
    data   = df_panels_f2 %>% filter(Status == "New this step"),
    aes(x = X, y = Y),
    fill = "#D55E00", shape = 23, size = 4.5, color = "white", stroke = 0.8
  ) +
  facet_wrap(~Panel, nrow = 1) +
  coord_sf(crs = 4326) +
  labs(
    title    = "BAS Spatial Trajectory \u2014 Local Dataset (N = 53)",
    subtitle = paste0(
      "Orange \u25c6: newly selected  |  Blue \u25cf: previously selected  ",
      "|  Black \u25cf: unsampled"
    ),
    x = "Longitude (\u00b0E)",
    y = "Latitude (\u00b0N)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    strip.text       = element_text(face = "bold", size = 13),
    panel.border     = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )

ggsave("BAS_SpatialTrajectory_Fig2.png", p_fig2,
       width = 15, height = 5.5, dpi = 300)
message("Saved: BAS_SpatialTrajectory_Fig2.png")

# --- Optional GIF animation for Fig2 ------------------------------------------
if (requireNamespace("gganimate", quietly = TRUE) &&
    requireNamespace("gifski",    quietly = TRUE)) {
  library(gganimate)
  library(gifski)

  # --- Unified data frame: unsampled = black, group = site_idx (no flying) ---
  df_anim_f2 <- bind_rows(lapply(seq_along(SNAPSHOTS_F2), function(s) {
    master_f2 %>%
      mutate(
        site_idx    = seq_len(n()),
        Status      = case_when(
          site_idx %in% traj_f2[[s]]$new     ~ "New this step",
          site_idx %in% traj_f2[[s]]$sampled ~ "Previously selected",
          TRUE                                ~ "Unsampled"
        ),
        Panel       = factor(step_labels_anim_f2[s], levels = step_labels_anim_f2),
        point_fill  = case_when(
          Status == "New this step"       ~ "#D55E00",
          Status == "Previously selected" ~ "#56B4E9",
          TRUE                            ~ "black"
        ),
        point_size  = case_when(
          Status == "New this step"       ~ 5.0,
          Status == "Previously selected" ~ 3.2,
          TRUE                            ~ 1.8
        ),
        point_alpha = ifelse(Status == "Unsampled", 0.6, 1.0)
      )
  }))

  p_anim_f2 <- ggplot(df_anim_f2,
                       aes(x = X, y = Y,
                           fill  = I(point_fill),
                           size  = I(point_size),
                           alpha = I(point_alpha),
                           group = site_idx)) +
    geom_sf(data = shp_boundary, fill = "grey96", color = "black",
            linewidth = 0.7, inherit.aes = FALSE) +
    geom_sf(data = shp_contour, color = "grey55", linewidth = 0.25,
            inherit.aes = FALSE) +
    geom_point(shape = 21, color = "white", stroke = 0.5) +
    coord_sf(crs = 4326) +
    labs(
      title    = "BAS Spatial Trajectory  \u2014  {current_frame}",
      subtitle = "Model selects sites with high predicted richness, novelty, and uncertainty",
      x = "Longitude (\u00b0E)", y = "Latitude (\u00b0N)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      panel.border  = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      legend.position = "none",
      aspect.ratio  = 1
    ) +
    transition_manual(Panel)

  # 6 steps × ~12 frames/step pause = 72 frames; fps=10 → ~7 sec
  anim_save("BAS_SpatialTrajectory_Fig2.gif",
            animate(p_anim_f2,
                    nframes = length(SNAPSHOTS_F2) * 12L,
                    fps = 10,
                    width = 720, height = 620,
                    renderer = gifski_renderer()))
  message("Saved: BAS_SpatialTrajectory_Fig2.gif")
} else {
  message("Tip: install 'gganimate' + 'gifski' to also generate an animated GIF.")
}

# ==============================================================================
# SECTION B — Fig3: GlobalFungi dataset (N = 5,000)
# Requires: master_data_GlobalFungi.rds, comm_pa_GlobalFungi.rds,
#           comm_hel_GlobalFungi.rds  (included in the repository; see README "Data")
# ==============================================================================
message("\n=== SECTION B: Fig3 GlobalFungi Dataset ===")

rds_needed <- c("master_data_GlobalFungi.rds",
                "comm_pa_GlobalFungi.rds",
                "comm_hel_GlobalFungi.rds")

if (!all(file.exists(rds_needed))) {
  message("One or more GlobalFungi RDS files not found:")
  message(paste(" ", rds_needed[!file.exists(rds_needed)], collapse = "\n"))
  message("These files are part of the repository (see README \"Data\").")
  message("Skipping Section B.")
} else {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("Package 'ranger' is required for Section B. Install with install.packages('ranger').")
  library(ranger)

  master_f3  <- readRDS("master_data_GlobalFungi.rds")
  comm_pa_f3 <- readRDS("comm_pa_GlobalFungi.rds")
  comm_hel_f3 <- readRDS("comm_hel_GlobalFungi.rds")

  N_F3 <- nrow(master_f3)
  message("N = ", N_F3, " sites loaded")

  # --- BAS trajectory for Fig3 (one run; follows Fig3 simulation settings) ---
  BATCH_F3_SIZE  <- 100
  N_INIT_F3      <- 1000
  SNAPSHOTS_F3   <- seq(N_INIT_F3, 3000, by = BATCH_F3_SIZE)  # every batch: 1000,1100,...,3000
  PANEL_STEPS_F3 <- c(1000, 1500, 2000, 3000)                  # milestone steps for static plot
  W_F3 <- as.numeric(strsplit(Sys.getenv("BAS_WEIGHTS_GLOBAL", "0.2,0.6,0.2"), ",")[[1]])   # same weights as 03

  set.seed(42)
  prob_w    <- ifelse(master_f3$latitude > 20, 1.0, 0.05)
  init_f3   <- sample(seq_len(N_F3), N_INIT_F3, prob = prob_w)
  rnd_f3    <- sample(setdiff(seq_len(N_F3), init_f3))

  message("Running BAS trajectory for Fig3 (this may take several minutes)...")
  traj_f3 <- list(); prev_f3 <- integer(0)
  rec_f3  <- function(S, b) { traj_f3[[b + 1L]] <<- list(n = length(S), sampled = S, new = setdiff(S, prev_f3)); prev_f3 <<- S
                              message("  Batch n = ", length(S), " recorded.") }
  invisible(run_campaign("BAS", master_f3, comm_pa_f3, init_f3, rnd_f3, batch = BATCH_F3_SIZE,
                         n_batches = (max(SNAPSHOTS_F3) - N_INIT_F3) / BATCH_F3_SIZE, engine = "qrf", w = W_F3, n_trees = 200, record = rec_f3))

  panel_idx_f3        <- which(SNAPSHOTS_F3 %in% PANEL_STEPS_F3)
  step_labels_f3      <- paste0("n = ", format(SNAPSHOTS_F3[panel_idx_f3], big.mark = ","))
  step_labels_anim_f3 <- paste0("n = ", format(SNAPSHOTS_F3, big.mark = ","))
  df_panels_f3 <- build_panel_df(master_f3, traj_f3[panel_idx_f3], step_labels_f3,
                                  lon_col = "longitude", lat_col = "latitude") %>%
    mutate(lon = master_f3$longitude[site_idx],
           lat = master_f3$latitude[site_idx])

  # --- 4-panel world map -------------------------------------------------------
  world <- get_world_layer()

  base_map <- if (!is.null(world)) {
    geom_sf(data = world, fill = "grey92", color = "grey65", linewidth = 0.15,
            inherit.aes = FALSE)
  } else {
    geom_blank()   # plain axes as fallback
  }

  p_fig3 <- ggplot() +
    base_map +
    # Unsampled — solid black
    geom_point(
      data   = df_panels_f3 %>% filter(Status == "Unsampled"),
      aes(x = lon, y = lat),
      fill = "black", shape = 21, size = 0.9, color = "black", stroke = 0.1, alpha = 0.45
    ) +
    # Previously selected
    geom_point(
      data   = df_panels_f3 %>% filter(Status == "Previously selected"),
      aes(x = lon, y = lat),
      fill = "#56B4E9", shape = 21, size = 1.6, color = "white", stroke = 0.25
    ) +
    # Newly selected this step
    geom_point(
      data   = df_panels_f3 %>% filter(Status == "New this step"),
      aes(x = lon, y = lat),
      fill = "#D55E00", shape = 23, size = 2.6, color = "white", stroke = 0.5
    ) +
    facet_wrap(~Panel, nrow = 2) +
    coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
    labs(
      title    = "BAS Spatial Trajectory \u2014 GlobalFungi Dataset (N = 5,000)",
      subtitle = paste0(
        "Orange \u25c6: newly selected  |  Blue \u25cf: previously selected  ",
        "|  Black \u25cf: unsampled"
      ),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(size = 9, color = "grey40"),
      strip.text      = element_text(face = "bold", size = 12),
      panel.border    = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      legend.position = "right",
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      panel.grid      = element_blank()
    )

  ggsave("BAS_SpatialTrajectory_Fig3.png", p_fig3,
         width = 14, height = 9, dpi = 300)
  message("Saved: BAS_SpatialTrajectory_Fig3.png")

  # --- Optional GIF animation for Fig3 ----------------------------------------
  if (requireNamespace("gganimate", quietly = TRUE) &&
      requireNamespace("gifski",    quietly = TRUE)) {
    library(gganimate)
    library(gifski)

    # --- Build from ALL batch steps (not just milestone panels) ---
    df_anim_f3 <- bind_rows(lapply(seq_along(SNAPSHOTS_F3), function(s) {
      master_f3 %>%
        mutate(
          site_idx    = seq_len(n()),
          Status      = case_when(
            site_idx %in% traj_f3[[s]]$new     ~ "New this step",
            site_idx %in% traj_f3[[s]]$sampled ~ "Previously selected",
            TRUE                                ~ "Unsampled"
          ),
          Panel       = factor(step_labels_anim_f3[s], levels = step_labels_anim_f3),
          lon         = longitude,
          lat         = latitude,
          point_fill  = case_when(
            Status == "New this step"       ~ "#D55E00",
            Status == "Previously selected" ~ "#56B4E9",
            TRUE                            ~ "black"
          ),
          point_size  = case_when(
            Status == "New this step"       ~ 2.8,
            Status == "Previously selected" ~ 1.6,
            TRUE                            ~ 0.9
          ),
          point_alpha = ifelse(Status == "Unsampled", 0.45, 1.0)
        )
    }))

    p_anim_f3 <- ggplot(df_anim_f3,
                         aes(x = lon, y = lat,
                             fill  = I(point_fill),
                             size  = I(point_size),
                             alpha = I(point_alpha),
                             group = site_idx)) +
      base_map +
      geom_point(shape = 21, color = "white", stroke = 0.3) +
      coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
      labs(
        title    = "BAS Spatial Trajectory  \u2014  {current_frame}",
        subtitle = "BAS targets globally underexplored, high-richness, high-uncertainty regions",
        x = NULL, y = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        panel.border  = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
        axis.text     = element_blank(),
        panel.grid    = element_blank()
      ) +
      transition_manual(Panel)

    # 21 steps × ~8 frames/step pause = 168 frames total; fps=12 → ~14 sec
    anim_save("BAS_SpatialTrajectory_Fig3.gif",
              animate(p_anim_f3,
                      nframes = length(SNAPSHOTS_F3) * 8L,
                      fps = 12,
                      width = 900, height = 540,
                      renderer = gifski_renderer()))
    message("Saved: BAS_SpatialTrajectory_Fig3.gif")
  } else {
    message("Tip: install 'gganimate' + 'gifski' to also generate an animated GIF.")
  }
}

message("\n=== All done ===")
