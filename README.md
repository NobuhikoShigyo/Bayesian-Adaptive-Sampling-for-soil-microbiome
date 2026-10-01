# Bayesian Adaptive Sampling for Soil Microbiomes

R code, data and figures for the paper
*"We're only in it for the information: A framework for Bayesian adaptive sampling in soil microbiomes."*

A data-driven approach for deciding where to sample soils next.

> [!IMPORTANT]
> **Work in Progress** — this repository accompanies a manuscript in submission/review. Files may still change.

---

## Overview

This project advocates a shift in soil microbiome sampling from static spatial coverage to a dynamic philosophy of *information discovery*. We propose a Bayesian Adaptive Sampling (BAS) framework that reframes sampling as a sequential decision-making process. Each iteration selects the next batch of sites by optimizing a weighted acquisition function

`a(x) = W_rich · μ_rich(x) + W_uniq · μ_uniq(x) + W_unc · σ(x)`

that combines

- predicted **richness**,
- predicted **uniqueness** — the *novelty* of a site, i.e. the number of taxa recorded there that have not been recorded at any other site sampled so far (the presence/absence counterpart of LCBD, recomputed after every batch from observed data only), and
- the prediction **uncertainty** of the richness model.

The surrogate model that produces these quantities is interchangeable. In the paper, regression kriging with an automatically fitted variogram is used at the local scale (Case 1) and Quantile Regression Forests at the global scale (Case 2).

| Strategy | Description |
|---|---|
| **BAS** | Adaptive sampling using environmental covariates + spatial coordinates |
| **Random** | Random (non-adaptive) baseline |
| **Oracle** | Greedy selection with full knowledge of the community table — a theoretical upper bound |

---

## Repository layout

```
BAS_publication/
├── analysis/          # working directory: scripts + every file they read or write
│   ├── bas_core.R     # shared engine: novelty, regression kriging, QRF, Oracle, pilot draws
│   ├── plot_curves.R  # shared figure code for Fig. 2 and Fig. 3
│   ├── 01–08_*.R      # main analyses (run in numeric order; see below)
│   ├── *.csv          # local dataset (Shigyo et al. 2022)
│   └── *.rds          # derived objects (regenerable; see "Data")
├── results/           # simulation outputs used for the figures
├── figures/
│   ├── main/          # Fig. 1–4 as used in the manuscript
│   └── supplementary/ # weight grid searches (S1, S2) and spatial trajectories (S3, S4)
├── LICENSE            # MIT (code)
└── DATA_LICENSE.md    # CC BY 4.0 (data; see "Data")
```

---

## Scripts

All scripts assume the working directory is `analysis/` and read/write files by bare name.
Weights are passed through environment variables (`BAS_WEIGHTS`, `BAS_WEIGHTS_LOCAL`, `BAS_WEIGHTS_GLOBAL`) and default to the values used in the paper: **(0.4, 0.2, 0.4)** locally and **(0.2, 0.6, 0.2)** globally, chosen by the grid searches (05, 07).

**Tuning and evaluation use disjoint pilot sets.** Pilot sets are generated from `seed = 123 + sim`. The grid searches use sims 1–50 (local) and 1–10 (global); the main simulations (01, 03, 04) use sims 51–150 and 11–110 (`SIM_OFFSET`). 01 and 03 also run a sensitivity arm with equal weights (1/3, 1/3, 1/3), plotted by 09 (Fig. S5).

With fewer than 10 training sites the kriging engine uses a nugget-only model (the covariate regression) instead of an automatically fitted variogram, because `automap::autofitVariogram` crashed the R session on some 5-point pilot sets (`MIN_AUTOFIT` in `bas_core.R`).

| Script | Produces | Figure |
|---|---|---|
| `01_BAS_simulation_Fig2AB.R` | Case 1 (local catchment, N = 53, regression kriging via `gstat` + `automap`) — BAS vs Random vs Oracle, 100 simulations → `simulation_results_Fig2AB.rds` | Fig. 2A–B |
| `02_BAS_Precompute_GlobalFungi.R` | `K_GlobalFungi.rds` (5000 × 5000 Hellinger kernel, not needed by v3 but kept for the LCBD comparison) and `comm_pa_sp_GlobalFungi.rds` (sparse presence/absence) — run once | — |
| `03_BAS_simulation_Fig3AB.R` | Case 2 (GlobalFungi, N = 5000, Quantile Regression Forest via `ranger`, 200 trees) — BAS vs Random vs Oracle, 100 simulations → `simulation_results_Fig3AB.rds` | Fig. 3A–B |
| `04_BAS_plot_Fig4.R` | Cumulative discovery of range-restricted taxa (local: ASVs at ≤ 2 of 53 sites; global: SHs at ≤ 5 of 5,000 sites) and the share of sampled sites south of 20°N, BAS vs Random vs Oracle, 100 simulations → `fig4_trajectories.rds` | Fig. 4A–C |
| `05_BAS_WeightGridSearch_local.R` | 66-point grid over (W_rich, W_uniq, W_unc), local dataset, 50 simulations each → `gridsearch_results_local.rds` | — |
| `06_plot_WeightGridSearch.R` | Ternary plots of 05 and 07 | Fig. S1, S2 |
| `07_BAS_WeightGridSearch_GlobalFungi.R` | Same grid on GlobalFungi (10 simulations, 100 trees — reduced for run time) → `gridsearch_results_global.rds` | — |
| `08_BAS_SpatialTrajectory_Fig2_Fig3.R` | Snapshots / animations of where BAS samples over time | Fig. S3, S4 |
| `09_plot_FigS5_equal_weights.R` | Tuned vs equal weights on the evaluation runs of 01 and 03 | Fig. S5 |

### Run order

```r
setwd("<path>/BAS_publication/analysis")
source("01_BAS_simulation_Fig2AB.R")        # ~2 min on 6 cores
source("02_BAS_Precompute_GlobalFungi.R")   # once; needs comm_hel_GlobalFungi.rds (included)
source("03_BAS_simulation_Fig3AB.R")        # ~12 min on 10 cores
source("04_BAS_plot_Fig4.R")                 # ~12 min on 10 cores
source("05_BAS_WeightGridSearch_local.R")   # ~15 min on 6 cores
source("07_BAS_WeightGridSearch_GlobalFungi.R")   # ~55 min on 8 cores
source("06_plot_WeightGridSearch.R")
source("08_BAS_SpatialTrajectory_Fig2_Fig3.R")
source("09_plot_FigS5_equal_weights.R")
```

> [!NOTE]
> If your R is linked to a multithreaded BLAS (OpenBLAS, MKL), set `OPENBLAS_NUM_THREADS=1` / `OMP_NUM_THREADS=1` in the environment **before** starting R. Scripts `01`, `03`, `05` and `07` spawn PSOCK worker processes, and each worker otherwise reserves several GB of virtual memory for BLAS thread buffers.

### Earlier versions

An earlier version of the analysis (git history before October 2026) used the abundance-weighted LCBD as the uniqueness term and a local regression/IDW predictor at the local scale. Replacing LCBD by novelty raised the share of the Random→Oracle gap closed at n = 3,000 on GlobalFungi from 43% to 76% under otherwise identical settings; the comparison is reported in the manuscript.

---

## Requirements

R 4.5.2 (Windows 11, 16 threads, 32 GB RAM was used). Package versions used for the results in `results/`:

| Package | Version | | Package | Version |
|---|---|---|---|---|
| gstat | 2.1.5 | | ranger | 0.18.0 |
| automap | 1.1.20 | | Matrix | 1.7.4 |
| sp | 2.2.1 | | dplyr | 1.2.0 |
| sf | 1.0.24 | | data.table | 1.18.2.1 |
| foreach | 1.5.2 | | terra | 1.9.27 |
| doParallel | 1.0.17 | | ggplot2 | 4.0.2 |
| vegan | 2.7.2 | | patchwork | 1.3.2 |
| ggtern | 4.0.0 | | gridExtra | 2.3 |
| gganimate / gifski | 1.0.11 / 1.32.0.2 | | | |

Also needed when rebuilding the GlobalFungi subset from raw files: `geodata`.

```r
install.packages(c("gstat","automap","sp","sf","vegan","dplyr","foreach","doParallel",
                   "ranger","Matrix","data.table","terra","geodata","ggplot2",
                   "patchwork","ggtern","gridExtra","gganimate","gifski"))
```

---

## Data

### Case 1 — local catchment (Shigyo et al. 2022) — included

| File | Description |
|---|---|
| `analysis/comm_data_Shigyo_et_al_2022.csv` | ASV table (rows = sites, columns = ASVs) |
| `analysis/site_data_Shigyo_et_al_2022.csv` | Site metadata (SiteID, lon, lat, pH, WC) |

### Case 2 — GlobalFungi v5 (Větrovský et al. 2020)

The raw GlobalFungi files (~12 GB) are **not** included. Download from https://globalfungi.com/:

| File needed | GlobalFungi filename |
|---|---|
| `GlobalFungi_5_sample_metadata.txt` | Sample metadata, GlobalFungi v5 |
| `GlobalFungi_5_SH_abundance_ITS1_ITS2.txt` | SH abundance table, GlobalFungi v5 |

The 5,000-sample subset was drawn (seed 123) from soil and topsoil samples with coordinates, measured pH and WorldClim values; it spans 341 studies, four sequencing platforms and ITS1/ITS2/both barcodes, and the SH table (ITS1 and ITS2 pooled, as released) was converted to presence/absence without rarefaction. The transformed objects derived from them (`comm_pa_GlobalFungi.rds`, `comm_pa_sp_GlobalFungi.rds`, `comm_hel_GlobalFungi.rds`, 12-column `master_data_GlobalFungi.rds`) are included under CC BY 4.0; see `DATA_LICENSE.md`.

---

## Results

| File | Content |
|---|---|
| `results/simulation_results_Fig2AB.rds` | Cumulative richness per (simulation, n, strategy) for Case 1: BAS / BAS (equal weights) / Random / Oracle × 100 evaluation simulations (sims 51–150) × n = 5…30 |
| `results/simulation_results_Fig3AB.rds` | Cumulative richness per (simulation, n, strategy) for Case 2: BAS / BAS (equal weights) / Random / Oracle × 100 evaluation simulations (sims 11–110) × n = 1000…3000 |
| `results/gridsearch_results_local.rds` | Local weight grid search (66 combinations × 50 simulations) |
| `results/gridsearch_results_global.rds` | Global weight grid search (66 × 10) |
| `results/fig4_trajectories.rds` | Range-restricted taxa found and southern share per (simulation, n, strategy), both scales (Fig. 4) |

---

## Author

**Nobuhiko Shigyo** (corresponding author)
Graduate School of Horticulture, Chiba University, Japan — shigyo@chiba-u.jp

## Funding

Japan Science and Technology Agency (JST) ACT-X (Grant No. JPMJAX25L5).

## License

Code: MIT (see `LICENSE`). Data: CC BY 4.0 (see `DATA_LICENSE.md`).
