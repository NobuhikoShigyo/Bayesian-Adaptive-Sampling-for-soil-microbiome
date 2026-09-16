# Bayesian Adaptive Sampling for Soil Microbiomes

R code, data and figures for the paper
*"We're only in it for the information: A framework for Bayesian adaptive sampling in soil microbiomes."*

A data-driven approach for deciding where to sample soils next.

> [!IMPORTANT]
> **Work in Progress** — this repository accompanies a manuscript in submission/review. Files may still change.

---

## Overview

This project advocates a shift in soil microbiome sampling from static spatial coverage to a dynamic philosophy of *information discovery*. We propose a Bayesian Adaptive Sampling (BAS) framework that reframes sampling as a sequential decision-making process. Each iteration selects the next batch of sites by optimizing a weighted acquisition function that combines

- predicted richness,
- predicted local contribution to beta diversity (LCBD; uniqueness), and
- model uncertainty.

The surrogate model that produces these three quantities is interchangeable. In the paper, Universal Kriging is used at the local scale (Case 1) and Quantile Regression Forests at the global scale (Case 2); a Gaussian-process variant of Case 2 is provided as a supplementary robustness check.

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
│   ├── 01–08_*.R      # main analyses (run in numeric order; see below)
│   ├── supplementary_GP/
│   ├── *.csv          # local dataset (Shigyo et al. 2022)
│   └── *.rds          # derived objects (regenerable; see "Data")
├── results/           # simulation outputs used for the figures
├── figures/
│   ├── main/          # Fig. 1–4 as used in the manuscript
│   └── supplementary/ # Oracle bounds, weight grid searches, spatial trajectories, GP variant
├── LICENSE            # MIT (code)
└── DATA_LICENSE.md    # CC BY 4.0 (data; see "Data")
```

---

## Scripts

All scripts assume the working directory is `analysis/` and read/write files by bare name.

| Script | Produces | Figure |
|---|---|---|
| `01_BAS_simulation_Fig2AB.R` | Case 1 (local catchment, N = 53, Universal Kriging via `gstat`) — BAS vs Random vs Oracle, 100 simulations → `simulation_results_Fig2AB.rds` | Fig. 2A–B, Fig. S1 |
| `02_BAS_Precompute_GlobalFungi.R` | `K_GlobalFungi.rds` (5000 × 5000 Hellinger kernel) and `comm_pa_sp_GlobalFungi.rds` (sparse presence/absence) — run once | — |
| `03_BAS_simulation_Fig3AB.R` | Case 2 (GlobalFungi, N = 5000, Quantile Regression Forest via `ranger`) — BAS vs Random vs Oracle, 100 simulations → `simulation_results_Fig3AB_v2.rds` | Fig. 3A–B, Fig. S2 |
| `04_BAS_plot_Fig4.R` | Rare-biosphere comparison: relative abundance of taxa unique to BAS vs unique to Random | Fig. 4A–B |
| `05_BAS_WeightGridSearch_local.R` | 66-point grid over (W_rich, W_uniq, W_unc), local dataset, 100 simulations each → `gridsearch_results.rds`, `gridsearch_summary.csv` | Fig. S3 |
| `06_BAS_WeightGridSearch_local_MultiN.R` | Post-processing of 05 at multiple evaluation depths | Fig. S4 |
| `07_BAS_WeightGridSearch_GlobalFungi.R` | Same grid on GlobalFungi (20 simulations, 50 trees — reduced for run time) → `gridsearch_results_GlobalFungi_light.rds` | Fig. S5 |
| `08_BAS_SpatialTrajectory_Fig2_Fig3.R` | Snapshots / animations of where BAS samples over time | Fig. S6–S7 |
| `supplementary_GP/S01_BAS_simulation_Fig2AB_GP.R` | Case 1 with a GP (kriging) engine for all three acquisition terms | Fig. S8 |
| `supplementary_GP/S02_BAS_simulation_Fig3AB_GP_RobustGaSP.R` | Case 2 with a full GP (`RobustGaSP`, 6-D Matérn 5/2 kernel on X, Y, pH, MAT, MAP, SOC). Shipped with `N_SIMULATIONS = 10`; see note below | — |

### Run order

```r
setwd("<path>/BAS_publication/analysis")
source("01_BAS_simulation_Fig2AB.R")        # ~minutes; writes master_data_local.rds etc.
source("02_BAS_Precompute_GlobalFungi.R")   # once; needs comm_hel_GlobalFungi.rds (included)
source("03_BAS_simulation_Fig3AB.R")        # ~15–35 min on 14 cores
source("04_BAS_plot_Fig4.R")                # needs objects from 01 and 03 in the session
```

`05`–`08` and the supplementary scripts are independent of one another once `01`–`03` have been run.

> [!NOTE]
> If your R is linked to a multithreaded BLAS (OpenBLAS, MKL), set `OPENBLAS_NUM_THREADS=1` / `OMP_NUM_THREADS=1` in the environment **before** starting R. Scripts `01`, `03`, `05` and `07` spawn up to `detectCores() - 1` worker processes, and each worker otherwise reserves several GB of virtual memory for BLAS thread buffers.

### Note on the GP variant of Case 2

`S02` is the script that produced `results/supplementary_GP/`. It is a 10-simulation run (log included) demonstrating that, with covariates in the kernel, the GP's uncertainty term stays discriminative across all batches (coefficient of variation 0.29–0.34) and BAS still beats Random. It starts from `N_INITIAL = 100` rather than the 1,000 used in `03`; set `N_INITIAL <- 1000` and `N_SIMULATIONS <- 30` for a like-for-like comparison (≈ 17 h at 33 min/simulation). The QRF engine outperforms it at n = 3000 (45,915 vs 44,975 cumulative SH; Random 43,325; Oracle 48,939), which is why QRF is the main analysis.

---

## Requirements

R 4.5.2 (Windows 11, 16 cores, 32 GB RAM was used). Package versions used for the results in `results/`:

| Package | Version | | Package | Version |
|---|---|---|---|---|
| vegan | 2.7.2 | | ranger | 0.18.0 |
| dplyr | 1.2.0 | | Matrix | 1.7.4 |
| gstat | 2.1.5 | | data.table | 1.18.2.1 |
| sp | 2.2.1 | | terra | 1.9.27 |
| sf | 1.0.24 | | ggplot2 | 4.0.2 |
| adespatial | 0.3.28 | | patchwork | 1.3.2 |
| foreach | 1.5.2 | | ggtern | 4.0.0 |
| doParallel | 1.0.17 | | gridExtra | 2.3 |
| RhpcBLASctl | 0.23.42 | | gganimate / gifski | 1.0.11 / 1.32.0.2 |
| automap (S01) | 1.1.20 | | RobustGaSP (S02) | 0.6.8 |

Also needed when rebuilding the GlobalFungi subset from raw files: `geodata`, `ggpubr` (Fig. 4).

```r
install.packages(c("vegan","dplyr","gstat","sp","sf","adespatial","foreach","doParallel",
                   "RhpcBLASctl","ranger","Matrix","data.table","terra","geodata","ggplot2",
                   "patchwork","ggpubr","ggtern","gridExtra","gganimate","gifski",
                   "automap","RobustGaSP"))
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


---

## Results

| File | Content |
|---|---|
| `results/simulation_results_Fig2AB.rds` | Cumulative richness per (simulation, n, strategy) for Case 1: BAS / Random / Oracle × 100 simulations × n = 5…50 |
| `results/simulation_results_Fig3AB_v2.rds` | Cumulative richness per (simulation, n, strategy) for Case 2: BAS / Random / Oracle × 100 simulations × n = 1000…3000 |
| `results/gridsearch_results.rds`, `gridsearch_summary.csv` | Local weight grid search (66 combinations × 100 simulations) |
| `results/gridsearch_results_GlobalFungi_light.rds` | Global weight grid search (66 × 20) |
| `results/supplementary_GP/` | GP variant of Case 2 (10 simulations) + run log |

---

## Author

**Nobuhiko Shigyo** (corresponding author)
Graduate School of Horticulture, Chiba University, Japan — shigyo@chiba-u.jp

## Funding

Japan Science and Technology Agency (JST) ACT-X (Grant No. JPMJAX25L5).

