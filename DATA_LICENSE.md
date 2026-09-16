# Data license

The **code** in this repository is released under the MIT License (see `LICENSE`).

The **data files** listed below are released under the
[Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/legalcode).
You may share and adapt them for any purpose, provided you give appropriate credit to the sources named here.

## Files covered

| File(s) | Content | Source to cite |
|---|---|---|
| `analysis/comm_data_Shigyo_et_al_2022.csv`, `analysis/site_data_Shigyo_et_al_2022.csv` | Local-catchment ASV table and site metadata (Case 1) | Shigyo et al. (2022) |
| `analysis/master_data_GlobalFungi.rds` | 5,000-site table: sample ID, source-paper ID/DOI, coordinates (WGS84 and Mollweide), pH, MAT, MAP, SOC, observed SH richness | GlobalFungi (pH, coordinates, sample and paper IDs); WorldClim 2.1 (MAT, MAP); SoilGrids (SOC) |
| `analysis/comm_pa_GlobalFungi.rds`, `analysis/comm_pa_sp_GlobalFungi.rds` | Presence/absence of 49,245 UNITE species hypotheses (SH) across the 5,000 sites | GlobalFungi |
| `analysis/comm_hel_GlobalFungi.rds` | Hellinger-transformed SH abundances across the 5,000 sites | GlobalFungi |
| `results/*.rds`, `results/*.csv` | Simulation outputs | this repository |

The GlobalFungi-derived files are a **filtered and transformed subset** of GlobalFungi Release 5
(soil samples with pH, random subsample of 5,000 sites; see `03_BAS_simulation_Fig3AB.R` Step 1).
The raw read-count table and the complete sample metadata are **not** redistributed here;
obtain them from https://globalfungi.com/ under the GlobalFungi citation terms.

## Required attribution

- **GlobalFungi:** Větrovský T, Morais D, Kohout P, et al. (2020). GlobalFungi, a global database of fungal occurrences from high-throughput-sequencing metabarcoding studies. *Scientific Data*, 7, 228. https://doi.org/10.1038/s41597-020-0567-7 — and the primary studies listed in `master_data_GlobalFungi.rds` (`paper_ID`, `paper_doi`) when individual samples are re-used.
- **WorldClim 2.1 (BIO1, BIO12):** Fick SE, Hijmans RJ (2017). WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. *International Journal of Climatology*, 37, 4302–4315. https://doi.org/10.1002/joc.5086 (CC BY 4.0)
- **SoilGrids (SOC):** Poggio L, de Sousa LM, Batjes NH, et al. (2021). SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty. *SOIL*, 7, 217–240. https://doi.org/10.5194/soil-7-217-2021 (CC BY 4.0)
- **Local catchment data:** Shigyo N, et al. (2022). [full reference to be inserted]
