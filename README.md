# Bayesian Adaptive Sampling for Soil Microbiome

R code and data for the paper *"We're only in it for the information: A framework for Bayesian adaptive sampling in soil microbiomes."* A data-driven approach for deciding where to sample soils next.

> [!IMPORTANT]
> **Work in Progress**
> This repository is currently under construction and accompanies a manuscript that is in the submission or review process. Files and documentation may be updated frequently.

---

## Overview

This project advocates for a shift in soil microbiome sampling from static spatial coverage to a dynamic philosophy of 'information discovery.' We propose a Bayesian Adaptive Sampling (BAS) framework that reframes sampling as a sequential decision-making process. By using probabilistic spatial models, BAS identifies the next optimal sampling locations to maximize scientific return while minimizing effort.

The framework compares three strategies for accumulating microbial richness across sampling iterations:

| Strategy | Description |
|---|---|
| **BAS (+Env) / BAS (Balanced)** | Adaptive sampling using environmental covariates (pH, climate, SOC) + spatial coordinates |
| **BAS (-Env)** | Adaptive sampling using spatial coordinates only |
| **Random selection** | Random (non-adaptive) baseline |

Each iteration selects the next batch of sites by optimizing a weighted acquisition function combining predicted richness, predicted local contribution to beta diversity (LCBD), and model uncertainty.

---

## Authors

**Nobuhiko Shigyo** (Corresponding author)
Graduate School of Horticulture, Chiba University, Japan
E-mail: shigyo@chiba-u.jp

---

## Funding

This work was supported by Japan Science and Technology Agency (JST) ACT-X (Grant No. JPMJAX25L5).

---
