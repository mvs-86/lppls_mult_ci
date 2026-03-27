# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This project implements the **Multi-Scale LPPLS Confidence Indicator** methodology from:

> Gabauer, Gupta, Karmakar & Nielsen (2022). *Stock Market Bubbles and the Forecastability of Gold Returns (and Volatility)*. University of Pretoria Working Paper 2022-28.

The reference paper is `Gabauer_Stock_2024.pdf`. The primary language is **R** (the paper explicitly uses `cv.glmnet` for LASSO). The R-btw MCP server is available for running R code interactively.

## Research Goal

Two objectives:
1. **Bubble detection**: Compute Multi-Scale LPPLS Confidence Indicators (LPPLS-CI) at short-, medium-, and long-term horizons for 12 stock markets (G7 + BRICS), applied to log price-dividend ratios.
2. **Gold forecasting**: Use those indicators (42 predictors for G7-only; 72 for G7+BRICS) to forecast weekly gold returns and volatility via LASSO + block-mean quantile regression, benchmarked against AR(1).

## Data

- **Gold prices**: Daily from LBMA, converted to weekly averages; compute log-returns.
- **Stock price-dividend ratios**: Daily from Refinitiv Datastream (local currencies, log-transformed).
- **G7 sample**: Week of Jan 7, 1973 – Sep 13, 2020 (2489 weekly obs). Countries: Canada, France, Germany, Italy, Japan, UK, US.
- **G7+BRICS sample**: Feb 14, 1999 – Sep 13, 2020 (1120 weekly obs). Adds: Brazil, Russia, India, China, South Africa.

## Methodology

Full documentation lives in `specs/`:
- **[`specs/model_desc.md`](specs/model_desc.md)** — narrative description of each model: what it does, why it is used, inputs and outputs.
- **[`specs/model_spec.md`](specs/model_spec.md)** — all equations with variable definitions.

Summary of the pipeline:

1. **LPPLS calibration** (spec §1): Fit the Filimonov–Sornette reparameterised LPPLS model to each rolling window of the log price-dividend ratio. Solve 4 linear parameters analytically via the 4×4 normal equations; solve 3 nonlinear parameters `{tc, m, ω}` via SLSQP (`nloptr` with `"NLOPT_LD_SLSQP"`).

2. **Multi-Scale LPPLS-CI** (spec §2): For each weekly observation sweep nested windows at short [30,90], medium [90,300], and long [300,745] trading-day ranges in steps of 2 days. Count qualified fits (5 filter conditions on `m`, `ω`, `tc`, `O`, `D`) separately for positive (`B < 0`) and negative (`B > 0`) bubbles. CI = qualified / total fits. Yields 6 indicators per country per week.

3. **LASSO + block-mean forecasting** (spec §3): Compare AR(1) benchmark against AR(1) + LPPLS-CI predictors (42 for G7; 72 for G7+BRICS) estimated via `cv.glmnet`. Point forecast = fitted values + mean of block-mean residuals over horizon `h`. Evaluate with Clark–West test.

4. **TVPGARCH volatility** (spec §4): Fit time-varying GARCH(1,1) via Epanechnikov kernel local likelihood (`bn = 0.25`). Use estimated `σ²` series as dependent variable in the same forecasting framework.

5. **TCI via TVP-VAR** (spec §5, Appendix): Optional aggregation — compute total connectedness index across the 6 bubble indicators per country group using a TVP-VAR + GFEVD approach (Antonakakis et al. 2020). Used as an alternative predictor set in the forecasting model.

## R Conventions

### Data manipulation
Use `data.table` for all data wrangling. Prefer `:=` for in-place column assignment, `.SD` / `.SDcols` for column subsets, and `on=` / `merge()` for joins. Avoid `dplyr` / `tidyr`.

### Plotting
Use `ggplot2` for all plots and `patchwork` to compose multi-panel figures (e.g., the per-country positive/negative bubble indicator panels in Fig. 2 of the paper). Avoid base R `plot()` and `lattice`.

### Parallel computing
Use the `futureverse` ecosystem (`future` + `furrr` or `future.apply`) for parallelism. Set the plan once at the top of the script (e.g., `plan(multisession)`). The LPPLS-CI computation—sweeping hundreds of nested windows across 12 countries—is the primary parallelization target: parallelize over countries or over the outer t2 loop. Avoid `parallel::mclapply` and `foreach`.

### Package management
The project uses `renv` for reproducible package snapshots. Run `renv::restore()` after cloning to install dependencies. Record new packages with `renv::snapshot()` before committing. The lockfile (`renv.lock`) is the source of truth for package versions.

## Key R Packages

- `nloptr` — SLSQP nonlinear optimization for LPPLS calibration (`nloptr()` with `"NLOPT_LD_SLSQP"` algorithm)
- `glmnet` — LASSO (`cv.glmnet` for cross-validated lambda, use `lambda.min`)
- `rugarch` or custom kernel code — TVPGARCH estimation
- `lmtest` / `sandwich` — Clark-West test for forecast comparison
- `data.table` — data manipulation
- `ggplot2` + `patchwork` — plotting
- `future` + `furrr` — parallel computing
- `renv` — package management
