# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Replication of the **Multi-Scale LPPLS Confidence Indicator** pipeline from:

> Gabauer, Gupta, Karmakar & Nielsen (2022). *Stock Market Bubbles and the Forecastability of Gold Returns (and Volatility)*. University of Pretoria Working Paper 2022-28.

The reference paper is `Gabauer_Stock_2024.pdf` (not committed). The primary language is **R**. The R-btw MCP server is available for running R code interactively.

## Research Goal

1. **Bubble detection**: Compute Multi-Scale LPPLS Confidence Indicators (LPPLS-CI) at short-, medium-, and long-term horizons for 12 stock markets (G7 + BRICS), applied to log price-dividend ratios.
2. **Gold forecasting**: Use those indicators (42 predictors for G7-only; 72 for G7+BRICS) to forecast weekly gold returns and volatility via LASSO + block-mean quantile regression, benchmarked against AR(1).

## Data

- **Gold prices**: Daily from LBMA, converted to weekly averages; compute log-returns.
- **Stock price-dividend ratios**: Daily from Refinitiv Datastream (local currencies, log-transformed). Free proxy: country ETFs via Yahoo Finance (`R/fetch_data.R`).
- **G7 sample**: Week of Jan 7, 1973 – Sep 13, 2020 (2489 weekly obs). Countries: Canada, France, Germany, Italy, Japan, UK, US.
- **G7+BRICS sample**: Feb 14, 1999 – Sep 13, 2020 (1120 weekly obs). Adds: Brazil, Russia, India, China, South Africa.

## Methodology

Full documentation lives in `specs/`:
- **[`specs/model_desc.md`](specs/model_desc.md)** — narrative description of each model.
- **[`specs/model_spec.md`](specs/model_spec.md)** — all equations with variable definitions.

Pipeline summary:

1. **LPPLS calibration** (spec §1): Filimonov–Sornette reparameterisation. Solve 4 linear params analytically via 4×4 normal equations; solve 3 nonlinear params `{tc, m, ω}` via SLSQP (`nloptr`, `"NLOPT_LD_SLSQP"`). Default initialisation: vectorised grid search over `(tc, m, ω)`.
2. **Multi-Scale LPPLS-CI** (spec §2): Sweep windows at short [30,90], medium [90,300], long [300,745] trading-day ranges (step 2). Five filter conditions on `m`, `ω`, `tc`, oscillation count `O`, damping ratio `D`. CI = qualified / total fits. Six indicators per country per week.
3. **LASSO + block-mean forecasting** (spec §3): AR(1) vs AR(1) + LPPLS-CI predictors via `cv.glmnet`. Block-mean residual correction for h-step forecasts. Evaluate with Clark–West test.
4. **TVPGARCH volatility** (spec §4): Time-varying GARCH(1,1) via Epanechnikov kernel local likelihood (`bn = 0.25`).
5. **TCI via TVP-VAR** (spec §5, Appendix): Optional connectedness aggregation via TVP-VAR + GFEVD.

## Module Map

| File | Status | Key exports |
|------|--------|-------------|
| `R/lppls_fit.R` | Complete | `lppls_basis`, `lppls_linear_solve`, `lppls_cost`, `lppls_gradient`, `lppls_grid_search`, `lppls_fit` |
| `R/lppls_ci.R` | Complete | `lppls_filter`, `compute_ci_t2`, `compute_ci_scale`, `compute_ci_asset`, `compute_ci_assets`, `compute_ci_country`, `compute_ci_all` |
| `R/ci_batch.R` | Complete | `compute_ci_batched`, `merge_ci_batches`, `ci_batch_status` |
| `R/data.R` | Complete | `load_gold`, `load_pd_ratios`, `align_to_weekly`, `build_predictor_matrix` |
| `R/fetch_data.R` | Complete | `fetch_gold_fred`, `fetch_pd_yahoo`, `fetch_all_data` |
| `R/tvpgarch.R` | Stub | `tvpgarch_series` |
| `R/forecast.R` | Stub | `ar1_forecast`, `block_mean_forecast`, `run_all_forecasts` |
| `R/evaluation.R` | Stub | `build_cw_table` |
| `R/plots.R` | Stub | `plot_gold_data`, `plot_ci_country` |

Stubs are sourced in `main.R` but not yet implemented.

## Key Function Signatures

```r
# Single asset (any log-price series)
compute_ci_asset(x, dates = NULL, value_col = NULL,
                 scales = SCALES, n_starts = 3L,
                 method = c("grid", "random", "grid_only"),
                 parallel_t2 = TRUE)

# Multiple assets in parallel (outer loop over assets)
compute_ci_assets(asset_list, value_col = NULL,
                  scales = SCALES, n_starts = 3L, method = "grid")

# Resumable batch execution
compute_ci_batched(x, dates = NULL, ...,
                   n_batches = 10L, cache_dir = "cache/ci_batches",
                   overwrite = FALSE)

# LPPLS fit with grid search initialisation
lppls_fit(t, lnp, n_starts = 3L,
          method = c("grid", "random", "grid_only"))
```

## Performance Notes

The CI sweep is expensive: ~113,400 LPPLS fits per asset for the full three-scale default. Key levers:

- **`method = "grid"` (default)**: vectorised 20×12×12 grid scan over `(tc, m, ω)` + SLSQP refinement from top-3 candidates; ~2–3× faster per fit than `"random"`.
- **`parallel_t2 = TRUE`** (default in `compute_ci_asset`): parallelises the inner t2 loop via `future.apply::future_lapply`. Set `plan(multisession, workers = N)` before calling. **Do not** combine with `compute_ci_assets` (which already parallelises over assets) — `compute_ci_assets` forces `parallel_t2 = FALSE`.
- **`compute_ci_batched`**: splits t2 into N chunks, saves each as RDS, resumes on restart.
- Typical runtime: ~60–90 s per asset on an 8-core machine.

## R Conventions

### Data manipulation
Use `data.table` throughout. Prefer `:=` for in-place assignment, `.SD`/`.SDcols` for column subsets, `on=`/`merge()` for joins. Avoid `dplyr`/`tidyr`.

### Plotting
`ggplot2` + `patchwork`. Avoid base R `plot()` and `lattice`.

### Parallel computing
`future` + `furrr` + `future.apply`. Set `plan(multisession)` once at the top of the script. Avoid nested future workers: when `compute_ci_assets` parallelises over assets, the inner t2 loop must be sequential (`parallel_t2 = FALSE`). Avoid `parallel::mclapply` and `foreach`.

### Package management
`renv`. Run `renv::restore()` after cloning. Record new packages with `renv::snapshot()` before committing.

## Running Tests

```r
Rscript --no-init-file tests/run_tests.R
```

66 tests: `tests/test_lppls_fit.R` (24) and `tests/test_lppls_ci.R` (42). All CI tests use `n_starts = 1L` and a micro-scale (`dt = [30,35]`) for speed; the full-scale `"short"` scale is used only for the qualitative bubble-vs-noise test with `method = "random"`.
