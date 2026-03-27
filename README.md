# lppls_mult_ci

R implementation of the **Multi-Scale LPPLS Confidence Indicator** pipeline from:

> Gabauer, D., Gupta, R., Karmakar, S., & Nielsen, J. (2022). *Stock Market Bubbles and the Forecastability of Gold Returns (and Volatility)*. University of Pretoria Working Paper 2022-28.

The paper shows that bubble signals extracted from G7 and BRICS stock markets — summarised by the LPPLS Confidence Indicator — significantly improve short-horizon forecasts of weekly gold returns and volatility over an AR(1) benchmark.

---

## What this project does

1. **Detects stock market bubbles** across 12 markets (G7 + BRICS) by fitting the Log-Periodic Power Law Singularity (LPPLS) model to rolling windows of log price-dividend ratios at three time scales (short, medium, long). Qualified fits are aggregated into a Confidence Indicator (CI) — a number in [0, 1] measuring how strongly a bubble pattern is present.

2. **Forecasts gold returns and volatility** using LASSO regression on the 42 (G7) or 72 (G7+BRICS) CI series as predictors, benchmarked against AR(1) via the Clark–West test.

The CI computation is also exposed as a **generic function** (`compute_ci_asset`) that works for any log-price series — equities, crypto, commodities, or FX — with user-configurable time scales.

---

## Repository layout

```
R/
  lppls_fit.R    — LPPLS single-window calibration (Filimonov-Sornette
                   reparameterisation; vectorised grid search + SLSQP)
  lppls_ci.R     — Multi-scale CI sweep; generic compute_ci_asset()
  ci_batch.R     — Resumable batch checkpointing for long CI runs
  data.R         — Data loading, weekly alignment, predictor matrix
  fetch_data.R   — Auto-download gold (FRED) and P/D proxies (Yahoo Finance)
  tvpgarch.R     — Time-varying GARCH(1,1) volatility (stub)
  forecast.R     — LASSO + block-mean forecasting, AR(1) benchmark (stub)
  evaluation.R   — Clark-West and Diebold-Mariano tests (stub)
  plots.R        — ggplot2 figures (stub)
specs/
  model_spec.md  — All equations with variable definitions
  model_desc.md  — Narrative descriptions of each model
tests/
  test_lppls_fit.R   — 24 unit tests for LPPLS calibration
  test_lppls_ci.R    — 42 unit tests for CI computation and generic API
main.R           — Orchestration script (runs the full pipeline)
```

---

## Quickstart

### 1. Install dependencies

```r
# install.packages("renv")
renv::restore()
```

Key packages: `nloptr`, `glmnet`, `data.table`, `ggplot2`, `patchwork`, `future`, `furrr`, `future.apply`.

### 2. Get data

The pipeline needs:
- Daily LBMA gold prices (`data/gold_prices.csv`)
- Daily log price-dividend ratios for each country (`data/pd_ratios/{country}.csv`)

**Option A — auto-download free proxies** (ETFs from Yahoo Finance + gold from FRED):

```r
source("R/fetch_data.R")
fetch_all_data()   # writes to data/
```

> BRICS ETFs have limited history (Russia from 2007, India/China from 2011-2012). For a faithful replication of the paper's results, replace the CSVs with Refinitiv Datastream exports.

**Option B — bring your own CSVs** in the expected format:

| File | Columns |
|------|---------|
| `data/gold_prices.csv` | `date` (YYYY-MM-DD), `price` (USD) |
| `data/pd_ratios/{country}.csv` | `date`, `pd_ratio` |

### 3. Run the pipeline

```r
Rscript main.R
```

Outputs are written to `output/tables/` (CSV) and `output/figures/` (PDF).

---

## Using the CI for other assets

`compute_ci_asset` works for any log-price series:

```r
source("R/lppls_fit.R")
source("R/lppls_ci.R")

# data.table input (auto-detects the value column)
ci <- compute_ci_asset(data.table(date = dates, log_price = log(btc_price)))

# Numeric vector input
ci <- compute_ci_asset(log(btc_price), dates = dates)

# Custom time scales (e.g., for a daily crypto series)
crypto_scales <- list(
  short  = list(dt_min = 20L,  dt_max =  60L, step = 2L),
  medium = list(dt_min = 60L,  dt_max = 200L, step = 2L),
  long   = list(dt_min = 200L, dt_max = 500L, step = 2L)
)
ci <- compute_ci_asset(log(btc_price), dates = dates, scales = crypto_scales)

# Multiple assets in parallel
ci_list <- compute_ci_assets(list(btc = btc_dt, eth = eth_dt))
```

Output is a `data.table(date, pos_short, neg_short, pos_med, neg_med, pos_long, neg_long)`.

### Batch checkpointing for long runs

```r
source("R/ci_batch.R")

# Split into 10 resumable chunks — safe to interrupt and restart
ci <- compute_ci_batched(log(btc_price), dates,
                         n_batches = 10, cache_dir = "cache/btc_ci")

# Check progress mid-run
ci_batch_status("cache/btc_ci")

# Merge completed batches manually
ci <- merge_ci_batches("cache/btc_ci")
```

---

## Performance

The CI sweep is computationally intensive: for 315 weekly observations and the default three-scale setup, ~113,400 LPPLS fits are needed per asset.

| Setting | Description |
|---------|-------------|
| `method = "grid"` (default) | Vectorised 20×12×12 grid scan over `(tc, m, ω)` + SLSQP refinement. ~2–3× faster than random starts per fit. |
| `method = "random"` | Original multi-start SLSQP (5 random starts, up to 500 iters each). |
| `method = "grid_only"` | Grid scan only, no SLSQP refinement. Fastest; less accurate. |
| `parallel_t2 = TRUE` (default in `compute_ci_asset`) | Parallelises the t2 loop via `future.apply`. Set `plan(multisession, workers = N)` before calling. |
| `parallel_t2 = FALSE` (default in `compute_ci_assets`) | Sequential t2 loop; outer parallelism is over assets instead. |
| `compute_ci_batched(..., n_batches = N)` | Splits computation into N resumable chunks. |

Typical runtime on an 8-core machine with `plan(multisession, workers = 8)` and `method = "grid"`: **~60–90 seconds per asset** for the full three-scale sweep.

---

## Tests

```r
Rscript --no-init-file tests/run_tests.R
# or from within R:
testthat::test_dir("tests/")
```

66 tests across two files (24 for LPPLS calibration, 42 for the CI and generic API), all passing.

---

## Methodology reference

Full mathematical specifications are in [`specs/model_spec.md`](specs/model_spec.md) and [`specs/model_desc.md`](specs/model_desc.md), covering:

- LPPLS model (Filimonov-Sornette 2013 reparameterisation)
- Multi-scale CI window construction and filter conditions
- LASSO + block-mean forecasting (Karmakar et al. 2021)
- Clark–West and Diebold–Mariano forecast evaluation tests
- TVPGARCH volatility (Karmakar & Roy 2021)
- TVP-VAR Total Connectedness Index (Antonakakis et al. 2020)
