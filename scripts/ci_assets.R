# scripts/ci_assets.R
#
# Fetch prices from Yahoo Finance, compute LPPLS Confidence Indicators in
# resumable batches, and save CI results + plots for a configurable asset list.
#
# Usage (from the project root):
#   Rscript scripts/ci_assets.R
#
# Resumable: re-running skips already-completed batch chunks (OVERWRITE=FALSE).
# Set OVERWRITE=TRUE to recompute from scratch.

# =============================================================================
# USER CONFIGURATION — edit this block to change assets / settings
# =============================================================================

start <- "2024-1-1"
end <- "2026-03-28"
ASSETS <- list(
  bova11  = list(ticker = "BOVA11.SA", start = start, end = end),
  petr4  = list(ticker = "PETR4.SA", start = start, end = end),
  vale3  = list(ticker = "VALE3.SA", start = start, end = end),
  bbdc4  = list(ticker = "BBDC4.SA", start = start, end = end),
  itsa4  = list(ticker = "ITSA4.SA", start = start, end = end),
  ivvb11  = list(ticker = "IVVB11.SA", start = start, end = end),
  small11  = list(ticker = "SMAL11.SA", start = start, end = end),
  wege3  = list(ticker = "WEGE3.SA", start = start, end = end),
  b3sa3  = list(ticker = "B3SA3.SA", start = start, end = end),
  ixic  = list(ticker = "^IXIC", start = "1997-1-1", end = "2001-12-31"),
  bvsp  = list(ticker = "^BVSP", start = start, end = end),
  rut  = list(ticker = "^RUT", start = start, end = end),
  oil  = list(ticker = "CL=F", start = start, end = end),
  silver  = list(ticker = "SL=F", start = start, end = end),
  btc  = list(ticker = "BTC-USD", start = start, end = end),
  mxn  = list(ticker = "MXN=X", start = start, end = end),
  brl  = list(ticker = "BRL=X", start = start, end = end),
  aud  = list(ticker = "AUDUSD=X", start = start, end = end)
)

# Number of chunks per asset — increase for very long series (>1000 obs)
N_BATCHES <- 10L

# Parallel workers for the inner t2 loop
N_WORKERS <- min(4L, parallelly::availableCores() - 1L)
N_WORKERS <- 12L

# Set TRUE to recompute even if batch cache already exists
OVERWRITE <- FALSE

# Custom LPPLS scales — NULL uses the default (short/medium/long trading-day windows)
# Example for shorter crypto horizons:
#   CUSTOM_SCALES <- list(
#     short  = list(dt_min = 20L, dt_max =  60L, step = 2L),
#     medium = list(dt_min = 60L, dt_max = 200L, step = 2L),
#     long   = list(dt_min = 200L, dt_max = 500L, step = 2L)
#   )
CUSTOM_SCALES <- NULL
CUSTOM_SCALES <- list(
  short  = list(dt_min = 20L,  dt_max = 60L,  step = 5L),
  medium = list(dt_min = 60L,  dt_max = 120L, step = 5L),
  long   = list(dt_min = 120L, dt_max = 250L, step = 5L)
)


# Output locations
PRICE_DIR <- "data/assets"
CI_DIR    <- "output/ci"
CACHE_DIR <- "cache/assets"
PLOT_DIR  <- "output/figures/assets"

# =============================================================================
# SETUP
# =============================================================================

# Ensure working directory is the project root when run via Rscript
if (!interactive()) {
  args        <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("--file=", "", args[grep("--file=", args)])
  if (length(script_file) == 1L)
    setwd(normalizePath(file.path(dirname(script_file), "..")))
}

library(future)
library(future.apply)
library(data.table)
library(ggplot2)

plan(multisession, workers = N_WORKERS)

source("R/fetch_data.R")
source("R/lppls_fit.R")
source("R/lppls_ci.R")
source("R/ci_batch.R")
source("R/plots.R")

for (d in c(PRICE_DIR, CI_DIR, PLOT_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

scales_to_use <- if (is.null(CUSTOM_SCALES)) SCALES else CUSTOM_SCALES

# =============================================================================
# PER-ASSET PIPELINE
# =============================================================================

results <- list()

for (name in names(ASSETS)) {
  cfg <- ASSETS[[name]]
  message(sprintf("\n=== %s (%s) ===", toupper(name), cfg$ticker))

  # ---- 1. Fetch or load prices ----------------------------------------
  price_path <- file.path(PRICE_DIR, paste0(name, ".csv"))

  if (!file.exists(price_path)) {
    price_dt <- fetch_asset_yahoo(
      name,
      cfg$ticker,
      output_dir = PRICE_DIR,
      start_date = cfg$start,
      end_date   = cfg$end
    )
  } else {
    message("  Loading cached prices from ", price_path)
    price_dt <- fread(price_path)
    price_dt[, date := as.Date(date)]
  }

  if (is.null(price_dt) || nrow(price_dt) == 0L) {
    warning(sprintf("  No price data for %s — skipping.", name))
    next
  }

  message(sprintf("  %d price rows (%s to %s)",
                  nrow(price_dt),
                  format(min(price_dt$date)),
                  format(max(price_dt$date))))

  # ---- 2. Compute CI in resumable batches -----------------------------
  ci_dt <- tryCatch(
    compute_ci_batched(
      price_dt,
      value_col   = "log_price",
      scales      = scales_to_use,
      n_batches   = N_BATCHES,
      cache_dir   = file.path(CACHE_DIR, name),
      overwrite   = OVERWRITE,
      parallel_t2 = TRUE
    ),
    error = function(e) {
      warning(sprintf("  CI computation failed for %s: %s",
                      name, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(ci_dt)) next

  # ---- 3. Save CI to CSV ----------------------------------------------
  ci_path <- file.path(CI_DIR, paste0(name, "_ci.csv"))
  fwrite(ci_dt, ci_path)
  message(sprintf("  CI saved: %d rows → %s", nrow(ci_dt), ci_path))

  # ---- 4. Plot --------------------------------------------------------
  lnp_aligned <- price_dt[date %in% ci_dt$date, log_price]
  # p <- plot_ci_asset(ci_dt,
  #                    lnp_series = lnp_aligned,
  #                    asset_name = toupper(name))
  p <- plot_ci_price_all(ci_dt,
                     price_series = exp(lnp_aligned),
                     dates_price = ci_dt$date,
                     asset_name = toupper(name))
  plot_path <- file.path(PLOT_DIR, paste0(name, "_ci.png"))
  ggsave(plot_path, p, width = 10, height = 7)
  message(sprintf("  Plot saved → %s", plot_path))

  results[[name]] <- list(price = price_dt, ci = ci_dt)
}

plan(sequential)   # release parallel workers

if (length(results) == 0L) {
  message("\nNo assets completed successfully.")
} else {
  message(sprintf("\nDone. Completed: %s",
                  paste(names(results), collapse = ", ")))
  message(sprintf("CI files: %s/", CI_DIR))
  message(sprintf("Plots   : %s/", PLOT_DIR))
}
