# main.R — Orchestration script for lppls_mult_ci
# Gabauer, Gupta, Karmakar & Nielsen (2022) replication
#
# Usage: Rscript main.R
# Prerequisite: renv::restore() to install all dependencies

library(future)
library(furrr)
library(future.apply)
library(data.table)

# ---------------------------------------------------------------------------
# Parallel plan — adjust workers to your machine
# ---------------------------------------------------------------------------
plan(multisession, workers = min(12L, parallelly::availableCores() - 1L))

# ---------------------------------------------------------------------------
# Source all modules
# ---------------------------------------------------------------------------
source("R/fetch_data.R")
source("R/data.R")
source("R/lppls_fit.R")
source("R/lppls_ci.R")
source("R/ci_batch.R")
source("R/tvpgarch.R")
source("R/forecast.R")
source("R/evaluation.R")
source("R/plots.R")

# ---------------------------------------------------------------------------
# Paths — edit to point at your data files
# ---------------------------------------------------------------------------
GOLD_PATH   <- "data/gold_prices.csv"        # daily LBMA gold prices (USD)
PD_DIR      <- "data/pd_ratios/"             # one CSV per country
OUTPUT_DIR  <- "output"

G7_COUNTRIES    <- c("canada", "france", "germany", "italy",
                     "japan",  "uk",     "us")
BRICS_COUNTRIES <- c("brazil", "russia", "india", "china", "south_africa")
ALL_COUNTRIES   <- c(G7_COUNTRIES, BRICS_COUNTRIES)

# ---------------------------------------------------------------------------
# Step 0 (optional): Download data if CSVs do not exist yet
# ---------------------------------------------------------------------------
if (!file.exists(GOLD_PATH)) {
  message("Data files not found — fetching from public sources...")
  fetch_all_data(
    countries  = ALL_COUNTRIES,
    output_dir = "data",
    start_date = "1970-01-01",
    end_date   = "2020-09-30"
  )
}

# ---------------------------------------------------------------------------
# Step 1: Load data
# ---------------------------------------------------------------------------
message("Loading gold prices...")
gold_dt <- load_gold(GOLD_PATH)

message("Loading price-dividend ratios...")
pd_list_g7 <- lapply(
  setNames(G7_COUNTRIES, G7_COUNTRIES),
  function(ctry) load_pd_ratios(file.path(PD_DIR, paste0(ctry, ".csv")), ctry)
)
pd_list_all <- lapply(
  setNames(ALL_COUNTRIES, ALL_COUNTRIES),
  function(ctry) load_pd_ratios(file.path(PD_DIR, paste0(ctry, ".csv")), ctry)
)

# Align weekly observations to gold series dates
dates_g7    <- gold_dt[date >= as.Date("1973-01-07") & date <= as.Date("2020-09-13"), date]
dates_brics <- gold_dt[date >= as.Date("1999-02-14") & date <= as.Date("2020-09-13"), date]

pd_input_g7 <- lapply(pd_list_g7, function(x) {
  align_to_weekly(x$lnp_daily, x$all_dates, dates_g7)
})
pd_input_all <- lapply(pd_list_all, function(x) {
  align_to_weekly(x$lnp_daily, x$all_dates, dates_brics)
})

# ---------------------------------------------------------------------------
# Step 2: Compute multi-scale LPPLS-CI
# ---------------------------------------------------------------------------
message("Computing LPPLS-CI for G7 (this is the slow step)...")
ci_g7  <- compute_ci_all(pd_input_g7)

message("Computing LPPLS-CI for G7 + BRICS...")
ci_all <- compute_ci_all(pd_input_all)

# ---------------------------------------------------------------------------
# Step 3: Build predictor matrices
# ---------------------------------------------------------------------------
X_g7  <- build_predictor_matrix(ci_g7,  dates_g7)    # 42 columns
X_all <- build_predictor_matrix(ci_all, dates_brics)  # 72 columns

# ---------------------------------------------------------------------------
# Step 4: Gold returns
# ---------------------------------------------------------------------------
y_g7    <- gold_dt[date %in% dates_g7,    log_return]
y_brics <- gold_dt[date %in% dates_brics, log_return]

# ---------------------------------------------------------------------------
# Step 5: TVPGARCH volatility
# ---------------------------------------------------------------------------
message("Fitting TVPGARCH...")
vol_g7    <- tvpgarch_series(y_g7)
vol_brics <- tvpgarch_series(y_brics)

# ---------------------------------------------------------------------------
# Step 6: Forecasting — returns and volatility
# ---------------------------------------------------------------------------
message("Running forecasting experiments...")
T_VEC <- c(100L, 250L, 500L, 1000L)  # 1000 for G7 only
H_VEC <- c(1L, 2L, 4L, 8L, 12L)

fc_ret_g7  <- run_all_forecasts(y_g7,    X_g7,  T_VEC,        H_VEC)
fc_ret_all <- run_all_forecasts(y_brics, X_all, T_VEC[-4L],   H_VEC)
fc_vol_g7  <- run_all_forecasts(vol_g7,  X_g7,  T_VEC,        H_VEC)
fc_vol_all <- run_all_forecasts(vol_brics, X_all, T_VEC[-4L], H_VEC)

# ---------------------------------------------------------------------------
# Step 7: Evaluation tables
# ---------------------------------------------------------------------------
message("Building CW/DM tables...")
tbl1 <- build_cw_table(fc_ret_g7,  fc_ret_all,  y_g7,    y_brics, "returns")
tbl2 <- build_cw_table(fc_vol_g7,  fc_vol_all,  vol_g7,  vol_brics, "volatility_tci")
tbl3 <- build_cw_table(fc_vol_g7,  fc_vol_all,  vol_g7,  vol_brics, "volatility")

fwrite(tbl1, file.path(OUTPUT_DIR, "tables", "table1_cw_returns.csv"))
fwrite(tbl2, file.path(OUTPUT_DIR, "tables", "table2_cw_vol_tci.csv"))
fwrite(tbl3, file.path(OUTPUT_DIR, "tables", "table3_cw_vol.csv"))

# ---------------------------------------------------------------------------
# Step 8: Plots
# ---------------------------------------------------------------------------
message("Generating figures...")
p1 <- plot_gold_data(gold_dt, data.table(date = dates_g7, vol = vol_g7))
ggplot2::ggsave(file.path(OUTPUT_DIR, "figures", "fig1_gold_data.pdf"), p1,
                width = 8, height = 6)

for (ctry in ALL_COUNTRIES) {
  ci_dt <- if (ctry %in% G7_COUNTRIES) ci_g7[[ctry]] else ci_all[[ctry]]
  pd_dt <- if (ctry %in% G7_COUNTRIES) pd_input_g7[[ctry]] else pd_input_all[[ctry]]
  p2 <- plot_ci_country(ci_dt, pd_dt$lnp_weekly, ctry)
  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "figures", paste0("fig2_ci_", ctry, ".pdf")),
    p2, width = 10, height = 8
  )
}

message("Done. Outputs written to ", OUTPUT_DIR)

# Restore sequential plan
plan(sequential)
