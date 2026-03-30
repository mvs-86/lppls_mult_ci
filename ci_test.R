library(future)
library(furrr)
library(data.table)

# ---------------------------------------------------------------------------
# Parallel plan — adjust workers to your machine
# ---------------------------------------------------------------------------
plan(multisession, workers = min(12L, parallelly::availableCores() - 1L))
# plan(multisession, workers = 12)

# ---------------------------------------------------------------------------
# Source all modules
# ---------------------------------------------------------------------------
# source("R/fetch_data.R")
# source("R/data.R")
source("R/lppls_fit.R")
source("R/lppls_ci.R")
source("R/ci_batch.R")

quantmod::getSymbols("^IXIC", from="1997-1-1", to="2002-1-1")

prices <- quantmod::Ad(IXIC["/2001-12"])
dates <- zoo::index(IXIC["/2001-12"])

# tictoc::tic("COMPUTE CI");ci <- compute_ci_asset(as.numeric(prices), dates); tictoc::toc()
tictoc::tic("COMPUTE CI");ci <- compute_ci_batched(as.numeric(prices), dates); tictoc::toc()
ci <- merge_ci_batches("cache/ci_batches")

source("R/plots.R")

gg <- plot_ci_asset(ci, prices, dates, asset_name = "IXIC")
gg
