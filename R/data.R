# R/data.R — Data loading and preprocessing
# Gabauer, Gupta, Karmakar & Nielsen (2022) replication

library(data.table)

# ---------------------------------------------------------------------------
# load_gold: Read daily LBMA gold prices, compute weekly averages and log-returns
# ---------------------------------------------------------------------------
# Expected CSV columns: date (YYYY-MM-DD), price (USD per troy oz)
# Returns data.table with columns: date (weekly Sunday/end-of-week), price, log_return
load_gold <- function(path) {
  dt <- fread(path)

  # Normalise column names to lowercase
  setnames(dt, tolower(names(dt)))

  # Parse date
  if (!inherits(dt$date, "Date")) {
    dt[, date := as.Date(date)]
  }

  setkey(dt, date)

  # Identify ISO week and year for weekly averaging
  dt[, week_id := format(date, "%G-%V")]  # ISO year-week

  # Weekly average price (use last trading day of week as representative date)
  weekly <- dt[, .(
    date  = max(date),
    price = mean(price, na.rm = TRUE)
  ), by = week_id]

  setkey(weekly, date)
  weekly[, week_id := NULL]

  # Log-return: r_t = log(P_t / P_{t-1})
  weekly[, log_return := c(NA_real_, diff(log(price)))]

  weekly[]
}

# ---------------------------------------------------------------------------
# load_pd_ratios: Read daily price-dividend ratio for one country
# ---------------------------------------------------------------------------
# Expected CSV columns: date, pd_ratio  (or price and dividend separately)
# Returns list with:
#   $lnp_daily  — numeric vector of log(P/D) values
#   $all_dates  — Date vector aligned to lnp_daily
load_pd_ratios <- function(path, country) {
  dt <- fread(path)
  setnames(dt, tolower(names(dt)))

  if (!inherits(dt$date, "Date")) {
    dt[, date := as.Date(date)]
  }
  setkey(dt, date)

  # Compute log P/D ratio
  # Accept either a 'pd_ratio' column directly, or 'price'+'dividend' columns
  if ("pd_ratio" %in% names(dt)) {
    dt[, lnpd := log(pd_ratio)]
  } else if (all(c("price", "dividend") %in% names(dt))) {
    dt[, lnpd := log(price / dividend)]
  } else {
    stop(sprintf(
      "load_pd_ratios: file '%s' must have columns 'pd_ratio' or 'price'+'dividend'",
      path
    ))
  }

  # Remove NAs/Inf introduced by log of zero or negative
  dt <- dt[is.finite(lnpd)]

  list(
    lnp_daily = dt$lnpd,
    all_dates  = dt$date,
    country    = country
  )
}

# ---------------------------------------------------------------------------
# align_to_weekly: Subsample a daily series to a set of weekly target dates
# ---------------------------------------------------------------------------
# lnp_daily : numeric vector (daily log P/D)
# all_dates  : Date vector, same length as lnp_daily
# target_dates: Date vector of weekly dates to align to
#
# For each target date, uses the most recent available daily observation
# (last-observation-carried-forward).
#
# Returns a data.table with columns: date, lnp_weekly
align_to_weekly <- function(lnp_daily, all_dates, target_dates) {
  # Build lookup: date -> lnp value
  daily_dt <- data.table(date = all_dates, lnp = lnp_daily)
  setkey(daily_dt, date)

  # For each target date, find last available observation <= target
  result <- data.table(date = target_dates)
  result[, lnp_weekly := {
    idx <- findInterval(date, daily_dt$date)
    ifelse(idx >= 1L, daily_dt$lnp[idx], NA_real_)
  }]

  result[]
}

# ---------------------------------------------------------------------------
# build_predictor_matrix: Stack CI time series into wide matrix
# ---------------------------------------------------------------------------
# ci_list  : named list (country -> data.table with columns date, pos_short,
#            neg_short, pos_med, neg_med, pos_long, neg_long)
# dates_ref: Date vector — rows of the output matrix
#
# Returns a matrix (length(dates_ref) x 6*length(ci_list)) with column names
# {country}_{sign}_{scale}, e.g. "us_pos_short", "canada_neg_long".
build_predictor_matrix <- function(ci_list, dates_ref) {
  ci_cols <- c("pos_short", "neg_short", "pos_med", "neg_med", "pos_long", "neg_long")

  parts <- lapply(names(ci_list), function(ctry) {
    ci_dt <- ci_list[[ctry]]
    setkey(ci_dt, date)

    # Align to dates_ref
    aligned <- ci_dt[.(dates_ref), on = "date"]

    mat <- as.matrix(aligned[, ..ci_cols])
    colnames(mat) <- paste(ctry, ci_cols, sep = "_")
    mat
  })

  X <- do.call(cbind, parts)

  # Replace NAs with 0 (no signal before series begins)
  X[is.na(X)] <- 0

  X
}
