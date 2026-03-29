# R/lppls_ci.R — Multi-scale LPPLS Confidence Indicator
# Gabauer, Gupta, Karmakar & Nielsen (2022) replication

library(data.table)
library(furrr)
library(future.apply)  # future_lapply for parallel t2 loop

# Window-length ranges (trading days) and step size (spec §2.1)
SCALES <- list(
  short  = list(dt_min = 30L,  dt_max = 90L,  step = 2L),
  medium = list(dt_min = 90L,  dt_max = 300L, step = 2L),
  long   = list(dt_min = 300L, dt_max = 745L, step = 2L)
)

# Sampling step within each window (spec: 5 trading days)
WINDOW_SAMPLE_STEP <- 5L

# ---------------------------------------------------------------------------
# lppls_filter: Apply all 5 filter conditions to a fitted LPPLS result
# ---------------------------------------------------------------------------
# fit  : list from lppls_fit() — must have fields tc, m, omega, B, C1, C2
# t1   : integer, start index of window
# t2   : integer, end index of window
#
# Returns TRUE if the fit is "qualified", FALSE otherwise.
lppls_filter <- function(fit, t1, t2) {
  if (is.null(fit)) return(FALSE)

  tc <- fit$tc
  m  <- fit$m
  om <- fit$omega
  B  <- fit$B
  C1 <- fit$C1
  C2 <- fit$C2
  dt <- t2 - t1

  # 1. m bounds
  if (m  < 0.01 || m  > 0.99) return(FALSE)

  # 2. omega bounds
  if (om < 2    || om > 15)   return(FALSE)

  # 3. tc bounds: tc in [max(t2-60, t2-0.5*dt), t2 + min(252, 0.5*dt)]
  tc_lb <- max(t2 - 60L, t2 - 0.5 * dt)
  tc_ub <- t2 + min(252L, 0.5 * dt)
  if (tc < tc_lb || tc > tc_ub) return(FALSE)

  # 4. Oscillation count O > 2.5
  # O = (omega / 2pi) * log( (tc - t1) / (tc - t2) )
  denom_O <- tc - t2
  if (denom_O <= 0) return(FALSE)
  O <- (om / (2 * pi)) * log((tc - t1) / denom_O)
  if (!is.finite(O) || O <= 2.5) return(FALSE)

  # 5. Damping ratio D > 0.5
  # D = m * |B| / (omega * |C|),  |C| = sqrt(C1^2 + C2^2)
  absC <- sqrt(C1^2 + C2^2)
  if (absC == 0) return(FALSE)
  D <- (m * abs(B)) / (om * absC)
  if (!is.finite(D) || D <= 0.5) return(FALSE)

  TRUE
}

# ---------------------------------------------------------------------------
# .tc_to_dates: Convert mean tc indices (into the full lnp_vec) to Dates
# ---------------------------------------------------------------------------
# tc_vec : numeric vector of mean tc index values (may exceed length(dates))
# dates  : Date vector aligned to lnp_vec (index 1 = dates[1])
#
# Indices within [1, n] are looked up directly; values outside are
# linearly extrapolated using the last date-step.
.tc_to_dates <- function(tc_vec, dates) {
  n <- length(dates)
  if (n < 2L) return(rep(as.Date(NA), length(tc_vec)))
  step <- as.numeric(dates[n] - dates[n - 1L])
  idx  <- round(tc_vec)
  as.Date(
    ifelse(is.na(idx),
      NA_real_,
      ifelse(idx < 1L,
        as.numeric(dates[1L]) + (idx - 1L) * step,
        ifelse(idx > n,
          as.numeric(dates[n]) + (idx - n) * step,
          as.numeric(dates[pmax(1L, pmin(n, idx))])
        )
      )
    ),
    origin = "1970-01-01"
  )
}

# ---------------------------------------------------------------------------
# compute_ci_t2: Compute CI (pos and neg) for one t2 at a given scale
# ---------------------------------------------------------------------------
# lnp_daily  : full numeric vector of daily log P/D (index 1..N)
# t2_idx     : integer, index in lnp_daily of the current week-end
# dt_min, dt_max, step: window length range and iteration step
#
# Returns list(pos, neg, tc_pos, tc_neg):
#   pos/neg    — fraction of qualified fits with B<0 (positive) / B>=0 (negative)
#   tc_pos/neg — mean tc index (in lnp_daily units) of those qualified fits;
#                NA when no qualified fits of that sign
compute_ci_t2 <- function(lnp_daily, t2_idx, dt_min, dt_max, step,
                          n_starts = 3L, method = "grid") {
  n_qualified_pos <- 0L
  n_qualified_neg <- 0L
  n_total         <- 0L
  tc_pos_vals     <- numeric(0)
  tc_neg_vals     <- numeric(0)

  dt_seq <- seq(dt_min, dt_max, by = step)

  for (dt in dt_seq) {
    t1_idx <- t2_idx - dt
    if (t1_idx < 1L) next

    # Subsample window at every WINDOW_SAMPLE_STEP trading days
    idx_window <- seq(t1_idx, t2_idx, by = WINDOW_SAMPLE_STEP)
    if (length(idx_window) < 5L) next  # need at least 5 points for a 4-param fit

    lnp_window <- lnp_daily[idx_window]
    if (any(!is.finite(lnp_window))) next

    # Normalise t to integers starting from 1 for numerical stability
    t_window <- seq_along(idx_window)
    t1_norm  <- 1L
    t2_norm  <- length(t_window)

    fit <- lppls_fit(t_window, lnp_window, n_starts = n_starts, method = method)

    n_total <- n_total + 1L

    if (lppls_filter(fit, t1_norm, t2_norm)) {
      # Map normalized tc back to original index space:
      #   tc_norm=1 → t1_idx,  tc_norm=t2_norm → t2_idx (approx)
      tc_actual <- t1_idx + (fit$tc - 1L) * WINDOW_SAMPLE_STEP

      if (fit$B < 0) {
        n_qualified_pos <- n_qualified_pos + 1L
        tc_pos_vals     <- c(tc_pos_vals, tc_actual)
      } else {
        n_qualified_neg <- n_qualified_neg + 1L
        tc_neg_vals     <- c(tc_neg_vals, tc_actual)
      }
    }
  }

  if (n_total == 0L) {
    return(list(pos = 0, neg = 0, tc_pos = NA_real_, tc_neg = NA_real_))
  }

  list(
    pos    = n_qualified_pos / n_total,
    neg    = n_qualified_neg / n_total,
    tc_pos = if (length(tc_pos_vals) > 0L) mean(tc_pos_vals) else NA_real_,
    tc_neg = if (length(tc_neg_vals) > 0L) mean(tc_neg_vals) else NA_real_
  )
}

# ---------------------------------------------------------------------------
# compute_ci_scale: Compute CI for all t2 values at one scale (sequential inner)
# ---------------------------------------------------------------------------
# lnp_daily  : full numeric vector
# t2_indices : integer vector of week-end indices into lnp_daily
# scale_name : character — "short", "medium", or "long"
#
# Returns data.table(t2_idx, ci_pos, ci_neg)
compute_ci_scale <- function(lnp_daily, t2_indices, scale_name,
                             n_starts = 3L, method = "grid") {
  sc   <- SCALES[[scale_name]]
  dt_min <- sc$dt_min
  dt_max <- sc$dt_max
  step   <- sc$step

  results <- lapply(t2_indices, function(t2_idx) {
    ci <- compute_ci_t2(lnp_daily, t2_idx, dt_min, dt_max, step,
                        n_starts = n_starts, method = method)
    list(t2_idx = t2_idx, ci_pos = ci$pos, ci_neg = ci$neg,
         tc_pos = ci$tc_pos, tc_neg = ci$tc_neg)
  })

  rbindlist(results)
}

# ---------------------------------------------------------------------------
# compute_ci_asset: Multi-scale CI for any log-price series
# ---------------------------------------------------------------------------
# Generic version — works for equities, commodities, crypto, etc.
#
# x         : numeric vector of log prices  OR  data.table with a date column
#             and exactly one numeric value column (or specify value_col).
# dates     : Date vector — required when x is a numeric vector; ignored otherwise.
# value_col : name of the value column when x is a data.table.
#             Auto-detected as the first non-date numeric column if NULL.
# scales    : named list of scale definitions, each with dt_min, dt_max, step.
#             Defaults to SCALES (spec §2.1: short/medium/long for daily equity data).
#             Pass custom scales for assets sampled at a different frequency, e.g.:
#               crypto_scales <- list(
#                 short  = list(dt_min = 20L,  dt_max =  60L, step = 2L),
#                 medium = list(dt_min = 60L,  dt_max = 200L, step = 2L),
#                 long   = list(dt_min = 200L, dt_max = 500L, step = 2L)
#               )
#
# Returns data.table(date, pos_{scale}, neg_{scale}) for every scale in `scales`.
# Column names follow the pattern pos_short / neg_short / pos_med / neg_med etc.,
# using the names() of the `scales` list (abbreviated to first 5 chars if needed).
compute_ci_asset <- function(x, dates = NULL, value_col = NULL,
                             scales = SCALES, n_starts = 3L,
                             method = "grid", parallel_t2 = TRUE) {
  # --- Resolve input to (lnp_vec, dates) ---
  if (is.data.table(x) || is.data.frame(x)) {
    dt_in <- as.data.table(x)

    # Detect value column
    if (is.null(value_col)) {
      num_cols <- setdiff(names(dt_in)[vapply(dt_in, is.numeric, logical(1))], "date")
      if (length(num_cols) == 0L)
        stop("compute_ci_asset: no numeric column found in input data.table")
      if (length(num_cols) > 1L)
        stop("compute_ci_asset: multiple numeric columns found — specify value_col: ",
             paste(num_cols, collapse = ", "))
      value_col <- num_cols[1L]
    }
    if (!value_col %in% names(dt_in))
      stop("compute_ci_asset: column '", value_col, "' not found")

    lnp_vec <- dt_in[[value_col]]
    dates   <- dt_in[["date"]]
    if (is.null(dates))
      stop("compute_ci_asset: data.table input must have a 'date' column")

  } else if (is.numeric(x)) {
    lnp_vec <- x
    if (is.null(dates))
      stop("compute_ci_asset: 'dates' must be supplied when x is a numeric vector")
    if (length(dates) != length(lnp_vec))
      stop("compute_ci_asset: length(dates) != length(x)")

  } else {
    stop("compute_ci_asset: x must be a numeric vector or a data.table")
  }

  # --- Validate scales ---
  if (length(scales) == 0L) stop("compute_ci_asset: 'scales' must be a non-empty list")
  for (sc_name in names(scales)) {
    sc <- scales[[sc_name]]
    if (!all(c("dt_min", "dt_max", "step") %in% names(sc)))
      stop("compute_ci_asset: scale '", sc_name,
           "' must have fields dt_min, dt_max, step")
  }

  # --- Compute CI for each scale ---
  t2_indices <- which(is.finite(lnp_vec))

  ci_list <- lapply(names(scales), function(sc_name) {
    sc <- scales[[sc_name]]
    compute_ci_scale_custom(lnp_vec, t2_indices,
                            dt_min      = sc$dt_min,
                            dt_max      = sc$dt_max,
                            step        = sc$step,
                            n_starts    = n_starts,
                            parallel_t2 = parallel_t2,
                            method      = method)
  })
  names(ci_list) <- names(scales)

  # --- Build output data.table ---
  # Column name: pos_{scale_name}, neg_{scale_name}
  # Abbreviate "medium" → "med" to stay consistent with existing convention
  scale_abbrev <- function(nm) {
    switch(nm, medium = "med", nm)
  }

  out <- data.table(date = dates[t2_indices])
  for (sc_name in names(scales)) {
    ab <- scale_abbrev(sc_name)
    out[, paste0("pos_",    ab) := ci_list[[sc_name]]$ci_pos]
    out[, paste0("neg_",    ab) := ci_list[[sc_name]]$ci_neg]
    out[, paste0("tc_pos_", ab) := .tc_to_dates(ci_list[[sc_name]]$tc_pos, dates)]
    out[, paste0("tc_neg_", ab) := .tc_to_dates(ci_list[[sc_name]]$tc_neg, dates)]
  }

  setkey(out, date)
  out
}

# ---------------------------------------------------------------------------
# compute_ci_scale_custom: Like compute_ci_scale but takes dt bounds directly
# ---------------------------------------------------------------------------
# parallel_t2 : if TRUE, use future_lapply to parallelise over t2 indices.
#               Set FALSE (default in compute_ci_assets) to avoid nested futures
#               when the outer loop already parallelises over assets.
#               Set TRUE (default in compute_ci_asset) for single-asset runs.
compute_ci_scale_custom <- function(lnp_vec, t2_indices, dt_min, dt_max, step,
                                    n_starts = 5L, parallel_t2 = FALSE,
                                    method = "grid") {
  worker <- function(t2_idx) {
    ci <- compute_ci_t2(lnp_vec, t2_idx, dt_min, dt_max, step,
                        n_starts = n_starts, method = method)
    list(t2_idx = t2_idx, ci_pos = ci$pos, ci_neg = ci$neg,
         tc_pos = ci$tc_pos, tc_neg = ci$tc_neg)
  }

  results <- if (parallel_t2) {
    future_lapply(t2_indices, worker, future.seed = TRUE)
  } else {
    lapply(t2_indices, worker)
  }

  rbindlist(results)
}

# ---------------------------------------------------------------------------
# compute_ci_assets: CI for a named list of assets (parallelised via furrr)
# ---------------------------------------------------------------------------
# asset_list : named list — each element is either:
#              (a) a data.table(date, <value_col>), or
#              (b) a list(x = numeric, dates = Date)
# value_col  : passed through to compute_ci_asset (NULL = auto-detect)
# scales     : passed through to compute_ci_asset
#
# Returns named list of data.tables, one per asset.
compute_ci_assets <- function(asset_list, value_col = NULL,
                              scales = SCALES, n_starts = 3L,
                              method = "grid") {
  # parallel_t2 = FALSE: outer loop already parallelises over assets,
  # so the inner t2 loop must stay sequential to avoid nested futures.
  future_map(
    asset_list,
    function(asset) {
      if (is.list(asset) && !is.data.table(asset) && !is.data.frame(asset)) {
        compute_ci_asset(asset$x, dates = asset$dates,
                         value_col = value_col, scales = scales,
                         n_starts = n_starts, method = method,
                         parallel_t2 = FALSE)
      } else {
        compute_ci_asset(asset, value_col = value_col, scales = scales,
                         n_starts = n_starts, method = method,
                         parallel_t2 = FALSE)
      }
    },
    .options = furrr_options(seed = TRUE)
  )
}

# ---------------------------------------------------------------------------
# compute_ci_country: backward-compatible wrapper around compute_ci_asset
# ---------------------------------------------------------------------------
# lnp_input : data.table(date, lnp_weekly)
compute_ci_country <- function(lnp_input) {
  compute_ci_asset(lnp_input, value_col = "lnp_weekly", scales = SCALES)
}

# ---------------------------------------------------------------------------
# compute_ci_all: backward-compatible wrapper around compute_ci_assets
# ---------------------------------------------------------------------------
# pd_input_list : named list of data.tables, each with columns date, lnp_weekly
compute_ci_all <- function(pd_input_list) {
  compute_ci_assets(pd_input_list, value_col = "lnp_weekly", scales = SCALES)
}
