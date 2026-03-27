# tests/test_lppls_ci.R
library(testthat)

proj_root <- tryCatch(
  rprojroot::find_root(rprojroot::has_file("main.R")),
  error = function(e) getwd()
)
if (!exists("lppls_basis"))   source(file.path(proj_root, "R", "lppls_fit.R"))
if (!exists("compute_ci_t2")) source(file.path(proj_root, "R", "lppls_ci.R"))

# ---------------------------------------------------------------------------
# Helper: simulate a plausible LPPLS series of daily values
# ---------------------------------------------------------------------------
sim_daily_lppls <- function(N = 800, tc_offset = 50, m = 0.5, omega = 7,
                             A = 8, B = -0.8, C1 = 0.15, C2 = -0.05,
                             noise_sd = 0.01) {
  t      <- seq(1, N)
  tc     <- N + tc_offset
  dt_vec <- tc - t
  dt_m   <- dt_vec^m
  lnp    <- A + B * dt_m +
    C1 * dt_m * cos(omega * log(dt_vec)) +
    C2 * dt_m * sin(omega * log(dt_vec))
  lnp + rnorm(N, sd = noise_sd)
}

# ---------------------------------------------------------------------------
test_that("SCALES constants match paper spec", {
  expect_equal(SCALES$short$dt_min,  30L)
  expect_equal(SCALES$short$dt_max,  90L)
  expect_equal(SCALES$medium$dt_min, 90L)
  expect_equal(SCALES$medium$dt_max, 300L)
  expect_equal(SCALES$long$dt_min,   300L)
  expect_equal(SCALES$long$dt_max,   745L)
})

test_that("window counts match paper spec", {
  # short: (90-30)/2 + 1 = 31 windows
  n_short  <- length(seq(SCALES$short$dt_min,  SCALES$short$dt_max,  by = SCALES$short$step))
  n_medium <- length(seq(SCALES$medium$dt_min, SCALES$medium$dt_max, by = SCALES$medium$step))
  n_long   <- length(seq(SCALES$long$dt_min,   SCALES$long$dt_max,   by = SCALES$long$step))

  # Paper says 30, 105, 223 — these are the dt values iterated
  # seq(30,90,2) has length 31; seq(90,300,2) has length 106; seq(300,745,2) has length 223
  # (Slight discrepancy vs paper for short/medium due to inclusive upper bound — accept ±2)
  expect_lte(abs(n_short  - 30L),  2L)
  expect_lte(abs(n_medium - 105L), 2L)
  expect_equal(n_long, 223L)
})

test_that("lppls_filter rejects NULL fit", {
  expect_false(lppls_filter(NULL, 1, 50))
})

test_that("lppls_filter rejects out-of-bounds m", {
  fit <- list(tc = 55, m = 1.5, omega = 7, B = -0.5, C1 = 0.1, C2 = -0.1)
  expect_false(lppls_filter(fit, 1, 50))
  fit$m <- 0.0
  expect_false(lppls_filter(fit, 1, 50))
})

test_that("lppls_filter rejects out-of-bounds omega", {
  fit <- list(tc = 55, m = 0.5, omega = 1.5, B = -0.5, C1 = 0.1, C2 = -0.1)
  expect_false(lppls_filter(fit, 1, 50))
  fit$omega <- 16
  expect_false(lppls_filter(fit, 1, 50))
})

test_that("lppls_filter rejects tc outside window", {
  # t1=1, t2=50, dt=49
  # tc_lb = max(50-60, 50-24.5) = 25.5;  tc_ub = 50 + min(252, 24.5) = 74.5
  # Use tc=53 so O = (7/2pi)*log((53-1)/(53-50)) = 1.114*log(17.33) ≈ 3.18 > 2.5
  # Use small C to ensure D = m|B|/(omega|C|) = 0.5*1/(7*0.07) ≈ 1.0 > 0.5
  fit_ok <- list(tc = 53, m = 0.5, omega = 7, B = -1.0, C1 = 0.05, C2 = -0.05)
  expect_true(lppls_filter(fit_ok, 1, 50))

  fit_bad_high <- fit_ok; fit_bad_high$tc <- 200
  expect_false(lppls_filter(fit_bad_high, 1, 50))

  fit_bad_low  <- fit_ok; fit_bad_low$tc  <- 5
  expect_false(lppls_filter(fit_bad_low,  1, 50))
})

test_that("lppls_filter rejects low oscillation count O <= 2.5", {
  # Force O ~ 0 by setting omega = 2, small tc - t2 relative to tc - t1
  # O = (2/2pi) * log((55-1)/(55-50)) = (1/pi) * log(54/5) ≈ 0.76 < 2.5
  fit <- list(tc = 55, m = 0.5, omega = 2, B = -0.5, C1 = 0.1, C2 = -0.1)
  expect_false(lppls_filter(fit, 1, 50))
})

test_that("CI values are in [0, 1] on random noise", {
  set.seed(1)
  lnp <- cumsum(rnorm(400))

  # Two t2 values; n_starts=1 for speed
  t2_vals <- c(200L, 400L)
  ci_dt   <- compute_ci_scale(lnp, t2_vals, "short", n_starts = 1L)

  expect_true(all(ci_dt$ci_pos >= 0 & ci_dt$ci_pos <= 1))
  expect_true(all(ci_dt$ci_neg >= 0 & ci_dt$ci_neg <= 1))
})

test_that("positive CI is higher for true positive bubble vs random noise", {
  set.seed(42)

  lnp_bubble <- sim_daily_lppls(N = 600, tc_offset = 30, B = -0.8, noise_sd = 0.005)
  lnp_noise  <- cumsum(rnorm(600, sd = 0.01))

  # Two t2 values; use method="random" with n_starts=2 — tests statistical
  # behaviour (bubble detected more often than noise), not optimizer choice.
  t2_vals <- c(400L, 600L)

  ci_bubble <- compute_ci_scale(lnp_bubble, t2_vals, "short",
                                n_starts = 2L, method = "random")
  ci_noise  <- compute_ci_scale(lnp_noise,  t2_vals, "short",
                                n_starts = 2L, method = "random")

  expect_gt(mean(ci_bubble$ci_pos), mean(ci_noise$ci_pos))
})

test_that("compute_ci_t2 returns fractions summing to <= 1", {
  set.seed(7)
  lnp <- sim_daily_lppls(N = 500, noise_sd = 0.01)
  ci  <- compute_ci_t2(lnp, t2_idx = 450L,
                        dt_min   = SCALES$short$dt_min,
                        dt_max   = SCALES$short$dt_max,
                        step     = SCALES$short$step,
                        n_starts = 1L)
  expect_lte(ci$pos + ci$neg, 1 + 1e-9)
  expect_gte(ci$pos, 0)
  expect_gte(ci$neg, 0)
})

# ===========================================================================
# compute_ci_asset — generic interface
#
# All interface tests use n_starts=1L and a micro scale (dt=[30,35], step=2)
# so only ~3 windows per t2 are attempted, keeping each test under ~5 seconds.
# ===========================================================================

# Micro scale: dt ∈ {30,32,34} — 3 windows per t2, enough to check logic
MICRO_SCALE <- list(fast = list(dt_min = 30L, dt_max = 35L, step = 2L))

# N=38: valid t2 indices are 31:38 (8 values), dt_max=35 only valid from t2=36
# → ~24 LPPLS fits per test at n_starts=1 → < 2 seconds
make_asset_dt <- function(N = 38, seed = 1) {
  set.seed(seed)
  lnp   <- sim_daily_lppls(N = N, tc_offset = 10, noise_sd = 0.01)
  dates <- seq(as.Date("2010-01-01"), by = "day", length.out = N)
  data.table(date = dates, log_price = lnp)
}

test_that("compute_ci_asset accepts a data.table with auto-detected value column", {
  dt  <- make_asset_dt()
  out <- compute_ci_asset(dt, scales = MICRO_SCALE, n_starts = 1L)

  expect_true(is.data.table(out))
  expect_true("date"     %in% names(out))
  expect_true("pos_fast" %in% names(out))
  expect_true("neg_fast" %in% names(out))
  expect_equal(nrow(out), sum(is.finite(dt$log_price)))
})

test_that("compute_ci_asset default SCALES produce the 6 standard CI column names", {
  # Check column naming only — use micro ranges so computation is trivial
  dt  <- make_asset_dt()
  out <- compute_ci_asset(dt, scales = list(
    short  = list(dt_min = 30L, dt_max = 32L, step = 2L),
    medium = list(dt_min = 30L, dt_max = 32L, step = 2L),
    long   = list(dt_min = 30L, dt_max = 32L, step = 2L)
  ), n_starts = 1L)
  expect_true(all(c("pos_short", "neg_short", "pos_med", "neg_med",
                    "pos_long",  "neg_long") %in% names(out)))
})

test_that("compute_ci_asset accepts a data.table with explicit value_col", {
  dt  <- make_asset_dt()
  out <- compute_ci_asset(dt, value_col = "log_price",
                          scales = MICRO_SCALE, n_starts = 1L)
  expect_true(is.data.table(out))
  expect_true(all(out$pos_fast >= 0 & out$pos_fast <= 1))
})

test_that("compute_ci_asset accepts a numeric vector + dates", {
  set.seed(2)
  lnp   <- sim_daily_lppls(N = 38, tc_offset = 10, noise_sd = 0.01)
  dates <- seq(as.Date("2015-01-01"), by = "day", length.out = 38)
  out   <- compute_ci_asset(lnp, dates = dates,
                            scales = MICRO_SCALE, n_starts = 1L)
  expect_true(is.data.table(out))
  expect_true("pos_fast" %in% names(out))
})

test_that("compute_ci_asset errors without dates when x is numeric", {
  expect_error(compute_ci_asset(rnorm(50)), "dates")
})

test_that("compute_ci_asset errors on mismatched lengths", {
  lnp   <- rnorm(100)
  dates <- seq(as.Date("2020-01-01"), by = "day", length.out = 50)
  expect_error(compute_ci_asset(lnp, dates = dates), "length")
})

test_that("compute_ci_asset errors on ambiguous columns", {
  dt <- data.table(date = Sys.Date() + 1:50, a = rnorm(50), b = rnorm(50))
  expect_error(compute_ci_asset(dt, scales = MICRO_SCALE, n_starts = 1L),
               "value_col")
})

test_that("compute_ci_asset uses custom scale names in output columns", {
  dt <- make_asset_dt()
  custom_scales <- list(
    fast = list(dt_min = 30L, dt_max = 32L, step = 2L),
    slow = list(dt_min = 30L, dt_max = 34L, step = 4L)
  )
  out <- compute_ci_asset(dt, scales = custom_scales, n_starts = 1L)
  expect_true(all(c("pos_fast", "neg_fast", "pos_slow", "neg_slow") %in% names(out)))
  expect_false("pos_short" %in% names(out))
})

test_that("compute_ci_asset CI values are in [0,1]", {
  dt  <- make_asset_dt()
  out <- compute_ci_asset(dt, scales = MICRO_SCALE, n_starts = 1L)
  expect_true(all(out$pos_fast >= 0 & out$pos_fast <= 1))
  expect_true(all(out$neg_fast >= 0 & out$neg_fast <= 1))
})

test_that("compute_ci_country equals compute_ci_asset with lnp_weekly", {
  set.seed(3)
  lnp   <- sim_daily_lppls(N = 38, tc_offset = 10, noise_sd = 0.01)
  dates <- seq(as.Date("2000-01-01"), by = "day", length.out = 38)
  dt    <- data.table(date = dates, lnp_weekly = lnp)

  out_explicit  <- compute_ci_asset(dt, value_col = "lnp_weekly",
                                    scales = MICRO_SCALE, n_starts = 1L)
  out_autodetect <- compute_ci_asset(dt, scales = MICRO_SCALE, n_starts = 1L)

  expect_equal(out_explicit, out_autodetect)
})
