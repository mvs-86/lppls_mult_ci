# tests/test_lppls_fit.R
library(testthat)

# Source relative to project root (works when run via run_tests.R or test_dir)
proj_root <- tryCatch(
  rprojroot::find_root(rprojroot::has_file("main.R")),
  error = function(e) getwd()
)
if (!exists("lppls_basis")) source(file.path(proj_root, "R", "lppls_fit.R"))

# ---------------------------------------------------------------------------
# Helper: simulate an LPPLS series with known parameters
# ---------------------------------------------------------------------------
sim_lppls <- function(N = 80, tc = 100, m = 0.5, omega = 7,
                      A = 10, B = -0.5, C1 = 0.1, C2 = -0.05) {
  t   <- seq(1, N)
  dt_vec <- tc - t
  stopifnot(all(dt_vec > 0))
  dt_m <- dt_vec^m
  lnp  <- A + B * dt_m +
    C1 * dt_m * cos(omega * log(dt_vec)) +
    C2 * dt_m * sin(omega * log(dt_vec))
  list(t = t, lnp = lnp, tc = tc, m = m, omega = omega,
       A = A, B = B, C1 = C1, C2 = C2)
}

# ---------------------------------------------------------------------------
test_that("lppls_basis returns correct dimensions", {
  t  <- 1:50
  tc <- 60; m <- 0.5; om <- 6
  b  <- lppls_basis(t, tc, m, om)
  expect_equal(length(b$f), 50)
  expect_equal(length(b$g), 50)
  expect_equal(length(b$h), 50)
})

test_that("lppls_basis returns NULL when tc <= max(t)", {
  t  <- 1:50
  expect_null(lppls_basis(t, tc = 50, m = 0.5, om = 6))
  expect_null(lppls_basis(t, tc = 30, m = 0.5, om = 6))
})

test_that("lppls_linear_solve recovers known linear params", {
  sim <- sim_lppls(N = 80, tc = 100, m = 0.5, omega = 7,
                   A = 10, B = -0.5, C1 = 0.1, C2 = -0.05)
  b   <- lppls_basis(sim$t, sim$tc, sim$m, sim$omega)
  lin <- lppls_linear_solve(b$f, b$g, b$h, sim$lnp)

  expect_false(is.null(lin))
  expect_equal(lin$A,  sim$A,  tolerance = 1e-6)
  expect_equal(lin$B,  sim$B,  tolerance = 1e-6)
  expect_equal(lin$C1, sim$C1, tolerance = 1e-6)
  expect_equal(lin$C2, sim$C2, tolerance = 1e-6)
  expect_lt(lin$ssr, 1e-8)  # noiseless simulation → near-zero SSR
})

test_that("lppls_linear_solve matches brute-force lm()", {
  sim <- sim_lppls(N = 80, tc = 100, m = 0.5, omega = 7)
  b   <- lppls_basis(sim$t, sim$tc, sim$m, sim$omega)

  # Brute-force via OLS
  ols  <- lm(sim$lnp ~ b$f + b$g + b$h)
  lin  <- lppls_linear_solve(b$f, b$g, b$h, sim$lnp)

  expect_equal(lin$A,  unname(coef(ols)[1]), tolerance = 1e-8)
  expect_equal(lin$B,  unname(coef(ols)[2]), tolerance = 1e-8)
  expect_equal(lin$C1, unname(coef(ols)[3]), tolerance = 1e-8)
  expect_equal(lin$C2, unname(coef(ols)[4]), tolerance = 1e-8)
})

test_that("lppls_cost returns large value for bad params", {
  sim  <- sim_lppls(N = 80, tc = 100)
  # Use the true params — should give near-zero cost
  true_cost <- lppls_cost(c(100, 0.5, 7), sim$t, sim$lnp)
  bad_cost  <- lppls_cost(c(110, 0.9, 4), sim$t, sim$lnp)
  expect_lt(true_cost, 1e-6)
  expect_gt(bad_cost, true_cost)
})

test_that("lppls_fit recovers parameters on noiseless simulation", {
  set.seed(42)
  sim <- sim_lppls(N = 80, tc = 100, m = 0.5, omega = 7,
                   A = 10, B = -0.5, C1 = 0.1, C2 = -0.05)
  fit <- lppls_fit(sim$t, sim$lnp, n_starts = 10L)

  expect_false(is.null(fit))
  expect_equal(fit$tc,    sim$tc,    tolerance = 0.5)
  expect_equal(fit$m,     sim$m,     tolerance = 0.05)
  expect_equal(fit$omega, sim$omega, tolerance = 0.5)
  expect_equal(fit$B,     sim$B,     tolerance = 0.05)
  expect_lt(fit$ssr, 1e-4)
})

test_that("lppls_fit returns NULL for a series that is too short", {
  # 4 points — fewer than 5 (hard minimum in compute_ci_t2)
  t   <- 1:4
  lnp <- c(1.0, 1.1, 1.2, 1.15)
  # nloptr will likely fail or converge poorly; NULL is acceptable
  # We only check it doesn't crash
  expect_error(lppls_fit(t, lnp, n_starts = 2L), NA)
})
