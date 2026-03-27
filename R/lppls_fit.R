# R/lppls_fit.R — LPPLS single-window fitting
# Filimonov-Sornette (2013) reparameterisation with SLSQP via nloptr
# Gabauer, Gupta, Karmakar & Nielsen (2022) replication

library(nloptr)

# ---------------------------------------------------------------------------
# lppls_basis: Compute basis vectors f, g, h for given nonlinear parameters
# ---------------------------------------------------------------------------
# t     : numeric vector of time indices (equally spaced integers or trading days)
# tc    : critical time
# m     : power-law exponent
# omega : log-periodic frequency
#
# Returns list(f, g, h) or NULL if any tc - t <= 0 (undefined power law)
lppls_basis <- function(t, tc, m, omega) {
  dt_vec <- tc - t
  if (any(dt_vec <= 0)) return(NULL)

  dt_m <- dt_vec^m
  log_dt <- log(dt_vec)

  list(
    f = dt_m,
    g = dt_m * cos(omega * log_dt),
    h = dt_m * sin(omega * log_dt)
  )
}

# ---------------------------------------------------------------------------
# lppls_linear_solve: Solve 4x4 normal equations for {A, B, C1, C2}
# ---------------------------------------------------------------------------
# f, g, h : basis vectors (numeric, same length as lnp)
# lnp     : log price-dividend ratio (numeric vector)
#
# Returns list(A, B, C1, C2, ssr) or NULL on singular matrix
lppls_linear_solve <- function(f, g, h, lnp) {
  N  <- length(lnp)
  sf <- sum(f);   sg <- sum(g);   sh <- sum(h)
  sf2 <- sum(f*f); sg2 <- sum(g*g); sh2 <- sum(h*h)
  sfg <- sum(f*g); sfh <- sum(f*h); sgh <- sum(g*h)
  slnp <- sum(lnp)
  sflnp <- sum(f*lnp); sglnp <- sum(g*lnp); shlnp <- sum(h*lnp)

  M <- matrix(c(
    N,   sf,  sg,  sh,
    sf,  sf2, sfg, sfh,
    sg,  sfg, sg2, sgh,
    sh,  sfh, sgh, sh2
  ), nrow = 4, byrow = TRUE)

  rhs <- c(slnp, sflnp, sglnp, shlnp)

  params <- tryCatch(
    solve(M, rhs),
    error = function(e) NULL
  )
  if (is.null(params)) return(NULL)

  A <- params[1]; B <- params[2]; C1 <- params[3]; C2 <- params[4]

  fitted <- A + B*f + C1*g + C2*h
  ssr    <- sum((lnp - fitted)^2)

  list(A = A, B = B, C1 = C1, C2 = C2, ssr = ssr)
}

# ---------------------------------------------------------------------------
# lppls_cost: Concentrated cost function F1(tc, m, omega) -> SSR
# ---------------------------------------------------------------------------
lppls_cost <- function(nonlin_params, t, lnp) {
  tc <- nonlin_params[1]
  m  <- nonlin_params[2]
  om <- nonlin_params[3]

  basis <- lppls_basis(t, tc, m, om)
  if (is.null(basis)) return(1e10)

  lin <- lppls_linear_solve(basis$f, basis$g, basis$h, lnp)
  if (is.null(lin)) return(1e10)

  lin$ssr
}

# ---------------------------------------------------------------------------
# lppls_gradient: Numerical gradient of concentrated cost (central differences)
# ---------------------------------------------------------------------------
lppls_gradient <- function(nonlin_params, t, lnp, eps = 1e-5) {
  grad <- numeric(3)
  for (i in seq_along(nonlin_params)) {
    p_up <- p_dn <- nonlin_params
    p_up[i] <- nonlin_params[i] + eps
    p_dn[i] <- nonlin_params[i] - eps
    grad[i] <- (lppls_cost(p_up, t, lnp) - lppls_cost(p_dn, t, lnp)) / (2 * eps)
  }
  grad
}

# ---------------------------------------------------------------------------
# lppls_grid_search: Vectorised grid scan over (tc, m, omega) → best starting point
# ---------------------------------------------------------------------------
# Sweeps a coarse grid, vectorising over omega for each (tc, m) pair so that
# basis computations and all dot-products are done with BLAS matrix ops rather
# than repeated scalar R calls.  Returns the top `top_k` candidate parameter
# vectors ordered by SSR (cheapest first).
#
# t, lnp  : window data (normalised, length N)
# tc_n    : number of tc grid points in [tc_lb, tc_ub]
# m_n     : number of m  grid points in (0.01, 0.99)
# omega_n : number of ω  grid points in (2, 15)
# top_k   : number of candidates to return (for multi-start refinement)
#
# Returns a matrix (top_k × 3) of [tc, m, omega] rows, best first.
lppls_grid_search <- function(t, lnp, tc_lb, tc_ub,
                              tc_n = 20L, m_n = 12L, omega_n = 12L,
                              top_k = 3L) {
  tc_grid    <- seq(tc_lb, tc_ub, length.out = tc_n)
  m_grid     <- seq(0.01,  0.99,  length.out = m_n)
  omega_grid <- seq(2,     15,    length.out = omega_n)

  N      <- length(lnp)
  ones_N <- rep(1.0, N)
  # Pre-compute constants (don't depend on omega)
  slnp  <- sum(lnp)
  N_val <- N

  best_ssrs   <- rep(Inf, top_k)
  best_params <- matrix(NA_real_, nrow = top_k, ncol = 3L)

  for (tc in tc_grid) {
    dt_vec <- tc - t
    if (any(dt_vec <= 0)) next

    for (m in m_grid) {
      f      <- dt_vec^m
      log_dt <- log(dt_vec)

      # Pre-compute constants that don't depend on omega
      sf    <- sum(f)
      sf2   <- sum(f * f)
      sflnp <- sum(f * lnp)

      # Vectorise over omega: G (omega_n × N), H (omega_n × N)
      # outer product: row j = omega_grid[j] * log_dt
      om_ld <- outer(omega_grid, log_dt)  # omega_n × N
      G <- cos(om_ld) * matrix(f, nrow = omega_n, ncol = N, byrow = TRUE)
      H <- sin(om_ld) * matrix(f, nrow = omega_n, ncol = N, byrow = TRUE)

      # Vectorised dot-products (each returns omega_n-length vector)
      sg    <- drop(G %*% ones_N)
      sh    <- drop(H %*% ones_N)
      sfg   <- drop(G %*% f)
      sfh   <- drop(H %*% f)
      sg2   <- rowSums(G * G)
      sh2   <- rowSums(H * H)
      sgh   <- rowSums(G * H)
      sglnp <- drop(G %*% lnp)
      shlnp <- drop(H %*% lnp)

      # Solve 4×4 for each omega and compute SSR
      for (j in seq_len(omega_n)) {
        M <- matrix(c(
          N_val,  sf,      sg[j],   sh[j],
          sf,     sf2,     sfg[j],  sfh[j],
          sg[j],  sfg[j],  sg2[j],  sgh[j],
          sh[j],  sfh[j],  sgh[j],  sh2[j]
        ), nrow = 4L, byrow = TRUE)
        rhs <- c(slnp, sflnp, sglnp[j], shlnp[j])

        params <- tryCatch(solve(M, rhs), error = function(e) NULL)
        if (is.null(params)) next

        # SSR
        fitted <- params[1] + params[2]*f + params[3]*G[j,] + params[4]*H[j,]
        ssr    <- sum((lnp - fitted)^2)

        if (!is.finite(ssr)) next

        # Keep top_k by insertion
        worst_idx <- which.max(best_ssrs)
        if (ssr < best_ssrs[worst_idx]) {
          best_ssrs[worst_idx]      <- ssr
          best_params[worst_idx, ]  <- c(tc, m, omega_grid[j])
        }
      }
    }
  }

  # Return sorted by SSR (best first), dropping NA rows
  valid <- which(!is.na(best_params[, 1]) & is.finite(best_ssrs))
  if (length(valid) == 0L) return(NULL)
  ord <- valid[order(best_ssrs[valid])]
  best_params[ord, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# lppls_fit: Fit LPPLS to a window [t1, t2] via SLSQP with multiple starts
# ---------------------------------------------------------------------------
# t         : integer vector of time indices (length N)
# lnp       : numeric vector of log P/D ratio (length N)
# n_starts  : number of starting points (random starts, or top-k grid candidates)
# method    : "random" (original multi-start SLSQP) |
#             "grid"   (vectorised grid scan → top-k refinement, default) |
#             "grid_only" (grid scan, no SLSQP refinement — fastest, least accurate)
#
# nloptr bounds (pre-filter; hard filter applied post-fit):
#   tc in [t2 - 0.5*dt, t2 + 0.5*dt]
#   m  in [0.01, 0.99]
#   om in [2, 15]
#
# Returns list(tc, m, omega, A, B, C1, C2, ssr, converged) or NULL on failure
lppls_fit <- function(t, lnp, n_starts = 3L,
                      method = c("grid", "random", "grid_only")) {
  method <- match.arg(method)

  t1 <- min(t)
  t2 <- max(t)
  dt <- t2 - t1

  tc_lb <- t2 - 0.5 * dt
  tc_ub <- t2 + 0.5 * dt
  lb <- c(tc_lb, 0.01,  2)
  ub <- c(tc_ub, 0.99, 15)

  # ---- Build candidate starting points ----
  if (method == "random") {
    set.seed(NULL)
    starts <- do.call(rbind, lapply(seq_len(n_starts), function(s)
      c(runif(1, tc_lb, tc_ub), runif(1, 0.01, 0.99), runif(1, 2, 15))
    ))
  } else {
    # Grid search: vectorised over omega; returns top n_starts candidates
    starts <- lppls_grid_search(t, lnp, tc_lb, tc_ub, top_k = n_starts)
    if (is.null(starts)) return(NULL)

    if (method == "grid_only") {
      # Accept the best grid point without SLSQP refinement
      x0    <- starts[1, ]
      basis <- lppls_basis(t, x0[1], x0[2], x0[3])
      if (is.null(basis)) return(NULL)
      lin <- lppls_linear_solve(basis$f, basis$g, basis$h, lnp)
      if (is.null(lin)) return(NULL)
      return(list(tc = x0[1], m = x0[2], omega = x0[3],
                  A = lin$A, B = lin$B, C1 = lin$C1, C2 = lin$C2,
                  ssr = lin$ssr, converged = FALSE))
    }
  }

  # ---- SLSQP refinement from each candidate ----
  best_ssr <- Inf
  best_fit <- NULL
  best_lin <- NULL

  for (s in seq_len(nrow(starts))) {
    x0  <- starts[s, ]

    res <- tryCatch(
      nloptr(
        x0     = x0,
        eval_f = function(x) list(
          objective = lppls_cost(x, t, lnp),
          gradient  = lppls_gradient(x, t, lnp)
        ),
        lb   = lb,
        ub   = ub,
        opts = list(
          algorithm = "NLOPT_LD_SLSQP",
          maxeval   = 200L,   # fewer iters needed from a good starting point
          xtol_rel  = 1e-6,
          ftol_rel  = 1e-8
        )
      ),
      error = function(e) NULL
    )

    if (is.null(res) || res$objective >= best_ssr) next

    basis <- lppls_basis(t, res$solution[1], res$solution[2], res$solution[3])
    if (is.null(basis)) next
    lin <- lppls_linear_solve(basis$f, basis$g, basis$h, lnp)
    if (is.null(lin)) next

    best_ssr <- res$objective
    best_fit <- res
    best_lin <- lin
  }

  if (is.null(best_fit)) return(NULL)

  list(
    tc        = best_fit$solution[1],
    m         = best_fit$solution[2],
    omega     = best_fit$solution[3],
    A         = best_lin$A,
    B         = best_lin$B,
    C1        = best_lin$C1,
    C2        = best_lin$C2,
    ssr       = best_ssr,
    converged = (best_fit$status > 0)
  )
}
