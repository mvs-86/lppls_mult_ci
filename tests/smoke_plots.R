setwd("C:/Users/marcus/Projetos/lppls_mult_ci")
source("R/plots.R")

set.seed(1)
n <- 100
dates <- seq(as.Date("2010-01-01"), by = "week", length.out = n)
lnp   <- cumsum(rnorm(n, sd = 0.02)) + 3

# --- Standard 3-scale CI ---
ci_dt <- data.table(
  date      = dates,
  pos_short = runif(n), neg_short = runif(n),
  pos_med   = runif(n), neg_med   = runif(n),
  pos_long  = runif(n), neg_long  = runif(n)
)

cat("Testing plot_ci_asset (with lnp)... ")
p1 <- plot_ci_asset(ci_dt, lnp_series = lnp, asset_name = "BTC")
ggplot2::ggsave("output/figures/smoke_ci_btc.pdf", p1, width = 10, height = 8)
cat("OK\n")

cat("Testing plot_ci_asset (no lnp, no title)... ")
p2 <- plot_ci_asset(ci_dt)
cat("OK\n")

cat("Testing plot_ci_country... ")
p3 <- plot_ci_country(ci_dt, lnp, "south_africa")
cat("OK (title: SOUTH AFRICA)\n")

# --- Custom single-scale CI ---
ci_custom <- data.table(date = dates, pos_fast = runif(n), neg_fast = runif(n))
cat("Testing plot_ci_asset (custom scale 'fast')... ")
p4 <- plot_ci_asset(ci_custom, lnp_series = lnp, asset_name = "ETH")
cat("OK\n")

# --- plot_gold_data ---
gold_dt <- data.table(
  date       = dates,
  price      = exp(lnp),
  log_return = c(NA_real_, diff(lnp))
)
vol_dt <- data.table(date = dates, vol = abs(rnorm(n, mean = 0.005, sd = 0.002)))

cat("Testing plot_gold_data... ")
p5 <- plot_gold_data(gold_dt, vol_dt)
ggplot2::ggsave("output/figures/smoke_gold.pdf", p5, width = 8, height = 6)
cat("OK\n")

# --- Error: no pos_* columns ---
cat("Testing error on missing pos_* columns... ")
bad_dt <- data.table(date = dates, ci = runif(n))
tryCatch(
  plot_ci_asset(bad_dt),
  error = function(e) cat("caught:", conditionMessage(e), "\n")
)

# --- plot_ci_price ---
cat("Testing plot_ci_price (explicit scale)... ")
p6 <- plot_ci_price(ci_dt, exp(lnp), dates_price = dates,
                    scale = "short", asset_name = "BTC")
ggplot2::ggsave("output/figures/smoke_ci_price_bars.pdf", p6, width = 10, height = 5)
cat("OK\n")

cat("Testing plot_ci_price (auto-scale)... ")
p7 <- plot_ci_price(ci_dt, exp(lnp))
cat("OK\n")

cat("Testing plot_ci_price error on unknown scale... ")
tryCatch(
  plot_ci_price(ci_dt, exp(lnp), scale = "xxx"),
  error = function(e) cat("caught:", conditionMessage(e), "\n")
)

cat("\nAll smoke tests passed.\n")
