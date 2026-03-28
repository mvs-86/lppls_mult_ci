# R/fetch_data.R — Download raw data from public sources
#
# Gold prices: LBMA AM fix from FRED (GOLDAMGBD228NLBM) — free, no API key.
#
# Price-dividend ratios: The paper uses Refinitiv Datastream (paid).
# This module provides a free proxy using NYSE/NASDAQ-listed country ETFs via
# Yahoo Finance. ETF prices + dividends → trailing-12-month P/D ratio.
#
# LIMITATIONS:
#   - ETF inception dates limit BRICS history (see COUNTRY_ETFS below).
#   - ETF prices are USD-denominated; the paper uses local currencies.
#   - ETF dividends are after withholding tax; raw index dividends are not.
#   - For a faithful replication, replace the ETF CSVs with Datastream exports
#     in the same format (date, pd_ratio).
#
# Output files (written to output_dir/):
#   gold_prices.csv          — columns: date, price
#   pd_ratios/{country}.csv  — columns: date, pd_ratio

library(quantmod)
library(data.table)

# ---------------------------------------------------------------------------
# ETF tickers used as country P/D proxies
# Inception dates noted where they limit coverage
# ---------------------------------------------------------------------------
COUNTRY_ETFS <- c(
  canada       = "EWC",   # iShares MSCI Canada, since 1996
  france       = "EWQ",   # iShares MSCI France, since 1996
  germany      = "EWG",   # iShares MSCI Germany, since 1996
  italy        = "EWI",   # iShares MSCI Italy, since 1996
  japan        = "EWJ",   # iShares MSCI Japan, since 1996
  uk           = "EWU",   # iShares MSCI United Kingdom, since 1996
  us           = "SPY",   # SPDR S&P 500 ETF, since 1993
  brazil       = "EWZ",   # iShares MSCI Brazil, since 2000 (< 1999 start)
  russia       = "RSX",   # VanEck Russia ETF, since 2007 (< 1999 start)
  india        = "INDA",  # iShares MSCI India, since 2012 (< 1999 start)
  china        = "MCHI",  # iShares MSCI China, since 2011 (< 1999 start)
  south_africa = "EZA"    # iShares MSCI South Africa, since 2003 (< 1999 start)
)

# ---------------------------------------------------------------------------
# fetch_gold_fred: Download daily LBMA gold AM fix from FRED
# ---------------------------------------------------------------------------
# Saves: {output_dir}/gold_prices.csv with columns date, price
# Returns the data.table invisibly.
fetch_gold_fred <- function(output_dir  = "data",
                            start_date  = "1970-01-01",
                            end_date    = "2020-09-30") {
  message("Fetching LBMA gold prices from FRED...")

  xts_gold <- tryCatch(
    getSymbols("GOLDAMGBD228NLBM", src = "FRED",
               from = start_date, to = end_date,
               auto.assign = FALSE),
    error = function(e) stop("Failed to download gold from FRED: ", conditionMessage(e))
  )

  dt <- data.table(
    date  = as.Date(index(xts_gold)),
    price = as.numeric(xts_gold[, 1])
  )

  # Drop NA rows (non-trading days returned as NA by FRED)
  dt <- dt[!is.na(price)]

  out_path <- file.path(output_dir, "gold_prices.csv")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(dt, out_path)
  message(sprintf("  Saved %d rows to %s", nrow(dt), out_path))

  invisible(dt)
}

# ---------------------------------------------------------------------------
# fetch_pd_yahoo: P/D ratio for one country via ETF price + dividends
# ---------------------------------------------------------------------------
# Computes pd_ratio = Adj.Close / TTM_dividends (trailing 12-month sum).
# Rows with zero TTM dividends are dropped.
#
# Saves: {output_dir}/{country}.csv with columns date, pd_ratio
# Returns the data.table invisibly.
fetch_pd_yahoo <- function(country,
                           output_dir  = "data/pd_ratios",
                           start_date  = "1970-01-01",
                           end_date    = "2020-09-30") {
  ticker <- COUNTRY_ETFS[country]
  if (is.na(ticker)) stop("Unknown country: ", country)

  message(sprintf("Fetching P/D for %s (%s) from Yahoo Finance...", country, ticker))

  # Price data
  prices <- tryCatch(
    getSymbols(ticker, src = "yahoo",
               from = start_date, to = end_date,
               auto.assign = FALSE, warnings = FALSE),
    error = function(e) {
      warning(sprintf("  Could not download %s: %s", ticker, conditionMessage(e)))
      return(NULL)
    }
  )
  if (is.null(prices)) return(invisible(NULL))

  price_dt <- data.table(
    date  = as.Date(index(prices)),
    price = as.numeric(Ad(prices))  # adjusted close (split- and dividend-adjusted)
  )
  price_dt <- price_dt[!is.na(price)]

  # Dividend data
  divs <- tryCatch(
    getDividends(ticker, src = "yahoo",
                 from = start_date, to = end_date,
                 auto.assign = FALSE),
    error = function(e) {
      warning(sprintf("  Could not download dividends for %s: %s", ticker, conditionMessage(e)))
      return(xts())
    }
  )

  if (nrow(divs) == 0) {
    warning(sprintf("  No dividend data for %s — skipping P/D computation.", ticker))
    return(invisible(NULL))
  }

  div_dt <- data.table(
    date     = as.Date(index(divs)),
    dividend = as.numeric(divs[, 1])
  )

  # Trailing 12-month dividends: rolling sum over the prior ~252 trading days
  # Merge dividend events onto the full price date grid (fill zero on non-ex-div days)
  merged <- merge(price_dt, div_dt, by = "date", all.x = TRUE)
  merged[is.na(dividend), dividend := 0]
  setorder(merged, date)

  # Rolling sum over 252-day window
  merged[, ttm_div := frollsum(dividend, n = 252L, align = "right", na.rm = TRUE)]

  # Drop rows where TTM dividends are zero (no dividend history yet) or price is NA
  merged <- merged[ttm_div > 0 & is.finite(price)]

  if (nrow(merged) == 0) {
    warning(sprintf("  No valid P/D observations for %s after filtering.", country))
    return(invisible(NULL))
  }

  merged[, pd_ratio := price / ttm_div]
  out_dt <- merged[, .(date, pd_ratio)]

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(output_dir, paste0(country, ".csv"))
  fwrite(out_dt, out_path)
  message(sprintf("  Saved %d rows to %s (ETF: %s)", nrow(out_dt), out_path, ticker))

  invisible(out_dt)
}

# ---------------------------------------------------------------------------
# fetch_all_data: Download gold + all country P/D ratios
# ---------------------------------------------------------------------------
# countries  : character vector, subset of names(COUNTRY_ETFS)
# output_dir : root directory; gold saved here, P/D in {output_dir}/pd_ratios/
# start_date, end_date: date range (character "YYYY-MM-DD")
#
# Returns named list with $gold and $pd (named list of data.tables or NULLs).
fetch_all_data <- function(countries  = names(COUNTRY_ETFS),
                           output_dir = "data",
                           start_date = "1970-01-01",
                           end_date   = "2020-09-30") {
  gold <- fetch_gold_fred(output_dir = output_dir,
                          start_date = start_date,
                          end_date   = end_date)

  pd_list <- lapply(
    setNames(countries, countries),
    function(ctry) {
      fetch_pd_yahoo(ctry,
                     output_dir = file.path(output_dir, "pd_ratios"),
                     start_date = start_date,
                     end_date   = end_date)
    }
  )

  # Report coverage gaps for BRICS
  brics <- c("brazil", "russia", "india", "china", "south_africa")
  brics_gaps <- vapply(intersect(countries, brics), function(ctry) {
    dt <- pd_list[[ctry]]
    if (is.null(dt)) return("NO DATA")
    min_yr <- format(min(dt$date), "%Y")
    if (as.integer(min_yr) > 1999L) {
      sprintf("starts %s (ETF inception; paper needs 1999)", min_yr)
    } else {
      "OK"
    }
  }, character(1))

  gaps <- brics_gaps[brics_gaps != "OK"]
  if (length(gaps) > 0) {
    message("\nCoverage gaps vs paper sample (G7+BRICS from 1999-02-14):")
    for (ctry in names(gaps)) message(sprintf("  %s: %s", ctry, gaps[ctry]))
    message("Replace affected CSVs with Refinitiv Datastream exports for full replication.")
  }

  invisible(list(gold = gold, pd = pd_list))
}

# ---------------------------------------------------------------------------
# fetch_asset_yahoo: Download adjusted close price for any Yahoo Finance ticker
# ---------------------------------------------------------------------------
# name       : string used as file stem and asset label (e.g. "btc", "gold")
# ticker     : Yahoo Finance ticker (e.g. "BTC-USD", "GLD", "CL=F")
# output_dir : directory; file saved as {output_dir}/{name}.csv
# start_date, end_date : character "YYYY-MM-DD"
#
# Saves: columns date, price, log_price
# Returns data.table invisibly, or NULL on failure.
fetch_asset_yahoo <- function(name, ticker,
                               output_dir = "data/assets",
                               start_date = "2015-01-01",
                               end_date   = format(Sys.Date())) {
  message(sprintf("Fetching %s (%s) from Yahoo Finance...", name, ticker))

  prices <- tryCatch(
    getSymbols(ticker, src = "yahoo",
               from = start_date, to = end_date,
               auto.assign = FALSE, warnings = FALSE),
    error = function(e) {
      warning(sprintf("  Could not download %s: %s", ticker, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(prices)) return(invisible(NULL))

  dt <- data.table(
    date  = as.Date(index(prices)),
    price = as.numeric(Ad(prices))
  )
  dt <- dt[is.finite(price) & price > 0]
  dt[, log_price := log(price)]

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(output_dir, paste0(name, ".csv"))
  fwrite(dt, out_path)
  message(sprintf("  Saved %d rows to %s", nrow(dt), out_path))

  invisible(dt)
}
