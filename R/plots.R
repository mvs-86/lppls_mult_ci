# R/plots.R — Visualisation for the lppls_mult_ci project
# Requires: ggplot2, patchwork, data.table

library(ggplot2)
library(patchwork)
library(data.table)

# ---------------------------------------------------------------------------
# .scale_label: human-readable panel title for a CI scale suffix
# ---------------------------------------------------------------------------
.scale_label <- function(suffix) {
  switch(suffix,
    short = "Short [30\u201390 days]",
    med   = "Medium [90\u2013300 days]",
    long  = "Long [300\u2013745 days]",
    tools::toTitleCase(suffix)   # fallback for custom scale names
  )
}

# ---------------------------------------------------------------------------
# .detect_scales: extract scale suffixes from pos_* columns of ci_dt
# ---------------------------------------------------------------------------
# Returns a named character vector: names = suffixes, values = display labels
.detect_scales <- function(ci_dt) {
  pos_cols <- grep("^pos_", names(ci_dt), value = TRUE)
  if (length(pos_cols) == 0L)
    stop("plot_ci_asset: ci_dt has no columns matching 'pos_*' — ",
         "run compute_ci_asset() first")
  suffixes <- sub("^pos_", "", pos_cols)
  setNames(vapply(suffixes, .scale_label, character(1L)), suffixes)
}

# ---------------------------------------------------------------------------
# plot_ci_asset: generic CI visualisation for any asset / scale set
# ---------------------------------------------------------------------------
# ci_dt      : data.table from compute_ci_asset() — must have date + pos_*/neg_* cols
# lnp_series : optional numeric vector of log prices / log P/D ratios to overlay
# dates_lnp  : Date vector for lnp_series; if NULL and lnp_series is given,
#              assumed aligned to ci_dt$date (same length)
# asset_name : string used as the overall plot title (e.g. "US", "BTC")
#
# Returns a patchwork object.
plot_ci_asset <- function(ci_dt,
                          lnp_series = NULL,
                          dates_lnp  = NULL,
                          asset_name = NULL) {

  if (!is.data.table(ci_dt))
    stop("plot_ci_asset: ci_dt must be a data.table")
  if (!"date" %in% names(ci_dt))
    stop("plot_ci_asset: ci_dt must have a 'date' column")

  # Resolve lnp dates
  if (!is.null(lnp_series)) {
    if (is.null(dates_lnp)) dates_lnp <- ci_dt$date
    if (length(lnp_series) != length(dates_lnp))
      stop("plot_ci_asset: lnp_series and dates_lnp must have the same length")
  }

  # ---- Pre-compute lnp scaling for secondary axis ----
  if (!is.null(lnp_series)) {
    lnp_dt    <- data.table(date = as.Date(dates_lnp), lnp = lnp_series)
    lnp_min   <- min(lnp_series, na.rm = TRUE)
    lnp_max   <- max(lnp_series, na.rm = TRUE)
    lnp_range <- lnp_max - lnp_min
    scale_lnp   <- function(x) -1 + 2 * (x - lnp_min) / lnp_range
    unscale_lnp <- function(x) lnp_min + (x + 1) / 2 * lnp_range
    y_label_sec <- if (!is.null(asset_name)) "Log P/D ratio" else "Log price"
  }

  scale_map <- .detect_scales(ci_dt)   # named label vector
  n_scales  <- length(scale_map)

  # ---- Build one ggplot per scale ----
  ci_panels <- lapply(seq_len(n_scales), function(i) {
    s         <- names(scale_map)[i]
    label     <- scale_map[[i]]
    is_bottom <- (i == n_scales)

    # Melt pos and neg to long; negate neg so it plots below zero
    wide <- ci_dt[, .(date,
                      pos = get(paste0("pos_", s)),
                      neg = get(paste0("neg_", s)))]
    long <- melt(wide, id.vars = "date",
                 measure.vars  = c("pos", "neg"),
                 variable.name = "direction",
                 value.name    = "raw")
    long[, value  := ifelse(direction == "neg", -raw, raw)]
    long[, signal := ifelse(direction == "pos", "Positive bubble", "Negative bubble")]
    long[, signal := factor(signal,
                            levels = c("Positive bubble", "Negative bubble"))]

    y_scale <- if (!is.null(lnp_series)) {
      scale_y_continuous(
        limits   = c(-1, 1),
        breaks   = c(-1, -0.5, 0, 0.5, 1),
        labels   = c("-1", "-0.5", "0", "0.5", "1"),
        sec.axis = sec_axis(~ unscale_lnp(.),
                            name = if (i == 1L) y_label_sec else NULL)
      )
    } else {
      scale_y_continuous(limits = c(-1, 1),
                         breaks = c(-1, -0.5, 0, 0.5, 1),
                         labels = c("-1", "-0.5", "0", "0.5", "1"))
    }

    p <- ggplot(long, aes(x = date, y = value, fill = signal)) +
      geom_area(position = "identity", alpha = 0.75) +
      geom_hline(yintercept = 0, linewidth = 0.35, colour = "black") +
      scale_fill_manual(
        values = c("Positive bubble" = "#d62728",
                   "Negative bubble" = "#1f77b4"),
        name   = NULL
      ) +
      y_scale +
      scale_x_date(date_labels = "%Y") +
      labs(title = label, x = NULL, y = "CI") +
      theme_bw(base_size = 10) +
      theme(
        legend.position   = if (i == 1L) "top" else "none",
        legend.key.size   = unit(0.4, "cm"),
        axis.text.x       = if (is_bottom) element_text() else element_blank(),
        axis.ticks.x      = if (is_bottom) element_line() else element_blank(),
        panel.grid.minor  = element_blank(),
        plot.title        = element_text(size = 9, face = "bold")
      )

    if (!is.null(lnp_series)) {
      p <- p + geom_line(
        data        = lnp_dt,
        aes(x = date, y = scale_lnp(lnp_series)),
        inherit.aes = FALSE,
        colour      = "black",
        linewidth   = 0.45
      )
    }

    p
  })

  combined <- Reduce(`/`, ci_panels) +
    plot_layout(heights = rep(1, n_scales))

  # ---- Overall title ----
  if (!is.null(asset_name)) {
    combined <- combined +
      plot_annotation(
        title = asset_name,
        theme = theme(plot.title = element_text(
          face = "bold", size = 12, hjust = 0.5
        ))
      )
  }

  combined
}

# ---------------------------------------------------------------------------
# plot_ci_country: paper-replication wrapper (matches main.R call signature)
# ---------------------------------------------------------------------------
# ci_dt        : data.table from compute_ci_all()[[country]]
# lnp_weekly   : numeric vector aligned to ci_dt$date (from pd_dt$lnp_weekly)
# country_name : e.g. "us", "south_africa"
plot_ci_country <- function(ci_dt, lnp_weekly, country_name) {
  clean_name <- toupper(gsub("_", " ", country_name))
  plot_ci_asset(
    ci_dt      = ci_dt,
    lnp_series = lnp_weekly,
    dates_lnp  = ci_dt$date,
    asset_name = clean_name
  )
}

# ---------------------------------------------------------------------------
# plot_gold_data: Figure 1 — gold log-returns + TVPGARCH volatility
# ---------------------------------------------------------------------------
# gold_dt : data.table from load_gold() — columns: date, price, log_return
# vol_dt  : data.table with columns date, vol (TVPGARCH σ²_t series)
plot_gold_data <- function(gold_dt, vol_dt) {
  ret_data <- gold_dt[is.finite(log_return)]

  p_ret <- ggplot(ret_data, aes(x = date, y = log_return)) +
    geom_col(fill = "grey60", colour = NA, width = 5) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    scale_x_date(date_labels = "%Y") +
    labs(title = "Gold log-return", x = NULL, y = "Log return") +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      panel.grid.minor = element_blank()
    )

  p_vol <- ggplot(vol_dt, aes(x = date, y = vol)) +
    geom_ribbon(aes(ymin = 0, ymax = vol), fill = "#aec7e8", alpha = 0.8) +
    geom_line(colour = "#1f77b4", linewidth = 0.45) +
    scale_x_date(date_labels = "%Y") +
    labs(
      title = expression("Conditional variance  " * sigma[t]^2 * "  (TVPGARCH)"),
      x     = "Date",
      y     = expression(sigma[t]^2)
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())

  (p_ret / p_vol) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      title = "Figure 1: Gold Returns and Conditional Variance",
      theme = theme(plot.title = element_text(
        face = "bold", size = 11, hjust = 0.5
      ))
    )
}

# ---------------------------------------------------------------------------
# plot_ci_price: price bars coloured by CI bubble signal
# ---------------------------------------------------------------------------
# Displays price as a bar chart where fill encodes the dominant bubble signal:
#   red  (light→dark) = positive bubble at increasing confidence
#   blue (light→dark) = negative bubble at increasing confidence
#   near-white        = no dominant signal
#
# signal = pos_{scale} − neg_{scale} ∈ [-1, 1] drives a diverging colour scale.
#
# ci_dt        : data.table from compute_ci_asset() — needs date + pos_*/neg_*
# price_series : numeric vector — prices or log prices (y-axis bar height)
# dates_price  : Date vector for price_series; NULL = assumed aligned to ci_dt$date
# scale        : CI scale suffix to colour by (e.g. "short", "med"); NULL = first detected
# asset_name   : plot title string
#
# Returns a ggplot object.
plot_ci_price <- function(ci_dt,
                           price_series,
                           dates_price = NULL,
                           scale       = NULL,
                           asset_name  = NULL) {

  if (!is.data.table(ci_dt) || !"date" %in% names(ci_dt))
    stop("plot_ci_price: ci_dt must be a data.table with a 'date' column")
  if (!is.numeric(price_series))
    stop("plot_ci_price: price_series must be a numeric vector")
  if (is.null(dates_price)) dates_price <- ci_dt$date
  if (length(price_series) != length(dates_price))
    stop("plot_ci_price: price_series and dates_price must have the same length")

  # Detect / validate scale
  scale_map <- .detect_scales(ci_dt)
  s <- if (is.null(scale)) names(scale_map)[1L] else scale
  if (!s %in% names(scale_map))
    stop("plot_ci_price: scale '", s, "' not found — available: ",
         paste(names(scale_map), collapse = ", "))

  # Inner-join price and CI by date
  price_dt <- data.table(date = as.Date(dates_price), price = as.numeric(price_series))
  dt <- merge(
    ci_dt[, .(date,
               pos = get(paste0("pos_", s)),
               neg = get(paste0("neg_", s)))],
    price_dt,
    by = "date"
  )
  dt[, signal := pos - neg]

  # Bar width = median inter-date gap (works for daily, weekly, or custom freq)
  bar_width <- if (nrow(dt) > 1L) as.numeric(median(diff(dt$date))) else 1

  p <- ggplot(dt, aes(x = date, y = price, fill = signal)) +
    geom_col(width = bar_width, colour = NA) +
    scale_fill_gradient2(
      low      = "#1f77b4",
      mid      = "grey92",
      high     = "#d62728",
      midpoint = 0,
      limits   = c(-1, 1),
      name     = paste0("CI signal\n(", .scale_label(s), ")")
    ) +
    scale_x_date(date_labels = "%Y") +
    labs(x = NULL, y = "Price") +
    theme_bw(base_size = 10) +
    theme(
      legend.position   = "right",
      legend.key.height = unit(1.5, "cm"),
      panel.grid.minor  = element_blank()
    )

  if (!is.null(asset_name))
    p <- p + labs(title = asset_name) +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))

  p
}
