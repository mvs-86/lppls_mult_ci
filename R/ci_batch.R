# R/ci_batch.R — Batch checkpointing for long LPPLS-CI computations
#
# Splits the t2 index range into chunks, computes each chunk independently,
# and saves results to RDS files.  Already-completed chunks are skipped on
# re-run, making the computation resumable after interruption.
#
# Typical workflow:
#
#   compute_ci_batched(btc_log_price, dates, n_batches = 10,
#                      cache_dir = "cache/btc_ci")
#
#   ci_dt <- merge_ci_batches("cache/btc_ci")
#
# Or check progress mid-run:
#   ci_batch_status("cache/btc_ci")

library(data.table)

# ---------------------------------------------------------------------------
# compute_ci_batched: Run compute_ci_asset in resumable batches
# ---------------------------------------------------------------------------
# x, dates, value_col, scales, n_starts, method, parallel_t2 :
#   passed through to compute_ci_asset — same semantics.
#
# n_batches  : split t2 indices into this many roughly equal chunks.
# cache_dir  : directory where per-batch RDS files are stored.
#              Created if it does not exist.
# overwrite  : if TRUE, recompute and overwrite existing batch files.
#              Default FALSE (skip completed batches).
#
# Returns the merged data.table (same shape as compute_ci_asset output),
# reading from cache when batches already exist.
compute_ci_batched <- function(x, dates = NULL, value_col = NULL,
                               scales    = SCALES,
                               n_starts  = 3L,
                               method    = "grid",
                               parallel_t2 = TRUE,
                               n_batches = 10L,
                               cache_dir = "cache/ci_batches",
                               overwrite = FALSE) {
  # --- Resolve lnp_vec + dates (same logic as compute_ci_asset) ---
  if (is.data.table(x) || is.data.frame(x)) {
    dt_in <- as.data.table(x)
    if (is.null(value_col)) {
      num_cols <- setdiff(names(dt_in)[vapply(dt_in, is.numeric, logical(1))], "date")
      if (length(num_cols) != 1L)
        stop("compute_ci_batched: cannot auto-detect value column — specify value_col")
      value_col <- num_cols[1L]
    }
    lnp_vec <- dt_in[[value_col]]
    dates   <- dt_in[["date"]]
  } else {
    lnp_vec <- x
  }
  if (is.null(dates) || length(dates) != length(lnp_vec))
    stop("compute_ci_batched: dates must be supplied and same length as x")

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  t2_all  <- which(is.finite(lnp_vec))
  n_t2    <- length(t2_all)
  if (n_t2 == 0L) stop("compute_ci_batched: no finite observations found")

  # Split into n_batches chunks
  batch_ids <- cut(seq_len(n_t2), breaks = n_batches, labels = FALSE)
  batches   <- split(t2_all, batch_ids)
  n_batches_actual <- length(batches)

  message(sprintf(
    "compute_ci_batched: %d t2 indices → %d batches of ~%d",
    n_t2, n_batches_actual, ceiling(n_t2 / n_batches_actual)
  ))

  # Process each batch
  for (b in seq_len(n_batches_actual)) {
    rds_path <- file.path(cache_dir, sprintf("batch_%04d.rds", b))

    if (!overwrite && file.exists(rds_path)) {
      message(sprintf("  Batch %d/%d — cached, skipping",
                      b, n_batches_actual))
      next
    }

    t2_batch <- batches[[b]]
    message(sprintf("  Batch %d/%d — %d t2 indices (t2 = %d..%d)",
                    b, n_batches_actual, length(t2_batch),
                    min(t2_batch), max(t2_batch)))

    batch_result <- .compute_ci_t2_subset(
      lnp_vec     = lnp_vec,
      dates       = dates,
      t2_indices  = t2_batch,
      scales      = scales,
      n_starts    = n_starts,
      method      = method,
      parallel_t2 = parallel_t2
    )

    saveRDS(batch_result, rds_path)
    message(sprintf("  Batch %d/%d — saved to %s", b, n_batches_actual, rds_path))
  }

  merge_ci_batches(cache_dir)
}

# ---------------------------------------------------------------------------
# merge_ci_batches: Combine all batch RDS files into one data.table
# ---------------------------------------------------------------------------
merge_ci_batches <- function(cache_dir) {
  rds_files <- sort(list.files(cache_dir, pattern = "^batch_[0-9]+\\.rds$",
                               full.names = TRUE))
  if (length(rds_files) == 0L)
    stop("merge_ci_batches: no batch files found in ", cache_dir)

  parts <- lapply(rds_files, readRDS)
  out   <- rbindlist(parts)
  setkey(out, date)
  out
}

# ---------------------------------------------------------------------------
# ci_batch_status: Report how many batches are done / pending / missing
# ---------------------------------------------------------------------------
ci_batch_status <- function(cache_dir) {
  rds_files <- list.files(cache_dir, pattern = "^batch_[0-9]+\\.rds$",
                           full.names = TRUE)
  if (length(rds_files) == 0L) {
    message("ci_batch_status: no batch files in ", cache_dir)
    return(invisible(NULL))
  }

  info <- data.table(
    file       = basename(rds_files),
    size_kb    = round(file.size(rds_files) / 1024, 1),
    modified   = file.mtime(rds_files)
  )
  setorder(info, file)

  message(sprintf("Cache dir : %s", cache_dir))
  message(sprintf("Batches   : %d complete", nrow(info)))
  message(sprintf("Total size: %.1f MB", sum(info$size_kb) / 1024))
  print(info)
  invisible(info)
}

# ---------------------------------------------------------------------------
# .compute_ci_t2_subset: Internal helper — run CI for a subset of t2 indices
# ---------------------------------------------------------------------------
# Same as compute_ci_asset but takes an explicit t2_indices vector rather than
# using all finite positions in lnp_vec.  Column layout is identical to
# compute_ci_asset output so batches can be rbindlist'd directly.
.compute_ci_t2_subset <- function(lnp_vec, dates, t2_indices,
                                  scales, n_starts, method, parallel_t2) {
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

  scale_abbrev <- function(nm) switch(nm, medium = "med", nm)

  out <- data.table(date = dates[t2_indices])
  for (sc_name in names(scales)) {
    ab <- scale_abbrev(sc_name)
    out[, paste0("pos_", ab) := ci_list[[sc_name]]$ci_pos]
    out[, paste0("neg_", ab) := ci_list[[sc_name]]$ci_neg]
  }
  setkey(out, date)
  out
}
