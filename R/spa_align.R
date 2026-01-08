#' Stratigraphic Plug Alignment (SPA)
#'
#' Linearly interpolates plug-based measurements (e.g., TOC, porosity, XRD)
#' onto a high-resolution reference depth grid (e.g., XRF).
#' The procedure uses base R's \code{approx()} with \code{rule = 1} (default)
#' to prevent extrapolation beyond the observed depth range, ensuring
#' stratigraphically consistent alignment of all datasets.
#'
#' SPA is intended for vertically ordered core or log data, where measurements
#' are indexed by depth along a stratigraphic profile.
#'
#' @param ref A data.frame containing the reference depth grid and (optionally)
#'   high-resolution variables (e.g., XRF). Must contain the depth column
#'   specified in \code{depth_col}.
#' @param ... One or more named data.frames containing plug-based measurements
#'   to be interpolated (e.g., \code{xrd = xrd_df}, \code{plugs = plug_df}).
#' @param depth_col A character string giving the name of the depth column
#'   shared by all input datasets. Defaults to \code{"Depth_m"}.
#'   The depth column may use any unit (e.g., meters, feet, centimeters);
#'   \code{"Depth_m"} is only a column label and does not require depths to be
#'   in meters. However, all input datasets must use the same depth unit for
#'   interpolation to be meaningful.
#' @param rule Integer passed to \code{approx()} (default \code{1}).
#'   \code{rule = 1} prevents extrapolation outside the observed depth range.
#'   \code{rule = 2} performs constant extrapolation (end-point hold).
#' @param add_suffix Logical; if \code{TRUE}, variable names are suffixed with
#'   the dataset name (e.g., \code{TOC_plugs}, \code{Quartz_xrd}).
#' @param trim Logical; if \code{TRUE}, drops rows where any interpolated column
#'   is \code{NA}. Default \code{FALSE} (preserves all reference depths).
#'
#' @return A data.frame containing the reference depth grid and interpolated
#'   variables aligned to the same resolution.
#'
#' @examples
#' ref <- data.frame(Depth_m = 0:10, Ca = runif(11, 100, 200))
#' xrd <- data.frame(Depth_m = c(2, 5, 7), Quartz = c(54, 60, 58))
#' plugs <- data.frame(Depth_m = c(3, 7, 9), TOC = c(3.0, 3.3, 3.5))
#'
#' # Default: preserves reference depth grid
#' aligned <- spa_align(ref, xrd = xrd, plugs = plugs)
#' head(aligned)
#'
#' # Overlap-only alignment
#' aligned_overlap <- spa_align(ref, xrd = xrd, plugs = plugs, trim = TRUE)
#'
#' @export
spa_align <- function(ref,
                      ...,
                      depth_col = "Depth_m",
                      rule = 1,
                      add_suffix = TRUE,
                      trim = FALSE) {

  plug_list <- list(...)
  if (length(plug_list) == 0L) {
    stop("Provide datasets to interpolate, e.g., xrd = xrd_df")
  }
  if (!is.data.frame(ref)) stop("'ref' must be a data.frame.")
  if (!depth_col %in% names(ref)) stop("Depth column not in ref dataset: ", depth_col)

  ref_depth <- ref[[depth_col]]
  if (!is.numeric(ref_depth)) stop("Reference depth column must be numeric: ", depth_col)

  # Sort reference by depth (does not remove rows)
  ref <- ref[order(ref_depth), , drop = FALSE]
  ref_depth <- ref[[depth_col]]

  out <- ref
  added_cols <- character(0)

  for (nm in names(plug_list)) {
    df <- plug_list[[nm]]
    if (!is.data.frame(df)) stop("Dataset '", nm, "' is not a data.frame.")
    if (!depth_col %in% names(df)) stop("Depth column missing in dataset '", nm, "'.")

    ddepth <- df[[depth_col]]
    if (!is.numeric(ddepth)) stop("Depth column must be numeric in dataset '", nm, "'.")

    # Sort plug dataset by depth
    df <- df[order(ddepth), , drop = FALSE]
    ddepth <- df[[depth_col]]

    # Numeric columns excluding depth
    numeric_cols <- vapply(df, is.numeric, logical(1))
    numeric_cols[depth_col] <- FALSE
    vars <- names(df)[numeric_cols]

    depth_ok <- is.finite(ddepth)

    for (v in vars) {
      y <- df[[v]]
      ok <- depth_ok & is.finite(y)

      new_name <- if (add_suffix) paste0(v, "_", nm) else v

      # Default NA output (preserves reference grid)
      out[[new_name]] <- NA_real_

      # Need at least 2 valid points for linear interpolation
      if (sum(ok) >= 2L) {
        out[[new_name]] <- stats::approx(
          x    = ddepth[ok],
          y    = y[ok],
          xout = ref_depth,
          rule = rule,
          ties = "ordered"
        )$y
      }

      added_cols <- c(added_cols, new_name)
    }
  }

  # Optional: overlap-only trimming
  if (isTRUE(trim) && length(added_cols) > 0L) {
    keep <- stats::complete.cases(out[, added_cols, drop = FALSE])
    out <- out[keep, , drop = FALSE]
  }

  out
}
