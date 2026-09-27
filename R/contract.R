SW_COLUMNS <- c("sample_id", "feature_type", "unit", "feature", "value",
                "n_cells", "estimable")

#' The per-sample long table
#'
#' Every derived feature passes through this shape before a model is fitted.
#' One row is one measurement of one feature in one sample.
#'
#' * `sample_id` - the biological replicate (mouse, donor, patient). Must
#'   match the sample column of the metadata given to [fit_features()].
#' * `feature_type` - the family of measurements: `"composition"`,
#'   `"communication"`, `"drug"`, ... By default one Benjamini-Hochberg
#'   correction is applied per `feature_type` and model term.
#' * `unit` - the biological context the value is measured in: a cell type,
#'   a `sender>receiver` pair.
#' * `feature` - what was measured within the unit: `"proportion_of_all"`,
#'   a signalling pathway, a ligand-receptor pair, a drug.
#' * `value` - the per-sample quantification, on the scale the builder
#'   documents. `NA` when not estimable.
#' * `n_cells` - the number of cells the value rests on in that sample.
#' * `estimable` - `TRUE` when the builder considers the value defined in
#'   that sample (enough cells, both partners present, ...). Rows that are
#'   not estimable are excluded from fitting. They are never zeros: a zero
#'   is a measured value, `NA` is the absence of a measurement.
#'
#' @param sample_id,feature_type,unit,feature Character vectors, recycled to
#'   the length of `value`.
#' @param value Numeric vector.
#' @param n_cells Integer vector, recycled.
#' @param estimable Logical vector, recycled. Defaults to `!is.na(value)`.
#' @return A validated `data.frame` with the seven columns above.
#' @examples
#' sw_long_table(c("m1", "m2"), "composition", "B cell", "proportion_of_all",
#'               c(0.41, 0.39), c(4100L, 3900L))
#' @export
sw_long_table <- function(sample_id, feature_type, unit, feature, value,
                          n_cells = NA_integer_, estimable = NULL) {
  n <- length(value)
  if (is.null(estimable)) estimable <- !is.na(value)
  x <- data.frame(
    sample_id = rep_len(as.character(sample_id), n),
    feature_type = rep_len(as.character(feature_type), n),
    unit = rep_len(as.character(unit), n),
    feature = rep_len(as.character(feature), n),
    value = as.numeric(value),
    n_cells = rep_len(as.integer(n_cells), n),
    estimable = rep_len(as.logical(estimable), n),
    stringsAsFactors = FALSE
  )
  x <- sw_validate(x)
  x
}

#' Validate a per-sample long table
#'
#' @param x A `data.frame`.
#' @return `x` invisibly, with the key columns coerced to character, or an
#'   error naming the first problem found.
#' @seealso [sw_long_table()] for the column definitions.
#' @export
sw_validate <- function(x) {
  if (!is.data.frame(x)) {
    stop("a long table must be a data.frame", call. = FALSE)
  }
  missing <- setdiff(SW_COLUMNS, names(x))
  if (length(missing)) {
    stop("long table is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  for (col in c("sample_id", "feature_type", "unit", "feature")) {
    if (!is.character(x[[col]])) x[[col]] <- as.character(x[[col]])
  }
  if (!is.numeric(x$value)) stop("`value` must be numeric", call. = FALSE)
  if (!is.logical(x$estimable)) stop("`estimable` must be logical", call. = FALSE)
  if (anyNA(x$estimable)) stop("`estimable` must not contain NA", call. = FALSE)
  if (any(x$estimable & is.na(x$value))) {
    stop("an estimable row must carry a value", call. = FALSE)
  }
  key <- paste(x$sample_id, x$feature_type, x$unit, x$feature, sep = "\r")
  if (anyDuplicated(key)) {
    stop("duplicate (sample_id, feature_type, unit, feature) rows", call. = FALSE)
  }
  invisible(x)
}

# cell-level metadata from the objects users actually have
.sw_cell_metadata <- function(cells) {
  if (is.data.frame(cells)) return(cells)
  if (inherits(cells, "Seurat")) return(methods::slot(cells, "meta.data"))
  if (inherits(cells, "SingleCellExperiment") &&
      requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    return(as.data.frame(SummarizedExperiment::colData(cells)))
  }
  stop("`cells` must be a data.frame of cell metadata, a Seurat object, ",
       "or a SingleCellExperiment", call. = FALSE)
}
