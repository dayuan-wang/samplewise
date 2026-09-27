#' Per-sample cell-type composition
#'
#' Counts cells per sample and cell type and turns them into proportions of
#' a chosen denominator, one row per sample and cell type, in the
#' [sw_long_table()] shape. The denominator is either every cell in the
#' sample or a named compartment (for example the lymphocyte types), so the
#' same call yields "proportion of all cells" and "proportion of
#' lymphocytes".
#'
#' A proportion is `estimable` when its denominator holds at least
#' `min_denominator` cells in that sample; otherwise `value` is `NA`. A cell
#' type absent from a sample is a measured zero, not a missing value.
#'
#' @param cells A cell-level `data.frame`, a Seurat object or a
#'   SingleCellExperiment.
#' @param sample_col,celltype_col Column names in the cell metadata.
#' @param denominator `NULL` for all cells, or a character vector of cell
#'   types that define the denominator. Only cell types inside the
#'   denominator become features.
#' @param denominator_name Label used in `feature`
#'   (`"proportion_of_<name>"`). Defaults to `"all"` or `"subset"`.
#' @param min_denominator Smallest denominator for a proportion to be
#'   estimable.
#' @return A long table with `feature_type = "composition"`, `unit` = the
#'   cell type and `feature = "proportion_of_<denominator_name>"`. The
#'   per-sample denominators are attached as `attr(, "denominator")`.
#' @examples
#' cells <- data.frame(
#'   sample_id = rep(c("m1", "m2"), each = 6),
#'   cell_type = c("B", "B", "T", "T", "NK", "Mono", "B", "T", "T", "T", "NK", "Mono")
#' )
#' build_composition(cells, min_denominator = 1)
#' build_composition(cells, denominator = c("B", "T", "NK"),
#'                   denominator_name = "lymphocytes", min_denominator = 1)
#' @export
build_composition <- function(cells, sample_col = "sample_id",
                              celltype_col = "cell_type", denominator = NULL,
                              denominator_name = NULL, min_denominator = 20L) {
  meta <- .sw_cell_metadata(cells)
  for (col in c(sample_col, celltype_col)) {
    if (!col %in% names(meta)) stop("column not found in cell metadata: ", col, call. = FALSE)
  }
  s <- as.character(meta[[sample_col]])
  ct <- as.character(meta[[celltype_col]])
  keep <- !is.na(s) & !is.na(ct)
  s <- s[keep]
  ct <- ct[keep]
  tab <- table(sample_id = s, cell_type = ct)
  types <- colnames(tab)
  if (is.null(denominator)) {
    den_types <- types
    den_label <- if (is.null(denominator_name)) "all" else denominator_name
  } else {
    absent <- setdiff(denominator, types)
    if (length(absent)) {
      stop("denominator cell type(s) not present in the data: ",
           paste(absent, collapse = ", "), call. = FALSE)
    }
    den_types <- denominator
    den_label <- if (is.null(denominator_name)) "subset" else denominator_name
  }
  den <- rowSums(tab[, den_types, drop = FALSE])
  counts <- as.data.frame(tab, stringsAsFactors = FALSE)
  counts <- counts[counts$cell_type %in% den_types, , drop = FALSE]
  d <- as.numeric(den[counts$sample_id])
  estimable <- d >= min_denominator
  value <- ifelse(estimable & d > 0, counts$Freq / d, NA_real_)
  out <- sw_long_table(
    sample_id = counts$sample_id,
    feature_type = "composition",
    unit = counts$cell_type,
    feature = paste0("proportion_of_", den_label),
    value = value,
    n_cells = as.integer(counts$Freq),
    estimable = estimable
  )
  attr(out, "denominator") <- data.frame(
    sample_id = names(den), n_denominator = as.integer(den),
    stringsAsFactors = FALSE
  )
  out
}
