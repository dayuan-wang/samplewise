#' Pseudobulk differential expression under the same design (planned)
#'
#' Interface reserved for v1. Aggregates counts per sample and cell type
#' (Seurat `AggregateExpression` semantics), drops pseudobulks built from
#' fewer than `min_cells` cells, applies the same per-group eligibility rule
#' as [fit_features()], and fits `formula` with the chosen engine. The
#' result table has the columns of [fit_features()] with `feature` = gene
#' and `unit` = cell type, so one report covers expression and derived
#' features alike.
#'
#' @param counts Genes x pseudobulk count matrix.
#' @param coldata `data.frame` with one row per column of `counts`, holding
#'   `sample_id`, `unit` (cell type) and `n_cells`.
#' @param metadata,formula,contrasts,min_per_group As in [fit_features()].
#' @param engine `"DESeq2"` (default, as published), `"edgeR"`
#'   (quasi-likelihood F-test) or `"dreamlet"` (precision-weighted linear
#'   mixed model, for designs with repeated samples).
#' @param min_cells Pseudobulks with fewer cells are dropped.
#' @export
fit_expression <- function(counts, coldata, metadata, formula, contrasts = NULL,
                           engine = c("DESeq2", "edgeR", "dreamlet"),
                           min_cells = 10L, min_per_group = 3L) {
  engine <- match.arg(engine)
  .NotYetImplemented()
}

#' Per-sample cell-cell communication strength from CellChat (planned)
#'
#' Interface reserved for v1. Takes CellChat objects built one per sample
#' and returns a long table at one aggregation level. A value is estimable
#' in a sample when sender and receiver both have more than `min_cells`
#' cells there; a zero in an estimable sample is CellChat's own result and
#' stays zero.
#'
#' @param objects Named list of CellChat objects, one per sample.
#' @param level One of `"pathway_edge"` (pathway x sender x receiver),
#'   `"lr_circuit"` (ligand-receptor pair x sender x receiver), `"edge"`
#'   (sender x receiver), `"pathway_total"`, `"lr_pair_total"`,
#'   `"sender_receiver"`.
#' @param min_cells Cells a sender or receiver needs in a sample for the
#'   value to be estimable.
#' @param senders,receivers Optional cell-type subsets.
#' @export
build_cellchat <- function(objects, level = c("pathway_edge", "lr_circuit", "edge",
                                              "pathway_total", "lr_pair_total",
                                              "sender_receiver"),
                           min_cells = 10L, senders = NULL, receivers = NULL) {
  level <- match.arg(level)
  .NotYetImplemented()
}
