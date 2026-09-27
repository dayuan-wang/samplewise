#' samplewise: sample-level inference for derived single-cell features
#'
#' The package moves every derived single-cell feature through one table
#' shape (see [sw_long_table()]), fits one user-supplied design formula to
#' every feature series with the same eligibility rules ([fit_features()]),
#' and reports the result in one format ([write_results()], [plot_forest()]).
#' The biological sample, not the cell, is the unit of inference throughout.
#'
#' Read `docs/design.md` in the source repository for the design and the
#' roadmap.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(".data")
