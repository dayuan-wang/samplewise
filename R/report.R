#' Write a results table
#'
#' A `.csv` path writes one flat table. An `.xlsx` path writes a workbook
#' with a `README` sheet first and one sheet per model term, rows sorted by
#' raw p, which is the layout collaborators in this project receive.
#'
#' @param res Output of [fit_features()].
#' @param path File path ending in `.csv` or `.xlsx`.
#' @param readme Optional character vector of extra lines for the README
#'   sheet (what was fitted, which design, which eligibility rules).
#' @return `path`, invisibly.
#' @export
write_results <- function(res, path, readme = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "xlsx") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("the openxlsx package is needed to write .xlsx", call. = FALSE)
    }
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "README")
    lines <- c(
      "samplewise results",
      paste("written", format(Sys.time(), "%Y-%m-%d %H:%M")),
      "one sheet per model term; rows sorted by raw p",
      "fdr = Benjamini-Hochberg within the family stated in the fit (default: feature_type x term)",
      "status other than ok / ok_with_warning means the series was not fitted; see ?fit_features",
      readme
    )
    openxlsx::writeData(wb, "README", data.frame(note = lines, stringsAsFactors = FALSE))
    terms <- unique(res$term[!is.na(res$term)])
    used <- character()
    for (tm in terms) {
      sheet <- substr(gsub("[][*/\\\\?:]", "_", tm), 1, 28)
      if (sheet %in% used) sheet <- paste0(sheet, "_", length(used))
      used <- c(used, sheet)
      openxlsx::addWorksheet(wb, sheet)
      d <- res[!is.na(res$term) & res$term == tm, , drop = FALSE]
      d <- d[order(d$p, na.last = TRUE), , drop = FALSE]
      openxlsx::writeData(wb, sheet, d)
    }
    skipped <- res[is.na(res$term), , drop = FALSE]
    if (nrow(skipped)) {
      openxlsx::addWorksheet(wb, "not_fitted")
      openxlsx::writeData(wb, "not_fitted", skipped)
    }
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  } else {
    utils::write.csv(res, path, row.names = FALSE)
  }
  invisible(path)
}

#' Forest plot of one model term
#'
#' @param res Output of [fit_features()].
#' @param term One value of `res$term`.
#' @param top_n How many series to show, chosen by smallest p.
#' @param level Confidence level of the intervals.
#' @param fdr_cutoff Series at or below this FDR are drawn filled.
#' @return A ggplot object.
#' @export
plot_forest <- function(res, term, top_n = 30L, level = 0.95, fdr_cutoff = 0.05) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("the ggplot2 package is needed for plot_forest()", call. = FALSE)
  }
  d <- res[!is.na(res$term) & res$term == term & !is.na(res$estimate), , drop = FALSE]
  if (!nrow(d)) stop("no fitted series for term `", term, "`", call. = FALSE)
  d <- d[order(d$p), , drop = FALSE]
  d <- d[seq_len(min(top_n, nrow(d))), , drop = FALSE]
  q <- ifelse(is.finite(d$df), stats::qt(1 - (1 - level) / 2, d$df),
              stats::qnorm(1 - (1 - level) / 2))
  d$lo <- d$estimate - q * d$se
  d$hi <- d$estimate + q * d$se
  d$label <- paste(d$unit, d$feature, sep = " | ")
  d$label <- factor(d$label, levels = rev(unique(d$label)))
  d$significant <- !is.na(d$fdr) & d$fdr <= fdr_cutoff
  ggplot2::ggplot(d, ggplot2::aes(x = .data$estimate, y = .data$label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_segment(ggplot2::aes(x = .data$lo, xend = .data$hi,
                                       y = .data$label, yend = .data$label)) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$significant), size = 2.5) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                                name = paste0("FDR <= ", fdr_cutoff)) +
    ggplot2::labs(x = paste0("estimate (", term, ") with ", round(level * 100), "% CI"),
                  y = NULL) +
    ggplot2::theme_minimal(base_size = 11)
}
