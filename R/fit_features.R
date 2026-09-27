.sw_rhs <- function(formula) {
  if (is.character(formula)) formula <- stats::as.formula(formula)
  if (!inherits(formula, "formula")) stop("`formula` must be a formula", call. = FALSE)
  if (length(formula) != 2L) {
    stop("`formula` must be one-sided, e.g. `~ age * sepsis + sex`; the ",
         "response is always the feature value", call. = FALSE)
  }
  formula
}

.sw_family <- function(family) {
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  if (!inherits(family, "family")) stop("`family` must be a stats family object", call. = FALSE)
  family
}

#' Coefficient names a design formula produces
#'
#' The names that `contrasts` in [fit_features()] refer to. Character and
#' logical columns are treated as factors, with the level order of the
#' metadata table (use [factor()] with an explicit `levels` argument to set
#' the reference level).
#'
#' @param metadata Sample-level `data.frame`.
#' @param formula One-sided design formula.
#' @return Character vector of coefficient names.
#' @examples
#' meta <- data.frame(age = c("young", "old"), sepsis = c("naive", "sepsis"))
#' sw_coef_names(meta, ~ age * sepsis)
#' @export
sw_coef_names <- function(metadata, formula) {
  formula <- .sw_rhs(formula)
  vars <- all.vars(formula)
  absent <- setdiff(vars, names(metadata))
  if (length(absent)) {
    stop("formula variable(s) not in metadata: ", paste(absent, collapse = ", "),
         call. = FALSE)
  }
  for (v in vars) {
    if (is.character(metadata[[v]]) || is.logical(metadata[[v]])) {
      metadata[[v]] <- factor(metadata[[v]])
    }
  }
  colnames(stats::model.matrix(formula, data = metadata))
}

# expand the user's `contrasts` argument into a named list of full-length
# contrast vectors over the fitted coefficients
.sw_terms <- function(coefs, contrasts) {
  cn <- names(coefs)
  unit <- function(nm) {
    v <- stats::setNames(numeric(length(cn)), cn)
    v[nm] <- 1
    v
  }
  if (is.null(contrasts)) {
    nms <- setdiff(cn, "(Intercept)")
    out <- lapply(nms, unit)
    names(out) <- nms
    return(out)
  }
  if (is.character(contrasts)) {
    bad <- setdiff(contrasts, cn)
    if (length(bad)) {
      stop("unknown coefficient(s): ", paste(bad, collapse = ", "),
           ". Available: ", paste(cn, collapse = ", "), call. = FALSE)
    }
    out <- lapply(contrasts, unit)
    names(out) <- contrasts
    return(out)
  }
  if (is.list(contrasts)) {
    if (is.null(names(contrasts)) || any(!nzchar(names(contrasts)))) {
      stop("`contrasts` given as a list must be fully named", call. = FALSE)
    }
    out <- lapply(contrasts, function(v) {
      if (!is.numeric(v)) stop("each contrast must be a numeric vector", call. = FALSE)
      if (!is.null(names(v))) {
        bad <- setdiff(names(v), cn)
        if (length(bad)) {
          stop("unknown coefficient(s) in contrast: ", paste(bad, collapse = ", "),
               ". Available: ", paste(cn, collapse = ", "), call. = FALSE)
        }
        full <- stats::setNames(numeric(length(cn)), cn)
        full[names(v)] <- v
        return(full)
      }
      if (length(v) != length(cn)) {
        stop("an unnamed contrast vector needs one entry per coefficient (",
             length(cn), "): ", paste(cn, collapse = ", "), call. = FALSE)
      }
      stats::setNames(v, cn)
    })
    return(out)
  }
  stop("`contrasts` must be NULL, a character vector of coefficient names, ",
       "or a named list of numeric contrast vectors", call. = FALSE)
}

.sw_fit_one <- function(s, metadata, sample_col, model_formula, fac_vars,
                        levels_full, family, transform, min_per_group,
                        min_detected, contrasts) {
  id <- s[1L, c("feature_type", "unit", "feature"), drop = FALSE]
  rownames(id) <- NULL
  d <- s[s$estimable & !is.na(s$value), c("sample_id", "value", "n_cells"), drop = FALSE]
  d <- merge(d, metadata, by.x = "sample_id", by.y = sample_col, sort = FALSE)
  n_fit <- nrow(d)
  n_detected <- sum(d$value > 0, na.rm = TRUE)
  row <- function(term, estimate = NA_real_, se = NA_real_, statistic = NA_real_,
                  df = NA_real_, p = NA_real_, status) {
    data.frame(id, term = term, estimate = estimate, se = se,
               statistic = statistic, df = df, p = p, n_fit = n_fit,
               n_detected = n_detected, status = status,
               stringsAsFactors = FALSE)
  }
  for (v in fac_vars) {
    d[[v]] <- factor(as.character(d[[v]]), levels = levels_full[[v]])
    if (any(table(d[[v]]) < min_per_group)) {
      return(row(NA_character_, status = paste0("too_few_per_group:", v)))
    }
  }
  if (min_detected > 0 && n_detected < min_detected) {
    return(row(NA_character_, status = "too_few_detected"))
  }
  d$.sw_y <- suppressWarnings(transform(d$value))
  if (any(!is.finite(d$.sw_y))) {
    return(row(NA_character_, status = "transform_not_finite"))
  }
  warn <- character()
  fit <- withCallingHandlers(
    tryCatch(stats::glm(model_formula, data = d, family = family),
             error = function(e) NULL),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (is.null(fit)) return(row(NA_character_, status = "fit_error"))
  b <- stats::coef(fit)
  V <- stats::vcov(fit)
  V[is.na(V)] <- 0
  t_based <- !(family$family %in% c("poisson", "binomial"))
  df_res <- fit$df.residual
  terms <- .sw_terms(b, contrasts)
  status <- if (length(warn)) "ok_with_warning" else "ok"
  out <- lapply(names(terms), function(nm) {
    cvec <- terms[[nm]]
    used <- cvec != 0
    if (any(is.na(b[used]))) return(row(nm, df = df_res, status = "term_not_estimable"))
    bb <- b
    bb[is.na(bb)] <- 0
    est <- sum(cvec * bb)
    se <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
    stat <- est / se
    p <- if (t_based) 2 * stats::pt(abs(stat), df = df_res, lower.tail = FALSE)
         else 2 * stats::pnorm(abs(stat), lower.tail = FALSE)
    row(nm, estimate = est, se = se, statistic = stat,
        df = if (t_based) df_res else Inf, p = p, status = status)
  })
  do.call(rbind, out)
}

#' Fit one design to every feature series
#'
#' Every (feature_type, unit, feature) series in the long table is joined to
#' the sample metadata and fitted with the same generalised linear model.
#' The sample is the unit: one row per sample enters each fit, so
#' pseudoreplication across cells cannot occur here by construction.
#'
#' A series is fitted only when it passes the eligibility rules; otherwise
#' the result row carries `NA` and a `status` naming the rule that failed:
#'
#' * every level of every factor in the formula (levels are taken from
#'   `metadata`, not from the series) has at least `min_per_group`
#'   estimable samples - `"too_few_per_group:<variable>"`;
#' * the value is above zero in at least `min_detected` samples -
#'   `"too_few_detected"` (set `min_detected = 0` to disable);
#' * `transform(value)` is finite - `"transform_not_finite"`;
#' * the fit converges - `"fit_error"`.
#'
#' Benjamini-Hochberg adjustment is applied within the groups defined by
#' `fdr_within` (default: per feature family and model term), across the
#' series that were fitted.
#'
#' @param x A long table, see [sw_long_table()].
#' @param metadata Sample-level `data.frame`, one row per sample, holding
#'   the sample column and every variable in `formula`. Character and
#'   logical columns become factors; set the reference level with
#'   [factor()] beforehand.
#' @param formula One-sided design formula, e.g. `~ age * sepsis + sex`.
#' @param contrasts What to report. `NULL` reports every non-intercept
#'   coefficient. A character vector names coefficients (see
#'   [sw_coef_names()]). A named list of numeric vectors defines linear
#'   contrasts of the coefficients, either named by coefficient or of full
#'   length.
#' @param family A [stats::family] (default `gaussian()`), or its name.
#' @param transform Function applied to `value` before fitting, e.g.
#'   `sqrt`, `log1p`, `stats::qlogis`.
#' @param min_per_group,min_detected Eligibility thresholds, see above.
#' @param fdr_within Columns of the result that define one
#'   Benjamini-Hochberg family.
#' @param sample_col Name of the sample column in `metadata`.
#' @return A `data.frame` with one row per fitted series and reported term:
#'   `feature_type`, `unit`, `feature`, `term`, `estimate`, `se`,
#'   `statistic`, `df`, `p`, `fdr`, `n_fit`, `n_detected`, `status`.
#'   Series that were not fitted keep one row with `term = NA`.
#' @examples
#' meta <- expand.grid(age = c("young", "old"), sepsis = c("naive", "sepsis"),
#'                     rep = 1:3, stringsAsFactors = FALSE)
#' meta$sample_id <- sprintf("m%02d", seq_len(nrow(meta)))
#' set.seed(1)
#' x <- sw_long_table(meta$sample_id, "composition", "B", "proportion_of_all",
#'                    plogis(-0.5 + 1 * (meta$sepsis == "sepsis") + rnorm(12, sd = 0.2)),
#'                    2000L)
#' fit_features(x, meta, ~ age * sepsis, transform = stats::qlogis)
#' @export
fit_features <- function(x, metadata, formula, contrasts = NULL,
                         family = stats::gaussian(), transform = identity,
                         min_per_group = 3L, min_detected = 5L,
                         fdr_within = c("feature_type", "term"),
                         sample_col = "sample_id") {
  x <- sw_validate(x)
  if (!is.data.frame(metadata)) stop("`metadata` must be a data.frame", call. = FALSE)
  if (!sample_col %in% names(metadata)) {
    stop("`metadata` has no column `", sample_col, "`", call. = FALSE)
  }
  metadata[[sample_col]] <- as.character(metadata[[sample_col]])
  if (anyDuplicated(metadata[[sample_col]])) {
    stop("`metadata` must have one row per sample", call. = FALSE)
  }
  unknown <- setdiff(unique(x$sample_id), metadata[[sample_col]])
  if (length(unknown)) {
    stop("sample(s) in the long table but not in `metadata`: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  formula <- .sw_rhs(formula)
  vars <- all.vars(formula)
  absent <- setdiff(vars, names(metadata))
  if (length(absent)) {
    stop("formula variable(s) not in metadata: ", paste(absent, collapse = ", "),
         call. = FALSE)
  }
  if (".sw_y" %in% names(metadata)) stop("`.sw_y` is reserved", call. = FALSE)
  family <- .sw_family(family)
  if (!is.function(transform)) stop("`transform` must be a function", call. = FALSE)
  fac_vars <- vars[vapply(vars, function(v) {
    is.character(metadata[[v]]) || is.factor(metadata[[v]]) || is.logical(metadata[[v]])
  }, logical(1))]
  for (v in fac_vars) if (!is.factor(metadata[[v]])) metadata[[v]] <- factor(metadata[[v]])
  levels_full <- lapply(fac_vars, function(v) levels(metadata[[v]]))
  names(levels_full) <- fac_vars
  model_formula <- stats::update(formula, .sw_y ~ .)
  key <- paste(x$feature_type, x$unit, x$feature, sep = "\r")
  parts <- split(x, factor(key, levels = unique(key)))
  res <- lapply(parts, .sw_fit_one, metadata = metadata, sample_col = sample_col,
                model_formula = model_formula, fac_vars = fac_vars,
                levels_full = levels_full, family = family, transform = transform,
                min_per_group = min_per_group, min_detected = min_detected,
                contrasts = contrasts)
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out$fdr <- NA_real_
  ok <- !is.na(out$p)
  if (any(ok)) {
    grp <- interaction(out[ok, fdr_within, drop = FALSE], drop = TRUE)
    out$fdr[ok] <- stats::ave(out$p[ok], grp, FUN = function(p) stats::p.adjust(p, "BH"))
  }
  out <- out[order(out$feature_type, out$term, out$p, na.last = TRUE), ]
  rownames(out) <- NULL
  out[, c("feature_type", "unit", "feature", "term", "estimate", "se",
          "statistic", "df", "p", "fdr", "n_fit", "n_detected", "status")]
}
