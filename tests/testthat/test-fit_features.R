meta <- expand.grid(age = c("young", "old"), sex = c("male", "female"),
                    sepsis = c("naive", "sepsis"), rep = 1:2,
                    stringsAsFactors = FALSE)
meta$sample_id <- sprintf("m%02d", seq_len(nrow(meta)))
set.seed(1)
valA <- 1 + 2 * (meta$sepsis == "sepsis") + rnorm(16, sd = 0.3)
valB <- 5 + rnorm(16, sd = 0.3)
valC <- c(rep(0, 13), 1, 2, 3)
x <- rbind(
  sw_long_table(meta$sample_id, "composition", "A", "proportion_of_all", valA, 100L),
  sw_long_table(meta$sample_id, "composition", "B", "proportion_of_all", valB, 100L),
  sw_long_table(meta$sample_id, "composition", "C", "proportion_of_all", valC, 100L)
)

test_that("sw_coef_names matches glm coefficient names", {
  expect_equal(sw_coef_names(meta, ~ age * sepsis + sex),
               c("(Intercept)", "ageyoung", "sepsissepsis", "sexmale",
                 "ageyoung:sepsissepsis"))
  expect_error(sw_coef_names(meta, ~ age + weight), "not in metadata")
  expect_error(sw_coef_names(meta, y ~ age), "one-sided")
})

test_that("the sepsis effect is recovered and BH is applied per term", {
  res <- fit_features(x, meta, ~ age * sepsis + sex)
  a <- res[res$unit == "A" & res$term == "sepsissepsis", ]
  expect_equal(a$status, "ok")
  expect_equal(a$estimate, 2, tolerance = 0.3)
  expect_lt(a$p, 1e-4)
  expect_equal(a$n_fit, 16L)
  b <- res[res$unit == "B" & res$term == "sepsissepsis", ]
  expect_lt(abs(b$estimate), 3 * b$se)
  expect_true(all(!is.na(res$fdr[res$status == "ok"])))
  fam <- res[res$term == "sepsissepsis" & !is.na(res$p), ]
  expect_equal(fam$fdr, p.adjust(fam$p, "BH"))
})

test_that("eligibility rules produce status rows instead of fits", {
  res <- fit_features(x, meta, ~ age * sepsis + sex)
  cc <- res[res$unit == "C", ]
  expect_equal(nrow(cc), 1)
  expect_true(is.na(cc$term))
  expect_equal(cc$status, "too_few_detected")
  expect_equal(cc$n_detected, 3L)
  res0 <- fit_features(x, meta, ~ age * sepsis + sex, min_detected = 0)
  expect_true(all(res0$status[res0$unit == "C"] %in% c("ok", "ok_with_warning")))
  x_young <- x[x$sample_id %in% meta$sample_id[meta$age == "young"], ]
  res_y <- fit_features(x_young, meta, ~ age * sepsis + sex)
  expect_true(all(res_y$status == "too_few_per_group:age"))
  expect_true(all(is.na(res_y$p)))
})

test_that("contrasts by name and by vector agree with the default fit", {
  res <- fit_features(x, meta, ~ age * sepsis + sex)
  by_name <- fit_features(x, meta, ~ age * sepsis + sex,
                          contrasts = "ageyoung:sepsissepsis")
  by_vec <- fit_features(x, meta, ~ age * sepsis + sex,
                         contrasts = list(interaction = c(`ageyoung:sepsissepsis` = 1)))
  ref <- res[res$term == "ageyoung:sepsissepsis" & res$unit == "A", ]
  expect_equal(by_name$estimate[by_name$unit == "A"], ref$estimate)
  expect_equal(by_vec$term[by_vec$unit == "A"], "interaction")
  expect_equal(by_vec$estimate[by_vec$unit == "A"], ref$estimate)
  expect_equal(by_vec$p[by_vec$unit == "A"], ref$p)
  # sepsis effect in the young = main effect + interaction
  young_sepsis <- fit_features(x, meta, ~ age * sepsis + sex,
                               contrasts = list(sepsis_in_young = c(sepsissepsis = 1, `ageyoung:sepsissepsis` = 1)))
  ys <- young_sepsis[young_sepsis$unit == "A", ]
  expect_equal(ys$estimate, 2, tolerance = 0.5)
  expect_error(fit_features(x, meta, ~ age, contrasts = "nope"), "unknown coefficient")
})

test_that("families, transforms and input checks behave", {
  res <- fit_features(x, meta, ~ sepsis, family = "gaussian", transform = sqrt)
  expect_true(all(res$status[res$unit == "A"] == "ok"))
  expect_true(is.finite(res$df[1]))
  expect_error(fit_features(x, meta[-1, ], ~ sepsis), "not in `metadata`")
  expect_error(fit_features(x, meta, ~ weight), "not in metadata")
  expect_error(fit_features(x, meta, ~ sepsis, family = "nonsense"))
  xx <- x
  xx$value[xx$unit == "A"][1] <- -1
  r <- fit_features(xx, meta, ~ sepsis, transform = log)
  expect_equal(r$status[r$unit == "A"], "transform_not_finite")
})

test_that("write_results and plot_forest run", {
  res <- fit_features(x, meta, ~ age * sepsis + sex)
  f <- tempfile(fileext = ".csv")
  write_results(res, f)
  expect_true(file.exists(f))
  skip_if_not_installed("openxlsx")
  g <- tempfile(fileext = ".xlsx")
  write_results(res, g, readme = "test")
  expect_true(file.exists(g))
  skip_if_not_installed("ggplot2")
  p <- plot_forest(res, "sepsissepsis")
  expect_s3_class(p, "ggplot")
})
