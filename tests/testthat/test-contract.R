test_that("sw_long_table builds and validates", {
  x <- sw_long_table(c("m1", "m2"), "composition", "B", "proportion_of_all",
                     c(0.4, 0.5), c(40L, 50L))
  expect_s3_class(x, "data.frame")
  expect_named(x, c("sample_id", "feature_type", "unit", "feature", "value",
                    "n_cells", "estimable"))
  expect_equal(x$estimable, c(TRUE, TRUE))
  expect_equal(nrow(x), 2)
})

test_that("sw_validate rejects broken tables", {
  x <- sw_long_table(c("m1", "m2"), "composition", "B", "proportion_of_all", c(0.4, 0.5))
  expect_error(sw_validate(x[, -1]), "missing column")
  y <- x
  y$estimable <- c(TRUE, NA)
  expect_error(sw_validate(y), "NA")
  z <- x
  z$value[1] <- NA
  expect_error(sw_validate(z), "estimable row must carry a value")
  w <- rbind(x, x)
  expect_error(sw_validate(w), "duplicate")
})

test_that("a non-estimable row may carry NA", {
  x <- sw_long_table(c("m1", "m2"), "communication", "B>T", "MHC-II",
                     c(NA, 0.2), c(3L, 40L), estimable = c(FALSE, TRUE))
  expect_true(is.na(x$value[1]))
  expect_false(x$estimable[1])
})
