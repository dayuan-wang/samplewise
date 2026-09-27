cells <- data.frame(
  sample_id = rep(c("m1", "m2"), each = 6),
  cell_type = c("B", "B", "T", "T", "NK", "Mono",
                "B", "T", "T", "T", "NK", "Mono"),
  stringsAsFactors = FALSE
)

test_that("proportions of all cells sum to one per sample", {
  x <- build_composition(cells, min_denominator = 1)
  expect_equal(unique(x$feature), "proportion_of_all")
  expect_equal(as.numeric(tapply(x$value, x$sample_id, sum)), c(1, 1))
  expect_equal(x$value[x$sample_id == "m1" & x$unit == "B"], 2 / 6)
  expect_equal(x$n_cells[x$sample_id == "m2" & x$unit == "B"], 1L)
  expect_equal(attr(x, "denominator")$n_denominator, c(6L, 6L))
})

test_that("a compartment denominator restricts features and rescales", {
  x <- build_composition(cells, denominator = c("B", "T", "NK"),
                         denominator_name = "lymphocytes", min_denominator = 1)
  expect_setequal(unique(x$unit), c("B", "T", "NK"))
  expect_equal(unique(x$feature), "proportion_of_lymphocytes")
  expect_equal(x$value[x$sample_id == "m2" & x$unit == "T"], 3 / 5)
  expect_error(build_composition(cells, denominator = c("B", "Plasma")), "not present")
})

test_that("a small denominator makes the proportion non-estimable, absent types are zeros", {
  x <- build_composition(cells, min_denominator = 10)
  expect_true(all(!x$estimable))
  expect_true(all(is.na(x$value)))
  cells2 <- rbind(cells, data.frame(sample_id = "m3", cell_type = rep("B", 30)))
  y <- build_composition(cells2, min_denominator = 10)
  expect_equal(y$value[y$sample_id == "m3" & y$unit == "T"], 0)
  expect_true(y$estimable[y$sample_id == "m3" & y$unit == "T"])
})

