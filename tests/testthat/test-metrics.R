test_that("combine_reps collapses replicates correctly", {
  raw <- data.frame(
    cell_line     = rep("CL1", 6),
    dose          = rep(c(1, 10), each = 3),
    response      = c(10, 20, 30, 40, 50, 60),
    ancestry      = "AFR",
    feature       = "none",
    stringsAsFactors = FALSE
  )
  agg <- combine_reps(raw)
  expect_equal(nrow(agg), 2)
  expect_equal(agg$mean_response[agg$dose == 1],  20)
  expect_equal(agg$mean_response[agg$dose == 10], 50)
})

test_that("load_data clamps negative responses to zero", {
  tmp <- tempfile(fileext = ".csv")
  write.csv(data.frame(
    dose = c(1, 10), response = c(-5, 50),
    cell_line = "CL1", ancestry = "AFR", feature = "none"
  ), tmp, row.names = FALSE)
  dat <- load_data(tmp)
  expect_equal(dat$response[dat$dose == 1], 0)
  unlink(tmp)
})

test_that("validate_columns errors clearly on missing columns", {
  df <- data.frame(a = 1, b = 2)
  expect_error(
    validate_columns(df, c("a", "c"), arg_name = "test_df"),
    regexp = "test_df.*missing.*c"
  )
})

test_that("standardize_ancestry maps correctly", {
  expect_equal(standardize_ancestry("African"), "AFR")
  expect_equal(standardize_ancestry("European"), "EUR")
  expect_equal(standardize_ancestry("East Asian"), "EAS")
})