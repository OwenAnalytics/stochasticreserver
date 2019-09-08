test_that("get_incremental_avg_matrix returns the correct matrix", {
  B0 <- stochasticreserver::B0
  A0 <- stochasticreserver::A0
  expect_equal(get_incremental_avg_matrix(B0), A0)
})
