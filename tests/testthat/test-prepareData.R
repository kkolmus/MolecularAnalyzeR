test_that("relevant output", {
  expect_s3_class(inputSLanalyzeR, "data.frame")
})

test_that("relevant output with respect to columns", {
  expect_equal(ncol(inputSLanalyzeR), 4)
})
