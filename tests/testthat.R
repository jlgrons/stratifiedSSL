library(testthat)
library(stratifiedSSL)

test_check("stratifiedSSL")

test_that("multiplication works", {
  expect_equal(AccuracyStdErrorEstimation, 4)
})