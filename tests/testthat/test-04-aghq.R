context("Quadrature")

test_that("Quadrature works",{
  expect_is(thequadrature,"aghq")
  expect_equal(names(thequadrature),c("normalized_posterior","marginals","optresults"))
  expect_is(summary(thequadrature),"aghqsummary")
})