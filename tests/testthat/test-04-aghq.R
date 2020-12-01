context("Quadrature")

test_that("Quadrature works",{
  # AGHQ
  expect_is(thequadrature,"aghq")
  expect_equal(names(thequadrature),c("normalized_posterior","marginals","optresults"))
  expect_is(summary(thequadrature),"aghqsummary")

  # Laplace approximation
  expect_is(thelaplace,"laplace")
  expect_equal(names(thelaplace),c("lognormconst","optresults"))
  expect_is(summary(thelaplace),"laplacesummary")

  # Marginal laplace approximation
  expect_is(themarginallaplace,"marginallaplace")
  expect_is(themarginallaplace,"aghq")
  expect_equal(names(themarginallaplace),c("normalized_posterior","marginals","optresults","modesandhessians"))
})
