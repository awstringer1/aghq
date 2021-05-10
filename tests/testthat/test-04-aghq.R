context("Quadrature")

test_that("Quadrature works",{
  # AGHQ
  expect_is(thequadrature,"aghq")
  expect_equal(names(thequadrature),c("normalized_posterior","marginals","optresults"))
  expect_is(summary(thequadrature),"aghqsummary")

  expect_is(thequadrature3d,"aghq")
  expect_equal(names(thequadrature3d),c("normalized_posterior","marginals","optresults"))
  expect_is(summary(thequadrature3d),"aghqsummary")



  # Laplace approximation
  expect_is(thelaplace,"laplace")
  expect_equal(names(thelaplace),c("lognormconst","optresults"))
  expect_is(summary(thelaplace),"laplacesummary")

  # Marginal laplace approximation
  expect_is(themarginallaplace,"marginallaplace")
  expect_is(themarginallaplace,"aghq")
  expect_equal(names(themarginallaplace),c("normalized_posterior","marginals","optresults","modesandhessians"))

  expect_is(themarginallaplace3d_1,"marginallaplace")
  expect_is(themarginallaplace3d_1,"aghq")
  expect_equal(names(themarginallaplace3d_1),c("normalized_posterior","marginals","optresults","modesandhessians"))

  expect_is(themarginallaplace3d_2,"marginallaplace")
  expect_is(themarginallaplace3d_2,"aghq")
  expect_equal(names(themarginallaplace3d_2),c("normalized_posterior","marginals","optresults","modesandhessians"))

  # Sampling from marginal Laplace approximation
  expect_is(themargsamps,"list")
  expect_length(themargsamps,3)
  expect_equal(names(themargsamps),c('samps','theta','thetasamples'))
  expect_equal(colnames(themargsamps$theta),'theta1')
  expect_equal(dim(themargsamps$theta),c(10,1))
  expect_equal(dim(themargsamps$samps),c(1,10))
  expect_equal(length(themargsamps$thetasamples),1)
  expect_equal(length(themargsamps$thetasamples[[1]]),10)

  expect_is(themargsamps3d_1,"list")
  expect_length(themargsamps3d_1,3)
  expect_equal(names(themargsamps3d_1),c('samps','theta','thetasamples'))
  expect_equal(colnames(themargsamps3d_1$theta),'theta1')
  expect_equal(dim(themargsamps3d_1$theta),c(10,1))
  expect_equal(dim(themargsamps3d_1$samps),c(2,10))
  expect_equal(length(themargsamps3d_1$thetasamples),1)
  expect_equal(length(themargsamps3d_1$thetasamples[[1]]),10)

  expect_is(themargsamps3d_2,"list")
  expect_length(themargsamps3d_2,3)
  expect_equal(names(themargsamps3d_2),c('samps','theta','thetasamples'))
  expect_equal(colnames(themargsamps3d_2$theta),c('theta1','theta2'))
  expect_equal(dim(themargsamps3d_2$theta),c(10,2))
  expect_equal(dim(themargsamps3d_2$samps),c(1,10))
  expect_equal(length(themargsamps3d_2$thetasamples),2)
  expect_equal(length(themargsamps3d_2$thetasamples[[1]]),10)
  expect_equal(length(themargsamps3d_2$thetasamples[[2]]),10)


  # Sparse grids!
  # UPDATE: this is not a supported featufre right now
  # expect_true(!any(is.na(sparsegrid_2d$normalized_posterior$nodesandweights$logpost_normalized)))
  # expect_true(!any(is.nan(sparsegrid_2d$normalized_posterior$nodesandweights$logpost_normalized)))
  # expect_true(all(is.numeric(sparsegrid_2d$normalized_posterior$nodesandweights$logpost_normalized)))
  #
  # expect_equal(sparsenormconst_2d,1)

})
