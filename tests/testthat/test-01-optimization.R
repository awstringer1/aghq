context("Optimization")

test_that("optimizations works", {
  # Mode
  expect_equal(round(opt_sparsetrust$mode,5),round(truemode,5))
  expect_equal(round(opt_sr1$mode,5),round(truemode,5))
  expect_equal(round(opt_trust$mode,5),round(truemode,5))
  expect_equal(round(opt_bfgs$mode,5),round(truemode,5))
  # Hessian
  expect_true(all(eigen(opt_sparsetrust$hessian)$values > 0))
  expect_true(all(eigen(opt_sr1$hessian)$values > 0))
  expect_true(all(eigen(opt_trust$hessian)$values > 0))
  expect_true(all(eigen(opt_bfgs$hessian)$values > 0))
  # Convergence
  expect_equal(opt_sparsetrust$convergence,"Success")
  expect_equal(opt_sr1$convergence,"Radius of trust region is less than stop.trust.radius")
  expect_true(opt_trust$convergence)
  expect_equal(opt_bfgs$convergence,0)
  # Methods
  expect_error(optimize_theta(funlist,1.5,control = list(method = "foo")))
})
