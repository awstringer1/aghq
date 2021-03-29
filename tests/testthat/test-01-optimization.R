context("Optimization")

test_that("optimizations works", {
  ## 1d

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
  expect_true(opt_sparsetrust$convergence %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_sr1$convergence %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_trust$convergence)
  expect_equal(opt_bfgs$convergence,0)
  # Methods
  expect_error(optimize_theta(funlist,1.5,control = list(method = "foo")))

  ## 2d

  # Mode
  expect_equal(round(opt_sparsetrust_2d$mode,5),round(truemode2d,5))
  expect_equal(round(opt_sr1_2d$mode,5),round(truemode2d,5))
  expect_equal(round(opt_trust_2d$mode,5),round(truemode2d,5))
  expect_equal(round(opt_bfgs_2d$mode,5),round(truemode2d,5))
  # Hessian
  expect_true(all(eigen(opt_sparsetrust_2d$hessian)$values > 0))
  expect_true(all(eigen(opt_sr1_2d$hessian)$values > 0))
  expect_true(all(eigen(opt_trust_2d$hessian)$values > 0))
  expect_true(all(eigen(opt_bfgs_2d$hessian)$values > 0))
  # Convergence
  expect_true(opt_sparsetrust_2d$convergence %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_sr1_2d$convergence %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_trust_2d$convergence)
  expect_equal(opt_bfgs_2d$convergence,0)

  ## 3d

  # Mode
  expect_equal(round(opt_sparsetrust_3d$mode,5),round(truemode3d,5))
  expect_equal(round(opt_sr1_3d$mode,5),round(truemode3d,5))
  expect_equal(round(opt_trust_3d$mode,5),round(truemode3d,5))
  expect_equal(round(opt_bfgs_3d$mode,5),round(truemode3d,5))
  # Hessian
  expect_true(all(eigen(opt_sparsetrust_3d$hessian)$values > 0))
  expect_true(all(eigen(opt_sr1_3d$hessian)$values > 0))
  expect_true(all(eigen(opt_trust_3d$hessian)$values > 0))
  expect_true(all(eigen(opt_bfgs_3d$hessian)$values > 0))
  # Convergence
  expect_true(opt_sparsetrust_3d$convergence %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_sr1_3d$convergence %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_trust_3d$convergence)
  expect_equal(opt_bfgs_3d$convergence,0)

})
