context("Quadrature")

quadnames <- c("normalized_posterior","optresults","control","transformation","marginals")
quadnamesmarg <- c("normalized_posterior","optresults","modesandhessians","control","transformation","marginals")


test_that("Quadrature works",{
  # AGHQ
  expect_is(thequadrature,"aghq")
  expect_true(setequal(names(thequadrature),quadnames))
  expect_is(summary(thequadrature),"aghqsummary")
  expect_is(thequadrature_correct,"aghq")
  expect_true(setequal(names(thequadrature_correct),quadnames))
  expect_is(summary(thequadrature_correct),"aghqsummary")

  expect_is(thequadrature3d,"aghq")
  expect_true(setequal(names(thequadrature3d),quadnames))
  expect_is(summary(thequadrature3d),"aghqsummary")

  # Laplace approximation
  expect_is(thelaplace,"laplace")
  expect_true(setequal(names(thelaplace),c("lognormconst","optresults")))
  expect_is(summary(thelaplace),"laplacesummary")

  # Marginal laplace approximation
  expect_is(themarginallaplace,"marginallaplace")
  expect_is(themarginallaplace,"aghq")
  expect_true(setequal(names(themarginallaplace),quadnamesmarg))

  expect_is(themarginallaplace3d_1,"marginallaplace")
  expect_is(themarginallaplace3d_1,"aghq")
  expect_true(setequal(names(themarginallaplace3d_1),quadnamesmarg))

  expect_is(themarginallaplace3d_2,"marginallaplace")
  expect_is(themarginallaplace3d_2,"aghq")
  expect_true(setequal(names(themarginallaplace3d_2),quadnamesmarg))

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


  # Sparse grids
  expect_true(!any(is.na(sparsegrid_2d$normalized_posterior$nodesandweights$logpost_normalized)))
  expect_true(!any(is.nan(sparsegrid_2d$normalized_posterior$nodesandweights$logpost_normalized)))
  expect_true(all(is.numeric(sparsegrid_2d$normalized_posterior$nodesandweights$logpost_normalized)))

  expect_equal(sparsenormconst_2d,1)

  # Test control params have correct options
  expect_true(all(c("method","negate","ndConstruction") %in% names(cntrl_base)))
  expect_true(all(c("method","inner_method","negate","ndConstruction") %in% names(cntrl_marg)))
  expect_true(all(c("method","numhessian","negate","ndConstruction") %in% names(cntrl_tmb)))

  # Laplace approximation
  expect_equal(la5,ls5)
  expect_equal(la10,ls10)
  expect_equal(la100,ls100)

  # Custom grid
  expect_equal(aghq_customgrid_gg1$normalized_posterior$lognormconst,aghq_customgrid_auto1$normalized_posterior$lognormconst)
  expect_equal(aghq_customgrid_gg2$normalized_posterior$lognormconst,aghq_customgrid_auto2$normalized_posterior$lognormconst)
  expect_equal(aghq_customgrid_gg2s$normalized_posterior$lognormconst,aghq_customgrid_auto2s$normalized_posterior$lognormconst)

  # Non-Gaussian custom grid- return an error
  expect_error(aghq(funlist2d,5,c(0,0),basegrid = gg4))


  # Control validation
  expect_error(validate_control(default_control_tmb(),type = "foo"))
  expect_true(validate_control(default_control()))
  expect_true(validate_control(default_control(),type = "aghq"))
  expect_true(validate_control(default_control_marglaplace(),type = "marglaplace"))
  expect_true(validate_control(default_control_tmb(),type = "tmb"))

  expect_error(validate_control(default_control(),type = 'marglaplace'))
  # expect_error(validate_control(default_control(),type = 'tmb')) # Currently, these have the same arguments
  expect_error(validate_control(default_control_marglaplace(),type = 'aghq'))
  expect_error(validate_control(default_control_marglaplace(),type = 'tmb'))
  # expect_error(validate_control(default_control_tmb(),type = 'aghq')) # Currently, these have the same arguments
  expect_error(validate_control(default_control_tmb(),type = 'marglaplace'))

  expect_error(validate_control(badcontrol1_aghq))
  expect_error(validate_control(badcontrol2_aghq))
  expect_error(validate_control(badcontrol3_aghq))

  expect_error(validate_control(badcontrol1_marglaplace,type='marglaplace'))
  expect_error(validate_control(badcontrol2_marglaplace,type='marglaplace'))
  expect_error(validate_control(badcontrol3_marglaplace,type='marglaplace'))

  expect_error(validate_control(badcontrol1_tmb,type='tmb'))
  expect_error(validate_control(badcontrol2_tmb,type='tmb'))
  expect_error(validate_control(badcontrol3_tmb,type='tmb'))

  # Test returning only normconst works
  expect_equal(class(aghq_normconst1),"numeric")
  expect_false(inherits(aghq_normconst1,'aghq'))
  expect_equal(aghq_normconst1,get_log_normconst(aghq_controlworks1))

  expect_equal(class(marglaplace_normconst1),"numeric")
  expect_false(inherits(marglaplace_normconst1,'aghq'))
  expect_false(inherits(marglaplace_normconst1,'marginallaplace'))
  expect_equal(marglaplace_normconst1,get_log_normconst(themarginallaplace))

  # Warnings for custom grids
  expect_warning(aghq(funlist2d,3,c(0,0),basegrid = gg5))

  # Setting k with custom grid works
  expect_equal(nrow(aghq_customgrid_gg6$normalized_posterior$nodesandweights),5^2)

  # Not modifying the grid
  expect_equal(gg1$features$move,gg7$features$move)

  # Test summary of marglaplace
  expect_equal(names(mlsumm1),c("aghqsummary","randomeffectsummary","info"))
  expect_equal(names(mlsumm2),c("aghqsummary","randomeffectsummary","info"))
  expect_equal(names(mlsumm3),c("aghqsummary","randomeffectsummary","info"))

  expect_is(mlsumm1,"marginallaplacesummary")
  expect_is(mlsumm2,"marginallaplacesummary")
  expect_is(mlsumm3,"marginallaplacesummary")
  expect_is(mlsumm4,"aghqsummary")
  expect_equal(nrow(mlsumm1$randomeffectsummary),1)
  expect_equal(nrow(mlsumm2$randomeffectsummary),2)
  expect_equal(nrow(mlsumm3$randomeffectsummary),1)

  expect_equal(mlsumm1$randomeffectsummary$variable,"W1")
  expect_equal(mlsumm2$randomeffectsummary$variable,c("W1","W2"))
  expect_equal(mlsumm3$randomeffectsummary$variable,"W1")

  expect_equal(mlsumm1$info["M"],c("M" = 1000))
  expect_equal(mlsumm2$info["M"],c("M" = 1000))
  expect_equal(mlsumm3$info["M"],c("M" = 100))

  expect_output(summary(themarginallaplace3d_1,max_print=1))

  # Test naming in marginallaplace
  expect_equal(names(themarginallaplace$modesandhessians$mode[[1]]),"W1")
  expect_equal(names(themarginallaplace3d_1$modesandhessians$mode[[1]]),c("W1","W2"))

  # Summaries with corrections
  expect_equal(thesummary_reuse$mode,thesummary_correct$mode)
  expect_equal(thesummary_reuse$lognormconst,thesummary_correct$lognormconst)
  expect_equal(thesummary_reuse$hessian,thesummary_correct$hessian)
  expect_equal(thesummary_reuse$covariance,thesummary_correct$covariance)
  expect_equal(thesummary_reuse$quadpoints,thesummary_correct$quadpoints)
  expect_equal(thesummary_reuse$dim,thesummary_correct$dim)
  expect_lt(sum(abs(thesummary_reuse$summarytable$mean - thesummary_correct$summarytable$mean)),1e-02)
  expect_lt(sum(abs(thesummary_reuse$summarytable$mean - truemean_correct)),1e-02)
  expect_lt(sum(abs(thesummary_correct$summarytable$mean - truemean_correct)),1e-02)
  expect_lt(sum(abs(thesummary_reuse$summarytable$sd - thesummary_correct$summarytable$sd)),1e-02)
  expect_lt(sum(abs(thesummary_reuse$summarytable$sd - truesd_correct)),1e-02)
  expect_lt(sum(abs(thesummary_correct$summarytable$sd - truesd_correct)),1e-02)
  expect_equal(thesummary_reuse$summarytable$mean == thesummary_correct$summarytable$mean,c(FALSE,FALSE))
  expect_equal(thesummary_reuse$summarytable$sd == thesummary_correct$summarytable$sd,c(FALSE,FALSE))
  expect_setequal(thequadrature_reuse$marginals[[1]]$theta1,thequadrature_correct$marginals[[1]]$theta1)
  expect_equal(thequadrature_reuse$marginals[[1]]$logmargpost==thequadrature_correct$marginals[[1]]$logmargpost,c(FALSE,FALSE,FALSE))
  expect_lt(sum(abs(thequadrature_reuse$marginals[[2]]$theta2 - thequadrature_correct$marginals[[2]]$theta2)),1e-15)
  expect_equal(thequadrature_reuse$marginals[[2]]$logmargpost==thequadrature_correct$marginals[[2]]$logmargpost,c(FALSE,FALSE,FALSE))



})
