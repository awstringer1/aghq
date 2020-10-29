context("Summary Statistics")

test_that("Marginal posteriors computed correctly",{
  # Formatting of output object
  expect_is(margpost_1d_1,"tbl_df")
  expect_is(margpost_2d_1,"tbl_df")
  expect_is(margpost_2d_2,"tbl_df")

  expect_equal(colnames(margpost_1d_1),c("theta1","logmargpost","w"))
  expect_equal(colnames(margpost_2d_1),c("theta1","logmargpost","w"))
  expect_equal(colnames(margpost_2d_2),c("theta2","logmargpost","w"))

  expect_equal(nrow(margpost_1d_1),3)
  expect_equal(nrow(margpost_2d_1),3)
  expect_equal(nrow(margpost_2d_2),3)
  expect_equal(nrow(margpost_2d_2_k7),7)

  # Integration
  expect_equal(sum(exp(margpost_1d_1$logmargpost) * margpost_1d_1$w),1)
  expect_equal(sum(exp(margpost_2d_1$logmargpost) * margpost_2d_1$w),1)
  expect_equal(sum(exp(margpost_2d_2$logmargpost) * margpost_2d_2$w),1)
  expect_equal(sum(exp(margpost_2d_2_k7$logmargpost) * margpost_2d_2_k7$w),1)

  # Moments
  expect_equal(aghqnormconst1d,c('theta1' = 1))
  expect_equal(aghqnormconst2d,c('theta1' = 1,'theta2' = 1))

  expect_lt(abs(aghqmean1d - truemean1d),.02)
  expect_lt(abs(aghqmean2d - truemean2d)[1],.03)
  expect_lt(abs(aghqmean2d - truemean2d)[2],.02)

  expect_lt(abs(aghqexpmean1d - trueexpmean1d),1e-05)
  expect_lt(abs(aghqexpmean2d - trueexpmean2d)[1],1e-04)
  expect_lt(abs(aghqexpmean2d - trueexpmean2d)[2],1e-04)

  expect_lt(abs(aghqexpsd1d - truesd1d),1e-04)
  expect_lt(abs(aghqexpsd2d_1 - truesd2d[1]),1e-03)
  expect_lt(abs(aghqexpsd2d_2 - truesd2d[2]),1e-04)

  # Interpolation
  inherits(margpostinterp,"function")
  expect_equal(margpost_2d_1$logmargpost,margpostinterp(margpost_2d_1$theta1))

  # PDF and CDF
  expect_is(thepdfandcdf,"tbl_df")
  expect_equal(colnames(thepdfandcdf),c("theta","pdf","cdf"))
  expect_true(all(thepdfandcdf$pdf > 0))
  expect_true(all(thepdfandcdf$cdf >= 0))
  expect_equal(thepdfandcdf$cdf[1],0)
  expect_equal(round(thepdfandcdf$cdf[1000],2),1)

  # Quantiles
  expect_equal(names(thequantiles),c("2.5%","97.5%"))
  expect_false(any(is.infinite(thequantiles)))
  expect_lt(thequantiles[1],thequantiles[2])
  expect_lt(abs(exp(thequantiles[1]) - qgamma(.025,1+sum(y1),1+n1)),.07)
  expect_lt(abs(exp(thequantiles[2]) - qgamma(.975,1+sum(y1),1+n1)),.09)

})

