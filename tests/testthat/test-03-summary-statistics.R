context("Summary Statistics")

test_that("Marginal posteriors computed correctly",{
  # Formatting of output object
  expect_is(margpost_1d_1,"data.frame")
  expect_is(margpost_2d_1,"data.frame")
  expect_is(margpost_2d_2,"data.frame")
  expect_is(margpost_3d_1,"data.frame")
  expect_is(margpost_3d_2,"data.frame")
  expect_is(margpost_3d_3,"data.frame")


  expect_equal(colnames(margpost_1d_1),c("theta1","logmargpost","w"))
  expect_equal(colnames(margpost_2d_1),c("theta1","logmargpost","w"))
  expect_equal(colnames(margpost_2d_2),c("theta2","logmargpost","w"))
  expect_equal(colnames(margpost_3d_1),c("theta1","logmargpost","w"))
  expect_equal(colnames(margpost_3d_2),c("theta2","logmargpost","w"))
  expect_equal(colnames(margpost_3d_3),c("theta3","logmargpost","w"))


  expect_equal(nrow(margpost_1d_1),3)
  expect_equal(nrow(margpost_2d_1),3)
  expect_equal(nrow(margpost_2d_2),3)
  expect_equal(nrow(margpost_2d_2_k7),7)
  expect_equal(nrow(margpost_3d_1),3)
  expect_equal(nrow(margpost_3d_2),3)
  expect_equal(nrow(margpost_3d_2),3)


  # Integration
  expect_equal(sum(exp(margpost_1d_1$logmargpost) * margpost_1d_1$w),1)
  expect_equal(sum(exp(margpost_2d_1$logmargpost) * margpost_2d_1$w),1)
  expect_equal(sum(exp(margpost_2d_2$logmargpost) * margpost_2d_2$w),1)
  expect_equal(sum(exp(margpost_2d_2_k7$logmargpost) * margpost_2d_2_k7$w),1)
  expect_equal(sum(exp(margpost_3d_1$logmargpost) * margpost_3d_1$w),1)
  expect_equal(sum(exp(margpost_3d_2$logmargpost) * margpost_3d_2$w),1)
  expect_equal(sum(exp(margpost_3d_3$logmargpost) * margpost_3d_3$w),1)

  # Moments
  expect_equal(aghqnormconst1d,1)
  expect_equal(aghqnormconst2d,1)
  expect_equal(aghqnormconst3d,1)


  expect_lt(abs(aghqmean1d - truemean1d),.02)
  expect_lt(abs(aghqmean2d - truemean2d)[1],.03)
  expect_lt(abs(aghqmean2d - truemean2d)[2],.02)
  expect_lt(abs(aghqmean3d - truemean3d)[1],.1)
  expect_lt(abs(aghqmean3d - truemean3d)[2],.1)
  expect_lt(abs(aghqmean3d - truemean3d)[3],.1)


  expect_lt(abs(aghqexpmean1d - trueexpmean1d),1e-05)
  expect_lt(abs(aghqexpmean2d - trueexpmean2d)[1],1e-04)
  expect_lt(abs(aghqexpmean2d - trueexpmean2d)[2],1e-04)

  expect_lt(abs(aghqexpmean3d - trueexpmean3d)[1],1e-04)
  expect_lt(abs(aghqexpmean3d - trueexpmean3d)[2],1e-04)
  expect_lt(abs(aghqexpmean3d - trueexpmean3d)[3],1e-04)

  expect_lt(abs(aghqexpsd1d - truesd1d),1e-04)
  expect_lt(abs(aghqexpsd2d_1 - truesd2d[1]),1e-03)
  expect_lt(abs(aghqexpsd2d_2 - truesd2d[2]),1e-04)
  expect_lt(abs(aghqexpsd3d_1 - truesd3d[1]),1e-03)
  expect_lt(abs(aghqexpsd3d_2 - truesd3d[2]),1e-03)
  expect_lt(abs(aghqexpsd3d_3 - truesd3d[3]),1e-03)

  expect_equal(aghqnormconst2d_new,1)
  expect_lt(abs(aghqmean2d - truemean2d)[1],.1)
  expect_lt(abs(aghqmean2d - truemean2d)[2],.1)

  # Interpolation
  expect_is(margpostinterp,"function")
  expect_is(margpostinterp_2,"function")
  expect_is(margpostinterp_3d_1,"function")
  expect_is(margpostinterp_3d_2,"function")
  expect_is(margpostinterp_3d_3,"function")

  expect_equal(margpost_2d_1$logmargpost,margpostinterp(margpost_2d_1$theta1))
  expect_equal(margpost_2d_2$logmargpost,margpostinterp_2(margpost_2d_2$theta2))
  expect_equal(margpost_3d_1$logmargpost,margpostinterp_3d_1(margpost_3d_1$theta1))
  expect_equal(margpost_3d_2$logmargpost,margpostinterp_3d_2(margpost_3d_2$theta2))
  expect_equal(margpost_3d_3$logmargpost,margpostinterp_3d_3(margpost_3d_3$theta3))

  # PDF and CDF
  expect_is(thepdfandcdf,"data.frame")
  expect_equal(colnames(thepdfandcdf),c("theta","pdf","cdf"))
  expect_true(all(thepdfandcdf$pdf > 0))
  expect_true(all(thepdfandcdf$cdf >= 0))
  expect_equal(thepdfandcdf$cdf[1],0)
  expect_equal(round(thepdfandcdf$cdf[1000],2),1)


  # PDF and CDF with transformation
  expect_is(pdfwithtrans,"data.frame")
  expect_equal(colnames(pdfwithtrans),c("theta","pdf","cdf","transparam","pdf_transparam"))
  expect_true(all(pdfwithtrans$pdf > 0))
  expect_true(all(pdfwithtrans$pdf_transparam > 0))
  expect_true(all(pdfwithtrans$cdf >= 0))
  expect_equal(pdfwithtrans$cdf[1],0)
  expect_equal(round(pdfwithtrans$cdf[1000],2),1)

  # 3d
  expect_is(thepdfandcdf3d_1,"data.frame")
  expect_equal(colnames(thepdfandcdf3d_1),c("theta","pdf","cdf","transparam","pdf_transparam"))
  expect_true(all(thepdfandcdf3d_1$pdf > 0))
  expect_true(all(thepdfandcdf3d_1$pdf_transparam > 0))
  expect_true(all(thepdfandcdf3d_1$cdf >= 0))
  expect_equal(thepdfandcdf3d_1$cdf[1],0)
  expect_equal(round(thepdfandcdf3d_1$cdf[1000],2),1)

  expect_is(thepdfandcdf3d_2,"data.frame")
  expect_equal(colnames(thepdfandcdf3d_2),c("theta","pdf","cdf","transparam","pdf_transparam"))
  expect_true(all(thepdfandcdf3d_2$pdf > 0))
  expect_true(all(thepdfandcdf3d_2$pdf_transparam > 0))
  expect_true(all(thepdfandcdf3d_2$cdf >= 0))
  expect_equal(thepdfandcdf3d_2$cdf[1],0)
  expect_equal(round(thepdfandcdf3d_2$cdf[1000],2),1)

  expect_is(thepdfandcdf3d_3,"data.frame")
  expect_equal(colnames(thepdfandcdf3d_3),c("theta","pdf","cdf","transparam","pdf_transparam"))
  expect_true(all(thepdfandcdf3d_3$pdf > 0))
  expect_true(all(thepdfandcdf3d_3$pdf_transparam > 0))
  expect_true(all(thepdfandcdf3d_3$cdf >= 0))
  expect_equal(thepdfandcdf3d_3$cdf[1],0)
  expect_equal(round(thepdfandcdf3d_3$cdf[1000],2),1)


  # Quantiles
  expect_equal(names(thequantiles),c("2.5%","97.5%"))
  expect_false(any(is.infinite(thequantiles)))
  expect_lt(thequantiles[1],thequantiles[2])
  expect_lt(abs(exp(thequantiles[1]) - qgamma(.025,1+sum(y1),1+n1)),.07)
  expect_lt(abs(exp(thequantiles[2]) - qgamma(.975,1+sum(y1),1+n1)),.09)

  expect_equal(names(thequantiles3d_1),c("2.5%","97.5%"))
  expect_false(any(is.infinite(thequantiles3d_1)))
  expect_lt(thequantiles3d_1[1],thequantiles3d_1[2])
  expect_lt(abs(exp(thequantiles3d_1[1]) - qgamma(.025,1+sum(y1),1+n1)),.07)
  expect_lt(abs(exp(thequantiles3d_1[2]) - qgamma(.975,1+sum(y1),1+n1)),.09)

  expect_equal(names(thequantiles3d_2),c("2.5%","97.5%"))
  expect_false(any(is.infinite(thequantiles3d_2)))
  expect_lt(thequantiles3d_2[1],thequantiles3d_2[2])
  expect_lt(abs(exp(thequantiles3d_2[1]) - qgamma(.025,1+sum(y2),1+n2)),.07)
  expect_lt(abs(exp(thequantiles3d_2[2]) - qgamma(.975,1+sum(y2),1+n2)),.09)

  expect_equal(names(thequantiles3d_3),c("2.5%","97.5%"))
  expect_false(any(is.infinite(thequantiles3d_3)))
  expect_lt(thequantiles3d_3[1],thequantiles3d_3[2])
  expect_lt(abs(exp(thequantiles3d_3[1]) - qgamma(.025,1+sum(y3),1+n3)),.07)
  expect_lt(abs(exp(thequantiles3d_3[2]) - qgamma(.975,1+sum(y3),1+n3)),.09)

  # Spline interp
  expect_warning(compute_pdf_and_cdf(thequadrature,transformation = list(totheta = log,fromtheta = exp),interpolation = 'spline')) # Too few points to use spline
  # Unstable
  expect_equal(max(pdf_poly_2d[[1]]$pdf_transparam),Inf)
  expect_equal(min(pdf_poly_2d[[1]]$cdf),0,tolerance = 1e-05)
  expect_equal(max(pdf_poly_2d[[1]]$cdf),Inf)
  # Stable
  expect_lt(max(pdf_spline_2d[[1]]$pdf_transparam),1)
  expect_equal(min(pdf_spline_2d[[1]]$cdf),0,tolerance = 1e-05)
  expect_equal(max(pdf_spline_2d[[1]]$cdf),1,tolerance = 1e-05)

  # Spline interp and sampling
  expect_lt(suppressWarnings(abs(ks.test(polysamps[[1]],splinesamps[[1]])$statistic)),.05)
  expect_lt(suppressWarnings(abs(ks.test(polysamps[[2]],splinesamps[[2]])$statistic)),.05)

})

