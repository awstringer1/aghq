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

  # New marginal posteriors
  expect_is(margpost_1d_1_correct,"data.frame")
  expect_is(margpost_2d_1_correct,"data.frame")
  expect_is(margpost_2d_2_correct,"data.frame")
  expect_is(margpost_3d_1_correct,"data.frame")
  expect_is(margpost_3d_2_correct,"data.frame")
  expect_is(margpost_3d_3_correct,"data.frame")


  expect_equal(colnames(margpost_1d_1_correct),c("theta1","logmargpost"))
  expect_equal(colnames(margpost_2d_1_correct),c("theta1","logmargpost"))
  expect_equal(colnames(margpost_2d_2_correct),c("theta2","logmargpost"))
  expect_equal(colnames(margpost_3d_1_correct),c("theta1","logmargpost"))
  expect_equal(colnames(margpost_3d_2_correct),c("theta2","logmargpost"))
  expect_equal(colnames(margpost_3d_3_correct),c("theta3","logmargpost"))


  expect_equal(nrow(margpost_1d_1_correct),3)
  expect_equal(nrow(margpost_2d_1_correct),3)
  expect_equal(nrow(margpost_2d_2_correct),3)
  expect_equal(nrow(margpost_2d_2_k7_correct),7)
  expect_equal(nrow(margpost_3d_1_correct),3)
  expect_equal(nrow(margpost_3d_2_correct),3)
  expect_equal(nrow(margpost_3d_2_correct),3)


  # Integration
  expect_lt(sum(abs(range(compute_pdf_and_cdf(margpost_1d_1_correct)$cdf) - c(0,1))),1e-03)
  expect_lt(sum(abs(range(compute_pdf_and_cdf(margpost_2d_1_correct)$cdf) - c(0,1))),1e-03)
  expect_lt(sum(abs(range(compute_pdf_and_cdf(margpost_2d_2_correct)$cdf) - c(0,1))),1e-03)
  expect_lt(sum(abs(range(compute_pdf_and_cdf(margpost_2d_2_k7_correct)$cdf) - c(0,1))),1e-03)
  expect_lt(sum(abs(range(compute_pdf_and_cdf(margpost_3d_1_correct)$cdf) - c(0,1))),1e-03)
  expect_lt(sum(abs(range(compute_pdf_and_cdf(margpost_3d_2_correct)$cdf) - c(0,1))),1e-03)
  expect_lt(sum(abs(range(compute_pdf_and_cdf(margpost_3d_3_correct)$cdf) - c(0,1))),1e-03)

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

  # Updated interface for moments
  expect_equal(compute_moment(thequadrature),1)
  expect_equal(compute_moment(thequadrature,method='correct'),1)
  expect_equal(compute_moment(thequadrature,ff = function(theta) rep(1,2)),rep(1,2))
  expect_equal(compute_moment(thequadrature,ff = function(theta) rep(1,10)),rep(1,10))
  expect_equal(compute_moment(thequadrature,ff = function(theta) 1:2),1:2)
  expect_equal(compute_moment(thequadrature,ff = function(theta) 1:10),1:10)

  expect_equal(compute_moment(thequadrature,ff = function(theta) rep(1,2),method='correct'),rep(1,2))
  expect_equal(compute_moment(thequadrature,ff = function(theta) rep(1,10),method='correct'),rep(1,10))
  expect_equal(compute_moment(thequadrature,ff = function(theta) 1:2,method='correct'),1:2)
  expect_equal(compute_moment(thequadrature,ff = function(theta) 1:10,method='correct'),1:10)

  expect_equal(compute_moment(thequadrature,gg = make_moment_function(function(x) log(2))),2)
  expect_equal(compute_moment(thequadrature,ff=NULL,gg = make_moment_function(function(x) log(2))),2)
  expect_equal(compute_moment(thequadrature,gg = make_moment_function(function(x) log(2)),method='correct'),2)
  expect_equal(compute_moment(thequadrature,ff=NULL,gg = make_moment_function(function(x) log(2)),method='correct'),2)
  expect_equal(compute_moment(thequadrature,gg = make_moment_function(function(x) 0),method='correct'),1)
  expect_equal(compute_moment(thequadrature,gg = make_moment_function(function(theta) log(rep(1,2)))),rep(1,2))
  expect_equal(compute_moment(thequadrature,gg = make_moment_function(function(theta) log(rep(1,10)))),rep(1,10))
  expect_equal(compute_moment(thequadrature,gg = make_moment_function(function(x) log(2)),method='correct'),2)

  expect_lt(abs(aghqexpmean1d_correct - trueexpmean1d),1e-05)
  expect_lt(abs(aghqexpmean2d_correct - trueexpmean2d)[1],1e-04)
  expect_lt(abs(aghqexpmean2d_correct - trueexpmean2d)[2],1e-04)

  expect_lt(abs(aghqexpmean3d_correct - trueexpmean3d)[1],1e-04)
  expect_lt(abs(aghqexpmean3d_correct - trueexpmean3d)[2],1e-04)
  expect_lt(abs(aghqexpmean3d_correct - trueexpmean3d)[3],1e-04)

  expect_lt(abs(aghqexpmean1d_correct2_1 - trueexpmean1d),1e-05)
  expect_lt(abs(aghqexpmean2d_correct2_1 - trueexpmean2d[1]),1e-05)
  expect_lt(abs(aghqexpmean2d_correct2_2 - trueexpmean2d[2]),1e-05)
  expect_lt(abs(aghqexpmean3d_correct2_1 - trueexpmean3d[1]),1e-05)
  expect_lt(abs(aghqexpmean3d_correct2_2 - trueexpmean3d[2]),1e-05)
  expect_lt(abs(aghqexpmean3d_correct2_3 - trueexpmean3d[3]),1e-05)


  expect_warning(compute_moment(thequadrature$normalized_posterior,method='correct'))




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
  expect_equal(colnames(thepdfandcdf),pdfandcdfnames)
  expect_true(all(thepdfandcdf$pdf > 0))
  expect_true(all(thepdfandcdf$cdf >= 0))
  expect_equal(thepdfandcdf$cdf[1],0)
  expect_equal(round(thepdfandcdf$cdf[1000],2),1)


  # PDF and CDF with transformation
  expect_is(pdfwithtrans,"data.frame")
  expect_equal(colnames(pdfwithtrans),pdfandcdfnames)
  expect_true(all(pdfwithtrans$pdf > 0))
  expect_true(all(pdfwithtrans$pdf_transparam > 0))
  expect_true(all(pdfwithtrans$cdf >= 0))
  expect_equal(pdfwithtrans$cdf[1],0)
  expect_equal(round(pdfwithtrans$cdf[1000],2),1)

  # 3d
  expect_is(thepdfandcdf3d_1,"data.frame")
  expect_equal(colnames(thepdfandcdf3d_1),pdfandcdfnames)
  expect_true(all(thepdfandcdf3d_1$pdf > 0))
  expect_true(all(thepdfandcdf3d_1$pdf_transparam > 0))
  expect_true(all(thepdfandcdf3d_1$cdf >= 0))
  expect_equal(thepdfandcdf3d_1$cdf[1],0)
  expect_equal(round(thepdfandcdf3d_1$cdf[1000],2),1)

  expect_is(thepdfandcdf3d_2,"data.frame")
  expect_equal(colnames(thepdfandcdf3d_2),pdfandcdfnames)
  expect_true(all(thepdfandcdf3d_2$pdf > 0))
  expect_true(all(thepdfandcdf3d_2$pdf_transparam > 0))
  expect_true(all(thepdfandcdf3d_2$cdf >= 0))
  expect_equal(thepdfandcdf3d_2$cdf[1],0)
  expect_equal(round(thepdfandcdf3d_2$cdf[1000],2),1)

  expect_is(thepdfandcdf3d_3,"data.frame")
  expect_equal(colnames(thepdfandcdf3d_3),pdfandcdfnames)
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
  # expect_gt(max(pdf_poly_2d[[1]]$pdf_transparam),1e10) # Arbitrary large threshold
  expect_equal(min(pdf_poly_2d[[1]]$cdf),0,tolerance = 1e-05)
  # expect_gt(max(pdf_poly_2d[[1]]$cdf),1e10)
  # Stable
  expect_lt(max(pdf_spline_2d[[1]]$pdf_transparam),1)
  expect_equal(min(pdf_spline_2d[[1]]$cdf),0,tolerance = 1e-05)
  expect_equal(max(pdf_spline_2d[[1]]$cdf),1,tolerance = 1e-05)

  # Spline interp and sampling
  expect_lte(suppressWarnings(abs(ks.test(polysamps[[1]],splinesamps[[1]])$statistic)),.05)
  # expect_lte(suppressWarnings(abs(ks.test(polysamps[[2]],splinesamps[[2]])$statistic)),.05)

  # Parallel sampling
  expect_true(all(themargsamps$samps == themargsamps_parallel$samps))
  expect_true(all(themargsamps$theta == themargsamps_parallel$theta))
  expect_true(all(themargsamps$thetasamples[[1]] == themargsamps_parallel$thetasamples[[1]]))

  # Transformations
  expect_equal(names(trans2),transnames)
  expect_equal(names(trans3),transnames)
  expect_equal(names(trans4),transnames)
  expect_is(trans2,"aghqtrans")
  expect_is(trans3,"aghqtrans")
  expect_is(trans4,"aghqtrans")

  expect_equal(trans1$totheta(tt),trans2$totheta(tt))
  expect_equal(trans2$totheta(tt),trans3$totheta(tt))
  expect_equal(trans3$totheta(tt),trans4$totheta(tt))
  expect_equal(trans1$fromtheta(tt),trans2$fromtheta(tt))
  expect_equal(trans2$fromtheta(tt),trans3$fromtheta(tt))
  expect_equal(trans3$fromtheta(tt),trans4$fromtheta(tt))
  expect_equal(trans2$jacobian(tt),trans3$jacobian(tt))
  expect_equal(trans3$jacobian(tt),trans4$jacobian(tt))

  expect_true(validate_transformation(default_transformation()))
  expect_true(validate_transformation(trans1))
  expect_true(validate_transformation(trans2))
  expect_true(validate_transformation(trans3))
  expect_true(validate_transformation(trans4))
  expect_error(validate_transformation(list(foo = 'bar')))
  expect_error(validate_transformation(t3,checkinverse = checkvals))

  # Moments
  expect_equal(names(mom1),momnames)
  expect_equal(names(mom2),momnames)
  expect_equal(names(mom3),momnames)
  expect_is(mom1,'aghqmoment')
  expect_is(mom2,'aghqmoment')
  expect_is(mom3,'aghqmoment')
  expect_true(validate_moment(mom1))
  expect_true(validate_moment(mom2))
  expect_true(validate_moment(mom3))
  expect_true(validate_moment(exp))
  expect_true(validate_moment('exp'))
  expect_true(validate_moment(list(fn=function(x) x,gr=function(x) 1,he = function(x) 0)))

  expect_error(make_moment_function(mombad1))
  expect_error(make_moment_function(mombad2))
  expect_error(validate_moment(mombad1))
  expect_error(validate_moment(mombad2))
  expect_error(validate_moment(mombad3,checkpositive = c(-1,0,1)))
  expect_error(validate_moment(NULL))

  expect_error(compute_moment(norm_sparse_2d_7,method = 'foo'))
  expect_warning(compute_moment(norm_sparse_2d_7,method = 'correct'))
  expect_error(compute_moment(norm_sparse_2d_7,method = 'foo'))


})

