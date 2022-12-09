# devtools::load_all()
# options(mc.cores = 2)
logfteta <- function(eta,y) {
  sum(y) * eta - (length(y) + 1) * exp(eta) - sum(lgamma(y+1)) + eta
}
truelogint <- function(y) lgamma(1 + sum(y)) - (1 + sum(y)) * log(length(y) + 1) - sum(lgamma(y+1))
set.seed(84343124)
y <- rpois(10,5) # Mode should be sum(y) / (10 + 1)
truemode <- log((sum(y) + 1)/(length(y) + 1))
truelognormconst <- truelogint(y)

objfunc <- function(x) logfteta(x,y)
funlist <- list(
  fn = objfunc,
  gr = function(x) numDeriv::grad(objfunc,x),
  he = function(x) numDeriv::hessian(objfunc,x)
)

opt_sparsetrust <- optimize_theta(funlist,1.5,control = default_control(method = "sparse_trust"))
opt_sr1 <- optimize_theta(funlist,1.5,control = default_control(method = "SR1"))
opt_trust <- optimize_theta(funlist,1.5,control = default_control(method = "trust"))
opt_bfgs <- optimize_theta(funlist,1.5,control = default_control(method = "BFGS"))

norm_sparse_1 <- normalize_logpost(opt_sparsetrust,1,1)
norm_sparse_3 <- normalize_logpost(opt_sparsetrust,3,1)
norm_sparse_5 <- normalize_logpost(opt_sparsetrust,5,1)
norm_sparse_7 <- normalize_logpost(opt_sparsetrust,7,1)

norm_trust_1 <- normalize_logpost(opt_trust,3,1)
norm_trust_3 <- normalize_logpost(opt_trust,3,1)
norm_trust_5 <- normalize_logpost(opt_trust,5,1)
norm_trust_7 <- normalize_logpost(opt_trust,7,1)

norm_bfgs_1 <- normalize_logpost(opt_bfgs,3,1)
norm_bfgs_3 <- normalize_logpost(opt_bfgs,3,1)
norm_bfgs_5 <- normalize_logpost(opt_bfgs,5,1)
norm_bfgs_7 <- normalize_logpost(opt_bfgs,7,1)

norm_sr1_1 <- normalize_logpost(opt_sr1,3,1)
norm_sr1_3 <- normalize_logpost(opt_sr1,3,1)
norm_sr1_5 <- normalize_logpost(opt_sr1,5,1)
norm_sr1_7 <- normalize_logpost(opt_sr1,7,1)

margpost_3 <- marginal_posterior(opt_sparsetrust,3,1)
pdfwithtrans <- compute_pdf_and_cdf(
  margpost_3,
  transformation = list(
    totheta = function(x) log(x),
    fromtheta = function(x) exp(x)
  ))






## A 2d example ##
logfteta2d <- function(eta,y) {
  # eta is now (eta1,eta2)
  # y is now (y1,y2)
  n <- length(y)
  n1 <- ceiling(n/2)
  n2 <- floor(n/2)
  y1 <- y[1:n1]
  y2 <- y[(n1+1):(n1+n2)]
  eta1 <- eta[1]
  eta2 <- eta[2]
  sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
    sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
}
set.seed(84343124)
n1 <- 5
n2 <- 5
n <- n1+n2
y1 <- rpois(n1,5)
y2 <- rpois(n2,5)
truemode2d <- c(log((sum(y1) + 1)/(length(y1) + 1)),log((sum(y2) + 1)/(length(y2) + 1)))
truelognormconst2d <- truelogint(y1) + truelogint(y2)

objfunc2d <- function(x) logfteta2d(x,c(y1,y2))
funlist2d <- list(
  fn = objfunc2d,
  gr = function(x) numDeriv::grad(objfunc2d,x),
  he = function(x) numDeriv::hessian(objfunc2d,x)
)

opt_sparsetrust_2d <- optimize_theta(funlist2d,c(1.5,1.5),control = default_control(method = "sparse_trust"))
opt_trust_2d <- optimize_theta(funlist2d,c(1.5,1.5),control = default_control(method = "trust"))
opt_sr1_2d <- optimize_theta(funlist2d,c(1.5,1.5),control = default_control(method = "SR1"))
opt_bfgs_2d <- optimize_theta(funlist2d,c(1.5,1.5),control = default_control(method = "BFGS"))

norm_sparse_2d_1 <- normalize_logpost(opt_sparsetrust_2d,1,1)
norm_sparse_2d_3 <- normalize_logpost(opt_sparsetrust_2d,3,1)
norm_sparse_2d_5 <- normalize_logpost(opt_sparsetrust_2d,5,1)
norm_sparse_2d_7 <- normalize_logpost(opt_sparsetrust_2d,7,1)

norm_trust_2d_1 <- normalize_logpost(opt_trust_2d,1,1)
norm_trust_2d_3 <- normalize_logpost(opt_trust_2d,3,1)
norm_trust_2d_5 <- normalize_logpost(opt_trust_2d,5,1)
norm_trust_2d_7 <- normalize_logpost(opt_trust_2d,7,1)

norm_bfgs_2d_1 <- normalize_logpost(opt_bfgs_2d,1,1)
norm_bfgs_2d_3 <- normalize_logpost(opt_bfgs_2d,3,1)
norm_bfgs_2d_5 <- normalize_logpost(opt_bfgs_2d,5,1)
norm_bfgs_2d_7 <- normalize_logpost(opt_bfgs_2d,7,1)

norm_sr1_2d_1 <- normalize_logpost(opt_sr1_2d,1,1)
norm_sr1_2d_3 <- normalize_logpost(opt_sr1_2d,3,1)
norm_sr1_2d_5 <- normalize_logpost(opt_sr1_2d,5,1)
norm_sr1_2d_7 <- normalize_logpost(opt_sr1_2d,7,1)

# Parameter vector reordering
norm_sparse_2d_reorder_3 <- normalize_logpost(opt_sparsetrust_2d,3,2)


# Marginal posteriors
margpost_1d_1 <- marginal_posterior(opt_sparsetrust,3,1)
margpost_2d_1 <- marginal_posterior(opt_sparsetrust_2d,3,1)
margpost_2d_2 <- marginal_posterior(opt_sparsetrust_2d,3,2)
margpost_2d_2_k7 <- marginal_posterior(opt_sparsetrust_2d,7,2)

# Moments
truemean1d <- truemode # note: this isn't exactly right...
truemean2d <- truemode2d
trueexpmean1d <- exp(truemean1d) # ...but this is
trueexpmean2d <- exp(truemean2d)
truesd1d <- sqrt(1+sum(y)) / (1 + length(y))
truesd2d <- c(sqrt(1+sum(y1)) / (1 + length(y1)),sqrt(1+sum(y2)) / (1 + length(y2)))

aghqnormconst1d <- compute_moment(norm_sparse_7)
aghqnormconst2d <- compute_moment(norm_sparse_2d_7)


aghqmean1d <- compute_moment(norm_sparse_7,function(x) x)
aghqmean2d <- compute_moment(norm_sparse_2d_7,function(x) x)

aghqexpmean1d <- compute_moment(norm_sparse_7,function(x) exp(x))
aghqexpmean2d <- compute_moment(norm_sparse_2d_7,function(x) exp(x))

aghqexpsd1d <- sqrt(compute_moment(norm_sparse_7,function(x) (exp(x) - trueexpmean1d)^2))
aghqexpsd2d_1 <- sqrt(compute_moment(norm_sparse_2d_7,function(x) (exp(x) - trueexpmean2d[1])^2))[1]
aghqexpsd2d_2 <- sqrt(compute_moment(norm_sparse_2d_7,function(x) (exp(x) - trueexpmean2d[2])^2))[2]




# Interpolation
margpostinterp <- interpolate_marginal_posterior(margpost_2d_1)
margpostinterp_2 <- interpolate_marginal_posterior(margpost_2d_2)

# pdf and cdf
pdfandcdfnames <- c("theta","pdf","cdf","transparam","pdf_transparam")
thepdfandcdf <- compute_pdf_and_cdf(margpost_2d_1)

# quantiles
thequantiles <- compute_quantiles(margpost_2d_1)
exp(thequantiles)
qgamma(.025,1+sum(y1),1+n1)
qgamma(.975,1+sum(y1),1+n1)

# Quadrature!
thequadrature <- aghq(funlist2d,3,c(0,0))
thesummary <- summary(thequadrature)

# Test new S3 interface for moments
aghqnormconst2d_new <- compute_moment(thequadrature)

aghqmean2d_new <- compute_moment(thequadrature,function(x) x)

# Test the two types of interpolation
thequadrature_k7 <- aghq(funlist2d,25,c(0,0)) # 25 quad points leads to bad poly interp
pdf_poly_2d <- compute_pdf_and_cdf(thequadrature_k7,transformation = list(totheta = log,fromtheta = exp),interpolation='polynomial')
pdf_spline_2d <- compute_pdf_and_cdf(thequadrature_k7,transformation = list(totheta = log,fromtheta = exp),interpolation = 'spline')
# Sampling
set.seed(708968)
polysamps <- sample_marginal(thequadrature_k7,1e03)
splinesamps <- sample_marginal(thequadrature_k7,1e03,interpolation = 'spline')
if (FALSE) {
  par(mfrow = c(1,2))
  hist(polysamps[[1]],breaks = 50,freq=FALSE)
  hist(splinesamps[[1]],breaks = 50,freq=FALSE)

}

# Laplace

thelaplace <- laplace_approximation(funlist2d,c(0,0))

# Marginal laplace...

objfunc2dmarg <- function(W,theta) objfunc2d(c(W,theta))
objfunc2dmarggr <- function(W,theta) {
  fn <- function(W) objfunc2dmarg(W,theta)
  numDeriv::grad(fn,W)
}
objfunc2dmarghe <- function(W,theta) {
  fn <- function(W) objfunc2dmarg(W,theta)
  numDeriv::hessian(fn,W)
}

funlist2dmarg <- list(
  fn = objfunc2dmarg,
  gr = objfunc2dmarggr,
  he = objfunc2dmarghe
)

themarginallaplace <- aghq::marginal_laplace(funlist2dmarg,3,list(W = 0,theta = 0))
set.seed(7809685)
themargsamps <- aghq::sample_marginal(themarginallaplace,10)
set.seed(7809685)
themargsamps_parallel <- aghq::sample_marginal(themarginallaplace,10,numcores = 2)

## A 3d example ##
# This is necessary because I want a 2d marginal example
logfteta3d <- function(eta,y) {
  # eta is now (eta1,eta2)
  # y is now (y1,y2)
  n <- length(y)
  n1 <- ceiling(n/3)
  n2 <- ceiling(n/3)
  n3 <- n - (n1+n2)

  y1 <- y[1:n1]
  y2 <- y[(n1+1):(n1+n2)]
  y3 <- y[(n1+n2+1):(n1+n2+n3)]

  eta1 <- eta[1]
  eta2 <- eta[2]
  eta3 <- eta[3]

  sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
    sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2 +
      sum(y3) * eta3 - (length(y3) + 1) * exp(eta3) - sum(lgamma(y3+1)) + eta3

}
set.seed(84343124)
n1 <- 5
n2 <- 5
n3 <- 5
n <- n1+n2+n3
y1 <- rpois(n1,5)
y2 <- rpois(n2,5)
y3 <- rpois(n3,5)

truemode3d <- c(log((sum(y1) + 1)/(length(y1) + 1)),log((sum(y2) + 1)/(length(y2) + 1)),log((sum(y3) + 1)/(length(y3) + 1)))
truelognormconst3d <- truelogint(y1)  + truelogint(y2) + truelogint(y3)

objfunc3d <- function(x) logfteta3d(x,c(y1,y2,y3))
funlist3d <- list(
  fn = objfunc3d,
  gr = function(x) numDeriv::grad(objfunc3d,x),
  he = function(x) numDeriv::hessian(objfunc3d,x)
)


opt_sparsetrust_3d <- optimize_theta(funlist3d,c(1.5,1.5,1.5),control = default_control(method = "sparse_trust"))
opt_trust_3d <- optimize_theta(funlist3d,c(1.5,1.5,1.5),control = default_control(method = "trust"))
opt_sr1_3d <- optimize_theta(funlist3d,c(1.5,1.5,1.5),control = default_control(method = "SR1"))
opt_bfgs_3d <- optimize_theta(funlist3d,c(1.5,1.5,1.5),control = default_control(method = "BFGS"))

norm_sparse_3d_1 <- normalize_logpost(opt_sparsetrust_3d,1,1)
norm_sparse_3d_3 <- normalize_logpost(opt_sparsetrust_3d,3,1)
norm_sparse_3d_5 <- normalize_logpost(opt_sparsetrust_3d,5,1)
norm_sparse_3d_7 <- normalize_logpost(opt_sparsetrust_3d,7,1)

norm_trust_3d_1 <- normalize_logpost(opt_trust_3d,1,1)
norm_trust_3d_3 <- normalize_logpost(opt_trust_3d,3,1)
norm_trust_3d_5 <- normalize_logpost(opt_trust_3d,5,1)
norm_trust_3d_7 <- normalize_logpost(opt_trust_3d,7,1)

norm_bfgs_3d_1 <- normalize_logpost(opt_bfgs_3d,1,1)
norm_bfgs_3d_3 <- normalize_logpost(opt_bfgs_3d,3,1)
norm_bfgs_3d_5 <- normalize_logpost(opt_bfgs_3d,5,1)
norm_bfgs_3d_7 <- normalize_logpost(opt_bfgs_3d,7,1)

norm_sr1_3d_1 <- normalize_logpost(opt_sr1_3d,1,1)
norm_sr1_3d_3 <- normalize_logpost(opt_sr1_3d,3,1)
norm_sr1_3d_5 <- normalize_logpost(opt_sr1_3d,5,1)
norm_sr1_3d_7 <- normalize_logpost(opt_sr1_3d,7,1)

# Parameter vector reordering
norm_sparse_3d_reorder_1 <- normalize_logpost(opt_sparsetrust_3d,3,1)
norm_sparse_3d_reorder_2 <- normalize_logpost(opt_sparsetrust_3d,3,2)
norm_sparse_3d_reorder_3 <- normalize_logpost(opt_sparsetrust_3d,3,3)

# Marginal posteriors
margpost_3d_1 <- marginal_posterior(opt_sparsetrust_3d,3,1)
margpost_3d_2 <- marginal_posterior(opt_sparsetrust_3d,3,2)
margpost_3d_3 <- marginal_posterior(opt_sparsetrust_3d,3,3)

margpostinterp_3d_1 <- interpolate_marginal_posterior(margpost_3d_1)
margpostinterp_3d_2 <- interpolate_marginal_posterior(margpost_3d_2)
margpostinterp_3d_3 <- interpolate_marginal_posterior(margpost_3d_3)


# Moments
truemean3d <- truemode3d
trueexpmean3d <- exp(truemean3d)
truesd3d <- c(sqrt(1+sum(y1)) / (1 + length(y1)),sqrt(1+sum(y2)) / (1 + length(y2)),sqrt(1+sum(y3)) / (1 + length(y3)))

aghqnormconst3d <- compute_moment(norm_sparse_3d_3)

aghqmean3d <- compute_moment(norm_sparse_3d_7,function(x) x)

aghqexpmean3d <- compute_moment(norm_sparse_3d_7,function(x) exp(x))

aghqexpsd3d_1 <- sqrt(compute_moment(norm_sparse_3d_7,function(x) (exp(x) - trueexpmean3d[1])^2))[1]
aghqexpsd3d_2 <- sqrt(compute_moment(norm_sparse_3d_7,function(x) (exp(x) - trueexpmean3d[2])^2))[2]
aghqexpsd3d_3 <- sqrt(compute_moment(norm_sparse_3d_7,function(x) (exp(x) - trueexpmean3d[3])^2))[3]

momquad <- aghq(funlist,9,0)
momquad2 <- aghq(funlist2d,9,c(0,0))
momquad3 <- aghq(funlist3d,9,c(0,0,0))

aghqexpmean1d_correct <- compute_moment(momquad,function(x) exp(x),method='correct')
aghqexpmean2d_correct <- compute_moment(momquad2,function(x) exp(x),method='correct')
aghqexpmean3d_correct <- compute_moment(momquad3,function(x) exp(x),method='correct')
aghqexpmean1d_correct2_1 <- compute_moment(momquad,gg = make_moment_function(function(x) x,method='correct'))
aghqexpmean2d_correct2_1 <- compute_moment(momquad2,gg = make_moment_function(function(x) x[1]),method='correct')
aghqexpmean2d_correct2_2 <- compute_moment(momquad2,gg = make_moment_function(function(x) x[2]),method='correct')
aghqexpmean3d_correct2_1 <- compute_moment(momquad3,gg = make_moment_function(function(x) x[1]),method='correct')
aghqexpmean3d_correct2_2 <- compute_moment(momquad3,gg = make_moment_function(function(x) x[2]),method='correct')
aghqexpmean3d_correct2_3 <- compute_moment(momquad3,gg = make_moment_function(function(x) x[3]),method='correct')

margquad <- aghq(funlist,3,0)
margquad2 <- aghq(funlist2d,3,c(0,0))
margquad2_k7 <- aghq(funlist2d,7,c(0,0))
margquad3 <- aghq(funlist3d,3,c(0,0,0))

margpost_1d_1_correct <- marginal_posterior(margquad,1,method='correct')
margpost_2d_1_correct <- marginal_posterior(margquad2,1,method='correct')
margpost_2d_2_correct <- marginal_posterior(margquad2,2,method='correct')
margpost_2d_2_k7_correct <- marginal_posterior(margquad2_k7,2,method='correct')
margpost_3d_1_correct <- marginal_posterior(margquad3,1,method='correct')
margpost_3d_2_correct <- marginal_posterior(margquad3,2,method='correct')
margpost_3d_3_correct <- marginal_posterior(margquad3,3,method='correct')

# check the correction of marginals
margpost_thequadrature_original <- thequadrature$marginals
thequadrature_correct <- correct_marginals(thequadrature)
margpost_thequadrature_correct <- thequadrature_correct$marginals


# Interpolation
margpostinterp3d_1 <- interpolate_marginal_posterior(margpost_3d_1)
margpostinterp3d_2 <- interpolate_marginal_posterior(margpost_3d_2)
margpostinterp3d_3 <- interpolate_marginal_posterior(margpost_3d_3)


# pdf and cdf
thepdfandcdf3d_1 <- compute_pdf_and_cdf(margpost_3d_1,transformation = list(totheta = log,fromtheta = exp))
thepdfandcdf3d_2 <- compute_pdf_and_cdf(margpost_3d_2,transformation = list(totheta = log,fromtheta = exp))
thepdfandcdf3d_3 <- compute_pdf_and_cdf(margpost_3d_3,transformation = list(totheta = log,fromtheta = exp))

# quantiles
thequantiles3d_1 <- compute_quantiles(margpost_3d_1)
thequantiles3d_2 <- compute_quantiles(margpost_3d_2)
thequantiles3d_3 <- compute_quantiles(margpost_3d_3)

# Quadrature!
thequadrature3d <- aghq(funlist3d,3,c(0,0,0))
thesummary3d <- summary(thequadrature3d)


# Laplace

thelaplace3d <- laplace_approximation(funlist3d,c(0,0,0))

# Marginal laplace...

objfunc3dmarg <- function(W,theta) objfunc3d(c(W,theta))
objfunc3dmarggr <- function(W,theta) {
  fn <- function(W) objfunc3dmarg(W,theta)
  numDeriv::grad(fn,W)
}
objfunc3dmarghe <- function(W,theta) {
  fn <- function(W) objfunc3dmarg(W,theta)
  numDeriv::hessian(fn,W)
}

funlist3dmarg <- list(
  fn = objfunc3dmarg,
  gr = objfunc3dmarggr,
  he = objfunc3dmarghe
)

# 1-d marglaplace
themarginallaplace3d_1 <- aghq::marginal_laplace(funlist3dmarg,3,list(W = c(0,0),theta = 0))
themargsamps3d_1 <- aghq::sample_marginal(themarginallaplace3d_1,10)

# 2d marglaplace
themarginallaplace3d_2 <- aghq::marginal_laplace(funlist3dmarg,3,list(W = 0,theta = c(0,0)))
themargsamps3d_2 <- aghq::sample_marginal(themarginallaplace3d_2,10)



## Sparse Grids ##
# doesn't make sense in 1d, do 2d
sparsegrid_2d <- aghq(funlist2d,5,c(0,0),control = default_control(ndConstruction = 'sparse'))
sparsenormconst_2d <- compute_moment(sparsegrid_2d) # Should be 1

## Control Params ##
cntrl_base <- default_control()
cntrl_marg <- default_control_marglaplace()
cntrl_tmb <- default_control_tmb()

# Laplace approx
logint1 <- function(x,n) n*log(x) - x
logaghq <- function(n) {
  ff <- list(
    fn = function(x) logint1(x,n),
    gr = function(x) numDeriv::grad(logint1,x,n=n),
    he = function(x) numDeriv::hessian(logint1,x,n=n)
  )
  aghq::laplace_approximation(ff,n)$lognormconst
}
logstirling <- function(n) (1/2)*log(2*pi) + (1/2)*log(n) + n*(log(n) - 1)
la5 <- round(logaghq(5),11)
ls5 <- round(logstirling(5),11)
la10 <- round(logaghq(10),11)
ls10 <- round(logstirling(10),11)
la100 <- round(logaghq(100),11)
ls100 <- round(logstirling(100),11)

## Custom grid ----
# GHQ, 1d
gg1 <- mvQuad::createNIGrid(1,'GHe',5)
aghq_customgrid_gg1 <- aghq(funlist,5,0,basegrid = gg1)
aghq_customgrid_auto1 <- aghq(funlist,5,0,basegrid = NULL)

# GHQ, 2d
gg2 <- mvQuad::createNIGrid(2,'GHe',5)
aghq_customgrid_gg2 <- aghq(funlist2d,5,c(0,0),basegrid = gg2)
aghq_customgrid_auto2 <- aghq(funlist2d,5,c(0,0),basegrid = NULL)

# GHQ, 2d, sparse grid
gg2s <- mvQuad::createNIGrid(2,'GHe',5,ndConstruction = 'sparse')
aghq_customgrid_gg2s <- aghq(funlist2d,5,c(0,0),basegrid = gg2s)
aghq_customgrid_auto2s <- aghq(funlist2d,5,c(0,0),basegrid = NULL,control = default_control(ndConstruction = 'sparse'))

# Non-GHQ, 2d, but still Gaussian
gg3 <- mvQuad::createNIGrid(2,'nHe',5)
aghq_customgrid_gg3 <- aghq(funlist2d,5,c(0,0),basegrid = gg2)

# Non-Gaussian kernel, should throw an error
gg4 <- mvQuad::createNIGrid(2,'GLe',5)

# Create a new grid for later check
gg5 <- mvQuad::createNIGrid(2,'GHe',5)
# Check that not providing k works
gg6 <- mvQuad::createNIGrid(2,'GHe',5)
aghq_customgrid_gg6 <- aghq(funlist2d,startingvalue = c(0,0),basegrid = gg6)

# Check that the grid isn't modified
gg7 <- mvQuad::createNIGrid(1,'GHe',5)



## Extraction of log normalizing constants
normconst1 <- get_log_normconst(thequadrature)
normconst2 <- get_log_normconst(thelaplace)
normconst3 <- get_log_normconst(themarginallaplace)

## Optimization: controls work

funlist3dneg <- with(funlist3d,list(
  fn = function(x) -1*fn(x),
  gr = function(x) -1*gr(x),
  he = function(x) -1*he(x)
))
opt_controlworks1 <- optimize_theta(funlist3d,c(0,0,0))
# Negate
opt_controlworks2 <- optimize_theta(funlist3dneg,c(0,0,0),control=default_control(negate=TRUE))
# Numeric hessian
funlist3dnohess <- funlist3dneg
funlist3dnohess$he <- NULL
opt_controlworks3 <- optimize_theta(funlist3d,c(0,0,0),control = default_control_tmb(negate=FALSE,numhessian = TRUE))

# but now, make sure aghq still works!
aghq_controlworks1 <- aghq(funlist3d,5,c(0,0,0))
aghq_controlworks2 <- aghq(funlist3dneg,5,c(0,0,0),control=default_control(negate=TRUE))
aghq_controlworks3 <- aghq(funlist3d,5,c(0,0,0),control = default_control(negate=FALSE,numhessian = TRUE))



## Control argument validation
goodcontrol_aghq <- default_control()
goodcontrol_marglaplace <- default_control_marglaplace()
goodcontrol_tmb <- default_control_tmb()

badcontrol1_aghq <- goodcontrol_aghq
badcontrol1_aghq$negate <- NULL
badcontrol2_aghq <- goodcontrol_aghq
badcontrol2_aghq$foo <- 'bar'
badcontrol3_aghq <- goodcontrol_aghq
badcontrol3_aghq$foo <- 'bar'
badcontrol3_aghq$negate <- NULL

badcontrol1_marglaplace <- goodcontrol_marglaplace
badcontrol1_marglaplace$negate <- NULL
badcontrol2_marglaplace <- goodcontrol_marglaplace
badcontrol2_marglaplace$foo <- 'bar'
badcontrol3_marglaplace <- goodcontrol_marglaplace
badcontrol3_marglaplace$foo <- 'bar'
badcontrol3_marglaplace$negate <- NULL

badcontrol1_tmb <- goodcontrol_tmb
badcontrol1_tmb$negate <- NULL
badcontrol2_tmb <- goodcontrol_tmb
badcontrol2_tmb$foo <- 'bar'
badcontrol3_tmb <- goodcontrol_tmb
badcontrol3_tmb$foo <- 'bar'
badcontrol3_tmb$negate <- NULL

## Test returning only the normconst
aghq_normconst1 <- aghq(funlist3d,5,c(0,0,0),control = default_control(onlynormconst = TRUE))
marglaplace_normconst1 <- aghq::marginal_laplace(funlist2dmarg,3,list(W = 0,theta = 0),control = default_control_marglaplace(onlynormconst = TRUE))

## Test summary for marglaplace
mlsumm1 <- summary(themarginallaplace)
mlsumm2 <- summary(themarginallaplace3d_1)
mlsumm3 <- summary(themarginallaplace,M=100)
mlsumm4 <- summary(themarginallaplace3d_1,max_print=1)

## Transformations ##
transnames <- c("totheta","fromtheta","jacobian")
tt <- exp(rnorm(10))
trans1 <- list(totheta = log,fromtheta = exp)
trans2 <- make_transformation(log,exp)
trans3 <- make_transformation("log","exp")
trans4 <- make_transformation(trans1)

t3 <- make_transformation(log,log)
checkvals <- exp(exp(rnorm(10)))

# thepdfandcdf_trans1 <- compute_pdf_and_cdf(margpost_3d_1,transformation = trans1)
# thepdfandcdf_trans2 <- compute_pdf_and_cdf(margpost_3d_1,transformation = trans2)
# thepdfandcdf_trans3 <- compute_pdf_and_cdf(margpost_3d_1,transformation = trans3)


## Moments ##
momnames <- c("fn","gr","he")
mom1 <- make_moment_function(exp)
mom2 <- make_moment_function('exp')
mom3 <- make_moment_function(list(fn=function(x) x,gr=function(x) 1,he = function(x) 0))
mombad1 <- list(exp,exp,exp) # No names
mombad2 <- list('exp','exp','exp') # List of not functions
mombad3 <- make_moment_function(function(x) NA)

## New integer moments ##
set.seed(4378)
n <- 10
lambda <- 2
y <- rpois(n,lambda)
# POSITIVE mode, with some negative quad points
momobjfunc1 <- function(eta) {
  sum(y) * eta - (length(y) + 1) * exp(eta) - sum(lgamma(y+1)) + eta
}
momfunlist1 <- list(
  fn = momobjfunc1,
  gr = function(x) numDeriv::grad(momobjfunc1,x),
  he = function(x) numDeriv::hessian(momobjfunc1,x)
)
momshiftquad1 <- aghq(momfunlist1,7,0)
truemoment <- digamma(sum(y) + 1) - log(length(y) + 1)
truesecondcentralmoment <- trigamma(sum(y)+1)
truesecondrawmoment <- truesecondcentralmoment + truemoment^2
# Make numeric moment function
ggmomnum_manual <- make_numeric_moment_function(1,1,momshiftquad1,shift = 20)
ggmomnum_auto <- make_numeric_moment_function(1,1,momshiftquad1)

# Old way
nummom_list1 <- compute_moment(momshiftquad1$normalized_posterior,nn=1)
nummom_aghq1 <- compute_moment(momshiftquad1,nn=1) # Should call compute_moment.list
nummom_aghq2 <- compute_moment(momshiftquad1,nn=2)
# Central
nummom_aghq_central1 <- compute_moment(momshiftquad1,nn=1,type='central')
nummom_aghq_central2 <- compute_moment(momshiftquad1,nn=2,type='central')
nummom_aghq2 - nummom_aghq1^2
# New way
nummom_aghq_correct1 <- compute_moment(momshiftquad1,nn=1,method = 'correct')
nummom_aghq_correct2 <- compute_moment(momshiftquad1,nn=2,method = 'correct')
nummom_aghq_correct_central1 <- compute_moment(momshiftquad1,nn=1,method = 'correct',type='central')
nummom_aghq_correct_central2 <- compute_moment(momshiftquad1,nn=2,method = 'correct',type='central')

# POSITIVE mode, with no negative quad points
# There should be NO shift
set.seed(4378)
n <- 100
lambda <- 20
y <- rpois(n,lambda)

momobjfunc2 <- function(eta) {
  sum(y) * eta - (length(y) + 1) * exp(eta) - sum(lgamma(y+1)) + eta
}
momfunlist2 <- list(
  fn = momobjfunc2,
  gr = function(x) numDeriv::grad(momobjfunc2,x),
  he = function(x) numDeriv::hessian(momobjfunc2,x)
)
momshiftquad2 <- aghq(momfunlist2,7,0)
truemoment2 <- digamma(sum(y) + 1) - log(length(y) + 1)
truesecondcentralmoment2 <- trigamma(sum(y)+1)
truesecondrawmoment2 <- truesecondcentralmoment2 + truemoment2^2

# Test that the moment functions are created without shift
ggmomshift2 <- make_numeric_moment_function(1,1,momshiftquad2,0,NULL)
ggmomshift2_withcentre <- make_numeric_moment_function(1,1,momshiftquad2,1,NULL)


# Moments
nummom_aghq_correct1_2 <- compute_moment(momshiftquad2,nn=1,method = 'correct')
nummom_aghq_correct2_2 <- compute_moment(momshiftquad2,nn=2,method = 'correct')
nummom_aghq_correct_central1_2 <- compute_moment(momshiftquad2,nn=1,method = 'correct',type='central')
nummom_aghq_correct_central2_2 <- compute_moment(momshiftquad2,nn=2,method = 'correct',type='central')

# NEGATIVE mode, with some negative quad points
set.seed(4378)
n <- 10
lambda <- 2
y <- rpois(n,lambda)

momobjfunc3 <- function(eta) {
  eta <- -1*eta
  sum(y) * eta - (length(y) + 1) * exp(eta) - sum(lgamma(y+1)) + eta
}
momfunlist3 <- list(
  fn = momobjfunc3,
  gr = function(x) numDeriv::grad(momobjfunc3,x),
  he = function(x) numDeriv::hessian(momobjfunc3,x)
)
momshiftquad3 <- aghq(momfunlist3,7,0)
truemoment3 <- -1*(digamma(sum(y) + 1) - log(length(y) + 1))
truesecondcentralmoment3 <- trigamma(sum(y)+1) # Still positive
truesecondrawmoment3 <- truesecondcentralmoment3 + truemoment3^2

# Moments
nummom_aghq_correct1_3 <- compute_moment(momshiftquad3,nn=1,method = 'correct')
nummom_aghq_correct2_3 <- compute_moment(momshiftquad3,nn=2,method = 'correct')
nummom_aghq_correct_central1_3 <- compute_moment(momshiftquad3,nn=1,method = 'correct',type='central')
nummom_aghq_correct_central2_3 <- compute_moment(momshiftquad3,nn=2,method = 'correct',type='central')

# EXTREMELY NEGATIVE mode, with no positive quad points
# This is a test for accuracy.
set.seed(4378)
n <- 100
lambda <- 10
y <- rpois(n,lambda)

momobjfunc4 <- function(eta) {
  eta <- -1*eta
  sum(y) * eta - (length(y) + 1) * exp(eta) - sum(lgamma(y+1)) + eta
}
momfunlist4 <- list(
  fn = momobjfunc4,
  gr = function(x) numDeriv::grad(momobjfunc4,x),
  he = function(x) numDeriv::hessian(momobjfunc4,x)
)
momshiftquad4 <- aghq(momfunlist4,7,0)
truemoment4 <- -1*(digamma(sum(y) + 1) - log(length(y) + 1))
truesecondcentralmoment4 <- trigamma(sum(y)+1) # Still positive
truesecondrawmoment4 <- truesecondcentralmoment4 + truemoment4^2

# Moments
nummom_aghq_correct1_4 <- compute_moment(momshiftquad4,nn=1,method = 'correct')
nummom_aghq_correct2_4 <- compute_moment(momshiftquad4,nn=2,method = 'correct')
nummom_aghq_correct_central1_4 <- compute_moment(momshiftquad4,nn=1,method = 'correct',type='central')
nummom_aghq_correct_central2_4 <- compute_moment(momshiftquad4,nn=2,method = 'correct',type='central')

## Indexing of moments, 2d example. ##
logfteta2d <- function(eta,y) {
  n <- length(y)
  n1 <- ceiling(n/2)
  n2 <- floor(n/2)
  y1 <- y[1:n1]
  y2 <- y[(n1+1):(n1+n2)]
  eta1 <- eta[1]
  eta2 <- eta[2]
  sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
    sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
}
set.seed(57380)
n11 <- 10
n22 <- 10
lambda11 <- 5
lambda22 <- 5
nn <- n11+n22
y11 <- rpois(n11,5)
y22 <- rpois(n22,5)

truemoment_2d <- digamma(c(sum(y11) + 1,sum(y22) + 1)) - log(c(n11 + 1,n22+1))
truesecondcentralmoment_2d <- trigamma(c(sum(y11)+1,sum(y22) + 1)) # Still positive
truesecondrawmoment_2d <- truesecondcentralmoment_2d + truemoment_2d^2

momobjfunc_2d <- function(x) logfteta2d(x,c(y11,y22))
momfunlist_2d <- list(
  fn = momobjfunc_2d,
  gr = function(x) numDeriv::grad(momobjfunc_2d,x),
  he = function(x) numDeriv::hessian(momobjfunc_2d,x)
)
momquad_2d <- aghq(momfunlist_2d,7,c(0,0))

nummom_aghq_correct1_2d <- compute_moment(momquad_2d,nn=1,method = 'correct')
nummom_aghq_correct2_2d <- compute_moment(momquad_2d,nn=2,method = 'correct')
nummom_aghq_correct_central1_2d <- compute_moment(momquad_2d,nn=1,method = 'correct',type='central')
nummom_aghq_correct_central2_2d <- compute_moment(momquad_2d,nn=2,method = 'correct',type='central')

## TODO: add "correct" option into aghq and control ##
set.seed(84343124)
y_correct <- rpois(10,5) # Mode should be sum(y) / (10 + 1)
objfunc_correct <- function(x) logfteta2d(x,y_correct)
funlist_correct <- list(
  fn = objfunc_correct,
  gr = function(x) numDeriv::grad(objfunc_correct,x),
  he = function(x) numDeriv::hessian(objfunc_correct,x)
)
thequadrature_reuse <- aghq(funlist_correct,3,c(0,0),control = default_control(method_summaries='reuse'))
thequadrature_correct <- aghq(funlist_correct,3,c(0,0),control = default_control(method_summaries='correct'))
thesummary_reuse <- summary(thequadrature_reuse)
thesummary_correct <- summary(thequadrature_correct)
truemean_correct <- digamma(c(sum(y_correct[1:5]) + 1,sum(y_correct[6:10]) + 1)) - log(c(5 + 1,5+1))
truesd_correct <- sqrt(trigamma(c(sum(y_correct[1:5])+1,sum(y_correct[6:10]) + 1)))

## Nested ##
quadtable_p1_k1_prod <- get_quadtable(1,1)
quadtable_p2_k1_prod <- get_quadtable(2,1)
quadtable_p3_k1_prod <- get_quadtable(3,1)
quadtable_p4_k1_prod <- get_quadtable(4,1)
quadtable_p5_k1_prod <- get_quadtable(5,1)

quadtable_p1_k3_prod <- get_quadtable(1,3)
quadtable_p2_k3_prod <- get_quadtable(2,3)
quadtable_p3_k3_prod <- get_quadtable(3,3)
quadtable_p4_k3_prod <- get_quadtable(4,3)
quadtable_p5_k3_prod <- get_quadtable(5,3)

quadtable_p1_k5_prod <- get_quadtable(1,5)
quadtable_p2_k5_prod <- get_quadtable(2,5)
quadtable_p3_k5_prod <- get_quadtable(3,5)
quadtable_p4_k5_prod <- get_quadtable(4,5)
quadtable_p5_k5_prod <- get_quadtable(5,5)

quadtable_p1_k1_sparse <- get_quadtable(1,1,'sparse')
quadtable_p2_k1_sparse <- get_quadtable(2,1,'sparse')
quadtable_p3_k1_sparse <- get_quadtable(3,1,'sparse')
quadtable_p4_k1_sparse <- get_quadtable(4,1,'sparse')
quadtable_p5_k1_sparse <- get_quadtable(5,1,'sparse')

quadtable_p1_k3_sparse <- get_quadtable(1,3,'sparse')
quadtable_p2_k3_sparse <- get_quadtable(2,3,'sparse')
quadtable_p3_k3_sparse <- get_quadtable(3,3,'sparse')
quadtable_p4_k3_sparse <- get_quadtable(4,3,'sparse')
quadtable_p5_k3_sparse <- get_quadtable(5,3,'sparse')

quadtable_p1_k5_sparse <- get_quadtable(1,5,'sparse')
quadtable_p2_k5_sparse <- get_quadtable(2,5,'sparse')
quadtable_p3_k5_sparse <- get_quadtable(3,5,'sparse')
quadtable_p4_k5_sparse <- get_quadtable(4,5,'sparse')
quadtable_p5_k5_sparse <- get_quadtable(5,5,'sparse')

nestedoptlist1 <- list(
  fn = function(x) sum(dnorm(x,log=TRUE)),
  mode = rep(0,1)
)
nestedoptlist2 <- list(
  fn = function(x) sum(dnorm(x,log=TRUE)),
  mode = rep(0,2)
)
nestedoptlist3 <- list(
  fn = function(x) sum(dnorm(x,log=TRUE)),
  mode = rep(0,3)
)
nestedoptlist4 <- list(
  fn = function(x) sum(dnorm(x,log=TRUE)),
  mode = rep(0,4)
)
nestedoptlist5 <- list(
  fn = function(x) sum(dnorm(x,log=TRUE)),
  mode = rep(0,5)
)

nq_p1_k1_prod <- nested_quadrature(nestedoptlist1,1,'product')
nq_p2_k1_prod <- nested_quadrature(nestedoptlist2,1,'product')
nq_p3_k1_prod <- nested_quadrature(nestedoptlist3,1,'product')
nq_p4_k1_prod <- nested_quadrature(nestedoptlist4,1,'product')
nq_p5_k1_prod <- nested_quadrature(nestedoptlist5,1,'product')

nq_p1_k3_prod <- nested_quadrature(nestedoptlist1,3,'product')
nq_p2_k3_prod <- nested_quadrature(nestedoptlist2,3,'product')
nq_p3_k3_prod <- nested_quadrature(nestedoptlist3,3,'product')
nq_p4_k3_prod <- nested_quadrature(nestedoptlist4,3,'product')
nq_p5_k3_prod <- nested_quadrature(nestedoptlist5,3,'product')

nq_p1_k5_prod <- nested_quadrature(nestedoptlist1,5,'product')
nq_p2_k5_prod <- nested_quadrature(nestedoptlist2,5,'product')
nq_p3_k5_prod <- nested_quadrature(nestedoptlist3,5,'product')
nq_p4_k5_prod <- nested_quadrature(nestedoptlist4,5,'product')
nq_p5_k5_prod <- nested_quadrature(nestedoptlist5,5,'product')

nq_p1_k1_sparse <- nested_quadrature(nestedoptlist1,1,'sparse')
nq_p2_k1_sparse <- nested_quadrature(nestedoptlist2,1,'sparse')
nq_p3_k1_sparse <- nested_quadrature(nestedoptlist3,1,'sparse')
nq_p4_k1_sparse <- nested_quadrature(nestedoptlist4,1,'sparse')
nq_p5_k1_sparse <- nested_quadrature(nestedoptlist5,1,'sparse')

nq_p1_k3_sparse <- nested_quadrature(nestedoptlist1,3,'sparse')
nq_p2_k3_sparse <- nested_quadrature(nestedoptlist2,3,'sparse')
nq_p3_k3_sparse <- nested_quadrature(nestedoptlist3,3,'sparse')
nq_p4_k3_sparse <- nested_quadrature(nestedoptlist4,3,'sparse')
nq_p5_k3_sparse <- nested_quadrature(nestedoptlist5,3,'sparse')

nq_p1_k5_sparse <- nested_quadrature(nestedoptlist1,5,'sparse')
nq_p2_k5_sparse <- nested_quadrature(nestedoptlist2,5,'sparse')
nq_p3_k5_sparse <- nested_quadrature(nestedoptlist3,5,'sparse')
nq_p4_k5_sparse <- nested_quadrature(nestedoptlist4,5,'sparse')
nq_p5_k5_sparse <- nested_quadrature(nestedoptlist5,5,'sparse')

# Adaptive nested quadrature
man1 <- 1
man2 <- 1:2
man3 <- 1:3
man4 <- 1:4
man5 <- 1:5

Han1 <- as.matrix(2)
Han2 <- Matrix::crossprod(matrix(rnorm(2^2),nrow=2))
Han3 <- Matrix::crossprod(matrix(rnorm(3^2),nrow=3))
Han4 <- Matrix::crossprod(matrix(rnorm(4^2),nrow=4))
Han5 <- Matrix::crossprod(matrix(rnorm(5^2),nrow=5))


dmnorm <- function(x,man,Han) { # multivariate normal
  p <- length(man)
  if (p == 1) return(dnorm(x,man,1/sqrt(Han),log = TRUE))
  up <- -0.5 * as.numeric(crossprod(x-man,crossprod(Han,x-man)))
  down <- -0.5*p*log(2*pi) + 0.5*as.numeric(determinant(Han,logarithm = TRUE)$modulus)
  up + down
}
adaptivenestedoptlist1 <- list(
  fn = function(x) dmnorm(x,man1,Han1),
  mode = man1,
  hessian = as.matrix(Han1)
)
adaptivenestedoptlist2 <- list(
  fn = function(x) dmnorm(x,man2,Han2),
  mode = man2,
  hessian = Han2
)
adaptivenestedoptlist3 <- list(
  fn = function(x) dmnorm(x,man3,Han3),
  mode = man3,
  hessian = Han3
)
adaptivenestedoptlist4 <- list(
  fn = function(x) dmnorm(x,man4,Han4),
  mode = man4,
  hessian = Han4
)
adaptivenestedoptlist5 <- list(
  fn = function(x) dmnorm(x,man5,Han5),
  mode = man5,
  hessian = Han5
)



anq_p1_k1_prod <- adaptive_nested_quadrature(adaptivenestedoptlist1,1,'product')
anq_p2_k1_prod <- adaptive_nested_quadrature(adaptivenestedoptlist2,1,'product')
anq_p3_k1_prod <- adaptive_nested_quadrature(adaptivenestedoptlist3,1,'product')
anq_p4_k1_prod <- adaptive_nested_quadrature(adaptivenestedoptlist4,1,'product')
anq_p5_k1_prod <- adaptive_nested_quadrature(adaptivenestedoptlist5,1,'product')

anq_p1_k3_prod <- adaptive_nested_quadrature(adaptivenestedoptlist1,3,'product')
anq_p2_k3_prod <- adaptive_nested_quadrature(adaptivenestedoptlist2,3,'product')
anq_p3_k3_prod <- adaptive_nested_quadrature(adaptivenestedoptlist3,3,'product')
anq_p4_k3_prod <- adaptive_nested_quadrature(adaptivenestedoptlist4,3,'product')
anq_p5_k3_prod <- adaptive_nested_quadrature(adaptivenestedoptlist5,3,'product')

anq_p1_k5_prod <- adaptive_nested_quadrature(adaptivenestedoptlist1,5,'product')
anq_p2_k5_prod <- adaptive_nested_quadrature(adaptivenestedoptlist2,5,'product')
anq_p3_k5_prod <- adaptive_nested_quadrature(adaptivenestedoptlist3,5,'product')
anq_p4_k5_prod <- adaptive_nested_quadrature(adaptivenestedoptlist4,5,'product')
anq_p5_k5_prod <- adaptive_nested_quadrature(adaptivenestedoptlist5,5,'product')

anq_p1_k1_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist1,1,'sparse')
anq_p2_k1_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist2,1,'sparse')
anq_p3_k1_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist3,1,'sparse')
anq_p4_k1_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist4,1,'sparse')
anq_p5_k1_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist5,1,'sparse')

anq_p1_k3_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist1,3,'sparse')
anq_p2_k3_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist2,3,'sparse')
anq_p3_k3_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist3,3,'sparse')
anq_p4_k3_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist4,3,'sparse')
anq_p5_k3_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist5,3,'sparse')

anq_p1_k5_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist1,5,'sparse')
anq_p2_k5_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist2,5,'sparse')
anq_p3_k5_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist3,5,'sparse')
anq_p4_k5_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist4,5,'sparse')
anq_p5_k5_sparse <- adaptive_nested_quadrature(adaptivenestedoptlist5,5,'sparse')





