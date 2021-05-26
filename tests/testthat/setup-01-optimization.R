# devtools::load_all()
options(mc.cores = 2)
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

opt_sparsetrust <- optimize_theta(funlist,1.5)
opt_sr1 <- optimize_theta(funlist,1.5,control = list(method = "SR1"))
opt_trust <- optimize_theta(funlist,1.5,control = list(method = "trust"))
opt_bfgs <- optimize_theta(funlist,1.5,control = list(method = "BFGS"))

norm_sparse_3 <- normalize_logpost(opt_sparsetrust,3,1)
norm_sparse_5 <- normalize_logpost(opt_sparsetrust,5,1)
norm_sparse_7 <- normalize_logpost(opt_sparsetrust,7,1)

norm_trust_3 <- normalize_logpost(opt_trust,3,1)
norm_trust_5 <- normalize_logpost(opt_trust,5,1)
norm_trust_7 <- normalize_logpost(opt_trust,7,1)

norm_bfgs_3 <- normalize_logpost(opt_bfgs,3,1)
norm_bfgs_5 <- normalize_logpost(opt_bfgs,5,1)
norm_bfgs_7 <- normalize_logpost(opt_bfgs,7,1)

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

opt_sparsetrust_2d <- optimize_theta(funlist2d,c(1.5,1.5))
opt_trust_2d <- optimize_theta(funlist2d,c(1.5,1.5),control = list(method = "trust"))
opt_sr1_2d <- optimize_theta(funlist2d,c(1.5,1.5),control = list(method = "SR1"))
opt_bfgs_2d <- optimize_theta(funlist2d,c(1.5,1.5),control = list(method = "BFGS"))

norm_sparse_2d_3 <- normalize_logpost(opt_sparsetrust_2d,3,1)
norm_sparse_2d_5 <- normalize_logpost(opt_sparsetrust_2d,5,1)
norm_sparse_2d_7 <- normalize_logpost(opt_sparsetrust_2d,7,1)

norm_trust_2d_3 <- normalize_logpost(opt_trust_2d,3,1)
norm_trust_2d_5 <- normalize_logpost(opt_trust_2d,5,1)
norm_trust_2d_7 <- normalize_logpost(opt_trust_2d,7,1)

norm_bfgs_2d_3 <- normalize_logpost(opt_bfgs_2d,3,1)
norm_bfgs_2d_5 <- normalize_logpost(opt_bfgs_2d,5,1)
norm_bfgs_2d_7 <- normalize_logpost(opt_bfgs_2d,7,1)

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
pdf_poly_2d <- compute_pdf_and_cdf(thequadrature_k7,transformation = list(totheta = log,fromtheta = exp))
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


opt_sparsetrust_3d <- optimize_theta(funlist3d,c(1.5,1.5,1.5))
opt_trust_3d <- optimize_theta(funlist3d,c(1.5,1.5,1.5),control = list(method = "trust"))
opt_sr1_3d <- optimize_theta(funlist3d,c(1.5,1.5,1.5),control = list(method = "SR1"))
opt_bfgs_3d <- optimize_theta(funlist3d,c(1.5,1.5,1.5),control = list(method = "BFGS"))

norm_sparse_3d_3 <- normalize_logpost(opt_sparsetrust_3d,3,1)
norm_sparse_3d_5 <- normalize_logpost(opt_sparsetrust_3d,5,1)
norm_sparse_3d_7 <- normalize_logpost(opt_sparsetrust_3d,7,1)

norm_trust_3d_3 <- normalize_logpost(opt_trust_3d,3,1)
norm_trust_3d_5 <- normalize_logpost(opt_trust_3d,5,1)
norm_trust_3d_7 <- normalize_logpost(opt_trust_3d,7,1)

norm_bfgs_3d_3 <- normalize_logpost(opt_bfgs_3d,3,1)
norm_bfgs_3d_5 <- normalize_logpost(opt_bfgs_3d,5,1)
norm_bfgs_3d_7 <- normalize_logpost(opt_bfgs_3d,7,1)

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
# sparsegrid_2d <- aghq(funlist2d,5,c(0,0),control = default_control(ndConstruction = 'sparse'))
# sparsenormconst_2d <- compute_moment(sparsegrid_2d) # Should be 1

## Control Params ##
cntrl_base <- default_control()
cntrl_marg <- default_control_marglaplace()
cntrl_tmb <- default_control_tmb()




