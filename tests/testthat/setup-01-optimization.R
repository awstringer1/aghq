library(aghq)
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

