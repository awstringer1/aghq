#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>
  
template<class Type>
Type objective_function<Type>::operator() ()
{
// Taken from glmmTMB: check if openmp configured, set parallelization accordingly.
#ifdef _OPENMP
  this -> max_parallel_regions = omp_get_max_threads();
  // std::cout << "OpenMP max_parallel_regions=" << this -> max_parallel_regions << "\n";
#else
  this -> max_parallel_regions = 1;
  // std::cout << "no OpenMP (max_parallel_regions=1)\n";
#endif

  // Bernoulli GLMM with nested random effects
  // Read in data
  DATA_VECTOR(y); // Response
  DATA_SPARSE_MATRIX(X); // Design matrix, g(E(Y)) = Xbeta + ZU
  DATA_SPARSE_MATRIX(Z); // Design matrix, g(E(Y)) = Xbeta + ZU
  
  // Parameters
  PARAMETER_VECTOR(beta); // Regression coefficients
  PARAMETER_VECTOR(theta); // sigma = SD(random intercept). theta = -2log(sigma), sigma = exp(-theta/2), sigma^2 = exp(-theta)
  int s = theta.size(); // Should be 2
  vector<Type> sigmasq(s);
  for (int i=0;i<s;i++) sigmasq(i) = exp(-theta(i));
  PARAMETER_VECTOR(Us); // Random intercepts for states
  PARAMETER_VECTOR(Ut); // Random intercepts for towns
  int m1 = Us.size();
  int m2 = Ut.size();
  double m1d;
  double m2d;
  m1d = m1;
  m2d = m2;
  int m = m1+m2;
  vector<Type> U(m);
  U << Us,Ut;
  REPORT(sigmasq);

  int n = y.size();
  vector<Type> eta = X*beta + Z*U;

  // log likelihood
  Type loglik = 0;
  for (int i=0;i<n;i++) PARALLEL_REGION {loglik -= y(i)*eta(i) - log(1 + exp(eta(i)));}
  // log priors 
  double pi = 3.141592653589793115998;
  // u
  PARALLEL_REGION loglik -= -0.5*log(2.0*pi*sigmasq(0))*m1d - (Us*Us).sum() / (2.0*sigmasq(0));
  PARALLEL_REGION loglik -= -0.5*log(2.0*pi*sigmasq(1))*m2d - (Ut*Ut).sum() / (2.0*sigmasq(1));
  // beta
  double pb = beta.size();
  PARALLEL_REGION loglik -= -0.5*log(2.0*pi*1000.0)*pb - (beta*beta).sum() / (2.0*1000.0);
  // log sigma
  Type lambda = -log(0.5)/0.5;
  for (int i=0;i<s;i++) PARALLEL_REGION {loglik -= log(lambda*0.5) - lambda * exp(-theta(i)*0.5) - theta(i)*0.5;}

  return loglik;

}
