#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Read in data
  // 
  // Distance matrix
  DATA_MATRIX(D);
  // Infection and removal times
  // Ordered by Infection time
  // Distance matrix appropriately sorted
  DATA_VECTOR(I);
  DATA_VECTOR(R);
  DATA_VECTOR(infected); // Indictor of infections

  int N = I.size();
  int m = infected.size();
  
  // Parameters
  // theta1 = log(alpha)
  // theta2 = log(beta)
  PARAMETER(theta1);
  PARAMETER(theta2);
  
  Type alpha = exp(theta1);
  Type beta = exp(theta2);
  
  // Infectivity
  matrix<Type> lambda(D);
  for (int i = 0;i<N;i++) {
    for (int j = 0;j<N;j++) {
      if (i == j) {
        lambda(i,j) = 0;
      } else {
        lambda(i,j) = alpha * pow(D(i,j),-beta);
      }
    }
  }
  
  Type t1 = 0;
  Type t1inner = 0;
  Type t2 = 0;
  
  // First term
  for (int j = 1;j<m;j++) {
    t1inner = 0;
    for (int i = 0;i<N;i++) {
      if (I(i) < I(j) & I(j) <= R(i)) {
        t1inner += lambda(i,j);
      } 
    }
    t1 += log(t1inner);
  }
  
  // Second term
  for (int i = 0;i<m;i++) {
    for (int j = 0;j<N;j++) {
      Type Rmin = R(i);
      if (R(i) > I(j)) Rmin = I(j);
      Type Imin = I(i);
      if (I(i) > I(j)) Imin = I(j);
      t2 += (Rmin - Imin) * lambda(i,j);
    }
  }
  
  Type logpost = t1 - t2;
  // Add the log-priors
  Type phi = 0.01;
  logpost += log(phi) + log(alpha) - phi*alpha;
  logpost += log(phi) + log(beta) - phi*beta;
  
  // Return NEGATED logpost, to facilitate MCMC sampling using tmbstan
  return -1.0 * logpost;
}
