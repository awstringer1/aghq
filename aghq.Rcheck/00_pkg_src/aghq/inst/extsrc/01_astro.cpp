#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW);

  DATA_MATRIX(y); //  Measurements
  // ADREPORT(W);
  DATA_MATRIX(B); // T x A_i, stacked rowwise; dimension (n x 3) x 3
  DATA_MATRIX(V); // The angle matrix, stacked rowwise; dimension (n x 3) x 3

  DATA_VECTOR(plx); // Parallax
  DATA_VECTOR(Xgc);
  DATA_VECTOR(Ygc);

  // Pass in the four nonlinear parameters
  // TRANSFORMED, because now we need theta derivatives.
  PARAMETER(theta1); // Psi0
  PARAMETER(theta2); // gamma
  PARAMETER(theta3); // alpha
  PARAMETER(theta4); // beta

  // Transform params
  // Hardcode in the bounds, easier
  Type p = (200.0 - 1.0) * exp(-1.0 * exp(theta1)) + 1.0;
  Type g = (0.7 - 0.3) * exp(-1.0 * exp(theta2)) + 0.3;
  Type a = exp(theta3) + 3.0;
  Type b = (1.0 + 0.5) * exp(-1.0 * exp(theta4)) - 0.5; // beta bounded below by NEGATIVE 0.5


  // Constants
  Type k = 4.74057;

  vector<Type> sunmotion(3);
  sunmotion(0) = 11.1;
  sunmotion(1) = 12.24;
  sunmotion(2) = 7.25;

  vector<Type> vLSR(3);
  vLSR(0) = 0.0;
  vLSR(1) = 220.0;
  vLSR(2) = 0.0;


  // Input data, with measurement errors added
  // y is 4 columns: Rgc, Vlos, PMra, PMdec
  int n = y.col(0).size();


  // Add the measurement errors to the data
  vector<Type> Rgc(n);
  vector<Type> Vlos(n);
  vector<Type> PMra(n);
  vector<Type> PMdec(n);
  for (int i = 0;i<n;i++) {
    Rgc(i) = y(i,0); // Galactocentric position
    Vlos(i) = y(i,1); // Heliocentric line of sight velocity
    PMra(i) = y(i,2); // Proper motion in right ascension (corrected for declanation)
    PMdec(i) = y(i,3); // Proper motion in declanation
  }

  // Transform data to galactocentric position/velocity
  matrix<Type> transdat(n,3);
  vector<Type> motion(3);
  vector<Type> UVW(3);
  matrix<Type> BB(3,3);
  matrix<Type> VV(3,3);
  Type pR;
  vector<Type> PTZ(3);
  vector<Type> vrvtvp(3);
  for (int i = 0;i<n;i++) {
    // Step 1
    motion(0) = Vlos(i); // Vlos
    motion(1) = k * PMra(i) / plx(i); // PMra
    motion(2) = k * PMdec(i) / plx(i); // PMdec

    BB(0,0) = B(3*i,0); BB(0,1) = B(3*i,1); BB(0,2) = B(3*i,2);
    BB(1,0) = B(3*i+1,0); BB(1,1) = B(3*i+1,1); BB(1,2) = B(3*i+1,2);
    BB(2,0) = B(3*i+2,0); BB(2,1) = B(3*i+2,1); BB(2,2) = B(3*i+2,2);

    UVW = BB * motion;

    UVW = UVW + sunmotion + vLSR;

    // Step 2
    pR = sqrt( pow(Xgc(i),2) + pow(Ygc(i),2) );
    PTZ(0) = UVW(0) * Xgc(i) / pR + UVW(1) * Ygc(i) / pR;
    PTZ(1) = -1 * UVW(0) * Ygc(i) / pR + UVW(1) * Xgc(i) / pR;
    PTZ(2) = UVW(2);

    // Step 3

    VV(0,0) = V(3*i,0); VV(0,1) = V(3*i,1); VV(0,2) = V(3*i,2);
    VV(1,0) = V(3*i+1,0); VV(1,1) = V(3*i+1,1); VV(1,2) = V(3*i+1,2);
    VV(2,0) = V(3*i+2,0); VV(2,1) = V(3*i+2,1); VV(2,2) = V(3*i+2,2);

    vrvtvp = VV * PTZ / 100;

    // Output
    transdat(i,0) = Rgc(i); // r
    transdat(i,1) = vrvtvp(0); // vr
    transdat(i,2) = sqrt( pow(vrvtvp(1),2) + pow(vrvtvp(2),2) ); // vt
  }

  // Log-likelihood (NOT NEGATIVE)
  Type ll;
  ll = 0;

  Type ppi;
  ppi = 3.141592653589793115998;

  Type normconst;
  //normconst = 0;

  normconst = (2.0*b/g - a/g)*log(p) - log( ppi*sqrt(ppi) * pow(2.0,-b+1.5) ) - lgamma(1.0-b) + lgamma(a/g - 2.0*b/g + 1.0) - lgamma(b*(g-2.0)/g + a/g - 0.5 );
  REPORT(normconst);


  // Energy
  Type Es;
  Type pot;

  // Nonlinear constraints

  Type const1;
  const1 = a - g;
  Type const2;
  const2 = a - b * (2.0 - g) - g/2.0;
  ADREPORT(const1);
  ADREPORT(const2);


  vector<Type> saveEs(n);

  for (int i = 0;i<n;i++) {
    // transdat(i,0) == r
    pot = p * pow(transdat(i,0),-g);
    // transdat(i,1) == vr; transdat(i,2) == vt
    Es = pot - 0.5 * ( pow(transdat(i,1),2) + pow(transdat(i,2),2) );

    ll += normconst -2*b*log( transdat(i,0) * transdat(i,2) ) + ( b*(g-2)/g + a/g - 1.5 )*log(Es);

    // Save quantities for constraints
    saveEs(i) = Es;
  }
  ADREPORT(saveEs); // For nonlinear constraint

  // Add on the log-prior
  Type logprior =
    (-1.0 * exp(theta1) + theta1) +
    (-1.0 * exp(theta2) + theta2) +
    (-1.0 * exp(theta4) + theta4) +
    ( theta3 + log(4.6) - 4.6 * exp(theta3) );

  ll += logprior;

  // NEGATED log-posterior, for compatibility with other TMB features
  return -1.0 * ll;
}
