#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW);
  
  DATA_MATRIX(y); //  Measurements
  PARAMETER_MATRIX(W); // Measurement errors
  // ADREPORT(W);
  DATA_MATRIX(B); // T x A_i, stacked rowwise; dimension (n x 3) x 3
  DATA_MATRIX(V); // The angle matrix, stacked rowwise; dimension (n x 3) x 3
  
  DATA_VECTOR(plx); // Parallax
  DATA_VECTOR(Xgc);
  DATA_VECTOR(Ygc);
  
  // Pass in the four nonlinear parameters
  // NOT transformed (transform them outside)
  // Pass as DATA. All the inference on the measurement errors is
  // conditional on the values of the nonlinear parameters.
  DATA_SCALAR(p);
  DATA_SCALAR(g);
  DATA_SCALAR(a);
  DATA_SCALAR(b);
  
  
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
    Rgc(i) = y(i,0) + W(i,0); // Rgc
    Vlos(i) = y(i,1) + W(i,1); // Vlos
    //PMra(i) = y(i,2) + W(i,2); // PMra
    //PMdec(i) = y(i,3) + W(i,3); // PMdec
    // For now, don't use measurement errors on the PM
    PMra(i) = y(i,2);
    PMdec(i) = y(i,3);
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
    
    //BBB = BB(i);
    //UVW = BBB * motion;
    
    
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
  REPORT(transdat);
  
  // Log-likelihood (NOT NEGATIVE)
  // Do NOT need the normalizing constant (with the loggammas)-- can compute that in R
  // Because we only need derivatives wrt W
  Type ll = 0;
  ll = 0;
  
  Type normconst;
  Type ppi;
  ppi = 3.141592653589793115998;
  
  // Energy
  Type Es;
  //Es = 0;
  Type pot;
  
  vector<Type> saveEs(n);
  //vector<Type> saveDatProd(n);
  //vector<Type> saveDat1(n);
  //vector<Type> saveF1(n);
  //vector<Type> saveF2(n);
  
  for (int i = 0;i<n;i++) {
    // transdat(i,0) == r
    pot = p * pow(transdat(i,0),-g);
    // transdat(i,1) == vr; transdat(i,2) == vt
    Es = pot - 0.5 * ( pow(transdat(i,1),2) + pow(transdat(i,2),2) );
    
    ll += -2*b*log( transdat(i,0) * transdat(i,2) ) + ( b*(g-2)/g + a/g - 1.5 )*log(Es);
    
    // Save quantities for constraints
    saveEs(i) = Es;
    //saveDatProd(i) = transdat(i,0) * transdat(i,2);
    //saveDat1(i) = transdat(i,0);
    //saveF1(i) = pow(saveDatProd(i),-2*b);
    //saveF2(i) = pow(Es,( b*(g-2)/g + a/g - 1.5 ));
  }
  ADREPORT(saveEs); // For nonlinear constraint
  //ADREPORT(saveDatProd);
  //ADREPORT(saveDat1);
  //ADREPORT(saveF1);
  //ADREPORT(saveF2);
  
  return ll;
}
