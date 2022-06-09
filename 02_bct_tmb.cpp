#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA
  DATA_VECTOR(y); // Response
  DATA_IVECTOR(id); // id(i) = n(i), the number of observations in group i
  DATA_SPARSE_MATRIX(Z); // random effect design matrix
  DATA_INTEGER(approx_type); // 0 = exact (but doesn't work at nu = 0), 1 = series, 2 = GLQ
  // Input the data for both approximations, easier to code
  DATA_INTEGER(m); // Number of series terms
  DATA_VECTOR(nn); // GLQ nodes
  DATA_VECTOR(ww); // GLQ weights
  int k = nn.size();
  // PARAMETERS
  PARAMETER_VECTOR(beta); // intercepts
  PARAMETER_VECTOR(theta); // log(alpha[1...4]); if re_dist = 2, then also log(df[1...4])
  PARAMETER_MATRIX(u); // random intercepts, 4 columns
  
  // CONSTANTS
  int d = u.rows();
  int n = y.size();
  int s = theta.size();
  vector<Type> alpha(s);
  for (int i=0;i<s;i++) alpha(i) = exp(theta(i));
  
  // NEGATIVE LOG LIKELIHOOD
  vector<Type> eta1 = Z*u.col(0);
  vector<Type> eta2 = Z*u.col(1);
  vector<Type> eta3 = Z*u.col(2);
  vector<Type> eta4 = Z*u.col(3);
  for (int i=0;i<n;i++) {
    eta1(i) += beta(0);
    eta2(i) += beta(1);
    eta3(i) += beta(2);
    eta4(i) += beta(3);
  }
  
  vector<Type> loglik(d);
  int l = 0; // Running index
  Type mu=0;
  Type sigma=0;
  Type nu=0;
  Type tau=0;
  Type z=0;
  for (int i=0;i<d;i++) {
    loglik(i) = 0;
    for (int j=0;j<id(i);j++) {
      mu = exp(eta1(l));
      sigma = exp(eta2(l));
      nu = eta3(l);
      tau = exp(eta4(l));
      // Approximate z
      if (approx_type == 0) {
        z = ( pow(y(l)/mu,nu) -1.0 ) / (nu*sigma); // Exact, but doesn't work at nu = 0
      } else if (approx_type == 1) {
        z = log(y(l)/mu);
        for (int b=2;b<=m;b++) z += pow(nu,b-1) * pow(log(y(l)/mu),b) * exp(-1.0*lfactorial(Type(b)));
        z = z/sigma;
      } else if (approx_type == 2) {
        z = 0;
        for (int b=0;b<k;b++) z += pow((y(l)/mu - 1)*nn(b) + 1,nu-1)*ww(b);
        z = z * (y(l)/mu - 1) / sigma;
      }
      loglik(i) -= dt(z,tau,true);
      loglik(i) -= (nu-1.0)*log(y(l)) - nu*log(mu) - log(sigma);
      l++;
    }
  }
  
  // RANDOM EFFECTS
  vector<Type> ranef(d);
  for (int i=0;i<d;i++) {
    ranef(i) = 0;
    ranef(i) -= dnorm(u(i,0), Type(0), alpha(0), true);
    ranef(i) -= dnorm(u(i,1), Type(0), alpha(1), true);
    ranef(i) -= dnorm(u(i,2), Type(0), alpha(2), true);
    ranef(i) -= dnorm(u(i,3), Type(0), alpha(3), true);
  }
  
  REPORT(loglik);
  REPORT(ranef);
  REPORT(u);
  REPORT(beta);
  REPORT(alpha);
  REPORT(eta1);
  REPORT(eta2);
  REPORT(eta3);
  REPORT(eta4);
  
  vector<Type> jnll = loglik + ranef;
  REPORT(jnll);
  
  return jnll.sum();
}
