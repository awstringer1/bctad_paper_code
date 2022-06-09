#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA
  DATA_VECTOR(y); // Response
  DATA_IVECTOR(id); // id(i) = n(i), the number of observations in group i
  DATA_SPARSE_MATRIX(Z); // random effect design matrix
  DATA_VECTOR(nn); // GLQ nodes
  DATA_VECTOR(ww); // GLQ weights
  int k = nn.size();
  // PARAMETERS
  PARAMETER_VECTOR(beta); // intercepts
  PARAMETER_VECTOR(theta); // log(alpha[1...2])
  PARAMETER_MATRIX(u); // random intercepts, 2 columns
  
  // CONSTANTS
  int d = u.rows();
  int n = y.size();
  int s = theta.size();
  vector<Type> alpha(s);
  for (int i=0;i<s;i++) alpha(i) = exp(theta(i));
  
  // NEGATIVE LOG LIKELIHOOD
  vector<Type> eta1 = Z*u.col(0);
  vector<Type> eta2 = Z*u.col(1);
  vector<Type> eta3(n);
  vector<Type> eta4(n);
  
  for (int i=0;i<n;i++) {
    eta1(i) += beta(0);
    eta2(i) += beta(1);
    eta3(i) = beta(2);
    eta4(i) = beta(3);
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
      z = 0;
      for (int b=0;b<k;b++) z += pow((y(l)/mu - 1)*nn(b) + 1,nu-1)*ww(b);
      z = z * (y(l)/mu - 1) / sigma;
      loglik(i) -= dt(z,tau,true);
      loglik(i) -= (nu-1.0)*log(y(l)) - nu*log(mu) - log(sigma);
      l++;
    }
  }
  
  // RANDOM EFFECTS
  vector<Type> ranef(d);
  for (int i=0;i<d;i++) {
    ranef(i) = 0;
    ranef(i) -= dt(u(i,0), alpha(0), true);
    ranef(i) -= dnorm(u(i,1), Type(0), alpha(1), true);
  }
  
  vector<Type> jnll = loglik + ranef;
  return jnll.sum();
}
