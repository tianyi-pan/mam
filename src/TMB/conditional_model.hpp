
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

#include <Eigen/Sparse>

template<class Type>
Type conditional_model(objective_function<Type>* obj){
  using namespace Eigen;
  using namespace tmbutils;
  using namespace density;

  // Bernoulli GLMM with random effects for smooth and iid terms
  // UPDATED to allow multiple smooth terms and random effects params
  // DATA
  DATA_VECTOR(y); // Response
  // eta = XF*betaF + XR*br + AU
  DATA_MATRIX(XF); // Unpenalized smooth term design matrix
  DATA_MATRIX(XR); // Penalized smooth term design matrix
  DATA_SPARSE_MATRIX(A); // Random intercept (etc) design matrix
  DATA_SPARSE_MATRIX(Lam); // cov(U) = Lam%*%Lam^T, sparse covariance factor
  DATA_IVECTOR(Lind); // Lam@x[] = theta[Lind]
  DATA_IVECTOR(diagind); // Identifies which elements of theta correspond to diagonal elements of Lambda

  // PARAMETERS
  PARAMETER_VECTOR(betaF); // Unpenalized smooth term regression coefficients
  PARAMETER_VECTOR(bR); // Penalized smooth term regression coefficients
  // PARAMETER_VECTOR(logresd); // Vector of log random effects standard deviations
  PARAMETER_VECTOR(theta); // vector of covariance factor parameters
  PARAMETER_VECTOR(logsmoothing); // Vector of log smoothing parameters

  PARAMETER_VECTOR(U) // Random effects

    // DIMENSIONS
    // Input the dimensions as data, integer types
    DATA_INTEGER(M); // Number of random effects vectors
  DATA_INTEGER(s); // Dimension of each random effects vector (same)
  int Udim = M*s; // Total number of random effects "parameters", i.e. length(U)
  DATA_INTEGER(p); // Number of smooth terms
  DATA_IVECTOR(r); // Vector of dimension of each smooth term. sum(r) = dim(bR), length(r) = p
  int bdim = r.sum();

  // TRANSFORMATIONS
  // Create the vector of smoothing params
  // vector<Type> reprec(logresd.size());
  vector<Type> smoothprec(logsmoothing.size());
  // for (int i=0;i<reprec.size();i++) reprec(i) = exp(-2.0*logresd(i));
  for (int i=0;i<smoothprec.size();i++) smoothprec(i) = exp(logsmoothing(i));
  // Exponentiate the elements of theta that correspond to diagonal elements of Lambda
  int ts = theta.size();
  vector<Type> covparam(ts);
  for (int i=0;i<ts;i++) {
    if (diagind(i) == 1) {
      covparam(i) = exp(theta(i));
    } else {
      covparam(i) = theta(i);
    }
  }
  // Create the random effects covariance factor, by updating the supplied Lam template.
  int tempitr=0;
  for (int k=0; k<Lam.outerSize(); ++k)
    for (typename SparseMatrix<Type,ColMajor>::InnerIterator it(Lam,k); it; ++it)
    {
      it.valueRef() = covparam(Lind(tempitr));
      tempitr++;
    }

    // Linear predictor
    vector<Type> eta = XF*betaF + XR*bR + A*U;
  int n = y.size();

  // NEGATIVE log likelihood
  Type loglik = 0;
  for (int i=0;i<n;i++) loglik -= y(i)*eta(i) - log(1 + exp(eta(i)));

  // NEGATIVE log priors
  double pi = 3.141592653589793115998;
  matrix<Type> Sigma = Lam * Lam.transpose(); // Kills sparsity; I don't see any other option though.
  Type normdens = MVNORM(Sigma)(U);
  REPORT(normdens);
  loglik += normdens;


  // bR
  // loglik -= -0.5*log(2.0*pi)*r +0.5*r*log(lambda) - lambda * (bR*bR).sum() / (2.0);
  double Rdub = double(bdim);
  vector<double> rdub(p);
  for (int j=0;j<p;j++) rdub(j) = double(r(j));
  loglik -= -0.5*Rdub*log(2.0*pi);
  int tmpidx2 = 0;
  for (int j=0;j<p;j++) {
    loglik -=0.5*rdub(j)*logsmoothing(j);
    for (int l=0;l<r(j);l++) {
      loglik -= -0.5*smoothprec(j) * bR(tmpidx2)*bR(tmpidx2);
      tmpidx2++;
    }
  }

  return loglik; // Actually minus loglik
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

