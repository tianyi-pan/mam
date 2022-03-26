#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

#include <Eigen/Sparse>

template<class Type>
Type marginal_mean(objective_function<Type>* obj){
  using namespace Eigen;

  // DATA
  // DATA_INTEGER(link); // See "valid_link" above for coding
  DATA_VECTOR(xf);
  DATA_VECTOR(xr);
  DATA_VECTOR(v);
  DATA_MATRIX(Q); // GHQ nodes, one node per row
  DATA_VECTOR(w); // GHQ weights, length(w) = nrow(Q)
  DATA_SPARSE_MATRIX(Lam); // cov(U) = Lam%*%Lam^T, sparse covariance factor, but now, JUST ONE group's, i.e. dim(Lam) = dim(U_i)
  DATA_IVECTOR(Lind); // Lam@x[] = theta[Lind]
  DATA_IVECTOR(diagind); // Identifies which elements of theta correspond to diagonal elements of Lambda
  // PARAMETERS
  // Ordered according to how we want the precision matrix ordered
  PARAMETER_VECTOR(betaF);
  PARAMETER_VECTOR(bR);
  PARAMETER_VECTOR(theta); // Covariance factor parameters
  // The LA margmean doesn't actually depend on lambda or U_hat
  // Pass them in so the Jacobian is computed correctly. It should have all zeroes
  PARAMETER_VECTOR(loglambda);
  PARAMETER_VECTOR(U_hat);

  // TRANSFORMATIONS
  Type xTb = (xf*betaF).sum() + (xr*bR).sum();
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

  // Compute output
  // double pi = 3.141592653589793115998;
  Type out = 0;
  vector<Type> tmp(Lam.rows());
  vector<Type> Qi(Q.cols());
  Type vTz = 0;
  for (int i=0;i<Q.rows();i++) {
    // tmp = 0;
    Qi = Q.row(i);
    tmp = Lam * Qi;
    vTz = (v*tmp).sum();
    // ONLY IMPLEMENTED FOR LOGIT CURRENTLY
    // out += (1.0 / (1.0 + exp(-1.0*(xTb+vTz))))*w(i);
    out += invlogit(xTb+vTz) * w(i);
    // out+=1;
  }
  // Return g(mu). ONLY IMPLEMENTED FOR LOGIT CURRENTLY
  Type gpix = log(out / (1-out));
  REPORT(out);
  return gpix;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

