#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

#include <Eigen/Sparse>

template<class Type>
Type marginal_mean2(objective_function<Type>* obj){
  using namespace Eigen;

  // DATA
  // DATA_INTEGER(link); // See "valid_link" above for coding
  // DATA_VECTOR(xf);
  // DATA_VECTOR(xr);
  // DATA_VECTOR(v);
  DATA_MATRIX(XF);
  DATA_MATRIX(XR);
  DATA_MATRIX(Z); // NOTE: collapse Z so it's a d-column matrix where d = dim(U)
  int n = XF.rows();
  int d = Z.cols();

  DATA_MATRIX(Q); // GHQ nodes, one node per row
  DATA_VECTOR(w); // GHQ weights, length(w) = nrow(Q)
  DATA_SPARSE_MATRIX(Lam); // cov(U) = Lam%*%Lam^T, sparse covariance factor, dim(Lam) = ncol(Q) = ncol(Z) = d
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

  int kd = Q.rows();
  // TRANSFORMATIONS
  // Type xTb = (xf*betaF).sum() + (xr*bR).sum();
  vector<Type> Xtb = XF*betaF + XR*bR;
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
  // REPORT(Lam);

  // Compute output
  // double pi = 3.141592653589793115998;
  vector<Type> out(n), tmp(d), Qi(d), v(d);
  Type vTz = 0,xTb = 0,tmpout;

  for (int i=0;i<n;i++) {
    xTb = Xtb(i);
    v = Z.row(i);
    out(i) = 0.;
    tmpout = 0.;
    for (int j=0;j<kd;j++) {
      Qi = Q.row(j);
      tmp = Lam * Qi;
      vTz = (v*tmp).sum();
      tmpout += invlogit(xTb+vTz) * w(j);
    }
    out(i) = log(tmpout / (1-tmpout));
    // out(i) = tmpout;
  }
  ADREPORT(out);

  // for (int i=0;i<Q.rows();i++) {
  //   // tmp = 0;
  //   Qi = Q.row(i);
  //   tmp = Lam * Qi;
  //   vTz = (Z.row(i)*tmp).sum();
  //   // ONLY IMPLEMENTED FOR LOGIT CURRENTLY
  //   // out += (1.0 / (1.0 + exp(-1.0*(xTb+vTz))))*w(i);
  //   out += invlogit(xtb+vTz) * w(i);
  //   // out+=1;
  // }
  // Return g(mu). ONLY IMPLEMENTED FOR LOGIT CURRENTLY
  // Type gpix = log(out / (1-out));
  // REPORT(out);
  return 0.;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

