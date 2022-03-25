### Control Options ###

#' Control Options for Marginal Additive Models
#'
#' Control options for the function \code{mam::mam()}.
#'
#' @param ... Used to provide all arguments. See details for arguments.
#'
#' @details Running with an empty argument list returns defaults. Setting any of the following
#' arguments overwrites the default:
#' \itemize{
#' \item{k}{Default 5, number of Gauss-Hermite Quadrature points to use. For \code{p}-dimensional
#' random effects, a product-rule extension is used, requiring \code{k^p} total points
#' per integral evaluation.}
#' \item{retcond}{Default FALSE. Retain estimates of the effects from the underlying conditional model?
#' Variance components are retained regardless, but turning this off doesn't compute uncertainty estimates
#' for the conditional model, at a potentially large computational benefit.}
#' \item{method}{Default \code{'BFGS'}, which optimization procedure to use? \code{'BFGS'} does
#' \code{stats::optim(...,method='BFGS')}, and \code{'trust'} does \code{trustOptim::trust.optim}. Will be
#' set to \code{'BFGS'} automatically, with a warning, if package \code{trustOptim} not installed.}
#' \item{varmethod}{Default \code{1}. Method to use to calculate the pointwise variances of the fitted curves. The
#' recommendation is \code{varmethod=1} which implements the "Kass and Steffey" method described in ..., Section ...
#' and Algorithm ...\code{varmethod = 2} uses the inverse of the negative Hessian of the joint penalized log-likelihood,
#' which may not be positive definite (in fact, may well be *indefinite*, having one or more negative eigenvalues) at
#' the point estimate used. Setting \code{varmethod=0} computes the \code{LDL} decomposition of this Hessian, and then sets
#' \code{varmethod=2} if all the \code{D} values are positive, and \code{varmethod=1} otherwise. We recommend just
#' using \code{varmethod=1}.}
#' \item{verbose}{Default \code{FALSE}. Print progress and diagnostic information as the function proceeds with
#' computations?}
#' }
#'
#' @examples
#' mam_control() # Return default settings
#' mam_control(k = 9) # Default settings, except changing \code{k} to \code{9}.
#'
#' @rawNamespace useDynLib(mam, .registration=TRUE); useDynLib(mam_TMBExports)
#'
#' @export
#'
mam_control <- function(...) {
  out <- list(
    k=5,
    retcond=FALSE,
    method = c("BFGS","trust"),
    varmethod = c(1,0,2),
    verbose = FALSE
  )
  specialargs <- list(...)
  for (arg in names(specialargs)) out[arg] <- specialargs[arg]
  out
}
