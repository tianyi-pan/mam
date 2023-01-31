### Misc Helper Functions ###

#' Create an empty sparse matrix
#'
#' Wrapper for new('dgCMatrix',...). Creates an empty sparse \code{dgCMatrix}
#' of specified dimension.
#'
#' @param n Number of rows
#' @param m Number of columns
#'
#' @return An object of class \code{dgCMatrix} having the appropriate structure,
#' and no data in teh \code{x} slot.
#' @examples
#' newsparsemat(10,5)
#'
#' @import Matrix
#' @export
#'
newsparsemat <- function(n,m) {
  methods::new('dgCMatrix',Dim=as.integer(c(n,m)),p=rep(0L,m+1),i=integer(0),x=numeric(0))
}

#' Flatten a block-diagonal sparse matrix
#'
#' Take a n x d block diagonal matrix with blocks of dimension m x p
#' and flatten to a (dense) matrix of dimension n x p
#'
#' @param Z Block diagonal matrix with equal-sized blocks (don't have to be square)
#' @param m Number of rows in each block
#' @param p Number of columns in each block
#'
#' @return A dense matrix with the appropriate structure
#'
#' @export
#'
flatten_bdmat <- function(Z,n,d) {
  # Z: sparse block diagonal matrix to flatten
  # n: block row dimension
  # d: block column dimension
  p <- ncol(Z)/d
  ZT <- methods::as(Z,'TsparseMatrix')
  ZT@j <- as.integer(rep(rep(0:(d-1),each=n),times=p))
  as.matrix(ZT[ ,1:d])
}

#' @name logit
#' @aliases logit
#' @aliases ilogit
#'
#' @title Logit and ilogit
#'
#' @description Quick implementations of the logit and inverse logit functions
#'
#' @param x Numeric vector
#'
#' @rdname logit
#' @returns \code{logit}: the logit of \code{x}, \code{log(x/(1-x))}
#' @examples
#' logit(seq(0,1,by=.1))
#' @export
#'
logit <- function(x) log(x/(1-x))
#' @rdname logit
#' @returns \code{ilogit}: the inverse logit of \code{x}, \code{1/(1+exp(-x))}
#' @examples
#' ilogit(seq(-5,5,by=1))
#' @export
#'
ilogit <- function(x) 1 / (1 + exp(-x))

# Not exported: splice
splice <- function(v,t,j) {
  # Insert t into v such that if vnew = splice(v,t,j) then vnew[j] == t
  if (j == 1) return(c(t,v))
  n <- length(v)
  if (j == (n+1)) return(c(v,t))
  if (j<=0) stop("j must be >0")
  if (j>(n+1)) stop("j must be <= n+1")
  c(v[1:(j-1)],t,v[j:n])
}

# Not exported: anyInf
anyInf <- function(x){
  sum(is.infinite(x))>0
}
