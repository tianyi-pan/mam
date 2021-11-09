### Marginal Additive Model ###

#' Estimate marginal coefficients from an additive model
#'
#' This function takes an object of class \code{gam} created using \code{mgcv::gam},
#' and returns a marginal additive model object.
#'
#' @param mod Object of class \code{gam} created using \code{mgcv::gam}.
#' @param ... not used.
#'
#' @returns Object of classes \code{c('mam','gam','glm','lm')} representing the marginal additive model.
#'
#' @export
#'
#' @examples
#' ## Not run
#'
mam <- function(mod,...) UseMethod('mam')
#' @rdname mam
#' @method mam default
#' @export
mam.default <- function(mod,...) 0
#' @rdname mam
#' @method mam gam
#' @export
mam.gam <- function(mod,...) 0

