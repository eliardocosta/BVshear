
#' Brutsaert-Veron shear driven distribution
#' 
#' Density, distribution function, quantile function and random generation for the Brutsaert-Veron shear driven distribution.
#' 
#' @aliases BVshear
#'
#' @param x vector of quantiles. 
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations for random generation.
#' @param u friction velocity.
#' @param ... additional arguments to be passed to the Pearson Type VI distribution of the PearsonDS package.
#'
#' @return \code{dBVshear} gives the density, \code{pBVshear} gives the distribution function, \code{qBVshear} gives the quantile function, and \code{rBVshear} provides the random generation.
#' @export
#' 
#' @name BVshear
#' @rdname bvshear
#' 
#' @note These functions are just wrappers of the Pearson Type VI distribution functions of the package \code{PearsonDS} applied to the context of the surface renewal distribution. See reference below. 
#' 
#' @references Solheid, B.S., Dias, N.L.C., and Costa, E.G. (2023). The equations for the nondimensional moments of the surface renewal distribution.
#'
#' @examples 
#' # density at 0.03 for u = 0.5
#' dBVshear(x = 0.03, u = 0.5) 
#'
#' # log-density at 0.03 for u = 0.5
#' dBVshear(x = 0.03, u = 0.5, log = TRUE) 
#'
#' # P(X <= 0.03) for u = 0.5
#' pBVshear(q = 0.03, u = 0.5) 
#'
#' # P(X > 0.03) for u = 0.5
#' pBVshear(q = 0.03, u = 0.5, lower.tail = FALSE) 
#'
#' # median of the distribution for u = 0.5
#' qBVshear(p = 0.5, u = 0.5)
#'
#' set.seed(1234)
#' tau <- rBVshear(n = 1E4, u = 0.5)
#' eps <- 0.005
#' hist(tau, freq = FALSE, breaks = seq(0, max(tau) + eps, eps), xlim = c(0, 0.2), 
#'      ylim = c(0, 40), xlab = "Contact time", main = "Histogram of contact times")
#' @aliases dBVshear
dBVshear <- function(x, u, ...) {
  u <- 0.0346930*u
  Re <- 0.452985*u^3
  s <- 0.003480316*Re^(-1/4)
  alpha <- (197.899543379*Re^(-1/12))/(197.899543379*Re^(-1/12) - 1)
  out <- PearsonDS::dpearsonVI(x, a = alpha, b = alpha, location = 0, scale = s, ...)
  return(out)
}
#' @rdname bvshear
#' @export
#' @aliases pBVshear
pBVshear <- function(q, u, ...) {
  u <- 0.0346930*u
  Re <- 0.452985*u^3
  s <- 0.003480316*Re^(-1/4)
  alpha <- (197.899543379*Re^(-1/12))/(197.899543379*Re^(-1/12) - 1)
  out <- PearsonDS::ppearsonVI(q, a = alpha, b = alpha, location = 0, scale = s, ...)
  return(out)
}
#' @rdname bvshear
#' @export
#' @aliases qBVshear
qBVshear <- function(p, u, ...) {
  u <- 0.0346930*u
  Re <- 0.452985*u^3
  s <- 0.003480316*Re^(-1/4)
  alpha <- (197.899543379*Re^(-1/12))/(197.899543379*Re^(-1/12) - 1)
  out <- PearsonDS::qpearsonVI(p, a = alpha, b = alpha, location = 0, scale = s, ...)
  return(out)
}
#' @rdname bvshear
#' @export
#' @aliases rBVshear
rBVshear <- function(n, u, ...) {
  u <- 0.0346930*u
  Re <- 0.452985*u^3
  s <- 0.003480316*Re^(-1/4)
  alpha <- (197.899543379*Re^(-1/12))/(197.899543379*Re^(-1/12) - 1)
  out <- PearsonDS::rpearsonVI(n, a = alpha, b = alpha, location = 0, scale = s, ...)
  return(out)
}
