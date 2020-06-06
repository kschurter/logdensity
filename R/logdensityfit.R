#' Estimate the derivatives of the log-density and the log-density itself
#'
#' @param data numeric vector of observations
#' @param x points at which to estimate (if shorter than h, recycled to length
#' of h). \code{logdensity.fit} only accepts scalar x.
#' @param h bandwidth (if shorter than x, recycled to the length of x). \code{logdensity.fit} only accepts scalar x.
#' @param g function of \code{u}, \code{zl}, and \code{zr} that must be equal to an \code{S}-length vector of
#' 0's at \code{zl = max((minx-x)/h,-1)} and \code{zr = min((maxx-x),1)}, where \code{S} is the order
#' of the local polynomial approximation to the log-density. Function must be
#' vectorized so that \code{g(u, zl, zr)} returns a matrix that is \code{length(u)} by \code{S}.
#' @param dg function that evaluates derivative of \code{g} with respect to \code{u}. Must be
#' vectorized and return a matrix that is \code{length(u)} by \code{S}.
#' @param m kernel function used to compute the density. Can be a function,
#' symbol, or character string that matches the name of a function or one of the
#' kernels allowed in evalkernel. If \code{m == "epanechnikov"} and the polynomial
#' order is 1, an exact solution is computed. Otherwise, the estimate of the 
#' log-density involves numerical integration.
#' @param minx lower bound of support of \code{x}
#' @param maxx upper bound of support of \code{x}
#' @param logf logical indicating whether the log-density should be compute,
#' in addition to its derivative(s).
#' @param exact logical indicating whether an exact solution should be used 
#' (if available) or numerical integration. Exact solution currently only
#' available with epanechnikov kernel and local linear approximation.
#' @param ... Further arguments supplied to \code{g} and \code{dg}
#'
#' @return \code{logdensity.fit} returns a numeric vector containing the log-density (if \code{logf == TRUE}) and its derivatives.
#' 
#' @rdname logdensity
#' 
#' @importFrom stats integrate
#' 
#' @export

logdensity.fit <- function(data, x, h, g, dg, m, minx, maxx, logf, exact, ...){
  zl <- min((x - minx) / h, 1)
  zr <- min((maxx - x) / h, 1)
  if(any(g(-zl, zl, zr,...) != 0) || anyNA(g(-zl,zl,zr,...))){
    stop(gettextf("g(%g) must equal 0 for estimating density at x=%g with bandwidth %g",-zl, x, h), domain = NA)
  }
  if(any(g(zr, zl, zr,...) != 0) || anyNA(g(-zl,zl,zr,...))){
    stop(gettextf("g(%g) must equal 0 for estimating density at x=%g with bandwidth %g",zr, x, h), domain = NA)
  }
  u <- (data - x) / h
  u <- u[(u <= zr) & (u >= -zl)]
  gu <- g(u, zl, zr, ...)
  S <- ncol(gu)
  if(S > length(u)){
    warning(gettextf("Can't do it. Only %d non-missing observations available to estimate %d derivatives at %g using bandwidth %g.", length(u), S, x, h))
    return(rep(NA, S + 1))
  }
  dgu <- dg(u, zl, zr, ...)
  R <- -1 * crossprod(gu, outer(u, 0:(S-1), `^`))
  if(!all(is.finite(R))){
    stop("Encountered a non-finite value of g", domain = NA)
  }
  l <-  colSums(dgu)
  if(!all(is.finite(l))){
    stop("Encountered a non-finite value of dg", domain = NA)
  }
  beta <- numeric(S + 1)
  beta[-1] <- drop(solve(R, l)) * factorial(0:(S-1)) / h^(1:S)
  if(logf){
    if(S == 1 & exact){ # Use an exact solution for denominator instead of numerical integration
      num <- 0.75 * sum(1-u^2) / length(data) 
      fn <- function(u) (exp(beta[-1] * h * u) * (2 + beta[-1] * h * u * (-2 + beta[-1] * h * u)))/(beta[-1] * h)^3
      den <- 0.75 * (fn(-1 * zl) - fn(zr) + (exp(zr * beta[-1] * h) - exp(-1*zl * beta[-1] * h))/(beta[-1] * h))
    }else{ # Use numerical integration
      num <- sum(m(u)) / length(data)
      den <-  tryCatch(integrate(function(w) m(w) * exp( outer(w * h, 1:S, `^`) %*% (beta[-1] / factorial(1:S))), lower = -1*zl, upper = zr)$value,
                       error = function(e) {warning(gettextf("Numerical integration failed at x=%g. Returning NA for log-density. Derivative estimates still valid. Error message from integrate(): %s", x, e$message)); return(NA)})
    }
    beta[1] <- log(num) - log(den) - log(h)
  }
  return(beta)
}
