#' A local polynomial log-density estimator
#'
#' @param data numeric vector of observations
#' @param x points at which to estimate (if shorter than h, recycled to length
#' of h)
#' @param h bandwidth (if shorter than x, recycled to the length of x)
#' @param g function u, zl, and zr that must be equal to an S-length vector of
#' 0's at zl = max((minx-x)/h,-1) and zr = min((maxx-x),1), where S is the order
#' of the local polynomial approximation to the log-density. Function must be
#' vectorized so that g(u, zl, zr) returns a matrix that is length(u) by S.
#' @param dg function that evaluates derivative of g with respect to u. Must be
#' vectorized and return a matrix that is length(u) by S.
#' @param m kernel function used to compute the density. Can be a function,
#' symbol, or character string that matches the name of a function or one of the
#' kernels allowed in evalkernel. If m == "epanechnikov" and the polynomial
#' order is 1, an exact solution is computed. Otherwise, the estimate of the 
#' log-density involves numerical integration.
#' @param minx lower bound of support of x
#' @param maxx upper bound of support of x
#' @param S degree of polynomial expansion of log-density to be used with
#' default g. If user supplies g and dg, this argument is ignored without
#' warning.
#' @param logf logical indicating whether the log-density should be compute,
#' in addition to its derivative(s).
#' @param mc.cores integer number of cores to use with mcmapply. If equal to 1
#' (default), mapply will be used to loop over x, instead.
#' @param ... further arguments supplied to g and dg
#'
#' @return an object of class "logdensity," which is a list containing
#' logdensity   matrix of estimated log-densities (if logf=TRUE) and derivatives
#' x            a vector of points at which the estimates were computed
#' n            the number of non-missing observations used in estimation
#' h            the bandwidth(s) used
#' call         the matched call
#'
#' @examples
#' dat <- rchisq(n = 100, df = 2)
#' x <- seq(from = 0, to = 2, length.out = 20)
#' logdensity(data = dat, x = x, h = 0.5, m = "epanechnikov", minx = 0, S = 1, logf = TRUE)
#' print(logdensity)
#' plot(logdensity)
#' 
#' @importFrom parallel mcmapply
#' 
#' @export

logdensity <- function(data, x, h, g, dg, m = "epanechnikov", minx = -Inf, maxx = Inf, S = 1, logf = TRUE, mc.cores = 1L, ...){
  if(any(!is.finite(h))){
    stop("h is not finite.")
  }
  if(any(h <= 0)){
    stop("h is not positive.")
  }
  if(any(is.na(data))){
    na_data <- is.na(data)
    data <- data[!na_data]
    warning(gettextf("Dropped %d missing observations", sum(na_data)))
  }
  if(any(data < minx) || any(data > maxx)){
    stop(gettextf("Data are not all between minx=%g and maxx=%g", minx, maxx))
  }
  if(missing(g) | missing(dg)){ # define default g and dg and S
    if(!missing(dg) | !missing(g)){
      stop("User must specify both g and dg or neither.")
    }
    if((S%%1 != 0) | (S < 1L)){
      stop("S must be a positive integer.")
    }
    g <- function(u, zl, zr,...){
      a <- outer(u + zl, 0:(S - 1), `^`)
      return(-1 *a * (u + zl) * (zr - u))
      
    }
    dg <- function(u, zl, zr,...){
      a <- outer(u + zl, 0:(S-1), `^`)
      return(a * (u+zl) - ((a * (zr-u)) %*% diag(1:S, nrow = S, ncol = S)))
    }
  }else{
    g <- match.fun(g)
    dg <- match.fun(dg)
  }
  kernchoices <- c("epanechnikov", "gaussian", "triweight", "uniform", "triangle", "cosinus", "quartic")
  kern <- pmatch(as.character(substitute(m)), choices, nomatch = 0L)
  if(logf){
    m <- if(kern){
      function(u) evalkernel(u, choices[kern])
    }else{
      match.fun(m)
    }
  }
  if(any(x < minx) || any(x > maxx)){
    warning(gettextf("Ignoring attempt(s) to estimate density outside the support [%g, %g].", minx, maxx))
    x <- x[(x >= minx) & (x<= maxx)]
    if(length(x)==0){
      stop("No values of x lie within the support.")
    }
  }
  cl <- match.call()
  if((length(x) > 1) || (length(h) > 1)){
    if(mc.cores > 1L){
        ld <- parallel::mcmapply(FUN = logdensity.fit, x = x, h = h, MoreArgs = list(data = data, g = g, dg = dg, m = m, minx = minx, maxx = maxx, logf = logf, exact = (kern==1L),...), SIMPLIFY = TRUE, mc.cores = mc.cores)
    }else{
      ld <- mapply(FUN = logdensity.fit, x = x, h = h, MoreArgs = list(data = data, g = g, dg = dg, m = m, minx = minx, maxx = maxx, logf = logf, exact = (kern==1L), ...), SIMPLIFY = TRUE)
    }
  }else{
    ld <- logdensity.fit(data = data, x = x, h = h, g = g, dg = dg, m = m, minx = minx, maxx = maxx, logf = logf, exact = (kern==1L), ...)
    ld <- as.matrix(ld)
  }
  z <- list(logdensity = ld, x = x, n = length(data), h = h, call = cl)
  class(z) <- "logdensity"
  return(z)
}