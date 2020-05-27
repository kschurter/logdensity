#' A local polynomial log-density estimator
#'
#' \code{logdensity.fit} is intended to be called from within
#' \code{logdensity}, after performing basic argument verification. Use
#' caution when calling \code{logdensity.fit} directly.
#'
#' @param data numeric vector of observations
#' @param x points at which to estimate (if shorter than h, recycled to length
#' of h)
#' @param h bandwidth (if shorter than x, recycled to the length of x)
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
#' @param S degree of polynomial expansion of log-density to be used with
#' default \code{g}. If user supplies \code{g} and \code{dg}, this argument is ignored without
#' warning.
#' @param logf logical indicating whether the log-density should be compute,
#' in addition to its derivative(s).
#' @param mc.cores integer number of cores to use with \code{mcmapply}. If equal to 1
#' (default), \code{mapply} will be used to loop over \code{x}, instead.
#' @param ... further arguments supplied to \code{g} and \code{dg}
#'
#' @return an object of class \code{logdensity} which inherits from \code{matrix}.
#'     The \code{S+1} by \code{length(x)} matrix of estimated log-densities (\code{NA} unless \code{logf}
#'     is \code{TRUE}) and derivatives has the following additional attributes:
#'     \code{x}       vector of points at which the estimates were computed
#'     \code{n}       number of non-missing observations used in estimation
#'     \code{h}       bandwidth(s) used
#'     \code{call}    matched call
#'
#' @examples
#' dat <- rchisq(n = 100, df = 2)
#' x <- seq(from = 0, to = 2, length.out = 20)
#' 
#' ## fixed bandwidth
#' ld <- logdensity(data = dat, x = x, h = 0.5, m = "epanechnikov", minx = 0, S = 1, logf = TRUE)
#' print(ld)
#' plot(ld)
#' 
#' ## variable bandwidth
#' h <- pmax(1-x, 0.5)
#' ld <- logdensity(data = dat, x = x, h = h, minx = 0, S = 2)
#' ld
#' plot(ld)
#' 
#' ## Faa di Bruno's formula for the density and its derivatives
#' deriv <- 0L  # integer between 0 and S (=2 for most recent estimation)
#' exp(ld[1, ]) * colSums(bellpoly(ld[-1, ], n = deriv)) # 
#' 
#' ### verify formula for deriv = 0, 1, and 2 
#' exp(ld[1, ]) * colSums(bellpoly(ld[-1, ], n = 0L))  # density (0th derivative)
#' exp(ld[1, ])  # equivalent
#' 
#' exp(ld[1, ]) * colSums(bellpoly(ld[-1, ], n = 1L))  # 1st derivative of density
#' exp(ld[1, ]) * ld[2,]  # equivalent
#' 
#' exp(ld[1, ]) * colSums(bellpoly(ld[-1, ], n = 2L))  # 2nd derivative of density
#' exp(ld[1, ]) * (ld[2,]^2 + ld[3, ])  # equivalent
#' 
#' @importFrom parallel mcmapply
#' @seealso \code{\link{mapply}}, \code{\link[parallel]{mcmapply}}, \code{\link[stats]{integrate}}, \code{\link[logdensity]{bellpoly}}
#' 
#' @references Pinkse, J. and Schurter, K. (2020) "Estimates of derivatives of (log) densities and related objects."
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
  kern <- pmatch(as.character(substitute(m)), kernchoices, nomatch = 0L)
  if(logf){
    m <- if(kern){
      function(u) evalkernel(u, kernchoices[kern])
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
  numCol <- max(length(x), length(h))
  if(numCol > 1){
    if(mc.cores > 1L){
        ld <- parallel::mcmapply(FUN = logdensity.fit, x = x, h = h, MoreArgs = list(data = data, g = g, dg = dg, m = m, minx = minx, maxx = maxx, logf = logf, exact = (kern==1L),...), SIMPLIFY = TRUE, mc.cores = mc.cores)
    }else{
      ld <- mapply(FUN = logdensity.fit, x = x, h = h, MoreArgs = list(data = data, g = g, dg = dg, m = m, minx = minx, maxx = maxx, logf = logf, exact = (kern==1L), ...), SIMPLIFY = TRUE)
    }
  }else{
    ld <- logdensity.fit(data = data, x = x, h = h, g = g, dg = dg, m = m, minx = minx, maxx = maxx, logf = logf, exact = (kern==1L), ...)
  }
  z <- structure(ld, x = x, n = length(data), h = h, call = cl, dim = c(S+1, numCol), class = c("logdensity", "matrix"))
  return(z)
}