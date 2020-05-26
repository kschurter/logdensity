#' Evaluate a kernel function
#'
#' @param u numeric array
#' @param kernel character string specifying name of kernel
#'
#' @return array with same dimensions as \code{u} containing kernel evaluated at \code{u}.
#' Note: \code{NA} values of \code{u} evaluate to 0.
#' 
#' @export
evalkernel <- function(u, kernel = c("epanechnikov","gaussian", "triweight", "uniform", "triangle", "cosinus", "quartic")){
  kernel <- match.arg(kernel)
  rval <- switch(kernel,
                 gaussian = 1 / sqrt(2 * pi) * exp(-(u^2) / 2),
                 epanechnikov = (0.75 * (1 - u^2)) * (abs(u) <= 1),
                 triweight = (1.09375 * (1 - u^2)^3) * (abs(u) <= 1),
                 uniform = rep(1/2, length(u)) * (abs(u) <= 1),
                 triangle = (1 - abs(u)) * (abs(u) <= 1),
                 cosinus = (pi / 4 * cos(pi / 2 * u)) * (abs(u) <= 1),
                 quartic = (0.9375 * (1 - u^2)^2) * (abs(u) <= 1))
  rval[!is.finite(rval)] <- 0
  rval
}