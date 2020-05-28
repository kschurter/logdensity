#' The logdensity package
#' 
#' A package for estimating (log) densities and their derivatives
#'
#' This package contains functions that estimate the logarithm of an unknown density function from iid data.
#' The estimation strategy is based on a local polynomial appoximation to the logarithm of the density. The
#' estimates behave well near the boundary of the support of the density and can be guaranteed to be nonnegative.
#' Details can be found in Pinske and Schurter (2020).
#' 
#' @aliases logdensity-package
#' @docType package
#' @name logdensity-package
#' @author Joris Pinkse and Karl Schurter (maintainer, \email{kschurter@@psu.edu})
#' 
#' @references Pinkse, J. and Schurter, K. (2020) ''Estimates of derivatives of (log) densities and related objects.''
NULL
