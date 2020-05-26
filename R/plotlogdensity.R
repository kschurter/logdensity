#' Plot density and log-density estimates
#'
#' @param z a logdensity object
#' @param main an overall title for the plots
#' @param xlab a title for the x axes
#' @param type what type of plot should be drawn
#' @param ... further arguments passed to plot
#' 
#' @export
plot.logdensity <- function(z, main = NULL, xlab = gettextf("n = %g", z$n), type = "l",...){
  par(mfrow = c(1,2))
  plot.default(z$x, z$logdensity[1,], main = main, xlab = xlab, ylab = "Log-density", type = type,...)
  plot.default(z$x, exp(z$logdensity[1,]), main = main, xlab = xlab, ylab = "Density", type = type,...)
  invisible(NULL)
}