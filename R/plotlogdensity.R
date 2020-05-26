#' Plot density and log-density estimates
#'
#' @param x a logdensity object
#' @param density logical indicating whether to plot density or log-density
#' @param type what type of plot should be drawn
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param ... further arguments passed to \code{plot.default}
#' 
#' @importFrom graphics plot.default
#' @method plot logdensity
#' @rdname plot.logdensity
#' @seealso \code{\link[graphics]{plot.default}}
#' @export
plot.logdensity <- function(x, density = FALSE, type = "l", xlab, ylab, ...){
  if(missing(xlab)){
    h <- attr(x, "h")
    xlab <- paste("n: ", attr(x,"n"), "  bandwidth: ", if(length(h) == 1L) h else "variable")
  }
  if(missing(ylab)){
    ylab <- if(density) "Density" else "Log-density"
  }
  y <- if(density) exp(x[1, ]) else x[1, ]
  plot.default(attr(x, "x"), y, xlab = xlab, ylab = ylab, type = type, ...)
  invisible(NULL)
}