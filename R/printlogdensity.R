#' Print log-density estimates
#'
#' @param x a logdensity object
#' @param digits minimal number of significant digits
#' @param max approximate max number of entries to be printed
#' @param quote logical, indicating whether strings should be printed with surrounding quotes
#' @param ... further arguments passed to print()
#'
#' @return invisibly returns x
#' @rdname print.logdensity
#' @method print logdensity
#' @seealso \code{\link{print.default}}
#' @export
print.logdensity <- function(x, digits = max(3L, getOption("digits")-3L), max = 100, quote = FALSE, ...){
  cat("Call:\n", paste(deparse(attr(x,"call")), sep = "\n", collapse = "\n"))
  cat("\nNumber of non-missing observations:\n", attr(x,"n"))
  S <- nrow(x) - 1L
  cat("\nOrder of polynomial approximation:\n", S)
  rownms <- c("Log-density","1st derivative", "2nd", "3rd", if(S>3) paste(4:S,"th",sep=""))
  m <- matrix(as.numeric(x), nrow = nrow(x), dimnames = list(rownms[1:(S+1)], format(attr(x,"x"), digits = digits)))
  h <- attr(x, "h")
  if(length(h)==1L){
    cat("\nBandwidth:\n", format(h, digits = digits))
  }else{
    m <- rbind("Bandwidth" = h, m)
  }
  cat("\n\n")
  print.default(m[,1:min(max/nrow(m),ncol(m))], digits = digits, quote = FALSE, right = FALSE, max = NULL, ...)
  invisible(x)
}