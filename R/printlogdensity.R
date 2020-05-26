#' Print log-density estimates
#'
#' @param z a logdensity object
#' @param digits minimal number of significant digits
#' @param max approximate max number of entries to be printed
#' @param quote logical, indicating whether strings should be printed with surrounding quotes
#' @param ... further arguments passed to print()
#'
#' @return invisibly returns z
#' @export
print.logdensity <- function(z, digits = max(3L, getOption("digits")-3L), max = 6 * ncol(m), quote = FALSE, ...){
  cat("Call:\n", paste(deparse(z$call), sep = "\n", collapse = "\n"))
  cat("\nNumber of non-missing observations:\n", z$n)
  S <- nrow(z$logdensity) - 1L
  cat("\nOrder of polynomial approximation:\n", S)
  m <- cbind(z$x, t(z$logdensity))
  colnms <- c("x","Log-density","1st derivative", "2nd", "3rd", if(S>3) paste(4:S,"th",sep=""))
  dimnames(m)[[2]]<- colnms[1:(S+2)]
  if(length(z$h)==1L){
    cat("\nBandwidth:\n", format(z$h, digits = digits))
  }else{
    m <- cbind("Bandwidth" = z$h, m)
  }
  cat("\n")
  print(m, digits = digits, quote = FALSE, max = max, ...)
  invisible(z)
}