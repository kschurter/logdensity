#' Evaluate partial Bell polynomials
#'
#' Evaluates the partial Bell polynomials, also known as the incomplete
#' Bell polynomials. The complete Bell polynomial is the sum of the 
#' partial Bell polynomials over \eqn{k=1,\dots,n}. The partial polynomials
#' are evaluted using the recurrence relation
#' \deqn{k B_{n,k}(x_1,\dots,x_{n-k+1}) = \sum_{r=k-1}^{n-1}{n \choose r}x_{n-r}B_{r,k-1}(x_1,\dots,x_{r-k+2})}
#' with \eqn{B_{0,0} = 1}, \eqn{B_{0,k} = 0} for \eqn{k \ge 1}, and \eqn{B_{n,0} = 0} for \eqn{n \ge 1}.
#'
#' @param x a matrix or object that can be coerced into one by \code{as.matrix}
#' @param n nonnegative integer representing the number of elements
#' @param k vector of integers no greater than \code{n} representing the sizes of partition
#'
#' @return \code{bellpoly} returns a vector or matrix matching the dimensions of \code{x} that 
#' contains the value of the partial Bell polynomial evaluated at columns of \code{as.matrix(x)}.
#' @export
#' @concept bell polynomial combinatoric
#'
#' @examples
#' n <- 5
#' k <- 1:n
#' 
#' ## Idempotent numbers
#' bellpoly(k)
#' choose(n, k) * k^(n - k) # equivalent
#' 
#' ## Stirling numbers (unsigned) of the first kind
#' bellpoly(x = factorial(k-1))
#' 
#' ## Stirling numbers of the second kind
#' bellpoly(x = rep(1, n))
#' 
bellpoly <- function(x, n = nrow(x), k = min(1,n):n){
  diminput <- dim(x)
  x <- as.matrix(x)
  if(length(n)>1L || (n < 0L)) stop("n must be a single nonnegative integer")
  nc <- ncol(x)
  if(n == 0){
    Bnk <- rep(ifelse(k==0,1,0), nc)
  }else{
    recur <- function(m, l) {
      if(l <= 1L){
        if(l < 1L) rep(0, nc) else x[m, ]
      }
      else{
        i <- (l-1L):(m-1L)
        Bnj <- t(sapply(i,FUN = recur, l = l-1L))
        dim(Bnj) <- c(length(i), nc)
        colSums(choose(m, i) * x[m-i, ] * Bnj ) 
      }
    }
    Bnk <- t(sapply(k, FUN = recur, m = n)) / factorial(k)
  }
  dim(Bnk) <- c(length(k), diminput[2])
  Bnk
}
