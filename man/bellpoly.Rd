% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bellpoly.R
\name{bellpoly}
\alias{bellpoly}
\title{Evaluate partial Bell polynomials}
\usage{
bellpoly(x, n = nrow(x), k = min(1, n):n)
}
\arguments{
\item{x}{a matrix or object that can be coerced into one by \code{as.matrix}}

\item{n}{nonnegative integer representing the number of elements}

\item{k}{vector of integers no greater than \code{n} representing the sizes of partition}
}
\value{
\code{bellpoly} returns a vector or matrix depending on the dimensions of \code{x} that 
contains the value of the partial Bell polynomial(s) evaluated at columns of \code{as.matrix(x)}.
}
\description{
Evaluates the partial Bell polynomials, also known as the incomplete
Bell polynomials. The complete Bell polynomial is the sum of the 
partial Bell polynomials over \eqn{k=1,\dots,n}. The partial polynomials
are evaluted using the recurrence relation
\deqn{k\, B_{n,k}(x_1,\dots,x_{n-k+1}) = \sum_{r=k-1}^{n-1}{n \choose r}x_{n-r}\,B_{r,k-1}(x_1,\dots,x_{r-k+2})}{k B_{n,k}(x_1,\dots,x_{n-k+1}) = \sum_{r=k-1}^{n-1} choose(n,r)x_{n-r}B_{r,k-1}(x_1,\dots,x_{r-k+2})}
with \eqn{B_{0,0} = 1}, \eqn{B_{0,k} = 0} for \eqn{k \ge 1}, and \eqn{B_{n,0} = 0} for \eqn{n \ge 1}.
}
\examples{
n <- 5
k <- 1:n

## Idempotent numbers
bellpoly(k)
choose(n, k) * k^(n - k) # equivalent

## Stirling numbers (unsigned) of the first kind
bellpoly(x = factorial(k-1))

## Stirling numbers of the second kind
bellpoly(x = rep(1, n))

}
\concept{bell polynomial combinatoric}
