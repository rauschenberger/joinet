% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{coef.joinet}
\alias{coef.joinet}
\title{Extract Coefficients}
\usage{
\method{coef}{joinet}(object, ...)
}
\arguments{
\item{object}{\link[joinet]{joinet} object}

\item{...}{further arguments (not applicable)}
}
\value{
This function returns the pooled coefficients.
The slot \code{alpha} contains the intercepts
in a vector of length \eqn{q},
and the slot \code{beta} contains the slopes
in a matrix with \eqn{p} rows (inputs) and \eqn{q} columns.
}
\description{
Extracts pooled coefficients.
(The meta learners linearly combines
the coefficients from the base learners.)
}
\examples{
\dontshow{
if(!grepl('SunOS',Sys.info()['sysname'])){
n <- 50; p <- 100; q <- 3
X <- matrix(rnorm(n*p),nrow=n,ncol=p)
Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
object <- joinet(Y=Y,X=X)
coef <- coef(object)}}
\dontrun{
n <- 50; p <- 100; q <- 3
X <- matrix(rnorm(n*p),nrow=n,ncol=p)
Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
object <- joinet(Y=Y,X=X)
coef <- coef(object)}

}
