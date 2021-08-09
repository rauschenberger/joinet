
#' @name joinet-package
#' @keywords documentation
#' @docType package
#' 
#' @title
#' Multivariate Elastic Net Regression
#' 
#' @description
#' The R package \code{joinet} implements multivariate
#' ridge and lasso regression using stacked generalisation.
#' This multivariate regression typically outperforms
#' univariate regression at predicting correlated outcomes.
#' It provides predictive and interpretable models
#' in high-dimensional settings.
#' 
#' @details
#' Use function \code{\link{joinet}} for model fitting.
#' Type \code{library(joinet)} and then \code{?joinet} or
#' \code{help("joinet)"} to open its help file.
#' 
#' See the vignette for further examples.
#' Type \code{vignette("joinet")} or \code{browseVignettes("joinet")}
#' to open the vignette.
#' 
#' @references
#' Armin Rauschenberger, Enrico Glaab (2021)
#' "Predicting correlated outcomes from molecular data"
#' \emph{Bioinformatics}. btab576
#' \doi{10.1093/bioinformatics/btab576}
#' 
#' \email{armin.rauschenberger@uni.lu}
#'
#' @examples
#' \dontshow{
#' if(!grepl('SunOS',Sys.info()['sysname'])){
#' #--- data simulation ---
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' # n samples, p inputs, q outputs
#' 
#' #--- model fitting ---
#' object <- joinet(Y=Y,X=X)
#' # slot "base": univariate
#' # slot "meta": multivariate
#' 
#' #--- make predictions ---
#' y_hat <- predict(object,newx=X)
#' # n x q matrix "base": univariate
#' # n x q matrix "meta": multivariate 
#' 
#' #--- extract coefficients ---
#' coef <- coef(object)
#' # effects of inputs on outputs
#' # q vector "alpha": intercepts
#' # p x q matrix "beta": slopes
#' 
#' #--- model comparison ---
#' loss <- cv.joinet(Y=Y,X=X)
#' # cross-validated loss
#' # row "base": univariate
#' # row "meta": multivariate
#' }}
#' \dontrun{
#' #--- data simulation ---
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' # n samples, p inputs, q outputs
#' 
#' #--- model fitting ---
#' object <- joinet(Y=Y,X=X)
#' # slot "base": univariate
#' # slot "meta": multivariate
#' 
#' #--- make predictions ---
#' y_hat <- predict(object,newx=X)
#' # n x q matrix "base": univariate
#' # n x q matrix "meta": multivariate 
#' 
#' #--- extract coefficients ---
#' coef <- coef(object)
#' # effects of inputs on outputs
#' # q vector "alpha": intercepts
#' # p x q matrix "beta": slopes
#' 
#' #--- model comparison ---
#' loss <- cv.joinet(Y=Y,X=X)
#' # cross-validated loss
#' # row "base": univariate
#' # row "meta": multivariate
#' }
#' 
NULL
