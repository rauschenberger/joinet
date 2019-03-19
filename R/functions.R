
#--- Workhorse function --------------------------------------------------------

#' @export
#' @aliases cornet-package
#' @title
#' Combined regression
#' 
#' @description
#' Implements lasso and ridge regression for dichotomised outcomes.
#' Such outcomes are not naturally but artificially binary.
#' They indicate whether an underlying measurement is greater than a threshold.
#'
#' @param y
#' continuous outcome\strong{:}
#' vector of length \eqn{n}
#' 
#' @param cutoff
#' cut-off point for dichotomising outcome into classes\strong{:}
#' \emph{meaningful} value between \code{min(y)} and \code{max(y)}
#' 
#' @param X
#' features\strong{:}
#' numeric matrix with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#' 
#' @param alpha
#' elastic net mixing parameter\strong{:}
#' numeric between \eqn{0} (ridge) and \eqn{1} (lasso)
#' 
#' @param foldid
#' fold identifiers\strong{:}
#' vector with entries between \eqn{1} and \code{nfolds};
#' or \code{NULL} (balance)
#' 
#' @param nfolds
#' number of folds\strong{:}
#' integer between \eqn{3} and \eqn{n}
#' 
#' @param type.measure
#' loss function for binary classification\strong{:}
#' character \code{"deviance"}, \code{"mse"}, \code{"mae"},
#' or \code{"class"} (see \code{\link[glmnet]{cv.glmnet}})
#' 
#' @param pi
#' pi sequence\strong{:}
#' vector of increasing values in the unit interval;
#' or \code{NULL} (default sequence)
#' 
#' @param npi
#' number of \code{pi} values (weighting)
#' 
#' @param sigma
#' sigma sequence\strong{:}
#' vector of increasing positive values;
#' or \code{NULL} (default sequence)
#' 
#' @param nsigma
#' number of \code{sigma} values (scaling)
#' 
#' @param ...
#' further arguments passed to \code{\link[glmnet]{glmnet}}
#' 
#' @details
#' The argument \code{family} is unavailable, because
#' this function fits a \emph{gaussian} model for the numeric response,
#' and a \emph{binomial} model for the binary response.
#' 
#' Linear regression uses the loss function \code{"deviance"} (or \code{"mse"}),
#' but the loss is incomparable between linear and logistic regression.
#' 
#' The loss function \code{"auc"} is unavailable for internal cross-validation.
#' If at all, use \code{"auc"} for external cross-validation only.
#' 
#' @return
#' Returns an object of class \code{cornet}, a list with multiple slots:
#' \itemize{
#'    \item \code{gaussian}: fitted linear model, class \code{glmnet}
#'    \item \code{binomial}: fitted logistic model, class \code{glmnet}
#'    \item \code{sigma}: scaling parameters \code{sigma},
#'           vector of length \code{nsigma}
#'    \item \code{pi}: weighting parameters \code{pi},
#'           vector of length \code{npi}
#'    \item \code{cvm}: evaluation loss,
#'           matrix with \code{nsigma} rows and \code{npi} columns
#'    \item \code{sigma.min}: optimal scaling parameter,
#'           positive scalar
#'    \item \code{pi.min}: optimal weighting parameter,
#'           scalar in unit interval
#'    \item \code{cutoff}: threshold for dichotomisation
#' }
#' 
#' @seealso
#' Methods for objects of class \code{cornet} include
#' \code{\link[=coef.cornet]{coef}} and
#' \code{\link[=predict.cornet]{predict}}.
#' 
#' @references 
#' Armin Rauschenberger and Enrico Glaab (2019).
#' "Lasso and ridge regression for dichotomised outcomes".
#' \emph{Manuscript in preparation}.
#' 
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' net <- cornet(y=y,cutoff=0,X=X)
#' net
#' 
cornet <- function(y,cutoff,X,alpha=1,npi=101,pi=NULL,nsigma=99,sigma=NULL,nfolds=10,foldid=NULL,type.measure="deviance",...){
  
  #--- temporary ---
  # cutoff <- 0; npi <- 101; pi <- NULL; nsigma <- 99; sigma <- NULL; nfolds <- 10;  foldid <- NULL; type.measure <- "deviance"; logistic <- TRUE
  test <- list()
  test$combined <- TRUE
  
  #--- checks ---
  n <- length(y)
  .check(x=y,type="vector")
  if(all(y %in% c(0,1))){warning("Binary response.",call.=FALSE)}
  .check(x=cutoff,type="scalar",min=min(y),max=max(y))
  if(length(y)!=nrow(X)){stop("Contradictory sample size.",call.=FALSE)}
  .check(x=X,type="matrix",dim=c(n,NA))
  .check(x=alpha,type="scalar",min=0,max=1)
  .check(x=npi,type="scalar",min=1)
  .check(x=pi,type="vector",min=0,max=1,null=TRUE)
  .check(x=nsigma,type="scalar",min=1)
  .check(x=sigma,type="vector",min=.Machine$double.eps,null=TRUE)
  .check(x=nfolds,type="scalar",min=3,max=n)
  .check(x=foldid,type="vector",dim=n,values=seq_len(nfolds),null=TRUE)
  .check(x=type.measure,type="string",values=c("deviance","mse","mae","class"))
  # never auc (min/max confusion)!
  if(!is.null(list(...)$family)){stop("Reserved argument \"family\".",call.=FALSE)}
  
  # binarisation
  z <- 1*(y > cutoff)
  
  # fold identifiers
  if(is.null(foldid)){
    foldid <- palasso:::.folds(y=z,nfolds=nfolds)
  }
  
  #--- model fitting ---
  fit <- list()
  fit$gaussian <- glmnet::glmnet(y=y,x=X,family="gaussian",alpha=alpha,...)
  fit$binomial <- glmnet::glmnet(y=z,x=X,family="binomial",alpha=alpha,...)
   
  #--- sigma sequence ---
  if(is.null(sigma)){
    fit$sigma <- exp(seq(from=log(0.05*stats::sd(y)),
                  to=log(10*stats::sd(y)),length.out=nsigma))
  } else {
    fit$sigma <- sigma
    nsigma <- length(sigma)
  }
  
  #--- pi sequence ---
  if(is.null(pi)){
    fit$pi <- seq(from=0,to=1,length.out=npi)
  } else {
    fit$pi <- pi
    npi <- length(pi)
  }
  
  #--- tuning parameters ---
  lab.pi <- paste0("pi",seq_len(npi))
  lab.sigma <- paste0("si",seq_len(nsigma))
  names(fit$sigma) <- lab.sigma
  names(fit$pi) <- lab.pi
  
  #--- cross-validation ---
  pred <- list()
  pred$y <- matrix(data=NA,nrow=n,ncol=length(fit$gaussian$lambda))
  pred$z <- matrix(data=NA,nrow=n,ncol=length(fit$binomial$lambda))

  if(test$combined){
    dimnames <- list(NULL,lab.sigma,lab.pi)
    pred$combined <- array(data=NA,dim=c(n,nsigma,npi),dimnames=dimnames)
  }

  for(k in seq_len(nfolds)){

    y0 <- y[foldid!=k]
    y1 <- y[foldid==k]
    z0 <- z[foldid!=k]
    z1 <- z[foldid==k]
    X0 <- X[foldid!=k,,drop=FALSE]
    X1 <- X[foldid==k,,drop=FALSE]
    
    # linear regression
    net <- glmnet::glmnet(y=y0,x=X0,family="gaussian",alpha=alpha,...)
    temp_y <- stats::predict(object=net,newx=X1,type="response",s=fit$gaussian$lambda)
    pred$y[foldid==k,seq_len(ncol(temp_y))] <- temp_y
    cvm <- palasso:::.loss(y=y1,fit=temp_y,family="gaussian",type.measure="deviance")[[1]]
    y_hat <- temp_y[,which.min(cvm)]
    
    # logistic regression
    net <- glmnet::glmnet(y=z0,x=X0,family="binomial",alpha=alpha,...)
    temp_z <- stats::predict(object=net,newx=X1,type="response",s=fit$binomial$lambda)
    pred$z[foldid==k,seq_len(ncol(temp_z))] <- temp_z
    cvm <- palasso:::.loss(y=z1,fit=temp_z,family="binomial",type.measure=type.measure)[[1]]
    z_hat <- temp_z[,which.min(cvm)]
    
    # combined regression
    if(test$combined){
      for(i in seq_along(fit$sigma)){
        for(j in seq_along(fit$pi)){
          cont <- stats::pnorm(q=y_hat,mean=cutoff,sd=fit$sigma[i])
          pred$combined[foldid==k,i,j] <- fit$pi[j]*cont + (1-fit$pi[j])*z_hat
        }
      }
    }
    
  }
  
  #--- evaluation ---
  
  # linear loss
  fit$gaussian$cvm <- palasso:::.loss(y=y,fit=pred$y,family="gaussian",type.measure="deviance")[[1]]
  fit$gaussian$lambda.min <- fit$gaussian$lambda[which.min(fit$gaussian$cvm)]
  
  # logistic loss
  fit$binomial$cvm <- palasso:::.loss(y=z,fit=pred$z,family="binomial",type.measure=type.measure)[[1]]
  fit$binomial$lambda.min <- fit$binomial$lambda[which.min(fit$binomial$cvm)]

  # combined loss
  if(test$combined){
    dimnames <- list(lab.sigma,lab.pi)
    fit$cvm <- matrix(data=NA,nrow=nsigma,ncol=npi,dimnames=dimnames)
    for(i in seq_len(nsigma)){
      for(j in seq_len(npi)){
        fit$cvm[i,j] <- palasso:::.loss(y=z,fit=pred$combined[,i,j],family="binomial",type.measure=type.measure)[[1]]
      }
    }
    temp <- which(fit$cvm==min(fit$cvm),arr.ind=TRUE,useNames=TRUE)
    if(nrow(temp)>1){warning("Multiple!",call.=FALSE);temp <- temp[1,,drop=FALSE]}
    fit$sigma.min <- fit$sigma[temp[1]]
    fit$pi.min <- fit$pi[temp[2]]
    if(fit$cvm[names(fit$sigma.min),names(fit$pi.min)]!=min(fit$cvm)){stop("Internal error.")}
  }
  
  #--- return ---
  fit$cutoff <- cutoff
  fit$info <- list(type.measure=type.measure,
                   sd.y=stats::sd(y),
                   "+"=sum(z==1),
                   "-"=sum(z==0),
                   table=table(z),
                   n=n,p=ncol(X),
                   test=as.data.frame(test))

  class(fit) <- "cornet"
  return(fit)
}

#--- Methods -------------------------------------------------------------------

#' @export
#' @title
#' Combined regression
#'
#' @description
#' Prints summary of cornet object.
#'
#' @param x
#' \link[cornet]{cornet} object
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @return
#' Returns sample size \eqn{n},
#' number of covariates \eqn{p},
#' information on dichotomisation,
#' tuned scaling parameter (sigma),
#' tuned weighting parameter (pi),
#' and corresponding loss.
#'
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' net <- cornet(y=y,cutoff=0,X=X)
#' print(net)
#' 
print.cornet <- function(x,...){
  cat("cornet object:\n")
  cat(paste0("n = ",x$info$n,","," p = ",x$info$p),"\n")
  cat(paste0("z = I(y > ",signif(x$cutoff,digits=2),"): "))
  cat(paste0(x$info$"+","+"," vs ",x$info$"-","-"),"\n")
  cat(paste0("sigma.min = ",signif(x$sigma.min,digits=1)),"\n")
  cat(paste0("pi.min = ",round(x$pi.min,digits=2)),"\n")
  type <- x$info$type.measure
  value <- signif(x$cvm[names(x$sigma.min),names(x$pi.min)],digits=2)
  cat(paste0(type," = ",value))
  return(invisible(NULL))
}

#' @export
#' @title
#' Extract estimated coefficients
#'
#' @description
#' Extracts estimated coefficients from linear and logistic regression,
#' under the penalty parameter that minimises the cross-validated loss.
#'
#' @param object
#' \link[cornet]{cornet} object
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @return
#' This function returns a matrix with \eqn{n} rows and two columns,
#' where \eqn{n} is the sample size. It includes the estimated coefficients
#' from linear regression (1st column: \code{"beta"})
#' and logistic regression (2nd column: \code{"gamma"}).
#'
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' net <- cornet(y=y,cutoff=0,X=X)
#' coef(net)
#' 
coef.cornet <- function(object,...){
  
  if(length(list(...))!=0){warning("Ignoring arguments.")}
  
  s <- object$gaussian$lambda.min
  beta <- glmnet::coef.glmnet(object=object$gaussian,s=s)
  
  s <- object$binomial$lambda.min
  gamma <- glmnet::coef.glmnet(object=object$binomial,s=s)
  
  coef <- cbind(beta,gamma)
  colnames(coef) <- c("beta","gamma")
  return(coef)
}

#' @export
#' @title
#' Plot loss matrix
#'
#' @description
#' Plots the loss for different combinations of
#' scaling (sigma) and weighting (pi) parameters.
#'
#' @param x
#' \link[cornet]{cornet} object
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @return
#' This function plots the evaluation loss (\code{cvm}).
#' Whereas the matrix has sigma in the rows, and pi in the columns,
#' the plot has sigma on the \eqn{x}-axis, and pi on the \eqn{y}-axis.
#' For all combinations of sigma and pi, the colour indicates the loss.
#' If the R package \code{RColorBrewer} is installed,
#' blue represents low. Otherwise, red represents low.
#' White always represents high.
#' 
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' net <- cornet(y=y,cutoff=0,X=X)
#' plot(net)
#' 
plot.cornet <- function(x,...){
  
  if(length(list(...))!=0){warning("Ignoring arguments.")}

  k <- 100
  levels <- stats::quantile(x$cvm,probs=seq(from=0,to=1,length.out=k+1))
  
  # colours
  if("RColorBrewer" %in% .packages(all.available=TRUE)){
    pal <- rev(c("white",RColorBrewer::brewer.pal(n=9,name="Blues")))
    col <- grDevices::colorRampPalette(colors=pal)(k)
  } else {
    col <- grDevices::heat.colors(n=k)
  }
  
  nsigma <- length(x$sigma)
  npi <- length(x$pi)
  
  graphics::plot.new()
  graphics::par(xaxs="i",yaxs="i")
  graphics::plot.window(xlim=c(1-0.5,nsigma+0.5),ylim=c(1-0.5,npi+0.5))
  
  graphics::title(xlab=expression(sigma),ylab=expression(pi),cex.lab=1)
  #graphics::.filled.contour(x=seq_along(x$sigma),y=seq_along(x$pi),z=x$cvm,levels=levels,col=col)
  graphics::image(x=seq_along(x$sigma),y=seq_along(x$pi),z=x$cvm,breaks=levels,col=col,add=TRUE)
  graphics::box()
  
  ssigma <- which(x$sigma %in% x$sigma.min)
  spi <- which(x$pi %in% x$pi.min)
  
  if(length(ssigma)==1 & length(spi)==1){
    # axes with labels for tuned parameters
    graphics::axis(side=1,at=c(1,ssigma,nsigma),labels=signif(x$sigma[c(1,ssigma,nsigma)],digits=2))
    graphics::axis(side=2,at=c(1,spi,npi),labels=signif(x$pi[c(1,spi,npi)],digits=2))
    # point for tuned parameters
    graphics::points(x=ssigma,y=spi,pch=4,col="white",cex=1)
  } else {
    # axes with standard labels
    at <- seq(from=1,to=nsigma,length.out=5)
    graphics::axis(side=1,at=at,labels=signif(x$sigma,digits=2)[at])
    at <- seq(from=1,to=nsigma,length.out=5)
    graphics::axis(side=2,at=at,labels=signif(x$pi,digits=2)[at])
    # points for selected parameters
    isigma <- sapply(x$sigma.min,function(y) which(x$sigma==y))
    ipi <- sapply(x$pi.min,function(y) which(x$pi==y))
    graphics::points(x=isigma,y=ipi,pch=4,col="white",cex=1)
  }
  
}

#' @export
#' @title
#' Predict binary outcome
#'
#' @description
#' Predicts the binary outcome with linear, logistic, and combined regression.
#' 
#' @details
#' For linear regression, this function tentatively transforms
#' the predicted values to predicted probabilities,
#' using a Gaussian distribution with a fixed mean (threshold)
#' and a fixed variance (estimated variance of the numeric outcome).
#' 
#' @param object
#' \link[cornet]{cornet} object
#' 
#' @param newx
#' covariates\strong{:}
#' numeric matrix with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#' 
#' @param type
#' \code{"probability"}, \code{"odds"}, \code{"log-odds"}
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' net <- cornet(y=y,cutoff=0,X=X)
#' predict(net,newx=X)
#' 
predict.cornet <- function(object,newx,type="probability",...){
  
  if(length(list(...))!=0){warning("Ignoring arguments.")}
  
  x <- object; rm(object)
  
  test <- x$info$test
  
  .check(x=newx,type="matrix")
  .check(x=type,type="string",values=c("probability","odds","log-odds"))
  
  # linear and logistic
  prob <- list()
  link <- as.numeric(stats::predict(object=x$gaussian,
                  newx=newx,s=x$gaussian$lambda.min,type="response"))
  prob$gaussian <- stats::pnorm(q=link,mean=x$cutoff,sd=x$info$sd.y)
  prob$binomial <- as.numeric(stats::predict(object=x$binomial,
                  newx=newx,s=x$binomial$lambda.min,type="response"))
  
  # combined
  if(test$combined){
    cont <- stats::pnorm(q=link,mean=x$cutoff,sd=x$sigma.min)
    prob$combined <- x$pi.min*cont + (1-x$pi.min)*prob$binomial
  }
  
  # consistency tests
  lapply(X=prob,FUN=function(p) .check(x=p,type="vector",min=0,max=1))
  .equal(link>x$cutoff,prob$gaussian>0.5)

  # transformation
  if(type=="probability"){
    frame <- prob
  } else if(type=="odds"){
    frame <- lapply(X=prob,FUN=function(x) x/(1-x))
  } else if(type=="log-odds"){
    frame <- lapply(X=prob,FUN=function(x) log(x/(1-x)))
  } else {
    stop("Invalid type.",call.=FALSE)
  }
  
  return(as.data.frame(frame))
}

#--- Internal functions --------------------------------------------------------

#' @title
#' Equality
#'
#' @description
#' Verifies whether two or more arguments are identical.
#'
#' @param ...
#' scalars, vectors, or matrices of equal dimensions
#' 
#' @param na.rm
#' remove missing values\strong{:}
#' logical
#' 
#' @examples 
#' cornet:::.equal(1,1,1)
#' 
.equal <- function(...,na.rm=FALSE){
  list <- list(...)
  cond <- vapply(X=list,
                 FUN=function(x) all(x==list[[1]],na.rm=na.rm),
                 FUN.VALUE=logical(1))
  if(any(!cond)){
    stop("Unequal elements.",call.=FALSE)
  }
  return(invisible(NULL))
}

#' @title
#' Arguments
#'
#' @description
#' Verifies whether an argument matches formal requirements.
#'
#' @param x
#' argument
#' 
#' @param type
#' character \code{"string"}, \code{"scalar"}, \code{"vector"}, \code{"matrix"}
#' 
#' @param dim
#' vector/matrix dimensionality\strong{:}
#' integer scalar/vector
#' 
#' @param miss
#' accept missing values\strong{:}
#' logical
#' 
#' @param min
#' lower limit\strong{:}
#' numeric
#' 
#' @param max
#' upper limit\strong{:}
#' numeric
#' 
#' @param values
#' only accept specific values\strong{:}
#' vector 
#' 
#' @param inf
#' accept infinite (\code{Inf} or \code{-Inf}) values\strong{:}
#' logical
#' 
#' @param null
#' accept \code{NULL}\strong{:}
#' logical
#'
#' @examples
#' cornet:::.check(0.5,type="scalar",min=0,max=1)
#' 
.check <- function(x,type,dim=NULL,miss=FALSE,min=NULL,max=NULL,values=NULL,inf=FALSE,null=FALSE){
  name <- deparse(substitute(x))
  if(null && is.null(x)){
    return(invisible(NULL)) 
  }
  if(type=="string"){
    cond <- is.vector(x) & is.character(x) & length(x)==1
  } else if(type=="scalar"){
    cond <- is.vector(x) & is.numeric(x) & length(x)==1
  } else if(type=="vector"){
    cond <- is.vector(x) & is.numeric(x)
  } else if(type=="matrix"){
    cond <- is.matrix(x) & is.numeric(x)
  } else {
    warning("Unknown type.")
  }
  if(!cond){
    stop(paste0("Argument \"",name,"\" does not match formal requirements."),call.=FALSE)
  }
  if(!is.null(dim) && length(dim)==1 && length(x)!=dim){
      stop(paste0("Argument \"",name,"\" has invalid length."),call.=FALSE)
  }
  if(!is.null(dim) && length(dim)>1 && any(dim(x)!=dim,na.rm=TRUE)){
      stop(paste0("Argument \"",name,"\" has invalid dimensions."),call.=FALSE)
  }
  
  #   } else if(length(dim)==2){
  #     if(!is.na(dim[1]) && nrow(x)!=dim[1]){
  #       stop(paste0("Argument \"",name,"\" has invalid row number."),call.=FALSE)
  #     }
  #     if(!is.na(dim[2]) && ncol(x)!=dim[2]){
  #       stop(paste0("Argument \"",name,"\" has invalid column number."),call.=FALSE)
  #     }
  #   } else {
  #     
  #   }
  # }
  
  # if(!is.null(length) && length(x)!=length){
  #   stop(paste0("Argument \"",name,"\" has invalid length."),call.=FALSE)
  # }
  # if(!is.null(nrow) && nrow(x)!=nrow){
  #   stop(paste0("Argument \"",name,"\" has invalid row number."),call.=FALSE)
  # }
  # if(!is.null(ncol) && ncol(x)!=ncol){
  #   stop(paste0("Argument \"",name,"\" has invalid column number."),call.=FALSE)
  # }
  
  
  if(!miss && any(is.na(x))){
    stop(paste0("Argument \"",name,"\" contains missing values."),call.=FALSE)
  }
  if(!is.null(min) && any(x<min)){
    stop(paste0("expecting ",name," >= ",min),call.=FALSE)
  }
  if(!is.null(max) && any(x>max)){
    stop(paste0("expecting ",name," <= ",max),call.=FALSE)
  }
  if(!is.null(values) && any(!x %in% values)){
    stop(paste0("Argument \"",name,"\" contains invalid values."),call.=FALSE)
  }
  if(!inf && any(is.infinite(values))){
    stop(paste0("Argument \"",name,"\ contains infinite values."),call.=FALSE)
  }
  return(invisible(NULL))
}

#--- Application ---------------------------------------------------------------

#' @title
#' Performance measurement
#'
#' @description
#' Compares models for a continuous response with a cut-off value.
#' 
#' @details
#' Uses k-fold cross-validation,
#' fits linear, logistic, and combined regression,
#' calculates different loss functions,
#' and examines squared deviance residuals.
#' 
#' @inheritParams  cornet
#' 
#' @examples
#' \dontshow{n <- 100; p <- 20
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' loss <- cornet:::.compare(y=y,cutoff=0,X=X,nfolds=2)
#' loss}
#' \donttest{n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' loss <- cornet:::.compare(y=y,cutoff=0,X=X)
#' loss}
#' 
.compare <- function(y,cutoff,X,alpha=1,nfolds=5,foldid=NULL,type.measure="deviance"){
  
  z <- 1*(y > cutoff)
  if(is.null(foldid)){
    fold <- palasso:::.folds(y=z,nfolds=nfolds)
  } else {
    fold <- foldid
  }
  
  #--- cross-validated loss ---
  
  cols <- c("gaussian","binomial","combined")
  pred <- matrix(data=NA,nrow=length(y),ncol=length(cols),
                 dimnames=list(NULL,cols))
  
  for(i in seq_len(nfolds)){
    fit <- cornet::cornet(y=y[fold!=i],cutoff=cutoff,X=X[fold!=i,],alpha=alpha,type.measure=type.measure)
    tryCatch(expr=plot.cornet(fit),error=function(x) NULL)
    temp <- predict.cornet(fit,newx=X[fold==i,])
    if(any(temp<0|temp>1)){stop("Outside unit interval.",call.=FALSE)}
    model <- colnames(pred)
    for(j in seq_along(model)){
      pred[fold==i,model[j]] <- temp[[model[j]]]
    }
  }
  
  type <- c("deviance","class","mse","mae","auc")
  loss <- lapply(X=type,FUN=function(x) palasso:::.loss(y=z,fit=pred,family="binomial",type.measure=x,foldid=fold)[[1]])
  names(loss) <- type
  
  #--- deviance residuals ---
  
  # squared deviance residuals
  limit <- 1e-05
  pred[pred < limit] <- limit
  pred[pred > 1 - limit] <- 1 - limit
  res <- -2 * (z * log(pred) + (1 - z) * log(1 - pred))
  rxs <- res[,"binomial"]
  rys <- res[,"combined"]
  
  # residual increase/decrease
  loss$resid.factor <- stats::median((rys-rxs)/rxs)
  
  # paired test for each fold
  loss$resid.pvalue <- numeric()
  for(i in seq_len(nfolds)){
    cond <- fold==i
    loss$resid.pvalue[i] <- stats::wilcox.test(x=rxs[cond],y=rys[cond],
                                               paired=TRUE,alternative="greater")$p.value
  }
  
  return(loss)
  
}

#' @title
#' Single-split test
#'
#' @description
#' Compares models for a continuous response with a cut-off value.
#' 
#' @details
#' Splits samples into 80% for training and 20% for testing,
#' calculates squared deviance residuals of logistic and combined regression,
#' conducts the paired one-sided Wilcoxon signed rank test,
#' and returns the p-value.
#' 
#' @inheritParams cornet
#' 
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' cornet:::.test(y=y,cutoff=0,X=X)
#' 
.test <- function(y,cutoff,X,alpha=1,type.measure="deviance"){
  
  z <- 1*(y > cutoff)
  fold <- palasso:::.folds(y=z,nfolds=5)
  fold <- ifelse(fold==1,1,0)
  
  fit <- cornet::cornet(y=y[fold==0],cutoff=cutoff,X=X[fold==0,],alpha=alpha)
  tryCatch(expr=plot.cornet(fit),error=function(x) NULL)
  pred <- predict.cornet(fit,newx=X[fold==1,])
  if(any(pred<0|pred>1)){stop("Outside unit interval.",call.=FALSE)}
  
  #res <- (pred-z[fold==1])^2 # MSE
  #pvalue <- wilcox.test(x=res[,"binomial"],y=res[,"combined"],paired=TRUE,alternative="greater")$p.value
  #colMeans(abs(pred-0.5)) # distance from 0.5
  
  limit <- 1e-05
  pred[pred < limit] <- limit
  pred[pred > 1 - limit] <- 1 - limit
  res <- -2 * (z[fold==1] * log(pred) + (1 - z[fold==1]) * log(1 - pred))
  pvalue <- stats::wilcox.test(x=res[,"binomial"],y=res[,"combined"],paired=TRUE,alternative="greater")$p.value
  
  return(pvalue)
}

#' @title
#' Data simulation
#'
#' @description
#' Simulates data for unit tests
#' 
#' @param n
#' sample size\strong{:}
#' positive integer
#' 
#' @param p
#' covariate space\strong{:}
#' positive integer
#' 
#' @param prob
#' (approximate) proportion of causal covariates\strong{:}
#' numeric between \eqn{0} and \eqn{1}
#' 
#' @param fac
#' noise factor\strong{:}
#' positive real number
#' 
#' @return
#' Returns invisible list with elements \code{y} and \code{X}.
#' 
#' @examples
#' data <- cornet:::.simulate(n=10,p=20,prob=0.2,fac=2)
#' names(data)
#' 
.simulate <- function(n,p,prob=0.2,fac=1){
  beta <- stats::rnorm(n=p)
  cond <- stats::rbinom(n=p,size=1,prob=prob)
  beta[cond==0] <- 0
  X <- matrix(stats::rnorm(n=n*p),nrow=n,ncol=p)
  mean <- X %*% beta
  y <- stats::rnorm(n=n,mean=mean,sd=fac*stats::sd(mean))
  return(invisible(list(y=y,X=X)))
}


#--- Legacy --------------------------------------------------------------------

# # Import this function from the palasso package.
# .loss <- function (y,fit,family,type.measure,foldid=NULL){
#   if (!is.list(fit)) {
#     fit <- list(fit)
#   }
#   loss <- list()
#   for (i in seq_along(fit)) {
#     if (is.vector(fit[[i]])) {
#       fit[[i]] <- as.matrix(fit[[i]])
#     }
#     if (is.null(foldid) & (family == "cox" | type.measure == 
#                            "auc")) {
#       stop("Missing foldid.", call. = FALSE)
#     }
#     if (family == "gaussian") {
#       if (type.measure %in% c("deviance", "mse")) {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean((x - y)^2))
#       }
#       else if (type.measure == "mae") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean(abs(x - y)))
#       }
#       else {
#         stop("Invalid type measure.", call. = FALSE)
#       }
#     }
#     else if (family == "binomial") {
#       if (type.measure == "deviance") {
#         limit <- 1e-05
#         fit[[i]][fit[[i]] < limit] <- limit
#         fit[[i]][fit[[i]] > 1 - limit] <- 1 - limit
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean(-2 * (y * log(x) + (1 - 
#                                                                         y) * log(1 - x))))
#       }
#       else if (type.measure == "mse") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) 2 * mean((x - y)^2))
#       }
#       else if (type.measure == "mae") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) 2 * mean(abs(x - y)))
#       }
#       else if (type.measure == "class") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean(abs(round(x) - y)))
#       }
#       else if (type.measure == "auc") {
#         weights <- table(foldid)
#         cvraw <- matrix(data = NA, nrow = length(weights), 
#                         ncol = ncol(fit[[i]])) # typo in palasso package !
#         for (k in seq_along(weights)) {
#           cvraw[k, ] <- apply(X = fit[[i]], MARGIN = 2, 
#                               FUN = function(x) glmnet::auc(y = y[foldid == 
#                                                                     k], prob = x[foldid == k]))
#         }
#         loss[[i]] <- apply(X = cvraw, MARGIN = 2, FUN = function(x) stats::weighted.mean(x = x, 
#                                                                                          w = weights, na.rm = TRUE))
#         names(loss[[i]]) <- colnames(fit[[i]]) # typo in palasso package!
#       }
#       else {
#         stop("Invalid type measure.", call. = FALSE)
#       }
#     }
#     else if (family == "poisson") {
#       if (type.measure == "deviance") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean(2 * (ifelse(y == 0, 
#                                                               0, y * log(y/x)) - y + x), na.rm = TRUE))
#       }
#       else if (type.measure == "mse") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean((x - y)^2))
#       }
#       else if (type.measure == "mae") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean(abs(x - y)))
#       }
#       else {
#         stop("Invalid type measure.", call. = FALSE)
#       }
#     }
#     else if (family == "cox") {
#       if (type.measure == "deviance") {
#         weights <- tapply(X = y[, "status"], INDEX = foldid, 
#                           FUN = sum)
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) stats::weighted.mean(x = x/weights, 
#                                                                   w = weights, na.rm = TRUE))
#       }
#       else {
#         stop("Invalid type measure.", call. = FALSE)
#       }
#     }
#     else {
#       stop("Invalid family.", call. = FALSE)
#     }
#     if (sum(diff(is.na(loss[[i]]))) == 1) {
#       loss[[i]] <- loss[[i]][!is.na(loss[[i]])]
#     }
#   }
#   return(loss)
# }
# 
# # Import this function from the palasso package.
# .folds <- function (y, nfolds, foldid = NULL){
#   if(!is.null(foldid)){
#     return(foldid)
#   }
#   #if (survival::is.Surv(y)){ # active in palasso
#   #  y <- y[, "status"] # active in palasso
#   #} # active in palasso
#   if(all(y %in% c(0, 1))){
#     foldid <- rep(x = NA, times = length(y))
#     foldid[y == 0] <- sample(x = rep(x = seq_len(nfolds), 
#                                      length.out = sum(y == 0)))
#     foldid[y == 1] <- sample(x = rep(x = seq_len(nfolds), 
#                                      length.out = sum(y == 1)))
#   } else {
#     foldid <- sample(x = rep(x = seq_len(nfolds), length.out = length(y)))
#   }
#   return(foldid)
# }
