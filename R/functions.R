
#--- Main function -------------------------------------------------------------

#' @export
#' @aliases joinet-package
#' @title
#' Multivariate Elastic Net Regression
#' 
#' @description
#' Implements multivariate elastic net regression.
#'  
#' @param Y
#' outputs\strong{:}
#' numeric matrix with \eqn{n} rows (samples)
#' and \eqn{q} columns (variables),
#' with positive correlation (see details)
#' 
#' @param X
#' inputs\strong{:}
#' numeric matrix with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#'
#' @param family
#' distribution\strong{:}
#' vector of length \eqn{1} or \eqn{q} with entries
#' \code{"gaussian"}, \code{"binomial"} or \code{"poisson"}
#'
#' @param nfolds
#' number of folds
#'
#' @param foldid
#' fold identifiers\strong{:}
#' vector of length \eqn{n} with entries between \eqn{1} and \code{nfolds};
#' or \code{NULL} (balance)
#' 
#' @param type.measure
#' loss function\strong{:}
#' vector of length \eqn{1} or \eqn{q} with entries
#' \code{"deviance"}, \code{"class"}, \code{"mse"} or \code{"mae"}
#' (see \code{\link[glmnet]{cv.glmnet}})
#'
#' @param alpha.base
#' elastic net mixing parameter for base learners\strong{:}
#' numeric between \eqn{0} (ridge) and \eqn{1} (lasso)
#' 
#' @param alpha.meta
#' elastic net mixing parameter for meta learner\strong{:}
#' numeric between \eqn{0} (ridge) and \eqn{1} (lasso)
#' 
#' @param ...
#' further arguments passed to \code{\link[glmnet]{glmnet}}
#' 
#' @references 
#' A Rauschenberger, E Glaab (2019)
#' "Multivariate elastic net regression through stacked generalisation"
#' \emph{Manuscript in preparation.}
#' 
#' @details
#' The \eqn{q} outcomes should be positively correlated.
#' Avoid negative correlations by changing the sign of the variable.
#' 
#' elastic net mixing parameters:
#' \code{alpha.base} controls input-output effects,
#' \code{alpha.meta} controls output-output effects;
#' ridge (\eqn{0}) renders dense models,
#' lasso (\eqn{1}) renders sparse models
#' 
#' @examples
#' n <- 30; q <- 2; p <- 20
#' Y <- matrix(rnorm(n*q),nrow=n,ncol=q)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- joinet(Y=Y,X=X)
#' 
joinet <- function(Y,X,family="gaussian",nfolds=10,foldid=NULL,type.measure="deviance",alpha.base=0,alpha.meta=0,...){
  
  #--- temporary ---
  # family <- "gaussian"; nfolds <- 10; foldid <- NULL; type.measure <- "deviance"
  
  #--- checks ---
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  
  cornet:::.check(x=Y,type="matrix",miss=TRUE)
  if(any(stats::cor(Y,use="pairwise.complete.obs")<0,na.rm=TRUE)){warning("Negative correlation!",call.=FALSE)}
  cornet:::.check(x=X,type="matrix")
  #cornet:::.check(x=family,type="string",values=c("gaussian","binomial","poisson"))
  if(nrow(Y)!=nrow(X)){stop("Contradictory sample size.",call.=FALSE)}
  cornet:::.check(x=nfolds,type="scalar",min=3)
  cornet:::.check(x=foldid,type="vector",values=seq_len(nfolds),null=TRUE)
  cornet:::.check(x=type.measure,type="string",values=c("deviance","class","mse","mae")) # not auc (min/max confusion)
  cornet:::.check(x=alpha.base,type="scalar",min=0,max=1)
  cornet:::.check(x=alpha.meta,type="scalar",min=0,max=1)
  if(!is.null(c(list(...)$lower.limits,list(...)$upper.limits))){
    stop("Reserved arguments \"lower.limits\" and \"upper.limits\".",call.=FALSE)
  }
  
  #--- dimensionality ---
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  #--- family ---
  if(length(family)==1){
    family <- rep(family,times=q)
  } else if(length(family)!=q){
    stop("Invalid argument family",call.=FALSE)
  }
  
  #--- fold identifiers ---
  # provide foldid as matrix?
  if(is.null(foldid)){
    foldid <- palasso:::.folds(y=Y[,1],nfolds=nfolds) # temporary Y[,1]
  } else {
    nfolds <- length(unique(foldid))
  }
  
  #--- full fit ---
  nlambda <- numeric()
  base <- lapply(seq_len(q),function(x) list())
  for(i in seq_len(q)){
    cond <- !is.na(Y[,i])
    #if(sum(cond)==0){nlambda[i] <- 0; next}
    base[[i]]$glmnet.fit <- glmnet::glmnet(y=Y[cond,i],x=X[cond,],family=family[i],alpha=alpha.base,...) # ellipsis
    base[[i]]$lambda <- base[[i]]$glmnet.fit$lambda
    nlambda[i] <- length(base[[i]]$glmnet.fit$lambda)
  }
  
  #--- predictions ---
  link <- list()
  for(i in seq_len(q)){
    link[[i]] <- matrix(data=NA,nrow=n,ncol=nlambda[i])
  }
  
  #--- base cross-validation ---
  for(k in seq_len(nfolds)){
    Y0 <- Y[foldid!=k,,drop=FALSE]
    Y1 <- Y[foldid==k,,drop=FALSE]
    X0 <- X[foldid!=k,,drop=FALSE]
    X1 <- X[foldid==k,,drop=FALSE]
    for(i in seq_len(q)){
      cond <- !is.na(Y0[,i])
      #if(sum(cond)==0){next}
      object <- glmnet::glmnet(y=Y0[cond,i],x=X0[cond,],family=family[i],alpha=alpha.base,...) # ellipsis
      temp <- stats::predict(object=object,newx=X1,type="link",
                             s=base[[i]]$glmnet.fit$lambda)
      link[[i]][foldid==k,seq_len(ncol(temp))] <- temp
    }
  }
  
  #--- tune base lambdas ---
  for(i in seq_len(q)){
    fit <- .mean.function(link[[i]],family=family[i])
    cond <- !is.na(Y[,i])
    base[[i]]$cvm <- palasso:::.loss(y=Y[cond,i],fit=fit[cond,],
                                     family=family[i],type.measure=type.measure)[[1]]
    base[[i]]$lambda.min <- base[[i]]$lambda[which.min(base[[i]]$cvm)]
  }
  
  #--- predictions ---
  hat <- matrix(NA,nrow=n,ncol=q)
  for(i in seq_len(q)){
    hat[,i] <- link[[i]][,base[[i]]$lambda==base[[i]]$lambda.min]
  }
  
  #--- meta cross-validation ---
  meta <- list()
  for(i in seq_len(q)){
    cond <- !is.na(Y[,i])
    meta[[i]] <- glmnet::cv.glmnet(y=Y[cond,i],x=hat[cond,],
                                   lower.limits=0, # important: 0
                                   upper.limits=Inf, # important: Inf
                                   foldid=foldid[cond],
                                   family=family[i],
                                   type.measure=type.measure,
                                   intercept=TRUE, # with intercept
                                   alpha=alpha.meta,...) # ellipsis
    # consider trying different alpha.meta and selecting best one
  }
  
  #--- return ---
  names(base) <- names(meta) <- paste0("y",seq_len(q))
  info <- data.frame(q=q,p=p,family=family,type.measure=type.measure)
  list <- list(base=base,meta=meta,info=info)
  class(list) <- "joinet"
  return(list)
}

.mean.function <- function(x,family){
  if(family=="gaussian"){
    return(x)
  } else if(family=="binomial"){
    return(1/(1+exp(-x)))
  } else if(family=="poisson"){
    return(exp(x))
  } else {
    stop("Family not implemented.",call.=FALSE)
  }
}

.link.function <- function(x,family){
  if(family=="gaussian"){
    return(x)
  } else if(family=="binomial"){
    return(log(1/(1-x)))
  } else if(family=="poisson"){
    return(log(x))
  } else {
    stop("Family not implemented.",call.=FALSE)
  }
}

#--- Methods for class "joinet" -----------------------------------------------

#' @export
#' @title
#' Make Predictions
#'
#' @description
#' Predicts outcome from features with stacked model.
#' 
#' @param object
#' \link[joinet]{joinet} object
#' 
#' @param newx
#' covariates\strong{:}
#' numeric matrix with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#' 
#' @param type
#' character "link" or "response"
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @examples
#' n <- 30; q <- 2; p <- 20
#' #Y <- matrix(rnorm(n*q),nrow=n,ncol=q)
#' Y <- matrix(rbinom(n=n*q,size=1,prob=0.5),nrow=n,ncol=q)
#' #Y <- matrix(rpois(n=n*q,lambda=4),nrow=n,ncol=q)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- joinet(Y=Y,X=X,family="binomial")
#' y_hat <- predict(object,newx=X)
#' 
predict.joinet <- function(object,newx,type="response",...){
  if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
  
  x <- object; rm(object)
  
  newx <- as.matrix(newx)
  cornet:::.check(x=newx,type="matrix")
  
  q <- length(x$base)
  n <- nrow(newx)
  base <- meta <- matrix(NA,nrow=n,ncol=q)
  
  # base learners
  for(i in seq_len(q)){
    #base[,i] <- as.numeric(stats::predict(object=x$base[[i]]$glmnet.fit,newx=newx,s=x$base[[i]]$lambda.min,type="link"))
    base[,i] <- as.numeric(glmnet::predict.cv.glmnet(object=x$base[[i]],newx=newx,s="lambda.min",type="link"))
    # check whether fine for "binomial" family
  }
  
  # meta learners
  for(i in seq_len(q)){
    meta[,i] <- as.numeric(stats::predict(object=x$meta[[i]],
                                          newx=base,s="lambda.min",type="link"))
  }
  
  list <- list(base=base,meta=meta)
  
  if(type=="response"){
    for(i in seq_len(q)){
      base[,i] <- .mean.function(x=base[,i],family=x$info$family[i])
      meta[,i] <- .mean.function(x=meta[,i],family=x$info$family[i])
    }
  } else if(type!="link"){
    stop("Invalid type.",call.=FALSE)
  }
  
  list <- list(base=base,meta=meta)
  return(list)
  
}

#' @export
#' @title
#' Extract Coefficients
#'
#' @description
#' Extracts pooled coefficients.
#' (The meta learners linearly combines
#' the coefficients from the base learners.)
#' 
#' @param object
#' \link[joinet]{joinet} object
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @examples
#' n <- 30; q <- 2; p <- 20
#' Y <- matrix(rnorm(n*q),nrow=n,ncol=q)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- joinet(Y=Y,X=X)
#' coef <- coef(object)
#' 
coef.joinet <- function(object,...){
  if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
  
  # base coefficients
  base <- list()
  coef <- sapply(object$base,function(x) glmnet::coef.glmnet(object=x$glmnet.fit,s=x$lambda.min))
  base$alpha <- sapply(coef,function(x) x[1,])
  base$beta <- sapply(coef,function(x) x[-1,])
  names(base$alpha) <- colnames(base$beta) <- names(object$base)
  
  # meta coefficients
  meta <- list()
  weights <- weights.joinet(object)
  meta$alpha <- weights[1,]
  meta$beta <- weights[-1,]
  
  # pooled coefficients
  pool <- list()
  pool$alpha <- meta$alpha + base$alpha %*% meta$beta
  pool$beta <- base$beta %*% meta$beta
  
  # q <- unique(object$info$q)
  # p <- unique(object$info$p)
  # pool$alpha <- rep(NA,times=q)
  # for(i in seq_len(q)){
  #   pool$alpha[i] <-  meta$alpha[i] + sum(meta$beta[,i] * base$alpha)
  # }
  # pool$beta <- matrix(NA,nrow=p,ncol=q)
  # for(i in seq_len(p)){
  #   for(j in seq_len(q)){
  #     pool$beta[i,j] <-  sum(meta$beta[,j] * base$beta[i,])
  #   }
  # }
  
  return(pool)
}

#' @export
#' @importFrom stats weights
#' @title
#' Extract Weights
#'
#' @description
#' Extracts coefficients from the meta learner,
#' i.e. the weights for the base learners.
#' 
#' @param object
#' \link[joinet]{joinet} object
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @examples
#' n <- 30; q <- 2; p <- 20
#' Y <- matrix(rnorm(n*q),nrow=n,ncol=q)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- joinet(Y=Y,X=X)
#' weights(object)
#' 
weights.joinet <- function(object,...){
  if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
  x <- object$meta
  coef <- lapply(object$meta,function(x) glmnet::coef.glmnet(object=x,s=x$lambda.min))
  coef <- do.call(what="cbind",args=coef)
  coef <- as.matrix(coef)
  colnames(coef) <- names(object$meta)
  return(coef)
}

print.joinet <- function(x,...){
  cat(paste0("joinet object"),"\n")
}

#--- Manuscript functions ------------------------------------------------------

#' @export
#' @title
#' Model comparison
#'
#' @description
#' Compares univariate and multivariate regression
#' 
#' @inheritParams joinet
#' 
#' @param nfolds.ext
#' number of external folds
#' 
#' @param nfolds.int
#' number of internal folds
#' 
#' @param foldid.ext
#' external fold identifiers\strong{:}
#' vector of length \eqn{n} with entries
#' between \eqn{1} and \code{nfolds.ext};
#' or \code{NULL}
#' 
#' @param foldid.int
#' internal fold identifiers\strong{:}
#' vector of length \eqn{n} with entries
#' between \eqn{1} and \code{nfolds.int};
#' or \code{NULL}
#' 
#' @param mnorm,spls,sier,mrce
#' experimental arguments\strong{:}
#' logical (install packages \code{spls}, \code{SiER}, or \code{MRCE})
#' 
#' @param ...
#' further arguments passed to \code{\link[glmnet]{glmnet}} and \code{\link[glmnet]{cv.glmnet}}
#' 
#' @examples
#' n <- 40; q <- 2; p <- 20
#' Y <- matrix(rnorm(n*q),nrow=n,ncol=q)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' cv.joinet(Y=Y,X=X)
#' 
cv.joinet <- function(Y,X,family="gaussian",nfolds.ext=5,nfolds.int=10,foldid.ext=NULL,foldid.int=NULL,type.measure="deviance",alpha.base=1,alpha.meta=0,mnorm=FALSE,spls=FALSE,sier=FALSE,mrce=FALSE,...){
  
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  #--- fold identifiers ---
  if(is.null(foldid.ext)){
    foldid.ext <- palasso:::.folds(y=Y[,1],nfolds=nfolds.ext) # temporary Y[,1]
  } else {
    nfolds.ext <- length(unique(foldid.ext))
  }
  
  #--- family ---
  if(length(family)==1){
    family <- rep(family,times=q)
  } else if(length(family)!=q){
    stop("Invalid argument family",call.=FALSE)
  }
  
  #--- cross-validated predictions ---
  
  models <- c("base","meta","spls"[spls],"mnorm"[mnorm],"sier"[sier],"mrce"[mrce],"none")
  pred <- lapply(X=models,function(x) matrix(data=NA,nrow=n,ncol=q))
  names(pred) <- models
  
  for(i in seq_len(nfolds.ext)){
    
    Y0 <- Y[foldid.ext!=i,,drop=FALSE]
    Y1 <- Y[foldid.ext==i,,drop=FALSE]
    X0 <- X[foldid.ext!=i,,drop=FALSE]
    X1 <- X[foldid.ext==i,,drop=FALSE]
    if(is.null(foldid.int)){
      foldid <- palasso:::.folds(y=Y0[,1],nfolds=nfolds.int) # temporary Y0[,1]
    } else {
      foldid <- foldid.int[foldid.ext!=i]
    }
    
    # base and meta learners
    fit <- joinet(Y=Y0,X=X0,family=family,type.measure=type.measure,alpha.base=alpha.base,alpha.meta=alpha.meta,foldid=foldid) # add ,...
    temp <- predict.joinet(fit,newx=X1)
    pred$base[foldid.ext==i,] <- temp$base
    pred$meta[foldid.ext==i,] <- temp$meta
    
    # other learners
    cond <- apply(X0,2,function(x) stats::sd(x)!=0)
    x0 <- X0[,cond]
    x1 <- X1[,cond]
    y0 <- apply(X=Y0,MARGIN=2,FUN=function(x) ifelse(is.na(x),sample(x[!is.na(x)],size=1),x))
    all(Y0==y0,na.rm=TRUE)
    
    if(mnorm){
      net <- glmnet::cv.glmnet(x=X0,y=y0,family="mgaussian",foldid=foldid,...) # ellipsis
      pred$mnorm[foldid.ext==i,] <- glmnet::predict.cv.glmnet(object=net,newx=X1,s="lambda.min",type="response")
    }
    if(spls){
      cv.spls <- spls::cv.spls(x=x0,y=y0,fold=nfolds.int,K=seq_len(10),
                               eta=seq(from=0.1,to=0.9,by=0.1),scale.x=FALSE,plot.it=FALSE)
      mspls <- spls::spls(x=x0,y=y0,K=cv.spls$K.opt,cv.spls$eta.opt,scale.x=FALSE)
      pred$spls[foldid.ext==i,] <- spls::predict.spls(object=mspls,newx=x1,type="fit")
    }
    if(sier){
      object <- SiER::cv.SiER(X=X0,Y=y0,K.cv=10)
      pred$sier[foldid.ext==i,] <- SiER::pred.SiER(cv.fit=object,X.new=X1)
    }
    if(mrce){
      lam1 <- rev(10^seq(from=-2,to=0,by=0.5))
      lam2 <- rev(10^seq(from=-2,to=0,by=0.5))
      object <- MRCE::mrce(X=x0,Y=y0,lam1.vec=lam1,lam2.vec=lam2,method="cv")
      pred$mrce[foldid.ext==i,] <- object$muhat + x1 %*% object$Bhat
    }
    
    pred$none[foldid.ext==i,] <- matrix(colMeans(Y0,na.rm=TRUE),nrow=sum(foldid.ext==i),ncol=ncol(Y0),byrow=TRUE)
    
  }
  
  #--- cross-validated loss ---
  loss <- matrix(data=NA,nrow=length(models),ncol=ncol(Y),
                 dimnames=list(models,colnames(Y)))
  for(j in models){
    for(i in seq_len(q)){
      cond <- !is.na(Y[,i])
      loss[j,i] <- palasso:::.loss(y=Y[cond,i],fit=pred[[j]][cond,i],family=family[i],type.measure=type.measure)[[1]]
    }
  }
  
  #--- model refit ---
  #fit <- joinet(Y=Y,X=X,family=family,type.measure=type.measure,alpha.base=alpha.base,alpha.meta=alpha.meta) # add ,...
  #list <- list(loss=loss,fit=fit)
  
  return(loss)
}
