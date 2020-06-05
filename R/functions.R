
# import unexported functions:
# FUNCTION <- get("FUNCTION",envir=asNamespace("PACKAGE"))

#--- Main function -------------------------------------------------------------

#' @export
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
#' Armin Rauschenberger, Enrico Glaab (2019)
#' "joinet: predicting correlated outcomes jointly
#' to improve clinical prognosis"
#' \emph{Manuscript in preparation}.
#' 
#' @details
#' \strong{correlation:}
#' The \eqn{q} outcomes should be positively correlated.
#' Avoid negative correlations by changing the sign of the variable.
#' 
#' \strong{elastic net:}
#' \code{alpha.base} controls input-output effects,
#' \code{alpha.meta} controls output-output effects;
#' lasso renders sparse models (\code{alpha}\eqn{=1}),
#' ridge renders dense models (\code{alpha}\eqn{=0})
#' 
#' @return
#' This function returns an object of class \code{joinet}.
#' Available methods include
#' \code{\link[=predict.joinet]{predict}},
#' \code{\link[=coef.joinet]{coef}},
#' and \code{\link[=weights.joinet]{weights}}.
#' The slots \code{base} and \code{meta} each contain
#' \eqn{q} \code{\link[glmnet]{cv.glmnet}}-like objects.
#' 
#' @seealso
#' \code{\link{cv.joinet}}, vignette
#' 
#' @examples
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' object <- joinet(Y=Y,X=X)
#' 
#' \dontrun{
#' browseVignettes("joinet") # further examples}
#' 
joinet <- function(Y,X,family="gaussian",nfolds=10,foldid=NULL,type.measure="deviance",alpha.base=1,alpha.meta=0,...){
  
  #--- temporary ---
  # family <- "gaussian"; nfolds <- 10; foldid <- NULL; type.measure <- "deviance"
  
  #--- checks ---
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  
  cornet:::.check(x=Y,type="matrix",miss=TRUE)
  if(any(stats::cor(Y,use="pairwise.complete.obs")<0,na.rm=TRUE)){warning("Negative correlation!",call.=FALSE)}
  cornet:::.check(x=X,type="matrix")
  #cornet:::.check(x=family,type="vector",values=c("gaussian","binomial","poisson"))
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
    class(base[[i]]) <- "cv.glmnet" # trial 2020-01-10
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
  if(family %in% c("gaussian","cox")){
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
  if(family %in% c("gaussian","cox")){
    return(x)
  } else if(family=="binomial"){
    if(any(x<0|x>1)){stop("Invalid!",call.=FALSE)}
    return(log(x/(1-x)))
  } else if(family=="poisson"){
    if(any(x<0)){stop("Invalid!",call.=FALSE)}
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
#' @return 
#' This function returns predictions from base and meta learners.
#' The slots \code{base} and \code{meta} each contain a matrix
#' with \eqn{n} rows (samples) and \eqn{q} columns (variables).
#' 
#' @examples
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' Y[,1] <- 1*(Y[,1]>median(Y[,1]))
#' object <- joinet(Y=Y,X=X,family=c("binomial","gaussian","gaussian"))
#' predict(object,newx=X)
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
    base[,i] <- as.numeric(stats::predict(object=x$base[[i]]$glmnet.fit,newx=newx,s=x$base[[i]]$lambda.min,type="link"))
    #base[,i] <- as.numeric(glmnet:::predict.cv.glmnet(object=x$base[[i]],newx=newx,s="lambda.min",type="link"))
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
#' @return
#' This function returns the pooled coefficients.
#' The slot \code{alpha} contains the intercepts
#' in a vector of length \eqn{q},
#' and the slot \code{beta} contains the slopes
#' in a matrix with \eqn{p} rows (inputs) and \eqn{q} columns.
#' 
#' @examples
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' object <- joinet(Y=Y,X=X)
#' coef <- coef(object)
#' 
coef.joinet <- function(object,...){
  if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
  
  # base coefficients
  base <- list()
  coef <- sapply(object$base,function(x) stats::coef(object=x$glmnet.fit,s=x$lambda.min))
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
#' @return
#' This function returns a matrix with
#' \eqn{1+q} rows and \eqn{q} columns.
#' The first row contains the intercepts,
#' and the other rows contain the slopes,
#' which are the effects of the outcomes
#' in the row on the outcomes in the column.
#' 
#' @examples
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' object <- joinet(Y=Y,X=X)
#' weights(object)
#' 
weights.joinet <- function(object,...){
  if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
  x <- object$meta
  coef <- lapply(object$meta,function(x) stats::coef(object=x,s=x$lambda.min))
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
#' Compares univariate and multivariate regression.
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
#' @param mnorm,spls,mrce,sier,mtps,rmtl,gpm
#' experimental arguments\strong{:}
#' logical
#' (\code{TRUE} requires packages \code{spls}, \code{MRCE}, \code{SiER}, \code{MTPS}, \code{RMTL} or \code{GPM})
#' 
#' @param mice
#' missing data imputation\strong{:}
#' logical (\code{mice=TRUE} requires package \code{mice})
#' 
#' @param cvpred
#' return cross-validated predicitions: logical
#' 
#' @param ...
#' further arguments passed to \code{\link[glmnet]{glmnet}}
#' and \code{\link[glmnet]{cv.glmnet}}
#' 
#' @return 
#' This function returns a matrix with \eqn{q} columns,
#' including the cross-validated loss from the univariate models
#' (\code{base}), the multivariate models (\code{meta}),
#' and the intercept-only models (\code{none}).
#' 
#' @examples
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' cv.joinet(Y=Y,X=X)
#' 
#' \dontrun{
#' # correlated features
#' n <- 50; p <- 100; q <- 3
#' mu <- rep(0,times=p)
#' Sigma <- 0.90^abs(col(diag(p))-row(diag(p)))
#' X <- MASS::mvrnorm(n=n,mu=mu,Sigma=Sigma)
#' mu <- rowSums(X[,sample(seq_len(p),size=5)])
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=mu))
#' #Y <- t(MASS::mvrnorm(n=q,mu=mu,Sigma=diag(n)))
#' cv.joinet(Y=Y,X=X)}
#' 
#' \dontrun{
#' # other distributions
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' eta <- rowSums(X[,1:5])
#' Y <- replicate(n=q,expr=rbinom(n=n,size=1,prob=1/(1+exp(-eta))))
#' cv.joinet(Y=Y,X=X,family="binomial")
#' Y <- replicate(n=q,expr=rpois(n=n,lambda=exp(scale(eta))))
#' cv.joinet(Y=Y,X=X,family="poisson")}
#' 
#' \dontrun{
#' # uncorrelated outcomes
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' y <- rnorm(n=n,mean=rowSums(X[,1:5]))
#' Y <- cbind(y,matrix(rnorm(n*(q-1)),nrow=n,ncol=q-1))
#' cv.joinet(Y=Y,X=X)}
#' 
#' \dontrun{
#' # sparse and dense models
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' set.seed(1) # fix folds
#' cv.joinet(Y=Y,X=X,alpha.base=1) # lasso
#' set.seed(1)
#' cv.joinet(Y=Y,X=X,alpha.base=0) # ridge}
#' 
cv.joinet <- function(Y,X,family="gaussian",nfolds.ext=5,nfolds.int=10,foldid.ext=NULL,foldid.int=NULL,type.measure="deviance",alpha.base=1,alpha.meta=0,mnorm=FALSE,spls=FALSE,mrce=FALSE,sier=FALSE,mtps=FALSE,rmtl=FALSE,gpm=FALSE,mice=FALSE,cvpred=FALSE,...){
  
  # family <- "gaussian"; nfolds.ext <- 5; nfolds.int <- 10; foldid.ext <- foldid.int <- NULL; type.measure <- "deviance"; alpha.base <- 1; alpha.meta <- 0; mnorm <- spls <- mrce <- sier <- mtps <- rmtl <- gpm <- mice <- cvpred <- TRUE
  # family <- "binomial"; nfolds.ext <- 1; foldid.ext <- fold; nfolds.int <- 10; foldid.int <- NULL; type.measure <- "deviance"; alpha.base <- alpha.meta <- 1; mnorm <- spls <- mrce <- sier <- mtps <- rmtl <- gpm <- mice <- cvpre <- TRUE
  
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  #--- fold identifiers ---
  if(is.null(foldid.ext)){
    foldid.ext <- palasso:::.folds(y=Y[,1],nfolds=nfolds.ext) # temporary Y[,1]
  } else {
    #nfolds.ext <- length(unique(foldid.ext))
    nfolds.ext <- max(foldid.ext)
  }
  
  #--- family ---
  if(length(family)==1){
    family <- rep(family,times=q)
  } else if(length(family)!=q){
    stop("Invalid argument family",call.=FALSE)
  }
  
  #--- checks ---
  if(any(mnorm,spls,mrce,sier,gpm) & any(family!="gaussian")){
    stop("\"mnorm\", \"spls\", \"mrce\" and \"sier\" require family \"gaussian\".",call.=FALSE)
  }
  if(any(mtps,rmtl) & any(!family %in% c("gaussian","binomial"))){
    stop("\"mtps\" and \"rmtl\" require family \"gaussian\" or \"binomial\".",call.=FALSE)
  }
  
  #--- cross-validated predictions ---
  
  models <- c("base","meta","mnorm"[mnorm],"spls"[spls],"mrce"[mrce],"sier"[sier],"mtps"[mtps],"rmtl"[rmtl],"gpm"[gpm],"none")
  pred <- lapply(X=models,function(x) matrix(data=NA,nrow=n,ncol=q))
  names(pred) <- models
  
  for(i in seq_len(nfolds.ext)){
    
    Y0 <- Y[foldid.ext!=i,,drop=FALSE]
    #Y1 <- Y[foldid.ext==i,,drop=FALSE]
    X0 <- X[foldid.ext!=i,,drop=FALSE]
    X1 <- X[foldid.ext==i,,drop=FALSE]
    # Remove constant features and impute missing values here?
    if(is.null(foldid.int)){
      foldid <- palasso:::.folds(y=Y0[,1],nfolds=nfolds.int) # temporary Y0[,1]
    } else {
      foldid <- foldid.int[foldid.ext!=i]
    }
    
    # base and meta learners
    fit <- joinet(Y=Y0,X=X0,family=family,type.measure=type.measure,alpha.base=alpha.base,alpha.meta=alpha.meta,foldid=foldid) # add ellipsis (...)
    temp <- predict.joinet(fit,newx=X1)
    pred$base[foldid.ext==i,] <- temp$base
    pred$meta[foldid.ext==i,] <- temp$meta
    
    # constant features, missing values
    # Consider moving this upwards!
    cond <- apply(X0,2,function(x) stats::sd(x)!=0)
    x0 <- X0[,cond]
    x1 <- X1[,cond]
    if(mice & any(is.na(Y0))){
      if(requireNamespace("mice",quietly=TRUE)){
        y0 <- as.matrix(mice::complete(data=mice::mice(Y0,m=1,method="pmm",seed=1,printFlag=FALSE),action="all")[[1]])
      } else {
        stop("mice=TRUE requires install.packages(\"mice\").",call.=FALSE)
      }
    } else {
      #y0 <- apply(X=Y0,MARGIN=2,FUN=function(x) ifelse(is.na(x),sample(x[!is.na(x)],size=1),x))
      y0 <- apply(X=Y0,MARGIN=2,FUN=function(x) ifelse(is.na(x),stats::median(x[!is.na(x)]),x))
    }
    all(Y0==y0,na.rm=TRUE)
    
    # other learners
    if(mnorm){
      net <- glmnet::cv.glmnet(x=X0,y=y0,family="mgaussian",foldid=foldid,alpha=alpha.base) # add ellipsis (...)
      pred$mnorm[foldid.ext==i,] <- stats::predict(object=net,newx=X1,s="lambda.min",type="response")
    }
    if(spls){
      cv.spls <- spls::cv.spls(x=x0,y=y0,fold=nfolds.int,K=seq_len(10),
                               eta=seq(from=0.1,to=0.9,by=0.1),plot.it=FALSE)
      spls <- spls::spls(x=x0,y=y0,K=cv.spls$K.opt,eta=cv.spls$eta.opt)
      pred$spls[foldid.ext==i,] <- spls::predict.spls(object=spls,newx=x1,type="fit")
    }
    if(mrce){
      stop("MRCE not yet implemented",call.=FALSE) # bug?
      lam1 <- rev(10^seq(from=-2,to=0,by=0.5))
      lam2 <- rev(10^seq(from=-2,to=0,by=0.5))
      #lam1 <- lam2 <- 10^seq(from=0,to=-0.7,length.out=5)
      #lam2 <- rev(10^seq(from=0,to=2,by=0.5))
      #lam2 <- c(2,1,0.5)
      #lam1 <- 10^seq(from=0,to=-5,length.out=11)
      #lam2 <- seq(from=1,to=0.1,by=-0.1)
      object <- MRCE::mrce(X=x0,Y=y0,lam1.vec=lam1,lam2.vec=lam2,method="cv",silent=FALSE,cov.tol=0.1,tol.out=1e-10)
      pred$mrce[foldid.ext==i,] <- matrix(object$muhat,nrow=nrow(x1),ncol=q,byrow=TRUE) + x1 %*% object$Bhat
    }
    if(sier){
      stop("SiER not yet implemented",call.=FALSE) # slow!
      object <- SiER::cv.SiER(X=X0,Y=y0,K.cv=5,thres=0.5)
      pred$sier[foldid.ext==i,] <- SiER::pred.SiER(cv.fit=object,X.new=X1)
    }
    if(mtps){
      if(alpha.base==0){
        step1 <- MTPS::glmnet.ridge
      } else if(alpha.base==1){
        step1 <- MTPS::glmnet.lasso
      } else {
        stop("MTPS requires alpha.base 0 or 1.",call.=FALSE)
      }
      step2 <- MTPS::rpart1
      object <- MTPS::MTPS(xmat=X0,ymat=y0,family=family,nfold=nfolds.int,
                           method.step1=step1,method.step2=step2)
      # nfold has no effect for cv=FALSE (default)
      pred$mtps[foldid.ext==i,] <- MTPS::predict.MTPS(object=object,newdata=X1,type="response")
    }
    if(rmtl){
      if(all(family=="gaussian")){
        type <- "Regression"
      } else if(all(family=="binomial")){
        type <- "Classification"
      } else {
        stop("RMTL requires either \"gaussian\" or \"binomial\".",call.=FALSE)
      }
      Y0l <- lapply(seq_len(ncol(Y0)),function(i) 2*Y0[,i,drop=FALSE]-1)
      X0l <- lapply(seq_len(ncol(Y0)),function(i) X0)
      X1l <- lapply(seq_len(ncol(Y0)),function(i) X1)
      #Lam2 <- 10^seq(from=0,to=-5,length.out=51)
      #cvm <- sapply(Lam2,function(x) min(RMTL::cvMTL(X=X0l,Y=Y0l,type=type,Lam2=x)$cvm))
      #Lam2 <- Lam2[which.min(cvm)]
      # Manually tune regularisation parameter lambda2!
      cvMTL <- RMTL::cvMTL(X=X0l,Y=Y0l,type=type) #,Lam2=Lam2)
      MTL <- RMTL::MTL(X=X0l,Y=Y0l,type=type,Lam1=cvMTL$Lam1.min)
      temp <- RMTL:::predict.MTL(object=MTL,newX=X1l)
      pred$rmtl[foldid.ext==i,] <- do.call(what="cbind",args=temp)
    }
    if(gpm){
      if(any(family!="gaussian")){
        stop("GPM requires \"gaussian\" family.",call.=FALSE)
      }
      object <- GPM::Fit(X=X0,Y=Y0)
      pred$gpm[foldid.ext==i,] <- GPM::Predict(XF=X1,Model=object)$YF
    }
    
    if(FALSE){ # bgsmtr (for SNPs only? error!)
      temp <- bgsmtr::bgsmtr(X=t(X0),Y=t(Y0),group=rep(1,times=ncol(X0)))
      # check dimensions of example, use group 1:p instead of rep(1,p)?
    }
    
    if(FALSE){ # MGLM (for multinomial data only?)
      MGLM::MGLMsparsereg.fit()
      MGLM::MGLMtune()
    }
    
    pred$none[foldid.ext==i,] <- matrix(colMeans(Y0,na.rm=TRUE),nrow=sum(foldid.ext==i),ncol=ncol(Y0),byrow=TRUE)
    
  }
  
  #--- cross-validated loss ---
  loss <- matrix(data=NA,nrow=length(models),ncol=ncol(Y),
                 dimnames=list(models,colnames(Y)))
  for(j in models){
    for(i in seq_len(q)){
      cond <- !is.na(Y[,i]) & foldid.ext!=0  # added foldid.ext!=0
      loss[j,i] <- palasso:::.loss(y=Y[cond,i],fit=pred[[j]][cond,i],family=family[i],type.measure=type.measure)[[1]]
    }
  }
  
  #--- model refit ---
  #fit <- joinet(Y=Y,X=X,family=family,type.measure=type.measure,alpha.base=alpha.base,alpha.meta=alpha.meta) # add ,...
  #list <- list(loss=loss,fit=fit)
  
  if(cvpred){
    loss <- list(loss=loss,cvpred=pred$meta)
  }
  
  return(loss)
}



plot.matrix <- function (X, margin = 0, labels = TRUE, las = 1, cex = 1, range = NULL, cutoff = 0, digits=2) {
  #margin <- 0; labels <- TRUE; las <- 1; cex <- 1; range <- NULL; cutoff <- 0; digits <- 2
  
  if(is.vector(X)){X <- as.matrix(X,ncol=1)}
  
  n <- nrow(X)
  p <- ncol(X)
  if(is.null(rownames(X))&n!=1){rownames(X) <- seq_len(n)}
  if(is.null(colnames(X))&p!=1){colnames(X) <- seq_len(p)}
  v <- ifelse(n==1,1,0.5/(n - 1))
  h <- ifelse(p==1,1,0.5/(p - 1))
  
  graphics::plot.new()
  graphics::plot.window(xlim = c(-h, 1 + h), ylim = c(-v, 1 + v))
  par_usr <- graphics::par()$usr
  graphics::par(usr = c(-h, 1 + h, -v, 1 + v))
  
  at <- (seq(nrow(X)) - 1)/(nrow(X) - 1)
  graphics::mtext(text = rev(rownames(X)), at = at, side = 2, las = las, cex = cex,line=0.1)
  at <- (seq(ncol(X)) - 1)/(ncol(X) - 1)
  graphics::mtext(text = colnames(X), at = at , side = 3, las = las, cex = cex,line=0.1)
  
  if(is.null(range)){range <- range(X,-X,na.rm=TRUE)}
  if(any(X<min(range),na.rm=TRUE)){stop("Invalid.")}
  if(any(X>max(range),na.rm=TRUE)){stop("Invalid.")}
  breaks <- c(seq(from=min(range),to=cutoff,length.out=101),
              seq(from=cutoff,to=max(range),length.out=101)[-1])
  col <- c(grDevices::colorRampPalette(c("darkblue","blue","white"))(100),
           grDevices::colorRampPalette(c("white","red","darkred"))(100))

  
  image <- t(X)[, seq(from = n, to = 1, by = -1),drop=FALSE]
  graphics::image(x = image, breaks=breaks, col=col, add = TRUE)
  
  if (any(margin==1)){
    graphics::segments(x0 = -h,
                       x1 = 1 + h,
                       y0 = seq(from = -v, to = 1 + v, by = 2 * v),
                       col = "white",
                       lwd = 3)
  }
  if (any(margin==2)){
    graphics::segments(x0 = seq(from = -h, to = 1 + h, by = 2 * h),
                       y0 = 1 + v,
                       y1 = 0 - v,
                       col = "white",
                       lwd = 3)
  }
  if (labels) {
    labels <- round(as.numeric(X), digits = digits)
    is.na <- is.na(labels)
    labels <- format(labels, digits = digits)
    labels[is.na] <- ""
    xs <- rep(seq_len(p), each = n)
    ys <- rep(seq_len(n), times = p)
    if(p==1){x <- 0.5}else{x <- (xs - 1)/(p - 1)}
    if(n==1){y <- 0.5}else{y <- (n - ys)/(n - 1)}
    graphics::text(x = x, y = y, labels = labels, col = "black",cex=cex)
  }
  graphics::par(usr = par_usr)
}
