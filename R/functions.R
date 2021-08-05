
#--- import unexported functions:
# FUNCTION <- get("FUNCTION",envir=asNamespace("PACKAGE"))

#--- deactivate on Solaris:
# if(!grepl('SunOS',Sys.info()['sysname'])){}

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
#' and \eqn{q} columns (outputs)
#' 
#' @param X
#' inputs\strong{:}
#' numeric matrix with \eqn{n} rows (samples)
#' and \eqn{p} columns (inputs)
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
#' elastic net mixing parameter for meta learners\strong{:}
#' numeric between \eqn{0} (ridge) and \eqn{1} (lasso)
#' 
#' @param weight
#' input-output relations\strong{:}
#' matrix with \eqn{p} rows (inputs) and \eqn{q} columns (outputs)
#' with entries \eqn{0} (exclude) and \eqn{1} (include),
#' or \code{NULL} (see details)
#' 
#' @param sign
#' output-output relations\strong{:}
#' matrix with \eqn{q} rows ("meta-inputs") and \eqn{q} columns (outputs), 
#' with entries \eqn{-1} (negative), \eqn{0} (none),
#' \eqn{1} (positive) and \eqn{NA} (any),
#' or \code{NULL} (see details)
#' 
#' @param ...
#' further arguments passed to \code{\link[glmnet]{glmnet}}
#' 
#' @references 
#' Armin Rauschenberger, Enrico Glaab (2021)
#' "Predicting correlated outcomes from molecular data"
#' \emph{Bioinformatics}. btab576
#' \url{https://doi.org/10.1093/bioinformatics/btab576}
#' 
#' @details
#' \strong{input-output relations:}
#' In this matrix with \eqn{p} rows and \eqn{q} columns,
#' the entry in the \eqn{j}th row and the \eqn{k}th column
#' indicates whether the \eqn{j}th input may be used for 
#' modelling the \eqn{k}th output
#' (where \eqn{0} means "exclude" and
#' \eqn{1} means "include").
#' By default (\code{sign=NULL}),
#' all entries are set to \eqn{1}.
#' 
#' \strong{output-output relations:}
#' In this matrix with \eqn{q} rows and \eqn{q} columns,
#' the entry in the \eqn{l}th row and the \eqn{k}th column
#' indicates how the \eqn{l}th output may be used for
#' modelling the \eqn{k}th output
#' (where \eqn{-1} means negative effect,
#' \eqn{0} means no effect,
#' \eqn{1} means positive effect,
#' and \eqn{NA} means any effect).
#' 
#' There are three short-cuts for filling up this matrix:
#' (1) \code{sign=1} sets all entries to \eqn{1} (non-negativity constraints).
#' This is useful if all pairs of outcomes
#' are assumed to be \emph{positively} correlated
#' (potentially after changing the sign of some outcomes).
#' (2) \code{code=NA} sets all diagonal entries to \eqn{1}
#' and all off-diagonal entries to \code{NA} (no constraints).
#' (3) \code{sign=NULL} uses Spearman correlation to determine the entries,
#' with \eqn{-1} for significant negative, \eqn{0} for insignificant,
#' \eqn{1} for significant positive correlations.
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
#' \dontshow{
#' if(!grepl('SunOS',Sys.info()['sysname'])){
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' object <- joinet(Y=Y,X=X)}}
#' \dontrun{
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' object <- joinet(Y=Y,X=X)}
#' 
#' \dontrun{
#' browseVignettes("joinet") # further examples}
#' 
joinet <- function(Y,X,family="gaussian",nfolds=10,foldid=NULL,type.measure="deviance",alpha.base=1,alpha.meta=1,weight=NULL,sign=NULL,...){
  # VERIFY CODE FOR WEIGHT AND SIGN!
  
  #--- temporary ---
  # family <- "gaussian"; nfolds <- 10; foldid <- NULL; type.measure <- "deviance"; alpha.base <- alpha.meta <- 1; weight <- NULL; sign <- NULL
  
  #--- checks ---
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  
  cornet:::.check(x=Y,type="matrix",miss=TRUE)
  #if(constraint){
  #  for(i in 1:(ncol(Y)-1)){
  #    for(j in i:ncol(Y)){
  #      cor <- stats::cor.test(Y[,i],Y[,j],use="pairwise.complete.obs")
  #      if(cor$estimate<0 & cor$p.value<0.05){
  #        warning(paste("Columns",i,"and",j,"are negatively correlated. Consider using constraint=FALSE."),call.=FALSE)
  #      }
  #    }
  #  }
  #}
  #if(any(stats::cor(Y,use="pairwise.complete.obs")<0,na.rm=TRUE)){warning("Negative correlation!",call.=FALSE)}
  cornet:::.check(x=X,type="matrix")
  #cornet:::.check(x=family,type="vector",values=c("gaussian","binomial","poisson"))
  if(nrow(Y)!=nrow(X)){stop("Contradictory sample size.",call.=FALSE)}
  n <- nrow(Y); q <- ncol(Y); p <- ncol(X)
  
  cornet:::.check(x=nfolds,type="scalar",min=3)
  cornet:::.check(x=foldid,type="vector",values=seq_len(nfolds),null=TRUE)
  cornet:::.check(x=type.measure,type="string",values=c("deviance","class","mse","mae")) # not auc (min/max confusion)
  cornet:::.check(x=alpha.base,type="scalar",min=0,max=1)
  cornet:::.check(x=alpha.meta,type="scalar",min=0,max=1)
  cornet:::.check(x=weight,type="matrix",min=0,max=1,null=TRUE,dim=c(p,q))
  if(!is.null(c(list(...)$lower.limits,list(...)$upper.limits))){
    stop("Reserved arguments \"lower.limits\" and \"upper.limits\".",call.=FALSE)
  }
  
  if(is.null(weight)){
    pf <- matrix(1,nrow=p,ncol=q)
  } else {
    pf <- 1/weight
  }
  
  #--- constraints --- 
  null <- is.null(sign)
  if(null){sign <- 0}
  if(!is.matrix(sign)){
    sign <- matrix(sign,nrow=q,ncol=q)
    diag(sign) <- 1
  }
  if(any(diag(sign)!=1)){
    warning("Matrix \"sign\" has entries other than one on the diagonal.",call.=FALSE)
  }
  temp <- sign[lower.tri(sign)|upper.tri(sign)]
  if(!null & all(!is.na(temp)&temp==0)){
    warning("Matrix \"sign\" only has zeros off the diagonal.",call.=FALSE)
  }
  for(i in seq(from=1,to=q,by=1)){
    for(j in seq(from=i,to=q,by=1)){
      cor <- stats::cor.test(Y[,i],Y[,j],use="pairwise.complete.obs",method="spearman",exact=FALSE)
      if(cor$p.value>0.05){next}
      if(null){
        sign[i,j] <- sign[j,i] <- sign(cor$estimate)
      } else if(!is.na(sign[i,j]) & sign[i,j]*sign(cor$estimate)==-1){
        warning(paste("Outputs",i,"and",j,"have a significant",ifelse(sign(cor$estimate)==1,"*positive*","*negative*"),"correlation. Consider changing argument \"sign\" (e.g. sign=NA or sign=NULL)."),call.=FALSE)
      }
    }
  }

  # # long version
  # 
  # if(is.null(sign)){
  #   sign <- matrix(0,nrow=q,ncol=q)
  #   for(i in seq(from=1,to=q,by=1)){
  #     for(j in seq(from=i,to=q,by=1)){
  #       cor <- stats::cor.test(Y[,i],Y[,j],use="pairwise.complete.obs",method="kendall")
  #       if(cor$p.value<0.05){
  #         sign[i,j] <- sign[j,i] <- sign(cor$estimate)
  #       }
  #     }
  #   }
  # } else {
  #   if(!is.matrix(sign)){
  #     sign <- matrix(sign,nrow=q,ncol=q) # accept 1 and NA
  #   }
  #   for(i in seq(from=1,to=q,by=1)){
  #     for(j in seq(from=i,to=q,by=1)){
  #       cor <- stats::cor.test(Y[,i],Y[,j],use="pairwise.complete.obs",method="kendall")
  #       if(!is.na(sign[i,j]) & cor$p.value<0.05 & sign[i,j]*sign(cor$estimate)==-1){
  #         warning(paste("Outputs",i,"and",j,"have an unexpected significant correlation. Consider using sign=NULL."),call.=FALSE)
  #       }
  #     }
  #   }
  # }
  
  
  ## old snippets
  
  #if(is.null(sign)){
  #  sign <- 0   # 0 or NA or 1
  #  check <- FALSE
  #} else {
  #  check <- TRUE
  #}
  #
  #if(is.matrix(sign)){
  #  # perform checks
  #} if else(!is.matrix(sign) & !is.null(sign)){
  #  sign <- matrix(sign,nrow=q,ncol=q)
  #  # perform checks
  #} else {  
  #  if(is.null(sign)){sign <- 0}  # 0 or NA or 1
  #  }
  
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
    base[[i]]$glmnet.fit <- glmnet::glmnet(y=Y[cond,i],x=X[cond,],family=family[i],alpha=alpha.base,penalty.factor=pf[,i],...) # ellipsis
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
      object <- glmnet::glmnet(y=Y0[cond,i],x=X0[cond,],family=family[i],alpha=alpha.base,penalty.factor=pf[,i],...) # ellipsis
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
    # trial start
    lower.limits <- rep(-Inf,times=q)
    upper.limits <- rep(Inf,times=q)
    lower.limits[sign[,i]>=0] <- 0
    upper.limits[sign[,i]<=0] <- 0
    # trial end
    cond <- !is.na(Y[,i])
    meta[[i]] <- glmnet::cv.glmnet(y=Y[cond,i],x=hat[cond,],
                                   lower.limits=lower.limits, # was first lower.limits=0 and later ifelse(constraint,0,-Inf) 
                                   upper.limits=upper.limits, # was first upper.limits=Inf
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
#' \dontshow{
#' if(!grepl('SunOS',Sys.info()['sysname'])){
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' Y[,1] <- 1*(Y[,1]>median(Y[,1]))
#' object <- joinet(Y=Y,X=X,family=c("binomial","gaussian","gaussian"))
#' predict(object,newx=X)}}
#' \dontrun{
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' Y[,1] <- 1*(Y[,1]>median(Y[,1]))
#' object <- joinet(Y=Y,X=X,family=c("binomial","gaussian","gaussian"))
#' predict(object,newx=X)}
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
#' \dontshow{
#' if(!grepl('SunOS',Sys.info()['sysname'])){
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' object <- joinet(Y=Y,X=X)
#' coef <- coef(object)}}
#' \dontrun{
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' object <- joinet(Y=Y,X=X)
#' coef <- coef(object)}
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
#' \dontshow{
#' if(!grepl('SunOS',Sys.info()['sysname'])){
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' object <- joinet(Y=Y,X=X)
#' weights(object)}}
#' \dontrun{
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' object <- joinet(Y=Y,X=X)
#' weights(object)}
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
#' @param compare
#' experimental arguments\strong{:}
#' character vector with entries "mnorm", "spls", "mrce",
#' "sier", "mtps", "rmtl", "gpm" and others
#' (requires packages \code{spls}, \code{MRCE}, \code{SiER}, \code{MTPS}, \code{RMTL} or \code{GPM})
#' 
#' @param mice
#' missing data imputation\strong{:}
#' logical (\code{mice=TRUE} requires package \code{mice})
#' 
#' @param cvpred
#' return cross-validated predictions: logical
#' 
#' @param times
#' measure computation time\strong{:}
#' logical
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
#' \dontshow{
#' if(!grepl('SunOS',Sys.info()['sysname'])){
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' cv.joinet(Y=Y,X=X)}}
#' \dontrun{
#' n <- 50; p <- 100; q <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- replicate(n=q,expr=rnorm(n=n,mean=rowSums(X[,1:5])))
#' cv.joinet(Y=Y,X=X)}
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
cv.joinet <- function(Y,X,family="gaussian",nfolds.ext=5,nfolds.int=10,foldid.ext=NULL,foldid.int=NULL,type.measure="deviance",alpha.base=1,alpha.meta=1,compare=FALSE,mice=FALSE,cvpred=FALSE,times=FALSE,...){
  
  if(FALSE){
  fold <- foldid.ext
  family <- "gaussian"; nfolds.ext <- 5; nfolds.int <- 10; foldid.ext <- foldid.int <- NULL; type.measure <- "deviance"; alpha.base <- alpha.meta <- 1; mice <- cvpred <- times <- FALSE
  foldid.ext <- fold; nfolds.ext <- 1
  #nfolds.ext <- 1; nfolds.int <- 10; foldid.int <- NULL; compare <- TRUE
  }
  
  if(length(compare)==1 && compare==TRUE){
    if(all(family=="gaussian")){
      compare <- c("mnorm","mars","spls","mrce","map","mrf","sier","mcen","gpm","rmtl","mtps")
    } else if(all(family=="binomial")){
      compare <- c("mars","mcen","rmtl","mtps")
    } else {
      stop("Comparison not implemented for mixed families.",call.=FALSE)
    }
  }
  
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
  
  # check packages
  pkgs <- .packages(all.available=TRUE)
  
  if(is.character(compare)){
    for(i in seq_along(compare)){
      pkg <- switch(compare[i],mnorm="glmnet",mars="earth",spls="spls",
                  mrce="MRCE",map="remMap",mrf="MultivariateRandomForest",
                  sier="SiER",mcen="mcen",gpm="GPM",rmtl="RMTL",mtps="MTPS",
                  stop("Invalid method.",call.=FALSE))
      if(!pkg %in% pkgs){
        stop("Method \"",compare[i],"\" requires package \"",pkg,"\".",call.=FALSE)
      }
    }
  }

  #--- checks ---
  #if(any( & any(family!="gaussian")){
  #  stop("\"mnorm\", \"spls\", \"mrce\" and \"sier\" require family \"gaussian\".",call.=FALSE)
  #}
  #if(any(mtps,rmtl) & any(!family %in% c("gaussian","binomial"))){
  #  stop("\"mtps\" and \"rmtl\" require family \"gaussian\" or \"binomial\".",call.=FALSE)
  #}
  
  #--- cross-validated predictions ---
  
  models <- c("base","meta",compare,"none")
  pred <- lapply(X=models,function(x) matrix(data=NA,nrow=n,ncol=q))
  time <- lapply(X=models,function(x) NA)
  names(pred) <- names(time) <- models
  
  for(i in seq_len(nfolds.ext)){
    
    Y0 <- Y[foldid.ext!=i,,drop=FALSE]
    #Y1 <- Y[foldid.ext==i,,drop=FALSE]
    X0 <- X[foldid.ext!=i,,drop=FALSE]
    X1 <- X[foldid.ext==i,,drop=FALSE]
    
    # standardise features (trial)
    #mu <- apply(X=X0,MARGIN=2,FUN=function(x) mean(x))
    #sd <- apply(X=X0,MARGIN=2,FUN=function(x) stats::sd(x))
    #X0 <- t((t(X0)-mu)/sd)[,sd!=0]
    #X1 <- t((t(X1)-mu)/sd)[,sd!=0]
    # or standardise once before cv?
    
    # remove constant features
    cond <- apply(X=X0,MARGIN=2,FUN=function(x) stats::sd(x)!=0)
    X0 <- X0[,cond]; X1 <- X1[,cond]
    
    if(is.null(foldid.int)){
      foldid <- palasso:::.folds(y=Y0[,1],nfolds=nfolds.int) # temporary Y0[,1]
    } else {
      foldid <- foldid.int[foldid.ext!=i]
    }
    
    # base and meta learners
    start <- Sys.time()
    fit <- joinet(Y=Y0,X=X0,family=family,type.measure=type.measure,alpha.base=alpha.base,alpha.meta=alpha.meta,foldid=foldid) # add ellipsis (...)
    # also do not standardise!
    temp <- predict.joinet(fit,newx=X1)
    pred$base[foldid.ext==i,] <- temp$base
    pred$meta[foldid.ext==i,] <- temp$meta
    end <- Sys.time()
    time$meta <- as.numeric(difftime(end,start,units="secs"))
    
    # missing values
    if(mice & any(is.na(Y0))){
      if(requireNamespace("mice",quietly=TRUE)){
        y0 <- as.matrix(mice::complete(data=mice::mice(Y0,m=1,method="pmm",seed=1,printFlag=FALSE),action="all")[[1]])
      } else {
        stop("Imputation by PMM requires install.packages(\"mice\").",call.=FALSE)
      }
    } else {
      y0 <- apply(X=Y0,MARGIN=2,FUN=function(x) ifelse(is.na(x),stats::median(x[!is.na(x)]),x))
    }
    all(Y0==y0,na.rm=TRUE)
    
    # other learners
    
    if("mnorm" %in% compare){
      cat("mnorm"," ")
      start <- Sys.time()
      if(any(family!="gaussian")){
        stop("mnorm requires \"gaussian\" family.",call.=FALSE)
      }
      net <- glmnet::cv.glmnet(x=X0,y=y0,family="mgaussian",foldid=foldid,alpha=alpha.base) # add ellipsis (...)
      pred$mnorm[foldid.ext==i,] <- stats::predict(object=net,newx=X1,s="lambda.min",type="response")
      end <- Sys.time()
      time$mnorm <- as.numeric(difftime(end,start,units="secs"))
    } else {
      net <- glmnet::glmnet(x=X0,y=y0,family="mgaussian")
    }
    
    if("mars" %in% compare){
      cat("mars"," ")
      start <- Sys.time()
      if(all(family=="gaussian")){
        object <- earth::earth(x=X0,y=y0)
        # equivalent: object <- mda::mars(x=X0,y=y0)
      } else if(all(family=="binomial")){
        object <- earth::earth(x=X0,y=y0,glm=list(family=stats::binomial))
      } else {
        stop("MARS requires either \"gaussian\" or \"binomial\" family.",call.=FALSE)
      }
      # pmethod="cv" not available for multivariate outputs
      
      ### start trial ###
      if(FALSE){
      #nk <- min(200, max(20, 2 * ncol(X0))) + 1
      #nprune <- round(seq(from=2,to=nk,length.out=10))
      #object <- list()
      #for(j in seq_along(nprune)){
      #  object[[j]] <- earth::earth(x=X0,y=y0,nprune=nprune[j],pmethod="cv",nfold=nfolds.int)
      #}
      #sapply(object,function(x) x$gcv)
      ## i.e. run earth/mars with tryCatch for each nprune
      ## and select run with best cvm (here gcv)
      # tune nprune (use default nk)!
      }
      
      #pred$mars[foldid.ext==i,] <- earth:::predict.earth(object=object,newdata=X1,type="response") # original
      pred$mars[foldid.ext==i,] <- stats::predict(object=object,newdata=X1,type="response") # trial
      end <- Sys.time()
      time$mars <- as.numeric(difftime(end,start,units="secs"))
    }
    
    if("spls" %in% compare){
      cat("spls"," ")
      start <- Sys.time()
      if(any(family!="gaussian")){
        stop("spls requires \"gaussian\" family.")
      }
      invisible(utils::capture.output(cv.spls <- spls::cv.spls(x=X0,y=y0,fold=nfolds.int,K=seq_len(min(ncol(X0),10)),
                               eta=seq(from=0.0,to=0.9,by=0.1),plot.it=FALSE))) #scale.x=FALSE
      object <- spls::spls(x=X0,y=y0,K=cv.spls$K.opt,eta=cv.spls$eta.opt) #scale.x=FALSE
      pred$spls[foldid.ext==i,] <- spls::predict.spls(object=object,newx=X1,type="fit")
      end <- Sys.time()
      time$spls <- as.numeric(difftime(end,start,units="secs"))
    }
    
    if("mrce" %in% compare){
      cat("mrce"," ")
      start <- Sys.time()
      if(any(family!="gaussian")){
        stop("MRCE requires \"gaussian\" family.",call.=FALSE)
      }
      lam1 <- lam2 <- 10^seq(from=1,to=-4,length.out=11)
      invisible(utils::capture.output(trials <- lapply(lam2,function(x) tryCatch(expr=MRCE::mrce(X=X0,Y=y0,lam1.vec=lam1,lam2.vec=x,method="cv",kfold=nfolds.int),error=function(x) NULL))))
      cv.err <- sapply(trials,function(x) ifelse(is.null(x),Inf,min(x$cv.err)))
      object <- trials[[which.min(cv.err)]]
      pred$mrce[foldid.ext==i,] <- matrix(object$muhat,nrow=nrow(X1),ncol=q,byrow=TRUE) + X1 %*% object$Bhat
      end <- Sys.time()
      time$mrce <- as.numeric(difftime(end,start,units="secs"))
    }
    
    if("map" %in% compare){
      cat("map"," ")
      start <- Sys.time()
      if(any(family!="gaussian")){
        stop("map requires \"gaussian\" family.")
      }
      mean <- colMeans(y0)
      y0s <- y0-matrix(data=mean,nrow=nrow(X0),ncol=ncol(y0),byrow=TRUE)
      lamL1.v <- lamL2.v <- exp(seq(from=0,to=5,length.out=11))
      cv <- remMap::remMap.CV(X=X0,Y=y0s,lamL1.v=lamL1.v,lamL2.v=lamL2.v,fold=nfolds.int)
      #graphics::plot(x=lamL1.v,y=log(as.numeric(cv$ols.cv[,3])))
      pick <- which.min(as.vector(cv$ols.cv))
      lamL1 <- cv$l.index[1,pick]
      lamL2 <- cv$l.index[2,pick]
      # index <- which(cv$ols.cv==min(cv$ols.cv),arr.ind=TRUE)[1,]
      # rev(lamL1.v)[index[1]]
      # rev(lamL2.v)[index[2]]
      ##cat("lam1:",lamL1,", lam2:",lamL2)
      object <- remMap::remMap(X.m=X0,Y.m=y0s,lamL1=lamL1,lamL2=lamL2)
      pred$map[foldid.ext==i,] <- matrix(data=mean,nrow=nrow(X1),ncol=ncol(y0),byrow=TRUE) + X1 %*% object$phi
      end <- Sys.time()
      time$map <- as.numeric(difftime(end,start,units="secs"))
    }
    
    if("mrf" %in% compare){
      cat("mrf"," ")
      start <- Sys.time()
      if(any(family!="gaussian")){
        stop("mrf requires \"gaussian\" family.")
      }
      pred$mrf[foldid.ext==i,] <- MultivariateRandomForest::build_forest_predict(trainX=X0,trainY=y0,
                                   n_tree=100,m_feature=min(ncol(X0)-1,5),min_leaf=min(nrow(X0),5),testX=X1)
      # use n_tree=500, m_feature=floor(ncol(X0)/3)
      # alternative: IntegratedMRF
      # Check why this does not work well!
      end <- Sys.time()
      time$mrf <- as.numeric(difftime(end,start,units="secs"))
    }
    
    if("sier" %in% compare){
      cat("sier"," ")
      start <- Sys.time()
      if(any(family!="gaussian")){
        stop("SiER requires \"gaussian\" family.",call.=FALSE)
      }
      invisible(utils::capture.output(object <- SiER::cv.SiER(X=X0,Y=y0,K.cv=3)))
      # trial with K.cv=3 (for spped-up)
      # use upper.comp=10 and thres=0.01  (changed for speed-up)
      pred$sier[foldid.ext==i,] <- SiER::pred.SiER(cv.fit=object,X.new=X1)
      end <- Sys.time()
      time$sier <- as.numeric(difftime(end,start,units="secs"))
    }
    
    if("mcen" %in% compare){
      cat("mcen"," ")
      start <- Sys.time()
      if(all(family=="gaussian")){
        type <- "mgaussian"
      } else if(all(family=="binomial")){
        type <- "mbinomial"
      } else {
        stop("mcen requires either \"gaussian\" or \"binomial\".",call.=FALSE)
      }
      object <- mcen::cv.mcen(x=X0,y=y0,family=type,folds=foldid,ky=1,
                              gamma_y=seq(from=0.1,to=5.1,by=1),ndelta=5)
      # TEMPORARY gamma_y=seq(from=0.1,to=5.1,length.out=3) and ndelta=3 (for speed-up)
      #temp <- mcen:::predict.cv.mcen(object=object,newx=X1) # original
      temp <- stats::predict(object=object,newx=X1) # trial
      pred$mcen[foldid.ext==i,] <- as.matrix(temp)
      # single cluster (ky=1) due to setting and error
      end <- Sys.time()
      time$mcen <- as.numeric(difftime(end,start,units="secs"))
    }
    
    if("gpm" %in% compare){
      cat("gpm"," ")
      start <- Sys.time()
      if(any(family!="gaussian")){
        stop("GPM requires \"gaussian\" family.",call.=FALSE)
      }
      object <- GPM::Fit(X=X0,Y=y0)
      pred$gpm[foldid.ext==i,] <- GPM::Predict(XF=X1,Model=object)$YF
      end <- Sys.time()
      time$gpm <- as.numeric(difftime(end,start,units="secs"))
    }
    
    if("rmtl" %in% compare){
      cat("rmtl"," ")
      start <- Sys.time()
      if(all(family=="gaussian")){
        type <- "Regression"
        y0l <- lapply(seq_len(ncol(y0)),function(i) y0[,i,drop=FALSE])
      } else if(all(family=="binomial")){
        type <- "Classification"
        y0l <- lapply(seq_len(ncol(y0)),function(i) 2*y0[,i,drop=FALSE]-1)
      } else {
        stop("RMTL requires either \"gaussian\" or \"binomial\".",call.=FALSE)
      }
      X0l <- lapply(seq_len(ncol(y0)),function(i) X0)
      X1l <- lapply(seq_len(ncol(y0)),function(i) X1)
      #---------------------------
      #--- manual tuning start ---
      Lam1_seq <- c(10^seq(from=1,to=-4,length.out=11),0)
      Lam2_seq <- c(10^seq(from=1,to=-4,length.out=11),0)
      cvMTL <- list()
      seed <- .Random.seed
      for(j in seq_along(Lam2_seq)){
        .Random.seed <- seed
        cvMTL[[j]] <- RMTL::cvMTL(X=X0l,Y=y0l,type=type,Lam1_seq=Lam1_seq,Lam2=Lam2_seq[j],nfolds=nfolds.int)
      }
      cvm <- vapply(X=cvMTL,FUN=function(x) min(x$cvm),FUN.VALUE=numeric(1))
      Lam1 <- cvMTL[[which.min(cvm)]]$Lam1.min
      Lam2 <- Lam2_seq[which.min(cvm)]
      #graphics::plot(x=Lam2_seq,y=cvm)
      #cat(Lam1," ",Lam2,"\n")
      #--- manual tuning end ----
      #--------------------------
      MTL <- RMTL::MTL(X=X0l,Y=y0l,type=type,Lam1=Lam1,Lam2=Lam2)
      #temp <- RMTL:::predict.MTL(object=MTL,newX=X1l) # original
      temp <- stats::predict(object=MTL,newX=X1l)
      pred$rmtl[foldid.ext==i,] <- do.call(what="cbind",args=temp)
      end <- Sys.time()
      time$rmtl <- as.numeric(difftime(end,start,units="secs"))
    }
    
    if("mtps" %in% compare){
      cat("mtps"," ")
      start <- Sys.time()
      
      if(alpha.base==0){
        step1 <- MTPS::glmnet.ridge
      } else if(alpha.base==1){
        step1 <- MTPS::glmnet.lasso
      } else {
        stop("MTPS requires alpha.base 0 or 1.",call.=FALSE)
      }
      
      #-------------------------------
      #--- manual lambda.min start ---
      #body <- body(step1)
      #body[[6]][[3]][[2]][[3]] <- "lambda.min"
      #body[[7]][[3]][[4]] <- "lambda.min"
      #body[[8]][[3]][[3]][[2]][[4]] <- "lambda.min"
      #body(step1) <- body
      #if(alpha.base==0){
      #  assignInNamespace(x="glmnet.ridge",value=step1,ns="MTPS")
      #} else if(alpha.base==1){
      #  assignInNamespace(x="glmnet.lasso",value=step1,ns="MTPS")
      #}
      #--- manual lambda.min end ---
      #-----------------------------
      
      #---------------------------
      #--- manual kmeans start ---
      #fun <- MTPS::cv.multiFit
      #body <- body(fun)
      #body[[11]][[3]][[3]] <- nrow(unique(y0))
      #body(fun) <- body
      #assignInNamespace(x="cv.multiFit",value=step1,ns="MTPS")
      #--- manual kmeans end ---
      #-------------------------
      
      if(all(family=="gaussian")){
        step2 <- MTPS::rpart1
      } else if(all(family=="binomial")){
        step2 <- MTPS::glm1
        #y0 <- as.data.frame(y0)
        #y0 <- as.data.frame(lapply(y0,function(x) factor(x)))
      } else {
        stop("MTPS requires family gaussian or binomial.",call.=FALSE)
      }
      
      object <- MTPS::MTPS(xmat=X0,ymat=y0,family=family,nfold=nfolds.int,
                           method.step1=step1,method.step2=step2,cv=TRUE,residual=TRUE)
      pred$mtps[foldid.ext==i,] <- MTPS::predict.MTPS(object=object,newdata=X1,type="response")
      end <- Sys.time()
      time$mtps <- as.numeric(difftime(end,start,units="secs"))
      # now using cross-validation residual stacking (CVRS) 
    }
    
    if(is.character(compare)){cat("\n")}
    
    # --- development ---
    
    if(FALSE){
      ## --- MLPUGS --- (binary outcome only)
      #X0f <- as.data.frame(X0)
      #y0f <- as.data.frame(y0)
      #X1f <- as.data.frame(X1)
      #object <- MLPUGS::ecc(x=X0f,y=y0f,.f=randomForest::randomForest)
      #pred_ecc <- predict(object,newdata=X1f,n.iters=300,burn.in=100,thin=2,
      #    .f = function(rF,newdata){randomForest:::predict.randomForest(rF, newdata, type = "prob")})
      #pred$ecc[foldid.ext==i,] <- summary(pred_ecc,type="prob")
      ## --- MSGLasso --- (many user inputs)
      #MSGLasso::MSGLasso.cv(X=X0,Y=Y0)
      ## --- PMA --- (not for prediction?)
      ## --- MSP --- (not for hd data?)
      #object <- MBSP::mbsp.tpbn(X=X0,Y=Y0)
      #X1 %*% object$B # adjust for intercept
      ## --- bgsmtr --- (for SNPs only?)
      #temp <- bgsmtr::bgsmtr(X=t(X0),Y=t(Y0),group=rep(1,times=ncol(X0)))
      ## --- MGLM --- (for multinomial data only?)
      #MGLM::MGLMsparsereg.fit()
      #MGLM::MGLMtune()
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
  if(times){
    loss <- list(loss=loss,time=time)
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

# 
# MTPS.MTPS <- function (xmat, ymat, family, cv = FALSE, residual = TRUE, nfold = 5, 
#           method.step1, method.step2, resid.type = c("deviance", 
#                                                      "pearson", "raw"), resid.std = FALSE) 
# {
#   resid.type <- match.arg(resid.type)
#   ny <- ncol(ymat)
#   if (length(family) == 1) {
#     if (!family %in% c("gaussian", "binomial")) {
#       stop("family must be gaussian or binomial")
#     }
#     if (family == "gaussian") {
#       family = rep("gaussian", ny)
#     }
#     else if (family == "binomial") {
#       family = rep("binomial", ny)
#     }
#   }
#   if (length(family) != ny) {
#     stop("length of family must be consistent with response")
#   }
#   if (sum(family %in% c("gaussian", "binomial")) != 
#       ny) {
#     stop("family must be gaussian or binomial or their combination")
#   }
#   if (length(method.step1) == 1) {
#     method.step1 <- rep(list(method.step1), ny)
#   }
#   if (length(method.step2) == 1) {
#     method.step2 <- rep(list(method.step2), ny)
#   }
#   if (length(method.step1) != ny) {
#     stop("length of method.step1 must be 1 or the same as response column")
#   }
#   if (length(method.step2) != ny) {
#     stop("length of method.step2 must be 1 or the same as response column")
#   }
#   #for (ii in 1:ny) {
#   #  if (!check.match(family[ii], FUN = method.step1[[ii]])) {
#   #    stop("method.step1 must be consistent with response category")
#   #  }
#   #}
#   if (!residual) {
#     for (ii in 1:ny) {
#       if (!MTPS::check.match(family[ii], FUN = method.step2[[ii]])) {
#         stop("method.step2 must be consistent with response category")
#       }
#     }
#   }
#   else {
#     for (ii in 1:ny) {
#       if (!MTPS::check.match("gaussian", FUN = method.step2[[ii]])) {
#         stop("residual stacking does not allow binary outcome model in second step")
#       }
#     }
#   }
#   if (cv) {
#     fit1 <- MTPS::cv.multiFit(xmat = xmat, ymat = ymat, nfold = nfold, 
#                         method = method.step1, family = family)
#   }
#   else {
#     fit1 <- MTPS.multiFit(xmat = xmat, ymat = ymat, method = method.step1, 
#                      family = family)
#   }
#   pred1 <- fit1$y.fitted
#   if (residual) {
#     fit2 <- MTPS::rs.multiFit(yhat = pred1, ymat = ymat, xmat = xmat, 
#                         family = family, resid.type = resid.type, resid.std = resid.std, 
#                         method = method.step2)
#   }
#   else {
#     fit2 <- MTPS.multiFit(xmat = pred1, ymat = ymat, method = method.step2, 
#                      family = family)
#   }
#   fit <- list(fit1 = fit1, fit2 = fit2, cv = cv, residual = residual)
#   class(fit) <- "MTPS"
#   return(fit)
# }
# 
# MTPS.multiFit <- function (xmat, ymat, method, family = family) 
# {
#   ny <- ncol(ymat)
#   nx <- ncol(xmat)
#   if (length(family) == 1) {
#     if (!family %in% c("gaussian", "binomial")) {
#       stop("family must be gaussian or binomial")
#     }
#     if (family == "gaussian") {
#       family = rep("gaussian", ny)
#     }
#     else if (family == "binomial") {
#       family = rep("binomial", ny)
#     }
#   }
#   if (length(family) != ny) {
#     stop("length of family must be consistent with response")
#   }
#   if (sum(family %in% c("gaussian", "binomial")) != 
#       ny) {
#     stop("each family must be gaussian or binomial")
#   }
#   if (length(method) == 1) {
#     method <- rep(list(method), ny)
#   }
#   if (length(method) != ny) {
#     stop("length of method.step1 must be 1 or the same as response column")
#   }
#   #for (ii in 1:ny) {
#   #  if (!check.match(family[ii], FUN = method[[ii]])) {
#   #    stop("method.step1 must be consistent with response category")
#   #  }
#   #}
#   y.fitted <- ymat
#   y.fitted[!is.na(y.fitted)] <- NA
#   models <- vector("list", ny)
#   colnames(y.fitted) <- names(models) <- colnames(ymat)
#   fit <- vector("list", ny)
#   colnames(xmat) <- paste0("X", 1:nx)
#   for (kk in 1:ny) {
#     fit[[kk]] <- method[[kk]](xmat, ymat[, kk], family = family[kk])
#     models[[kk]] <- fit[[kk]]$model
#     y.fitted[, kk] <- fit[[kk]]$y.fitted
#   }
#   multiFit.fits <- list(fit = fit, y.fitted = y.fitted, model = models, 
#                         method = method, family = family)
#   class(multiFit.fits) <- "multiFit"
#   return(multiFit.fits)
# }
# 
# 
# MTPS.glmnet.ridge <- function (xmat, ymat, family, alpha = 0, ...)
# {
#   tmp0 <- data.frame(yy = ymat, xmat)
#   if (family == "binomial")
#     ymat <- factor(ymat)
#   foldid <- MTPS::createFolds(ymat, k = 5, list = F)
#   model <- MTPS::cv.glmnet2(xmat, ymat, alpha = alpha, foldid = foldid,
#                       family = family, ...)
#   coef.mat <- as.numeric(coef(model, s = "lambda.min"))
#   y.fitted <- predict(model, xmat, s = "lambda.min",
#                       type = "response")
#   predFun <- function(model, xnew) {
#     predict(model, newx = as.matrix(xnew), s = "lambda.min",
#             type = "response")
#   }
#   return(list(model = model, y.fitted = y.fitted, predFun = predFun))
# }
# 
# MTPS.glmnet.lasso <- function (xmat, ymat, family, alpha = 1, ...)
# {
#   tmp0 <- data.frame(yy = ymat, xmat)
#   if (family == "binomial")
#     ymat <- factor(ymat)
#   foldid <- MTPS::createFolds(ymat, k = 5, list = F)
#   model <- MTPS::cv.glmnet2(xmat, ymat, alpha = alpha, foldid = foldid,
#                       family = family, ...)
#   coef.mat <- as.numeric(coef(model, s = "lambda.min"))
#   y.fitted <- predict(model, xmat, s = "lambda.min",
#                       type = "response")
#   predFun <- function(model, xnew) {
#     predict(model, newx = as.matrix(xnew), s = "lambda.min",
#             type = "response")
#   }
#   return(list(model = model, y.fitted = y.fitted, predFun = predFun))
# }
# 
# 
# 
