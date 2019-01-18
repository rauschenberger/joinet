
#--- Workhorse function --------------------------------------------------------

#' @export
#' @title
#' Logistic regression with a continuous response
#' 
#' @description
#' Implements logistic regression with a continuous response.
#'  
#' @param y
#' continuous response\strong{:}
#' vector of length \eqn{n}
#' 
#' @param cutoff
#' cutoff point for dichotomising response into classes\strong{:}
#' value between \code{min(y)} and \code{max(y)}
#' 
#' @param X
#' covariates\strong{:}
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
#' number of folds
#' 
#' @param type.measure
#' loss function for binary classification
#' (linear regression uses the deviance)
#' 
#' @param pi
#' pi sequence\strong{:}
#' vector of values in the unit interval;
#' or \code{NULL} (default sequence)
#' 
#' @param npi
#' number of \code{pi} values
#' 
#' @param sigma
#' sigma sequence\strong{:}
#' vector of increasing positive values;
#' or \code{NULL} (default sequence)
#' 
#' @param nsigma
#' number of \code{sigma} values
#' 
#' @param ...
#' further arguments passed to \code{\link[glmnet]{glmnet}}
#' 
#' @details
#' - INCLUDE note on deviance (not comparable between lin and log models)
#' - alpha: elastic net parameter\strong{:}
#' numeric between \eqn{0} (ridge) and \eqn{1} (lasso)
#' - do not use "family"
#' 
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' net <- cornet(y=y,cutoff=0,X=X)
#' ### Add ... to all glmnet::glmnet calls !!! ###
cornet <- function(y,cutoff,X,alpha=1,npi=101,pi=NULL,nsigma=99,sigma=NULL,nfolds=10,foldid=NULL,type.measure="deviance",...){
  
  #--- temporary ---
  # cutoff <- 0; npi <- 101; pi <- NULL; nsigma <- 99; sigma <- NULL; nfolds <- 10;  foldid <- NULL; type.measure <- "deviance"; logistic <- TRUE
  test <- list()
  test$sigma <- test$pi <- FALSE
  test$grid <- TRUE
  
  #--- checks ---
  cornet:::.check(x=y,type="vector")
  if(all(y %in% c(0,1))){warning("Binary response.",call.=FALSE)}
  cornet:::.check(x=cutoff,type="scalar",min=min(y),max=max(y))
  cornet:::.check(x=X,type="matrix")
  cornet:::.check(x=alpha,type="scalar",min=0,max=1)
  if(length(y)!=nrow(X)){stop("Contradictory sample size.",call.=FALSE)}
  cornet:::.check(x=npi,type="scalar",min=1)
  cornet:::.check(x=pi,type="vector",min=0,max=1,null=TRUE)
  cornet:::.check(x=nsigma,type="scalar",min=1)
  cornet:::.check(x=sigma,type="vector",min=.Machine$double.eps,null=TRUE)
  cornet:::.check(x=nfolds,type="scalar",min=3)
  cornet:::.check(x=foldid,type="vector",values=seq_len(nfolds),null=TRUE)
  cornet:::.check(x=type.measure,type="string",values=c("deviance","class","mse","mae")) # not auc (min/max confusion)
  if(!is.null(list(...)$family)){stop("Reserved argument \"family\".",call.=FALSE)}
  n <- length(y)
  
  # binarisation
  z <- 1*(y > cutoff)
  
  # fold identifiers
  if(is.null(foldid)){
    foldid <- palasso:::.folds(y=z,nfolds=nfolds)
  }
  
  #--- model fitting ---
  fit <- list()
  fit$gaussian <- glmnet::glmnet(y=y,x=X,family="gaussian",alpha=alpha)
  fit$binomial <- glmnet::glmnet(y=z,x=X,family="binomial",alpha=alpha)
   
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

  if(test$sigma){
    pred$sigma <- matrix(data=NA,nrow=n,ncol=nsigma)
  }
  if(test$pi){
    pred$pi <- matrix(data=NA,nrow=n,ncol=npi)
  }
  if(test$grid){
    dimnames <- list(NULL,lab.sigma,lab.pi)
    pred$grid <- array(data=NA,dim=c(n,nsigma,npi),dimnames=dimnames)
  }

  for(k in seq_len(nfolds)){

    y0 <- y[foldid!=k]
    y1 <- y[foldid==k]
    z0 <- z[foldid!=k]
    z1 <- z[foldid==k]
    X0 <- X[foldid!=k,,drop=FALSE]
    X1 <- X[foldid==k,,drop=FALSE]
    
    # linear regression
    net <- glmnet::glmnet(y=y0,x=X0,family="gaussian",alpha=alpha)
    temp_y <- stats::predict(object=net,newx=X1,type="response",s=fit$gaussian$lambda)
    pred$y[foldid==k,seq_len(ncol(temp_y))] <- temp_y
    cvm <- cornet:::.loss(y=y1,fit=temp_y,family="gaussian",type.measure="deviance")[[1]]
    y_hat <- temp_y[,which.min(cvm)]
    
    # logistic regression
    net <- glmnet::glmnet(y=z0,x=X0,family="binomial",alpha=alpha)
    temp_z <- stats::predict(object=net,newx=X1,type="response",s=fit$binomial$lambda)
    pred$z[foldid==k,seq_len(ncol(temp_z))] <- temp_z
    cvm <- cornet:::.loss(y=z1,fit=temp_z,family="binomial",type.measure=type.measure)[[1]]
    z_hat <- temp_z[,which.min(cvm)]
    
    # fusion (sigma)
    if(test$sigma){
      for(i in seq_along(fit$sigma)){
        #pred$sigma[foldid==k,i] <- stats::pnorm(q=y_hat,mean=cutoff,sd=fit$sigma[i])
      }
    }
    
    # fusion (pi)
    if(test$pi){
      for(i in seq_along(fit$pi)){
        #cont <- stats::pnorm(q=y_hat,mean=cutoff,sd=stats::sd(y))
        #pred$pi[foldid==k,i] <- fit$pi[i]*cont + (1-fit$pi[i])*z_hat
      }
    }

    # fusion (grid)
    if(test$grid){
      for(i in seq_along(fit$sigma)){
        for(j in seq_along(fit$pi)){
          cont <- stats::pnorm(q=y_hat,mean=cutoff,sd=fit$sigma[i])
          pred$grid[foldid==k,i,j] <- fit$pi[j]*cont + (1-fit$pi[j])*z_hat
        }
      }
    }
  }
  
  #--- evaluation ---
  
  # deviance (not comparable between Gaussian and binomial families)
  fit$gaussian$cvm <- cornet:::.loss(y=y,fit=pred$y,family="gaussian",type.measure="deviance")[[1]]
  fit$gaussian$lambda.min <- fit$gaussian$lambda[which.min(fit$gaussian$cvm)]
  
  fit$binomial$cvm <- cornet:::.loss(y=z,fit=pred$z,family="binomial",type.measure=type.measure)[[1]]
  fit$binomial$lambda.min <- fit$binomial$lambda[which.min(fit$binomial$cvm)]

  if(test$sigma){
    #fit$sigma.cvm <- cornet:::.loss(y=z,fit=pred$sigma,family="binomial",type.measure=type.measure)[[1]]
    #fit$sigma.min1 <- fit$sigma[which.min(fit$sigma.cvm)]
  }

  if(test$pi){
    #fit$pi.cvm <- cornet:::.loss(y=z,fit=pred$pi,family="binomial",type.measure=type.measure)[[1]] # trial
    #fit$pi.min1 <- fit$pi[which.min(fit$pi.cvm)]
  }

  if(test$grid){
    dimnames <- list(lab.sigma,lab.pi)
    fit$cvm <- matrix(data=NA,nrow=nsigma,ncol=npi,dimnames=dimnames)
    for(i in seq_len(nsigma)){
      for(j in seq_len(npi)){
        fit$cvm[i,j] <- cornet:::.loss(y=z,fit=pred$grid[,i,j],family="binomial",type.measure=type.measure)[[1]]
      }
    }
    temp <- which(fit$cvm==min(fit$cvm),arr.ind=TRUE,useNames=TRUE)
    if(nrow(temp)>1){warning("MULTIPLE!",call.=FALSE);temp <- temp[1,,drop=FALSE]}
    fit$sigma.min <- fit$sigma[temp[1]]
    fit$pi.min <- fit$pi[temp[2]]
    if(fit$cvm[names(fit$sigma.min),names(fit$pi.min)]!=min(fit$cvm)){stop("Internal error.")}
  }
  
  #--- return ---
  fit$cutoff <- cutoff
  fit$info <- list(type.measure=type.measure,
                   sd.y=stats::sd(y),
                   table=table(z),
                   test=as.data.frame(test))

  class(fit) <- "cornet"
  return(fit)
}

#--- Methods -------------------------------------------------------------------

#' @export
#' @title
#' to do
#'
#' @description
#' to do
#'
#' @param object
#' cornet object
#' 
#' @param ...
#' to do
#'
#' @examples
#' NA
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
#' to do
#'
#' @description
#' to do
#'
#' @param x
#' to do
#' 
#' @param ...
#' further arguments
#'
#' @examples
#' NA
#' 
plot.cornet <- function(x,...){
  
  if(length(list(...))!=0){warning("Ignoring arguments.")}

  k <- 100
  levels <- stats::quantile(x$cvm,probs=seq(from=0,to=1,length.out=k+1))
  col <- colorspace::diverge_hsv(n=k)
  nsigma <- length(x$sigma)
  npi <- length(x$pi)
  
  graphics::plot.new()
  graphics::par(xaxs="i",yaxs="i")
  graphics::plot.window(xlim=c(1-0.5,nsigma+0.5),ylim=c(1-0.5,npi+0.5))
  
  ssigma <- which(x$sigma==x$sigma.min)
  graphics::axis(side=1,at=c(1,ssigma,nsigma),labels=signif(x$sigma[c(1,ssigma,nsigma)],digits=2))
  
  spi <- which(x$pi==x$pi.min)
  graphics::axis(side=2,at=c(1,spi,npi),labels=signif(x$pi[c(1,spi,npi)],digits=2))
  
  graphics::title(xlab=expression(sigma),ylab=expression(pi))
  #graphics::.filled.contour(x=seq_along(x$sigma),y=seq_along(x$pi),z=x$cvm,levels=levels,col=col)
  graphics::image(x=seq_along(x$sigma),y=seq_along(x$pi),z=x$cvm,breaks=levels,col=col,add=TRUE)
  graphics::box()
  
  #graphics::abline(v=ssigma,lty=2,col="grey")
  #graphics::abline(h=spi,lty=2,col="grey")
  
  graphics::points(x=ssigma,y=spi,pch=4,col="black",cex=1)
  
}

#' @export
#' @title
#' to do
#'
#' @description
#' to do
#'
#' @param object
#' cornet object
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
#' to do
#' 
#' @examples
#' NA
#' 
predict.cornet <- function(object,newx,type="probability",...){
  
  if(length(list(...))!=0){warning("Ignoring arguments.")}
  
  x <- object; rm(object)
  
  test <- x$info$test
  
  .check(x=newx,type="matrix")
  .check(x=type,type="string",values=c("probability","odds","log-odds"))
  
  # linear, logistic and mixed
  prob <- list()
  link <- as.numeric(stats::predict(object=x$gaussian,
                  newx=newx,s=x$gaussian$lambda.min,type="response"))
  prob$gaussian <- stats::pnorm(q=link,mean=x$cutoff,sd=x$info$sd.y)
  prob$binomial <- as.numeric(stats::predict(object=x$binomial,
                  newx=newx,s=x$binomial$lambda.min,type="response"))
  
  if(test$sigma){
    #prob$sigma <- stats::pnorm(q=link,mean=x$cutoff,sd=x$sigma.min)
  }
  
  if(test$pi){
     #prob$pi <- x$pi.min*prob$gaussian + (1-x$pi.min)*prob$binomial
  }
 
  if(test$grid){
    cont <- stats::pnorm(q=link,mean=x$cutoff,sd=x$sigma.min)
    prob$grid <- x$pi.min*cont + (1-x$pi.min)*prob$binomial
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

#' @export
#' @title
#' Comparison
#'
#' @description
#' Compares models for a continuous response with a cutoff value
#'
#' @inheritParams  cornet
#'
#' @examples
#' NA
#' 
.compare <- function(y,cutoff,X,alpha=1,nfolds=5,foldid=NULL,type.measure="deviance"){
  
  z <- 1*(y > cutoff)
  if(is.null(foldid)){
    fold <- palasso:::.folds(y=z,nfolds=nfolds)
  } else {
    fold <- foldid
  }
  
  #cols <- c("gaussian","binomial","sigma","pi","trial","grid","unit")
  cols <- c("gaussian","binomial","grid")
  pred <- matrix(data=NA,nrow=length(y),ncol=length(cols),
                 dimnames=list(NULL,cols))
  
  select <- list()
  for(i in seq_len(nfolds)){
    fit <- cornet::cornet(y=y[fold!=i],cutoff=cutoff,X=X[fold!=i,],alpha=alpha,logistic=TRUE,type.measure=type.measure)
    tryCatch(expr=cornet:::plot.cornet(fit),error=function(x) NULL)
    #cornet:::plot.cornet(fit)
    temp <- cornet:::predict.cornet(fit,newx=X[fold==i,])
    if(any(temp<0|temp>1)){stop("Outside unit interval.",call.=FALSE)}
    model <- colnames(pred)
    for(j in seq_along(model)){
      pred[fold==i,model[j]] <- temp[[model[j]]]
    }
  }
  
  type <- c("deviance","class","mse","mae","auc")
  loss <- lapply(X=type,FUN=function(x) cornet:::.loss(y=z,fit=pred,family="binomial",type.measure=x,foldid=fold)[[1]])
  names(loss) <- type

  return(loss)
}

#' @title
#' Arguments
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
#' NA
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
#' NA
#' 
.check <- function(x,type,miss=FALSE,min=NULL,max=NULL,values=NULL,inf=FALSE,null=FALSE){
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
  if(!miss && any(is.na(x))){
    stop(paste0("Argument \"",name,"\" contains missing values."),call.=FALSE)
  }
  if(!is.null(min) && any(x<min)){
    stop(paste0("expecting ",name," >= ",min),call.=FALSE)
  }
  if(!is.null(max) && any(x>max)){
    stop(paste0("expecting ",name," <= ",max),call.=FALSE)
  }
  if(!is.null(values) && !(x %in% values)){
    stop(paste0("Argument \"",name,"\" contains invalid values."),call.=FALSE)
  }
  if(!inf && any(is.infinite(values))){
    stop(paste0("Argument \"",name,"\ contains infinite values."),call.=FALSE)
  }
  return(invisible(NULL))
}

# Simulates y and X.
.simulate <- function(n,p,prob=0.2,fac=1){
  beta <- stats::rnorm(n=p)
  cond <- stats::rbinom(n=p,size=1,prob=prob)
  beta[cond==0] <- 0
  X <- matrix(stats::rnorm(n=n*p),nrow=n,ncol=p)
  mean <- X %*% beta
  y <- stats::rnorm(n=n,mean=mean,sd=fac*stats::sd(mean))
  return(list(y=y,X=X))
}

#--- start trial ---
if(FALSE){
  n <- 1000
  y_hat <- runif(n)
  y <- y_hat > 0.9
  y <- rbinom(n=n,size=1,prob=0.5)
  foldid <- rep(1:10,length.out=n)
  .loss(y=y,fit=y_hat,family="binomial",type.measure="auc",foldid=foldid)
}
#--- end trial ---



# Correct this function in the palasso package (search twice for "# typo").
.loss <- function (y,fit,family,type.measure,foldid=NULL){
  if (!is.list(fit)) {
    fit <- list(fit)
  }
  loss <- list()
  for (i in seq_along(fit)) {
    if (is.vector(fit[[i]])) {
      fit[[i]] <- as.matrix(fit[[i]])
    }
    if (is.null(foldid) & (family == "cox" | type.measure == 
                           "auc")) {
      stop("Missing foldid.", call. = FALSE)
    }
    if (family == "gaussian") {
      if (type.measure %in% c("deviance", "mse")) {
        loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
                           FUN = function(x) mean((x - y)^2))
      }
      else if (type.measure == "mae") {
        loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
                           FUN = function(x) mean(abs(x - y)))
      }
      else {
        stop("Invalid type measure.", call. = FALSE)
      }
    }
    else if (family == "binomial") {
      if (type.measure == "deviance") {
        limit <- 1e-05
        fit[[i]][fit[[i]] < limit] <- limit
        fit[[i]][fit[[i]] > 1 - limit] <- 1 - limit
        loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
                           FUN = function(x) mean(-2 * (y * log(x) + (1 - 
                                                                        y) * log(1 - x))))
      }
      else if (type.measure == "mse") {
        loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
                           FUN = function(x) 2 * mean((x - y)^2))
      }
      else if (type.measure == "mae") {
        loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
                           FUN = function(x) 2 * mean(abs(x - y)))
      }
      else if (type.measure == "class") {
        loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
                           FUN = function(x) mean(abs(round(x) - y)))
      }
      else if (type.measure == "auc") {
        weights <- table(foldid)
        cvraw <- matrix(data = NA, nrow = length(weights), 
                        ncol = ncol(fit[[i]])) # typo in palasso package !
        for (k in seq_along(weights)) {
          cvraw[k, ] <- apply(X = fit[[i]], MARGIN = 2, 
                              FUN = function(x) glmnet::auc(y = y[foldid == 
                                                                    k], prob = x[foldid == k]))
        }
        loss[[i]] <- apply(X = cvraw, MARGIN = 2, FUN = function(x) stats::weighted.mean(x = x, 
                                                                                         w = weights, na.rm = TRUE))
        names(loss[[i]]) <- colnames(fit[[i]]) # typo in palasso package!
      }
      else {
        stop("Invalid type measure.", call. = FALSE)
      }
    }
    else if (family == "poisson") {
      if (type.measure == "deviance") {
        loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
                           FUN = function(x) mean(2 * (ifelse(y == 0, 
                                                              0, y * log(y/x)) - y + x), na.rm = TRUE))
      }
      else if (type.measure == "mse") {
        loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
                           FUN = function(x) mean((x - y)^2))
      }
      else if (type.measure == "mae") {
        loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
                           FUN = function(x) mean(abs(x - y)))
      }
      else {
        stop("Invalid type measure.", call. = FALSE)
      }
    }
    else if (family == "cox") {
      if (type.measure == "deviance") {
        weights <- tapply(X = y[, "status"], INDEX = foldid, 
                          FUN = sum)
        loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
                           FUN = function(x) stats::weighted.mean(x = x/weights, 
                                                                  w = weights, na.rm = TRUE))
      }
      else {
        stop("Invalid type measure.", call. = FALSE)
      }
    }
    else {
      stop("Invalid family.", call. = FALSE)
    }
    if (sum(diff(is.na(loss[[i]]))) == 1) {
      loss[[i]] <- loss[[i]][!is.na(loss[[i]])]
    }
  }
  return(loss)
}

#--- Lost and found ------------------------------------------------------------

# calibrate (for cornet)
#if(test$calibrate){
#  fit$calibrate <- CalibratR::calibrate(actual=z,predicted=pred$y[,which.min(fit$gaussian$cvm)],nCores=1,model_idx=5)$calibration_models
#}
# calibrate (for predict.cornet)
#if(test$calibrate){
#  prob$calibrate <- CalibratR::predict_calibratR(calibration_models=x$calibrate,new=link,nCores=1)$GUESS_2
#}

.args <- function(...){
  args <- list(...)
  names <- names(formals(glmnet::glmnet))
  if(!is.null(args$family)){
    warning("Unexpected argument \"family\".",call.=FALSE) 
  }
  if(any(!names(args) %in% names)){
    stop("Unexpected argument.",call.=FALSE)
  }
  if(is.null(args$alpha)) {
    args$alpha <- 1
  }
  if(is.null(args$nlambda)){
    args$nlambda <- 100
  }
  if(is.null(args$lambda)){
    if(is.null(args$nlambda)){
      args$nlambda <- 100
    }
  } else {
    args$nlambda <- length(args$lambda)
  }
  return(args)
}





