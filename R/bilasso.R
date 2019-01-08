

# transforms to unit interval
.prob <- function(x,cutoff,shape1,shape2,plot=FALSE){
  q <- exp(x-cutoff)/(1+exp(x-cutoff))
  p <- stats::pbeta(q=q,shape1=shape1,shape2=shape2)
  xlim <- range(c(x,-x))
  mu <- stats::pbeta(q=0.5,shape1=shape1,shape2=shape2)
  a <- rep(NA,times=length(p))
  a[p<mu] <- 0.5*(p/mu)[p<mu]
  a[p>mu] <- 1-0.5*((1-p)/(1-mu))[p>mu]
  if(plot){
    graphics::plot(x=x,y=a,ylim=c(0,1),xlim=xlim)
    graphics::abline(v=cutoff,lty=2,col="red")
    graphics::abline(h=0.5,lty=2,col="grey")
  }
  return(a)
}

#.prob(x=rnorm(100),cutoff=2,shape1=2,shape2=1,plot=TRUE)


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
#' @param sigma
#' sigma sequence\strong{:}
#' vector of increasing positive values;
#' or \code{NULL} (default sequence)
#' 
#' @param nsigma
#' number of \code{sigma} values
#' 
#' @param logistic
#' fit logistic regression for comparison\strong{:}
#' logical
#' (currently all methods require \code{TRUE})
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
#' net <- bilasso(y=y,cutoff=0,X=X)
#' ### Add ... to all glmnet::glmnet calls !!! ###
bilasso <- function(y,cutoff,X,npi=100,nsigma=100,sigma=NULL,nfolds=10,foldid=NULL,type.measure="deviance",logistic=TRUE,...){
  
  #--- temporary ---
  # cutoff <- 0; npi <- 100; nsigma <- 99; sigma <- NULL; nfolds <- 10;  foldid <- NULL; type.measure <- "deviance"; logistic <- TRUE
  test <- list()
  
  test$sigma <- TRUE
  test$pi <- TRUE
  test$trial <- TRUE
  test$unit <- TRUE

  test$grid <- TRUE
  test$max <- FALSE
  test$grid2 <- FALSE
  
  test$calibrate <- FALSE
  
  #--- checks ---
  colasso:::.check(x=y,type="vector")
  if(all(y %in% c(0,1))){warning("Binary response.",call.=FALSE)}
  colasso:::.check(x=cutoff,type="scalar",min=min(y),max=max(y))
  colasso:::.check(x=X,type="matrix")
  if(length(y)!=nrow(X)){stop("Contradictory sample size.",call.=FALSE)}
  colasso:::.check(x=npi,type="scalar",min=1)
  colasso:::.check(x=nsigma,type="scalar",min=1)
  colasso:::.check(x=sigma,type="vector",min=.Machine$double.eps,null=TRUE)
  colasso:::.check(x=nfolds,type="scalar",min=3)
  colasso:::.check(x=foldid,type="vector",values=seq_len(nfolds),null=TRUE)
  colasso:::.check(x=type.measure,type="string",values=c("deviance","class","mse","mae")) # not auc (min/max confusion)
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
  fit$gaussian <- glmnet::glmnet(y=y,x=X,family="gaussian")
  if(logistic){
    fit$binomial <- glmnet::glmnet(y=z,x=X,family="binomial")
  }
  
  #--- sigma, nsigma ---
  if(is.null(sigma)){
    fit$sigma <- exp(seq(from=log(0.05*stats::sd(y)),
                  to=log(10*stats::sd(y)),length.out=nsigma))
  } else {
    fit$sigma <- sigma
    nsigma <- length(sigma)
  }

  #--- tuning parameters ---
  fit$lambda <- fit$gaussian$lambda
  nlambda <- length(fit$gaussian$lambda)
  lab.sigma <- paste0("si",seq_len(nsigma))
  lab.lambda <- paste0("la",seq_len(nlambda))
  names(fit$sigma) <- lab.sigma
  names(fit$lambda) <- lab.lambda
  
  if(test$pi){
    fit$pi <- seq(from=0,to=1,length.out=npi)
    lab.pi <- paste0("pi",seq_len(npi))
  }
  if(test$max){
    fit$max <- exp(seq(from=log(0.05*max(abs(y-cutoff))),
                to=log(max(abs(y-cutoff))),length.out=100))
  }
  
  if(test$unit){
    fit$shape1 <- seq(from=0.01,to=10,length.out=100)
    fit$shape2 <- seq(from=0.01,to=10,length.out=100)
    lab.s1 <- paste0("a",1:100)
    lab.s2 <- paste0("b",1:100)
  }
  
  #--- cross-validation ---
  pred <- list()
  pred$y  <- matrix(data=NA,nrow=n,ncol=nlambda)
  if(logistic){
    pred$z  <- matrix(data=NA,nrow=n,ncol=length(fit$binomial$lambda))
  }
  if(test$sigma){
    pred$sigma <- matrix(data=NA,nrow=n,ncol=nsigma)
  }
  if(test$pi){
    pred$pi <- matrix(data=NA,nrow=n,ncol=length(fit$pi))
  }
  if(test$max){
    pred$max <- matrix(data=NA,nrow=n,ncol=length(fit$max))
  }
  if(test$grid){
    dimnames <- list(NULL,lab.sigma,lab.lambda)
    pred$grid <- array(data=NA,dim=c(n,nsigma,nlambda),dimnames=dimnames)
  }
  if(test$grid2){
    dimnames <- list(NULL,lab.pi,lab.lambda)
    pred$grid2 <- array(data=NA,dim=c(n,npi,nlambda),dimnames=dimnames)
  }
  if(test$trial){
    dimnames <- list(NULL,lab.sigma,lab.pi)
    pred$trial <- array(data=NA,dim=c(n,nsigma,npi),dimnames=dimnames)
  }
  if(test$unit){
    dimnames <- list(NULL,lab.s1,lab.s2)
    pred$unit <- array(data=NA,dim=c(n,100,100),dimnames=dimnames)
  }

  for(k in seq_len(nfolds)){

    y0 <- y[foldid!=k]
    y1 <- y[foldid==k]
    z0 <- z[foldid!=k]
    z1 <- z[foldid==k]
    X0 <- X[foldid!=k,,drop=FALSE]
    X1 <- X[foldid==k,,drop=FALSE]
    
    # linear regression
    net <- glmnet::glmnet(y=y0,x=X0,family="gaussian")
    temp_y <- stats::predict(object=net,newx=X1,type="response",s=fit$gaussian$lambda)
    pred$y[foldid==k,seq_len(ncol(temp_y))] <- temp_y
    
    # logistic regression
    if(logistic){
      net <- glmnet::glmnet(y=z0,x=X0,family="binomial")
      temp_z <- stats::predict(object=net,newx=X1,type="response",s=fit$binomial$lambda)
      pred$z[foldid==k,seq_len(ncol(temp_z))] <- temp_z
    }

    # fusion (sigma)
    if(test$sigma){
      cvm <- colasso:::.loss(y=y1,fit=temp_y,family="gaussian",type.measure="deviance")[[1]]
      y_hat <- temp_y[,which.min(cvm)]
      for(i in seq_along(fit$sigma)){
        pred$sigma[foldid==k,i] <- stats::pnorm(q=y_hat,mean=cutoff,sd=fit$sigma[i])
      }
    }

    # fusion (grid)
    if(test$grid){
      for(i in seq_along(fit$sigma)){
        pred$grid[foldid==k,i,] <- stats::pnorm(q=temp_y,mean=cutoff,sd=fit$sigma[i])
      }
    }
    
    # fusion (grid2)
    if(test$grid2){
      for(i in seq_along(fit$sigma)){
        cont <- stats::pnorm(q=temp_y,mean=cutoff,sd=stats::sd(y))
        pred$grid2[foldid==k,i,] <- fit$pi[i]*cont + (1-fit$pi[i])*temp_z
      }
    }
    
    # fusion (pi)
    if(test$pi){
      cvm <- colasso:::.loss(y=z1,fit=temp_z,family="binomial",type.measure=type.measure)[[1]]
      z_hat <- temp_z[,which.min(cvm)]
      for(i in seq_along(fit$pi)){
        cont <- stats::pnorm(q=y_hat,mean=cutoff,sd=stats::sd(y))
        pred$pi[foldid==k,i] <- fit$pi[i]*cont + (1-fit$pi[i])*z_hat
      }
    }
    
    # fusion (max)
    if(test$max){
        for(i in seq_along(fit$max)){
        temp <- ((y_hat-cutoff)/fit$max[i] + 1)/2
        pred$max[foldid==k,i] <- pmax(0,pmin(temp,1))
      }
    }
    
    # fusion (trial)
    if(test$trial){
      for(i in seq_along(fit$sigma)){
        for(j in seq_along(fit$pi)){
          cont <- stats::pnorm(q=y_hat,mean=cutoff,sd=fit$sigma[i])
          pred$trial[foldid==k,i,j] <- fit$pi[j]*cont + (1-fit$pi[j])*z_hat
        }
      }
    }
    
    # fusion (unit)
    if(test$unit){
      for(i in seq_len(100)){
        for(j in seq_len(100)){
          pred$unit[foldid==k,i,j] <- .prob(x=y_hat,cutoff=cutoff,shape1=fit$shape1[i],shape2=fit$shape2[j])
        }
      }
    }
    
  }
  
  #--- evaluation ---
  
  # deviance (not comparable between Gaussian and binomial families)
  fit$gaussian$cvm <- colasso:::.loss(y=y,fit=pred$y,family="gaussian",type.measure="deviance")[[1]]
  fit$gaussian$lambda.min <- fit$gaussian$lambda[which.min(fit$gaussian$cvm)]
  
  if(logistic){
    fit$binomial$cvm <- colasso:::.loss(y=z,fit=pred$z,family="binomial",type.measure=type.measure)[[1]]
    fit$binomial$lambda.min <- fit$binomial$lambda[which.min(fit$binomial$cvm)]
  }

  if(test$sigma){
    fit$sigma.cvm <- colasso:::.loss(y=z,fit=pred$sigma,family="binomial",type.measure=type.measure)[[1]]
    fit$sigma.min <- fit$sigma[which.min(fit$sigma.cvm)]
  }

  #graphics::plot(x=fit$sigma,y=fit$sigma.cvm)
  #graphics::abline(v=fit$sigma.min,col="red",lty=2)
  #graphics::abline(v=stats::sd(y),col="grey",lty=2)
  
  if(test$pi){
    fit$pi.cvm <- colasso:::.loss(y=z,fit=pred$pi,family="binomial",type.measure=type.measure)[[1]] # trial
    fit$pi.min <- fit$pi[which.min(fit$pi.cvm)] # trial
  }

  if(test$max){
    fit$max.cvm <- colasso:::.loss(y=z,fit=pred$max,family="binomial",type.measure=type.measure)[[1]] # trial
    fit$max.min <- fit$max[which.min(fit$max.cvm)] # trial
  }
  
  if(test$grid){
    dimnames <- list(lab.sigma,lab.lambda)
    fit$cvm <- matrix(data=NA,nrow=nsigma,ncol=nlambda,dimnames=dimnames)
    for(i in seq_len(nsigma)){
      for(j in seq_len(nlambda)){
        fit$cvm[i,j] <- colasso:::.loss(y=z,fit=pred$grid[,i,j],family="binomial",type.measure=type.measure)[[1]]
      }
    }
    temp <- which(fit$cvm==min(fit$cvm),arr.ind=TRUE)
    if(nrow(temp)>1){warning("MULTIPLE!",call.=FALSE);temp <- temp[1,,drop=FALSE]}
    fit$grid.min <- data.frame(sigma=fit$sigma[temp[,"row"]],lambda=fit$gaussian$lambda[temp[,"col"]])
  }
  
  if(test$grid2){
    dimnames <- list(lab.sigma,lab.lambda)
    fit$cvm2 <- matrix(data=NA,nrow=nsigma,ncol=nlambda,dimnames=dimnames)
    for(i in seq_len(nsigma)){
      for(j in seq_len(nlambda)){
        fit$cvm2[i,j] <- colasso:::.loss(y=z,fit=pred$grid2[,i,j],family="binomial",type.measure=type.measure)[[1]]
      }
    }
    temp <- which(fit$cvm2==min(fit$cvm2),arr.ind=TRUE)
    if(nrow(temp)>1){warning("MULTIPLE!",call.=FALSE);temp <- temp[1,]}
    fit$grid2.min <- data.frame(sigma=fit$sigma[temp[,"row"]],lambda=fit$gaussian$lambda[temp[,"col"]])
  }
  
  if(test$trial){
    dimnames <- list(lab.sigma,lab.pi)
    fit$trial.cvm <- matrix(data=NA,nrow=nsigma,ncol=npi,dimnames=dimnames)
    for(i in seq_len(nsigma)){
      for(j in seq_len(npi)){
        fit$trial.cvm[i,j] <- colasso:::.loss(y=z,fit=pred$trial[,i,j],family="binomial",type.measure=type.measure)[[1]]
      }
    }
    temp <- which(fit$trial.cvm==min(fit$trial.cvm),arr.ind=TRUE)
    if(nrow(temp)>1){warning("MULTIPLE!",call.=FALSE);temp <- temp[1,,drop=FALSE]}
    fit$trial.min <- data.frame(sigma=fit$sigma[temp[,"row"]],pi=fit$pi[temp[,"col"]])
  }
  
  if(test$unit){
    dimnames <- list(lab.s1,lab.s2)
    fit$unit.cvm <- matrix(data=NA,nrow=100,ncol=100,dimnames=dimnames)
    for(i in seq_len(100)){
      for(j in seq_len(100)){
        fit$unit.cvm[i,j] <- colasso:::.loss(y=z,fit=pred$unit[,i,j],family="binomial",type.measure=type.measure)[[1]]
      }
    }
    temp <- which(fit$unit.cvm==min(fit$unit.cvm),arr.ind=TRUE)
    if(nrow(temp)>1){warning("MULTIPLE!",call.=FALSE);temp <- temp[1,,drop=FALSE]}
    fit$unit.min <- data.frame(shape1=fit$shape1[temp[,"row"]],shape2=fit$shape2[temp[,"col"]])
    cond <- fit$unit.cvm[lab.s1[temp[,"row"]],lab.s2[temp[,"col"]]]==min(fit$unit.cvm) # check
    if(!cond){stop("internal mistake",call.=FALSE)}
  }
  
  # calibrate
  if(test$calibrate){
    fit$calibrate <- CalibratR::calibrate(actual=z,predicted=pred$y[,which.min(fit$gaussian$cvm)],nCores=1,model_idx=5)$calibration_models
  }
  
  #--- return ---
  fit$cutoff <- cutoff
  fit$info <- list(type.measure=type.measure,
                   sd.y=stats::sd(y),
                   table=table(z),
                   test=test)

  class(fit) <- "bilasso"
  return(fit)
}


coef.bilasso <- function(x){
  s <- x$gaussian$lambda.min
  beta <- glmnet::coef.glmnet(object=x$gaussian,s=s)
  
  s <- x$binomial$lambda.min
  gamma <- glmnet::coef.glmnet(object=x$binomial,s=s)
  
  coef <- cbind(beta,gamma)
  colnames(coef) <- c("beta","gamma")
  return(coef)
}

# plot.bilasso <- function(x){
#   #graphics::plot(x=x$sigma,y=x$sigma.cvm,xlab=expression(sigma),ylab="deviance")
#   #graphics::abline(v=x$sigma.min,lty=2,col="red")
#   #graphics::abline(v=x$info$sd.y,lty=2,col="grey")
#   
#   ### original ###
#   #x$grid.cvm[40,40] <- -100
#   #k <- 100
#   #levels <- quantile(x$grid.cvm,probs=seq(from=0,to=1,length.out=k+1))
#   #col <- colorspace::diverge_hsv(n=k)
#   #graphics::filled.contour(x$grid.cvm,xlab="",ylab="",levels=levels,col=col,)
#   
#   ### trial ###
#   k <- 100
#   levels <- stats::quantile(x$cvm,probs=seq(from=0,to=1,length.out=k+1))
#   col <- colorspace::diverge_hsv(n=k)
#   nsigma <- length(x$sigma)
#   nlambda <- length(x$gaussi$lambda)
#   
#   sigma.min <- x$grid.min$sigma
#   lambda.min <- x$grid.min$lambda
#   
#   graphics::plot.new()
#   graphics::par(xaxs="i",yaxs="i")
#   graphics::plot.window(xlim=c(1-0.5,nsigma+0.5),ylim=c(1-0.5,nlambda+0.5))
#   
#   sel <- which(x$sigma==sigma.min)
#   graphics::axis(side=1,at=c(1,sel,nsigma),labels=signif(x$sigma[c(1,sel,nsigma)],digits=2))
#   
#   sel <- which(x$gaussian$lambda==lambda.min)
#   graphics::axis(side=2,at=c(1,sel,nlambda),labels=signif(x$gaussian$lambda[c(1,sel,nlambda)],digits=2))
#   
#   graphics::title(xlab=expression(sigma),ylab=expression(lambda))
#   #graphics::.filled.contour(x=seq_along(x$sigma),y=seq_along(x$gaussian$lambda),z=x$cvm,levels=levels,col=col)
#   graphics::image(x=seq_along(x$sigma),y=seq_along(x$gaussian$lambda),z=x$cvm,breaks=levels,col=col,add=TRUE)
#   graphics::box()
# 
# }


plot.bilasso <- function(x){

  k <- 100
  levels <- stats::quantile(x$trial.cvm,probs=seq(from=0,to=1,length.out=k+1))
  col <- colorspace::diverge_hsv(n=k)
  nsigma <- length(x$sigma)
  npi <- length(x$pi)
  
  sigma.min <- x$trial.min$sigma
  pi.min <- x$trial.min$pi
  
  graphics::plot.new()
  graphics::par(xaxs="i",yaxs="i")
  graphics::plot.window(xlim=c(1-0.5,nsigma+0.5),ylim=c(1-0.5,npi+0.5))
  
  sel <- which(x$sigma==sigma.min)
  graphics::axis(side=1,at=c(1,sel,nsigma),labels=signif(x$sigma[c(1,sel,nsigma)],digits=2))
  
  sel <- which(x$pi==pi.min)
  graphics::axis(side=2,at=c(1,sel,npi),labels=signif(x$pi[c(1,sel,npi)],digits=2))
  
  graphics::title(xlab=expression(sigma),ylab=expression(pi))
  #graphics::.filled.contour(x=seq_along(x$sigma),y=seq_along(x$pi),z=x$cvm,levels=levels,col=col)
  graphics::image(x=seq_along(x$sigma),y=seq_along(x$pi),z=x$trial.cvm,breaks=levels,col=col,add=TRUE)
  graphics::box()
  
}


predict.bilasso <- function(x,newx,type="probability"){
  
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
    prob$sigma <- stats::pnorm(q=link,mean=x$cutoff,sd=x$sigma.min) # original
  }
  
  ### DELETE THE FOLLOWING LINE ###
  ###prob$sigma <- stats::pnorm(q=link,mean=x$cutoff,sd=max(x$sigma.min,x$info$sd.y)) # delete
  ### DELETE THE PREVIOUS LINE ###
  
  if(test$pi){
     prob$pi <- x$pi.min*prob$gaussian + (1-x$pi.min)*prob$binomial # trial pi
  }
 
  if(test$max){
    temp <- ((link-x$cutoff)/x$max.min + 1)/2 # trial max
    prob$max <- pmax(0,pmin(temp,1)) # trial max
  }  
  
  if(test$grid){
    temp <- as.numeric(stats::predict(object=x$gaussian,
                                      newx=newx,s=x$grid.min$lambda,type="response"))
    prob$grid <- stats::pnorm(q=temp,mean=x$cutoff,sd=x$grid.min$sigma)
  }
  
  if(test$grid2){
    temp <- as.numeric(stats::predict(object=x$gaussian,
                                      newx=newx,s=x$grid2.min$lambda,type="response"))
    prob$grid2 <- x$grid2.min$pi*temp + (1-x$grid2.min$pi)*prob$binomial
  }
  
  if(test$trial){
    cont <- stats::pnorm(q=link,mean=x$cutoff,sd=x$trial.min$sigma)
    prob$trial <- x$trial.min$pi*cont + (1-x$trial.min$pi)*prob$binomial
  }
  
  if(test$unit){
    prob$unit <- .prob(x=link,cutoff=x$cutoff,shape1=x$unit.min$shape1,shape2=x$unit.min$shape2)
  }
  
  if(test$calibrate){
    prob$calibrate <- CalibratR::predict_calibratR(calibration_models=x$calibrate,new=link,nCores=1)$GUESS_2
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

#' @export
#' @title
#' Comparison
#'
#' @description
#' Compares models for a continuous response with a cutoff value
#'
#' @inheritParams  bilasso
#'
#' @examples
#' NA
#' 
bilasso_compare <- function(y,cutoff,X,nfolds=5,type.measure="deviance"){
  
  z <- 1*(y > cutoff)
  fold <- palasso:::.folds(y=z,nfolds=nfolds)
  
  cols <- c("gaussian","binomial","sigma","pi","trial","grid","unit")
  pred <- matrix(data=NA,nrow=length(y),ncol=length(cols),
                 dimnames=list(NULL,cols))
  
  select <- list()
  for(i in seq_len(nfolds)){
    fit <- colasso::bilasso(y=y[fold!=i],cutoff=cutoff,X=X[fold!=i,],logistic=TRUE,type.measure=type.measure)
    tryCatch(expr=colasso:::plot.bilasso(fit),error=function(x) NULL)
    #colasso:::plot.bilasso(fit)
    temp <- colasso:::predict.bilasso(fit,newx=X[fold==i,])
    if(any(temp<0|temp>1)){stop("Outside unit interval.",call.=FALSE)}
    model <- colnames(pred)
    for(j in seq_along(model)){
      pred[fold==i,model[j]] <- temp[[model[j]]]
    }
  }
  
  type <- c("deviance","class","mse","mae","auc")
  loss <- lapply(X=type,FUN=function(x) colasso:::.loss(y=z,fit=pred,family="binomial",type.measure=x,foldid=fold)[[1]])
  names(loss) <- type

  return(loss)
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

# Verifies whether two or more arguments are identical.
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


#.glmnet <- function(y,x,family,args){
#  args$y <- y
#  args$x <- X
#  args$family <- "gaussian"
#  do.call(what=glmnet::glmnet,args=args)
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


