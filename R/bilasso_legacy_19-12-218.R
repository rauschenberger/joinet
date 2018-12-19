if(FALSE){

#' @export
#' @title
#' Logistic regression with a continuous response
#' 
#' @description
#' Implements penalised logistic regression
#' with both a binary and a continuous response.
#' 
#' @details
#' Finds a compromise between binomial (\eqn{pi=0})
#' and linear (\eqn{pi=1}) regression.
#'  
#' @param y
#' continuous response\strong{:}
#' vector of length \eqn{n}
#' 
#' @param cutoff
#' value between \code{min(y)} and \code{max(y)}
#' 
#' @param X
#' covariates\strong{:}
#' matrix with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#' 
#' @param alpha
#' elastic net parameter\strong{:}
#' numeric between \eqn{0} (ridge)
#' and \eqn{1} (lasso)
#' 
#' @param nfolds
#' number of folds
#' 
#' @param type.measure
#' loss function for logistic regression
#' (linear regression uses the deviance)
#' 
#' @param sigma
#' sigma sequence\strong{:}
#' 
#' @param nsigma
#' number of \code{sigma} values
#' 
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' fit <- bilasso(y,cutoff=0,X)
#' 
bilasso <- function(y,cutoff,X,alpha=1,nfolds=10,type.measure="deviance",res=100){

  # checks
  .check(x=y,type="vector")
  if(all(y %in% c(0,1))){stop("Binary response.",call.=FALSE)}
  .check(x=cutoff,type="scalar",min=min(y),max=max(y))
  .check(x=X,type="matrix")
  .check(x=alpha,type="scalar",min=0,max=1)
  .check(x=nfolds,type="scalar",min=3)
  .check(x=type.measure,type="string",values=c("deviance","class","mse","mae","auc"))
  .check(x=res,type="scalar",min=10)
  if(length(y)!=nrow(X)){stop("Contradictory sample size.",call.=FALSE)}
  
  # binarisation
  z <- 1*(y > cutoff)
  
  # fold identifiers
  foldid <- palasso:::.folds(y=z,nfolds=nfolds)

  # model fitting
  fit <- list()
  fit$gaussian <- glmnet::glmnet(y=y,x=X,family="gaussian",alpha=alpha)
  fit$binomial <- glmnet::glmnet(y=z,x=X,family="binomial",alpha=alpha)
  
  # weights
  fit$pi <- seq(from=0,to=1,length.out=res)
  #fit$base <- exp(seq(from=log(1),to=log(100),length.out=100)) # old base
  fit$base <- exp(seq(from=log(0.05*stats::sd(y)),to=log(2*stats::sd(y)),length.out=res)) # new base
  #fit$grid <- expand.grid(pi=fit$pi,base=fit$base) # temporary
  #fit$grid <- expand.grid(sd0=fit$base,sd1=fit$base) # trial
  fit$max <- exp(seq(from=log(0.05*max(abs(y-cutoff))),
                     to=log(max(abs(y-cutoff))),
                     length.out=res))
  
  # cross-validation
  pred <- list()
  pred$y  <- matrix(data=NA,nrow=length(y),ncol=length(fit$gaussian$lambda))
  pred$z  <- matrix(data=NA,nrow=length(y),ncol=length(fit$binomial$lambda))
  pred$pi <- matrix(data=NA,nrow=length(y),ncol=length(fit$pi))
  pred$base <- matrix(data=NA,nrow=length(y),ncol=length(fit$base))
  pred$max <- matrix(data=NA,nrow=length(y),ncol=length(fit$max))
  #pred$grid <- matrix(data=NA,nrow=length(y),ncol=nrow(fit$grid)) # trial
  
  for(k in unique(foldid)){

    y0 <- y[foldid!=k]
    y1 <- y[foldid==k]
    z0 <- z[foldid!=k]
    z1 <- z[foldid==k]
    X0 <- X[foldid!=k,,drop=FALSE]
    X1 <- X[foldid==k,,drop=FALSE]
    
    # linear regression
    net <- glmnet::glmnet(y=y0,x=X0,family="gaussian",alpha=alpha)
    temp <- stats::predict(object=net,newx=X1,type="response",s=fit$gaussian$lambda)
    cvm <- .loss(y=y1,fit=temp,family="gaussian",type.measure="deviance")[[1]]
    pred$y[foldid==k,seq_len(ncol(temp))] <- temp
    y_hat <- temp[,which.min(cvm)]
    
    # logistic regression
    net <- glmnet::glmnet(y=z0,x=X0,family="binomial",alpha=alpha)
    temp <- stats::predict(object=net,newx=X1,type="response",s=fit$binomial$lambda)
    cvm <- .loss(y=z1,fit=temp,family="binomial",type.measure=type.measure)[[1]]
    pred$z[foldid==k,seq_len(ncol(temp))] <- temp
    z_hat <- temp[,which.min(cvm)]
    
    # fusion (pi)
    for(i in seq_along(fit$pi)){
      #pred$pi[foldid==k,i] <- fit$pi[i]*(y_hat > cutoff) + (1-fit$pi[i])*z_hat # original
      cont <- stats::pnorm(q=y_hat,mean=cutoff,sd=stats::sd(y)) # trial
      pred$pi[foldid==k,i] <- fit$pi[i]*cont + (1-fit$pi[i])*z_hat #trial
    }
    
    # fusion (base)
    for(i in seq_along(fit$base)){
      #pred$base[foldid==k,i] <- 1/(1+fit$base[i]^(cutoff-y_hat)) # old trial
      pred$base[foldid==k,i] <- stats::pnorm(q=y_hat,mean=cutoff,sd=fit$base[i]) # new trial
    }
    
    # fusion (max)
    for(i in seq_along(fit$max)){
      pred$max[foldid==k,i] <- ((y_hat-cutoff)/fit$max[i] + 1)/2
    }
    
    # fusion (pi and base)
    #for(i in seq_len(nrow(fit$grid))){
    #  cont <- stats::pnorm(q=y_hat,mean=cutoff,sd=fit$grid$base[i])
    #  temp <- fit$grid$pi[i]*cont + (1-fit$grid$pi[i])*z_hat
    #  pred$grid[foldid==k,i] <- temp
    #}
    
    ## fusion (trial, two bases)
    #for(i in seq_len(nrow(fit$grid))){
    #  p0 <- stats::pnorm(q=y_hat,mean=cutoff,sd=fit$grid$sd0[i])
    #  p1 <- stats::pnorm(q=y_hat,mean=cutoff,sd=fit$grid$sd1[i])
    #  pred$grid[foldid==k,i] <- ifelse(y_hat<cutoff,p0,p1)
    #}
    
  }
  
  # deviance (not comparable between Gaussian and binomial families)
  fit$gaussian$cvm <- .loss(y=y,fit=pred$y,family="gaussian",type.measure="deviance")[[1]]
  fit$gaussian$lambda.min <- fit$gaussian$lambda[which.min(fit$gaussian$cvm)]
  fit$binomial$cvm <- .loss(y=z,fit=pred$z,family="binomial",type.measure=type.measure)[[1]]
  fit$binomial$lambda.min <- fit$binomial$lambda[which.min(fit$binomial$cvm)]
  fit$pi.cvm <- .loss(y=z,fit=pred$pi,family="binomial",type.measure=type.measure)[[1]]
  fit$pi.min <- fit$pi[which.min(fit$pi.cvm)]
  fit$base.cvm <- .loss(y=z,fit=pred$base,family="binomial",type.measure=type.measure)[[1]]
  fit$base.min <- fit$base[which.min(fit$base.cvm)]
  #fit$grid.cvm <- .loss(y=z,fit=pred$grid,family="binomial",type.measure=type.measure)[[1]] # trial
  #fit$grid.min <- fit$grid[which.min(fit$grid.cvm),] # trial
  
  ## start trial ##
  pred$max[pred$max < 0] <- 0
  pred$max[pred$max > 1] <- 1
  fit$max.cvm <- .loss(y=z,fit=pred$max,family="binomial",type.measure=type.measure)[[1]]
  fit$max.min <- fit$max[which.min(fit$max.cvm)]
  #graphics::plot(x=fit$max,y=fit$max.cvm)
  ## end trial ##
  
  fit$cutoff <- cutoff
  fit$sd.y <- stats::sd(y)

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


# Consider predicting: linear predictor, probability, odds, log(odds)
predict.bilasso <- function(x,newx,type="response"){
  
  if(type!="response"){stop("Invalid type.",call.=FALSE)}
  
  # predicted values - gaussian
  s <- x$gaussian$lambda.min
  pred_y <- as.numeric(stats::predict(object=x$gaussian,newx=newx,s=s,type=type))
  
  # predicted values - binomial
  s <- x$binomial$lambda.min
  pred_z <- as.numeric(stats::predict(object=x$binomial,newx=newx,s=s,type=type))
  
  # gaussian
  #gaussian <- ((pred_y-x$cutoff)/max(abs(pred_y-x$cutoff))+1)/2 # old
  gaussian <- stats::pnorm(q=pred_y,mean=x$cutoff,sd=x$sd.y)
  if(any((pred_y>=x$cutoff)!=(gaussian>=0.5))){
    stop("Wrong check sum.",call.=FALSE)
  }
  if(any(round(gaussian)!=1*(pred_y > x$cutoff))){
    stop("Not compatible.",call.=FALSE)
  }
  if(any(gaussian<0|gaussian>1)){
    stop("unit interval",call.=FALSE)
  }
  
  # binomial
  binomial <- pred_z
  
  # pi-model
  pi <- x$pi.min*gaussian + (1-x$pi.min)*binomial
  if(any((gaussian <= pi) != (pi < binomial))){ # check this
    warning("consistency",call.=FALSE) # check why this happens
  }
  
  # base-model
  #base <- 1/(1+x$base.min^(x$cutoff-pred_y)) # old trial
  base <- stats::pnorm(q=pred_y,mean=x$cutoff,sd=x$base.min) # new trial
  
  # # grid
  #cont <- stats::pnorm(q=pred_y,mean=x$cutoff,sd=x$grid.min$base)
  #grid <- x$grid.min$pi*cont + (1-x$grid.min$pi)*pred_z
  
  # # trial
  #cont <- stats::pnorm(q=pred_y,mean=x$cutoff,sd=x$base.min)
  #trial <- x$pi.min*cont + (1-x$pi.min)*pred_z
  
  ## trial
  #p0 <- stats::pnorm(q=pred_y,mean=x$cutoff,sd=x$grid.min$sd0)
  #p1 <- stats::pnorm(q=pred_y,mean=x$cutoff,sd=x$grid.min$sd1)
  #trial <- ifelse(pred_y<x$cutoff,p0,p1)
  
  trial <- ((pred_y-x$cutoff)/x$max.min + 1)/2
  
  frame <- data.frame(gaussian=gaussian,binomial=binomial,pi=pi,base=base,trial=trial)
  return(frame)
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
bilasso_compare <- function(y,cutoff,X,type.measure="deviance",res=100){
  
  z <- 1*(y > cutoff)
  fold <- palasso:::.folds(y=z,nfolds=5)
  
  cols <- c("gaussian","binomial","pi","base","trial")
  pred <- matrix(data=NA,nrow=length(y),ncol=length(cols),
                 dimnames=list(NULL,cols))
  
  select <- list()
  for(i in sort(unique(fold))){
    fit <- bilasso(y=y[fold!=i],cutoff=cutoff,X=X[fold!=i,],type.measure=type.measure,res=res)
    
    #gaussian <- 1*(stats::predict(object=fit$gaussian,
    #                          newx=X[fold==i,],
    #                          s=fit$gaussian$lambda.min,
    #                          type="response") > cutoff)
    #binomial <- stats::predict(object=fit$binomial,
    #                          newx=X[fold==i,],
    #                          s=fit$binomial$lambda.min,
    #                          type="response")
  #  
  #  pred[fold==i,"gaussian"] <- gaussian
  #  pred[fold==i,"binomial"] <- binomial
  #  pred[fold==i,"mixed"] <- fit$pi.min*pred[fold==i,"gaussian"] + (1-fit$pi.min)*pred[fold==i,"binomial"]
   
    temp <- colasso:::predict.bilasso(fit,newx=X[fold==i,])
    model <- colnames(pred)
    for(j in seq_along(model)){
      pred[fold==i,model[j]] <- temp[[model[j]]]
    }
    
    #pred[fold==i,"gaussian"] <- temp$gaussian
    #pred[fold==i,"binomial"] <- temp$binomial
    #pred[fold==i,"mixed"] <- temp$mixed
    #pred[fold==i,"extra"] <- temp$extra
    #pred[fold==i,"grid"] <- temp$grid
     
  }
  
  type <- c("deviance","class","mse","mae","auc")
  loss <- lapply(X=type,FUN=function(x) .loss(y=z,fit=pred,family="binomial",type.measure=x,foldid=fold)[[1]])
  names(loss) <- type
  
  #loss <- list()
  #loss$deviance <- .loss(y=z,fit=pred,family="binomial",type.measure="deviance")[[1]]
  #loss$class <- .loss(y=z,fit=pred,family="binomial",type.measure="class")[[1]]
  #loss$mse <- .loss(y=z,fit=pred,family="binomial",type.measure="mse")[[1]]
  #loss$mae <- .loss(y=z,fit=pred,family="binomial",type.measure="mae")[[1]]
  #loss$auc <- .loss(y=z,fit=pred,family="binomial",type.measure="auc",foldid=fold)[[1]]
  
  return(loss)
}


#' @export
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
#' @examples
#' NA
#' 
.check <- function(x,type,miss=FALSE,min=NULL,max=NULL,values=NULL,inf=FALSE){
  name <- deparse(substitute(x))
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

}
