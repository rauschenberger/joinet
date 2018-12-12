
#' @export
#' @title
#' bilasso
#' 
#' @description
#' Implements penalised regression with response duality.
#' \code{pi=0} represents binomial regression,
#' \code{pi=1} represents linear regression
#'  
#' @param y
#' continuous response\strong{:}
#' vector of length \eqn{n}
#' 
#' @param z
#' binary response\strong{:}
#' vector of length \eqn{n}
#' 
#' @param cutoff
#' value between \code{min(y)} and \code{max(y)}
#' 
#' @param X
#' covariates\strong{:}
#' matrix with \eqn{n} rows (samples) and \eqn{p} columns (variables)
#' 
#' @param alpha
#' elastic net parameter\strong{:}
#' numeric between \eqn{0} and \eqn{1};
#' \eqn{alpha=1} for lasso,
#' \eqn{alpha=0} for ridge
#' 
#' @param nfolds
#' number of folds
#' 
#' @examples
#' NA
#' 
bilasso <- function(y,cutoff,X,alpha=1,nfolds=10){
  
  z <- 1*(y > cutoff)
  # alpha <- 1; nfolds <- 10
  
  # properties
  n <- nrow(X); p <- ncol(X)
  if(length(y)!=n){stop("sample size")}
  foldid <- palasso:::.folds(y=z,nfolds=nfolds)
  if(cutoff < min(y) | max(y) < cutoff){stop("Cutoff outside.")}
  
  # model fitting
  fit <- list()
  fit$gaussian <- glmnet::glmnet(y=y,x=X,family="gaussian",alpha=alpha)
  fit$binomial <- glmnet::glmnet(y=z,x=X,family="binomial",alpha=alpha)
  
  # weights
  fit$pi <- seq(from=0,to=1,length.out=101) # adapt this
  
  # inner cross-validation
  pred_y <- pred_z <- matrix(data=NA,nrow=length(y),ncol=100)
  pred <- matrix(data=NA,nrow=length(y),ncol=length(fit$pi))
  for(k in unique(foldid)){
    y0 <- y[foldid!=k]
    y1 <- y[foldid==k]
    z0 <- z[foldid!=k]
    z1 <- z[foldid==k]
    X0 <- X[foldid!=k,,drop=FALSE]
    X1 <- X[foldid==k,,drop=FALSE]
    
    foldid_int <- palasso:::.folds(y=z0,nfolds=nfolds)
    
    net_y <- glmnet::glmnet(y=y0,x=X0,family="gaussian",alpha=alpha)
    net_z <- glmnet::glmnet(y=z0,x=X0,family="binomial",alpha=alpha)
    
    temp_y <- stats::predict(object=net_y,newx=X1,type="response",s=fit$gaussian$lambda)
    cvm_y <- .loss(y=y1,fit=temp_y,family="gaussian",type.measure="deviance")[[1]]
    sel_y <- which.min(cvm_y)
    pred_y[foldid==k,seq_len(ncol(temp_y))] <- temp_y
    
    temp_z <- stats::predict(object=net_z,newx=X1,type="response",s=fit$binomial$lambda)
    cvm_z <- .loss(y=z1,fit=temp_z,family="binomial",type.measure="deviance")[[1]]
    sel_z <- which.min(cvm_z)
    pred_z[foldid==k,seq_len(ncol(temp_z))] <- temp_z
    
    for(i in seq_along(fit$pi)){
      pred[foldid==k,i] <- fit$pi[i]*(temp_y[,sel_y] > cutoff) + (1-fit$pi[i])*temp_z[,sel_z]
    }
  
  }
  
  fit$gaussian$cvm <- .loss(y=y,fit=pred_y,family="gaussian",type.measure="deviance")[[1]]
  fit$gaussian$lambda.min <- fit$gaussian$lambda[which.min(fit$gaussian$cvm)]
  fit$binomial$cvm <- .loss(y=z,fit=pred_z,family="binomial",type.measure="deviance")[[1]]
  fit$binomial$lambda.min <- fit$binomial$lambda[which.min(fit$binomial$cvm)]
  
  fit$cvm <- .loss(y=z,fit=pred,family="binomial",type.measure="deviance")[[1]]
  sel <- which.min(fit$cvm)
  fit$pi.min <- fit$pi[sel]
  
  class(fit) <- "bilasso"
  return(fit)
}

bilasso_compare <- function(y,cutoff,X){
  
  z <- 1*(y > cutoff)
  
  fold <- palasso:::.folds(y=z,nfolds=5)
  pred <- matrix(data=NA,nrow=length(y),ncol=3,
                 dimnames=list(NULL,c("gaussian","binomial","mixed")))
  
  select <- list()
  for(i in sort(unique(fold))){
    cat("i =",i,"\n")
    fit <- bilasso(y=y[fold!=i],X=X[fold!=i,],cutoff=cutoff)
    
    gaussian <- 1*(stats::predict(object=fit$gaussian,
                              newx=X[fold==i,],
                              s=fit$gaussian$lambda.min,
                              type="response") > cutoff)
    binomial <- stats::predict(object=fit$binomial,
                              newx=X[fold==i,],
                              s=fit$binomial$lambda.min,
                              type="response")
    
    pred[fold==i,"gaussian"] <- gaussian
    pred[fold==i,"binomial"] <- binomial
    pred[fold==i,"mixed"] <- fit$pi.min*pred[fold==i,"gaussian"] + (1-fit$pi.min)*pred[fold==i,"binomial"]
      
  }
  
  loss <- list()
  loss$deviance <- .loss(y=z,fit=pred,family="binomial",type.measure="deviance")[[1]]
  loss$class <- .loss(y=z,fit=pred,family="binomial",type.measure="class")[[1]]
  loss$mse <- .loss(y=z,fit=pred,family="binomial",type.measure="mse")[[1]]
  loss$mae <- .loss(y=z,fit=pred,family="binomial",type.measure="mae")[[1]]
  loss$auc <- .loss(y=z,fit=pred,family="binomial",type.measure="auc",foldid=fold)[[1]]
  
  return(loss)
}

# Correct this function in the palasso package (search for "# typo").

.loss <- function (y, fit, family, type.measure, foldid = NULL) 
{
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

