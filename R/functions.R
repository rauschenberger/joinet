
#' @export
#' @aliases colasso-package
#' @title
#' colasso
#' 
#' @description
#' Implements penalised regression with response moderation.
#'  
#' @param y
#' response\strong{:}
#' vector of length \eqn{n}
#' 
#' @param Y
#' response\strong{:}
#' matrix with \eqn{n} rows and \eqn{p} columns,
#' or vector of length \eqn{n} (see details)
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
#' @param family
#' see glmnet
#' 
#' @param type.measure
#' see glmnet
#' 
#' @examples
#' n <- 100; p <- 20; q <- 10
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Y <- matrix(rnorm(n*q),nrow=n,ncol=q)
#' #y <- rbinom(n=n,size=1,prob=0.2)
#' y <- rnorm(n=n)
#' test <- colasso(y=y,Y=Y,X=X)
#' 
colasso <- function(y,Y,X,alpha=1,nfolds=10,family="gaussian",type.measure="deviance"){
  
  # properties
  n <- nrow(X); p <- ncol(X)
  if(!family %in% c("gaussian","poisson","binomial")){
    stop("Family not implemented.")
  }
  if(length(y)!=n){stop("sample size")}
  #foldid <- sample(x=rep(x=seq_len(nfolds),length.out=n))
  foldid <- palasso:::.folds(y=y,nfolds=nfolds)
  
  # weights
  pi <- seq(from=0,to=1,by=0.2) # adapt this

  # model fitting
  fit <- list()
  ym <- colasso::colasso_moderate(Y=Y) # trial
  for(i in seq_along(pi)){
    weights <- rep(c(1-pi[[i]],pi[[i]]),each=n)
    fit[[i]] <- glmnet::glmnet(y=c(y,ym),x=rbind(X,X),weights=weights,family=family,alpha=alpha)
  }
  names(fit) <- paste0("pi",pi)
  
  # inner cross-validation
  pred <- lapply(pi,function(x) matrix(data=NA,nrow=length(y),ncol=100))
  for(k in unique(foldid)){
    y0 <- y[foldid!=k]
    y1 <- y[foldid==k]
    Y0 <- Y[foldid!=k,,drop=FALSE]
    Y1 <- Y[foldid==k,,drop=FALSE]
    X0 <- X[foldid!=k,,drop=FALSE]
    X1 <- X[foldid==k,,drop=FALSE]
    
    y0m <- colasso_moderate(Y=Y0)
    for(i in seq_along(pi)){
      weights <- rep(c(1-pi[[i]],pi[[i]]),each=sum(foldid!=k)) # trial
      glmnet <- glmnet::glmnet(y=c(y0,y0m),x=rbind(X0,X0),weights=weights,family=family,alpha=alpha)
      temp <- stats::predict(object=glmnet,newx=X1,type="response",s=fit[[i]]$lambda)
      pred[[i]][foldid==k,seq_len(ncol(temp))] <- temp
    }
  }
  
  # loss sequence 
  for(i in seq_along(pi)){
    # WATCH OUT: adapt to all loss fuctions
    #fit[[i]]$cvm <- apply(X=pred[[i]],MARGIN=2,FUN=function(x) mean((x-y)^2))
    fit[[i]]$cvm <- palasso:::.loss(y=y,fit=pred[[i]],family=family,type.measure=type.measure,foldid=foldid)[[1]]
    # WATCH OUT: minimise or maximise
    if(type.measure=="AUC"){
      fit[[i]]$lambda.min <- fit[[i]]$lambda[which.max(fit[[i]]$cvm)]
    } else {
      fit[[i]]$lambda.min <- fit[[i]]$lambda[which.min(fit[[i]]$cvm)]
    }
  }
  
  # loss sequence
  #cvm <- palasso:::.loss(y=y,fit=pred,family=family,type.measure=type.measure,foldid=foldid)
  
  # optimisation
  #model <- .extract(fit=fit.full,lambda=lambda,cvm=cvm,type.measure=args$type.measure)
  
  # selection
  cvm <- sapply(fit,function(x) x$cvm[which(x$lambda==x$lambda.min)])
  if(type.measure=="AUC"){
    sel <- which.max(cvm)
  } else {
    sel <- which.min(cvm)
  }
  fit[[length(pi)+1]] <- fit[[sel]]
  
  #graphics::plot(cvm); graphics::abline(v=sel,lty=2)
  names(fit) <- c("glmnet",paste0("pi",pi[-1]),"conet")
  return(fit)
}




#' @export
#' @title
#' moderated response
#'
#' @description
#' This function ...
#'
#' @inheritParams colasso
#' 
#' @param ...
#' further arguments (currently not implemented)
#' vector with entries between \eqn{0} and \eqn{1} (rename argument)
#'
#'
#'
#' @examples
#' NA
colasso_moderate <- function(Y,...){
  # (most basic version possible)
  y <- rowMeans(Y)
  y <- apply(Y,1,stats::median)
  if(all(y %in% c(0,0.5,1))){
    y[y==0.5] <- 1
    warning("Invalid unless binomial family.")
  }
  return(y)
}

#' @export
#' @title
#' simulate data
#' 
#' @description
#' This function ...
#'  
#' @param n
#' sample size
#' 
#' @param p
#' number of covariates
#' 
#' @param cor
#' correlation structure
#' 
#' @param plot
#' logical
#' 
#' @param family
#' character
#' 
#' @examples
#' # CONTINUE HERE
#' 
colasso_simulate <- function(n=100,p=500,cor="constant",family="gaussian",plot=TRUE){
    # correlation matrix
    if(cor=="none"){
        Sigma <- matrix(data=0,nrow=p,ncol=p)
    } else if(cor=="constant"){
        Sigma <- matrix(data=0.05,nrow=p,ncol=p)
    } else if(cor=="autoregressive"){
        # adjust 0.9 to p, such that mean(Sigma)=0.05
        # sum(2*(p-seq_len(p)+1)*0.9^seq_len(p))/(p*p)
        Sigma <- 0.9^abs(col(diag(p))-row(diag(p)))
    } else if(cor=="unstructured"){
        Sigma <- matrix(data=stats::rbeta(n=p,shape1=0.05,shape2=1),nrow=p,ncol=p)
    }
    diag(Sigma) <- 1
    
    X <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma)
    stats::median(abs(as.numeric(cor(X))))

    # non-sparse effects
    #beta <- stats::rnorm(n=p,mean=0,sd=1)
    #mu <- X %*% beta
    #y <- stats::rnorm(n=n,mean=mu)

    # sparse effects
    beta <- rep(x=1,times=p) 
    beta[stats::rbinom(n=p,size=1,prob=0.95)==1] <- 0
    mu <- X %*% beta
    
    if(family=="gaussian"){
      Y <- replicate(n=10,expr=stats::rnorm(n=n,mean=mu,sd=10))
    } else if(family=="binomial"){
      prob <- exp(mu)/(1+exp(mu))
      Y <- replicate(n=10,expr=stats::rbinom(n=n,size=1,prob=prob))
    } else if(family=="poisson"){
      lambda <- exp(mu)
      Y <- replicate(n=10,expr=stats::rpois(n=n,lambda=lambda))
    } else if(family=="cox"){
      warning("Cox regression not implemented!")
    }
    
    #y <- stats::rnorm(n=n,mean=mu,sd=10)
    y <- Y[,1]

    # predictivity -----------------------------------------------------------------
    if(plot){
        graphics::par(mar=c(3,3,1,1))
        test <- glmnet::cv.glmnet(x=X,y=y)
        graphics::plot(x=log(test$lambda),y=test$cvm)
        graphics::abline(h=test$cvm[test$lambda==max(test$lambda)],lty=2)
    }
        
    return(list(y=y,Y=Y,X=X))
}

#' @export
#' @title
#' External cross-validation
#' 
#' @description
#' This function ...
#'  
#' @param y
#' response vector
#' 
#' @param Y
#' response moderation matrix
#' 
#' @param X
#' covariates
#' 
#' @param plot
#' logical
#' 
#' @param nfolds.int
#' internal folds
#' 
#' @inheritParams colasso
#' 
#' @examples
#' NA
#'
colasso_compare <- function(y,Y,X,plot=TRUE,nfolds.int=10,family="gaussian",type.measure="deviance"){
    
    fold <- sample(x=rep(x=seq_len(5),length.out=length(y)))
    pred <- matrix(data=NA,nrow=length(y),ncol=8)
    select <- list()
    for(i in sort(unique(fold))){
        cat("i =",i,"\n")
        fit <- colasso(y=y[fold!=i],Y=Y[fold!=i,],X=X[fold!=i,],alpha=1,nfolds=nfolds.int,type.measure=type.measure)
        for(j in seq_along(fit)){
            pred[fold==i,j] <- glmnet::predict.glmnet(object=fit[[j]],
                                                  newx=X[fold==i,],
                                                  s=fit[[j]]$lambda.min,
                                                  type="response")
        }
        select[[i]] <- lapply(fit,function(x) which(x$beta[,x$lambda==x$lambda.min]!=0))
        pred[fold==i,8] <- mean(y[fold!=i]) # intercept-only model
    }
    colnames(pred) <- c(names(fit),"intercept")
    
    if(family=="gaussian" & type.measure=="deviance" | family=="binomial" & type.measure=="mse"){
      loss <- apply(X=pred,MARGIN=2,FUN=function(x) sum((y-x)^2))
    } else if(family=="binomial" & type.measure=="class"){
      loss <- apply(X=pred,MARGIN=2,FUN=function(x) mean(y!=x))
    } else {
      warning("Implement other loss functions!")
    }
    
    ### start temporary ###
    #stability <- numeric()
    #for(k in seq_along(fit)){
    #    matrix <- matrix(data=NA,nrow=5,ncol=5)
    #    for(i in seq_len(5)){
    #        for(j in seq_len(5)){
    #            a <- select[[i]][[k]]
    #            b <- select[[j]][[k]]
    #            matrix[i,j] <- length(intersect(a,b))/length(union(a,b))
    #        }
    #    }
    #    diag(matrix) <- NA
    #    stability[k] <- mean(matrix,na.rm=TRUE)
    #}
    #cat(stability)
    ### end temporary ###
    
    if(plot){
        graphics::par(mar=c(3,3,1,1))
        col <- rep(x=0,times=length(loss)-1)
        col[1] <- col[length(col)] <- 1
        graphics::plot.new()
        graphics::plot.window(xlim=c(1,length(loss)-1),ylim=range(loss))
        graphics::axis(side=1,at=seq_len(length(loss)-1),labels=names(loss)[-length(loss)])
        graphics::axis(side=2)
        graphics::box()
        graphics::points(y=loss[-length(loss)],
                       x=seq_len(length(loss)-1),
                       col=col+1,pch=col)
        graphics::abline(v=c(1.5,length(loss)-1.5),lty=2)
        graphics::grid()
        graphics::abline(h=loss[length(loss)],lty=2,col="red")
    }
    
    return(loss)
}
