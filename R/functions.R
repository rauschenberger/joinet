

#' @export
#' @aliases colasso-package
#' @title
#' colasso
#' 
#' @description
#' This function ...
#'  
#' @param y
#' response\strong{:}
#' vector of length \eqn{n}
#' 
#' @param X
#' covariates\strong{:}
#' matrix with \eqn{n} rows (samples) and \eqn{p} columns (variables)
#' 
#' @param nfolds
#' number of folds
#' 
#' @param alpha
#' elastic net parameter
#' 
#' @examples
#' n <- 100; p <- 20
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' #y <- rbinom(n=n,size=1,prob=0.2)
#' y <- rnorm(n=n)
#' #y[1] <- 0.5
#' #a <- glmnet::glmnet(y=y,x=x,family="binomial")
#' #b <- stats::glm(y~x,family="binomial")
colasso <- function(y,X,nfold=10,alpha=1,nfolds=10){
    
    # properties
    n <- nrow(X); p <- ncol(X)
    if(length(y)!=n){stop("sample size")}
    foldid <- sample(x=rep(x=seq_len(nfolds),length.out=n))
    pi <- seq(from=0,to=0.5,by=0.1) # adapt this

    # model fitting
    fit <- list()
    Y <- colasso_moderate(y=y,X=X,pi=pi)
    for(j in seq_along(pi)){
        fit[[j]] <- glmnet::glmnet(y=Y[,j],x=X,alpha=alpha)
    }
    
    # inner cross-validation
    pred <- lapply(pi,function(x) matrix(data=NA,nrow=length(y),ncol=100))
    for(k in sort(unique(foldid))){
        y0 <- y[foldid!=k]
        y1 <- y[foldid==k]
        X0 <- X[foldid!=k,,drop=FALSE]
        X1 <- X[foldid==k,,drop=FALSE]
        
        Y0 <- colasso_moderate(y=y0,X=X0,pi=pi) 
        for(j in seq_along(pi)){
            glmnet <- glmnet::glmnet(y=Y0[,j],x=X0,alpha=alpha) # was lambda=lambda
            temp <- stats::predict(object=glmnet,newx=X1,type="response",s=fit[[j]]$lambda)
            pred[[j]][foldid==k,seq_len(ncol(temp))] <- temp
        }
    }
    
    # loss sequence
    for(i in seq_along(pi)){
        fit[[i]]$cvm <- apply(X=pred[[i]],MARGIN=2,FUN=function(x) mean((x-y)^2))
        fit[[i]]$lambda.min <- fit[[i]]$lambda[which.min(fit[[i]]$cvm)]
    }
    
    # selection
    cvm <- sapply(fit,function(x) x$cvm[which(x$lambda==x$lambda.min)])
    sel <- which.min(cvm) # BUT "sel <- which.max(cvm)" FOR AUC!
    fit[[length(pi)+1]] <- fit[[sel]]
    
    #graphics::plot(cvm); graphics::abline(v=sel,lty=2)
    names(fit) <- c("standard",paste0("pi",pi[-1]),"select")
    return(fit)
}




# colasso <- function(y,X,alpha=1,nfolds=10){
#     n <- nrow(X); p <- ncol(X)
#     if(length(y)!=n){stop("sample size")}
#     foldid <- sample(x=rep(x=seq_len(nfolds),length.out=n))
#     pi <- seq(from=0,to=0.5,by=0.1) # adapt this
#     lambda <- glmnet::glmnet(y=y,x=X,alpha=alpha)$lambda
#     
#     pred <- lapply(pi,function(x) matrix(data=NA,nrow=length(y),ncol=100))
#     for(k in sort(unique(foldid))){
#         y0 <- y[foldid!=k]
#         y1 <- y[foldid==k]
#         X0 <- X[foldid!=k,,drop=FALSE]
#         X1 <- X[foldid==k,,drop=FALSE]
#         
#         Y0 <- colasso_moderate(y=y0,X=X0,pi=pi) 
#         for(j in seq_along(pi)){
#             glmnet <- glmnet::glmnet(y=Y0[,j],x=X0,alpha=alpha) # was lambda=lambda
#             temp <- stats::predict(object=glmnet,newx=X1,type="response",s=lambda)
#             pred[[j]][foldid==k,seq_len(ncol(temp))] <- temp
#         }
#     }
#     
#     Y <- colasso_moderate(y=y,X=X,pi=pi)
#     fit <- list()
#     for(j in seq_along(pi)){
#         fit[[j]] <- glmnet::glmnet(y=Y[,j],x=X,alpha=alpha) # was lambda=lambda
#         #fit[[j]]$cvm <- palasso::loss.trial(y=y,fit=pred[[j]],family=family,type.measure=type.measure)[[1]]
#         fit[[j]]$cvm <- apply(X=pred[[j]],MARGIN=2,FUN=function(x) mean((x-y)^2))
#         fit[[j]]$lambda.min <- lambda[which.min(fit[[j]]$cvm)]
#     }
#     
#     cvm <- sapply(fit,function(x) x$cvm[which.min(abs(x$lambda-x$lambda.min))[1]])
#     
#     #sel <- which.max(cvm) # only for AUC
#     sel <- which.min(cvm)
#     
#     fit[[length(pi)+1]] <- fit[[sel]]
#     #graphics::plot(cvm); graphics::abline(v=sel,lty=2)
#     names(fit) <- c("standard",paste0("pi",pi[-1]),"select")
#     return(fit)
# }

#' @export
#' @title
#' marginal significance
#' 
#' @description
#' This function ...
#'  
#' @inheritParams colasso
#' 
#' @examples
#' NA
#' 
# obtain p-values --------------------------------------------------------------
# input: y, X; output: p-value
colasso_marginal_significance <- function(y,X){
    x <- vector()
    for(i in seq_len(ncol(X))){
        if(stats::var(X[,i])==0){
            x[i] <- 0
        } else {
            x[i] <- summary(stats::lm(y~X[,i]))$coefficients["X[, i]","Pr(>|t|)"]
        }
    }
    return(x)
}

#' @export
#' @title
#' covariate weights
#' 
#' @description
#' This function ...
#'  
#' @param x
#' \eqn{p}-values\strong{:}
#' vector with entries between \eqn{0} and \eqn{1}
#' 
#' @param max
#' maximum \eqn{p}-value to receive weight one
#' 
#' @param min
#' minimum \eqn{p}-value to receive weight zero
#' 
#' @param version
#' (temporary argument)
#' 
#' @examples
#' # alternative weights: weight <- abs(cor(y,X))
#' 
colasso_covariate_weights <- function(x,max=0.05/length(x),min=0.05,version=1){ # was min=0.2
    if(version==1){
        weight <- pmin(-log10(max),pmax(-log10(x)+log10(min),0))
        weight <- weight/(-log10(max))
        # more compact: (log10(x)-log10(min))/(log10(max)) # (and then truncate)
    } else {
        weight <- pmax(1-(1/min)*x,0) 
    }
    return(weight)
}

#' @export
#' @title
#' weighted correlation
#' 
#' @description
#' This function ...
#'  
#' @inheritParams colasso
#' 
#' @param w
#' covariate weights\strong{:}
#' vector of length \eqn{p}, with entries between \eqn{0} and \eqn{1}
#' 
#' @examples
#' NA
#' 
# weighted covariate correlation -----------------------------------------------
# input: X, weights; output: correlation
colasso_weighted_correlation <- function(X,w=NULL){
    n <- nrow(X); p <- ncol(X)
    if(is.null(w)){w <- rep(x=1,times=p)}
    mx <- rep(x=NA,times=p)
    for(i in seq_len(p)){
        mx[i] <- sum(w*X[,i])/sum(w)
    }
    sx <- cor <- matrix(NA,nrow=p,ncol=p)
    diag(cor) <- 1
    for(i in seq(from=1,to=p,by=1)){
        sx[i,i] <- sum(w*(X[,i]-mx[i])^2)/sum(w)
    }
    for(i in seq(from=1,to=p,by=1)){
        for(j in seq(from=i,to=p,by=1)){
            sx[i,j] <- sx[j,i] <- sum(w*(X[,i]-mx[i])*(X[,j]-mx[j]))/sum(w)
            cor[i,j] <- cor[j,i] <- sx[i,j]/sqrt(sx[i,i]*sx[j,j])
        }
    }
    return(cor)
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
#' @param pi
#' vector with entries between \eqn{0} and \eqn{1} (rename argument)
#' 
#' @param plot
#' logical
#' 
#' 
#' @examples
#' NA
#' 
colasso_moderate <- function(y,X,pi,plot=FALSE){
    pvalue <- colasso_marginal_significance(y=y,X=X) 
    weight <- colasso_covariate_weights(x=pvalue)
    #cor <- abs(colasso_weighted_correlation(t(X),w=weight)) # robust
    cor <- abs(weights::wtd.cors(t(X),weight=weight)) # fast
    Y <- matrix(data=NA,nrow=length(y),ncol=length(pi),dimnames=list(NULL,pi))
    for(i in seq_along(y)){
        for(j in seq_along(pi)){
            Y[i,j] <- (1-pi[j])*y[i]+pi[j]*sum(cor[,i]*y)/sum(cor[,i])
        }
    }
    if(plot){
        #plot(ecdf(x)); abline(a=0,b=1,lty=2)
        #plot(x=x,y=weight,ylim=c(0,1)); abline(v=0.05,lty=2)
        #graphics::image(cor)
        #graphics::image(Y)
    }
    return(Y)
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
#' @examples
#' # CONTINUE HERE
#' 
colasso_simulate <- function(n=100,p=500,cor="constant",plot=TRUE){
    # correlation matrix
    if(cor=="none"){
        Sigma <- matrix(data=0,nrow=p,ncol=p)
    } else if(cor=="constant"){
        Sigma <- matrix(data=0.05,nrow=p,ncol=p)
    } else if(cor=="autoregressive"){
        # adjust 0.9 to p, such that mean(Sigma)=0.05
        #sum(2*(p-seq_len(p)+1)*0.9^seq_len(p))/(p*p)
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
    beta <- rep(1,times=p) 
    beta[stats::rbinom(n=p,size=1,prob=0.97)==1] <- 0
    mu <- X %*% beta
    y <- stats::rnorm(n=n,mean=mu,sd=10)

    # predictivity -----------------------------------------------------------------
    if(plot){
        test <- glmnet::cv.glmnet(x=X,y=y)
        graphics::plot(x=log(test$lambda),y=test$cvm)
        graphics::abline(h=test$cvm[test$lambda==max(test$lambda)],lty=2)
    }
        
    return(list(y=y,X=X))
}




#' @export
#' @title
#' External cross-validation
#' 
#' @description
#' This function ...
#'  
#' @param y
#' response
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
#' @examples
#' NA
#'
colasso_compare <- function(y,X,plot=TRUE,nfolds.int=10){
    
    fold <- sample(x=rep(x=seq_len(5),length.out=length(y)))
    pred <- matrix(data=NA,nrow=length(y),ncol=8)
    for(i in sort(unique(fold))){
        cat("i =",i,"\n")
        fit <- colasso(y=y[fold!=i],X=X[fold!=i,],alpha=1,nfolds=nfolds.int)
        for(j in seq_along(fit)){
            pred[fold==i,j] <- glmnet::predict.glmnet(object=fit[[j]],
                                                  newx=X[fold==i,],
                                                  s=fit[[j]]$lambda.min,
                                                  type="response")
        }
        pred[fold==i,8] <- mean(y[fold!=i]) # intercept-only model
    }
    colnames(pred) <- c(names(fit),"intercept")
    loss <- apply(X=pred,MARGIN=2,FUN=function(x) sum((y-x)^2))
    
    if(plot){
        graphics::par(mar=c(3,3,1,1))
        col <- rep(x=0,times=length(loss)-1)
        col[1] <- col[length(col)] <- 1
        graphics::plot(y=loss[-length(loss)],
                       x=seq_len(length(loss)-1),
                       col=col+1,ylim=range(loss),pch=col)
        graphics::abline(v=c(1.5,length(loss)-1.5),lty=2)
        graphics::grid()
        graphics::abline(h=loss[length(loss)],lty=2,col="red")
    }
    
    return(loss)
}


