
if(!grepl('SunOS',Sys.info()['sysname'])){

#--- Initialisation ---

set.seed(1)
n <- 30; q <- 3; p <- 20
Y <- matrix(c(rnorm(n),rbinom(n,size=1,prob=0.5),rpois(n,lambda=4)),
            nrow=n,ncol=q)
X <- matrix(rnorm(n*p),nrow=n,ncol=p)
family <- c("gaussian","binomial","poisson")
foldid <- palasso:::.folds(y=Y[,2],nfolds=5)

testthat::test_that("link-mean",{
  x <- stats::rnorm(n=100)
  for(family in c("gaussian","binomial","poisson")){
    mean <- joinet:::.mean.function(x,family=family)
    link <- joinet:::.link.function(mean,family=family)
    cond <- all(abs(x-link)<1e-06)
    testthat::expect_true(cond)
  }
})

for(alpha in c(0.05,0.95)){
  
  object <- joinet::joinet(Y=Y,X=X,family=family,alpha.base=alpha,foldid=foldid)
  
  glmnet <- list()
  for(i in seq_len(q)){
    glmnet[[i]] <- glmnet::cv.glmnet(x=X,y=Y[,i],family=family[i],alpha=alpha,foldid=foldid,)
  }
  
  #--- Equality glmnet and joinet ---
  
  testthat::test_that("lambda: glmnet = joinet",{
    for(i in seq_len(q)){
      a <- glmnet[[i]]$lambda
      b <- object$base[[i]]$lambda
      max <- min(length(a),length(b))
      cond <- all(a[seq_len(max)]==b[seq_len(max)])
      testthat::expect_true(cond)
    }
  })
  
  testthat::test_that("lambda.min: glmnet = joinet",{
    for(i in seq_len(q)){
      a <- glmnet[[i]]$lambda.min
      b <- object$base[[i]]$lambda.min
      cond <- (a==b)
      testthat::expect_true(cond)
    }
  })
  
  testthat::test_that("cvm: glmnet = joinet",{
    for(i in seq_len(q)){
      a <- glmnet[[i]]$cvm
      b <- object$base[[i]]$cvm
      max <- min(length(a),length(b))
      cond <- all(abs(a[seq_len(max)]-b[seq_len(max)])<1e-06)
      testthat::expect_true(cond)
    }
  })
  
  testthat::test_that("glmnet.fit: glmnet = joinet",{
    for(i in seq_len(q)){
      a <- glmnet[[i]]$glmnet.fit
      b <- object$base[[i]]$glmnet.fit
      names <- setdiff(x=names(a),y="call")
      for(j in names){
        cond <- all(a[[j]]==b[[j]])
        testthat::expect_true(cond)
      }
    }
  })
  
  #--- Coherence joinet ---
  
  testthat::test_that("predict: glmnet = joinet",{
    a <- joinet:::predict.joinet(object=object,newx=X)$base
    for(i in seq_len(q)){
      b <- stats::predict(object=glmnet[[i]],newx=X,type="response",s="lambda.min")
      cond <- all(a[,i]==b)
      testthat::expect_true(cond)
    }
  })
  
  #--- Equivalence stacking and pooling ---
  
  testthat::test_that("stacking = pooling",{
    pred0 <- joinet:::predict.joinet(object,newx=X)$meta
    coef <- joinet:::coef.joinet(object)
    pred1 <- matrix(data=NA,nrow=n,ncol=q)
    for(i in seq_len(q)){
      pred1[,i] <- joinet:::.mean.function(coef$alpha[i] + X %*% coef$beta[,i],family=family[i])
    }
    cond <- all(abs(pred0-pred1)<1e-06)
    testthat::expect_true(cond)
  })
  
}

}

testthat::test_that("argument weight",{
  weight <- matrix(stats::rbinom(n=p*q,size=1,prob=0.3),nrow=p,ncol=q)
  object <- joinet(Y=Y,X=X,weight=weight,alpha.base=0,family=family)
  coef <- sapply(object$base, function(x) as.numeric(stats::coef(object = x$glmnet.fit, s = x$lambda.min))[-1])
  cond <- all((weight==0)==(coef==0))
  testthat::expect_true(cond)
})

testthat::test_that("argument sign",{
  sign <- matrix(stats::rbinom(n=q*q,size=2,prob=0.3)-1,nrow=q,ncol=q) # asymmetric
  diag(sign) <- 1
  object <- joinet(Y=Y,X=X,sign=sign,alpha.meta=0,family=family)
  weights <- weights(object)[-1,]
  cond <- all(weights[sign==-1]<=0) & all(weights[sign==0]==0) & all(weights[sign==1]>=0)
  testthat::expect_true(cond)
})
