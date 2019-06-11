
#--- Initialisation ---

set.seed(1)
n <- 30; q <- 3; p <- 20
Y <- matrix(c(rnorm(n),rbinom(n,size=1,prob=0.5),rpois(n,lambda=4)),
            nrow=n,ncol=q)
X <- matrix(rnorm(n*p),nrow=n,ncol=p)
family <- c("gaussian","binomial","poisson")
foldid <- palasso:::.folds(y=Y[,2],nfolds=5)

for(alpha in c(0.05,0.95)){
  
  object <- mixnet::mixnet(Y=Y,X=X,family=family,alpha.base=alpha,foldid=foldid)
  
  glmnet <- list()
  for(i in seq_len(q)){
    glmnet[[i]] <- glmnet::cv.glmnet(x=X,y=Y[,i],family=family[i],alpha=alpha,foldid=foldid,)
  }
  
  #--- Equality glmnet and mixnet ---
  
  testthat::test_that("lambda: glmnet = mixnet",{
    for(i in seq_len(q)){
      a <- glmnet[[i]]$lambda
      b <- object$base[[i]]$lambda
      max <- min(length(a),length(b))
      cond <- all(a[seq_len(max)]==b[seq_len(max)])
      testthat::expect_true(cond)
    }
  })
  
  testthat::test_that("lambda.min: glmnet = mixnet",{
    for(i in seq_len(q)){
      a <- glmnet[[i]]$lambda.min
      b <- object$base[[i]]$lambda.min
      cond <- (a==b)
      testthat::expect_true(cond)
    }
  })
  
  testthat::test_that("cvm: glmnet = mixnet",{
    for(i in seq_len(q)){
      a <- glmnet[[i]]$cvm
      b <- object$base[[i]]$cvm
      max <- min(length(a),length(b))
      cond <- all(abs(a[seq_len(max)]-b[seq_len(max)])<1e-06)
      testthat::expect_true(cond)
    }
  })
  
  testthat::test_that("glmnet.fit: glmnet = mixnet",{
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
  
  #--- Coherence mixnet ---
  
  testthat::test_that("predict: glmnet = mixnet",{
    a <- mixnet:::predict.mixnet(object=object,newx=X)$base
    for(i in seq_len(q)){
      b <- stats::predict(object=glmnet[[i]],newx=X,type="response",s="lambda.min")
      cond <- all(a[,i]==b)
      testthat::expect_true(cond)
    }
  })
  
  #--- Equivalence stacking and pooling ---
  
  testthat::test_that("stacking = pooling",{
    pred0 <- mixnet:::predict.mixnet(object,newx=X)$meta
    coef <- mixnet:::coef.mixnet(object)
    pred1 <- matrix(data=NA,nrow=n,ncol=q)
    for(i in seq_len(q)){
      pred1[,i] <- mixnet:::.mean.function(coef$alpha[i] + X %*% coef$beta[,i],family=family[i])
    }
    cond <- all(abs(pred0-pred1)<1e-06)
    testthat::expect_true(cond)
  })
  
}
