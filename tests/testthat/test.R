
set.seed(1)
n <- 30; q <- 3; p <- 20
Y <- matrix(c(rnorm(n),rbinom(n,size=1,prob=0.5),rpois(n,lambda=4)),
            nrow=n,ncol=q)
X <- matrix(rnorm(n*p),nrow=n,ncol=p)
family <- c("gaussian","binomial","poisson")
foldid <- palasso:::.folds(y=Y[,2],nfolds=5)
object <- mixnet::mixnet(Y=Y,X=X,family=family,foldid=foldid)

glmnet <- list()
for(i in seq_len(q)){
  glmnet[[i]] <- glmnet::cv.glmnet(x=X,y=Y[,i],family=family[i],foldid=foldid)
}

testthat::test_that("lambda: glmnet = stacknet",{
  for(i in seq_len(q)){
    a <- glmnet[[i]]$lambda
    b <- object$base[[i]]$lambda
    max <- min(length(a),length(b))
    cond <- all(a[seq_len(max)]==b[seq_len(max)])
    testthat::expect_true(cond)
  }
})

testthat::test_that("cvm: glmnet = stacknet",{
  for(i in seq_len(q)){
    a <- glmnet[[i]]$cvm
    b <- object$base[[i]]$cvm
    max <- min(length(a),length(b))
    cond <- all(abs(a[seq_len(max)]-b[seq_len(max)])<1e-06)
    testthat::expect_true(cond)
  }
})

testthat::test_that("glmnet.fit: glmnet = stacknet",{
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
