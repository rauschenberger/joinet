
# data simulation
list <- cornet:::.simulate(n=100,p=200)
y <- list$y; X <- list$X

# penalised regression
cutoff <- 1
foldid <- cornet:::.folds(y=y>cutoff,nfolds=10)
fit <- cornet::cornet(y=y,cutoff=cutoff,X=X,foldid=foldid)
net <- list()
net$gaussian <- glmnet::cv.glmnet(y=y,x=X,family="gaussian",foldid=foldid)
net$binomial <- glmnet::cv.glmnet(y=y>cutoff,x=X,family="binomial",foldid=foldid)

for(dist in c("gaussian","binomial")){
  
  testthat::test_that("cross-validated loss",{
    a <- fit[[dist]]$sigma.cvm
    b <- net[[dist]]$cvm
    diff <- abs(a[seq_along(b)]-b)
    testthat::expect_true(all(diff<1e-06))
  })
  
  testthat::test_that("optimal lambda",{
    a <- fit[[dist]]$lambda.min
    b <- net[[dist]]$lambda.min
    testthat::expect_true(a==b)
  })
  
  testthat::test_that("lambda sequence",{
    a <- fit[[dist]]$lambda
    b <- net[[dist]]$lambda
    testthat::expect_true(all(a[seq_along(b)]==b))
  })
  
  testthat::test_that("predicted values",{
    a <- stats::predict(object=fit[[dist]],newx=X)
    b <- stats::predict(object=net[[dist]]$glmnet.fit,newx=X)
    testthat::expect_true(all(a==b))
  })
  
}

testthat::test_that("predicted values (logistic)",{
  a <- cornet:::predict.cornet(object=fit,newx=X)$binomial
  b <- as.numeric(stats::predict(object=net$binomial,newx=X,s="lambda.min",type="response"))
  testthat::expect_true(all(a==b))
})


