#
# libmlts.R - A small library of machine learning methods for time series, v01.08 (2024-04-20)
#
# Pedro Afonso Fernandes, UCP, CLSBE, Lisbon, Portugal (paf@ucp.pt)
#
# This library is free software; you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation. The library is
# distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details: https://www.gnu.org/licenses/.
#


# Load auxiliary packages

library(MASS)
library(textir)
library(gamlr)
library(tree)
library(kernlab)
library(h2o)
library(glmnet)

h2o.init()

### CROSS VALIDATION ###

# cv.ols(x, y, k, z)
#
# Simple k-fold cross validation for a linear regression of time series.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.ols <- function(x, y, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[i:(m+i-1),1]
    xtrain <- x[i:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- lm(ytrain ~ xtrain)
    
    pred[i]  <- xtest %*% coefficients(mod)[-1] + coefficients(mod)[1]
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k

  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.logit(x, y, k, z)
#
# Simple k-fold cross validation for a logistic regression of time series.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.logit <- function(x, y, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[i:(m+i-1),1]
    xtrain <- x[i:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- glm(ytrain ~ xtrain, family = binomial(link = "logit"))
    
    pred[i]  <- 1 / (1 + exp(xtest %*% coefficients(mod)[-1] + coefficients(mod)[1]) )
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.pls(x, y, k)
#
# Simple k-fold cross validation for a boosted Partial Least Squares model of time series.
#
# Returns a list with k OOS predictions, the associated errors (MAE and MSE), AIC (Greene approach) and 
# the number of PLS directions that min(MAE).
#

cv.pls <- function(x, y, k=4){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  q <- ncol(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)

  mae_aux <- 1000000000
  
  for (d in 1:5){
  
    for (i in 1:k){
      
        ytrain <- y[i:(m+i-1),1]
        xtrain <- x[i:(m+i-1),]
        xtest  <- x[m+i,]
    
        mod <- pls(xtrain, ytrain, K=d)
    
        pred[i]  <- predict(mod, newdata=xtest)
    
        error[i] <- y[m+i] - pred[i]
    
    }
  
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k

    aic <- log(sum(error*error)/k) + 2*q/k

    if (mae <= mae_aux){
        
        pred_aux <- pred
        mae_aux  <- mae
        mse_aux  <- mse
        aic_aux  <- aic
        d_aux    <- d
    }
  }
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred_aux, "MAE" = mae_aux, "MSE" = mse_aux, "AIC" = aic_aux, "directions" = d_aux, "cum_error" = cumerr)
  
}


# cv.pls2(x, y, k, d, z)
#
# Simple k-fold cross validation for a boosted Partial Least Squares model of time series with d directions.
#
# Returns a list with k OOS predictions, the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.pls2 <- function(x, y, k=4, d=3, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  q <- ncol(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)

  for (i in 1:k){
  	      
    ytrain <- y[i:(m+i-1),1]
    xtrain <- x[i:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- pls(xtrain, ytrain, K=d)
    
    pred[i]  <- predict(mod, newdata=xtest)
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k

  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}



# cv.lasso(x, y, k, z)
#
# Simple k-fold cross validation for a linear LASSO regression of time series.
#
# Returns a list with k OOS predictions that min(AICc) and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.lasso <- function(x, y, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[i:(m+i-1),1]
    xtrain <- x[i:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- gamlr(xtrain, ytrain, gamma=0)
    
    pred[i]  <- predict(mod, newdata=t(as.matrix(xtest)))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k

  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.lasso.logit(x, y, k, z)
#
# Simple k-fold cross validation for a logistic LASSO regression of time series.
#
# Returns a list with k OOS predictions that min(AICc) and the associated errors (MAE and MSE) and AIC (Greene approach)..
#

cv.lasso.logit <- function(x, y, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[i:(m+i-1),1]
    xtrain <- x[i:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- gamlr(xtrain, ytrain, family = "binomial", gamma=0)
    
    pred[i]  <- predict(mod, newdata=t(as.matrix(xtest)))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.svm(x, y, k, z)
#
# Simple k-fold cross validation for a support vector machine (with polynomial kernel) of time series.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.svm <- function(x, y, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[i:(m+i-1),1]
    xtrain <- x[i:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- ksvm(ytrain ~ xtrain, kernel ="polydot")
    
    pred[i]  <- predict(mod, newdata=t(as.matrix(xtest)))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.rtree(formula, df, md, k, z)
#
# Simple k-fold cross validation for a regression tree of time series. 
#
# The input must be a formula y ~ x, a data frame (df) and the min deviance improvement (md).
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.rtree <- function(formula, df, md=0.01, k=4, z=2){
  
  n <- nrow(df)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    dtrain <- df[i:(m+i-1),]
    dtest  <- df[m+i,]
    
    mod <- tree(formula, data = dtrain, mindev = md)
    
    pred[i]  <- predict(mod, newdata = dtest)
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(df) - 1
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.forest(formula, df, md, b, k, z)
#
# Simple k-fold cross validation for a forest of b regression trees of time series.
#
# The input must be a formula y ~ x, a data frame (df), the min deviance improvement (md) and the number of trees (b) in the forest.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.forest <- function(formula, df, md=0.01, b=5, k=4, z=2){
  
  n <- nrow(df)
  m <- n - k
  
  pred  <- rep(0,k)
  pred2 <- matrix(0,k,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    dtrain <- df[i:(m+i-1),]
    dtest  <- df[(m+1):n,]
    
    mod <- tree(formula, data = dtrain, mindev = md)
    
    pred2[i,]  <- predict(mod, newdata = dtest)
    
  }
  
  for (i in 1:k){
  
    a <- max(1,i-b+1)
    
    pred[i] <- mean(pred2[a:i,i])
  
  }
  
  error <- y[(m+1):n] - pred
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(df) - 1
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.bag.ols(x, y, b, k, z)
#
# Simple k-fold cross validation for a bag of b linear regressions of time series.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.bag.ols <- function(x, y, b=5, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)

  pred  <- rep(0,k)
  pred2 <- matrix(0,k,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[i:(m+i-1),1]
    xtrain <- x[i:(m+i-1),]
    xtest  <- x[(m+1):n,]
    
    mod <- lm(ytrain ~ xtrain)
    
    pred2[i,]  <- t(xtest %*% coefficients(mod)[-1] + coefficients(mod)[1])
    
  }
  
  for (i in 1:k){
  
    a <- max(1,i-b+1)
    
    pred[i] <- mean(pred2[a:i,i])
  
  }
  
  error <- y[(m+1):n] - pred
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.drf(y, df, k, z)
#
# Simple k-fold cross validation for a Distributed Random Forest (DRF). 
#
# The input must be a "string" with the name of the response, a data frame (df) and the number of folds (k).
#
# NB: all columns of df except y are used as covariates. 
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.drf <- function(response, df, k=4, z=2){
  
  n <- nrow(df)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    dtrain <- as.h2o(df[i:(m+i-1),])
    dtest  <- as.h2o(df[m+i,])
    
    mod <- h2o.randomForest(y = response, training_frame = dtrain)
    
    pred[i]  <- as.vector(h2o.predict(mod, newdata = dtest))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(df) - 1
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.gbm(y, df, k, z)
#
# Simple k-fold cross validation for a Gradient Boosting Machine (GBM) of regression trees. 
#
# The input must be a "string" with the name of the response, a data frame (df) and the number of folds (k).
#
# NB: all columns of df except y are used as covariates. 
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.gbm <- function(response, df, k=4, z=2){
  
  n <- nrow(df)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    dtrain <- as.h2o(df[i:(m+i-1),])
    dtest  <- as.h2o(df[m+i,])
    
    mod <- h2o.gbm(y = response, training_frame = dtrain)
    
    pred[i]  <- as.vector(h2o.predict(mod, newdata = dtest))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(df) - 1
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.rulefit(y, df, k, z)
#
# Simple k-fold cross validation for a RuleFit regression of time series. 
#
# The general algorithm fits a tree ensemble to the data, builds a rule ensemble by traversing each tree, 
# evaluates the rules on the data to build a rule feature set, and fits a LASSO model to the rule feature 
# set joined with the original feature set.
#
# The input must be a "string" with the name of the response, a data frame (df) and the number of folds (k).
#
# NB: all columns of df except y are used as covariates. 
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.rulefit <- function(response, df, k=4, z=2){
  
  n <- nrow(df)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    dtrain <- as.h2o(df[i:(m+i-1),])
    dtest  <- as.h2o(df[m+i,])
    
    mod <- h2o.rulefit(y = response, training_frame = dtrain)
        
    pred[i]  <- as.vector(h2o.predict(mod, newdata = dtest))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(df) - 1
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.dnn(y, df, k, z)
#
# Simple k-fold cross validation for a Deep Neural Network (multi-layer perceptron) of time series. 
#
# The input must be a "string" with the name of the response, a data frame (df) and the number of folds (k).
#
# NB: all columns of df except y are used as covariates. 
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.dnn <- function(response, df, k=4, z=2){
  
  n <- nrow(df)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    dtrain <- as.h2o(df[i:(m+i-1),])
    dtest  <- as.h2o(df[m+i,])
    
    mod <- h2o.deeplearning(y = response, training_frame = dtrain)
    
    pred[i]  <- as.vector(h2o.predict(mod, newdata = dtest))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(df) - 1
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.ridge(x, y, k, z)
#
# Simple k-fold cross validation for a ridge regression of time series.
#
# Returns a list with k OOS predictions for lambda/s=0.001 and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.ridge <- function(x, y, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[i:(m+i-1),1]
    xtrain <- x[i:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- glmnet(xtrain, ytrain, alpha=0)
    
    pred[i]  <- predict(mod, newx=t(as.matrix(xtest)), s=0.001)
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}




### CROSS VALIDATION WITH GROWING TRAINING WINDOW ###

# cv1.ols(x, y, k, z)
#
# Simple k-fold cross validation for a linear regression of time series with growing training window.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.ols <- function(x, y, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    xtrain <- x[1:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- lm(ytrain ~ xtrain)
    
    pred[i]  <- xtest %*% coefficients(mod)[-1] + coefficients(mod)[1]
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.logit(x, y, k, z)
#
# Simple k-fold cross validation for a logistic regression of time series with growing training window.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.logit <- function(x, y, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    xtrain <- x[1:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- glm(ytrain ~ xtrain, family = binomial(link = "logit"))
    
    pred[i]  <- 1 / (1 + exp(xtest %*% coefficients(mod)[-1] + coefficients(mod)[1]) )
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.pls(x, y, k)
#
# Simple k-fold cross validation for a boosted Partial Least Squares model of time series with growing training window.
#
# Returns a list with k OOS predictions, the associated errors (MAE and MSE), AIC (Greene approach) and 
# the number of PLS directions that min(MAE).
#

cv1.pls <- function(x, y, k=4){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  q <- ncol(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  mae_aux <- 1000000000
  
  for (d in 1:5){
    
    for (i in 1:k){
      
      ytrain <- y[1:(m+i-1),1]
      xtrain <- x[1:(m+i-1),]
      xtest  <- x[m+i,]
      
      mod <- pls(xtrain, ytrain, K=d)
      
      pred[i]  <- predict(mod, newdata=xtest)
      
      error[i] <- y[m+i] - pred[i]
      
    }
    
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k
    
    aic <- log(sum(error*error)/k) + 2*q/k
    
    if (mae <= mae_aux){
      
      pred_aux <- pred
      mae_aux  <- mae
      mse_aux  <- mse
      aic_aux  <- aic
      d_aux    <- d
    }
  }
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred_aux, "MAE" = mae_aux, "MSE" = mse_aux, "AIC" = aic_aux, "directions" = d_aux, "cum_error" = cumerr)
  
}


# cv1.pls2(x, y, k, d, z)
#
# Simple k-fold cross validation for a boosted Partial Least Squares model of time series with d directions and growing training window.
#
# Returns a list with k OOS predictions, the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.pls2 <- function(x, y, k=4, d=3, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  q <- ncol(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)

  for (i in 1:k){
  	      
    ytrain <- y[1:(m+i-1),1]
    xtrain <- x[1:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- pls(xtrain, ytrain, K=d)
    
    pred[i]  <- predict(mod, newdata=xtest)
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k

  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.lasso(x, y, k, z)
#
# Simple k-fold cross validation for a linear LASSO regression of time series with growing training window.
#
# Returns a list with k OOS predictions that min(AICc) and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.lasso <- function(x, y, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    xtrain <- x[1:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- gamlr(xtrain, ytrain, gamma=0)
    
    pred[i]  <- predict(mod, newdata=t(as.matrix(xtest)))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.lasso.logit(x, y, k, z)
#
# Simple k-fold cross validation for a logistic LASSO regression of time series with growing training window.
#
# Returns a list with k OOS predictions that min(AICc) and the associated errors (MAE and MSE) and AIC (Greene approach)..
#

cv1.lasso.logit <- function(x, y, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    xtrain <- x[1:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- gamlr(xtrain, ytrain, family = "binomial", gamma=0)
    
    pred[i]  <- predict(mod, newdata=t(as.matrix(xtest)))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.svm(x, y, k, z)
#
# Simple k-fold cross validation for a support vector machine (with polynomial kernel) of time series with growing training window.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.svm <- function(x, y, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    xtrain <- x[1:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- ksvm(ytrain ~ xtrain, kernel ="polydot")
    
    pred[i]  <- predict(mod, newdata=t(as.matrix(xtest)))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.rtree(formula, df, md, k, z)
#
# Simple k-fold cross validation for a regression tree of time series with growing training window. 
#
# The input must be a formula y ~ x, a data frame (df) and the min deviance improvement (md).
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.rtree <- function(formula, df, md=0.01, k=4, z=2){
  
  n <- nrow(df)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    dtrain <- df[1:(m+i-1),]
    dtest  <- df[m+i,]
    
    mod <- tree(formula, data = dtrain, mindev = md)
    
    pred[i]  <- predict(mod, newdata = dtest)
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(df) - 1
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.forest(formula, df, md, b, k, z)
#
# Simple k-fold cross validation for a forest of b regression trees of time series with growing training window.
#
# The input must be a formula y ~ x, a data frame (df), the min deviance improvement (md) and the number of trees (b) in the forest.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.forest <- function(formula, df, md=0.01, b=5, k=4, z=2){
  
  n <- nrow(df)
  m <- n - k
  
  pred  <- rep(0,k)
  pred2 <- matrix(0,k,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    dtrain <- df[1:(m+i-1),]
    dtest  <- df[(m+1):n,]
    
    mod <- tree(formula, data = dtrain, mindev = md)
    
    pred2[i,]  <- predict(mod, newdata = dtest)
    
  }
  
  for (i in 1:k){
    
    a <- max(1,i-b+1)
    
    pred[i] <- mean(pred2[a:i,i])
    
  }
  
  error <- y[(m+1):n] - pred
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(df) - 1
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.bag.ols(x, y, b, k, z)
#
# Simple k-fold cross validation for a bag of b linear regressions of time series with growing training window.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.bag.ols <- function(x, y, b=5, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  pred  <- rep(0,k)
  pred2 <- matrix(0,k,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    xtrain <- x[1:(m+i-1),]
    xtest  <- x[(m+1):n,]
    
    mod <- lm(ytrain ~ xtrain)
    
    pred2[i,]  <- t(xtest %*% coefficients(mod)[-1] + coefficients(mod)[1])
    
  }
  
  for (i in 1:k){
    
    a <- max(1,i-b+1)
    
    pred[i] <- mean(pred2[a:i,i])
    
  }
  
  error <- y[(m+1):n] - pred
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.drf(y, df, k, z)
#
# Simple k-fold cross validation for a Distributed Random Forest (DRF) with growing training window. 
#
# The input must be a "string" with the name of the response, a data frame (df) and the number of folds (k).
#
# NB: all columns of df except y are used as covariates. 
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.drf <- function(response, df, k=4, z=2){
  
  n <- nrow(df)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    dtrain <- as.h2o(df[1:(m+i-1),])
    dtest  <- as.h2o(df[m+i,])
    
    mod <- h2o.randomForest(y = response, training_frame = dtrain)
    
    pred[i]  <- as.vector(h2o.predict(mod, newdata = dtest))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(df) - 1
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.gbm(y, df, k, z)
#
# Simple k-fold cross validation for a Gradient Boosting Machine (GBM) of regression trees with growing training window. 
#
# The input must be a "string" with the name of the response, a data frame (df) and the number of folds (k).
#
# NB: all columns of df except y are used as covariates. 
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.gbm <- function(response, df, k=4, z=2){
  
  n <- nrow(df)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    dtrain <- as.h2o(df[1:(m+i-1),])
    dtest  <- as.h2o(df[m+i,])
    
    mod <- h2o.gbm(y = response, training_frame = dtrain)
    
    pred[i]  <- as.vector(h2o.predict(mod, newdata = dtest))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(df) - 1
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.rulefit(y, df, k, z)
#
# Simple k-fold cross validation for a RuleFit regression of time series with growing training window. 
#
# The general algorithm fits a tree ensemble to the data, builds a rule ensemble by traversing each tree, 
# evaluates the rules on the data to build a rule feature set, and fits a LASSO model to the rule feature 
# set joined with the original feature set.
#
# The input must be a "string" with the name of the response, a data frame (df) and the number of folds (k).
#
# NB: all columns of df except y are used as covariates. 
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.rulefit <- function(response, df, k=4, z=2){
  
  n <- nrow(df)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    dtrain <- as.h2o(df[1:(m+i-1),])
    dtest  <- as.h2o(df[m+i,])
    
    mod <- h2o.rulefit(y = response, training_frame = dtrain)
    
    pred[i]  <- as.vector(h2o.predict(mod, newdata = dtest))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(df) - 1
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.dnn(y, df, k, z)
#
# Simple k-fold cross validation for a Deep Neural Network (multi-layer perceptron) of time series with growing training window. 
#
# The input must be a "string" with the name of the response, a data frame (df) and the number of folds (k).
#
# NB: all columns of df except y are used as covariates. 
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.dnn <- function(response, df, k=4, z=2){
  
  n <- nrow(df)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    dtrain <- as.h2o(df[1:(m+i-1),])
    dtest  <- as.h2o(df[m+i,])
    
    mod <- h2o.deeplearning(y = response, training_frame = dtrain)
    
    pred[i]  <- as.vector(h2o.predict(mod, newdata = dtest))
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(df) - 1
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.ridge(x, y, k, z)
#
# Simple k-fold cross validation for a ridge regression of time series and growwing training window.
#
# Returns a list with k OOS predictions for lambda/s=0.001 and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.ridge <- function(x, y, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    xtrain <- x[1:(m+i-1),]
    xtest  <- x[m+i,]
    
    mod <- glmnet(xtrain, ytrain, alpha=0)
    
    pred[i]  <- predict(mod, newx=t(as.matrix(xtest)), s=0.001)
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x)
  
  aic <- log(sum(error*error)/k) + 2*q/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}




