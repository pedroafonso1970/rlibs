#
# libpls.R - A small library of Partial Least Squares methods in R, v01.03 (2025-05-15)
#
# Pedro Afonso Fernandes, UCP, CLSBE, Lisbon, Portugal (paf@ucp.pt)
# 
# This library is free software; you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation. The library is
# distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details: https://www.gnu.org/licenses/.
#


# Load libraries

library(robustbase)


# Set random number generator for reproducibility of covMcd() function from robustbase

set.seed(1970)


### UTILITIES ###

# zsm(X)
#
# Center and normalize a matrix of data by subtracting the mean of each column and dividing by its standard deviation (z-scores)
#

zsm <- function(X){
  
  m   <- apply(X,2,mean)
  s   <- apply(X,2,sd)
  
  R <- t((t(X) - m) / s) 
  
  return(R)
}


# zsv(X)
#
# Center and normalize a vector of data by subtracting its mean and dividing by its standard deviation (z-scores)
#

zsv <- function(x){
  
  m   <- mean(x)
  s   <- sd(x)
  
  y <- (x - m) / s 
  
  return(y)
}


# rzsm(X)
#
# Center and normalize a matrix of data by subtracting the mean of each column and dividing by its ROBUST standard deviation (z-scores)
#
# Uses the fast covMcd() function from the package robustbase.
#

rzsm <- function(X){
  
  mcd <- covMcd(X, seed = .Random.seed)
  
  m <- mcd$center
  
  robust_cov <- mcd$cov
  
  s <- sqrt(diag(robust_cov))
  
  R <- t((t(X) - m) / s) 
  
  return(R)
}


# rzsv(X)
#
# Center and normalize a vector of data by subtracting its mean and dividing by its ROBUST standard deviation (z-scores)
#
# Uses the fast covMcd() function from the package robustbase.
#

rzsv <- function(x){
  
  mcd <- covMcd(x, seed = .Random.seed)
  
  m <- as.vector(mcd$center)
  
  robust_cov <- as.vector(mcd$cov)
  
  s   <- sqrt(robust_cov)
  
  y <- (x - m) / s 
  
  return(y)
}


# norm(X)
#
# Normalize a matrix of data by dividing each column by its standard deviation
#

norm <- function(X){
  
  s   <- apply(X,2,sd)
  
  R <- t(t(X) / s) 
  
  return(R)
}


# norm2(X)
#
# Normalize a matrix of data by dividing each column by its square root of the sum of squares (Abdi, 2003).
#

norm2 <- function(X){
  
  p <- ncol(X)
  
  s   <- as.vector(matrix(0,1,p))
  
  for (i in 1:p){
    x    <- X[,i]
    s[i] <- sqrt(sum(x*x))
  }
  
  R <- t(t(X) / s) 
  
  return(R)
}


# rnorm(X)
#
# Normalize a matrix of data by dividing each column by its ROBUST standard deviation
#
# Uses the fast covMcd() function from the package robustbase.
#

rnorm <- function(X){
  
  mcd <- covMcd(X, seed = .Random.seed)
  robust_cov <- mcd$cov
  
  s <- sqrt(diag(robust_cov))
  
  R <- t(t(X) / s) 
  
  return(R)
}



# norv(X)
#
# Normalize a vector of data by dividing it by its standard deviation
#

norv <- function(x){
  
  s   <- sd(x)
  
  y <- x / s 
  
  return(y)
}


# norv2(X)
#
# Normalize a vector of data by dividing it by the square root of the sum of squares (Abdi, 2003).
#

norv2 <- function(x){
  
  s   <- sqrt(sum(x*x))
  
  y <- x / s 
  
  return(y)
}


# rnorv(X)
#
# Normalize a vector of data by dividing it by its ROBUST standard deviation
#
# Uses the fast covMcd() function from the package robustbase.
#

rnorv <- function(x){
  
  mcd <- covMcd(x, seed = .Random.seed)
  robust_cov <- mcd$cov
  
  s   <- sqrt(robust_cov)
  
  y <- x / s 
  
  return(y)
}



# rmean(X)
#
# Computes the ROBUST mean of a matrix or vector X.
#
# Uses the fast covMcd() function from the package robustbase.
#

rmean <- function(X){
  
  mcd <- covMcd(X, seed = .Random.seed)
  rm  <- mcd$center
  
  return(rm)
}


# rmean2(X)
#
# Computes the ROBUST mean of a matrix or vector X.
#
# Uses the covOGK() function from the package robustbase.
#

rmean2 <- function(X){
  
  mcd <- covOGK(X, sigmamu = scaleTau2)
  rm  <- mcd$center
  
  return(rm)
}



# rstd(X)
#
# Computes the ROBUST standard deviation of a matrix or vector X.
#
# Uses the fast covMcd() function from the package robustbase.
#

rstd <- function(X){
  
  mcd <- covMcd(X, seed = .Random.seed)
  rs <- sqrt(diag(mcd$cov))
  
  return(rs)
}


# rstd2(X)
#
# Computes the ROBUST standard deviation of a matrix or vector X.
#
# Uses the covOGK() function from the package robustbase.
#

rstd2 <- function(X){
  
  mcd <- covOGK(X, sigmamu = scaleTau2)
  rs <- sqrt(diag(mcd$cov))
  
  return(rs)
}



# rcov(X)
#
# Computes the ROBUST covariance (q x q) matrix of the (n x q) matrix or vector X.
#
# Uses the fast covMcd() function from the package robustbase.
#

rcov <- function(X){
  
  mcd <- covMcd(X, seed = .Random.seed)
  robust_cov <- mcd$cov
  
  return(robust_cov)
}


# rcov2(X)
#
# Alternative computation of the ROBUST covariance (q x q) matrix of the (n x q) matrix or vector X.
#
# Uses the covOGK() function from the package robustbase.
#

rcov2 <- function(X){
  
  mcd <- covOGK(X, sigmamu = scaleTau2)
  robust_cov <- mcd$cov
  
  return(robust_cov)
}



# rcrosscov(X,y)
#
# Computes the ROBUST cross-covariance (vector) between the (n x q) matrix X and the (n x 1) vector y.
#
# Uses the fast covMcd() function from the package robustbase.
#

rcrosscov <- function(X,y){
  
  Z <- cbind(X,y)
  
  mcd <- covMcd(Z, seed = .Random.seed)
  robust_cov <- mcd$cov
  
  q <- ncol(X)
  q1 <- q+1

  R <- robust_cov[1:q,q1]
  
  return(R)
}


# rcrosscov2(X,y)
#
# Alternative computation of the ROBUST cross-covariance (vector) between the (n x q) matrix X and the (n x 1) vector y.
#
# Uses the cov.rob(method = "classical") function from the package MASS (product-moment method).
#

rcrosscov2 <- function(X,y){
  
  Z <- cbind(X,y)
  
  mcd <- covOGK(Z, sigmamu = scaleTau2)
  
  robust_cov <- mcd$cov
  
  q <- ncol(X)
  q1 <- q+1
  
  R <- robust_cov[1:q,q1]
  
  return(R)
}



## MULTIVARIATE PARTIAL LEAST SQUARES (PLS2) ##

# mpls(X,Y,d)
#
# Estimates a multivariate partial least squares (PLS2) regression following the strategy of Abdi (2003).
#
# X and Y are the matrices of independent and dependent variables, respectively, and d is the number of PLS directions.
#
# It returns the X weights (W), X factor scores/directions (T), Y weights (C), Y scores (U), loading matrix (P),
# forward coefficients to predict Y (b), estimated X (X_hat), estimated Y (Y_hat) and error's metrics.
#

mpls <- function(X, Y, d=3){
  
  # Precision for convergence and max of iterations
  
  epsilon <- 2.2204e-16
  maxit   <- 100
  
  # Dimensions
  
  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(X)
  
  # Initialization
  
  W <- matrix(0,q,d)  # X weights
  T <- matrix(0,n,d)  # X factor
  
  C <- matrix(0,p,d)  # Y weights
  U <- matrix(0,n,d)  # Y factor
  
  P <- matrix(0,q,d)  # X loadings
  
  b <- matrix(0,1,d)  # b coefficient
  
  # Z-scores
  
  E <- zsm(X)
  F <- zsm(Y)
  
  # Descriptive statistics
  
  m_X   <- apply(X,2,mean)
  s_X   <- apply(X,2,sd)
  
  m_Y   <- apply(Y,2,mean)
  s_Y   <- apply(Y,2,sd)
  
  # Main cycle
  
  for (i in 1:d){
    
    t <- norv2(F[,1])
    u <- t
    
    dif <- 9999
    j    <- 1
    
    while (dif > epsilon & j < maxit) {
      t0 <- t
      
      w <- norv2(t(E) %*% u)
      t <- norv2(E %*% w)
      c <- norv2(t(F) %*% t)
      u <- F %*% c
      
      dif <- t(t0-t) %*% (t0-t)
      
      j <- j+1
    }
    
    pi  <- t(E) %*% t
    bi <- t(u) %*% t
    
    W[,i] <- w
    T[,i] <- t
    C[,i] <- c
    U[,i] <- u
    P[,i] <- pi
    b[1,i]  <- bi
    
    E <- E - t %*% t(pi)
    F <- F - (as.vector(bi) * (t %*% t(c)))
    
  }
  
  # Prediction of independent and dependent variables
  
  if (d==1){
    B <- b
  }
  else {
    B <- diag(as.vector(b))  
  }
  
  
  E_hat <- T %*% t(P)
  
  if (q==1){
    X_hat <- E_hat * s_X + m_X
  }
  else {
    X_hat <- E_hat %*% diag(as.vector(s_X)) + matrix(1,n,1) %*% m_X
  }
  
  F_hat <- T %*% B %*% t(C)
  
  if (p==1){
    Y_hat <- F_hat * s_Y + m_Y
  }
  else {
    Y_hat <- F_hat %*% diag(as.vector(s_Y)) + matrix(1,n,1) %*% m_Y
  }
  

  # In-sample errors for the dependents variables
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/(n*p)
  mse <- sum(error*error)/(n*p)
  
  aic <- log(sum(error*error)/(n*p)) + 2*q/n
  
  cumerr <- cumsum(abs(apply(error,1,sum)))

  # Output
  
  list("W" = W, "T" = T, "C" = C, "U" = U, "P" = P, "b" = b, "X_hat" = X_hat, "Y_hat" = Y_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr, "m_X" = m_X, "s_X" = s_X, "m_Y" = m_Y, "s_Y" = s_Y)
  
}  


# zmpls(X,Y,Z,d)
#
# Variant of mpls() with out-of-sample estimate YZ_hat given Z.
#

zmpls <- function(X, Y, Z, d=3){
  
  # Precision for convergence and max of iterations
  
  epsilon <- 2.2204e-16
  maxit   <- 100
  
  # Dimensions
  
  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(X)
  
  nz <- nrow(Z) # OOS dimension
  
  # Initialization
  
  W <- matrix(0,q,d)  # X weights
  T <- matrix(0,n,d)  # X factor
  
  C <- matrix(0,p,d)  # Y weights
  U <- matrix(0,n,d)  # Y factor
  
  P <- matrix(0,q,d)  # X loadings
  
  b <- matrix(0,1,d)  # b coefficient
  
  TZ <- matrix(0,nz,d) # OOS factor
  PZ <- matrix(0,q,d)  # OOS loadings
  
  # Z-scores
  
  E <- zsm(X)
  F <- zsm(Y)
  
  # Descriptive statistics
  
  m_X   <- apply(X,2,mean)
  s_X   <- apply(X,2,sd)
  
  m_Y   <- apply(Y,2,mean)
  s_Y   <- apply(Y,2,sd)
  
  # OOS scores
  
  XZ <- rbind(X,Z)
  
  m_XZ   <- apply(XZ,2,mean)
  s_XZ   <- apply(XZ,2,sd)
  
  G <- t((t(Z) - m_XZ) / s_XZ)
  

  # Main cycle
  
  for (i in 1:d){
    
    t <- norv2(F[,1])
    u <- t
    
    dif <- 9999
    j    <- 1
    
    while (dif > epsilon & j < maxit) {
      t0 <- t
      
      w <- norv2(t(E) %*% u)
      t <- norv2(E %*% w)
      c <- norv2(t(F) %*% t)
      u <- F %*% c
      
      tz <- norv2(G %*% w)  # OOS factor
      
      dif <- t(t0-t) %*% (t0-t)
      
      j <- j+1
    }
    
    pi  <- t(E) %*% t
    bi <- t(u) %*% t
    
    W[,i] <- w
    T[,i] <- t
    C[,i] <- c
    U[,i] <- u
    P[,i] <- pi
    b[1,i]  <- bi
    
    E <- E - t %*% t(pi)
    F <- F - (as.vector(bi) * (t %*% t(c)))

    # OOS loadings and deflated score
    
    TZ[,i] <- tz
    piz  <- t(G) %*% tz
    G <- G - tz %*% t(piz)
    
  }
  
  # Prediction of independent and dependent variables
  
  if (d==1){
    B <- b
  }
  else {
    B <- diag(as.vector(b))  
  }
  
  
  E_hat <- T %*% t(P)
  
  if (q==1){
    X_hat <- E_hat * s_X + m_X
  }
  else {
    X_hat <- E_hat %*% diag(as.vector(s_X)) + matrix(1,n,1) %*% m_X
  }
  
  F_hat <- T %*% B %*% t(C)
  
  if (p==1){
    Y_hat <- F_hat * s_Y + m_Y
  }
  else {
    Y_hat <- F_hat %*% diag(as.vector(s_Y)) + matrix(1,n,1) %*% m_Y
  }
  
  # OOS prediction
  
  H_hat <- TZ %*% B %*% t(C)
  
  if (p==1){
    YZ_hat <- H_hat * s_Y + m_Y
  }
  else {
    YZ_hat <- H_hat %*% diag(as.vector(s_Y)) + matrix(1,nz,1) %*% m_Y
  }
  
  
  # In-sample errors for the dependents variables
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/(n*p)
  mse <- sum(error*error)/(n*p)
  
  aic <- log(sum(error*error)/(n*p)) + 2*q/n
  
  cumerr <- cumsum(abs(apply(error,1,sum)))
  
  # Output
  
  list("W" = W, "T" = T, "C" = C, "U" = U, "P" = P, "b" = b, "X_hat" = X_hat, "Y_hat" = Y_hat, "YZ_hat" = YZ_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr, "m_X" = m_X, "s_X" = s_X, "m_Y" = m_Y, "s_Y" = s_Y)
  
}  


# rmpls(X,Y,d)
#
# Estimates a ROBUST multivariate partial least squares (PLS2) regression following the strategy of Abdi (2003).
#
# Uses the fast covMcd() function from the package robustbase.
#
# X and Y are the matrices of independent and dependent variables, respectively, and d is the number of PLS directions.
#
# It returns the X weights (W), X factor scores/directions (T), Y weights (C), Y scores (U), loading matrix (P),
# forward coefficients to predict Y (b), estimated X (X_hat), estimated Y (Y_hat) and error's metrics.
#

rmpls <- function(X, Y, d=3){
  
  # Precision for convergence and max of iterations
  
  epsilon <- 2.2204e-16
  maxit   <- 100
  
  # Dimensions
  
  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(X)
  
  # Initialization
  
  W <- matrix(0,q,d)  # X weights
  T <- matrix(0,n,d)  # X factor
  
  C <- matrix(0,p,d)  # Y weights
  U <- matrix(0,n,d)  # Y factor
  
  P <- matrix(0,q,d)  # X loadings
  
  b <- matrix(0,1,d)  # b coefficient
  
  # ROBUST Descriptive statistics
  
  m_X   <- rmean(X)
  s_X   <- rstd(X)

  m_Y   <- rmean(Y)
  s_Y   <- rstd(Y)
 
  # ROBUST Z-scores
  
  E <- t((t(X) - m_X) / s_X)  # rzsm(X)
  F <- t((t(Y) - m_Y) / s_Y)  # rzsm(Y)
  
  # Main cycle
  
  for (i in 1:d){
    
    t <- norv2(F[,1])
    u <- t
    
    dif <- 9999
    j    <- 1
    
    while (dif > epsilon & j < maxit) {
      t0 <- t
      
      w <- norv2(rcrosscov(E,u)*n)
      t <- norv2(E %*% w)
      c <- norv2(t(F) %*% t)
      u <- F %*% c
      
      dif <- t(t0-t) %*% (t0-t)
      
      j <- j+1
    }
    
    pi  <- t(E) %*% t
    bi <- t(u) %*% t
    
    W[,i] <- w
    T[,i] <- t
    C[,i] <- c
    U[,i] <- u
    P[,i] <- pi
    b[1,i]  <- bi
    
    E <- E - t %*% t(pi)
    F <- F - (as.vector(bi) * (t %*% t(c)))
    
  }
  
  # Prediction of independent and dependent variables
  
  if (d==1){
    B <- b
  }
  else {
    B <- diag(as.vector(b))  
  }
  
  
  E_hat <- T %*% t(P)
  
  if (q==1){
    X_hat <- E_hat * s_X + m_X
  }
  else {
    X_hat <- E_hat %*% diag(as.vector(s_X)) + matrix(1,n,1) %*% m_X
  }
  
  F_hat <- T %*% B %*% t(C)
  
  if (p==1){
    Y_hat <- F_hat * s_Y + m_Y
  }
  else {
    Y_hat <- F_hat %*% diag(as.vector(s_Y)) + matrix(1,n,1) %*% m_Y
  }
  
  
  # In-sample errors for the dependents variables
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/(n*p)
  mse <- sum(error*error)/(n*p)
  
  aic <- log(sum(error*error)/(n*p)) + 2*q/n
  
  cumerr <- cumsum(abs(apply(error,1,sum)))
  
  # Output
  
  list("W" = W, "T" = T, "C" = C, "U" = U, "P" = P, "b" = b, "X_hat" = X_hat, "Y_hat" = Y_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr, "m_X" = m_X, "s_X" = s_X, "m_Y" = m_Y, "s_Y" = s_Y)
  
}  


# rmpls2(X,Y,d)
#
# Estimates a ROBUST multivariate partial least squares (PLS2) regression following the strategy of Abdi (2003).
#
# Uses the covOGK() function from the package robustbase.
#
# X and Y are the matrices of independent and dependent variables, respectively, and d is the number of PLS directions.
#
# It returns the X weights (W), X factor scores/directions (T), Y weights (C), Y scores (U), loading matrix (P),
# forward coefficients to predict Y (b), estimated X (X_hat), estimated Y (Y_hat) and error's metrics.
#

rmpls2 <- function(X, Y, d=3){
  
  # Precision for convergence and max of iterations
  
  epsilon <- 2.2204e-16
  maxit   <- 100
  
  # Dimensions
  
  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(X)
  
  # Initialization
  
  W <- matrix(0,q,d)  # X weights
  T <- matrix(0,n,d)  # X factor
  
  C <- matrix(0,p,d)  # Y weights
  U <- matrix(0,n,d)  # Y factor
  
  P <- matrix(0,q,d)  # X loadings
  
  b <- matrix(0,1,d)  # b coefficient
  
  # ROBUST Descriptive statistics
  
  m_X   <- rmean2(X)
  s_X   <- rstd2(X)
  
  m_Y   <- rmean2(Y)
  s_Y   <- rstd2(Y)
  
  # ROBUST Z-scores
  
  E <- t((t(X) - m_X) / s_X)  # rzsm(X)
  F <- t((t(Y) - m_Y) / s_Y)  # rzsm(Y)
  
  # Main cycle
  
  for (i in 1:d){
    
    t <- norv2(F[,1])
    u <- t
    
    dif <- 9999
    j    <- 1
    
    while (dif > epsilon & j < maxit) {
      t0 <- t
      
      w <- norv2(rcrosscov2(E,u)*n)
      t <- norv2(E %*% w)
      c <- norv2(t(F) %*% t)
      u <- F %*% c
      
      dif <- t(t0-t) %*% (t0-t)
      
      j <- j+1
    }
    
    pi  <- t(E) %*% t
    bi <- t(u) %*% t
    
    W[,i] <- w
    T[,i] <- t
    C[,i] <- c
    U[,i] <- u
    P[,i] <- pi
    b[1,i]  <- bi
    
    E <- E - t %*% t(pi)
    F <- F - (as.vector(bi) * (t %*% t(c)))
    
  }
  
  # Prediction of independent and dependent variables
  
  if (d==1){
    B <- b
  }
  else {
    B <- diag(as.vector(b))  
  }
  
  
  E_hat <- T %*% t(P)
  
  if (q==1){
    X_hat <- E_hat * s_X + m_X
  }
  else {
    X_hat <- E_hat %*% diag(as.vector(s_X)) + matrix(1,n,1) %*% m_X
  }
  
  F_hat <- T %*% B %*% t(C)
  
  if (p==1){
    Y_hat <- F_hat * s_Y + m_Y
  }
  else {
    Y_hat <- F_hat %*% diag(as.vector(s_Y)) + matrix(1,n,1) %*% m_Y
  }
  
  
  # In-sample errors for the dependents variables
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/(n*p)
  mse <- sum(error*error)/(n*p)
  
  aic <- log(sum(error*error)/(n*p)) + 2*q/n
  
  cumerr <- cumsum(abs(apply(error,1,sum)))
  
  # Output
  
  list("W" = W, "T" = T, "C" = C, "U" = U, "P" = P, "b" = b, "X_hat" = X_hat, "Y_hat" = Y_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr, "m_X" = m_X, "s_X" = s_X, "m_Y" = m_Y, "s_Y" = s_Y)
  
}  


# zrmpls(X,Y,Z,d)
#
# Variant of rmpls() with out-of-sample estimate YZ_hat given Z.
#

zrmpls <- function(X, Y, Z, d=3){
  
  # Precision for convergence and max of iterations
  
  epsilon <- 2.2204e-16
  maxit   <- 100
  
  # Dimensions
  
  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(X)
  
  nz <- nrow(Z) # OOS dimension
  
  # Initialization
  
  W <- matrix(0,q,d)  # X weights
  T <- matrix(0,n,d)  # X factor
  
  C <- matrix(0,p,d)  # Y weights
  U <- matrix(0,n,d)  # Y factor
  
  P <- matrix(0,q,d)  # X loadings
  
  b <- matrix(0,1,d)  # b coefficient
  
  TZ <- matrix(0,nz,d) # OOS factor
  PZ <- matrix(0,q,d)  # OOS loadings
  
  # ROBUST Descriptive statistics
  
  m_X   <- rmean(X)
  s_X   <- rstd(X)
  
  m_Y   <- rmean(Y)
  s_Y   <- rstd(Y)
  
  # ROBUST Z-scores
  
  E <- t((t(X) - m_X) / s_X)  # rzsm(X)
  F <- t((t(Y) - m_Y) / s_Y)  # rzsm(Y)
  
  # ROBUST OOS scores
  
  XZ <- rbind(X,Z)
  
  m_XZ   <- rmean(XZ)
  s_XZ   <- rstd(XZ)
  
  G <- t((t(Z) - m_XZ) / s_XZ)
  
  # Main cycle
  
  for (i in 1:d){
    
    t <- norv2(F[,1])
    u <- t
    
    dif <- 9999
    j    <- 1
    
    while (dif > epsilon & j < maxit) {
      t0 <- t
      
      w <- norv2(rcrosscov(E,u)*n)
      t <- norv2(E %*% w)
      c <- norv2(t(F) %*% t)
      u <- F %*% c
      
      tz <- norv2(G %*% w)  # OOS factor
      
      dif <- t(t0-t) %*% (t0-t)
      
      j <- j+1
    }
    
    pi  <- t(E) %*% t
    bi <- t(u) %*% t
    
    W[,i] <- w
    T[,i] <- t
    C[,i] <- c
    U[,i] <- u
    P[,i] <- pi
    b[1,i]  <- bi
    
    E <- E - t %*% t(pi)
    F <- F - (as.vector(bi) * (t %*% t(c)))
    
    # OOS loadings and deflated score
    
    TZ[,i] <- tz
    piz  <- t(G) %*% tz
    G <- G - tz %*% t(piz)
    
  }
  
  # Prediction of independent and dependent variables
  
  if (d==1){
    B <- b
  }
  else {
    B <- diag(as.vector(b))  
  }
  
  
  E_hat <- T %*% t(P)
  
  if (q==1){
    X_hat <- E_hat * s_X + m_X
  }
  else {
    X_hat <- E_hat %*% diag(as.vector(s_X)) + matrix(1,n,1) %*% m_X
  }
  
  F_hat <- T %*% B %*% t(C)
  
  if (p==1){
    Y_hat <- F_hat * s_Y + m_Y
  }
  else {
    Y_hat <- F_hat %*% diag(as.vector(s_Y)) + matrix(1,n,1) %*% m_Y
  }
  
  # OOS prediction: multivariate
  
  H_hat <- TZ %*% B %*% t(C)
  
  if (p==1){
    YZ_hat <- H_hat * s_Y + m_Y
  }
  else {
    YZ_hat <- H_hat %*% diag(as.vector(s_Y)) + matrix(1,nz,1) %*% m_Y
  }
  
  # In-sample errors for the dependents variables
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/(n*p)
  mse <- sum(error*error)/(n*p)
  
  aic <- log(sum(error*error)/(n*p)) + 2*q/n
  
  cumerr <- cumsum(abs(apply(error,1,sum)))
  
  # Output
  
  list("W" = W, "T" = T, "C" = C, "U" = U, "P" = P, "b" = b, "X_hat" = X_hat, "Y_hat" = Y_hat, "YZ_hat" = YZ_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr, "m_X" = m_X, "s_X" = s_X, "m_Y" = m_Y, "s_Y" = s_Y)
  
}  


# zrmpls2(X,Y,Z,d)
#
# Variant of rmpls2() with out-of-sample estimate YZ_hat given Z.
#

zrmpls2 <- function(X, Y, Z, d=3){
  
  # Precision for convergence and max of iterations
  
  epsilon <- 2.2204e-16
  maxit   <- 100
  
  # Dimensions
  
  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(X)
  
  nz <- nrow(Z) # OOS dimension
  
  # Initialization
  
  W <- matrix(0,q,d)  # X weights
  T <- matrix(0,n,d)  # X factor
  
  C <- matrix(0,p,d)  # Y weights
  U <- matrix(0,n,d)  # Y factor
  
  P <- matrix(0,q,d)  # X loadings
  
  b <- matrix(0,1,d)  # b coefficient
  
  TZ <- matrix(0,nz,d) # OOS factor
  PZ <- matrix(0,q,d)  # OOS loadings
  
  # ROBUST Descriptive statistics
  
  m_X   <- rmean2(X)
  s_X   <- rstd2(X)
  
  m_Y   <- rmean2(Y)
  s_Y   <- rstd2(Y)
  
  # ROBUST Z-scores
  
  E <- t((t(X) - m_X) / s_X)  # rzsm(X)
  F <- t((t(Y) - m_Y) / s_Y)  # rzsm(Y)
  
  # ROBUST OOS scores
  
  XZ <- rbind(X,Z)
  
  m_XZ   <- rmean2(XZ)
  s_XZ   <- rstd2(XZ)
  
  G <- t((t(Z) - m_XZ) / s_XZ)
  
  # Main cycle
  
  for (i in 1:d){
    
    t <- norv2(F[,1])
    u <- t
    
    dif <- 9999
    j    <- 1
    
    while (dif > epsilon & j < maxit) {
      t0 <- t
      
      w <- norv2(rcrosscov2(E,u)*n)
      t <- norv2(E %*% w)
      c <- norv2(t(F) %*% t)
      u <- F %*% c
      
      tz <- norv2(G %*% w)  # OOS factor
      
      dif <- t(t0-t) %*% (t0-t)
      
      j <- j+1
    }
    
    pi  <- t(E) %*% t
    bi <- t(u) %*% t
    
    W[,i] <- w
    T[,i] <- t
    C[,i] <- c
    U[,i] <- u
    P[,i] <- pi
    b[1,i]  <- bi
    
    E <- E - t %*% t(pi)
    F <- F - (as.vector(bi) * (t %*% t(c)))
    
    # OOS loadings and deflated score
    
    TZ[,i] <- tz
    piz  <- t(G) %*% tz
    G <- G - tz %*% t(piz)
    
  }
  
  # Prediction of independent and dependent variables
  
  if (d==1){
    B <- b
  }
  else {
    B <- diag(as.vector(b))  
  }
  
  
  E_hat <- T %*% t(P)
  
  if (q==1){
    X_hat <- E_hat * s_X + m_X
  }
  else {
    X_hat <- E_hat %*% diag(as.vector(s_X)) + matrix(1,n,1) %*% m_X
  }
  
  F_hat <- T %*% B %*% t(C)
  
  if (p==1){
    Y_hat <- F_hat * s_Y + m_Y
  }
  else {
    Y_hat <- F_hat %*% diag(as.vector(s_Y)) + matrix(1,n,1) %*% m_Y
  }
  
  # OOS prediction: multivariate
  
  H_hat <- TZ %*% B %*% t(C)
  
  if (p==1){
    YZ_hat <- H_hat * s_Y + m_Y
  }
  else {
    YZ_hat <- H_hat %*% diag(as.vector(s_Y)) + matrix(1,nz,1) %*% m_Y
  }
  
  # In-sample errors for the dependents variables
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/(n*p)
  mse <- sum(error*error)/(n*p)
  
  aic <- log(sum(error*error)/(n*p)) + 2*q/n
  
  cumerr <- cumsum(abs(apply(error,1,sum)))
  
  # Output
  
  list("W" = W, "T" = T, "C" = C, "U" = U, "P" = P, "b" = b, "X_hat" = X_hat, "Y_hat" = Y_hat, "YZ_hat" = YZ_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr, "m_X" = m_X, "s_X" = s_X, "m_Y" = m_Y, "s_Y" = s_Y)
  
}  



## UNIVARIATE PARTIAL LEAST SQUARES (PLS1) ##

# pls(X,Y,d)
#
# Estimates a univariate partial least squares (PLS1) regression following the strategy of Abdi (2003).
#
# X is the matrix of independent variables, Y the vector of the dependent variable, and d is the number of PLS directions.
#
# It returns the X weights (W), X factor scores/directions (T), loading matrix (P),
# forward coefficients to predict Y (b), estimated X (X_hat), estimated Y (Y_hat) and error's metrics.
#

pls <- function(X, Y, d=3){
  
  # Dimensions
  
  n <- length(Y)
  q <- ncol(X)
  
  # Initialization
  
  W <- matrix(0,q,d)  # X weights
  T <- matrix(0,n,d)  # X factor
  
  P <- matrix(0,q,d)  # X loadings
  
  b <- matrix(0,1,d)  # b coefficient
  
  # Z-scores
  
  E <- zsm(X)
  F <- zsv(Y)
  
  # Descriptive statistics
  
  m_X   <- apply(X,2,mean)
  s_X   <- apply(X,2,sd)
  
  m_Y   <- mean(Y)
  s_Y   <- sd(Y)
  
  # Main cycle
  
  for (i in 1:d){
    
    w <- norv2(t(E) %*% F)
    t <- norv2(E %*% w)
    
    pi  <- t(E) %*% t
    bi <- t(F) %*% t
    
    W[,i] <- w
    T[,i] <- t
    P[,i] <- pi
    b[1,i]  <- bi
    
    E <- E - t %*% t(pi)
    F <- F - (as.vector(bi) * t)
    
  }
  
  # Prediction of independent and dependent variables
  
  if (d==1){
    B <- b
  }
  else {
    B <- diag(as.vector(b))  
  }
  
  
  E_hat <- T %*% t(P)
  
  if (q==1){
    X_hat <- E_hat * s_X + m_X
  }
  else {
    X_hat <- E_hat %*% diag(as.vector(s_X)) + matrix(1,n,1) %*% m_X
  }
  
  F_hat <- T %*% B %*% matrix(1,d,1)
  
  Y_hat <- F_hat * s_Y + m_Y
  
  
  # In-sample errors for the dependents variables
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/n
  mse <- sum(error*error)/n
  
  aic <- log(sum(error*error)/n) + 2*q/n
  
  cumerr <- cumsum(abs(error))
  
  # Output
  
  list("W" = W, "T" = T, "P" = P, "b" = b, "X_hat" = X_hat, "Y_hat" = Y_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr, "m_X" = m_X, "s_X" = s_X, "m_Y" = m_Y, "s_Y" = s_Y)
  
}  


# zpls(X,Y,Z,d)
#
# Variant of pls() with out-of-sample estimate YZ_hat given Z.
#

zpls <- function(X, Y, Z, d=3){
  
  # Dimensions
  
  n <- length(Y)
  q <- ncol(X)
  
  nz <- nrow(Z) # OOS dimension
  
  # Initialization
  
  W <- matrix(0,q,d)  # X weights
  T <- matrix(0,n,d)  # X factor
  
  P <- matrix(0,q,d)  # X loadings
  
  b <- matrix(0,1,d)  # b coefficient
  
  TZ <- matrix(0,nz,d) # OOS factor
  PZ <- matrix(0,q,d)  # OOS loadings
  
  # Z-scores
  
  E <- zsm(X)
  F <- zsv(Y)
  
  # Descriptive statistics
  
  m_X   <- apply(X,2,mean)
  s_X   <- apply(X,2,sd)
  
  m_Y   <- mean(Y)
  s_Y   <- sd(Y)
  
  # OOS scores
  
  XZ <- rbind(X,Z)
  
  m_XZ   <- apply(XZ,2,mean)
  s_XZ   <- apply(XZ,2,sd)
  
  G <- t((t(Z) - m_XZ) / s_XZ)
  
  # Main cycle
  
  for (i in 1:d){
    
    w <- norv2(t(E) %*% F)
    t <- norv2(E %*% w)
    
    tz <- norv2(G %*% w)  # OOS factor
    
    pi  <- t(E) %*% t
    bi <- t(F) %*% t
    
    W[,i] <- w
    T[,i] <- t
    P[,i] <- pi
    b[1,i]  <- bi
    
    E <- E - t %*% t(pi)
    F <- F - (as.vector(bi) * t)
    
    # OOS loadings and deflated score
    
    TZ[,i] <- tz
    piz  <- t(G) %*% tz
    G <- G - tz %*% t(piz)
    
  }
  
  # Prediction of independent and dependent variables
  
  if (d==1){
    B <- b
  }
  else {
    B <- diag(as.vector(b))  
  }
  
  
  E_hat <- T %*% t(P)
  
  if (q==1){
    X_hat <- E_hat * s_X + m_X
  }
  else {
    X_hat <- E_hat %*% diag(as.vector(s_X)) + matrix(1,n,1) %*% m_X
  }
  
  F_hat <- T %*% B %*% matrix(1,d,1)
  
  Y_hat <- F_hat * s_Y + m_Y

  # OOS prediction
  
  H_hat <- TZ %*% B %*% matrix(1,d,1)
  
  YZ_hat <- H_hat * s_Y + m_Y
  

  # In-sample errors for the dependents variables
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/n
  mse <- sum(error*error)/n
  
  aic <- log(sum(error*error)/n) + 2*q/n
  
  cumerr <- cumsum(abs(error))
  
  # Output
  
  list("W" = W, "T" = T, "P" = P, "b" = b, "X_hat" = X_hat, "Y_hat" = Y_hat, "YZ_hat" = YZ_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr, "m_X" = m_X, "s_X" = s_X, "m_Y" = m_Y, "s_Y" = s_Y)
  
}  


# rpls(X,Y,d)
#
# Estimates a ROBUST univariate partial least squares (PLS1) regression following the strategy of Abdi (2003).
#
# Uses the fast covMcd() function from the package robustbase.
# 
# X is the matrix of independent variables, Y the vector of the dependent variable, and d is the number of PLS directions.
#
# It returns the X weights (W), X factor scores/directions (T), loading matrix (P),
# forward coefficients to predict Y (b), estimated X (X_hat), estimated Y (Y_hat) and error's metrics.
#

rpls <- function(X, Y, d=3){
  
  # Dimensions
  
  n <- length(Y)
  q <- ncol(X)
  
  # Initialization
  
  W <- matrix(0,q,d)  # X weights
  T <- matrix(0,n,d)  # X factor
  
  P <- matrix(0,q,d)  # X loadings
  
  b <- matrix(0,1,d)  # b coefficient
  
  # ROBUST descriptive statistics
  
  m_X   <- rmean(X)
  s_X   <- rstd(X)
  
  m_Y   <- rmean(Y)
  s_Y   <- rstd(Y)
  
  # ROBUST Z-scores
  
  E <- t((t(X) - m_X) / s_X)  # rzsm(X)
  F <- (Y - m_Y) / s_Y        # rzsV(Y)
  
  # Main cycle
  
  for (i in 1:d){
    
    w <- norv2(rcrosscov(E,F)*n)
    t <- norv2(E %*% w)
    
    pi  <- t(E) %*% t
    bi <- t(F) %*% t
    
    W[,i] <- w
    T[,i] <- t
    P[,i] <- pi
    b[1,i]  <- bi
    
    E <- E - t %*% t(pi)
    F <- F - (as.vector(bi) * t)
    
  }
  
  # Prediction of independent and dependent variables
  
  if (d==1){
    B <- b
  }
  else {
    B <- diag(as.vector(b))  
  }
  
  
  E_hat <- T %*% t(P)
  
  if (q==1){
    X_hat <- E_hat * s_X + m_X
  }
  else {
    X_hat <- E_hat %*% diag(as.vector(s_X)) + matrix(1,n,1) %*% m_X
  }
  
  F_hat <- T %*% B %*% matrix(1,d,1)
  
  Y_hat <- F_hat * s_Y + m_Y
  
  
  # In-sample errors for the dependents variables
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/n
  mse <- sum(error*error)/n
  
  aic <- log(sum(error*error)/n) + 2*q/n
  
  cumerr <- cumsum(abs(error))
  
  # Output
  
  list("W" = W, "T" = T, "P" = P, "b" = b, "X_hat" = X_hat, "Y_hat" = Y_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr, "m_X" = m_X, "s_X" = s_X, "m_Y" = m_Y, "s_Y" = s_Y)
  
}  


# rpls2(X,Y,d)
#
# Estimates a ROBUST univariate partial least squares (PLS1) regression following the strategy of Abdi (2003).
#
# Uses the fast covOGK() function from the package robustbase.
# 
# X is the matrix of independent variables, Y the vector of the dependent variable, and d is the number of PLS directions.
#
# It returns the X weights (W), X factor scores/directions (T), loading matrix (P),
# forward coefficients to predict Y (b), estimated X (X_hat), estimated Y (Y_hat) and error's metrics.
#

rpls2 <- function(X, Y, d=3){
  
  # Dimensions
  
  n <- length(Y)
  q <- ncol(X)
  
  # Initialization
  
  W <- matrix(0,q,d)  # X weights
  T <- matrix(0,n,d)  # X factor
  
  P <- matrix(0,q,d)  # X loadings
  
  b <- matrix(0,1,d)  # b coefficient
  
  # ROBUST descriptive statistics
  
  m_X   <- rmean2(X)
  s_X   <- rstd2(X)
  
  m_Y   <- rmean(Y)
  s_Y   <- rstd(Y)
  
  # ROBUST Z-scores
  
  E <- t((t(X) - m_X) / s_X)  # rzsm(X)
  F <- (Y - m_Y) / s_Y        # rzsV(Y)
  
  # Main cycle
  
  for (i in 1:d){
    
    w <- norv2(rcrosscov2(E,F)*n)
    t <- norv2(E %*% w)
    
    pi  <- t(E) %*% t
    bi <- t(F) %*% t
    
    W[,i] <- w
    T[,i] <- t
    P[,i] <- pi
    b[1,i]  <- bi
    
    E <- E - t %*% t(pi)
    F <- F - (as.vector(bi) * t)
    
  }
  
  # Prediction of independent and dependent variables
  
  if (d==1){
    B <- b
  }
  else {
    B <- diag(as.vector(b))  
  }
  
  
  E_hat <- T %*% t(P)
  
  if (q==1){
    X_hat <- E_hat * s_X + m_X
  }
  else {
    X_hat <- E_hat %*% diag(as.vector(s_X)) + matrix(1,n,1) %*% m_X
  }
  
  F_hat <- T %*% B %*% matrix(1,d,1)
  
  Y_hat <- F_hat * s_Y + m_Y
  
  
  # In-sample errors for the dependents variables
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/n
  mse <- sum(error*error)/n
  
  aic <- log(sum(error*error)/n) + 2*q/n
  
  cumerr <- cumsum(abs(error))
  
  # Output
  
  list("W" = W, "T" = T, "P" = P, "b" = b, "X_hat" = X_hat, "Y_hat" = Y_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr, "m_X" = m_X, "s_X" = s_X, "m_Y" = m_Y, "s_Y" = s_Y)
  
}  



# zrpls(X,Y,Z,d)
#
# Variant of rpls() with out-of-sample estimate YZ_hat given Z.
#

zrpls <- function(X, Y, Z, d=3){
  
  # Dimensions
  
  n <- length(Y)
  q <- ncol(X)
  
  nz <- nrow(Z) # OOS dimension
  
  # Initialization
  
  W <- matrix(0,q,d)  # X weights
  T <- matrix(0,n,d)  # X factor
  
  P <- matrix(0,q,d)  # X loadings
  
  b <- matrix(0,1,d)  # b coefficient
  
  TZ <- matrix(0,nz,d) # OOS factor
  PZ <- matrix(0,q,d)  # OOS loadings
  
  # ROBUST descriptive statistics
  
  m_X   <- rmean(X)
  s_X   <- rstd(X)
  
  m_Y   <- rmean(Y)
  s_Y   <- rstd(Y)
  
  # ROBUST Z-scores
  
  E <- t((t(X) - m_X) / s_X)  # rzsm(X)
  F <- (Y - m_Y) / s_Y        # rzsV(Y)
  
  # ROBUST OOS scores
  
  XZ <- rbind(X,Z)
  
  m_XZ   <- rmean(XZ)
  s_XZ   <- rstd(XZ)
  
  G <- t((t(Z) - m_XZ) / s_XZ)

  # Main cycle
  
  for (i in 1:d){
    
    w <- norv2(rcrosscov(E,F)*n)
    t <- norv2(E %*% w)
    
    tz <- norv2(G %*% w)  # OOS factor
    
    pi  <- t(E) %*% t
    bi <- t(F) %*% t
    
    W[,i] <- w
    T[,i] <- t
    P[,i] <- pi
    b[1,i]  <- bi
    
    E <- E - t %*% t(pi)
    F <- F - (as.vector(bi) * t)
    
    # OOS loadings and deflated score
    
    TZ[,i] <- tz
    piz  <- t(G) %*% tz
    G <- G - tz %*% t(piz)
    
  }
  
  # Prediction of independent and dependent variables
  
  if (d==1){
    B <- b
  }
  else {
    B <- diag(as.vector(b))  
  }
  
  
  E_hat <- T %*% t(P)
  
  if (q==1){
    X_hat <- E_hat * s_X + m_X
  }
  else {
    X_hat <- E_hat %*% diag(as.vector(s_X)) + matrix(1,n,1) %*% m_X
  }
  
  F_hat <- T %*% B %*% matrix(1,d,1)
  
  Y_hat <- F_hat * s_Y + m_Y
  
  # OOS prediction: univariate
  
  H_hat <- TZ %*% B %*% matrix(1,d,1)
  
  YZ_hat <- H_hat * s_Y + m_Y
  
  
  # In-sample errors for the dependents variables
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/n
  mse <- sum(error*error)/n
  
  aic <- log(sum(error*error)/n) + 2*q/n
  
  cumerr <- cumsum(abs(error))
  
  # Output
  
  list("W" = W, "T" = T, "P" = P, "b" = b, "X_hat" = X_hat, "Y_hat" = Y_hat, "YZ_hat" = YZ_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr, "m_X" = m_X, "s_X" = s_X, "m_Y" = m_Y, "s_Y" = s_Y)
  
}  


# zrpls2(X,Y,Z,d)
#
# Variant of rpls2() with out-of-sample estimate YZ_hat given Z.
#

zrpls2 <- function(X, Y, Z, d=3){
  
  # Dimensions
  
  n <- length(Y)
  q <- ncol(X)
  
  nz <- nrow(Z) # OOS dimension
  
  # Initialization
  
  W <- matrix(0,q,d)  # X weights
  T <- matrix(0,n,d)  # X factor
  
  P <- matrix(0,q,d)  # X loadings
  
  b <- matrix(0,1,d)  # b coefficient
  
  TZ <- matrix(0,nz,d) # OOS factor
  PZ <- matrix(0,q,d)  # OOS loadings
  
  # ROBUST descriptive statistics
  
  m_X   <- rmean2(X)
  s_X   <- rstd2(X)
  
  m_Y   <- rmean(Y)
  s_Y   <- rstd(Y)
  
  # ROBUST Z-scores
  
  E <- t((t(X) - m_X) / s_X)  # rzsm(X)
  F <- (Y - m_Y) / s_Y        # rzsV(Y)
  
  # ROBUST OOS scores
  
  XZ <- rbind(X,Z)
  
  m_XZ   <- rmean2(XZ)
  s_XZ   <- rstd2(XZ)
  
  G <- t((t(Z) - m_XZ) / s_XZ)
  
  # Main cycle
  
  for (i in 1:d){
    
    w <- norv2(rcrosscov2(E,F)*n)
    t <- norv2(E %*% w)
    
    tz <- norv2(G %*% w)  # OOS factor
    
    pi  <- t(E) %*% t
    bi <- t(F) %*% t
    
    W[,i] <- w
    T[,i] <- t
    P[,i] <- pi
    b[1,i]  <- bi
    
    E <- E - t %*% t(pi)
    F <- F - (as.vector(bi) * t)
    
    # OOS loadings and deflated score
    
    TZ[,i] <- tz
    piz  <- t(G) %*% tz
    G <- G - tz %*% t(piz)
    
  }
  
  # Prediction of independent and dependent variables
  
  if (d==1){
    B <- b
  }
  else {
    B <- diag(as.vector(b))  
  }
  
  
  E_hat <- T %*% t(P)
  
  if (q==1){
    X_hat <- E_hat * s_X + m_X
  }
  else {
    X_hat <- E_hat %*% diag(as.vector(s_X)) + matrix(1,n,1) %*% m_X
  }
  
  F_hat <- T %*% B %*% matrix(1,d,1)
  
  Y_hat <- F_hat * s_Y + m_Y
  
  # OOS prediction: univariate
  
  H_hat <- TZ %*% B %*% matrix(1,d,1)
  
  YZ_hat <- H_hat * s_Y + m_Y
  
  
  # In-sample errors for the dependents variables
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/n
  mse <- sum(error*error)/n
  
  aic <- log(sum(error*error)/n) + 2*q/n
  
  cumerr <- cumsum(abs(error))
  
  # Output
  
  list("W" = W, "T" = T, "P" = P, "b" = b, "X_hat" = X_hat, "Y_hat" = Y_hat, "YZ_hat" = YZ_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr, "m_X" = m_X, "s_X" = s_X, "m_Y" = m_Y, "s_Y" = s_Y)
  
}  



# npls(X,Y,d)
#
# Estimates a nonlinear partial least squares regression with a logistic relationship between the latent variables
# of Y (vector u) and X (vector t): u = G(t) = Asym/(1+exp((xmid-t)/scal))
#
# where Asym is the asymptote of the logistic curve, xmid is its inflexion point and scal is a scale parameter.
#
# X and Y are the matrices of independent and dependent variables, respectively, and d is the number of PLS directions.
#
# It returns the X weights (W), X factor scores/directions (T), Y weights (C), Y scores (U), loading matrix (P),
# NLS/logistic coefficients to predict Y (B), estimated X (X_hat), estimated Y (Y_hat) and error's metrics.


npls <- function(X, Y, d=3){
  
  # Precision for convergence and max of iterations
  
  epsilon <- 2.2204e-16
  maxit   <- 100
  
  # Dimensions
  
  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(X)
  
  # Initialization
  
  W <- matrix(0,q,d)  # X weights
  T <- matrix(0,n,d)  # X factor
  
  C <- matrix(0,p,d)  # Y weights
  U <- matrix(0,n,d)  # Y factor
  
  P <- matrix(0,q,d)  # X loadings
  
  B <- matrix(0,4,d)  # Logistic coefficients (Asym, xmid and scal) and shift factor (alpha)
  
  # Z-scores
  
  E <- zsm(X)
  F <- zsm(Y)
  
  # Descriptive statistics
  
  m_X   <- apply(X,2,mean)
  s_X   <- apply(X,2,sd)
  
  m_Y   <- apply(Y,2,mean)
  s_Y   <- apply(Y,2,sd)
  
  # Main cycle
  
  for (i in 1:d){
    
    t <- norv2(F[,1])
    u <- t
    
    dif <- 9999
    j    <- 1
    
    while (dif > epsilon & j < maxit) {
      t0 <- t
      
      w <- norv2(t(E) %*% u)
      t <- norv2(E %*% w)
      c <- norv2(t(F) %*% t)
      u <- F %*% c
      
      dif <- t(t0-t) %*% (t0-t)
      
      j <- j+1
    }
    
    pi  <- t(E) %*% t
    
    # Nonlinear model u = G(t)
    
    alpha <- min(u)  # Shift factor, to be subtracted to u before the logistic adjustment
    
    dt <- as.data.frame(cbind(t = t, u = u-alpha))
    
    mod <- nls(u ~ SSlogis(t, Asym, xmid, scal), data=dt)
    
    ai <- coefficients(mod)[1]
    bi <- coefficients(mod)[2]
    ci <- coefficients(mod)[3]
    di <- alpha
    
    W[,i] <- w
    T[,i] <- t
    C[,i] <- c
    U[,i] <- u
    P[,i] <- pi
    
    B[1,i]  <- ai
    B[2,i]  <- bi
    B[3,i]  <- ci
    B[4,i]  <- di
    
    E <- E - t %*% t(pi)
    F <- F - (predict(mod) + alpha) %*% t(c)
    
  }
  
  # Prediction of independent and dependent variables
  
  E_hat <- T %*% t(P)
  
  if (q==1){
    X_hat <- E_hat * s_X + m_X
  }
  else {
    X_hat <- E_hat %*% diag(as.vector(s_X)) + matrix(1,n,1) %*% m_X
  }
  
  
  TB <- matrix(0,n,d)
  
  for (i in 1:d){
    TB[,i] <- B[1,i]/(1+exp((B[2,i]-T[,i])/B[3,i])) + B[4,i]
  }
  
  F_hat <- TB %*% t(C)
  
  
  if (p==1){
    Y_hat <- F_hat * s_Y + m_Y
  }
  else {
    Y_hat <- F_hat %*% diag(as.vector(s_Y)) + matrix(1,n,1) %*% m_Y
  }
  
  
  # In-sample errors for the dependents variables
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/(n*p)
  mse <- sum(error*error)/(n*p)
  
  aic <- log(sum(error*error)/(n*p)) + 2*q/n
  
  cumerr <- cumsum(abs(error))
  
  # Output
  
  list("W" = W, "T" = T, "C" = C, "U" = U, "P" = P, "B" = B, "X_hat" = X_hat, "Y_hat" = Y_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr, "m_X" = m_X, "s_X" = s_X, "m_Y" = m_Y, "s_Y" = s_Y)
  
}  



## UNIVARIATE ORDINARY LEAST SQUARES (OLS) ##

# ols(X,Y)
#
# Estimates an ordinary least squares (OLS) regression where X is the matrix of independent variables and
# Y the vector of the dependent variable.
#
# It returns the estimated Y (Y_hat) and in-sample error's metrics.
#

ols <- function(X, Y){
  
  n <- length(Y)
  q <- ncol(X)
  
  mod <- lm(Y~X)
  
  Y_hat <- mod$fitted.values
  
  error <- Y - Y_hat
  
  mae <- sum(abs(error))/n
  mse <- sum(error*error)/n
  
  aic <- log(sum(error*error)/n) + 2*q/n
  
  cumerr <- cumsum(abs(error))
  
  # Output
  
  list("Y_hat" = Y_hat, "Y_error" = error, "MAE" = mae, "MSE" = mse, "AIC" = aic, "cum_error" = cumerr)
  
}  


# fols(X,Y,Z)
#
# Forecasts the dependent variable with OLS given the out-of-sample attribute Z. 
#
# X is the matrix of independent variables, Y the vector of the dependent variable.
#

fols <- function(X, Y, Z){
  
  mod <- lm(Y~X)
  
  Y_hat <- Z %*% coefficients(mod)[-1] + coefficients(mod)[1]
  
  return(Y_hat)
}  





