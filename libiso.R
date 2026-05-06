#
# libiso.R - A small library of R functions for isotonic regression, v01.01 (2022-05-06)
#
# Pedro Afonso Fernandes, UCP, CLSBE, Lisbon, Portugal (paf@ucp.pt)
# 
# This library is free software; you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation. The library is
# distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details: https://www.gnu.org/licenses/.
#

# isofit(x,z)
#
# Monotone/isotone regression of one independent variable x on z (Gebhardt, 1970).
#

isofit <- function(x,z){
 
  n <- length(x)      # Number of observations
  
  if (n != length(z)) return("ERROR: x and z have different lengths!")
   
  dt <- cbind(x,z)    # Matrix with x and z binden in column
  
  o <- order(dt[,1])  # Index with the order of x
  
  xz <- dt[o,]        # Data matrix sorted by x
  
  xo <- xz[,1]        # x sorted 
  zo <- xz[,2]        # z sorted by x 
  
  u <- zo             # Initializes the regression values u = z sorted (default)
  m <- 1              # Initializes the size of the current (first) block
  s <- zo[1]          # Value of the current (first) block
  
  d     <- min(diff(u))  # Min slope of the sorted z series (should be positive)
  i     <- 1             # Number of iterations of the main algorithm
  maxit <- 1000          # Max number of iterations
  
  while (d < 0) {
  
    for (k in 2:n){
      if(u[k] < u[k-1]){
        # Non-monotone value of z, thus ...
        m <- m+1             # Increments the size of the current block
        s <- s + (u[k]-s)/m  # Updates the value (average) of the current block
        
        for (j in 1:m){
          u[k-j+1] <- s      # Update the regression values of the current block
        }
      }
      else{
        m <- 1               # Initializes the size of the current (k) block
        s <- u[k]            # Value of the current (k) block
      }
    }
  
    d <- min(diff(u))  # Updates the min slope of the fitted values (should be positive)
    
    i <- i+1
    if (i > maxit) break
      
  }
  
  aux <- list("yf" = u, "ord" = o)  # Follows the isoreg() notation from "stats"
  
  return(aux)
  
}


# concord(X,w)
#
# Concordance matrix of a decision matrix X with n actions (rows) and 
# m criteria (columns), given the weights vector w.
#
# Uses true-criteria like ELECTRE I (Roy, 1968).
#

concord <- function(X,w){

  n <- nrow(X)  # Number of actions
  m <- ncol(X)  # Number of criteria
  
  if (m != length(w)) return("ERROR: length(w) w is different from ncol(X)!")
  
  w <- w / sum(w)  # Normalizes the weights in the interval [0,1]
  
  C <- matrix(0,n,n)  # Initializes the concordance matrix
  
  for (i in 1:n){
    for (j in 1:n){
      for (k in 1:m){
        if(X[i,k] >= X[j,k]){
          C[i,j] <- C[i,j] + w[k]  # Updates the weights' sum for the pair of actions (i,j)
        }
      }
    }
  }
  
  diag(C) <- 0
  
  return(C)
  
}


# electre(X,w)
#
# Lightweight version of the ELECTRE III method with true-criteria like ELECTRE I,
# concordance analysis and simple ranking using the non-domination degree proposed by 
# Siskos and Hubert (1983).
#
# The inputs are a decision matrix X with n actions (rows) and m criteria 
# (columns), and the weights vector w.
#
# It returns the non-domination degree. The higher it is, the more preferable are
# the alternatives; thus, it can be used to order() the last ones.

electre <- function(X,w){
  
  n <- nrow(X)  # Number of actions
  
  C <- concord(X,w)  # Concordance matrix
  
  D <- matrix(0,n,n)  # Initializes the domination matrix
  
  for (i in 1:n){
    for (j in 1:n){
      D[i,j] <- C[j,i] - C[i,j]
    }
  }
  
  nd <- 1 - apply(D,2,max)  # Degree of non-domination
  
  return(nd)
  
}


# mfit(X,z)
#
# Simple multivariate Monotone/isotone regression of the columns of X on z.
#
# The columns of X are ordered with the method ELECTRE and then a univariate
# isotone regression is computed using that order as independent variable.
#

mfit <- function(X,z){
  
  n <- nrow(X)  # Number of observations ("actions")
  m <- ncol(X)  # Number of independent variables ("criteria")
  
  if (n != length(z)) return("ERROR: X and z have different lengths!")
 
  w <- rep(1,m)       # Vector of "criteria" weights
  
  nd <- electre(X,w)  # Degree of non-dominance of the observations ("actions")
  
  fit <- isofit(nd,z)  # Isotonic regression with nd as independent variable
  
  aux <- list("yf" = fit$yf, "ord" = fit$ord)  # Follows the isoreg() notation from "stats"
  
  return(aux)
  
}


# biv(x,y,z,a,b)
#
# Creates a bivariate grid matrix with the mean of z ordered on variables x and y.
#
# The ordered classes are defined as sequences between the min and max of x and y
# using a and b, respectively, as increments or steps.
#
# The output can be used as input for the function biviso() from package "Iso".
#

biv <- function(x,y,z,a=1,b=1){
  
  n <- length(x)      # Number of observations
  
  if (n != length(y)) return("ERROR: x and y have different lengths!")
  if (n != length(z)) return("ERROR: x and z have different lengths!")
  
  df <- as.data.frame(cbind(x,y,z))  # Data frame with all variables
  
  x_min = round(min(x))
  x_max = round(max(x))
  
  y_min = round(min(y))
  y_max = round(max(y))
  
  df$x_class <- cut(df$x, breaks = seq(x_min,x_max,a), include.lowest = TRUE)
  df$y_class <- cut(df$y, breaks = seq(y_min,y_max,b), include.lowest = TRUE)
  
  M <- tapply(df$z, list(df$x_class, df$y_class), FUN = mean, default = 0)
  
  return(M)

}

