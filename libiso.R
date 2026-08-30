#
# libiso.R - A small library of R functions for isotonic regression, v01.04 (2026-07-24)
#
# Pedro Afonso Fernandes, UCP, CLSBE, Lisbon, Portugal (paf@ucp.pt)
# 
# This library is free software; you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation. The library is
# distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details: https://www.gnu.org/licenses/.
#


# isofit(x,y)
#
# Monotone/isotonic regression of one independent variable x on the 
# the dependent variable y (Gebhardt, 1970).
#

isofit <- function(x,y){
 
  n <- length(x)      # Number of observations
  
  if (n != length(y)) return("ERROR: x and y have different lengths!")
   
  dt <- cbind(x,y)    # Matrix with x and y in column
  
  o <- order(dt[,1])  # Index with the order of x
  
  xy <- dt[o,]        # Data matrix sorted by x
  
  u <- xy[,2]         # Initializes the regression values u = y sorted (default)
  m <- 1              # Initializes the size of the current (first) block
  s <- u[1]           # Value of the current (first) block
  
  d     <- min(diff(u))  # Min slope of the sorted z series (should be positive)
  i     <- 1             # Number of iterations of the main algorithm
  maxit <- 1000          # Max number of iterations
  
  while (d < 0) {
  
    for (k in 2:n){
      if(u[k] < u[k-1]){
        # Non-monotone value of y, thus ...
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


# isodp(x,y)
#
# Unweighted isotonic L1 regression of one independent variable x on the 
# dependent variable y by dynamic programming (Rote, 2019).
#

isodp <- function(x,y){
  
  n <- length(x)      # Number of observations
  
  if (n != length(y)) return("ERROR: x and y have different lengths!")
  
  dt <- cbind(x,y)    # Matrix with x and y in column
  
  o <- order(dt[,1])  # Index with the order of x
  
  xy <- dt[o,]        # Data matrix sorted by x
  
  yo <- xy[,2]        # y sorted by x 
  
  p <- rep(0,n)       # arg min of each cost-to-go function
  z <- rep(0,n)       # Optimal solution
  
  n1 <- n-1
  
  Q <- matrix(-999999999,n+1,2) # Priority queue of breakpoints (position and value in column)
  
  for (k in 1:n){
    
    # Insert a new breakpoint with position yo[k] and value 2
    Q[k,1] <- yo[k]
    Q[k,2] <- 2
    
    # Sort the queue by position/priority (descending order)  
    Q <- Q[order(-Q[,1]),]
    
    # find/peek max i.e. the breakpoint with highest position/priority (Q head)
    B <- Q[1,2]
    
    if (B==1){
      # The max is from a previous interaction, so "delete" it (pop)
      Q[1,1] <- -999999999
      Q[1,2] <- -999999999
      
      # Sort the queue by position/priority (descending order) once again
      Q <- Q[order(-Q[,1]),]
    }
    else{
      # update the value with max priority (Q head)
      Q[1,2] <- 1
    }
    
    p[k] <- Q[1,1]  # Save the highest position (arg min of k cost-to-go function)
    
  }
  
  # Compute the optimal solution
  
  z[n] <- p[n]  # Terminal solution
  
  for (i in n1:1){
    z[i] <- min(z[i+1], p[i])
  }
  
  aux <- list("yf" = z, "ord" = o, "min" = p)  # Follows the isoreg() notation from "stats"
  
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
# The output can be used as input for bivdp() or biviso() from package "Iso".
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


# bivc(y)
#
# Function to correct the output of biv() before applying bivdp() when we get 
# empty values in some cells of the grid.
#
# It replaces the empty values of the bivariate grid y with the last values in
# row or column, following a carry-on appending strategy (Wen and Zeng, 1999).
#

bivc <- function(y){

  m <- nrow(y)        # Number of rows
  n <- ncol(y)        # Number of columns

  for (j in 1:n){
    for (i in 2:m){
      if (y[i,j] == 0 & y[i-1,j] != 0) y[i,j] <- y[i-1,j]
    }
  }
  
  for (i in 1:m){
    for (j in 2:n){
      if (y[i,j] == 0 & y[i,j-1] != 0) y[i,j] <- y[i,j-1]
    }
  }

  return(y)
}


# bivdp(y)
#
# Bivariate unweighted isotonic L1 regression with respect to a simple (increasing)
# order on both variables (Fernandes, 2026). 
# 
# y is the matrix (2D grid) of observations to be isotonized.
#

bivdp <- function(y){
  
  m <- nrow(y)        # Number of rows
  n <- ncol(y)        # Number of columns
  
  P <- matrix(0,m,n)  # arg min of each cost-to-go function
  Z <- matrix(0,m,n)  # Optimal solution
  
  Q <- matrix(-999999999,m*n+1,2) # Priority queue of breakpoints (position and value in column)
  
  # Anti-diagonal traversal of matrix P from d = 0 up to (m+n-2)  
  
  k <- 1
  
  for (d in 0:(m+n-2)){
    
    i_min <- max(0, d-(n-1))
    i_max <- min(d, m-1)
    
    for (i in i_min:i_max){
      
      j <- d-i  # NB: in anti-diagonals, we have i + j = d 
    
      # Insert a new breakpoint with position y[i,j] and value 2
      Q[k,1] <- y[i+1,j+1]  # convert to 1-indexed for R
      Q[k,2] <- 2
      
      # Sort the queue by position/priority (descending order)  
      Q <- Q[order(-Q[,1]),]
      
      # find/peek max i.e. the breakpoint with highest position/priority (Q head)
      B <- Q[1,2]
      
      if (B==1){
        # The max is from a previous interaction, so "delete" it (pop)
        Q[1,1] <- -999999999
        Q[1,2] <- -999999999
        
        # Sort the queue by position/priority (descending order) once again
        Q <- Q[order(-Q[,1]),]
      }
      else{
        # update the value with max priority (Q head)
        Q[1,2] <- 1
      }
      
      P[i+1,j+1] <- Q[1,1]  # Save the highest position (arg min of k cost-to-go function)
      
      k <- k+1  # Increments the queue
    }
  }
  
  # Compute the optimal solution Z
  
  Z[m,n] <- P[m,n]  # Terminal solution
  
  for (i in (m-1):1){
    Z[i,n] <- min(Z[i+1,n], P[i,n])  # Last column
  }
  
  for (j in (n-1):1){
    Z[m,j] <- min(Z[m,j+1], P[m,j])  # Last row
  }
  
  # Anti-diagonal traversal of matrix Z's "middle" from d = (m+n-4) down to 0  
  
  for (d in (m+n-4):0){
    
    i_min <- max(0, d-(n-2))
    i_max <- min(d, m-2)
    
    for (i in i_max:i_min){
      
      j <- d-i  # NB: in anti-diagonals, we have i + j = d
    
      Z[i+1,j+1] <- min(Z[i+2,j+1], Z[i+1,j+2], P[i+1,j+1])
    }
  
  }
  
  aux <- list("yf" = Z, "min" = P)
  
  return(aux)
  
}


# slogis(x)
#
# Computes the logistic scale parameter given a vector x using the formula: 
# 
# scale = sqrt(3) * StdDev / pi
#

slogis <- function(x){

  scale <- sd(x) * sqrt(3) / pi
  
}


# alogis(x)
#
# Computes the logistic correlation given two vectors x and y using the formula: 
# 
# alpha = cor(x,y) * pi*pi / 3
#

alogis <- function(x,y){
  
  alpha <- cor(x,y) * pi*pi / 3
  
}


# pbivlogis(x, y, loc_x, scale_x, loc_y, scale_y)
#
# Standard Type I bivariate logistic distribution (Gumbel, 1961): 
# 
# H(x,y) = 1 / (1 + exp(-x) + exp(-y))
#

pbivlogis <- function(x, y, loc_x=0, scale_x=1, loc_y=0, scale_y=1){
  
  s_x <- (x - loc_x) / scale_x 
  s_y <- (y - loc_y) / scale_y
  
  bivdist <- 1 / (1 + exp(-s_x) + exp(-s_y))

  }


# pbivlogisFGM(x, y, loc_x, scale_x, loc_y, scale_y, alpha)
#
# Standard Type II bivariate logistic distribution (Gumbel, 1961): 
# 
# H(x,y) = F(x) F(y) [1 + alpha (1-F(x)) (1-F(y))]
#
# where F(.) is the standard logistic distribution.
#

pbivlogisFGM <- function(x, y, loc_x=0, scale_x=1, loc_y=0, scale_y=1, alpha=0){
  
  s_x <- (x - loc_x) / scale_x 
  s_y <- (y - loc_y) / scale_y
  
  F_x <- 1 / (1 + exp(-s_x))
  F_y <- 1 / (1 + exp(-s_y))
  
  bivdist <- F_x * F_y * (1 + alpha * (1-F_x) *  (1-F_y))
  
}


# fbivlogis(x, y, loc_x, scale_x, loc_y, scale_y, loc_z, scale_z)
#
# Simple forecasting of variable z conditional to the standard Type I bivariate
# logistic distribution of x and y.
# 

fbivlogis <- function(x, y, loc_x=0, scale_x=1, loc_y=0, scale_y=1, loc_z=0, scale_z=1){
  
  h <- pbivlogis(x, y, loc_x, scale_x, loc_y, scale_y)
  
  f <- qlogis(h, location = loc_z, scale = scale_z)
  
  f <- f + loc_z - mean(f)  # To assure that mean(forecast) = mean(target)
  
}


# fbivlogisFGM(x, y, loc_x, scale_x, loc_y, scale_y, loc_z, scale_z)
#
# Simple forecasting of variable z conditional to the standard Type II bivariate
# logistic distribution of x and y.
# 

fbivlogisFGM <- function(x, y, loc_x=0, scale_x=1, loc_y=0, scale_y=1, loc_z=0, scale_z=1){
  
  h <- pbivlogisFGM(x, y, loc_x, scale_x, loc_y, scale_y)
  
  f <- qlogis(h, location = loc_z, scale = scale_z)
  
  f <- f + loc_z - mean(f)  # To assure that mean(forecast) = mean(target)
  
}






