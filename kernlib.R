#
# kernlib.R - A small library of kernel methods in R, v01.01 (2026-04-09)
#
# Pedro Afonso Fernandes, UCP, CLSBE, Lisbon, Portugal (paf@ucp.pt)
# 
# This library is free software; you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation. The library is
# distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details: https://www.gnu.org/licenses/.
#

# kgauss(x)
#
# Gaussian kernel indicated by Friedman (2006).

kgauss <- function(x){
  
  y <- exp(-x*x/2)
  
  return(y)
}


# lcrg(x, y, k)
#
# Computes the local constant regression of series y on x with k lags and leads.
#
# This is the classic Nadaraya-Watson estimator with Gaussian kernel.
#
# Uses the 'first and last values carry-on appending strategy' (Wen and Zeng, 1999).
#
# It returns a time series object with the original data y filtered.

lcrg <- function(x, y, k = 11){
  
  n <- length(y)
  m <- n + 2 * k
  
  a <- rep(0, m)
  b <- rep(0, m)
  
  for (i in 1:m) {
    
    if (i <= k){
      a[i] <- x[1]
      b[i] <- y[1]
    }
    else if(i > (n+k)){
      a[i] <- x[n]
      b[i] <- y[n]
    }
    else{
      a[i] <- x[i-k]
      b[i] <- y[i-k]
    }
  }
  
  for (i in 1:n) {
    
    h <- max(a[i:(i+2*k)]) - min(a[i:(i+2*k)])  # Local bandwidth (range)
    
    w <- kgauss( (a[i:(i+2*k)]-x[i])/h )  # Local weights
    
    s <- sum(w)  # Sum of local weights
    
    x[i] <- (b[i:(i+2*k)] %*% w) / s
  }
  
  out <- ts(data = x, start = 1, end = n)
}
