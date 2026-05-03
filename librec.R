#
# librec.R - A small library of recursive methods in R, v01.22 (2026-04-03)
#
# Adaptation to R with extensions of the Thomas Sargent's MatLab toolkits available at:
#
#  http://www.tomsargent.com/source_code/mitbook.zip 
#  http://www.tomsargent.com/source_code/hansarprograms.zip
#
# Pedro Afonso Fernandes, UCP, CLSBE, Lisbon, Portugal (paf@ucp.pt)
# 
# This library is free software; you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation. The library is
# distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details: https://www.gnu.org/licenses/.
#


# Load other (auxiliary) libraries

library(MASS)



### THOMAS SARGENT'S FUNCTIONS ADAPTED TO R ###


# compn(a)
#
# Creates companion matrix B for the vector a such that
#
#          | a(1) a(2) ... a(n-1) a(n) |
#          |   1   0         0      0  |
#      B = |   0   1         0      0  |
#          |   .                       |
#          |   .                       |
#          |   .                       |
#          |   0   0         1      0  |
#
# The matrix B is returned

compn <- function(a){

    n <- length(a)
    m <- n-1
    
    B <- rbind(a, cbind(diag(m), rep(0,m)))

    return(B)
}


# cmean(A,G,x0,N)
#
# In the framework of a discrete-time linear state-space (lss) system, computes the conditional mean:
#
#   E[y(t) | x(0)] = G * A^t * x0 for t = 1,2,...,N
#
# where:
#
#   A is a n*n transition matrix of the lss system
#   G is a m*n output matrix of the lss system
#   x0 is a n*1 vector with the initial value of the state (condition for the expected value)
#   N is the time horizon (integer)
#
# NB: this function can be used to predict the expected values of the observations in a lss system
#     or to simulate the control in the optimal linear regulator problem (OLRP) by making G = -F
#     and A = A0 = (A-BF)
#
# NB: this function can be used to simulate the state in the optimal linear regulator problem (OLRP) 
#     by making G = diag(n) and A = A0 = (A-BF)
#
# It returns a m*(1+N) matrix B with first column B[,1] = G * x0

cmean <- function(A,G,x0,N){

    m <- nrow(G)

    B1 <- matrix(0,m,N)

    for (i in 1:N){

        B1[,i] <- G %*% mpower(A,i) %*% x0
    }

    B <- cbind(G %*% x0, B1)

    return(B)
}


# dlsim(A,C,G,D,w,x0)
#
# Simulate the response of the discrete-time linear state-space (lss) system:
#
#   x(t+1) = A * x(t) + C * w(t+1)
#   y(t)   = G * x(t) + D * w(t+1)
#
# to a N*p matrix of perturbations/shocks w (given) where:
#
#   N is the number of simulations (obs)
#   x(t) is a n*1 vector of the states with an initial condition x0 (given)
#   A is the n*n transition matrix of the lss system
#   G is the m*n output matrix of the lss system
#   C is a n*p state's volatility matrix
#   D is a m*p output's volatility matrix
#
# It returns the m time series y(t) for t = 0, ..., N with y(0) = t(G * x0)

dlsim <- function(A,C,G,D,w,x0){

    n <- nrow(A)
    m <- nrow(G)

    N <- nrow(w)

    y1 <- matrix(0,m,N)
    x1 <- matrix(0,n,N)

    x_0 <- x0

    for (i in 1:N){

        x_1 <- A %*% x_0 + C %*% t(w[i,])

        y1[,i] <- G %*% x_1 + D %*% t(w[i,])

        x1[,i] <- x_1

        x_0 <- x_1

    }

    y <- cbind(G %*% x0, y1)

    ty <- ts(data = t(y), start = 0, end = N)

}


# dimpulse(A,C,G,D,wj,N)
#
# Impulse response of the discrete-time linear state-space (lss) system:
#
#   x(t+1) = A * x(t) + C * w(t+1)
#   y(t)   = G * x(t) + D * w(t+1)
#
# to an unit shock applied to the wj perturbation with 1 <= wj <= p where p
# is the number of columns of the matrices C and D
#
# Integer N specifies how many points of the impulse response to find
#
# It returns a time series object y(t), t = 0, ..., N with m columns

dimpulse <- function(A,C,G,D,wj,N){

	n <- nrow(A)

	x0 <- as.matrix(rep(0,n))	

	w  <- as.matrix(rep(0,N))
	w[1] <- 1

	I <- dlsim(A,as.matrix(C[,wj]),G,as.matrix(D[,wj]),w,x0)

}


# doublej(A,B)
#
# Solves the discrete Lyapunov equation V = A * V * A' + B using the Sargent's "doubling algorithm"
# implemented with the routine doublej.m from http://www.tomsargent.com/source_code/mitbook.zip
#
# This algorithm computes V = SUM (A^j) * B * t(A^j) from j = 0 to j = infinity

doublej <- function(A,B){

    A0 <- A
    V0 <- B

    dif   <- 9999
    i    <- 1
    maxit <- 1000

    while (dif > 1e-15) {

        A1 <- A0 %*% A0

        V1 <- V0 + A0 %*% V0 %*% t(A0)

        dif <- max(abs(V1-V0))

        V0 <- V1
        A0 <- A1
        i  <- i + 1

        if (i > maxit) break

    }

    return(V1)
}


# kfilter(A,Q,G,R)
#
# Calculates the time-invariant gain, K, and the stationary covariance matrix, S, using the 
# Kalman filter for the linear state space (lss) system:
#
#   x(t+1)  = A x(t) + C w(t+1)
#   y(t)    = G x(t) + v(t)
#
# where:
#
#   x(t) is the n*1 vector of unobserved states at time t = 0,1,2,...
#   A is the n*n transition matrix of the lss system
#   w(t+1) is an iid sequence of p*1 Gaussian random numbers ~ N(O,I); typically, p=1
#   C is the n*p state's volatility matrix
#
#   y(t) is the m*1 vector of observations at time t; typically, m=1
#   G is the m*n output matrix
#   v(t+1) is an iid sequence of m*1 Gaussian random numbers ~ N(O,R)
#   R is the m*m observation's volatility matrix; typically, a single number
#
# NB: the argument C should be passed as Q = C %*% t(C), a matrix n*n; i.e. C*w[t+1] ~ N(0,Q=CC')
#
# Based on the Sargent's routine kfilter.m from http://www.tomsargent.com/source_code/mitbook.zip
# that iterates on the Riccati diference equation:
#
# S1 = Q + A*S0*A' - A*S0*G' / (G*S0*G'+ R) * G*S0*A' = Q + A*S0*A' - K0 * G*S0*A' 
#
# starting from a initial covariance matrix S0 = 0.01 * I(n) and noting that the Kalman gain is
#
# K0 = A * S0 * G' / (G * S0 * G'+ R)
#
# It returns a bundle with the time-invariant gain K, the stationary covariance matrix S and 
# the number of steps i until convergence

kfilter <- function(A,Q,G,R){

    n <- nrow(A)

    S0 <- 0.01 * diag(n)

    dif <- 9999
    i    <- 1
    maxit <- 1000

    while (dif > 1e-8) {

        K0 <- A %*% S0 %*% t(G) %*% solve(G %*% S0 %*% t(G) + R)    # Kalman gain

        S1 <- Q + A %*% S0 %*% t(A) - K0 %*% G %*% S0 %*% t(A)      # Riccati dif equation

        K1 <- A %*% S1 %*% t(G) %*% solve(G %*% S1 %*% t(G) + R)    # Kalman gain (next iteration)

        dif <- max(abs(K1-K0))

        S0 <- S1

        i  <- i + 1

        if (i > maxit) break

    }

    aux <- list("K" = K1, "S" = S1, "i" = i)

    return(aux)

}


# ikalman(A,Q,G,R)
#
# Improved version of the Kalman filter that uses a Howard policy improvement algorithm as suggested by
# Ljungqvist and Sargent (2018, p. 103), Exercice 2.29, part e)
#
# Instead of iterating on the Riccati diference equation, this algorithm starts from a guess K0 to the 
# Kalman gain (the "policy function" in this context) and then solves the following Lyapunov equation for S:
#
# S = (A - K0 * G) * S * (A - K0 * G)' + (Q + K0 * R * K0')
#
# Given this time-invariant covariance matrix S, it recomputes the Kalman gain:
#
# K1 = A * S * G' / (G * S * G'+ R)
#
# and continues by solving again the Lyapunov equation given K1 until convergence between K1 and K0 
#
# It returns a bundle with the time-invariant gain K, the stationary covariance matrix S and 
# the number of steps i until convergence


ikalman <- function(A,Q,G,R){

    n <- nrow(A)
    m <- nrow(G)

    K0 <- matrix(0.01,n,m)

    dif <- 9999
    i    <- 1
    maxit <- 1000

    while (dif > 1e-8) {

        S <- doublej(A - K0 %*% G, Q + K0 %*% R %*% t(K0))      # Lyapunov equation

        K1 <- A %*% S %*% t(G) %*% solve(G %*% S %*% t(G) + R)  # Kalman gain (next iteration)

        dif <- max(abs(K1-K0))

        K0 <- K1

        i  <- i + 1

        if (i > maxit) break

    }

    aux <- list("K" = K1, "S" = S, "i" = i)

    return(aux)

}


# krec(A,Q,G,R,y,x0)
#
# Computes the Kalman recursions:
#
# xhat(t+1) = A * xhat(t) + K * a(t) with the inovation a(t) = y(t) - G * xhat(t)
#
# where:
#
#   A is the n*n transition matrix of the linear state space (lss) system
#   Q = CC' is the n*n state's volatility matrix
#   G is the m*n output matrix
#   R is the m*m observation's volatility matrix
#   K is the Kalman gain computed with the aux function kfilter(A,Q,G,R)
#   y is a N*m matrix of observations
#   x0 is a n*1 vector with the initial value of the state
#
# It returns a time series with the original N*m series filtered, i.e., yhat(t) = G * xhat(t)

krec <- function(A,Q,G,R,y,x0){

    RS <- kfilter(A,Q,G,R)

    n <- nrow(A)
    N <- nrow(y)

    B1 <- matrix(0,n,N)

    ty <- t(y)

    xhat_0 <- x0

    for (i in 1:N){

        a_0 <- ty[,i] - G %*% xhat_0

        xhat_1 <- A %*% xhat_0 + RS$K %*% a_0

        B1[,i] <- xhat_1

        xhat_0 <- xhat_1

    }

    B <- G %*% B1

    T <- ts(data = t(B), start = 1, end = N)

    return(T)

}


# markov(T,n,s0,v)
#
# Generates a Markov chain where:
#
#   T is a transition matrix with size m*m
#   n is the number of periods to simulate + initial period (n > 1)
#   s0 is the initial state (integer number from 1 to m)
#   v is a row vector with the quantity (or value) corresponding to each state i = 1,...,m
#
# It returns a list with:
#
#   Chain of values with size 1*n
#   Matrix of states with size m*n

markov <- function(T,n,s0,v){

    m <- nrow(T)
    X <- runif(n-1)

    s <- matrix(0,m,1)
    s[s0] <- 1

    S <- matrix(0,m,n)
    S[,1] <- s

    A <- matrix(1,m,m)
    A[lower.tri(A)] <- 0

    cum <- T %*% A  # Cumulative probabilities of transition from state/row i to state/column j

    for (k in 1:(n-1)){

        aux <- t(s) %*% cum

        ppi <- cbind(0, aux)

        s <- as.numeric(X[k] <= ppi[2:(m+1)]) * as.numeric(X[k] > ppi[1:m])

        S[,k+1] <- s
    }

    c <- v %*% S

    list("chain" = c, "states" = S)

}


# sims(A,G,C,x0,N)
#
# In the framework of a discrete-time linear state-space (lss) system, simulates a sequence {y(t)} of observations such that:
#
#   y(t) = G * x(t) with x(t) = A * x(t-1) + C * w(t) for t = 1,2,...,N
#
# where:
#
#   x(t) is a n*1 vector of unobserved states with an initial condition x(0) = x0
#   A is the n*n transition matrix of the lss system
#   G is the m*n output matrix of the lss system
#   C is a n*1 state's volatility matrix
#   w(t+1) is an iid sequence of Gaussian random numbers ~ N(O,I), generated by the function
#   N is the time horizon (integer)
#
# It returns a m*(1+N) matrix B with B[,1] = G*x0

sims <- function(A,G,C,x0,N){

    m <- nrow(G)

    B1 <- matrix(0,m,N)

    w <- rnorm(N, mean=0, sd=1)

    x_0 <- x0

    for (i in 1:N){

        x_1 <- A %*% x_0 + C * w[i]

        B1[,i] <- G %*% x_1

        x_0 <- x_1

    }

    B <- cbind(G %*% x0, B1)

    return(B)
}


# olrp(beta,A,B,R,Q,H)
#
# Solves the discounted Optimal Linear Regulator Problem (OLRP):
#
# Maximize { sum [beta^t (x'Rx + u'Qu + 2u'Hx)] }
#
# subject to
#
#   x(t+1) = A x(t) + B u(t)
#
# where:
#
#   x(t) is a n*1 vector of states at time t = 0,1,2,... with x(0) given
#   u(t) is a k*1 vector of controls
#   A is the n*n transition matrix associated with the states
#   B is the n*k transition matrix associated with the controls
#   R is a n*n positive semidefinite symmetric matrix 
#   Q is a k*k positive definite symmetric matrix
#   H is a k*n cross-product matrix
#   beta is the discount factor (0.96 by default)
#
# by iterating on the following Riccati diference equation:
#
# P1 = R + beta * A'*P0*A - (beta * A'*P0*B + H') / (Q + beta * B'*P0*B) * (beta * B'*P0*A + H)
#
# starting from an initial matrix P0 = 0
#
# The optimal value function will be x'Px associated with the optimal policy function u = -Fx where 
#
# F = (Q + beta * B' * P * B) \ (beta * B' * P * A + H)
#
# It returns a bundle with the optimal policy matrix F, the steady-state solution P of the Riccati dif 
# equation, the optimal closed loop transition matrix A0 = A - BF and the number of steps i until convergence

olrp <- function(beta = 0.96,A,B,R,Q,H){ 

    n <- ncol(A)

    P0 <- matrix(0.01,n,n)

    F0 <- solve(Q + beta * t(B) %*% P0 %*% B) %*% (beta * t(B) %*% P0 %*% A + H)

    dif <- 9999
    i    <- 1
    maxit <- 1000

    while (dif > 1e-6) {

        # Ricatti dif equation:

        P1 <- R + beta * t(A) %*% P0 %*% A - (beta * t(A) %*% P0 %*% B + t(H)) %*% F0

        F1 <- solve(Q + beta * t(B) %*% P1 %*% B) %*% (beta * t(B) %*% P1 %*% A + H)

        dif <- max(abs(F1-F0))

        P0 <- P1
        F0 <- F1

        i  <- i + 1

        if (i > maxit) break

    }

    A0 <- A - B %*% F1

    list("F" = F1, "P" = P1, "A0" = A0, "i" = i)

}


# policyi(beta,A,B,R,Q,H)
#
# Howard policy improvement algorithm that computes the matrix F from the feedback rule u = -Fx for the
# Optimal Linear Regulator Problem (OLRP):
#
# Maximize { sum [beta^t (x'Rx + u'Qu + 2u'Hx)] }
#
# subject to
#
#   x(t+1) = A x(t) + B u(t)
#
# where:
#
#   x(t) is a n*1 vector of states at time t = 0,1,2,... with x(0) given
#   u(t) is a k*1 vector of controls
#   A is the n*n transition matrix associated with the states
#   B is the n*k transition matrix associated with the controls
#   R is a n*n positive semidefinite symmetric matrix 
#   Q is a k*k positive definite symmetric matrix
#   H is a k*n cross-product matrix
#   beta is the discount factor (0.96 by default)
#
# Instead of iterating on the Riccati diference equation, this algorithm starts from a guess F0 to the 
# policy function and then solves the following Lyapunov equation for P:
#
# P = (R + F0' * Q * F0 - 2 * F0' * H) + beta * (A - B * F0)' * P * (A - B * F0)
#
# Given this steady-state solution P, it recomputes the policy function:
#
# F1 = (Q + beta * B' * P * B) \ (beta * B' * P * A + H)
#
# and continues by solving again the Lyapunov equation given F1 until convergence between F1 and F0
#
# The optimal value function will be x'Px associated with the optimal policy function u = -Fx 
#
# It returns a bundle with the optimal policy matrix F, the steady-state solution P of the Riccati dif 
# equation, the optimal closed loop transition matrix A0 = A - BF and the number of steps i until convergence


policyi <- function(beta = 0.06,A,B,R,Q,H){

    k <- ncol(B)
    n <- ncol(A)

    s <- rnorm(k*n,0,0.01)

    F0 <- matrix(s,k,n)

    dif <- 9999
    i    <- 1
    maxit <- 1000

    while (dif > 1e-6) {

        U <- t(A - B %*% F0) * sqrt(beta) 
        V <- R + t(F0) %*% Q %*% F0 - 2 * t(F0) %*% H

        P <- doublej(U, V)      # Lyapunov equation

        F1 <- solve(Q + beta * t(B) %*% P %*% B) %*% (beta * t(B) %*% P %*% A + H)  # Improved policy (next iteration)

        dif <- max(abs(F1-F0))

        F0 <- F1

        i  <- i + 1

        if (i > maxit) break

    }

    A0 <- A - B %*% F1

    list("F" = F1, "P" = P, "A0" = A0, "steps" = i)

}


# nash(beta,A,B1,B2,R1,R2,Q1,Q2,S1,S2,W1,W2,M1,M2)
#
# Computes the limit of a Nash-Markov linear quadratic two-player game by iterating on a pair of Ricatti equations.
# This is a complex version of the algorithm, with cross terms between states and controls
#
# Each player i:
#
# Maximize { sum [beta^t (x' Ri x + ui' Qi ui + uj' Si uj + 2 x' Wi ui + 2 uj' Mi ui)] }
#
# subject to
#
#   x(t+1) = A x(t) + Bi ui(t) + Bj uj(t)
#
# and a perceived control uj(t) = - Fj x(t) for the other played j.
#
# where:
#
#   x(t) is a n*1 vector of states at time t = 0,1,2,... with x(0) given
#   ui(t) is a ki*1 vector of controls
#   A is the n*n transition matrix associated with the states
#   Bi is the n*ki transition matrix associated with the controls of player i
#   Ri is a n*n positive semidefinite symmetric matrix 
#   Qi is a ki*ki positive definite symmetric matrix
#   Si is a kj*kj positive definite symmetric matrix
#   Wi is n x ki
#   Mi is kj x ki
#   beta is the discount factor (0.96 by default)
#

nash <- function(beta=0.96,A,B1,B2,R1,R2,Q1,Q2,S1,S2,W1,W2,M1,M2){

    n <- ncol(A)    # Number of states

    k1 <- ncol(Q1)  # Number of controls of player 1
    k2 <- ncol(Q2)  # Number of controls of player 2

    P1 <- matrix(0,n,n)
    P2 <- matrix(0,n,n)

    v1 <- rnorm(k1*n,0,0.01)
    v2 <- rnorm(k2*n,0,0.01)

    F1 <- matrix(v1,k1,n)
    F2 <- matrix(v2,k2,n)

    H1 <- t(W1) - t(M1) %*% F2
    H2 <- t(W2) - t(M2) %*% F1

    dif <- 9999
    i    <- 1
    maxit <- 1000

    while (dif > 1e-6) {

        F10 <- F1
        F20 <- F2

        F1 <- solve(Q1 + beta * t(B1) %*% P1 %*% B1) %*% (beta * t(B1) %*% P1 %*% (A - B2 %*% F20) + H1)
        F2 <- solve(Q2 + beta * t(B2) %*% P2 %*% B2) %*% (beta * t(B2) %*% P2 %*% (A - B1 %*% F10) + H2)

        A1 <- A - B2 %*% F2
	      A2 <- A - B1 %*% F1

        H1 <- t(W1) - t(M1) %*% F2
        H2 <- t(W2) - t(M2) %*% F1

        # Ricatti dif equation:

        P1 <- (R1 + t(F2) %*% S1 %*% F2) + beta * t(A1) %*% P1 %*% A1 - (beta * t(A1) %*% P1 %*% B1 + t(H1)) %*% F1
        P2 <- (R2 + t(F1) %*% S2 %*% F1) + beta * t(A2) %*% P2 %*% A2 - (beta * t(A2) %*% P2 %*% B2 + t(H2)) %*% F2

        dif <- max(abs(F10-F1))+max(abs(F20-F2))

        i  <- i + 1

        if (i > maxit) break

    }

    A0 <- A - B1 %*% F1 - B2 %*% F2

    list("F1" = F1, "F2" = F2, "P1" = P1, "P2" = P2, "A0" = A0, "i" = i)

}


# nash2(beta,A,B1,B2,R1,R2,Q1,Q2,S1,S2)
#
# Computes the limit of a Nash-Markov linear quadratic two-player game by iterating on a pair of Ricatti equations.
# This is a simple version of the algorithm, without cross terms between states and controls
#
# Each player i:
#
# Maximize { sum [beta^t (x' Ri x + ui' Qi ui + uj' Si uj)] }
#
# subject to
#
#   x(t+1) = A x(t) + Bi ui(t) + Bj uj(t)
#
# and a perceived control uj(t) = - Fj x(t) for the other played j.
#
# where:
#
#   x(t) is a n*1 vector of states at time t = 0,1,2,... with x(0) given
#   ui(t) is a ki*1 vector of controls
#   A is the n*n transition matrix associated with the states
#   Bi is the n*ki transition matrix associated with the controls of player i
#   Ri is a n*n positive semidefinite symmetric matrix 
#   Qi is a ki*ki positive definite symmetric matrix
#   Si is a kj*kj positive definite symmetric matrix
#   beta is the discount factor (0.96 by default)
#

nash2 <- function(beta=0.96,A,B1,B2,R1,R2,Q1,Q2,S1,S2){

    n <- ncol(A)    # Number of states

    k1 <- ncol(Q1)  # Number of controls of player 1
    k2 <- ncol(Q2)  # Number of controls of player 2

    P1 <- matrix(0,n,n)
    P2 <- matrix(0,n,n)

    v1 <- rnorm(k1*n,0,0.01)
    v2 <- rnorm(k2*n,0,0.01)

    F1 <- matrix(v1,k1,n)
    F2 <- matrix(v2,k2,n)

    dif <- 9999
    i    <- 1
    maxit <- 1000

    while (dif > 1e-6) {

        F10 <- F1
        F20 <- F2

        F1 <- solve(Q1 + beta * t(B1) %*% P1 %*% B1) %*% (beta * t(B1) %*% P1 %*% (A - B2 %*% F20))
        F2 <- solve(Q2 + beta * t(B2) %*% P2 %*% B2) %*% (beta * t(B2) %*% P2 %*% (A - B1 %*% F10))

        A1 <- A - B2 %*% F2
		    A2 <- A - B1 %*% F1

        # Ricatti dif equation:

        P1 <- (R1 + t(F2) %*% S1 %*% F2) + beta * t(A1) %*% P1 %*% A1 - (beta * t(A1) %*% P1 %*% B1) %*% F1
        P2 <- (R2 + t(F1) %*% S2 %*% F1) + beta * t(A2) %*% P2 %*% A2 - (beta * t(A2) %*% P2 %*% B2) %*% F2

        dif <- max(abs(F10-F1))+max(abs(F20-F2))

        i  <- i + 1

        if (i > maxit) break

    }

    A0 <- A - B1 %*% F1 - B2 %*% F2

    list("F1" = F1, "F2" = F2, "P1" = P1, "P2" = P2, "A0" = A0, "i" = i)

}


# nash2s(beta,A,B1,B2,R1,R2,Q1,Q2)
#
# Slim version of the nash2() function without matrices Si, that is, each player i = 1,2 simply: 
#
# Maximize { sum [beta^t (x' Ri x + ui' Qi ui)] }
#
# subject to
#
#   x(t+1) = A x(t) + Bi ui(t) + Bj uj(t)
#

nash2s <- function(beta=0.96,A,B1,B2,R1,R2,Q1,Q2){

    n <- ncol(A)    # Number of states

    k1 <- ncol(Q1)  # Number of controls of player 1
    k2 <- ncol(Q2)  # Number of controls of player 2

    P1 <- matrix(0,n,n)
    P2 <- matrix(0,n,n)

    v1 <- rnorm(k1*n,0,0.01)
    v2 <- rnorm(k2*n,0,0.01)

    F1 <- matrix(v1,k1,n)
    F2 <- matrix(v2,k2,n)

    dif <- 9999
    i    <- 1
    maxit <- 1000

    while (dif > 1e-6) {

        F10 <- F1
        F20 <- F2

        F1 <- solve(Q1 + beta * t(B1) %*% P1 %*% B1) %*% (beta * t(B1) %*% P1 %*% (A - B2 %*% F20))
        F2 <- solve(Q2 + beta * t(B2) %*% P2 %*% B2) %*% (beta * t(B2) %*% P2 %*% (A - B1 %*% F10))

        A1 <- A - B2 %*% F2
		A2 <- A - B1 %*% F1

        # Ricatti dif equation:

        P1 <- R1 + beta * t(A1) %*% P1 %*% A1 - (beta * t(A1) %*% P1 %*% B1) %*% F1
        P2 <- R2 + beta * t(A2) %*% P2 %*% A2 - (beta * t(A2) %*% P2 %*% B2) %*% F2

        dif <- max(abs(F10-F1))+max(abs(F20-F2))

        i  <- i + 1

        if (i > maxit) break

    }

    A0 <- A - B1 %*% F1 - B2 %*% F2

    list("F1" = F1, "F2" = F2, "P1" = P1, "P2" = P2, "A0" = A0, "i" = i)

}


# nash8s(beta,A,B1,B2,B3,B4,B5,B6,B7,B8,R1,R2,R3,R4,R5,R6,R7,R8,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8)
#
# Version of the nash2s() algorithm for 8 players where each player i
#
# Maximize { sum [beta^t (x' Ri x + ui' Qi ui)] }
#
# subject to
#
#   x(t+1) = A x(t) + B1 u1(t) + ... + Bi ui(t) + ... + B8 u8(t)
#

nash8s <- function(beta=0.96,A,B1,B2,B3,B4,B5,B6,B7,B8,R1,R2,R3,R4,R5,R6,R7,R8,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8){

    n <- ncol(A)    # Number of states

    k1 <- ncol(Q1)  # Number of controls of player 1
    k2 <- ncol(Q2)  # Number of controls of player 2
  	k3 <- ncol(Q3)  # Number of controls of player 3
	  k4 <- ncol(Q4)  # Number of controls of player 4
	  k5 <- ncol(Q5)  # Number of controls of player 5
	  k6 <- ncol(Q6)  # Number of controls of player 6
	  k7 <- ncol(Q7)  # Number of controls of player 7
	  k8 <- ncol(Q8)  # Number of controls of player 8

    P1 <- matrix(0,n,n)
    P2 <- matrix(0,n,n)
	  P3 <- matrix(0,n,n)
	  P4 <- matrix(0,n,n)
	  P5 <- matrix(0,n,n)
	  P6 <- matrix(0,n,n)
	  P7 <- matrix(0,n,n)
	  P8 <- matrix(0,n,n)

    v1 <- rnorm(k1*n,0,0.01)
    v2 <- rnorm(k2*n,0,0.01)
	  v3 <- rnorm(k3*n,0,0.01)
	  v4 <- rnorm(k4*n,0,0.01)
	  v5 <- rnorm(k5*n,0,0.01)
	  v6 <- rnorm(k6*n,0,0.01)
	  v7 <- rnorm(k7*n,0,0.01)
	  v8 <- rnorm(k8*n,0,0.01)

    F1 <- matrix(v1,k1,n)
    F2 <- matrix(v2,k2,n)
	  F3 <- matrix(v3,k3,n)
	  F4 <- matrix(v4,k4,n)
	  F5 <- matrix(v5,k5,n)
	  F6 <- matrix(v6,k6,n)
	  F7 <- matrix(v7,k7,n)
	  F8 <- matrix(v8,k8,n)

    dif <- 9999
    i    <- 1
    maxit <- 1000

    while (dif > 1e-6) {

        F10 <- F1
        F20 <- F2
		    F30 <- F3
		    F40 <- F4
		    F50 <- F5
		    F60 <- F6
		    F70 <- F7
		    F80 <- F8

        F1 <- solve(Q1 + beta * t(B1) %*% P1 %*% B1) %*% (beta * t(B1) %*% P1 %*% (A              - B2 %*% F20 - B3 %*% F30 - B4 %*% F40 - B5 %*% F50 - B6 %*% F60 - B7 %*% F70 - B8 %*% F80))
        F2 <- solve(Q2 + beta * t(B2) %*% P2 %*% B2) %*% (beta * t(B2) %*% P2 %*% (A - B1 %*% F10              - B3 %*% F30 - B4 %*% F40 - B5 %*% F50 - B6 %*% F60 - B7 %*% F70 - B8 %*% F80))
        F3 <- solve(Q3 + beta * t(B3) %*% P3 %*% B3) %*% (beta * t(B3) %*% P3 %*% (A - B1 %*% F10 - B2 %*% F20              - B4 %*% F40 - B5 %*% F50 - B6 %*% F60 - B7 %*% F70 - B8 %*% F80))
        F4 <- solve(Q4 + beta * t(B4) %*% P4 %*% B4) %*% (beta * t(B4) %*% P4 %*% (A - B1 %*% F10 - B2 %*% F20 - B3 %*% F30              - B5 %*% F50 - B6 %*% F60 - B7 %*% F70 - B8 %*% F80))
        F5 <- solve(Q5 + beta * t(B5) %*% P5 %*% B5) %*% (beta * t(B5) %*% P5 %*% (A - B1 %*% F10 - B2 %*% F20 - B3 %*% F30 - B4 %*% F40              - B6 %*% F60 - B7 %*% F70 - B8 %*% F80))
        F6 <- solve(Q6 + beta * t(B6) %*% P6 %*% B6) %*% (beta * t(B6) %*% P6 %*% (A - B1 %*% F10 - B2 %*% F20 - B3 %*% F30 - B4 %*% F40 - B5 %*% F50              - B7 %*% F70 - B8 %*% F80))
        F7 <- solve(Q7 + beta * t(B7) %*% P7 %*% B7) %*% (beta * t(B7) %*% P7 %*% (A - B1 %*% F10 - B2 %*% F20 - B3 %*% F30 - B4 %*% F40 - B5 %*% F50 - B6 %*% F60              - B8 %*% F80))
        F8 <- solve(Q8 + beta * t(B8) %*% P8 %*% B8) %*% (beta * t(B8) %*% P8 %*% (A - B1 %*% F10 - B2 %*% F20 - B3 %*% F30 - B4 %*% F40 - B5 %*% F50 - B6 %*% F60 - B7 %*% F70             ))

        A1 <- A             - B2 %*% F2 - B3 %*% F3 - B4 %*% F4 - B5 %*% F5 - B6 %*% F6 - B7 %*% F7 - B8 %*% F8
        A2 <- A - B1 %*% F1             - B3 %*% F3 - B4 %*% F4 - B5 %*% F5 - B6 %*% F6 - B7 %*% F7 - B8 %*% F8
        A3 <- A - B1 %*% F1 - B2 %*% F2             - B4 %*% F4 - B5 %*% F5 - B6 %*% F6 - B7 %*% F7 - B8 %*% F8
        A4 <- A - B1 %*% F1 - B2 %*% F2 - B3 %*% F3             - B5 %*% F5 - B6 %*% F6 - B7 %*% F7 - B8 %*% F8
        A5 <- A - B1 %*% F1 - B2 %*% F2 - B3 %*% F3 - B4 %*% F4             - B6 %*% F6 - B7 %*% F7 - B8 %*% F8
        A6 <- A - B1 %*% F1 - B2 %*% F2 - B3 %*% F3 - B4 %*% F4 - B5 %*% F5             - B7 %*% F7 - B8 %*% F8
        A7 <- A - B1 %*% F1 - B2 %*% F2 - B3 %*% F3 - B4 %*% F4 - B5 %*% F5 - B6 %*% F6             - B8 %*% F8
        A8 <- A - B1 %*% F1 - B2 %*% F2 - B3 %*% F3 - B4 %*% F4 - B5 %*% F5 - B6 %*% F6 - B7 %*% F7

        # Ricatti dif equation:

        P1 <- R1 + beta * t(A1) %*% P1 %*% A1 - (beta * t(A1) %*% P1 %*% B1) %*% F1
        P2 <- R2 + beta * t(A2) %*% P2 %*% A2 - (beta * t(A2) %*% P2 %*% B2) %*% F2
        P3 <- R3 + beta * t(A3) %*% P3 %*% A3 - (beta * t(A3) %*% P3 %*% B3) %*% F3
        P4 <- R4 + beta * t(A4) %*% P4 %*% A4 - (beta * t(A4) %*% P4 %*% B4) %*% F4
        P5 <- R5 + beta * t(A5) %*% P5 %*% A5 - (beta * t(A5) %*% P5 %*% B5) %*% F5
        P6 <- R6 + beta * t(A6) %*% P6 %*% A6 - (beta * t(A6) %*% P6 %*% B6) %*% F6
        P7 <- R7 + beta * t(A7) %*% P7 %*% A7 - (beta * t(A7) %*% P7 %*% B7) %*% F7
        P8 <- R8 + beta * t(A8) %*% P8 %*% A8 - (beta * t(A8) %*% P8 %*% B8) %*% F8

        dif <- max(abs(F10-F1)) + max(abs(F20-F2)) + max(abs(F30-F3)) + max(abs(F40-F4)) + max(abs(F50-F5)) + max(abs(F60-F6)) + max(abs(F70-F7)) + max(abs(F80-F8))

        i  <- i + 1

        if (i > maxit) break

    }

    A0 <- A - B1 %*% F1 - B2 %*% F2 - B3 %*% F3 - B4 %*% F4 - B5 %*% F5 - B6 %*% F6 - B7 %*% F7 - B8 %*% F8

    list("F1" = F1, "F2" = F2, "F3" = F3, "F4" = F4, "F5" = F5, "F6" = F6, "F7" = F7, "F8" = F8, "P1" = P1, "P2" = P2, "P3" = P3, "P4" = P4, "P5" = P5, "P6" = P6, "P7" = P7, "P8" = P8, "A0" = A0, "i" = i)

}



### TIME SERIES FILTERS AND FORECASTING METHODS ###

# hpfilt(y, lambda)
#
# Computes the Hodrick-Prescott (HP) trend component of series y for some lambda value, namely:
#
#   lambda = 6.25 or 100 for yearly data;
#   lambda = 1600 for quarterly data (default value);
#   lambda = 14400 for monthly data.
#
# ... following the computational strategy of Kim et al. (2009).
#
# It returns a time series object with the original data y filtered.

hpfilt <- function(y, lambda = 1600){

    n <- length(y)

    V <- rep(0, n)
    V[1] <- 1
    V[2] <- -2
    V[3] <- 1

    T <- toeplitz(V)
    T[lower.tri(T)] <- 0
    m <- n-2
    D <- T[1:m,]
    D2 <- crossprod(D)

    A <- solve( diag(n) + lambda * D2 )
    
    x <- A %*% y
    
    out <- ts(data = x, start = 1, end = n)
}


# bhpfilt(y, lambda, m)
#
# Computes the Boosted Hodrick-Prescott (HP) trend component of series y for some lambda value and 
# m iterations, following the computational strategy of Phillips and Shi (2021).
#
# It returns a time series object with the original data y filtered.

bhpfilt <- function(y, lambda = 1600, m = 1){

    n <- length(y)

    V <- rep(0, n)
    V[1] <- 1
    V[2] <- -2
    V[3] <- 1

    T <- toeplitz(V)
    T[lower.tri(T)] <- 0
    k <- n-2
    D <- T[1:k,]
    D2 <- crossprod(D)

    S <- solve( diag(n) + lambda * D2 )

    I <- diag(n)
    C <- I
    
    for (i in 1:m) {
      C <- C %*% (I - S)
    }
    
    x <- (I - C) %*% y
    
    out <- ts(data = x, start = 1, end = n)
}


# mafilt(y, k)
#
# Computes the symetric Moving Average (MA) of series y with k lags and leads.
#
# It returns a time series object with the original data y filtered.

mafilt <- function(y, k = 11){

    n <- length(y)
    m <- n + 2 * k

    a <- rep(0, m)
    x <- y

    lavg <- mean(x[1:(k+1)])
    uavg <- mean(x[(n-k):n])

    for (i in 1:m) {
	
	    if (i <= k){
	    	a[i] <- lavg
    	}
	    else if(i > (n+k)){
	    	a[i] <- uavg
	    }
	    else
	    	a[i] <- x[i-k]
    }

    for (i in 1:n) {
    	x[i] <- mean(a[i:(i+2*k)])
    }

    out <- ts(data = x, start = 1, end = n)
}


# mafilt2(y, k)
#
# Computes the symetric Moving Average (MA) of series y with k lags and leads.
#
# Uses the 'first and last values carry-on appending strategy' (Wen and Zeng, 1999)
#
# It returns a time series object with the original data y filtered.

mafilt2 <- function(y, k = 11){
  
  n <- length(y)
  m <- n + 2 * k
  
  a <- rep(0, m)
  x <- y
  
  for (i in 1:m) {
    
    if (i <= k){
      a[i] <- x[1]
    }
    else if(i > (n+k)){
      a[i] <- x[n]
    }
    else
      a[i] <- x[i-k]
  }
  
  for (i in 1:n) {
    x[i] <- mean(a[i:(i+2*k)])
  }
  
  out <- ts(data = x, start = 1, end = n)
}


# medfilt(y, k)
#
# Computes the Median (MED) filter of series y with k lags and leads, following the
# 'first and last values carry-on appending strategy' proposed by Wen and Zeng (1999).
#
# It returns a time series object with the original data y filtered.

medfilt <- function(y, k = 15){

    n <- length(y)
    m <- n + 2 * k

    a <- rep(0, m)
    x <- y

    for (i in 1:m) {
	
	    if (i <= k){
	    	a[i] <- x[1]
    	}
	    else if(i > (n+k)){
	    	a[i] <- x[n]
	    }
	    else
	    	a[i] <- x[i-k]
    }

    for (i in 1:n) {
    	x[i] <- median(a[i:(i+2*k)])
    }

    out <- ts(data = x, start = 1, end = n)
}


# bmedfilt(y, k, m)
#
# Computes the Boosted Median (MED) filter of series y with k lags and leads and 
# m iterations.
#
# It returns a time series object with the original data y filtered.

bmedfilt <- function(y, k = 15, m = 1){
  
  n <- length(y)
  
  D <- y
  
  for (i in 1:m) {
    T <- medfilt(D,k)
    C <- D - T
    D <- C
  }
  
  x <- y - C
  
  out <- ts(data = x, start = 1, end = n)
}


# mrfilt(y, theta)
#
# Computes the Mosheiov-Raveh (1997) filter using the LAD Meketon's algorithm, given the smooth parameter theta, namely:
#
#   theta = sqrt(100) = 10 for yearly data;
#   theta = srt(1600) = 40 for quarterly data (default value);
#   theta = sqrt(14400) = 120 for monthly data.
#
# It returns a time series object with the original data y filtered.

mrfilt <- function (y, theta = 40){

    n <- length(y)

    V <- rep(0, n)
    V[1] <- 1
    V[2] <- -2
    V[3] <- 1

    T <- toeplitz(V)
    T[lower.tri(T)] <- 0

    m <- n-2

    D <- T[1:m,]

    X <- rbind(diag(n), theta * D)

    y_tilde <- c(y, rep(0,m))

    aux <- meketon(y_tilde,X)

    out <- ts(data = aux$coef, start = 1, end = n)
}


# bmrfilt(y, theta, m)
#
# Computes the Boosted MR trend component of series y for some theta value and 
# m iterations.
#
# It returns a time series object with the original data y filtered.

bmrfilt <- function(y, theta = 40, m = 1){
  
  n <- length(y)
  
  D <- y
  
  for (i in 1:m) {
    T <- mrfilt(D,theta)
    C <- D - T
    D <- C
  }
  
  x <- y - C
  
  out <- ts(data = x, start = 1, end = n)
}


# hamfilt(y, p)
#
# Computes the Hamilton's (2017) cyclical component of series y with periodicity p = 4 by default (quarterly data).

hamfilt <- function(y, p = 4) {

    n <- length(y)
    h <- 2 * p
    m <- h + 4

    fm <- lm(y[m:n] ~ y[4:(n-h)] +  y[3:(n-h-1)] +  y[2:(n-h-2)] +  y[1:(n-h-3)])

    uhat <- as.matrix(residuals(fm))

    x <- c(rep(NA,times=(m-1)), uhat)

return(x)

}


# hampel(y, k, h)
#
# Detect the values outside the interval formed by the median plus or minus h median absolute deviations 
# (MAD) and replace them by those lower and upper bounds. The default threshold is 2.5 (Leys et al. 2013).
#
# The outliers' detection use a rolling window with k lags and leads, and follows the 'first and last 
# values carry-on appending strategy'. The default is k = 9 i.e. a windows of 19 periods (approx. 5 years
# with quarterly data -> good for COVID times).
#
# It returns a vector with the original data y corrected.

hampel <- function(y, k = 9, h = 2.5){
  
  n <- length(y)
  m <- n + 2 * k
  
  a <- rep(0, m)
  x <- y
  
  for (i in 1:m) {
    
    if (i <= k){
      a[i] <- x[1]
    }
    else if(i > (n+k)){
      a[i] <- x[n]
    }
    else
      a[i] <- x[i-k]
  }
  
  for (i in 1:n) {
    med <- median(a[i:(i+2*k)])
    mad <- mad(a[i:(i+2*k)], constant = 1)
    lb  <- med - h*mad
    ub  <- med + h*mad
    x[i] <- min(max(y[i],lb),ub)
    
  }
  
  out <- as.vector(x)
}



# theta(y,h,gamma)
#
# Returns the theta forecast for the time horizon h>0 given the time series y and the exponential smoothing parameter gamma.
#

theta <- function(y, h=1, gamma=0.3){
  
  n <- length(y)
  
  t <- 0:(n-1)
  
  mod0 <- lm(y ~ t)
  mod2 <- lm(-y ~ t)
  
  y0 <- mod0$coefficients[1] + mod0$coefficients[2]*t
  y2 <- mod2$coefficients[1] + mod2$coefficients[2]*t + 2*y
  
  t1 <- n:(n+h-1)
  
  yhat0 <- mod0$coefficients[1] + mod0$coefficients[2]*t1
  
  yhat2 <- rep(0,n)
  
  yhat2[1] <- y2[1]
  
  for (i in 2:n) yhat2[i] <- gamma * y2[i] + (1-gamma) * yhat2[i-1]
  
  yhat <- (yhat0 + yhat2[n]) * 0.5    

  }


# AR1(y)
#
# One-step look ahead forecast computed with a simple auto-regressive model of order 1.
#

AR1 <- function(y){
  
  n <- length(y)
  
  y0 <- y[1:(n-1)]
  y1 <- y[2:n]
  
  mod <- lm(y1 ~ y0)
  
  yhat <- mod$coefficients[1] + mod$coefficients[2]*y[n]

}


# AR1MAX(x,y)
#
# One-step look ahead forecast computed with a simple auto-regressive model of order 1 on y using the x covariates.
#
# NB: The vector y and the matrix x must have the same number of observations / rows
#

AR1MAX <- function(x,y){
  
  n <- length(y)
  
  y0 <- y[1:(n-1)]
  y1 <- y[2:n]
  yn <- y[n]
  
  x0 <- cbind(y0,x[1:(n-1),])
  xn <- c(yn,x[n,])
  
  mod <- lm(y1 ~ x0)
  
  yhat <- xn %*% coefficients(mod)[-1] + mod$coefficients[1]
  
}



### UTILITIES ###


# meketon(y,x,eps,beta)
#
# Estimates the least absolute deviations (LAD) coefficients of y on x using the Affine Scalling Algorithm
# proposed by Meketon (1986) and implemented by Koenker (2008)

meketon <- function (y, x, eps = 1e-04, beta = 0.97){

    f <- lm.fit(x,y)
    n <- length(y)
    w <- rep(0, n)
    d <- rep(1, n)
    its <- 0
    while(sum(abs(f$resid)) - crossprod(y, w) > eps)
    {
        its <- its + 1
        s <- f$resid * d
        alpha <- max(pmax(s/(1 - w), -s/(1 + w)))
        w <- w + (beta/alpha) * s
        d <- pmin(1 - w, 1 + w)^2
        f <- lm.wfit(x,y,d)
    }

    list(coef = f$coef, steps = its)
}


# mpower(A,k)
#
# Multiplicates a square A matrix k times

mpower <- function(A,k){

    n <- nrow(A)
    P <- diag(n)

    for (i in 1:k) P <- P %*% A

    return(P)
}


# nodes(X,R)
#
# Computes the k nodes of a neural network layer given two k*k matrices with tensor values (X) and 
# weights (R), using the ReLU activation function.
#

nodes <- function(X,R){
  
  h <- apply(X*R,1,sum)
  
  h[which(h<0)] <- 0
  
  return(h)
}



# dnodes(X,R)
#
# Computes the partial derivatives of the k nodes with respect to the weights of a neural network given two 
# k*k matrices with tensor values (X) and weights (R), using the ReLU activation function.
#

dnodes <- function(X,R){
  
  D <- X
  
  h <- as.vector(nodes(X,R))
  
  D[which(h==0),] <- 0 
  
  return(D)
}


# regm(X)
#
# Regularize a matrix of data by subtracting the min of each column and dividing by its range (max-min)
#

regm <- function(X){
  
  max   <- apply(X,2,max)
  min   <- apply(X,2,min)
  range <- max - min 
  
  R <- t((t(X) - min) / range) 
  
  return(R)
}


# regv(X)
#
# Regularize a vector of data by subtracting its min and dividing by its range (max-min)
#

regv <- function(x){
  
  max   <- max(x)
  min   <- min(x)
  range <- max - min 
  
  y <- (x - min) / range 
  
  return(y)
}



### REINFORCEMENT LEARNING ALGORITHMS ###

# td0(ts,p,d,w,beta,c)
#
# Simple version of the tabular temporal difference TD(0) algorithm proposed by Sutton (1988) and described by Szepesvari (2010, pp. 11-14).
#
# For a single state x with m discrete categories (values), it estimates the value or reward-to-go function V(x) from data, given:
#
#   - a sequence ts of n integer numbers from 1 to m, which one represents the state of the system (e.g. assets detained) at time t = 1, ..., n;
#   - an immediate reward p associated with the last state x (p = 10 by default);
#   - an adjustment cost d associated with the transition from the last state x to the next state y, that is, with u = y - x  (d = 1 by default);
#   - a sequence w of n random "disturbances", N(0,1) by default, or the given values for an exogenous (not controlled) state (e.g. wages at time t = 1, ..., n);
#   - a discount factor beta (= 0.96 by default);
#   - a parameter c > 0 for the step-size sequence coefficient alpha = c/t (c = 1 by default).
#
# For each transition between x and y, the algorithm computes:
#
#   u     <- y - x
#   delta <- p * x + w - d * u + beta * V(y) - V(x)
#   V(x)  <- V(x) + alpha * delta
#
# starting from the initial state x = ts[1].
#
# It returns a bundle (list) with: 
#
#   - the estimated value function V(x);
#   - the one-step lookahead policy mu(x) with the adjustment u that maximizes the immediate and future reward from each state x, using the estimated value V(x+u) as terminal.
#

td0 <- function(ts,p=10,d=1,w=1,beta=0.96,c=1){

    m <- max(s)
    n <- length(s)

    V <- as.matrix(rep(0,m))

    if (length(w) == 1) w <- rnorm(n,0,w)


    # Value function V(x)

    for (i in 1:(n-1)){

        x <- ts[i]
        y <- ts[i+1]
        u <- y - x

        delta <- p * x + w[i] - d * u + beta * V[y] - V[x]

        alpha <- c / i

        V[x] <- V[x] + alpha * delta

    }


    # Policy function mu(x)

    z <- as.matrix(rep(0,m))

    ones <- as.matrix(rep(1,m))

    X <- as.matrix(1:m) %*% t(ones)
    Y <- ones %*% t(as.matrix(1:m))

    A <- (p + d) * X + w[n] - d * Y + beta * ones %*% t(V)

    W <- apply(A,1,max)

    for (i in 1:m){
        for (j in 1:m){
            if (A[i,j] == W[i]) z[i] <- j
        }
    }

    mu <- z - as.matrix(1:m)

    list("V" = V, "mu" = mu)

}


# td0LFA(ts,p,d,w,beta,c)
#
# Simple temporal difference TD(0) algorithm with linear function approximation and quadratic polynomial basis (Szepesvari, 2010, pp. 12, 18).
#
# For a single real-valued state x, it estimates from data the linear value or reward-to-go function V(x) with quadratic polynomial basis phi(x), given:
#
#   - a time series ts of n real-valued numbers, where each observation represents the state of the system (e.g. stock of capital) at time t = 1, ..., n;
#   - an immediate reward p associated with the last state x (p = 1 by default);
#   - a quadratic adjustment cost d associated with the transition from the last state x to the next state y, that is, with u = y - x (e.g. investment);
#   - a sequence w of n random "disturbances", N(0,1) by default, or the given values for an exogenous (not controlled) state (e.g. depreciation with w < 0);
#   - a discount factor beta (= 0.96 by default);
#   - a parameter c > 0 for the step-size sequence coefficient alpha = c/t (c = 1 by default).
#
# For each transition between x and y, the algorithm computes:
#
#   u     <- y - x
#   delta <- p * x + w - 0.5 d * u^2 + r1 * (beta * y - x) + 0.5 r2 * (beta * y^2 - x^2)
#   r0 <- r0 + alpha * delta
#   r1 <- r1 + alpha * delta * x
#   r2 <- r2 + alpha * delta * x^2
#
# starting from the initial state x = ts[1].
#
# It returns a bundle (list) with: 
#
#   - the estimated parameters r0, r1 and r2 of the value function V(x);
#   - the linear coefficients a and b of the one-step look-ahead policy mu(x) = a + b * x.
#

td0LFA <- function(ts,p=1,d=1,w=0,beta=0.96,c=0.0000001){
  
  n <- length(ts)
  
  r0 <- 0
  r1 <- 0
  r2 <- 0
  
  if (length(w) == 1) w <- rnorm(n,0,w)
  
  
  # Value function V(x)
  
  for (i in 1:(n-1)){
    
    x <- ts[i]
    y <- ts[i+1]
    u <- y - x
    
    delta <- p * x + w[i] - d/2 * u^2 + r1 * (beta * y - x) + r2/2 * (beta * y^2 - x^2)
    
    alpha <- c / i
    
    r0 <- r0 + alpha * delta
    r1 <- r1 + alpha * delta * x
    r2 <- r2 + alpha * delta * x^2
    
  }
  
 
  # Policy function mu(x)
  
  a <- beta * r1 / (d - beta * r2)
  b <- beta * r2 / (d - beta * r2)
  
  list("r0" = r0, "r1" = r1, "r2" = r2, "a" = a, "b" = b)
  
}


# td0LFA2(ts,p,d,w,beta,c)
#
# Version 2 of the simple temporal difference TD(0) algorithm with the immediate reward as a function of u = y - x instead of x:
#
#   delta <- p * x + w - 0.5 d * u^2 + r1 * (beta * y - x) + 0.5 r2 * (beta * y^2 - x^2)
#

td0LFA2 <- function(ts,p=1,d=1,w=0,beta=0.96,c=0.0000001){
  
  n <- length(ts)
  
  r0 <- 0
  r1 <- 0
  r2 <- 0
  
  if (length(w) == 1) w <- rnorm(n,0,w)
  
  
  # Value function V(x)
  
  for (i in 1:(n-1)){
    
    x <- ts[i]
    y <- ts[i+1]
    u <- y - x
    
    delta <- p * u + w[i] - d/2 * u^2 + r1 * (beta * y - x) + r2/2 * (beta * y^2 - x^2)
    
    alpha <- c / i
    
    r0 <- r0 + alpha * delta
    r1 <- r1 + alpha * delta * x
    r2 <- r2 + alpha * delta * x^2
    
  }
  
 
  # Policy function mu(x)
  
  a <- (p + beta * r1) / (d - beta * r2)
  b <- beta * r2 / (d - beta * r2)
  
  list("r0" = r0, "r1" = r1, "r2" = r2, "a" = a, "b" = b)
  
}


# td0LFA3(z,ts,p,d,w,beta,c)
#
# Version 3 of the TD(0) algorithm with cross-effects between ts (x&y) and a matrix n*s of non-controllable states z:
#
#   delta <- p * u + w - 0.5 d * u^2 + r1 * (beta * y - x) + 0.5 r2 * (beta * y^2 - x^2) + (beta*y*z1 - x*z0) * r3
#
# where r3 is a column vector with s elements and z0/z1 are rows of z with s elements too.
#

td0LFA3 <- function(z,ts,p=1,d=1,w=0,beta=0.96,c=0.0000001){
  
  n <- length(ts)
  s <- ncol(z)
  
  r0 <- 0
  r1 <- 0
  r2 <- 0
  
  r3 <- rep(0,s) #matrix(0,s,1)
  
  if (length(w) == 1) w <- rnorm(n,0,w)
  
  # Value function V(x)
  
  for (i in 1:(n-1)){
    
    x <- ts[i]
    y <- ts[i+1]
    u <- y - x
    
    z0 <- z[i,]
    z1 <- z[i+1,]
    
    delta <- p * u + w[i] - d/2 * u^2 + r1 * (beta * y - x) + r2/2 * (beta * y^2 - x^2) + (beta*y*z1 - x*z0) %*% r3
    
    alpha <- c / i
    
    r0 <- r0 + alpha * delta
    r1 <- r1 + alpha * delta * x
    r2 <- r2 + alpha * delta * x^2
    
    r3 <- r3 + alpha * as.vector(delta) * as.vector(x) * as.vector(z0)
    
  }
  
  
  # Policy function mu(x)
  
  a <- beta * r1 / (d - beta * r2)
  b <- beta * r2 / (d - beta * r2)
  q <- beta * r3 / (d - beta * as.vector(2))
  
  list("r0" = r0, "r1" = r1, "r2" = r2, "r3" = r3, "a" = a, "b" = b, "q" = q)
  
}


# tdLambdaLFA(ts,p,d,w,beta,c,lambda)
#
# Simple temporal difference TD(lambda) algorithm with linear function approximation and quadratic polynomial basis (Szepesvari, 2010, pp. 16-18, 22).
#
# For each transition between x and y, the algorithm computes:
#
#   u     <- y - x
#   delta <- p * x + w - 0.5 d * u^2 + r1 * (beta * y - x) + 0.5 r2 * (beta * y^2 - x^2)
#   z  <- x + beta * lambda * z
#   r0 <- r0 + alpha * delta
#   r1 <- r1 + alpha * delta * z
#   r2 <- r2 + alpha * delta * z^2
#
# starting from the initial state x = ts[1].
#
# It returns a bundle (list) with: 
#
#   - the estimated parameters r0, r1 and r2 of the value function V(x);
#   - the linear coefficients a and b of the one-step look-ahead policy mu(x) = a + b * x.
#

tdLambdaLFA <- function(ts,p=1,d=1,w=0,beta=0.96,c=0.0000001,lambda=0.5){
  
  n <- length(ts)
  
  r0 <- 0
  r1 <- 0
  r2 <- 0
  
  z <- 0
  
  if (length(w) == 1) w <- rnorm(n,0,w)
  
  
  # Value function V(x)
  
  for (i in 1:(n-1)){
    
    x <- ts[i]
    y <- ts[i+1]
    u <- y - x
    
    delta <- p * x + w[i] - d/2 * u^2 + r1 * (beta * y - x) + r2/2 * (beta * y^2 - x^2)
    
    z <- x + beta * lambda * z
    
    alpha <- c / i
    
    r0 <- r0 + alpha * delta
    r1 <- r1 + alpha * delta * z
    r2 <- r2 + alpha * delta * z^2
    
  }
  
 
  # Policy function mu(x)
  
  a <- beta * r1 / (d - beta * r2)
  b <- beta * r2 / (d - beta * r2)
  
  list("r0" = r0, "r1" = r1, "r2" = r2, "a" = a, "b" = b)
  
}


# GTD2(ts,p,d,w,beta,c)
#
# Simple "gradient temporal difference learning, version 2" algorithm with linear function approximation and quadratic polynomial basis (Szepesvari, 2010, p. 26).
#
# For each transition between x and y, the algorithm computes:
#
#   u     <- y - x
#   delta <- p * x + w - 0.5 d * u^2 + r1 * (beta * y - x) + 0.5 r2 * (beta * y^2 - x^2)
#   a0    <- w0
#   a1    <- w1 * x
#   a2    <- w2 * x^2
#   r0    <- r0 + alpha * a0
#   r1    <- r1 + alpha * a1 * (x - beta * y)
#   r2    <- r2 + alpha * a2 * (x^2 - beta * y^2)
#   w0    <- w0 + alpha * (delta - a0)
#   w1    <- w1 + alpha * (delta - a1) * x
#   w2    <- w2 + alpha * (delta - a2) * x^2
#
# starting from the initial state x = ts[1].
#
# It returns a bundle (list) with: 
#
#   - the estimated parameters r0, r1 and r2 of the value function V(x);
#   - the linear coefficients a and b of the one-step look-ahead policy mu(x) = a + b * x.
#

GTD2 <- function(ts,p=1,d=1,w=0,beta=0.96,c=0.0000001){
  
  n <- length(ts)
  
  r0 <- 0
  r1 <- 0
  r2 <- 0
  
  w0 <- 0
  w1 <- 0
  w2 <- 0
  
  if (length(w) == 1) w <- rnorm(n,0,w)
  
  
  # Value function V(x)
  
  for (i in 1:(n-1)){
    
    x <- ts[i]
    y <- ts[i+1]
    u <- y - x
    
    delta <- p * x + w[i] - d/2 * u^2 + r1 * (beta * y - x) + r2/2 * (beta * y^2 - x^2)
    
    a0    <- w0
    a1    <- w1 * x
    a2    <- w2 * x^2
    
    alpha <- c / i
    
    r0 <- r0 + alpha * a0
    r1 <- r1 + alpha * a1 * (x - beta * y)
    r2 <- r2 + alpha * a2 * (x^2 - beta * y^2)
    
    w0    <- w0 + alpha * (delta - a0)
    w1    <- w1 + alpha * (delta - a1) * x
    w2    <- w2 + alpha * (delta - a2) * x^2
  }
  
 
  # Policy function mu(x)
  
  a <- beta * r1 / (d - beta * r2)
  b <- beta * r2 / (d - beta * r2)
  
  list("r0" = r0, "r1" = r1, "r2" = r2, "a" = a, "b" = b)
  
}


# td0NN(x,y,beta,c,h)
#
# Simple temporal difference TD(0) algorithm with nonlinear function approximation (Bertsekas and Tsitsiklis, 1996).
#
# Given a set of states x and a response y, it estimates a cost-to-go function J(x,r,R) using a neural network with ReLU activation function.
#
# Beta is the discount factor beta and c > 0 the step-size sequence coefficient alpha = c/t.
#
# h is the initial value of the neural network weights r and R.
#
# It returns a bundle (list) with: 
#
#   - the estimated weights of the final linear combination of the nodes of neural network (r);
#   - the estimated weights of the hidden layer of the neural network (R).
#

td0NN <- function(x,y,beta=0.96,c=0.0001,h=0.1){
  
  x <- as.matrix(x)
  
  k <- ncol(x)
  n <- length(y)
  
  r <- rep(h,k)
  R <- matrix(h,k,k)
  
  # Cost-to-go function J(x,r,R)
  
  for (i in 1:(n-1)){
    
    X0 <- t(t(x[i,])) %*% x[i,]     # Tensor product for the initial state
    X1 <- t(t(x[i+1,])) %*% x[i+1,] # Tensor product for the next state
    
    y0 <- y[i]
    y1 <- y[i+1]
    
    u <- y1 - y0
    
    g <- u^2     # Immediate cost-to-go
    
    delta <- nodes(X0,R) %*% r - beta * nodes(X1,R) %*% r - g
      
    gamma <- c / i
    
    r <- r - as.vector(gamma) * as.vector(delta) * nodes(X0,R)
    
    R <- R - as.vector(gamma) * as.vector(delta) * dnodes(X0,R) * r
    
  }
  
  
  list("r" = r, "R" = R)
  
}



# ptd0NN(x0,x1,y0,y1,beta,c,h)
#
# TD(0) algorithm with nonlinear function approximation (Bertsekas and Tsitsiklis, 1996) for panel data.
#
# Given the initial state x0, the next state x1, the initial response y0 and the next response y1,
# it estimates a cost-to-go function J(x,r,R) using a neural network with ReLU activation function.
#
# Beta is the discount factor beta and c > 0 the step-size sequence coefficient alpha = c/t.
#
# h is the initial value of the neural network weights r and R.
#
# It returns a bundle (list) with: 
#
#   - the estimated weights of the final linear combination of the nodes of neural network (r);
#   - the estimated weights of the hidden layer of the neural network (R).
#

ptd0NN <- function(x0,x1,y0,y1,beta=0.96,c=0.0001,h=0.1){
  
  x0 <- as.matrix(x0)
  x1 <- as.matrix(x1)
  
  k <- ncol(x0)
  n <- length(y0)
  
  r <- rep(h,k)
  R <- matrix(h,k,k)
  
  # Cost-to-go function J(x,r,R)
  
  for (i in 1:n){
    
    X0 <- t(t(x0[i,])) %*% x0[i,]  # Tensor product for the initial state
    X1 <- t(t(x1[i,])) %*% x1[i,]  # Tensor product for the next state
    
    u <- y1[i] - y0[i]
    
    g <- u^2  # Immediate cost-to-go
    
    delta <- nodes(X0,R) %*% r - beta * nodes(X1,R) %*% r - g
    
    gamma <- c / i
    
    r <- r - as.vector(gamma) * as.vector(delta) * nodes(X0,R)
    
    R <- R - as.vector(gamma) * as.vector(delta) * dnodes(X0,R) * r
    
  }
  
  
  list("r" = r, "R" = R)
  
}


# td0LIN(x,y,beta,c,h)
#
# Simple temporal difference TD(0) algorithm with linear function approximation (Bertsekas and Tsitsiklis, 1996).
#
# Given a set of states x and a response y, it estimates a cost-to-go function J(x,r,R) using a linear architecture.
#
# Beta is the discount factor beta and c > 0 the step-size sequence coefficient alpha = c/t.
#
# h is the initial value of the weights r.
#
# It returns the weights r.
#

td0LIN <- function(x,y,beta=0.96,c=0.0001,h=0.1){
  
  x <- as.matrix(x)
  
  k <- ncol(x)
  n <- length(y)
  
  r <- rep(h,k)
  
  # Cost-to-go function J(x,r,R)
  
  for (i in 1:(n-1)){
    
    x0 <- x[i,] 
    x1 <- x[i+1,]
    
    y0 <- y[i]
    y1 <- y[i+1]
    
    u <- y1 - y0
    
    g <- u^2     # Immediate cost-to-go
    
    delta <- x0 %*% r - beta * x1 %*% r - g
    
    gamma <- c / i
    
    r <- r - as.vector(gamma) * as.vector(delta) * x0
    
  }
  
  return(r)
  
}


# ptd0LIN(x0,x1,y0,y1,beta,c,h)
#
# TD(0) algorithm with linear function approximation (Bertsekas and Tsitsiklis, 1996) for panel data.
#
# Given the initial state x0, the next state x1, the initial response y0 and the next response y1,
# it estimates a cost-to-go function J(x,r,R) using a linear architecture.
#
# Beta is the discount factor beta and c > 0 the step-size sequence coefficient alpha = c/t.
#
# h is the initial value of the weights r.
#
# It returns the weights r.
#

ptd0LIN <- function(x0,x1,y0,y1,beta=0.96,c=0.0001,h=0.1){
  
  x0 <- as.matrix(x0)
  x1 <- as.matrix(x1)
  
  k <- ncol(x0)
  n <- length(y0)
  
  r <- rep(h,k)
  
  # Cost-to-go function J(x,r,R)
  
  for (i in 1:n){
    
    X0 <- x0[i,]  
    X1 <- x1[i,]  
    
    u <- y1[i] - y0[i]
    
    g <- u^2  # Immediate cost-to-go
    
    delta <- X0 %*% r - beta * X1 %*% r - g
    
    gamma <- c / i
    
    r <- r - as.vector(gamma) * as.vector(delta) * X0
    
  }
  
  return(r)
  
}


# LSTD(x,y,beta,h)
#
# Least squares temporal differences (LSTD) algorithm described in (Bertsekas, 2011, p. 367-369).
#
# Given a set of states x and a response y, it estimates the weights r of a linear cost-to-go function J(x,r)
# with immediate cost given by (y1-y0)^2.
#
# h is the initial value of the weights r.
#
# It returns the weights r.
#

LSTD <- function(x,y,beta=0.96,h=0.0001){
  
  x <- as.matrix(x)
  
  s <- ncol(x)
  n <- length(y)
  
  r <- rep(h,s)
  
  
  C <- matrix(0,s,s)
  
  d <- rep(0,s)
  
  # Cost-to-go function J(x,r,R)
  
  for (i in 1:(n-1)){
    
    C <- C + x[i,] %*% t(x[i,] - beta * x[i+1,])
    
    y0 <- y[i]
    y1 <- y[i+1]
    
    u <- y1 - y0
    
    g <- u^2     # Immediate cost-to-go
    
    d <- d + x[i,] * as.vector(g)
    
  }
  
  C <- C / (n-1)
  d <- d / (n-1)
  
  r <- solve(C,d)  
  
  return(r)
  
}


# pLSTD(x0,x1,y0,y1,beta,h)
#
# Least squares temporal differences (LSTD) algorithm described in (Bertsekas, 2011, p. 367-369) for panel data.
#
# Given the initial state x0, the next state x1, the initial response y0 and the next response y1,
# it estimates the weights r of a linear cost-to-go function J(x,r) with immediate cost given by (y1-y0)^2.
#
# h is the initial value of the weights r.
#
# It returns the weights r.
#

pLSTD <- function(x0,x1,y0,y1,beta=0.96,h=0.0001){
  
  x0 <- as.matrix(x0)
  x1 <- as.matrix(x1)
  
  s <- ncol(x0)
  n <- length(y0)
  
  r <- rep(h,s)
  
  
  C <- matrix(0,s,s)
  
  d <- rep(0,s)
  
  # Cost-to-go function J(x,r,R)
  
  for (i in 1:n){
    
    C <- C + x0[i,] %*% t(x0[i,] - beta * x1[i,])
    
    u <- y1[i] - y0[i]
    
    g <- u^2   # Immediate cost-to-go
    
    d <- d + x0[i,] * as.vector(g)
    
  }
  
  C <- C / n
  d <- d / n
  
  r <- solve(C,d)  
  
  return(r)
  
}



### CROSS-VALIDATIONS FOR RL ALGORITHMS ###


# cv.td0LFA(y,p,w,beta,c,k)
#
# Simple k-fold cross validation for a temporal difference TD(0) algorithm with linear function approximation and quadratic polynomial basis.
#
# TIP: calibrate c to give a reasonable min(MAE)).

cv.td0LFA <- function(y,p=1,w=0,beta=0.96,c=0.00001,k=5){

  n <- length(y)
  m <- n - k
  
  nw <- length(w)
  
  y <- as.matrix(y)
 
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  mae_aux <- 1000000000
  
  for (d in 1:20){
  
    for (i in 1:k){
    
      ytrain <- y[i:(m+i-1),1]
      wtrain <- ifelse(nw == 1, w, w[i:(m+i-1),1])
    
      mod <- td0LFA(ytrain,p,d/10,wtrain,beta,c)
    
      pred[i]  <- mod$a + (1 + mod$b) * y[m+i-1]
    
      error[i] <- y[m+i] - pred[i]
    
    }
  
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k

    if (mae <= mae_aux){
        
        pred_aux <- pred
        mae_aux  <- mae
        mse_aux  <- mse
        d_aux    <- d/10
    }
  }
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred_aux, "MAE" = mae_aux, "MSE" = mse_aux, "cost" = d_aux, "cum_error" = cumerr)

}


# cv.td0LFA2(y,p,w,beta,c,k)
#
# Simple k-fold cross validation for a temporal difference TD(0) algorithm with linear function approximation and quadratic polynomial basis, version 2.
#
# TIP: calibrate c to give a reasonable min(MAE).

cv.td0LFA2 <- function(y,p=1,w=0,beta=0.96,c=0.00001,k=5){

  n <- length(y)
  m <- n - k
  
  nw <- length(w)
  
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  mae_aux <- 1000000000
  
  for (d in 1:20){
  
    for (i in 1:k){
    
      ytrain <- y[i:(m+i-1),1]
      wtrain <- ifelse(nw == 1, w, w[i:(m+i-1),1])
    
      mod <- td0LFA2(ytrain,p,d/10,wtrain,beta,c)
    
      pred[i]  <- mod$a + (1 + mod$b) * y[m+i-1]
    
      error[i] <- y[m+i] - pred[i]
    
    }
  
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k

    if (mae <= mae_aux){
        
        pred_aux <- pred
        mae_aux  <- mae
        mse_aux  <- mse
        d_aux    <- d/10
    }
  }
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred_aux, "MAE" = mae_aux, "MSE" = mse_aux, "cost" = d_aux, "cum_error" = cumerr)

}


# cv.td0LFA3(z,y,p,w,beta,c,k)
#
# Simple k-fold cross validation for a temporal difference TD(0) algorithm with:
#
# - Linear value function approximation;
# - Quadratic polynomial basis;
# - Cross-effects between the controlled state y and the non-controllable states z.
#
# TIP: calibrate c to give a reasonable min(MAE)).

cv.td0LFA3 <- function(z,y,p=1,w=0,beta=0.96,c=0.00001,k=5){
  
  n <- length(y)
  m <- n - k
  
  nw <- length(w)
  
  y <- as.matrix(y)
  z <- as.matrix(z)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  mae_aux <- 1000000000
  
  for (d in 1:20){
    
    for (i in 1:k){
      
      ytrain <- y[i:(m+i-1),1]
      ztrain <- z[i:(m+i-1),]
      
      wtrain <- ifelse(nw == 1, w, w[i:(m+i-1),1])
      
      mod <- td0LFA3(ztrain,ytrain,p,d/10,wtrain,beta,c)
      
      pred[i]  <- mod$a + (1 + mod$b) * y[m+i-1] + z[m+i,] %*% mod$q
      
      error[i] <- y[m+i] - pred[i]
      
    }
    
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k
    
    if (mae <= mae_aux){
      
      pred_aux <- pred
      mae_aux  <- mae
      mse_aux  <- mse
      d_aux    <- d/10
    }
  }
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred_aux, "MAE" = mae_aux, "MSE" = mse_aux, "cost" = d_aux, "cum_error" = cumerr)
  
}


# cv.tdLambdaLFA(y,p,w,beta,c,lambda,k)
#
# Simple k-fold cross validation for a temporal difference TD(lambda) algorithm with linear function approximation and quadratic polynomial basis.
#
# TIP: firstly, with lambda = 0, calibrate c to give a reasonable min(MAE); then, calibrate lambda.

cv.tdLambdaLFA <- function(y,p=1,w=0,beta=0.96,c=0.0000001,lambda=0.5,k=5){

  n <- length(y)
  m <- n - k
  
  nw <- length(w)
  
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  mae_aux <- 1000000000
  
  for (d in 1:20){
  
    for (i in 1:k){
    
      ytrain <- y[i:(m+i-1),1]
      wtrain <- ifelse(nw == 1, w, w[i:(m+i-1),1])
    
      mod <- tdLambdaLFA(ytrain,p,d/10,wtrain,beta,c,lambda)
    
      pred[i]  <- mod$a + (1 + mod$b) * y[m+i-1]
    
      error[i] <- y[m+i] - pred[i]
    
    }
  
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k
    
    if (mae <= mae_aux){
        
        pred_aux <- pred
        mae_aux  <- mae
        mse_aux  <- mse
        d_aux    <- d/10
    }
  }
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred_aux, "MAE" = mae_aux, "MSE" = mse_aux, "cost" = d_aux, "cum_error" = cumerr)

}


# cv.GTD2(y,p,w,beta,c,k)
#
# Simple k-fold cross validation for a "gradient temporal difference learning, version 2" algorithm with linear function approximation and quadratic polynomial basis.
#
# TIP: calibrate c to give a reasonable min(MAE).

cv.GTD2 <- function(y,p=1,w=0,beta=0.96,c=0.00001,k=5){

  n <- length(y)
  m <- n - k
  
  nw <- length(w)
  
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  mae_aux <- 1000000000
  
  for (d in 1:20){
  
    for (i in 1:k){
    
      ytrain <- y[i:(m+i-1),1]
      wtrain <- ifelse(nw == 1, w, w[i:(m+i-1),1])
    
      mod <- GTD2(ytrain,p,d/10,wtrain,beta,c)
    
      pred[i]  <- mod$a + (1 + mod$b) * y[m+i-1]
    
      error[i] <- y[m+i] - pred[i]
    
    }
  
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k
  
    if (mae <= mae_aux){
  
        pred_aux <- pred
        mae_aux  <- mae
        mse_aux  <- mse
        d_aux    <- d/10
    }
  }
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred_aux, "MAE" = mae_aux, "MSE" = mse_aux, "cost" = d_aux, "cum_error" = cumerr)

}


# cv.td0ols(z,y,p,w,beta,c,k)
#
# Simple k-fold cross validation for a temporal difference TD(0) algorithm with:
#
# - Linear value function approximation;
# - Quadratic polynomial basis;
# - State ytilde estimated from the attributes z with linear regression (OLS) given the response y.
# - One-step lookahead policy as function of ytilde(t+1) instead of x = y(t) -> z(t+1) is incorporated in the prediction (richer info set).
#
# TIP: calibrate c to give a reasonable min(MAE)).

cv.td0ols <- function(z,y,p=1,w=0,beta=0.96,c=0.00001,k=5){

  n <- length(y)
  m <- n - k
  
  nw <- length(w)
  
  y <- as.matrix(y)
  z <- as.matrix(z)

  pred  <- rep(0,k)
  error <- rep(0,k)
  
  mae_aux <- 1000000000
  
  for (d in 1:20){
  
    for (i in 1:k){
    
      ytrain <- y[i:(m+i-1),1]
      ztrain <- z[i:(m+i-1),]
      ztest  <- z[m+i,]
      
      wtrain <- ifelse(nw == 1, w, w[i:(m+i-1),1])
    
      ols <- lm(ytrain ~ ztrain)
    
      ytilde   <- predict(ols)
      ytilde1  <- ztest %*% coefficients(ols)[-1] + coefficients(ols)[1]
    
      mod <- td0LFA(ytilde,p,d/10,wtrain,beta,c)
    
      pred[i]  <- y[m+i-1] + beta*mod$r1*10/d + beta*mod$r2*10/d * ytilde1
    
      error[i] <- y[m+i] - pred[i]
    
    }
  
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k

    if (mae <= mae_aux){
        
        pred_aux <- pred
        mae_aux  <- mae
        mse_aux  <- mse
        d_aux    <- d/10
    }
  }
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred_aux, "MAE" = mae_aux, "MSE" = mse_aux, "cost" = d_aux, "cum_error" = cumerr)

}


# cv.td0NN(x,y,beta,c,k,h)
#
# Simple k-fold cross validation for a temporal difference TD(0) algorithm with nonlinear function approximation (neural network).
#
# TIP: calibrate c to give a reasonable min(MAE)).

cv.td0NN <- function(x,y,beta=0.96,c=0.00001,k=5,h=0.1){
  
  n <- length(y)
  m <- n - k
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
      
      ytrain <- y[i:(m+i-1),1]
      
      xtrain <- x[i:(m+i-1),]
      xtest  <- x[m+i,]
      
      mod <- td0NN(xtrain,ytrain,beta,c,h)
      
      X1 <- t(t(xtest)) %*% xtest     # Tensor product for the next state
      
      pred[i]  <- nodes(X1,mod$R) %*% mod$r
      
      error[i] <- y[m+i] - pred[i]
      
    }
    
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k
    
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "cum_error" = cumerr)
  
}


# cv.ptd0NN(x0,x1,y0,y1,z1,beta,c,k,h)
#
# Simple k-fold validation for the TD(0) algorithm with nonlinear function approximation (neural network) and panel data.
#
# TIP: calibrate c to give a reasonable min(MAE)).

cv.ptd0NN <- function(x0,x1,y0,y1,z,w,beta=0.96,c=0.00001,k=5,h=0.1){
  
  n <- length(w)
  m <- n - k
  
  x0 <- as.matrix(x0)
  y0 <- as.matrix(y0)
  
  x1 <- as.matrix(x1)
  y1 <- as.matrix(y1)
  
  z <- as.matrix(z)
  w <- as.matrix(w)
  
  mod <- ptd0NN(x0,x1,y0,y1,beta,c,h)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    xtest  <- z[m+i,]
    
    X1 <- t(t(xtest)) %*% xtest     # Tensor product for the next state
    
    pred[i]  <- nodes(X1,mod$R) %*% mod$r
    
    error[i] <- w[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "cum_error" = cumerr)
  
}


# cv.td0LIN(x,y,beta,c,k,h)
#
# Simple k-fold cross validation for a temporal difference TD(0) algorithm with linear function approximation.
#
# TIP: calibrate c to give a reasonable min(MAE)).

cv.td0LIN <- function(x,y,beta=0.96,c=0.00001,k=5,h=0.1){
  
  n <- length(y)
  m <- n - k
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[i:(m+i-1),1]
    
    xtrain <- x[i:(m+i-1),]
    xtest  <- x[m+i,]
    
    r <- td0LIN(xtrain,ytrain,beta,c,h)
    
    pred[i]  <- x[m+i,] %*% r
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "cum_error" = cumerr)
  
}


# cv.ptd0LIN(x0,x1,y0,y1,z1,beta,c,k,h)
#
# Simple k-fold validation for the TD(0) algorithm with linear function approximation and panel data.
#
# TIP: calibrate c to give a reasonable min(MAE)).

cv.ptd0LIN <- function(x0,x1,y0,y1,z,w,beta=0.96,c=0.00001,k=5,h=0.1){
  
  n <- length(w)
  m <- n - k
  
  x0 <- as.matrix(x0)
  y0 <- as.matrix(y0)
  
  x1 <- as.matrix(x1)
  y1 <- as.matrix(y1)
  
  z <- as.matrix(z)
  w <- as.matrix(w)
  
  r <- ptd0LIN(x0,x1,y0,y1,beta,c,h)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    pred[i]  <- z[m+i,] %*% r
    
    error[i] <- w[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "cum_error" = cumerr)
  
}


# cv.LSTD(x,y,beta,k,z)
#
# Simple k-fold cross validation for the least squares temporal differences (LSTD) algorithm.
#

cv.LSTD <- function(x,y,beta=0.96,k=5,z=0.0001){
  
  n <- length(y)
  m <- n - k
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[i:(m+i-1),1]
    
    xtrain <- x[i:(m+i-1),]
    
    r <- LSTD(xtrain,ytrain,beta,z)
    
    pred[i]  <- x[m+i,] %*% r
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "cum_error" = cumerr)
  
}


# cv.pLSTD(x0,x1,y0,y1,z,w,beta,k,h)
#
# Simple k-fold validation for the least squares temporal differences (LSTD) algorithm with panel data.
#

cv.pLSTD <- function(x0,x1,y0,y1,z,w,beta=0.96,k=5,h=0.0001){
  
  n <- length(w)
  m <- n - k
  
  x0 <- as.matrix(x0)
  y0 <- as.matrix(y0)
  
  x1 <- as.matrix(x1)
  y1 <- as.matrix(y1)
  
  z <- as.matrix(z)
  w <- as.matrix(w)
  
  r <- pLSTD(x0,x1,y0,y1,beta,h)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    pred[i]  <- z[m+i,] %*% r
    
    error[i] <- w[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "cum_error" = cumerr)
  
}



### CROSS-VALIDATIONS FOR FILTERS, OLS AND THETA MODEL ###


# cv.hp(y, lambda, k, z)
#
# Simple k-fold cross validation for an univariate time series model with the trend computed by the HP filter.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.hp <- function(y, lambda=1600, k=4, z=2){

    n <- length(y)
    m <- n - k

    y <- as.matrix(y)

    pred  <- rep(0,k)
    error <- rep(0,k)

    for (i in 1:k){
    
        ytrain <- y[i:(m+i-1),1]
    
        trend <- hpfilt(ytrain, lambda)
    
        pred[i]  <- y[m+i-1] + trend[m] - trend[m-1]
    
        error[i] <- y[m+i] - pred[i]
    
    }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  aic <- log(sum(error*error)/k) + 2/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error)) 
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.ma(y, q, k, z)
#
# Simple k-fold cross validation for an univariate time series model with the trend computed by the MA filter with q lags and leads.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.ma <- function(y, q=11, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[i:(m+i-1),1]
    
    trend <- mafilt(ytrain, q)
    
    pred[i]  <- y[m+i-1] + trend[m] - trend[m-1]
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  aic <- log(sum(error*error)/k) + 2/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.med(y, p, k, z)
#
# Simple k-fold cross validation for an univariate time series model with the trend computed by the median filter with p lags and leads.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.med <- function(y, p=11, k=4, z=2){

    n <- length(y)
    m <- n - k

    y <- as.matrix(y)

    pred  <- rep(0,k)
    error <- rep(0,k)

    for (i in 1:k){
    
        ytrain <- y[i:(m+i-1),1]
    
        trend <- medfilt(ytrain, p)
    
        pred[i]  <- y[m+i-1] + trend[m] - trend[m-1]
    
        error[i] <- y[m+i] - pred[i]
    
    }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  aic <- log(sum(error*error)/k) + 2/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.mr(y, theta, k, z)
#
# Simple k-fold cross validation for an univariate time series model with the trend computed by the MR filter.
#
# Returns a list with k OOS predictions and associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.mr <- function(y, theta=40, k=4, z=2){

    n <- length(y)
    m <- n - k

    y <- as.matrix(y)

    pred  <- rep(0,k)
    error <- rep(0,k)

    for (i in 1:k){
    
        ytrain <- y[i:(m+i-1),1]
    
        trend <- mrfilt(ytrain, theta)
    
        pred[i]  <- y[m+i-1] + trend[m] - trend[m-1]
    
        error[i] <- y[m+i] - pred[i]
    
    }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  aic <- log(sum(error*error)/k) + 2/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
 
}


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


# cv.theta(y, k, z)
#
# Simple k-fold cross validation with a sliding window for the Theta model.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.theta <- function(y, k=4, gamma=0.3, z=2){
  
  y <- as.matrix(y)
  
  n <- length(y)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[i:(m+i-1),1]
    
    pred[i]  <- theta(ytrain, h=1, gamma)
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  aic <- log(sum(error*error)/k) + 2/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.AR1(y, k, z)
#
# Simple k-fold cross validation with a sliding window for the AR(1) model.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.AR1 <- function(y, k=4, z=2){
  
  y <- as.matrix(y)
  
  n <- length(y)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[i:(m+i-1),1]
    
    pred[i]  <- AR1(ytrain)
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  aic <- log(sum(error*error)/k) + 2/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv.AR1MAX(x, y, k, z)
#
# Simple k-fold cross validation with a sliding window for the AR(1) model with covariates x.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv.AR1MAX <- function(x, y, k=4, z=2){
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  n <- length(y)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[i:(m+i-1),1]
    xtrain <- x[i:(m+i-1),]
    
    pred[i]  <- AR1MAX(xtrain, ytrain)
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x) + 1
  
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




### CROSS-VALIDATIONS FOR RL ALGORITHMS WITH GROWING TRAINING WINDOW ###


# cv1.td0LFA(y,p,w,beta,c,k)
#
# Simple k-fold cross validation for a temporal difference TD(0) algorithm with linear function approximation,
# quadratic polynomial basis and growing training window.
#
# TIP: calibrate c to give a reasonable min(MAE)).

cv1.td0LFA <- function(y,p=1,w=0,beta=0.96,c=0.00001,k=5){
  
  n <- length(y)
  m <- n - k
  
  nw <- length(w)
  
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  mae_aux <- 1000000000
  
  for (d in 1:20){
    
    for (i in 1:k){
      
      ytrain <- y[1:(m+i-1),1]
      wtrain <- ifelse(nw == 1, w, w[1:(m+i-1),1])
      
      mod <- td0LFA(ytrain,p,d/10,wtrain,beta,c)
      
      pred[i]  <- mod$a + (1 + mod$b) * y[m+i-1]
      
      error[i] <- y[m+i] - pred[i]
      
    }
    
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k
    
    if (mae <= mae_aux){
      
      pred_aux <- pred
      mae_aux  <- mae
      mse_aux  <- mse
      d_aux    <- d/10
    }
  }
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred_aux, "MAE" = mae_aux, "MSE" = mse_aux, "cost" = d_aux, "cum_error" = cumerr)
  
}


# cv1.td0LFA2(y,p,w,beta,c,k)
#
# Simple k-fold cross validation for a temporal difference TD(0) algorithm with linear function approximation,
# quadratic polynomial basis, version 2, and growing training window.
#
# TIP: calibrate c to give a reasonable min(MAE).

cv1.td0LFA2 <- function(y,p=1,w=0,beta=0.96,c=0.00001,k=5){
  
  n <- length(y)
  m <- n - k
  
  nw <- length(w)
  
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  mae_aux <- 1000000000
  
  for (d in 1:20){
    
    for (i in 1:k){
      
      ytrain <- y[1:(m+i-1),1]
      wtrain <- ifelse(nw == 1, w, w[1:(m+i-1),1])
      
      mod <- td0LFA2(ytrain,p,d/10,wtrain,beta,c)
      
      pred[i]  <- mod$a + (1 + mod$b) * y[m+i-1]
      
      error[i] <- y[m+i] - pred[i]
      
    }
    
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k
    
    if (mae <= mae_aux){
      
      pred_aux <- pred
      mae_aux  <- mae
      mse_aux  <- mse
      d_aux    <- d/10
    }
  }
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred_aux, "MAE" = mae_aux, "MSE" = mse_aux, "cost" = d_aux, "cum_error" = cumerr)
  
}


# cv1.tdLambdaLFA(y,p,w,beta,c,lambda,k)
#
# Simple k-fold cross validation for a temporal difference TD(lambda) algorithm with linear function approximation,
# quadratic polynomial basis and growing training window.
#
# TIP: firstly, with lambda = 0, calibrate c to give a reasonable min(MAE); then, calibrate lambda.

cv1.tdLambdaLFA <- function(y,p=1,w=0,beta=0.96,c=0.0000001,lambda=0.5,k=5){
  
  n <- length(y)
  m <- n - k
  
  nw <- length(w)
  
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  mae_aux <- 1000000000
  
  for (d in 1:20){
    
    for (i in 1:k){
      
      ytrain <- y[1:(m+i-1),1]
      wtrain <- ifelse(nw == 1, w, w[1:(m+i-1),1])
      
      mod <- tdLambdaLFA(ytrain,p,d/10,wtrain,beta,c,lambda)
      
      pred[i]  <- mod$a + (1 + mod$b) * y[m+i-1]
      
      error[i] <- y[m+i] - pred[i]
      
    }
    
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k
    
    if (mae <= mae_aux){
      
      pred_aux <- pred
      mae_aux  <- mae
      mse_aux  <- mse
      d_aux    <- d/10
    }
  }
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred_aux, "MAE" = mae_aux, "MSE" = mse_aux, "cost" = d_aux, "cum_error" = cumerr)
  
}


# cv1.GTD2(y,p,w,beta,c,k)
#
# Simple k-fold cross validation for a "gradient temporal difference learning, version 2" algorithm with linear 
# function approximation, quadratic polynomial basis and growing training window.
#
# TIP: calibrate c to give a reasonable min(MAE).

cv1.GTD2 <- function(y,p=1,w=0,beta=0.96,c=0.00001,k=5){
  
  n <- length(y)
  m <- n - k
  
  nw <- length(w)
  
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  mae_aux <- 1000000000
  
  for (d in 1:20){
    
    for (i in 1:k){
      
      ytrain <- y[1:(m+i-1),1]
      wtrain <- ifelse(nw == 1, w, w[1:(m+i-1),1])
      
      mod <- GTD2(ytrain,p,d/10,wtrain,beta,c)
      
      pred[i]  <- mod$a + (1 + mod$b) * y[m+i-1]
      
      error[i] <- y[m+i] - pred[i]
      
    }
    
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k
    
    if (mae <= mae_aux){
      
      pred_aux <- pred
      mae_aux  <- mae
      mse_aux  <- mse
      d_aux    <- d/10
    }
  }
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred_aux, "MAE" = mae_aux, "MSE" = mse_aux, "cost" = d_aux, "cum_error" = cumerr)
  
}


# cv1.td0ols(z,y,p,w,beta,c,k)
#
# Simple k-fold cross validation for a temporal difference TD(0) algorithm with:
#
# - Linear value function approximation;
# - Quadratic polynomial basis;
# - State ytilde estimated from the attributes z with linear regression (OLS) given the response y;
# - One-step lookahead policy as function of ytilde(t+1) instead of x = y(t) -> z(t+1) is incorporated in the prediction (richer info set);
# - Growing training window.
#
# TIP: calibrate c to give a reasonable min(MAE)).

cv1.td0ols <- function(z,y,p=1,w=0,beta=0.96,c=0.00001,k=5){
  
  n <- length(y)
  m <- n - k
  
  nw <- length(w)
  
  y <- as.matrix(y)
  z <- as.matrix(z)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  mae_aux <- 1000000000
  
  for (d in 1:20){
    
    for (i in 1:k){
      
      ytrain <- y[1:(m+i-1),1]
      ztrain <- z[1:(m+i-1),]
      ztest  <- z[m+i,]
      
      wtrain <- ifelse(nw == 1, w, w[1:(m+i-1),1])
      
      ols <- lm(ytrain ~ ztrain)
      
      ytilde   <- predict(ols)
      ytilde1  <- ztest %*% coefficients(ols)[-1] + coefficients(ols)[1]
      
      mod <- td0LFA(ytilde,p,d/10,wtrain,beta,c)
      
      pred[i]  <- y[m+i-1] + beta*mod$r1*10/d + beta*mod$r2*10/d * ytilde1
      
      error[i] <- y[m+i] - pred[i]
      
    }
    
    mae <- sum(abs(error))/k
    mse <- sum(error*error)/k
    
    if (mae <= mae_aux){
      
      pred_aux <- pred
      mae_aux  <- mae
      mse_aux  <- mse
      d_aux    <- d/10
    }
  }
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred_aux, "MAE" = mae_aux, "MSE" = mse_aux, "cost" = d_aux, "cum_error" = cumerr)
  
}



### CROSS-VALIDATIONS FOR FILTERS, OLS AND THETA MODEL WITH GROWING TRAINING WINDOW ###


# cv1.hp(y, lambda, k, z)
#
# Simple k-fold cross validation for an univariate time series model with the trend computed by the HP filter
# and growing training window.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.hp <- function(y, lambda=1600, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    
    trend <- hpfilt(ytrain, lambda)
    
    pred[i]  <- y[m+i-1] + trend[m] - trend[m-1]
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  aic <- log(sum(error*error)/k) + 2/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.ma(y, q, k, z)
#
# Simple k-fold cross validation for an univariate time series model with the trend computed by the MA filter
# with q lags and leads and growing training window.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.ma <- function(y, q=11, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    
    trend <- mafilt(ytrain, q)
    
    pred[i]  <- y[m+i-1] + trend[m] - trend[m-1]
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  aic <- log(sum(error*error)/k) + 2/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.med(y, p, k, z)
#
# Simple k-fold cross validation for an univariate time series model with the trend computed by the median filter
# with p lags and leads and growing training window.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.med <- function(y, p=11, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    
    trend <- medfilt(ytrain, p)
    
    pred[i]  <- y[m+i-1] + trend[m] - trend[m-1]
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  aic <- log(sum(error*error)/k) + 2/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.mr(y, theta, k, z)
#
# Simple k-fold cross validation for an univariate time series model with the trend computed by the MR filter
# and growing training window.
#
# Returns a list with k OOS predictions and associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.mr <- function(y, theta=40, k=4, z=2){
  
  n <- length(y)
  m <- n - k
  
  y <- as.matrix(y)
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    
    trend <- mrfilt(ytrain, theta)
    
    pred[i]  <- y[m+i-1] + trend[m] - trend[m-1]
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  aic <- log(sum(error*error)/k) + 2/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.ols(x, y, k, z)
#
# Simple k-fold cross validation for a linear regression of time series and growing training window.
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


# cv1.theta(y, k, z)
#
# Simple k-fold cross validation with a growing training window for the Theta model.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.theta <- function(y, k=4, gamma=0.3, z=2){
  
  y <- as.matrix(y)
  
  n <- length(y)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    
    pred[i]  <- theta(ytrain, h=1, gamma)
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  aic <- log(sum(error*error)/k) + 2/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.AR1(y, k, z)
#
# Simple k-fold cross validation with a growing training window for the AR(1) model.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.AR1 <- function(y, k=4, z=2){
  
  y <- as.matrix(y)
  
  n <- length(y)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    
    pred[i]  <- AR1(ytrain)
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  aic <- log(sum(error*error)/k) + 2/k
  
  sd <- sd(error)
  a <- abs(error)
  b <- ifelse(a<z*sd,a,0)
  c <- ifelse(a<z*sd,1,0)
  mae2 <- sum(b)/sum(c)
  mse2 <- sum(b*b)/sum(c)
  
  cumerr <- cumsum(abs(error))
  
  list("prediction" = pred, "MAE" = mae, "MSE" = mse, "AIC" = aic, "MAE2" = mae2, "MSE2" = mse2, "cum_error" = cumerr)
  
}


# cv1.AR1MAX(x, y, k, z)
#
# Simple k-fold cross validation with a growing window for the AR(1) model with covariates x.
#
# Returns a list with k OOS predictions and the associated errors (MAE and MSE) and AIC (Greene approach).
#

cv1.AR1MAX <- function(x, y, k=4, z=2){
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  n <- length(y)
  m <- n - k
  
  pred  <- rep(0,k)
  error <- rep(0,k)
  
  for (i in 1:k){
    
    ytrain <- y[1:(m+i-1),1]
    xtrain <- x[1:(m+i-1),]
    
    pred[i]  <- AR1MAX(xtrain, ytrain)
    
    error[i] <- y[m+i] - pred[i]
    
  }
  
  mae <- sum(abs(error))/k
  mse <- sum(error*error)/k
  
  q <- ncol(x) + 1
  
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



### IN-SAMPLE ERRORS ###


# err.hp(y, lambda)
#
# Computes the in-sample mean absolute and squared fitting error (= cycle) of the HP filter.
#

err.hp <- function(y, lambda=1600){
  
  n <- length(y)  
  
  error <- y - hpfilt(y, lambda)
  
  mae <- sum(abs(error))/n
  mse <- sum(error*error)/n
  sd  <- sd(error)
  
  list("MAE" = mae, "MSE" = mse, "sd" = sd)
  
}


# err.mr(y, theta)
#
# Computes the in-sample mean absolute and squared fitting error (= cycle) of the MR filter.
#

err.mr <- function(y, theta=40){
  
  n <- length(y)  
  
  error <- y - mrfilt(y, theta)
  
  mae <- sum(abs(error))/n
  mse <- sum(error*error)/n
  sd  <- sd(error)
  
  list("MAE" = mae, "MSE" = mse, "sd" = sd)
  
}


