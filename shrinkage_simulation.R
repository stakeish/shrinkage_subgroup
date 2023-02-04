rm(list=ls())

Rcpp::sourceCpp("src/cppshrinkage.cpp")
source("R/shrinkage_func.R")
source("R/shrinkage_estimation.R")

library(Matrix)
library(quadprog)
library(nloptr)
library(lbfgs)

##Compute a quantile of half-chi-square distribution
critical_value <- qnorm(0.95)^2

##Simulation setting
set.seed(10)

nrep <- 2000
n <- 500
p.all <- seq(from = 9.25, to = 13, by = 0.25)
p.length <- length(p.all)
rejection <- matrix(double(nrep * p.length), nrow = nrep, ncol = p.length)

alpha <- c(1, 2)
beta <- 1
lambda <- 0
gamma <- rep(0.5, 10)
sigma <- 1

alpha.dim <- length(alpha)
gamma.dim <- length(gamma)

#x1big <- rep(1, n*nrep)
##x2big <- rnorm(n*nrep, -1, 1)
#xbig <- cbind(x1big, x2big)
#dbig <- rbinom(n*nrep, 1, 0.5)
#zbig <- matrix(rnorm(n*nrep*(gamma.dim-1), mean = 0, sd = 1), ncol=gamma.dim-1)
#zbig <- cbind(rep(1, n*nrep), zbig)
#ybig <- double(n*nrep)
#for(iter in 1:nrep){
#  x <- xbig[((iter-1)*n + 1):(iter*n), ]
#  d <- dbig[((iter-1)*n + 1):(iter*n)]
#  z <- zbig[((iter-1)*n + 1):(iter*n), ]
#  ybig[((iter-1)*n + 1):(iter*n)] <- rlognormal(n, x, d, z, alpha, beta, lambda, gamma, sigma)
#}

b.time <- proc.time()
for(rep in 1:nrep){
  x1 <- rep(1, n)
  x2 <- rnorm(n, -1, 1)
  x <- cbind(x1, x2)
  d <- rbinom(n, 1, 0.5)
  z <- matrix(rnorm(n*(gamma.dim-1), mean = 0, sd = 1), nrow=n, ncol= gamma.dim - 1)
  z <- cbind(rep(1, n), z)
  y <- rlognormal(n, x, d, z, alpha, beta, lambda, gamma, sigma)
  ### We compute MLE under the null hypothesis.
  out.null <- homogeneousMLE(x, d, y)
  ll.null <- out.null$ll
  
  ##We compute MLE under mixture model with equal weight
  out.equal <- heterogeneousMLE(x, d, y, 10)
  ll.equal <- out.equal$ll
  alpha.initial <- out.equal$par$alpha
  beta.initial <- out.equal$par$beta
  lambda.initial <- out.equal$par$lambda
  sigma.initial <- out.equal$par$sigma
  theta.initial <- c(alpha.initial, beta.initial, lambda.initial, sigma.initial, rep(0, gamma.dim))
  
  ### We compute MLE for logistic normal mixture model
  ### Calculate the initial value for parameter
  #y.nontreat <- y[d == 0]
  #x.nontreat <- x[d == 0, ]
  #alpha.initial <- solve(t(x.nontreat) %*% x.nontreat) %*% t(x.nontreat) %*% y.nontreat
  #sigma.initial <- sqrt(mean((y.nontreat - x.nontreat %*% alpha.initial)^2))
  
  for(piter in 1:p.length){
    p <- p.all[piter]
    out.alt <- shrinkageMLE(x, d, z, y, p, ninits=10, theta.initial)
    gamma.alt <- out.alt$par$gamma
    ll.alt <- out.alt$ploglik + p*norm(as.matrix(gamma.alt))
    SLR <- 2*(ll.alt - ll.null)
    rejection[rep, piter] <- SLR > critical_value
  }
  if(floor(rep/10)==rep/10){
    print(rep)
  }
}
e.time <- proc.time()
colMeans(rejection)