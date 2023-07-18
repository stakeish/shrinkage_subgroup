rm(list=ls())

library(parallel)
library(abind)
library(Matrix)
library(quadprog)
library(nloptr)
library(lbfgs)

load("coefficients.RData")
source("R/shrinkage_func.R")
source("R/shrinkage_estimation.R")
Rcpp::sourceCpp("src/cppshrinkage.cpp")
intercept <- coefficients[1]
coef <- coefficients[2]
nreps <- 2000
n <- 750
d <- 25
p <- intercept + coef * n^(7/8) * sqrt(log(d))
p_all <- seq(from = 11.3, to = 12, by = 0.1)

critical_value <- qnorm(0.95)^2

alpha <- c(1, 2)
beta <- 1
lambda <- 0
gamma <-  rep(1, 25)
sigma <- 1

alpha.dim <- length(alpha)
gamma.dim <- length(gamma)

rfreq_all <- double(length(p_all))

# Set up a cluster and load cpp files
nworker <- min(nreps,detectCores()-1)
cl <- makePSOCKcluster(rep('localhost',nworker),master='localhost')
clusterSetRNGStream(cl, 123456)
# workerごとに異なる乱数を発生させる場合に使う
clusterEvalQ(cl, {
  library(Matrix)
  library(quadprog)
  library(nloptr)
  library(lbfgs)
  source("R/shrinkage_func.R")
  source("R/shrinkage_estimation.R")
  Rcpp::sourceCpp("src/cppshrinkage.cpp")
})
# end of setting up a cluster

b.time <- proc.time()
set.seed(10)

intercept_set <- array(rep(1, n*nreps), dim = c(n, 1, nreps))
x_set <- array(rnorm(n*nreps*(alpha.dim - 1)), dim = c(n, alpha.dim - 1, nreps))
x_set <- abind(intercept_set, x_set, along = 2)
d_set <- array(rbinom(n*nreps, 1, 0.5), dim = c(n, 1, nreps))
z_set <- array(rnorm(n*nreps*(gamma.dim - 1)), dim = c(n, gamma.dim - 1, nreps))
z_set <- abind(intercept_set, z_set, along = 2)
y_set <- lapply(1:nreps, function(j, x_set, d_set, z_set, alpha, beta, lambda, gamma, sigma){
  x <- x_set[, , j]
  d <- d_set[, , j]
  z <- z_set[, , j]
  n <- nrow(x)
  rlognormal(n, x, d, z, alpha, beta, lambda, gamma, sigma)
}, x_set = x_set, d_set = d_set, z_set = z_set, alpha = alpha, beta = beta, lambda = lambda, gamma = gamma, sigma = sigma)

for(i in 1:(length(p_all) + 1)){
  
  if(i == 1) {
    ll_null_all <- parSapply(cl, 1:nreps, function(j, x_set, d_set, y_set){
      x <- x_set[, , j]
      d <- d_set[, , j]
      y <- y_set[[j]]
      n <- nrow(x)
      out <- homogeneousMLE(x, d, y)
      out$ll
    }, x_set = x_set, d_set = d_set, y_set = y_set)
    
    out_equal_all <- parLapply(cl, 1:nreps, function(j, x_set, d_set, y_set, gamma.dim) {
      x <- x_set[, , j]
      d <- d_set[, , j]
      y <- y_set[[j]]
      out <- heterogeneousMLE(x, d, y, 10)
      ll <- out$ll
      alpha <- out$par$alpha
      beta <- out$par$beta
      lambda <- out$par$lambda
      sigma <- out$par$sigma
      theta <- c(alpha, beta, lambda, sigma, rep(0, gamma.dim))
      list(ll = ll, theta_initial = theta)
    }, x_set = x_set, d_set = d_set, y_set = y_set, gamma.dim = gamma.dim)
    
    ll_equal_alt <- sapply(out_equal_all, "[[", "ll")
    theta_initial_set <- lapply(out_equal_all, "[[", "theta_initial")
    LRT_equal <- 2*(ll_equal_alt - ll_null_all)
    rfreq_equal <- mean(LRT_equal > critical_value)
  } else {
    p <- p_all[i - 1]
    ll_alt <- parSapply(cl, 1:nreps, function(j, x_set, d_set, z_set, y_set, theta_initial_set, p) {
      x <- x_set[, , j]
      d <- d_set[, , j]
      z <- z_set[, , j]
      y <- y_set[[j]]
      theta_initial <- theta_initial_set[[j]]
      out <- shrinkageMLE(x, d, z, y, p, ninits = 10, theta_initial)
      gamma <- out$par$gamma
      out$ploglik + p*norm(as.matrix(gamma))
    }, x_set = x_set, d_set = d_set, z_set = z_set, y_set = y_set, theta_initial_set = theta_initial_set, p = p)
    SLR <- 2*(ll_alt - ll_null_all)
    rfreq <- mean(SLR > critical_value)
    rfreq_all[i - 1] <- rfreq
  }
  
  if(i == 1){
    print(rfreq_equal)
  } else {
    print(rfreq_all)
  }
}
e.time <- proc.time()

print(rfreq_equal)
print(rfreq_all)

stopCluster(cl)