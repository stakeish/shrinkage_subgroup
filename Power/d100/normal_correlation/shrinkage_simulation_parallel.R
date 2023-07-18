rm(list=ls())

#setwd("power_d_100_lambda_1_correlation")

library(parallel)
library(abind)
library(Matrix)
library(quadprog)
library(nloptr)
library(lbfgs)
library(MASS)

load("coefficients.RData")
load("summary_null.RData")
source("R/shrinkage_func.R")
source("R/shrinkage_estimation.R")
Rcpp::sourceCpp("src/cppshrinkage.cpp")
intercept <- coefficients[1]
coef <- coefficients[2]
nreps <- 5000
gamma.dim <- 100
sample_size <- c(100, 250, 500, 750, 1000)

adj_crt_vec_equal <- summary_matrix[, 3]
adj_crt_vec_SLR <- summary_matrix[, 4]
critical_value <- qnorm(0.95)^2

alpha <- c(1, 2)
beta <- 1
lambda <- 1
gamma <-  rep(1, gamma.dim)
sigma <- 1

set.seed(1234)
Sigma_mat_root <- matrix(runif((gamma.dim - 1)^2, min = -1, max = 1), ncol = gamma.dim - 1)
Sigma_mat_square <- Sigma_mat_root^2
Sigma_mat_sum <- apply(Sigma_mat_square, 1, sum)
Sigma_mat_root <- Sigma_mat_root / sqrt(Sigma_mat_sum)
Sigma_mat <- Sigma_mat_root %*% t(Sigma_mat_root)

alpha.dim <- length(alpha)

summary_matrix <- matrix(double(length(sample_size)*4), nrow = length(sample_size), ncol = 4)
colnames(summary_matrix) <- c("power equal", "power SLR", "adj. power equal", "adj. power SLR")

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

for(i in 5:length(sample_size)) {
  n <- sample_size[i]
  p <- intercept + coef * n^(7/8) * sqrt(log(gamma.dim))
  adj_crt_equal <- adj_crt_vec_equal[i]
  adj_crt_SLR <- adj_crt_vec_SLR[i]
  
  set.seed(20)
  
  intercept_set <- array(rep(1, n*nreps), dim = c(n, 1, nreps))
  x_set <- array(rnorm(n*nreps*(alpha.dim - 1)), dim = c(n, alpha.dim - 1, nreps))
  x_set <- abind(intercept_set, x_set, along = 2)
  d_set <- array(rbinom(n*nreps, 1, 0.5), dim = c(n, 1, nreps))
  z_set <- lapply(1:nreps, function(j, Sigma_mat, gamma.dim, n){
    norm_part <- mvrnorm(n, mu = rep(0, gamma.dim - 1), Sigma = Sigma_mat)
    cbind(rep(1, n), norm_part)
  }, Sigma_mat = Sigma_mat, gamma.dim = gamma.dim, n = n)
  y_set <- lapply(1:nreps, function(j, x_set, d_set, z_set, alpha, beta, lambda, gamma, sigma){
    x <- x_set[, , j]
    d <- d_set[, , j]
    z <- z_set[[j]]
    n <- nrow(x)
    rlognormal(n, x, d, z, alpha, beta, lambda, gamma, sigma)
  }, x_set = x_set, d_set = d_set, z_set = z_set, alpha = alpha, beta = beta, lambda = lambda, gamma = gamma, sigma = sigma)
  
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
  rfreq_equal_adj <- mean(LRT_equal > adj_crt_equal)
  save(ll_null_all, file = "ll_null_all.RData")
  save(out_equal_all, file = "out_equal_all.RData")
  
  rm(list=c("intercept_set", "out_equal_all"))
  
  ll_alt <- parSapply(cl, 1:nreps, function(j, x_set, d_set, z_set, y_set, theta_initial_set, p) {
    x <- x_set[, , j]
    d <- d_set[, , j]
    z <- z_set[[j]]
    y <- y_set[[j]]
    theta_initial <- theta_initial_set[[j]]
    out <- shrinkageMLE(x, d, z, y, p, ninits = 10, theta_initial)
    gamma <- out$par$gamma
    out$ploglik + p*norm(as.matrix(gamma))
  }, x_set = x_set, d_set = d_set, z_set = z_set, y_set = y_set, theta_initial_set = theta_initial_set, p = p)
  SLR <- 2*(ll_alt - ll_null_all)
  rfreq <- mean(SLR > critical_value)
  rfreq_adj <- mean(SLR > adj_crt_SLR)
  
  summary_matrix[i, ] <- c(rfreq_equal, rfreq, rfreq_equal_adj, rfreq_adj)
  save(summary_matrix, file = "summary.RData")
  print(summary_matrix)
  
}

summary_matrix

stopCluster(cl)