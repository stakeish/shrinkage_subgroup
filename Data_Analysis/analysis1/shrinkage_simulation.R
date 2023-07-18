rm(list=ls())

library(parallel)
library(abind)
library(Matrix)
library(quadprog)
library(nloptr)
library(lbfgs)
library(speff2trial)
load("coefficients.RData")
source("R/shrinkage_func.R")
source("R/shrinkage_estimation.R")
Rcpp::sourceCpp("src/cppshrinkage.cpp")

data("ACTG175")
intercept <- coefficients[1]
coef <- coefficients[2]

data <- ACTG175
x <- cbind(rep(1, nrow(data)), data$age, data$wtkg, data$karnof, data$cd40, data$cd80, data$hemo, data$homo, data$drugs, data$race, data$gender, data$str2, data$symptom)
d <- as.numeric(data$arms > 0)
y <- data$cd420
z <- cbind(rep(1, nrow(data)), scale(data$age), scale(data$wtkg), scale(data$karnof), scale(data$cd40), scale(data$cd80), scale(data$hemo), scale(data$homo), scale(data$drugs), scale(data$race), scale(data$gender), scale(data$str2), scale(data$symptom))

n <- nrow(data)
gamma.dim <- ncol(z)
p <- intercept + coef * n^(7/8) * sqrt(log(gamma.dim))

critical_value <- qnorm(0.95)^2

ll_null <- homogeneousMLE(x, d, y)$ll

out_equal <- heterogeneousMLE(x, d, y, 100)
ll_equal <- out_equal$ll
alpha_equal <- out_equal$par$alpha
beta_equal <- out_equal$par$beta
lambda_equal <- out_equal$par$lambda
sigma_equal <- out_equal$par$sigma
theta_equal <- c(alpha_equal, beta_equal, lambda_equal, sigma_equal, rep(0, gamma.dim))

out <- shrinkageMLE(x, d, z, y, p, ninits = 100, theta_equal)
gamma <- out$par$gamma
beta <- out$par$beta
lambda <- out$par$lambda
ll_alt <- out$ploglik + p*norm(as.matrix(gamma))

SLR <- 2*(ll_alt - ll_null)
pval <- 1-pnorm(sqrt(pmax(0,SLR)))

save(pval, file = "pval.RData")
#out_nonpenalty <- shrinkageMLE(x, d, z, y, 0, ninits = 50, theta_equal)
#gamma_nonpenalty <- out_nonpenalty$par$gamma