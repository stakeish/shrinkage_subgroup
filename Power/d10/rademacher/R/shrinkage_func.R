##Updated on January 30th
shrinkageMLEinit <- function(ninits, gamma.dim){
  #alpha <- matrix(runif(alpha.dim*ninits, min = -10, max = 10), nrow=alpha.dim)
  beta <- matrix(runif(ninits, min = -5, max=5), nrow = 1)
  lambda <- matrix(runif(ninits, min = 0, max=5), nrow=1)
  #sigma <- matrix(runif(ninits, min = 0.1, max=10), nrow=1)
  gamma <- matrix(runif(ninits*gamma.dim, min=-5, max=5), nrow=gamma.dim)
  
  list(beta=beta, lambda=lambda, gamma=gamma)
}

rlognormal <- function(n, x, d, z, alpha, beta, lambda, gamma, sigma){
  membership <- rbinom(n, 1, plogis(z %*% gamma))
  subgroup_mean <- ifelse(membership==1, x %*% alpha + d * (beta + lambda), x %*% alpha + d * beta)
  rnorm(n, subgroup_mean, sigma)
}