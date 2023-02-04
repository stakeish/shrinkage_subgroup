homogeneousMLE <- function(x, d, y){
  x.dim <- ncol(x)
  X <- cbind(x, d)
  estimator <- solve(t(X) %*% X) %*% t(X) %*% y
  alpha.hat <- estimator[1:x.dim]
  beta.hat <- estimator[x.dim+1]
  sigma.hat <- sqrt(mean((y - X %*% estimator)^2))
  ll <- sum(log(dnorm(y - X %*% estimator, mean=0, sd=sigma.hat)))
  par <- list(alpha = alpha.hat, beta = beta.hat, sigma = sigma.hat)
  
  list(ll = ll, par=par)
}

heterogeneousMLE <- function(x, d, y, ninits){
  n <- length(y)
  alpha.dim <- ncol(x)
  x.null <- x[d == 0, ]
  y.null <- y[d == 0]
  alpha.null <- solve(t(x.null) %*% x.null) %*% t(x.null) %*% y.null
  sigma.null <- sqrt(mean((y.null - x.null %*% alpha.null)^2))
  
  alpha.set <- matrix(0, nrow = alpha.dim, ncol = ninits)
  beta.set <- double(ninits)
  lambda.set <- double(ninits)
  sigma.set <- double(ninits)
  ll.set <- double(ninits)
  
  for(i in 1:ninits){
    
    ##Initialize EM iteration
    pi <- 1/2
    alpha <- alpha.null
    beta <- runif(1, min = -5, max=5)
    lambda <- runif(1, min=0, max=5)
    sigma <- sigma.null
    diff <- 1
    oldll <- -Inf
    for(j in 1:500){
      
      mean.1 <- y - x %*% alpha - d*(beta + lambda)
      mean.2 <- y - x %*% alpha - d*beta
      w1 <- pi*dnorm(mean.1, 0, sigma)
      w2 <- (1-pi)*dnorm(mean.2, 0, sigma)
      w.sum <- w1+w2
      w <- w1 / w.sum
      
      ll <- sum(log(w.sum))
      
      diff <- ll - oldll 
      oldll <- ll
      #print(diff)
      if(diff < 1e-08){
        break
      }
      
      
      coef_vec <- c(t(x) %*% y, sum(y * d), sum(w * y * d))
      data_mat_w <- cbind(x, d, w * d)
      coef_mat <- t(data_mat_w) %*% data_mat_w
      coef_mat[alpha.dim + 2, alpha.dim + 2] <- sum(w * d)
      const_mat <- matrix(0, nrow = alpha.dim + 2, ncol = alpha.dim + 2)
      const_mat[alpha.dim + 2, alpha.dim + 2] <- 1
      result <- tryCatch({solve.QP(coef_mat, coef_vec, const_mat)}
                         , error = function(err) {coef_mat <- nearPD(coef_mat)$mat
                         return(solve.QP(coef_mat, coef_vec, const_mat))})
      
      alpha <- result$solution[1:alpha.dim]
      beta <- result$solution[alpha.dim + 1]
      lambda <- result$solution[alpha.dim + 2]
      sigma <- sqrt(mean(w * (y - x %*% alpha - d * (beta + lambda))^2) + mean((1 - w) * (y - x %*% alpha - d * beta)^2))
    }
    
    alpha.set[, i] <- alpha
    beta.set[i] <- beta
    lambda.set[i] <- lambda
    sigma.set[i] <- sigma
    ll.set[i] <- ll
  }
  
  index <- which.max(ll.set)
  alpha <- alpha.set[, index]
  beta <- beta.set[index]
  lambda <- lambda.set[index]
  sigma <- sigma.set[index]
  parm <- list(alpha=alpha, beta=beta, lambda=lambda, sigma=sigma)
  ll <- ll.set[index]
  
  list(ll=ll, parm=parm)
}

shrinkageMLE <- function(x, d, z, y, p, ninits=5, theta.initial, epsilon = 1e-08, maxit = 500, 
                         epsilon.short = 1e-02, maxit.short = 100){
  ##Updated on January 30th
  n <- length(y)
  alpha.dim <- ncol(x)
  gamma.dim <- ncol(z)
  ninits.short <- ninits*gamma.dim
  alpha.initial <- theta.initial[1:alpha.dim]
  sigma.initial <- theta.initial[alpha.dim+3]
  
  tmp <- shrinkageMLEinit(ninits=ninits.short, gamma.dim=gamma.dim)
  alpha.tmp <- matrix(rep(alpha.initial, ninits.short), nrow = alpha.dim, ncol = ninits.short)
  sigma.tmp <- rep(sigma.initial, nrow = 1, ncol = ninits.short)
  
  theta0 <- rbind(alpha.tmp, tmp$beta, tmp$lambda, sigma.tmp, tmp$gamma)
  theta0[, 1] <- theta.initial
  ploglikset <- double(ninits.short)
  
  w.mat <- matrix(nrow=n, ncol=2)
  
  for(ii in 1:ninits.short){
    alpha <- theta0[1:alpha.dim, ii]
    beta <- theta0[alpha.dim+1, ii]
    lambda <- theta0[alpha.dim+2, ii]
    sigma <- theta0[alpha.dim+3, ii]
    gamma <- theta0[(alpha.dim+4):(alpha.dim+gamma.dim+3), ii]
    diff <- 1.0
    oldpenloglik <- -Inf
    for(jj in 1:maxit.short){
      
      ##E-step
      mean.2 <- y - x%*%alpha - d*beta
      mean.1 <- mean.2 - d*lambda
      log.prop <- plogis(z %*% gamma)
      w.mat[, 1] <- log.prop * dnorm(mean.1, 0, sigma)
      w.mat[, 2] <- (1-log.prop) * dnorm(mean.2, 0, sigma)
      l.vec <- rowSums(w.mat)
      w.mat <- w.mat / l.vec
      
      ##Exit loop or not
      penloglik <- sum(log(l.vec)) - p*norm(as.matrix(gamma))
      diff <- penloglik - oldpenloglik
      oldpenloglik <- penloglik
      
      
      if(diff < epsilon.short){
        break
      }
      
      ##M-step
      w <- w.mat[, 1]
      coef_vec <- c(t(x) %*% y, sum(y * d), sum(w * y * d))
      data_mat_w <- cbind(x, d, w * d)
      coef_mat <- t(data_mat_w) %*% data_mat_w
      coef_mat[alpha.dim + 2, alpha.dim + 2] <- sum(w * d)
      const_mat <- matrix(0, nrow = alpha.dim + 2, ncol = alpha.dim + 2)
      const_mat[alpha.dim + 2, alpha.dim + 2] <- 1
      result <- tryCatch({solve.QP(coef_mat, coef_vec, const_mat)}
                         , error = function(err) {coef_mat <- nearPD(coef_mat)$mat
                         return(solve.QP(coef_mat, coef_vec, const_mat))})
      alpha <- result$solution[1:alpha.dim]
      beta <- result$solution[alpha.dim + 1]
      lambda <- result$solution[alpha.dim + 2]
      mean.2 <- y - x %*% alpha - d * beta
      mean.1 <- mean.2 - d*lambda
      sigma <- sqrt(mean(w * (mean.1)^2) + mean((1 - w) * (mean.2)^2))
      
      out.logit <- lbfgs(llogit, dlogit, gamma, invisible=1, orthantwise_c = p, z=z, w=w)
      gamma <- out.logit$par
      
    }
    
    ploglikset[ii] <- penloglik
    theta0[1:alpha.dim, ii] <- alpha
    theta0[alpha.dim+1, ii] <- beta
    theta0[alpha.dim+2, ii] <- lambda
    theta0[alpha.dim+3, ii] <- sigma
    theta0[(alpha.dim+4):(alpha.dim+gamma.dim+3), ii] <- gamma
  }
  
  components <- order(ploglikset, decreasing = TRUE)[1:ninits]
  theta0 <- theta0[, components]
  ploglikset <- double(ninits)
  w.mat <- matrix(nrow=n, ncol=2)
  
  for(ii in 1:ninits){
    alpha <- theta0[1:alpha.dim, ii]
    beta <- theta0[alpha.dim+1, ii]
    lambda <- theta0[alpha.dim+2, ii]
    sigma <- theta0[alpha.dim+3, ii]
    gamma <- theta0[(alpha.dim+4):(alpha.dim+gamma.dim+3), ii]
    diff <- 1.0
    oldpenloglik <- -Inf
    
    for(jj in 1:maxit){
      
      ##E-step
      mean.2 <- y - x%*%alpha - d*beta
      mean.1 <- mean.2 - d*lambda
      log.prop <- plogis(z %*% gamma)
      w.mat[, 1] <- log.prop * dnorm(mean.1, 0, sigma)
      w.mat[, 2] <- (1-log.prop) * dnorm(mean.2, 0, sigma)
      l.vec <- rowSums(w.mat)
      w.mat <- w.mat / l.vec
      
      ##Exit loop or not
      penloglik <- sum(log(l.vec)) - p*norm(as.matrix(gamma))
      diff <- penloglik - oldpenloglik
      oldpenloglik <- penloglik
      
      if(diff < epsilon){
        break
      }
      
      ##M-step
      w <- w.mat[, 1]
      coef_vec <- c(t(x) %*% y, sum(y * d), sum(w * y * d))
      data_mat_w <- cbind(x, d, w * d)
      coef_mat <- t(data_mat_w) %*% data_mat_w
      coef_mat[alpha.dim + 2, alpha.dim + 2] <- sum(w * d)
      const_mat <- matrix(0, nrow = alpha.dim + 2, ncol = alpha.dim + 2)
      const_mat[alpha.dim + 2, alpha.dim + 2] <- 1
      result <- tryCatch({solve.QP(coef_mat, coef_vec, const_mat)}
                         , error = function(err) {coef_mat <- nearPD(coef_mat)$mat
                         return(solve.QP(coef_mat, coef_vec, const_mat))})
      alpha <- result$solution[1:alpha.dim]
      beta <- result$solution[alpha.dim + 1]
      lambda <- result$solution[alpha.dim + 2]
      mean.2 <- y - x %*% alpha - d * beta
      mean.1 <- mean.2 - d*lambda
      sigma <- sqrt(mean(w * (mean.1)^2) + mean((1 - w) * (mean.2)^2))
      
      out.logit <- lbfgs(llogit, dlogit, gamma, invisible=1, orthantwise_c = p, z=z, w=w)
      gamma <- out.logit$par
      
    }
    
    ploglikset[ii] <- penloglik
    theta0[1:alpha.dim, ii] <- alpha
    theta0[alpha.dim+1, ii] <- beta
    theta0[alpha.dim+2, ii] <- lambda
    theta0[alpha.dim+3, ii] <- sigma
    theta0[(alpha.dim+4):(alpha.dim+gamma.dim+3), ii] <- gamma
  }
  
  
  index <- which.max(ploglikset)
  alpha <- theta0[1:alpha.dim, index]
  beta <- theta0[alpha.dim+1, index]
  lambda <- theta0[alpha.dim+2, index]
  sigma <- theta0[alpha.dim+3, index]
  gamma <- theta0[(alpha.dim+4):(alpha.dim+gamma.dim+3), index]
  ploglik <- ploglikset[index]
  
  
  #out.update <- lbfgs(obj, grad, gamma, invisible=1, orthantwise_c = p, z=z, x=x,
  #                    y=y, alpha=alpha, beta=beta, lambda=lambda, d=d, sigma=sigma)
  #value.update <- out.update$value
  #ll.update <- -value.update
  
  par <- list(alpha=alpha, beta=beta, lambda=lambda, sigma=sigma, gamma=gamma)
  list(ploglik=ploglik, par=par) 
}