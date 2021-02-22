library(ggplot2)
library(truncnorm)
library(tidyverse)


# Functions for constraints (Estevão did)
square_root_matrix <- function(x){
  
  # When Q = 1, x will be a scalar
  if (nrow(x) == 1) {return(sqrt(x))}
  
  # When Q > 1, then x will be a matrix
  if (nrow(x) > 1) {
    # Jordan normal form
    X = eigen(x)
    P = X$vectors
    A = diag(X$values)
    
    A_sqrt = diag(sqrt(X$values))
    P_inv = solve(P)
    x_sqrt = P %*% A_sqrt %*%  P_inv
    return(x_sqrt)
  }
}
generate_gamma_delta <- function(INDEX, Q) {
  
  first_row = TRUE
  
  while(first_row) {
    raw_par = matrix(rnorm(INDEX*Q), ncol=Q)
    par_mean  = matrix(rep(apply(raw_par,2,mean), each = nrow(raw_par)), ncol=Q)
    par_aux  = raw_par - par_mean
    
    # Constraints ----
    # apply(par_aux,2,sum)
    parTpar = solve(t(par_aux)%*%(par_aux))
    A = square_root_matrix(parTpar)
    samples = par_aux%*%A
    
    # Force the first to be positive
    for (i in 1:nrow(samples)){
      row1 = samples[1, ]
      if (all(samples[i, ] > 0)) {
        aux = samples[i, ]
        samples[1,] = aux
        samples[i,] = row1
        return(samples)
      }
    }
    # t(samples)%*%samples == 0
    # apply(samples,2,sum) == diag(Q)
  }
}



sim.data <- function(I, J, Q, m, s.mu, s.g, s.e, s.lambda, sigma.E, seed = 123){
  
  # Number of observations
  N <- I*J
  
  # Generate mu
  mu <- rnorm(1, mean = m, sd = s.mu)
  
  # Generate g
  g <- rnorm(I,0,s.g)
  
  # Generate e
  e <- rnorm(J,0,s.e)
  
  # Generate lambda
  lambda <- matrix(sort(rtruncnorm(Q, a=0,mean = 0, sd = s.lambda)),ncol = 1) %>% sort()
  
  
  gamma <- generate_gamma_delta(I, Q)
  delta <- generate_gamma_delta(J, Q)

  # Generate gamma
 # gamma <- matrix(NA, nrow = I ,ncol = Q)
 # gamma[1,] <- rtruncnorm(Q, a=0)
 # for(k in 2:nrow(gamma)){
 #   gamma[k,] <- rnorm(Q)
 # }
 # 
 # # Generate delta
 # delta <- matrix(rnorm(J*Q), nrow = J ,ncol = Q)
 # 
  
  # Now simulate the values
  x <- expand.grid(1:I, 1:J)
  names(x) <- c("gen", "env")
  
  blin <- rep(0, I*J)
  for (k in 1:length(lambda)) {
    blin <- blin + lambda[k]*gamma[x[,1],k]*delta[x[,2],k]
  }
  
  
  mu.ij <- mu + g[x[,1]] + e[x[,2]] + blin
  Y <- rnorm(N, mu.ij, sigma.E)
  
  hyper.par <- list(m = m,
                    s.mu = s.mu, 
                    s.g = s.g, 
                    s.e = s.e, 
                    s.lambda = s.lambda, 
                    sigma.E = sigma.E)

  return(list(Y = Y, x = x,I = I, J = J, Q = Q, mu = mu, g = g, e = e, lambda = lambda, gamma = gamma, delta = delta, blin = blin, hyper.par = hyper.par))
  
}


# Example

sim.test <- sim.data(I = 5,
                     J = 4,
                     Q = 1,
                     m = 90,
                     s.mu = 10,
                     s.g = 10,
                     s.e = 10,
                     s.lambda = 10,
                     sigma.E = 2)

sim.test
# Can create some plots
p.gen <- qplot(x = sim.test$x[,1], y = sim.test$Y, geom = 'boxplot', group = sim.test$x[,1], xlab = 'Genotype', ylab = 'Y')
p.env <- qplot(x = sim.test$x[,2], y = sim.test$Y, geom = 'boxplot', group = sim.test$x[,2], xlab = 'Environment', ylab = 'Y')

gridExtra::grid.arrange(p.gen,p.env)

