library(ggplot2)
library(truncnorm)

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
  lambda <- matrix(sort(rtruncnorm(Q, a=0,mean = 0, sd = s.lambda)),ncol = 1)

  # Generate gamma
  gamma <- matrix(NA, nrow = I ,ncol = Q)
  gamma[1,] <- rtruncnorm(Q, a=0)
  for(k in 2:nrow(gamma)){
    gamma[k,] <- rnorm(Q)
  }
  
  # Generate delta
  delta <- matrix(rnorm(J*Q), nrow = J ,ncol = Q)
  
  
  # Now simulate the values
  x <- expand.grid(1:I, 1:J)
  names(x) <- c("gen", "env")
  
  blin <- rep(0, I*J)
  for (k in 1:length(lambda)) {
    blin <- blin + lambda[k]*gamma[x[,1],k]*delta[x[,2],k]
  }
  
  
  mu.ij <- mu + g[x[,1]] + e[x[,2]] + blin
  Y <- rnorm(N, mu.ij, sigma.E)
  

  return(list(Y = Y, x = x,I = I, J = J, Q = Q, mu = mu, g = g, e = e, lambda = lambda, gamma = gamma, delta = delta))
  
}


# Example

sim.test <- sim.data(I = 2,
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

p.gen
p.env


