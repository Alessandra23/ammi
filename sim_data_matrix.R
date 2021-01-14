# Simulate data - matrix form

library(ggplot2)
library(truncnorm)
library(gridExtra)

sim.data.matrix <- function(values){
  
  seed <- values$seed
  I <- values$I
  J <- values$J
  Q <- values$Q
  N <- I*J
  
  m <- values$m
  s.mu <- values$s.mu
  s.alpha <- values$s.alpha
  s.beta <- values$s.beta
  s.lambda <- values$s.lambda
  sigma.E <- values$sigma.E
  
  # Generate mu
  mu <- rnorm(1, mean = m, sd = s.mu)
  
  # Generate alpha
  alpha <- rnorm(I,0,s.alpha)
  
  # Generate beta
  beta <- rnorm(J,0,s.beta)
  
  # Generate lambda
  lambda <- matrix(sort(rtruncnorm(Q, a=0,mean = 0, sd = s.lambda)),ncol = 1)
  #lambda <- (sort(rtruncnorm(Q, a=0,mean = 0, sd = s.lambda)))
  # Generate gamma
  gamma <- matrix(NA, nrow = I ,ncol = Q)
  gamma[1,] <- rtruncnorm(Q, a=0)
  for(k in 2:nrow(gamma)){
    gamma[k,] <- rnorm(Q)
  }
  
  # Generate delta
  delta <- matrix(rnorm(J*Q), nrow = J ,ncol = Q)
  
  
  # Now simulate the values
  G_by_E = expand.grid(1:I, 1:J)
  # consertar isso
  for (k in 1:length(lambda)) {
    blin <- lambda[k]*gamma[G_by_E[,1],k]*delta[G_by_E[,2],k]
  }
  
  mu_ij = mu + alpha[G_by_E[,1]] + beta[G_by_E[,2]] + blin
  Y2 = rnorm(N, mu_ij, sigma.E)
  
  I1 <- rep(1,I)
  J1 <- rep(1,J)
  
  mu.Y <- mu*I1%*%t(J1) + kronecker(alpha,t(J1)) + kronecker(t(beta), (I1)) + gamma%*%diag(lambda)%*%t(delta)
  Y <- rnorm(N,c(t(mu.Y)),sigma.E)
  
  # Can create some plots
  p.gen <- qplot(x = G_by_E[,1], y = Y, geom = 'boxplot', group = G_by_E[,1], xlab = 'Genotype')
  p.env <- qplot(x = G_by_E[,2], y = Y, geom = 'boxplot', group = G_by_E[,2], xlab = 'Environment')
  
  p.gen2 <- qplot(x = G_by_E[,1], y = Y2, geom = 'boxplot', group = G_by_E[,1], xlab = 'Genotype')
  p.env2 <- qplot(x = G_by_E[,2], y = Y2, geom = 'boxplot', group = G_by_E[,2], xlab = 'Environment')
  
  comp.gen <- grid.arrange(p.gen, p.gen2)
  comp.env <- grid.arrange(p.env, p.env2)
  
  grid.arrange(p.gen,p.env)
  
  return(list(Y = Y2, Y2 = Y, p.gen = p.gen, p.env = p.env, p.gen2 = p.gen2, p.env2 = p.env2, comp.gen = comp.gen, comp.env = comp.env))
  
}


values <- list(seed = set.seed(123),
               I = 2,
               J = 4,
               Q = 2,
               m = 90,
               s.mu = 10,
               s.alpha = 10,
               s.beta = 10,
               s.lambda = 10,
               sigma.E = 2)


sim.data(values = values)
