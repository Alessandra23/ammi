## Complete Model

rm(list = ls())
library(R2jags)
library(ggplot2)
library(MCMCvis)
require(reshape2)
library(tidyverse)
library(gridExtra)
theme_set(theme_light())

# Specify fixed values
Q = 1 # Number of components
I = 5 # Number of genotypes
J = 9 # Number of environments
N = I*J # Total number of obs
m = 90
s_mu = 20
s_alpha = 10
s_beta = 10
s_lambda = 10
#S_ME = 2
a = 0.1
b = 0.1

# Some further fixed values
mu = 100
sigma_E = 3/2 # Not sure why S_ME was specified if they're also giving sigma_E
alpha = c(-1, -1, 0, 1, 1)
beta = -4:4
lambda = 1
gamma = seq(2, -2)/sqrt(10)
delta = c(0.5, 0.5, rep(0, 5), -0.5, -0.5)
# delta = c(1, 0.5, rep(0, 5), -0.5, -1) 

# Now simulate the values
set.seed(123)
G_by_E = expand.grid(1:I, 1:J)
mu_ij = mu + alpha[G_by_E[,1]] + beta[G_by_E[,2]] + lambda * gamma[G_by_E[,1]] * delta[G_by_E[,2]]
Y = rnorm(N, mu_ij, sigma_E)

# Can create some plots
qplot(x = G_by_E[,1], y = Y, geom = 'boxplot', group = G_by_E[,1], xlab = 'Genotype')
qplot(x = G_by_E[,2], y = Y, geom = 'boxplot', group = G_by_E[,2], xlab = 'Environment')


# Q = 1 model -------------------------------------------------------------

# This is the simple Q = 1 model - only here for understanding
# Jags code to fit the model to the simulated data
model_code = '
model
{
  # Likelihood
  for (k in 1:N) {
    Y[k] ~ dnorm(mu[k], sigma_E^-2)
    mu[k] = mu_all + alpha[genotype[k]] + beta[environment[k]] + lambda * gamma[genotype[k]] * delta[environment[k]]
  }
  # Priors
  mu_all ~ dnorm(m, s_mu^-2) # Prior on grand mean
  for(i in 1:I) {
    alpha[i] ~ dnorm(0, s_alpha^-2) # Prior on genotype effect
  }
  gamma[1] ~ dnorm(0, 1)T(0,)
  for(i in 2:I){
  gamma[i] ~ dnorm(0, 1) # Prior on genotype interactions
  }
  for(j in 1:J) {
    beta[j] ~ dnorm(0, s_beta^-2) # Prior on environment effect
    delta[j] ~ dnorm(0, 1)
  }
  # Prior on first (and only) eigenvalue
  lambda ~ dnorm(0, s_lambda^-2)T(0,)
  # Prior on residual standard deviation
  sigma_E ~ dgamma(a, b)
}
'