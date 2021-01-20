library(R2jags)

bayesian.ammi <- function(data,
                          m = 90,
                          s.mu = 20,
                          s.lambda = 10,
                          s.g_hyperpar = 10,
                          s.e_hyperpar = 10,
                          s.me = 2,
                          n.thin = 1,
                          n.burnin = 2000){
  
  # Specify the Bayesian AMMI model similar to Josse et al (JABES, 2014)
  model_code = '
  model
  {
  # Likelihood
  for (k in 1:N) {
  Y[k] ~ dnorm(mu[k], sigma_E^-2)
  mu[k] = mu_all + g[genotype[k]] + e[environment[k]] + sum(lambda*gamma[genotype[k],1:Q] * delta[environment[k],1:Q])
  }
  # Priors
  mu_all ~ dnorm(m, s.mu^-2) # Prior on grand mean
  for(i in 1:I) {
  g[i] ~ dnorm(0, s.g^-2) # Prior on genotype effect
  }
  for(j in 1:J) {
  e[j] ~ dnorm(0, s.e^-2) # Prior on environment effect
  }
  # Priors on gamma
  for(q in 1:Q) {
  gamma[1, q] ~ dnorm(0, 1)T(0,) # First one is restriced to be positive
  for(i in 2:I) {
  gamma[i, q] ~ dnorm(0, 1) # Prior on genotype interactions
  }
  }
  # Priors on delta
  for(q in 1:Q) {
  for(j in 1:J) {
  delta[j, q] ~ dnorm(0, 1) # Prior on environment interactions
  }
  }
  # Prior on eigenvalues
  for(q in 1:Q) {
  lambda_raw[q] ~ dnorm(0, s.lambda^-2)T(0,)
  }
  lambda = sort(lambda_raw)
  # Prior on residual standard deviation
  sigma_E ~ dunif(0, s.me)
  }
  '
  
  # Get the quantities needed in the JAGS model list
  N = data$I * data$J
  Y = data$Y
  I = data$I
  J = data$J
  Q = data$Q
  genotype = data$x[,'gen']
  environment = data$x[,'env']
  
  # Set up the data
  model_data = list(N = N,
                    Y = Y,
                    I = I,
                    J = J,
                    Q = Q,
                    genotype = genotype,
                    environment = environment,
                    m = m,
                    s.mu = s.mu,
                    s.g = s.g_hyperpar,
                    s.e = s.e_hyperpar,
                    s.lambda = s.lambda,
                    s.me = s.me)
  
  # Choose the parameters to watch
  model_parameters =  c("g", "e", "lambda", "gamma", "delta",
                        'sigma_E', 'mu_all')
  
  # Run the model
  model_run = jags(data = model_data,
                   parameters.to.save = model_parameters,
                   model.file=textConnection(model_code),
                   progress.bar = 'none',
                   n.thin = n.thin,
                   n.burnin = n.burnin,
                   n.iter=2000*n.thin+n.burnin)
  
  return(model_run)
  
}


# Example
data = sim.test
bayesian.ammi(sim.test)

