library(rstan)

data <- sim.test


ammi_stan_model <- stan_model(file = "ammi_stan.stan")

ammi_stan_data <- list(
  N = data$I*data$J,
  Y = data$Y,
  I = data$I,
  J = data$J,
  Q = data$Q,
  m = data$hyper.par$m,
  genotype = data$x[,'gen'],
  environment = data$x[,'env'],
  s_mu = data$hyper.par$s.mu,
  s_g = data$hyper.par$s.g,
  s_e = data$hyper.par$s.e,
  s_lambda = data$hyper.par$s.lambda,
  a_sigma = 2,
  b_sigma = 0.5
)



ammi_stan_fitted_model <- sampling(ammi_stan_model,
                                   data = ammi_stan_data,
                                   chains = 1,
                                   warmup = 3000,
                                   iter = 10000,
                                   thin = 1)


#ammi_fit_stan <- stan(
#  file = "ammi_stan.stan",
#  data = ammi_stan_data,
#  chains = 3
#)
#

mu.hat.stan     <- extract(ammi_stan_fitted_model, "mu_all", permuted = F) %>% data.frame()
g.hat.stan      <- extract(ammi_stan_fitted_model, 'g', permuted = F) %>% data.frame()
e.hat.stan      <- extract(ammi_stan_fitted_model, 'e', permuted = F) %>% data.frame()
gamma.hat.stan  <- extract(ammi_stan_fitted_model, 'gamma', permuted = F) %>% data.frame()
delta.hat.stan  <- extract(ammi_stan_fitted_model, 'delta', permuted = F) %>% data.frame()
lambda.hat.stan <- extract(ammi_stan_fitted_model, 'lambda', permuted = F) %>% data.frame()
sigma.hat.stan  <- extract(ammi_stan_fitted_model, 'sigma2', permuted = F) %>% data.frame()




# Function to create plots

ggplot_gd <- function(obj1, obj2){
  p <- ggplot(data = melt(obj1), aes(x=variable, y=value)) 
  p <- p + geom_boxplot(aes(fill=variable) ,fill = "white", color = "black") + 
    geom_point(data = data.frame(x = names(obj1), y = as.vector(obj2)), aes(x=x, y = y), color = 'red') +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title=element_text(size=14)) + #,face="bold"))+
    ylab(" ")
  options(warn=-1)
  return(p)
}


p.mu.stan <- ggplot_gd(mu.hat.stan, data$mu)  + xlab(expression(paste(~mu)))
p.g.stan <- ggplot_gd(g.hat.stan, data$g) + xlab("Genotype")
p.e.stan <- ggplot_gd(e.hat.stan, data$e) + xlab("Environment")
p.gamma.stan <- ggplot_gd(gamma.hat.stan, data$gamma) + xlab(expression(paste(~gamma))) 
p.delta.stan <- ggplot_gd(delta.hat.stan, data$delta) + xlab(expression(paste(~delta))) 
p.lambda.stan <- ggplot_gd(lambda.hat.stan, data$lambda) + xlab(expression(paste(~lambda))) 
p.sigma.stan <- ggplot_gd(sigma.hat.stan, data$hyper.par$sigma.E) + xlab(expression(paste(~sigma^2)))


if(data$Q == 1){
  all.plots.stan <- (p.g.stan / p.e.stan / p.gamma.stan / p.delta.stan) | (p.mu.stan / p.lambda.stan / p.sigma.stan)
}else{
  all.plots.stan <- (p.g.stan / p.e.stan / p.gamma.stan / p.delta.stan / p.lambda.stan) | (p.mu.stan / p.sigma.stan)
}

all.plots.stan





# mu_pos <- rstan::extract(ammi_fit_stan, "mu_all", permuted = F) %>% data.frame()
# g_pos <- rstan::extract(ammi_fit_stan, 'g', permuted = F) %>% data.frame()
# beta_pos <- rstan::extract(ammi_fit_stan, 'e', permuted = F) %>% data.frame()
# lambda_pos <- rstan::extract(ammi_fit_stan, 'lambda', permuted = F) %>% data.frame()
# gamma_pos <- rstan::extract(ammi_fit_stan, 'gamma', permuted = F) %>% data.frame()
# delta_pos <- rstan::extract(ammi_fit_stan, 'delta', permuted = F) %>% data.frame()
# sigma_E_pos <- rstan::extract(ammi_fit_stan, 'sigma2', permuted = F) %>% data.frame()




