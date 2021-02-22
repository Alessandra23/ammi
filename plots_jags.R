## Plots bayesian_ammi - Jags code

library(ggplot2)
library(gridExtra)
library(reshape2)
library(patchwork)
theme_set(theme_light())


## Data

data <- sim.test


# nburn    <- object$BUGSoutput$n.burnin
# seq_burn <-  seq(1, nburn, by=1)
object     <- b.ammi
mu.hat     <- object$mu.hat %>% data.frame()
g.hat      <- object$g.hat %>% data.frame()
e.hat      <- object$e.hat %>% data.frame()
gamma.hat  <- object$gamma.hat %>% data.frame()
delta.hat  <- object$delta.hat %>% data.frame()
lambda.hat <- object$lambda.hat %>% data.frame()
sigma.hat  <- object$sigma.hat %>% data.frame()


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


p.mu <- ggplot_gd(mu.hat, data$mu)  + xlab(expression(paste(~mu)))
p.mu
p.g <- ggplot_gd(g.hat, data$g) + xlab("Genotype")
p.e <- ggplot_gd(e.hat, data$e) + xlab("Environment")
p.gamma <- ggplot_gd(gamma.hat, data$gamma) + xlab(expression(paste(~gamma))) 
p.delta <- ggplot_gd(delta.hat, data$delta) + xlab(expression(paste(~delta))) 
p.lambda <- ggplot_gd(lambda.hat, data$lambda) + xlab(expression(paste(~lambda))) 
p.sigma <- ggplot_gd(sigma.hat, data$hyper.par$sigma.E) + xlab(expression(paste(~sigma^2)))


if(data$Q == 1){
  all.plots <- (p.g / p.e / p.gamma / p.delta) | (p.mu / p.lambda / p.sigma)
}else{
  all.plots <- (p.g / p.e / p.gamma / p.delta / p.lambda) | (p.mu / p.sigma)
}

all.plots



## See the bilinear term

p.blin <- function(object, data, i, j){
  Q = data$Q
  lambda.hat <- object$lambda.hat
  gamma.hat  <- object$gamma.hat 
  delta.hat  <- object$delta.hat

  if(Q == 1){
    prod.lgd <- lambda.hat*gamma.hat[,i]*delta.hat[,j] %>% data.frame()
    p.lgd <- ggplot_gd(prod.lgd, data$blin[i]*blin[j])  + xlab(expression(paste(~lambda~gamma[i]~delta[j])))
  } else{
    prod.lgd = c(0,Q)
    for (k in 1:Q) {
      prod.lgd <- prod.lgd + lambda.hat[,k]*gamma.hat[,paste0("gamma[", i, ',', k,"]")]*delta.hat[,paste0("delta[", j, ',', k,"]")] %>% data.frame()
    }
    p.lgd <- ggplot_gd(prod.lgd, data$blin[i]*data$blin[j])  + xlab(bquote(~Sigma[q]~lambda[q]~gamma[.(i)~q]~delta[.(j)~q]))  #xlab(expression(paste(~lambda~gamma[.(i)]~delta["j"])))
  }
  
  options(warn=-1)
  return(p.lgd)
}

p.blin(object, data, 2,3)

blin.gamma <- lapply(1:data$J, p.blin, object = object, data = data, i = 1)
blin.delta <- lapply(1:data$I, p.blin, object = object, data = data, j = 1)
do.call(grid.arrange, c(blin.gamma, ncol = length(blin.gamma)))
do.call(grid.arrange, c(blin.delta, ncol = length(blin.delta)))

