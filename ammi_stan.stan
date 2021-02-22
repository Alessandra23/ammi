data {
  int<lower=0> N;
  int<lower=0> I;
  int<lower=0> J;
  int<lower=0> Q;
  vector[N] Y;
  real m;
  int genotype[N];
  int environment[N];
  real<lower=0> s_mu;
  real<lower=0> s_g;
  real<lower=0> s_e;
  real<lower=0> s_lambda;
  real<lower=0> a_sigma;
  real<lower=0> b_sigma;
}

parameters {
  real<lower=0> sigma2;
  matrix[I,Q] gamma;
  matrix[J,Q] delta;
  vector[I] z_g;
  vector[J] z_e;
  vector[Q] z_lambda;
  real z_mu;
  real<lower=0> sigma_E;
}

transformed parameters{
  real mu_all = (s_mu)*z_mu + m;
  vector[I] g = (s_g)*z_g;
  vector[J] e = (s_e)*z_e;
  //vector[Q] lambda_no_sort = (s_lambda)*z_lambda;
  vector[Q] lambda = sort_asc(z_lambda);
}

model {
  
  for(n in 1:N){
    real mu = 0;
    vector[Q] blin ;
    
    // Compute the bilinear term 
    for (k in 1:Q) {
       blin[k] = lambda[k]*gamma[genotype[n],k]*delta[environment[n],k];
    } 
    
    mu = mu_all + g[genotype[n]] + e[environment[n]] + sum(blin);
    Y[n] ~ normal(mu, sqrt(sigma2));
  }
  
  // Prior on gamma
  
  for(q in 1:Q) {
   gamma[1, q] ~ normal(0,1)T[0, ]; // First one is restriced to be positive
   for(i in 2:I) {
      gamma[i, q] ~ normal(0,1); 
   }
  }
  
  // Prior on delta
  
  for(q in 1:Q) {
    for(j in 1:J){
      delta[j, q] ~ normal(0, 1); 
    }
  }
  
  // Prion on environment effect
  
  for(j in 1:J){
    z_e[j] ~ normal(0,1);
  }
  
  // Prior on genotype effect
  
  for(i in 1:I){
    z_g[i] ~ normal(0,1);
  }
  
  
  // Prior on lambda
  
  for(q in 1:Q) {
    z_lambda[q] ~ normal(0,s_lambda)T[0, ];
  }
  
  // Prion on grand mean
  
  z_mu ~ normal(0,1);
  
  // Prion on variance
  
  sigma2 ~ inv_gamma(a_sigma,b_sigma);
}

