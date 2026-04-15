data{
  int<lower=1> R;   //total num of person period rows 
  int<lower=1> N;   //num of patients
  int<lower=1> K;   //num of intervals
  int<lower=1> P;   //num of covariates
  array[R] int<lower=0, upper=1> y; // binary outcome
  matrix[R, P] X;   // covariate matrix (person-period)
  array[R] int<lower=1, upper=K> interval;  // interval index
  array[R] int<lower=1, upper=N> patient;  // patient index
}

parameters {
  vector[K] alpha;  // baseline log-odds per interval
  vector[P] beta;   // covariate effects
  vector[N] b;      // patient frailties
  real<lower=0> sigma_b;  // frailty SD
  real<lower=0> sigma_alpha; // RW prior SD
}


model {
  // -- priors 
  // RW(1) baseline hazard
  alpha[1] ~ normal(0, 3); 
  
  for (t in 1:K)
    alpha[t] ~ normal( alpha[t-1], sigma_alpha);
  sigma_alpha ~ cauchy(0, 1);
  
  // covariate effects
  beta ~ normal(0, sqrt(10)); 
  
  // frailties
  b ~ normal(0, sigma_b); 
  sigma_b ~ cauchy(0, 1);
  
  // likelihood
  {
    vector[R] eta; 
    for (r in 1:R)
      eta[r] = alpha[interval[r]] + X[r] * beta + b[patient[r]];
    y ~ bernoulli_logit(eta); 
  }
}






