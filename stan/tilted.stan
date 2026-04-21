data {
  int<lower=1> R;
  int<lower=1> K;
  int<lower=1> P;
  array[R] int<lower=0, upper=1> y;
  matrix[R, P] X;
  array[R] int<lower=1, upper=K> interval;
  vector[R] offset_vec;
}
parameters {
  vector[K] alpha;
  vector[P] beta;
  real<lower=0> sigma_alpha;
}
model {
  alpha[1] ~ normal(0, 3);
  for (t in 2:K) alpha[t] ~ normal(alpha[t-1], sigma_alpha);
  sigma_alpha ~ student_t(3, 0, 1);
  beta ~ normal(0, sqrt(10));

  {
    vector[R] eta;
    for (r in 1:R) eta[r] = alpha[interval[r]] + X[r] * beta + offset_vec[r];
    y ~ bernoulli_logit(eta);
  }
}
generated quantities {
  vector[P] beta_out = beta;
}




