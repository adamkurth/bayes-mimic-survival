data {
  int<lower=1> R;
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> P;
  array[R] int<lower=0, upper=1> y;
  matrix[R, P] X;
  array[R] int<lower=1, upper=K> interval;
  array[R] int<lower=1, upper=N> patient;
}
parameters {
  vector[K] alpha;
  vector[P] beta;
  vector[N] b_raw;
  real<lower=0> sigma_b;
  real<lower=0> sigma_alpha;
}
transformed parameters {
  vector[N] b = sigma_b * b_raw;   // NON-CENTERED
}
model {
  alpha[1] ~ normal(0, 3);
  for (t in 2:K) alpha[t] ~ normal(alpha[t-1], sigma_alpha);
  sigma_alpha ~ student_t(3, 0, 1);

  beta    ~ normal(0, sqrt(10));
  b_raw   ~ std_normal();
  sigma_b ~ exponential(1);        // PC prior: shrinks toward sigma_b = 0

  {
    vector[R] eta;
    for (r in 1:R) eta[r] = alpha[interval[r]] + X[r] * beta + b[patient[r]];
    y ~ bernoulli_logit(eta);
  }
}
