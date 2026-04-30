// mnar_sensitivity.stan
// ============================================================================
// MNAR model for discrete-time survival with informative censoring.
//
// selection factorization:
//   f(m, y, x | omega) = f(m | y, x; xi) * f(y | x; eta) * f(x; theta)
//
// Here:
//   - Y_{i,t} is the person-period mortality outcome
//   - M_i = 1 if outcome fully observed (died or survived past window)
//   - M_i = 0 if censored (discharged alive)
//   - gamma_sens plays the role of "xi_2": the sensitivity parameter
//     encoding the relationship between M and Y conditional on X
//
// gamma_sens = 0, -> reduces to MAR (ignorable censoring).
// gamma_sens != 0,  -> censoring is informative (MNAR).
//
// The prior on gamma_sens encodes beliefs about the degree of MNAR:
//   gamma_sens ~ N(0, sigma_gamma)  centers on MAR with uncertainty
// ============================================================================

data {
  int<lower=1> R;                              // person-period rows
  int<lower=1> K;                              // number of intervals
  int<lower=1> P;                              // number of covariates
  array[R] int<lower=0, upper=1> y;            // event indicator (0/1)
  matrix[R, P] X;                              // covariate matrix
  array[R] int<lower=1, upper=K> interval;     // interval index for each row
  vector[R] censored_terminal;                 // 1 = terminal row of a censored patient, 0 otherwise
  real<lower=0> sigma_gamma;                   // prior SD for gamma_sens (default: 1.0)
}

parameters {
  vector[K] alpha;                             // baseline log-odds per interval
  vector[P] beta;                              // covariate effects
  real gamma_sens;                             // sensitivity parameter ("xi_2")
  real<lower=0> sigma_alpha;                   // RW innovation SD
}

model {
  // --- Priors ---
  alpha[1] ~ normal(0, 3);
  for (t in 2:K) alpha[t] ~ normal(alpha[t-1], sigma_alpha);
  sigma_alpha ~ student_t(3, 0, 1);
  beta ~ normal(0, sqrt(10));

  // Sensitivity parameter: centered on MAR (gamma=0)
  // sigma_gamma controls how much departure from MAR is plausible a priori
  gamma_sens ~ normal(0, sigma_gamma);

  // --- Likelihood: selection factorization ---
  // The hazard for each person-period row includes the MNAR offset (violation of MAR)
  // at the terminal interval of censored patients
  {
    vector[R] eta;
    for (r in 1:R)      // r = person period rows, X[r] fixed, o_it = I(C_i < T_i) in linear predictor
      eta[r] = alpha[interval[r]] + X[r] * beta + gamma_sens * censored_terminal[r];
    y ~ bernoulli_logit(eta);
  }
}

generated quantities {
  // Export key parameters for post-processing
  vector[P] beta_out = beta;
  real gamma_out = gamma_sens;
}
