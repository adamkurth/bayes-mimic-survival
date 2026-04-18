#---------------------------------------------------------
# core.r (core functions)
#---------------------------------------------------------
# sourced by pg_demo.r and stan_demo.r

# contains FUNCTIONS ONLY -- no executable code
# contents:
#   * generate_synthetic_icu()   – synthetic ICU cohort with informative censoring
#   * expand_person_period()     – person-period (long-form) design matrix
#   * pg_gibbs_standard()        – PG Gibbs sampler, no frailty
#   * pg_gibbs_frailty()         – PG Gibbs sampler, shared frailty (centered)
#   * pg_gibbs_tilted()          – PG Gibbs sampler with exponential-tilting offset
#   * make_offset_vector()       – helper for sensitivity analysis offsets
#
# Usage:
#   source("core.r")
#---------------------------------------------------------

library(BayesLogit); library(Matrix); library(coda)

# 1. synthetic data generation --------------------------------------------------
# generate N patients with static covariates and discrete-time survival outcomes
# embed informative censoring: discharge probability depends on latent (b_i)


generate_synthetic_icu <- function( N = 500, K = 30, delta = 1, sigma_b_true = 0.6,  cens_frailty_coef = -0.8) {
    # N                 : number of patients
    # K                 : max discrete intervals (e.g., 30 days or ~4 weeks)
    # delta             : interval width in days (1 = daily, 7 = weekly)
    # sigma_b_true      : SD of the latent frailty
    # cens_frailty_coef : coefficient linking frailty to censoring hazard
    #                     (negative => healthier patients discharged earlier)

    # static X_i
    age                 <- rnorm(N, mean = 65, sd = 15)
    age                 <- pmin( pmax(age, 18), 90) # truncate to [18,90]
    male                <- rbinom(n=N, size=1, prob=0.5)
    emergency           <- rbinom(n=N, size=1, prob = 0.65)
    # icd codes (comorbidities)
    icd_sepsis          <- rbinom(n=N, size=1, prob = 0.25) # sepsis
    icd_heart_failure   <- rbinom(n=N, size=1, prob = 0.20) # heart failure
    icd_resp_failure    <- rbinom(n=N, size=1, prob = 0.20) # respiratory failure
    icd_aki             <- rbinom(n=N, size=1, prob = 0.15) # acute kidney injury
    gcs_impaired        <- rbinom(n=N, size=1, prob = 0.20) # impaired consciousness (GCS <= 12)
    on_vasopressors     <- rbinom(n=N, size=1, prob = 0.30) # on vasopressors at admission
    age_std <- (age - mean(age)) / sd(age) # standardize age

    X <- cbind(
        age_std        = age_std,
        male           = male,
        emergency      = emergency,
        icd_sepsis     = icd_sepsis,
        icd_hf         = icd_heart_failure,
        icd_resp       = icd_resp_failure,
        icd_aki        = icd_aki,
        gcs_impaired   = gcs_impaired,
        vasopressor    = on_vasopressors
    )
    P <- ncol(X)

    # ---- true parameters
    # baseline hazard: smooth increasing then plateau
    alpha_true <- -3.5 + 0.8 * log(1:K) - 0.02 * (1:K)

    # covariate effects (log odds ratios)
    beta_true <- c(
        age_std = 0.25,         # older increases risk
        male = 0.10,            # slight higher risk
        emergency = 0.30,       # emergency much worse
        icd_sepsis = 0.60,      # sepsis large effect
        icd_resp = 0.35,        # respiratory moderate
        icd_hf = 0.40,          # heart failure moderate
        icd_aki = 0.25,         # acute kidney injury
        gcs_impaired = 0.50,    # impaired consciousness large effect
        vasopressor = 0.45      # on vasopressors large effect
    )

    # ----- latent frailty (b_i)
    b_true <- rnorm(n=N, mean=0, sd=sigma_b_true) # patient-specific frailty (log hazard multiplier)

    # ----- generate latent survival
    T_latent <- rep(K + 1, N)

    for ( i in 1:N) {

        # linear predictor: eta_{i,t} = alpha_t + X_i * beta + b_i
        eta_i <- alpha_true + as.numeric(X[i, ] %*% beta_true) + b_true[i]
        h_i <- plogis(eta_i) # hazard probability at each time t

        for (t in 1:K){
            if (runif(1) < h_i[t]){ # event occurs at time t
                T_latent[i] <- t
                break
            }
        }
    }

    # ----- generate informative censoring
    # censoring depends on frailty:
    #   - healthier (b_i < 0) discharged earlier
    #   - sicker (b_i > 0) discharged later

    C_latent <- rep(K + 1, N)

    for (i in 1:N) {
        # censoring hazard: depends on b_i (informative!)
        # FIX: use cens_frailty_coef parameter (was hardcoded as -0.8)
        cens_logit <- -2.5 - 0.3 * log(1:K) + cens_frailty_coef * b_true[i]
        # when cens_frailty_coef = -0.8 (default): healthier patients (b_i < 0)
        # have HIGHER censoring hazard (more likely discharged early)
        # when cens_frailty_coef = 0: censoring is INDEPENDENT of frailty

        cens_h <- plogis(cens_logit)

        for (t in 1:K){
            if (runif(1) < cens_h[t]){ # censored at time t
                C_latent[i] <- t
                break
            }
        }

    }

    y_obs <- pmin(T_latent, C_latent) # observed time is min of event or censoring
    delta_i <- as.numeric(T_latent <= C_latent) # event indicator: 1 if event observed, 0 if censored


    list(
        N = N, K = K, P = P, delta_width = delta,
        X = X, y_obs = y_obs, delta = delta_i,
        T_latent = T_latent, C_latent = C_latent,
        b_true = b_true,
        alpha_true = alpha_true, beta_true = beta_true,
        sigma_b_true = sigma_b_true
    )
}

expand_person_period <- function(dat){
    # dat: list output from generate_synthetic_icu()
    # returns list with:
    #   pp: person-period data frame with columns (i, t, y_it, X_i)
    #   X.ast: design matrix for logistic regression with interval dummies and covariates
    #   y: vector of binary outcomes (event indicator for each person-period)
    #   patient_id: vector of patient IDs for each row
    #   interval: vector of interval numbers for each row
    #   R: total number of person-period rows
    #   K: number of intervals
    #   P: number of covariates

    # Initialize vectors to store results
    rows_i <- NULL      # patient IDs
    rows_t <- NULL      # time intervals
    rows_y <- NULL      # outcome indicators
    rows_X <- NULL      # covariates

    for (i in 1:dat$N){
        Ki <- min(dat$y_obs[i], dat$K) # cap at K; number of total rows for patient i

        for (t in 1:Ki){
            rows_i <- c(rows_i, i)
            rows_t <- c(rows_t, t)
            y_it <- as.integer(t == Ki && dat$delta[i] == 1)
            rows_y <- c(rows_y, y_it)
            rows_X <- rbind(rows_X, dat$X[i, ])
        }
    }

    # Create person-period data frame
    pp <- data.frame(
        i = rows_i,
        t = rows_t,
        y_it = rows_y,
        rows_X
    )

    # build design matrix X_ast = [D_t, X_i], D_t interval dummies
    K <- dat$K
    P <- dat$P
    R <- nrow(pp)

    D <- matrix(0, nrow=R, ncol=K) # interval dummies
    colnames(D) <- paste0("alpha_", 1:K)

    for (r in 1:R) D[r, pp$t[r]] <- 1 # set D_{r,t} = 1 if row r corresponds to interval t

    # Extract covariates and ensure proper matrix form
    X_covariates <- as.matrix(pp[, 4:ncol(pp), drop = FALSE])
    X_ast <- cbind(D, X_covariates)
    colnames(X_ast) <- c(colnames(D), colnames(X_covariates))

    list(
        pp = pp,
        X.ast = X_ast,
        y = pp$y_it,
        patient_id = pp$i,
        interval = pp$t,
        R = R, K = K, P = P, N = dat$N
    )
}


# 2. PG Gibbs sampler (no frailty) --------------------------------------------------

# Full conditionals (Polson, Scott & Windle 2013):
#   omega_r | theta   ~  PG(1, S_r theta)
#   theta   | omega,y ~  N( (Sigma_0^{-1} + S'OmegaS)^{-1}(Sigma_0^{-1}mu_0 + S'kappa),
#                           (Sigma_0^{-1} + S'OmegaS)^{-1} )
# where kappa_r = y_r - 1/2.

pg_gibbs_standard <- function(  pp_data, n_iter = 2000, burn_in = 500,
                                tau2 = 10, sigma_alpha2 = 1, verbose = TRUE) {
    # pp_data: list output from expand_person_period()
    # n_iter: total number of MCMC iterations
    # burn_in: number of burn-in iterations to discard
    # tau2: prior variance for betas (covariate effects)
    # sigma_alpha2: prior variance for alphas (baseline hazard)

    X.ast <- pp_data$X.ast
    y <- pp_data$y
    R <- pp_data$R
    K <- pp_data$K
    P <- pp_data$P
    dim_theta <- K + P # number of parameters: K alphas + P betas

    # priors
    mu_0 <- rep(0, dim_theta) # prior mean, all 0
    # simple use diag prior
    Sigma_0_inv <- diag( c( rep( 1/sigma_alpha2, K), rep(1/tau2, P) ) )

    # storage
    theta_samples <- matrix(0, nrow=n_iter, ncol=dim_theta)
    theta <- rep(0, dim_theta)                  # start at 0 (log-odds = 0 => p=0.5)
    kappa <- y - 0.5                            # kappa_i = y_i - 1/2

    if (verbose) {
        cat("Running standard PG Gibbs sampler...\n")
        pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
    }

    for (iter in 1:n_iter){

        # 1: sample PG augmentation
        eta <- as.numeric(X.ast %*% theta) # lin pred
        omega <- rpg(n=R, h=1, z=eta) # omega_i ~ PG(1, eta_i)

        # 2: sample: theta | omega, kappa ~ N(mu_theta, Sigma_theta)
        # Sigma_theta = (X'*Omega*X + Sigma_0_inv)^(-1)
        # mu_theta = Sigma_theta * (X'*kappa + Sigma_0_inv*mu_0)
        # posterior: N(mu_theta, Sigma_theta)

        X0X <- crossprod(X.ast * omega, X.ast) # X'*Omega*X
        Sigma_theta_inv <- X0X + Sigma_0_inv
        Sigma_theta <- solve(Sigma_theta_inv) # posterior covariance
        mu_theta <- Sigma_theta %*% ( Sigma_0_inv %*% mu_0 + crossprod(X.ast, kappa) ) # posterior mean

        theta <- as.numeric( mu_theta + chol(Sigma_theta) %*% rnorm(dim_theta) ) # sample from N(mu_theta, Sigma_theta)
        theta_samples[iter, ] <- theta

        if (verbose) setTxtProgressBar(pb, iter)
    }
    if (verbose) close(pb)

    # discard burn-in
    keep <- (burn_in + 1):n_iter
    samples <- theta_samples[keep, ]
    colnames(samples) <- colnames(pp_data$X.ast)

    list(
        samples = samples,
        alpha   = samples[, 1:K, drop = FALSE],
        beta    = samples[, (K+1):(K+P), drop = FALSE]
    )
}


# 3. PG Gibbs sampler with shared frailty (centered) --------------------------------------------------

# Adds patient-level random intercept b_i ~ N(0, sigma_b^2) with conjugate
# hyperprior sigma_b^2 ~ InvGamma(a0, b0).

# known limitation: centered parameterization (b, sigma_b^2) exhibits Neal's funnel geometry, that mixes very slowly in Gibbs.
# This is an intentional demo to illustrate the problem; stan_demo.r provides a non-centered parameterization that fixes this issue.

pg_gibbs_frailty <- function(pp_data, dat, n_iter = 2000, burn_in = 500,
                            tau2 = 10, sigma_alpha2 = 1, a0 = 0.01, b0 = 0.01,
                            verbose = TRUE) {
    X.ast <- pp_data$X.ast
    y <- pp_data$y
    R <- pp_data$R
    K <- pp_data$K
    P <- pp_data$P
    N <- pp_data$N
    pid <- pp_data$patient_id
    dim_theta <- K + P

    # priors
    mu_0 <- rep(0, dim_theta) # prior mean, all 0
    Sigma_0_inv <- diag( c( rep(1/sigma_alpha2, K), rep(1/tau2, P) ) )

    # hyperprior for sigma_b^2: Inverse-Gamma(a,b)
    a0 <- 0.01; b0 <- 0.01

    # storage
    theta_samples <- matrix(data=0, nrow=n_iter, ncol=dim_theta)
    sigma_b_samples <- numeric(n_iter)
    b_samples <- matrix(data=0, nrow=n_iter, ncol=N)

    # initialize
    theta <- rep(0, dim_theta)      # fixed effects start at 0
    b <- rep(0, N)                  # frailty starts at 0
    sigma_b2 <- 1                   # initial value for sigma_b^2
    kappa <- y - 0.5                # kappa_i = y_i - 1/2

    # precompute patient membership: which rows belong to patient i
    patient_rows <- split(1:R, pid)
    if (verbose) {
        cat("Running frailty PG Gibbs sampler (centered)...\n")
        pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
    }

    for( iter in 1:n_iter){

        #  1. sample PG augmentation
        b_expanded <- b[pid] # expand b_i to row-level
        eta <- as.numeric(X.ast %*% theta + b_expanded) # lin pred with frailty
        omega <- rpg(n=R, h=1, z=eta) # omega_i ~ PG(1, eta_i)

        #  2. sample theta | omega, kappa, b ~ N(mu_theta, Sigma_theta)
        kappa_tilde <- kappa - omega * b_expanded # adjust kappa for frailty
        X0X <- crossprod(X.ast * omega, X.ast) # X'*Omega*X
        Sigma_theta_inv <- X0X + Sigma_0_inv
        Sigma_theta <- solve(Sigma_theta_inv) # posterior covariance
        mu_theta <- Sigma_theta %*% ( Sigma_0_inv %*% mu_0 + crossprod(X.ast, kappa_tilde) ) # posterior mean
        theta <- as.numeric( mu_theta + chol(Sigma_theta) %*% rnorm(dim_theta) )


        # 3. sample b_i | omega, kappa, theta, sigma_b^2 ~ N(mu_bi, sigma_bi^2)
        eta_fixed <- as.numeric(X.ast %*% theta) # fixed part
        for (i in 1:N) {
            rows_i <- patient_rows[[i]]
            omega_i <- omega[rows_i]
            kappa_i <- kappa[rows_i]
            eta_fix_i <- eta_fixed[rows_i]
            # posterior precision and mean for b_i
            prec_bi <- 1/sigma_b2 + sum(omega_i) # precision from prior + likelihood
            var_bi <- 1 / prec_bi
            mu_bi <- var_bi * sum(kappa_i - omega_i * eta_fix_i)
            b[i] <- rnorm(1, mean = mu_bi, sd = sqrt(var_bi))

        }

        # 4. sample sigma_b^2 | b ~ Inverse-Gamma(a', b')
        a_post <- a0 + N / 2
        b_post <- b0 + sum((b)^2) / 2
        sigma_b2 <- 1 / rgamma(1, shape = a_post, rate = b_post)

        # store
        theta_samples[iter, ]   <- theta
        sigma_b_samples[iter]   <- sigma_b2
        b_samples[iter, ]       <- b
        if (verbose) setTxtProgressBar(pb, iter)
    }
    if (verbose) close(pb)

    keep <- (burn_in + 1):n_iter
    samples <- theta_samples[keep, ]
    colnames(samples) <- colnames(pp_data$X.ast)

    list(
        samples  = samples,
        alpha    = samples[, 1:K, drop = FALSE],
        beta     = samples[, (K+1):(K+P), drop = FALSE],
        sigma_b2 = sigma_b_samples[keep],
        b        = b_samples[keep, ]
    )
}


# 4. PG Gibbs sampler with exponential-tilting offset --------------------------------------------------

# Exponential tilting: adds offset o_{i,t} = gamma * I(C_i = t, delta_i = 0)
# to the linear predictor, only on the terminal censoring interval.

make_offset_vector <- function(pp_data, dat, gamma_val) {
    # pp_data: list output from expand_person_period()
    # dat: original data list from generate_synthetic_icu()
    # gamma_val: value of gamma for the exponential tilting offset

  R <- pp_data$R
  N <- dat$N
  pid <- pp_data$patient_id
  offset <- rep(0, R)

  for (i in 1:N) {

    if (dat$delta[i] == 0) {

      rows_i   <- which(pid == i) # rows corresponding to patient i
      terminal <- rows_i[length(rows_i)] # censoring interval
      offset[terminal] <- gamma_val # set offset
    }

  }
  offset
}


pg_gibbs_tilted <- function(pp_data, dat, gamma_val, n_iter = 2000, burn_in = 500,
                            tau2 = 10, sigma_alpha2 = 1,
                            a0 = 0.01, b0 = 0.01, verbose = FALSE) {
    # pp_data: list output from expand_person_period()
    # dat: original data list from generate_synthetic_icu()
    # gamma_val: value of gamma for the exponential tilting offset
    # n_iter: total number of MCMC iterations
    # burn_in: number of burn-in iterations to discard
    # tau2: prior variance for betas (covariate effects)
    # sigma_alpha2: prior variance for alphas (baseline hazard)
    # a0, b0: hyperparameters for sigma_b^2 ~ InvGamma(a0, b0)

    X.ast <- pp_data$X.ast
    y <- pp_data$y
    pid <- pp_data$patient_id
    R <- pp_data$R
    K <- pp_data$K
    P <- pp_data$P
    N <- dat$N
    dim_theta <- K + P

    offset <- make_offset_vector(pp_data, dat, gamma_val) # compute offset vector

    # priors
    mu_0 <- rep(0, dim_theta) # prior mean, all 0
    Sigma_0_inv <- diag( c( rep(1/sigma_alpha2, K), rep(1/tau2, P) ) )

    # hyperprior for sigma_b^2: Inverse-Gamma(a,b)
    a0 <- 0.01; b0 <- 0.01

    # initialize
    theta <- rep(0, dim_theta)      # fixed effects start at 0
    b <- rep(0, N)                  # frailty starts at 0
    sigma_b2 <- 1                   # initial value for sigma_b^2
    kappa <- y - 0.5                # kappa_i = y_i - 1/2

    # precompute patient membership: which rows belong to patient i
    patient_rows <- split(1:R, pid)

    # storage
    beta_samples <- matrix(data=0, nrow=n_iter - burn_in, ncol=P)

    for( iter in 1:n_iter){

        #  1. sample PG augmentation
        b_expanded <- b[pid] # expand b_i to row-level
        eta <- as.numeric(X.ast %*% theta + b_expanded + offset) # lin pred with frailty and offset
        omega <- rpg(n=R, h=1, z=eta) # omega_i ~ PG(1, eta_i)

        #  2. sample theta | omega, kappa, b ~ N(mu_theta, Sigma_theta)
        kappa_tilde <- kappa - omega * (b_expanded + offset) # adjust kappa for frailty and offset
        X0X <- crossprod(X.ast * omega, X.ast) # X'*Omega*X
        Sigma_theta_inv <- X0X + Sigma_0_inv
        Sigma_theta <- solve(Sigma_theta_inv) # posterior covariance
        mu_theta <- Sigma_theta %*% ( Sigma_0_inv %*% mu_0 + crossprod(X.ast, kappa_tilde) ) # posterior mean
        theta <- as.numeric( mu_theta + chol(Sigma_theta) %*% rnorm(dim_theta) )

        for (i in 1:N) {
            rows_i <- patient_rows[[i]]
            eta_fix_i <- as.numeric(X.ast[rows_i, ] %*% theta) # fixed part for patient i
            omega_i <- omega[rows_i]
            kappa_i <- kappa[rows_i]
            offset_i <- offset[rows_i]
            # posterior precision and mean for b_i
            prec_bi <- 1/sigma_b2 + sum(omega_i) # precision from prior + likelihood
            var_bi <- 1 / prec_bi
            mu_bi <- var_bi * sum(kappa_i - omega_i * (eta_fix_i + offset_i))
            b[i] <- rnorm(1, mean = mu_bi, sd = sqrt(var_bi))
        }

        a_post <- a0 + N / 2
        b_post <- b0 + sum((b)^2) / 2
        sigma_b2 <- 1 / rgamma(1, shape = a_post, rate = b_post)

        if (iter > burn_in) beta_samples[iter - burn_in, ] <- theta[(K+1):(K+P)]
    }

    colnames(beta_samples) <- colnames(pp_data$X.ast)[(K+1):(K+P)]
    beta_samples
}

# end of core.r
