# ==============================================================================
# pg_demo.R
# Bayesian Discrete-Time Survival with Polya-Gamma Augmentation
# ==============================================================================
# demonstrates full modeling procesure from write-up, including:
# 1. synthetic ICU data generation
# 2. persion-period expansion
# 3. PG-augmented Gibbs sampling: standard model
# 4. PG-augmented Gibbs sampling: shared frailty model
# 5. sensitivity analysis
# 6. stan comparison
# 
# covariates that mirror data-prep.R "X" matrix:
#   age, gender, admission_type, comorbidity flags, GCS, HR, etc.
# ==============================================================================
rm(list = ls())
library(BayesLogit);library(coda);library(Matrix);library(ggplot2)
set.seed(2530)

# 1. synthetic data generation --------------------------------------------------
# generate N pateints with static covariates and discrete-time survival outcomes
# embed informative censoring: discharge probability depends on latent (b_i)

gen.synthetic.icu <- function(N = 500, K = 30, delta = 1, beta = rep(0.1, 9)){
    # N    : number of patients
    # K    : max discrete time periods (e.g. 30 days, ~4 weeks)
    # delta: interval widths in days (e.g. 1 = daily, 7 = weekly)

    # staic X_i
    age         <- rnorm(N, mean = 65, sd = 15)
    age         <- pmin( pmax(age, 18), 90) # truncate to [18,90] 
    male        <- rbinom(n=N, size=1, prob=0.5)
    emergency   <- rbinom(n=N, size=1, prob = 0.65)
    # icd codes (comorbidities)
    icd_sepsis  <- rbinom(n=N, size=1, prob = 0.25) # sepsis
    icd_heart_failure   <- rbinom(n=N, size=1, prob = 0.20) # heart failure
    icd_resp_failure    <- rbinom(n=N, size=1, prob = 0.20) # respiratory failure
    icu_aki     <- rbinom(n=N, size=1, prob = 0.15) # acute kidney injury
    gcs_impaired <- rbinom(n=N, size=1, prob = 0.20) # impaired consciousness (GCS <= 12)
    on_vasopressors <- rbinom(n=N, size=1, prob = 0.30) # on vasopressors at admission
    # standardize contuinous covariates (age)
    age_std <- (age - mean(age)) / sd(age)

    X <- cbind(
        age_std = age_std,
        male = male,
        emergency = emergency,
        icd_sepsis = icd_sepsis,
        icd_heart_failure = icd_heart_failure,
        icd_resp_failure = icd_resp_failure,
        icu_aki = icu_aki,
        gcs_impaired = gcs_impaired,
        on_vasopressors = on_vasopressors
    )
    P <- ncol(X)

    # ---- true parameters 
    # baseline hazard: smooth increasing then plateau
    alpha_true <- -3.5 + 0.8 * (1:K)  - 0.02 * (1:K)^2

    # covariate effects (log odds ratios)
    beta_true <- c(
        age_std = 0.25, # older increases risk
        male = 0.10, # slight higher risk
        emergency = 0.30, # emergency much worse 
        icd_sepsis = 0.60, # sepsis large effect
        icd_resp = 0.35, # respiratory moderate
        icd_hf = 0.40, # heart failure moderate
        icu_aki = 0.25, # acute kidney 
        gcs_impaired = 0.45, # impaired consciousness large effect
        vasopressor = 0.45 # on vasopressors large effect
    )

    # ----- latent frailty (b_i) 
    sigma_b_true <- 0.5 # frailty variance
    b_true <- rnorm(n=N, mean=0, sd=sigma_b_true)

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
        cens_logit <- -2.5 - 0.3 * log(1:K) - 0.8 * b_true[i] 
        # "0.8 * b_i" healthier patients (b_i < 0) have HIGHER censoring hazard
        # (more likely to be discharged early) T_i ~indep C_i | X_i
        
        cens_h <- plogis(cens_logit)

        for (t in 1:K){ 
            if (runif(1) < cens_h[t]){ # censored at time t
                C_latent[i] <- t
                break
            }
        }

    }

    # ----- observed data
    y_obs <- pmin(T_latent, C_latent) # observed time is min of event and censoring
    delta_i <- as.integer(T_latent <= C_latent) # event indicator: 1 if event observed, 0 if censored

    list(
        N = N,
        K = K,
        P = P,
        X = X,
        y_obs = y_obs,
        delta_i = delta_i,
        T_lantent = T_latent,
        C_latent = C_latent,
        b_true = b_true,
        alpha_true = alpha_true,
        beta_true = beta_true,
        sigma_b_true = sigma_b_true
    )
}





# choose discretization 
# delta = 1 : daily intervals, K = 30 (30 days)
# delta = 7 : weekly intervals, K = 4 (4 weeks)
DELTA <- 1
K_MAX <- ceiling(30 / DELTA)

dat <- gen.synthetic.icu(N = 500, K = K_MAX, delta = DELTA, beta = c(0.5, -0.3, rep(0,7)))

cat(sprintf("Discretization: delta = %d day(s), K = %d intervals\n", DELTA, K_MAX))
cat(sprintf("Generated %d patients, %d events (%.1f%% mortality)\n",
            dat$N, sum(dat$delta), 100 * mean(dat$delta)))
cat(sprintf("Median follow-up: %.1f intervals\n", median(dat$y_obs)))


# 2. person-period expansion --------------------------------------------------
# expand data to long format: one row per patient-period
# for i=1:N, t=1:y_obs[i], create row with:
#   - covariates X_i (static)
#   - time t (discrete time period)
#   - event indicator: 1 if t == y_obs[i] and delta_i == 1, else 0

expand.person.period <- function(dat){

    rows <- list(); idx <- 0 
    
    for (i in 1:dat$N){ 

        Ki <- dat$y_obs[i] # number of total rows for patient i

        for (t in 1:Ki){
            idx <- idx + 1
            y_it <- as.integer( t == Ki && dat$delta_i[i] == 1)
            rows[[idx]] <- c(
                i = i,
                t = t,
                y_it = y_it,
                dat$X[i, ]
            )
        }
    }

    pp <- as.data.frame(do.call(rbind, rows))

    # build design matrix X_ast = [D_t, X_i], D_t interval dummies
    K <- dat$K; P <- dat$P; R <- nrow(pp)

    D <- matrix(data= 0 , nrow=R, ncol=K) # interval dummies
    colnames(D) <- paste0("alpha_", 1:K)
    
    for (r in 1:R) D [r, pp$t[r]] <- 1 # set D_{r,t} = 1 if row r corresponds to interval t

    X_ast <- cbind( D, as.matrix(pp[, 4:ncol(pp)] ))

    list(
        pp = pp,
        X.ast = X_ast, 
        y = pp$y_it,
        patient_id = pp$i,
        interval = pp$t, 
        R = R, K = K, P = P, N = dat$N
    )

}


pp_dat <- expand.person.period(dat=dat)

cat(sprintf("Person-period dataset: %d rows (%d patients x avg %.1f intervals)\n",
            pp_dat$R, dat$N, pp_dat$R / dat$N))

# 3. PG-augmented Gibbs sampling: standard model --------------------------------

pg.gibbs.standard <- function( pp_data, n_iter = 2000, burn_in = 500){ 
    X.ast <- pp_data$X.ast; y <- pp_data$y; R <- pp_data$R; K <- pp_data$K; P <- pp_data$P
    # pid <- pp_data$patient_id
    dim_theta <- K + P # number of parameters: K alphas + P betas

    # priors
    mu_0 <- rep(0, dim_theta) # prior mean, all 0
    
    # block diag prior: precision (alpha, tau^2I)
    tau2 <- 10 # prior var. for betas
    sigma_alpha2 <- 1 # RW variance for alphas (baseline hazard)
    
    # simple use diag prior
    Sigma_0_inv <- diag( c( rep( 1/sigma_alpha2, K), rep(1/tau2, P) ) )

    # storage 
    theta_samples <- matrix(data=0, nrow=n_iter, ncol=dim_theta)

    # initialize 
    theta <- rep(0, dim_theta) # start at 0 (log-odds = 0 => p=0.5)
    kappa <- y - 0.5 # kappa_i = y_i - 1/2

    cat("Running standard PG Gibbs sampler...\n")
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)


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

        setTxtProgressBar(pb, iter)
    }
    close(pb)

    # discard burn-in
    keep <- (burn_in + 1):n_iter
    samples <- theta_samples[keep, ]
    colnames(samples) <- colnames(pp_data$X.ast)

    list(
        samples = samples,
        alpha = samples[, 1:K],
        beta = samples[, (K+1):(K+P) ]
    )
}


fit_standard <- pg.gibbs.standard(pp_dat, n_iter = 3000, burn_in = 1000)


# 4. PG-augmented Gibbs sampling: shared frailty model --------------------------------
pg.gibbs.frailty <- function( pp_data, n_iter = 2000, burn_in = 500){ 
    X.ast <- pp_data$X.ast; y <- pp_data$y; R <- pp_data$R; K <- pp_data$K; P <- pp_data$P; N <- pp_data$N
    patient_id <- pp_data$patient_id
    dim_theta <- K + P # number of fixed effects

    # priors
    mu_0 <- rep(0, dim_theta) # prior mean, all 0
    tau2 <- 10 # prior var. for betas
    sigma_alpha2 <- 1 # RW variance for alphas 
    Sigma_0_inv <- diag( c( rep(1/sigma_alpha2, K), rep(1/tau2, P) ) )

    # hyperprior for sigma_b^2: Inverse-Gamma(a,b)
    a0 <- 0.01; b0 <- 0.01

    # storage
    theta_samples <- matrix(data=0, nrow=n_iter, ncol=dim_theta)
    sigma_b_samples <- numeric(n_iter)
    b_samples <- matrix(data=0, nrow=n_iter, ncol=N)

    # initialize
    theta <- rep(0, dim_theta) # fixed effects start at 0
    b <- rep(0, N) # frailty starts at 0
    sigma_b2 <- 1 # initial value for sigma_b^2
    kappa <- y - 0.5 # kappa_i = y_i - 1/2

    # precompute patient membership: which rows belong to patient i
    patient_rows <- split(1:R, patient_id)

    cat("Running frailty PG Gibbs sampler...\n")
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)

    for( iter in 1:n_iter){
        
        #  1. sample PG augmentation
        b_expanded <- b[patient_id] # expand b_i to row-level
        eta <- as.numeric(X.ast %*% theta + b_expanded) # lin pred
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
        theta_samples[iter, ] <- theta
        sigma_b_samples[iter] <- sigma_b2
        b_samples[iter, ] <- b

        setTxtProgressBar(pb, iter)
    }
    close(pb)

    keep <- (burn_in + 1):n_iter
    samples <- theta_samples[keep, ]
    colnames(samples) <- colnames(pp_data$X.ast)

    list(
        samples = samples,
        alpha = samples[, 1:K],
        beta = samples[, (K+1):(K+P) ],
        sigma_b2 = sigma_b_samples[keep],
        b = b_samples[keep, ]
    )
}

fit_frailty <- pg.gibbs.frailty(pp_dat, n_iter = 3000, burn_in = 1000)




# 5. sensitivity analysis --------------------------------------------------
pg.gibbs.tilted <- function(pp_data, dat, gamma_val, n_iter = 2000, burn_in = 500){
    
    # same frailty sampler but with offset o_{i,t} = gamma * I(C_i = t, delta_i = 0) to adjust for informative censoring
    X.ast <- pp_data$X.ast; y <- pp_data$y; R <- pp_data$R; K <- pp_data$K; P <- pp_data$P; N <- pp_data$N
    patient_id <- pp_data$patient_id
    patient_rows <- split(1:R, patient_id) # list of row indices per patient
    # 1:R = row indices in pp_data, 
    # patient_id = vector of patient IDs for each row, length R
    dim_theta <- K + P # number of fixed effects
    
    # build offset vector: gamma at terminal interval for censored patients
    offset <- rep(0, R)

    for( i in 1:N){ 
        if (dat$delta[i] == 0){ # censored patient
            rows_i <- which(patient_id == i)
            terminal <- rows_i[length(rows_i)] # last row for patient i
            offset[terminal] <- gamma_val # add offset at term.
        }
    }

    # priors (same as frailty)
    mu_0 <- rep(0, dim_theta) # prior mean, all 0
    tau2 <- 10 # prior var. for betas
    sigma_alpha2 <- 1 # RW variance for alphas
    Sigma_0_inv <- diag( c( rep(1/sigma_alpha2, K), rep(1/tau2, P) ) ) 
    a0 <- 0.01; b0 <- 0.01

    # 
    theta <- rep(0, dim_theta)
    b <- rep(0, N)
    sigma_b2 <- 1
    kappa <- y - 0.5
    patient_rows <- lapply(1:N, function(i) which(patient_id == i))

    # storage
    beta_samples <- matrix(data=0, nrow=n_iter - burn_in, ncol=P)

    for (iter in 1:n_iter) {
        b_expanded <- b[patient_id] # expand b_i to row-level
        eta <- as.numeric(X.ast %*% theta) + b_expanded + offset
        omega <- rpg(R, h = 1, z = eta)

        kappa_tilde <- kappa - omega * (b_expanded + offset)
        X0X <- crossprod( X.ast * omega, X.ast)
        Sigma_theta <- solve(Sigma_0_inv + X0X)
        mu_theta <- Sigma_theta %*% (Sigma_0_inv %*% mu_0 + crossprod(X.ast, kappa_tilde))
        theta <- as.numeric(mu_theta + chol(Sigma_theta) %*% rnorm(dim_theta))

        eta_fixed <- as.numeric(X.ast %*% theta) + offset
        
    for (i in 1:N) {
        rows_i <- patient_rows[[i]]
        prec_bi <- 1/sigma_b2 + sum(omega[rows_i])
        var_bi  <- 1 / prec_bi
        mu_bi   <- var_bi * sum(kappa[rows_i] - omega[rows_i] * eta_fixed[rows_i])
        b[i] <- rnorm(1, mu_bi, sqrt(var_bi))
    }

    a_post <- a0 + N/2
    b_post <- b0 + sum(b^2)/2
    sigma_b2 <- 1/rgamma(1, shape = a_post, rate = b_post)

    if (iter > burn_in) {
      beta_samples[iter - burn_in, ] <- theta[(K+1):(K+P)]
    }
  }

  colnames(beta_samples) <- colnames(pp_data$X.ast)[(K+1):(K+P)]
  beta_samples
}

gamma_grid <- c(-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5) # test gamma vals 
sensitivity_results <- list()

for (g in gamma_grid) {
    cat(sprintf("  gamma = %.1f ... ", g))
    beta_g <- pg.gibbs.tilted(pp_dat, dat, gamma_val = g,
                                n_iter = 2000, burn_in = 500)
    sensitivity_results[[as.character(g)]] <- beta_g
    cat("done\n")
}


# 6. diagnostics ----------------------------------------------------------------

# posterior summaries
# beta comparison (standard vs frailty)
beta_names <- colnames(fit_standard$beta)
cat("\nCovariate effects (posterior mean [95% CI]):\n")
cat(sprintf("%-15s %8s %20s %20s\n", "Parameter", "True",
            "Standard", "Frailty"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (j in 1:length(beta_names)){ 
    true_val <- dat$beta_true[j]
    std_mean <- mean(fit_standard$beta[, j])
    std_ci <- quantile(fit_standard$beta[, j], probs = c(0.025, 0.975))
    fra_mean <- mean(fit_frailty$beta[, j])
    fra_ci <- quantile(fit_frailty$beta[, j], probs = c(0.025, 0.975))

  cat(sprintf("%-15s %8.3f   %6.3f [%5.2f,%5.2f]   %6.3f [%5.2f,%5.2f]\n",
              beta_names[j], true_val,
              std_mean, std_ci[1], std_ci[2],
              fra_mean, fra_ci[1], fra_ci[2]))
}


# EES diagnostics (coda)
cat("\nEffective Sample Size (frailty model, beta):\n")
for (j in 1:length(beta_names)){
    ess <- effectiveSize(fit_frailty$beta[, j])
    cat(sprintf("%-15s %8.1f\n", beta_names[j], ess))
}

cat(sprintf("  %-15s ESS = %.0f\n", "sigma_b2", effectiveSize(fit_frailty$sigma_b2)))

# frailty variance 
cat(sprintf("\nsigma_b^2: true = %.3f, posterior mean = %.3f [%.3f, %.3f]\n",
            dat$sigma_b_true^2,
            mean(fit_frailty$sigma_b2),
            quantile(fit_frailty$sigma_b2, 0.025),
            quantile(fit_frailty$sigma_b2, 0.975)))

# sensitivity: emergency admission across gamma
cat("\nSensitivity: 'emergency' coefficient across gamma:\n")
cat(sprintf("%-8s %10s %20s\n", "gamma", "mean", "95% CI"))

for (g in names(sensitivity_results)) {
    beta_emerg <- sensitivity_results[[g]][, "emergency"]
    cat(sprintf("%-8s %10.3f   [%.3f, %.3f]\n",
                g, mean(beta_emerg),
                quantile(beta_emerg, 0.025),
                quantile(beta_emerg, 0.975)))
}





# 7. plots -------------------------------------------------------------------------

# baseline hazard comparison
alpha_std_mean <- colMeans(fit_standard$alpha)
alpha_fra_mean <- colMeans(fit_frailty$alpha)
K <- dat$K
pdf("pg_demo_plots.pdf", width = 10, height = 8)
par(mfrow=c(2,2))


# Plot 1: Baseline hazard
plot(1:K, plogis(dat$alpha_true), type = "l", lwd = 2, col = "black",
     xlab = "Interval", ylab = "h(t)", main = "Baseline Hazard",
     ylim = c(0, max(plogis(alpha_std_mean)) * 1.3))
lines(1:K, plogis(alpha_std_mean), col = "blue", lwd = 2, lty = 2)
lines(1:K, plogis(alpha_fra_mean), col = "red", lwd = 2, lty = 3)
legend("topleft", c("True", "Standard", "Frailty"),
       col = c("black", "blue", "red"), lty = 1:3, lwd = 2, cex = 0.8)

# Plot 2: Frailty variance trace
plot(fit_frailty$sigma_b2, type = "l", col = "darkgrey",
     xlab = "Iteration", ylab = expression(sigma[b]^2),
     main = "Frailty Variance Trace")
abline(h = dat$sigma_b_true^2, col = "red", lwd = 2)
abline(h = mean(fit_frailty$sigma_b2), col = "blue", lwd = 2, lty = 2)

# Plot 3: Beta comparison
beta_std <- colMeans(fit_standard$beta)
beta_fra <- colMeans(fit_frailty$beta)
idx <- 1:length(beta_names)
plot(idx, dat$beta_true, pch = 16, cex = 1.5, col = "black",
     xlab = "", ylab = "Coefficient", main = "Covariate Effects",
     xaxt = "n", ylim = range(c(dat$beta_true, beta_std, beta_fra)) * 1.2)
points(idx - 0.1, beta_std, pch = 17, col = "blue", cex = 1.3)
points(idx + 0.1, beta_fra, pch = 15, col = "red", cex = 1.3)
axis(1, at = idx, labels = beta_names, las = 2, cex.axis = 0.7)
abline(h = 0, lty = 3)
legend("topleft", c("True", "Standard", "Frailty"),
       pch = c(16, 17, 15), col = c("black", "blue", "red"), cex = 0.8)

# Plot 4: Sensitivity analysis
emerg_means <- sapply(sensitivity_results, function(x) mean(x[, "emergency"]))
emerg_lo    <- sapply(sensitivity_results, function(x) quantile(x[, "emergency"], 0.025))
emerg_hi    <- sapply(sensitivity_results, function(x) quantile(x[, "emergency"], 0.975))
gamma_vals  <- as.numeric(names(sensitivity_results))

plot(gamma_vals, emerg_means, type = "b", pch = 16, lwd = 2,
     xlab = expression(gamma), ylab = "Emergency coefficient",
     main = "Sensitivity Analysis",
     ylim = range(c(emerg_lo, emerg_hi)))
arrows(gamma_vals, emerg_lo, gamma_vals, emerg_hi,
       angle = 90, code = 3, length = 0.05, col = "darkgrey")
abline(h = 0, lty = 3, col = "red")
abline(v = 0, lty = 3, col = "grey")

dev.off()
cat("\nPlots saved to pg_demo_plots.pdf\n")





# stan model code in pg_demo.stan 
if(requireNamespace("cmdstanr", quietly = TRUE)) {
    cat("\nFitting Stan model for comparison...\n")
    library(cmdstanr)

    mod <- cmdstan_model("pg_demo.stan")
    stan_data <- list(
        R = pp_dat$R,
        N = pp_dat$N,
        K = pp_dat$K,
        P = pp_dat$P,
        y = as.integer(pp_dat$y),
        X = pp_dat$X.ast[ , (dat$K + 1):(dat$K + dat$P) ], # covariate part of design matrix
        interval = as.integer(pp_dat$interval),
        patient_id = as.integer(pp_dat$patient_id)
    )

    fit_stan <- mod$sample(
        data = stan_data,
        chains = 2, 
        parallel_chains = 2,
        iter_warmup = 500,
        iter_sampling = 1000,
        refresh = 500
    )

    cat("\n=== Stan vs PG Gibbs: Beta Comparison ===\n")
    stan_summary <- fit_stan$summary(variables = paste0("beta[", 1:dat$P, "]"))
    cat(sprintf("%-15s %8s %12s %12s\n", "Parameter", "True", "PG Gibbs", "Stan"))
    cat(paste(rep("-", 50), collapse = ""), "\n")
    for (j in 1:dat$P) {
        cat(sprintf("%-15s %8.3f %12.3f %12.3f\n",
                    beta_names[j], dat$beta_true[j],
                    mean(fit_frailty$beta[, j]),
                    stan_summary$mean[j]))
    }

    cat(sprintf("\nsigma_b: PG = %.3f, Stan = %.3f (true = %.3f)\n",
                sqrt(mean(fit_frailty$sigma_b2)),
                fit_stan$summary("sigma_b")$mean,
                dat$sigma_b_true))

    } else {
    cat("\ncmdstanr not available. Install with:\n")
    cat("  install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))\n")
    cat("  cmdstanr::install_cmdstan()\n")
    cat("Then rerun this section.\n")
}















