## ------------------------------------------------------------------
##  Polson–Scott–Windle Polya‑Gamma Gibbs sampler for logistic regression
##  Re‑implementation of the Python example in pure R
#   https://gregorygundersen.com/blog/2019/09/20/polya-gamma/#a1-negative-binomial-likelihood
## ------------------------------------------------------------------
rm(list = ls())
# install.packages("BayesLogit") 
library(BayesLogit);library(MASS); library(scales); library(latex2exp)
set.seed(123)

expit <- function(x) ifelse(x > 0, 1 / (1 + exp(-x)), exp(x) / (1 + exp(x)))

# gen binomial data
gen.binom <- function(n, p=0.3){
    y <- rbinom(n=n, size=1, prob=p)
    X <- numeric(n)
    X[y==1] <- rnorm( n=sum(y==1), mean=0, sd=1) # component 1
    X[y==0] <- rnorm( n=sum(y==0), mean=4, sd=1.2) # component 2 
    return(list(X=X, y=y))
}

# place priors 
N.train <- 1000; N.test <- 1000
b       <- c(0, 0) # prior mean
B       <- diag(2) # precision

d.train     <- gen.binom(n=N.train, p=0.3); d.test <- gen.binom(n=N.test, p=0.3)

X.train     <- cbind( rep(1, N.train), d.train$X) # N x 2 (intercept + 1 covariate)
y.train     <- d.train$y
X.test      <- cbind( rep(1, N.test), d.test$X)
y.test      <- d.test$y


# Gibbs sampler
T           <- 100 
# w = (w_1, ..., w_n) ~ PG(n_i, x_i' beta) where n_i = # trials for obs. i (n_i=1 for binary data)
Omega.diag  <- rep(1, N.train)   # initialize omega
# S^{-1} = B + X' Omega X = B + X' diag(omega) X
Sigma.inv   <- solve(B + t(X.train) %*% diag(Omega.diag) %*% X.train) # precompute for beta update
beta.hat    <- mvrnorm(n=1, mu=b, Sigma=Sigma.inv) # initialize beta
kappa       <- y.train - 0.5          # kappa_i = y_i - n_i/2 = y_i - 1/2


for(i in seq_len(T)){ 

    # 3.a) sample w_i ~ PG(n_i, x_i' beta) for i=1,...,n
    # rpg(num = N, h = 1, z = X*beta )(BayesLogit) samples from PG(b, c) (b = n_i, c = x_i' beta)
    Omega.diag <- rpg(num = N.train, h = 1, z = X.train %*% beta.hat)   # sample omega_i

    # 3.b) sample beta | w ~ N(m, S) where S^{-1} = B + X' Omega X and m = S (B b + X' kappa)
    V <- solve( t(X.train) %*% diag(Omega.diag) %*% X.train + B )   # (post covariance)
    m <- V %*% (B %*% b + t(X.train) %*% kappa)                     # post mean
    beta.hat <- mvrnorm(n=1, mu=m, Sigma=V) # sample beta

}

# Predictions on test 
y.pred.prob <- expit(X.test %*% beta.hat)
y.pred      <- rbinom(n=N.test, size=1, prob=y.pred.prob)


# bins like np.linspace
bins <- seq(min(X.test[,2]) - 3,
            max(X.test[,2]) + 3,
            length.out = 100)

# use base R rgb() for transparency (alpha = 0.4)
hist(X.test[y.pred == 0, 2], breaks = bins,
     col   = rgb(1, 0, 0, 0.4),   # red, 40 % opacity
     main  = "Predicted y = 0 vs y = 1",
     xlab  = "Feature value (X[2])",
     ylab  = "Count",
     xlim  = range(bins))

hist(X.test[y.pred == 1, 2], breaks = bins,
     col   = rgb(0, 0, 1, 0.4),   # blue, 40 % opacity
     add   = TRUE)

legend("topright",
       legend = c("Predicted y = 0", "Predicted y = 1"),
       fill   = c(rgb(1,0,0,0.4), rgb(0,0,1,0.4)))

x.seq <- seq(min(X.test[,2]) - 3, max(X.test[,2]) + 3, length.out = 400)
lines(x.seq, expit(x.seq), col = "darkgreen", lwd = 2)

# add latex annotation for the logistic function
legend("topleft",
       legend = TeX("$P(y=1|x) = \\frac{1}{1 + e^{-x}}$"),
       bty = "n",  # no box around the legend
       cex = 1.2,  # increase text size
       text.col = "darkgreen")  
