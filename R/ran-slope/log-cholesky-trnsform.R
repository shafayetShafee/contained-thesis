# glmmadaptive reprex -----------------------------------------------------

library(GLMMadaptive)

set.seed(1234)
n <- 100 # number of subjects
K <- 8 # number of measurements per subject
t_max <- 15 # maximum follow-up time

# we construct a data frame with the design: 
# everyone has a baseline measurement, and then measurements at random follow-up times
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects
X <- model.matrix(~ sex * time, data = DF)
Z <- model.matrix(~ time, data = DF)

betas <- c(-2.13, -0.25, 0.24, -0.05) # fixed effects coefficients
D11 <- 0.48 # variance of random intercepts
D22 <- 0.1 # variance of random slopes

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
# linear predictor
eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
# we simulate binary longitudinal data
DF$y <- rbinom(n * K, 1, plogis(eta_y))

fm <- mixed_model(fixed = y ~ sex * time, random = ~ time | id, data = DF,
                  family = binomial())


D <- vcov(fm, parm = "var-cov")
D

# estimated covariance matrix of random effects
D <- fm$D

# transform from covariance matrix to entries of cholesky factor with 
# log-transformed main diagonal
D_chol_entries <- GLMMadaptive:::chol_transf(D)

D_chol_to_D <- function(x) {
  
  # transform from entries of cholesky factor with log-transformed main diagonal
  # to covariance matrix
  D <- GLMMadaptive:::chol_transf(x)
  
  D[upper.tri(D, diag = TRUE)]
}

J <- numDeriv::jacobian(D_chol_to_D, D_chol_entries)

# estimated covariance matrix of D_chol_entries
V_chol <- vcov(fm, parm = "var-cov")

# estimated covariance matrix of entries of D
V <- J %*% V_chol %*% t(J)

se <- sqrt(diag(V))

cat("--- D_11 ---\nEstimate:", D[1, 1], "\nStd. Error:", se[1], fill = TRUE)
# --- D_11 ---
# Estimate: 0.4639974 
# Std. Error: 0.5992368
cat("--- D_22 ---\nEstimate:", D[2, 2], "\nStd. Error:", se[3], fill = TRUE)
# --- D_22 ---
# Estimate: 0.06304059 
# Std. Error: 0.02843676

