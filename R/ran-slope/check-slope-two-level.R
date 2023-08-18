m = 10 # number of cluster
n = 10 # size of each cluster
N = m*n 

set.seed(1083)

x1c <- rnorm(N)
x2 <- rnorm(N)
x2b <- ifelse(x2 <= 0.5, 0, 1)

sigma2_u0 <- 2.5
sigma2_u1 <- 2.1
sigma2_01 <- 0.75
sigma_mat <- matrix(c(sigma2_u0, sigma2_01, sigma2_01, sigma2_u1), 
                    byrow = TRUE, nrow = 2, ncol = 2)

u <- MASS::mvrnorm(m, mu = c(0, 0), Sigma = sigma_mat)
u0j <- rep(u[, 1], each = n)
u1j <- rep(u[, 2], each = n)

beta0 <- 2
beta1 <- 1.75
beta2 <- 0.67

eta_ij <- beta0 + u0j + beta1*x1c + u1j*x1c + beta2*x2b 
pi_ij <- exp(eta_ij) / (1 + exp(eta_ij))

yij <- rbinom(N, size = 1, prob = pi_ij)

multi_data <- data.frame(
  cluster = rep(seq(m), each = n),
  id = rep(seq(n), times = m),
  u0j = u0j,
  u1j = u1j,
  X1c = x1c,
  X2b= x2b,
  Yij = yij
)

library(GLMMadaptive)

multi_model <- mixed_model(fixed = Yij ~ X1c + X2b, random = ~ X1c | cluster, 
                           data = multi_data, nAGQ = 11, family = binomial())

summary(multi_model)
multi_model$D

sigma_u_hat <- sym_mat_to_vec(multi_model$D) # var-cov matrix (estimate of random effects)

rand_var <- vcov_orig_scale(multi_model)

v <- diag(rand_var)
v_mat <- diag(c(v[1], v[3], v[2]))

library(numDeriv)

log_mor_expr <- function(x) {
  sqrt((2*x[1]) + (2*0.5^2*x[2]) + (4*0.5*x[3])) * qnorm(0.75)
}

x1_val <- 0.5
true_mor <- exp(sqrt(2*sigma2_u0 + 2*x1_val^2*sigma2_u1 + 4*x1_val*sigma2_01) * qnorm(0.75))

mor_hat <- exp(
      sqrt(2*sigma_u_hat[1] + 2*x1_val^2*sigma_u_hat[3] + 4*x1_val*sigma_u_hat[2]) 
      * qnorm(0.75)
  )

J <- numDeriv::jacobian(log_mor_expr, c(sigma_u_hat[1], sigma_u_hat[3], sigma_u_hat[2])) # success

var_log_mor <- as.numeric(J %*% v_mat %*% t(J))



# attempt with lme4 -------------------------------------------------------

# library(lme4)
# library(broom.mixed)

# multi_model <- glmer(Yij ~ X1c + X2b + (X1c| cluster),
#                      family = binomial("logit"), nAGQ = 1,
#                      data = multi_data)
# vc_slope <- VarCorr(multi_model)
# print(vc_slope, comp=c("Variance","Std.Dev."), corr = F)
# 
# sigma_u_hat <- (output_df[output_df$effect == "ran_pars", ]$estimate)
#
# library(merDeriv)
#
# vv <- merDeriv::vcov.glmerMod(multi_model, full = TRUE, ranpar = "var")
# colnames(vv)
# 
# v <- diag(vv)[4:6]

