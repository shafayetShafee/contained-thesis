m = 50 # number of cluster
n = 50 # size of each cluster
N = m*n 

set.seed(1083)

x1c <- rnorm(N)
x2 <- rnorm(N)
x2b <- ifelse(x2 <= 0.5, 0, 1)

sigma_u_sq <- 2.5
beta0 <- 2
beta1 <- 1.75
beta2 <- 0.67

uj <- rep(rnorm(m, 0, sqrt(sigma_u_sq)), each = n)
eta_ij <- beta0 + beta1*x1c + beta2*x2b + uj
pi_ij <- exp(eta_ij) / (1 + exp(eta_ij))

yij <- rbinom(N, size = 1, prob = pi_ij)

multi_data_int <- data.frame(
  cluster = rep(seq(m), each = n),
  id = rep(seq(n), times = m),
  uj = uj,
  X1c = x1c,
  X2b= x2b,
  Yij = yij
)

library(GLMMadaptive)

multi_model_int <- mixed_model(fixed = Yij ~ X1c + X2b, random = ~ 1 | cluster,
                              nAGQ = 20, family = binomial("logit"), 
                              data = multi_data_int)

summary(multi_model_int)
multi_model_int$coefficients[1]

sigma_u_sq_hat <- multi_model_int$D[[1]] # sigma_u squared

MOR <- exp(sqrt(2 * sigma_u_sq) * qnorm(0.75))
MOR_hat <- exp(sqrt(2 * sigma_u_sq_hat) * qnorm(0.75))

var_sigma_u_sq_hat <- vcov_orig_scale(multi_model_int)

log_mor_int_expr <- function(x) {
  # x => variance term
  sqrt(2 * x) * qnorm(0.75)
}

J <- numDeriv::jacobian(log_mor_int_expr, x = sigma_u_sq_hat)

# delta method 
log_se_mor_hat <- as.numeric(sqrt(t(J) %*% var_sigma_u_sq_hat %*% J))


# ------ prevalence

m_data = gen_two_level_int_data(10, 5, 2.5, sample(1:10000, 1, F))
mean(m_data$Yij)

# library(lme4)
# library(broom.mixed)
# 
# multi_model <- glmer(Yij ~ X1c + X2b + (1 | cluster),
#                      family = binomial("logit"), nAGQ = 20,
#                      data = multi_data)
# 
# summary(multi_model)
# output_df <- tidy(multi_model)
# 
# sigma_b2_hat <- as.numeric(output_df[output_df$effect == "ran_pars", "estimate"]) ^ 2

# library(merDeriv)
# vv <- merDeriv::vcov.glmerMod(multi_model, full = TRUE, ranpar = "var")
# diag(vv)[4] # variance of random effect estimate (i.e. Var(sigma_u_sq))