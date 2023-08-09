m = 100 # number of cluster
n = 50 # size of each cluster
N = m*n 

set.seed(1083)

x1c <- rnorm(N)
x2 <- rnorm(N)
x2b <- ifelse(x2 <= 0.5, 0, 1)

sigma_u0 <- 2.5
sigma_u1 <- 1.75
beta0 <- 2
beta1 <- 1.75
beta2 <- 0.67

u0j <- rep(rnorm(m, 0, sqrt(sigma_u0)), each = n)
u1j <- rep(rnorm(m, 0, sqrt(sigma_u1)), each = n)
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

library(lme4)
library(broom.mixed)

multi_model <- glmer(Yij ~ X1c + X2b + (X1c| cluster), 
                     family = binomial("logit"), 
                     data = multi_data)

summary(multi_model)
output_df <- tidy(multi_model)

vc_slope <- VarCorr(multi_model)

print(vc_slope, comp=c("Variance","Std.Dev."), corr = F)

sigma_u_hat <- as.numeric(output_df[output_df$effect == "ran_pars", ] |> dplyr::pull(estimate)) ^ 2

library(merDeriv)


vv <- merDeriv::vcov.glmerMod(multi_model, full = TRUE, ranpar = "var")
colnames(vv)

v <- diag(vv)[4:6]

v_mat <- diag(c(v[1], v[3], v[2]))

library(numDeriv)

log_mor_expr <- function(x) {
  sqrt((2*x[1]) + (2*0.5^2*x[2]) + (4*0.5*x[3])) * qnorm(0.75)
}

J <- jacobian(log_mor_expr, c(2.58, 1.80, -0.164)) # success

var_log_mor <- J %*% v_mat %*% t(J)  

# MOR <- exp(sqrt(2 * sigma_b2) * qnorm(0.75))
# MOR_hat <- exp(sqrt(2 * sigma_b2_hat) * qnorm(0.75))
