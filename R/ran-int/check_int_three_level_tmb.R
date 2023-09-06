library(dplyr)
library(purrr)
library(glmmTMB)

source(here::here("R/sim_utils.R"))
source(here::here("R/three_level_int_sim_funcs.R"))

l = 20 # number of EA
m = 10 # number of HH
n = 15 # size of each HH
N = l*m*n 

set.seed(1083)

x1c <- rnorm(N)
x2 <- rnorm(N)
x2b <- ifelse(x2 <= 0.5, 0, 1)

sigma_u_sq <- 2.5
sigma_v_sq <- 1.5

beta0 <- -4.1
beta1 <- 1.75
beta2 <- 0.67

ujk <- rep(rnorm(l*m, 0, sqrt(sigma_u_sq)), each = n)
vk <- rep(rnorm(l, 0, sqrt(sigma_v_sq)), each = m*n)

eta_ijk <- beta0 + beta1*x1c + beta2*x2b + ujk + vk
pi_ijk <- exp(eta_ijk) / (1 + exp(eta_ijk))
yijk <- rbinom(N, size = 1, prob = pi_ijk)

multi_data_int <- data.frame(
  ea = rep(seq(l), each = m*n),
  hh = rep(rep(seq(m), each = n), times = l),
  id = rep(seq(n), times = l*m),
  ujk = ujk,
  vk = vk,
  X1c = x1c,
  X2b = x2b,
  Yijk = yijk
)


model = glmmTMB(Yijk ~ X1c + X2b + (1 | ea) + (1 | ea:hh),
                data = multi_data_int, family = "binomial")

# glmmTMB gives the var comp in log_sd scale

se_logsd <- sqrt(diag(vcov(model,full=TRUE)))[4:5] # sd(log_sd)
logsd <- model$sdr$par.fixed[4:5] # log_sd

sigma_sq_hat <- unname(exp(2*logsd))
# returns in order of (sigma_vk, sigma_ujk)

# by delta method
# se(sigma_sq_hat) = sqrt((2 * exp(2*logsd))^2) * se(logsd)
se_sigma_sq_hat <- diag(unname(se_logsd*2*exp(2*logsd))) 
# se of (sigma_vk, sigma_ujk)

var_sigma_sq_hat <- se_sigma_sq_hat^2

true_mor1 <- exp(sqrt(2 * sigma_u_sq) * qnorm(0.75))
mor1_hat <- exp(sqrt(2 * sigma_sq_hat[2]) * qnorm(0.75))
log_mor1_hat <- log(mor1_hat)

log_mor1_expr <- function(x) {
  # x is ran-effect parameterized as sd
  sqrt(2 * x) * qnorm(0.75)
}

J1 <- numDeriv::jacobian(log_mor1_expr, x = sigma_sq_hat[2])
log_se_mor1_hat <- as.numeric(sqrt(t(J1) %*% var_sigma_sq_hat[2, 2] %*% J1))
se_mor_hat <- exp(log_se_mor1_hat)


true_mor2 <- exp(sqrt(2 * (sigma_u_sq + sigma_v_sq)) * qnorm(0.75))
mor2_hat <- exp(sqrt(2 * sum(sigma_sq_hat)) * qnorm(0.75))
log_mor2_hat <- log(mor2_hat)

log_mor2_expr <- function(x) {
  sqrt(2 * sum(x)) * qnorm(0.75)
}

J2 <- numDeriv::jacobian(log_mor2_expr, x = sigma_sq_hat)
log_se_mor2_hat <- as.numeric(sqrt(J2 %*% var_sigma_sq_hat %*% t(J2)))
se_mor2_hat <- exp(log_se_mor2_hat)

ci_1 <- log_mor1_hat + c(-1, 1) * 1.96 * log_se_mor1_hat
ci_1_exp <- exp(ci_1)
coverage_1 <- as.numeric(ci_1_exp[1] <= true_mor1 && ci_1_exp[2] >= true_mor1)

ci_2 <- log_mor2_hat + c(-1, 1) * 1.96 * log_se_mor2_hat
ci_2_exp <- exp(ci_2)
coverage_2 <- as.numeric(ci_2_exp[1] <= true_mor2 && ci_2_exp[2] >= true_mor2)


