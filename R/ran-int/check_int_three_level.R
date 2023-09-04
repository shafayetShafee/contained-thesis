library(lme4)
library(dplyr)
library(purrr)
library(numDeriv)
library(broom.mixed)

source(here::here("R/sim_utils.R"))
source(here::here("R/three_level_int_sim_funcs.R"))

l = 50 # number of EA
m = 30 # number of HH
n = 5 # size of each HH
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

t_data <- gen_three_level_int_data(l, m, n, c(-4.1, 1.75, 0.67), c(2.5, 1.5), 1083)

identical(multi_data_int, t_data)


model <- glmer(Yijk ~ X1c + X2b + (1 | ea) + (1 | ea:hh),
               data = multi_data_int, family = "binomial")

s = summary(model)

sigma_hat <- broom.mixed::tidy(model) %>% filter(effect == "ran_pars") %>% pull(estimate)
# returns in order of (sigma_ujk, sigma_vk)

sd_sigma_hat = diag(diag(vcov(VarCorr(model), model)))

true_mor1 <- exp(sqrt(2 * sigma_u_sq) * qnorm(0.75))

mor1_hat <- exp(sqrt(2 * sigma_hat[1]^2) * qnorm(0.75))

log_mor1_hat <- log(mor1_hat)

log_mor1_expr <- function(x) {
  # x is ran-effect parameterized as sd
  sqrt(2 * x^2) * qnorm(0.75)
}

J <- numDeriv::jacobian(log_mor1_expr, x = sigma_hat[1])
log_se_mor_hat <- as.numeric(sqrt(t(J) %*% sd_sigma_hat[1, 1] %*% J))
se_mor_hat <- exp(log_se_mor_hat)


true_mor2 <- exp(sqrt(2 * (sigma_u_sq + sigma_v_sq)) * qnorm(0.75))
 
mor2_hat <- exp(sqrt(2 * sum(sigma_hat^2)) * qnorm(0.75))

log_mor2_hat <- log(mor2_hat)

log_mor2_expr <- function(x) {
  sqrt(2 * sum(x^2)) * qnorm(0.75)
}

J2 <- numDeriv::jacobian(log_mor2_expr, x = sigma_hat)
log_se_mor2_hat <- as.numeric(sqrt(J2 %*% sd_sigma_hat %*% t(J2)))
se_mor2_hat <- exp(log_se_mor2_hat)

t_mor_info = est_three_lvl_int_mor(l, m, n, c(beta0, beta1, beta2), 
                                   sigma_sq = c(2.5, 1.5), 1083)

mor_info <- c(true_mor1, mor1_hat, se_mor_hat, true_mor2, mor2_hat, se_mor2_hat, 
              sigma_hat^2, unname(s$coefficients[, "Estimate"]))


sim_test = simulate_three_lvl_int(l, m, n, fixed_coeff = c(beta0, beta1, beta2), 
                                  sigma_sq = c(2.5, 1.5), nsims = 10, 
                                  log_file = "test.txt", seed = 1083, 
                                  more_iter = 1)
