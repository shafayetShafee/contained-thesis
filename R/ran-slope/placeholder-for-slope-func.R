source(here::here("R/sim_utils.R"))

# -------------------------------------


fixed_coeff <- c(2, 1.75, 0.67)

sigma2_u0 <- 1
sigma2_u1 <- 2
sigma2_01 <- 0 # 0.75
sigma_mat <- matrix(c(sigma2_u0, sigma2_01, sigma2_01, sigma2_u1), 
                    byrow = TRUE, nrow = 2, ncol = 2)

# test <- simulate_two_lvl_slope(10, 5, fixed_coeff, sigma_mat, nsims = 10, "test.txt", 1083)

m_data <- gen_two_level_slope_data(10, 10, fixed_coeff = fixed_coeff,
                                   sigma_mat, 1083)

tictoc::tic()
multi_model <- mixed_model(fixed = Yij ~ X1c + X2b, random = ~ X1c | cluster, 
                           data = m_data, nAGQ = 11, family = binomial())

tictoc::toc()

summary(multi_model)
multi_model$D

sigma_u_hat <- sym_mat_to_vec(multi_model$D) # var-cov matrix (estimate of random effects)
var_sigma_u_sq_hat <- var_slope_rand_effect(multi_model)
var_sigma_u_sq_hat

x1_val <- mean(m_data$X1c)
true_mor <- exp(sqrt(2*sigma2_u0 + 2*x1_val^2*sigma2_u1 + 4*x1_val*sigma2_01) * qnorm(0.75))
true_mor

mor_hat <- exp(
  sqrt(2*sigma_u_hat[1] + 2*x1_val^2*sigma_u_hat[3] + 4*x1_val*sigma_u_hat[2]) 
  * qnorm(0.75)
)
mor_hat

log_mor_expr <- function(x, x1_val) {
  sqrt((2*x[1]) + (2*x1_val^2*x[3]) + (4*x1_val*x[2])) * qnorm(0.75)
}

J <- numDeriv::jacobian(log_mor_expr, sigma_u_hat, 
                        x1_val = x1_val) # success

log_se_mor_hat <- sqrt(as.numeric(J %*% var_sigma_u_sq_hat %*% t(J)))
exp(log_se_mor_hat)
