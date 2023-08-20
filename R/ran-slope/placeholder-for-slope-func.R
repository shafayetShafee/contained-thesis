source(here::here("R/sim_utils.R"))

gen_two_level_slope_data <- function(m, n, sigma_mat, data_seed) {
  m = 10 # number of cluster
  n = 10 # size of each cluster
  N = m*n 
  
  set.seed(data_seed)
  
  x1c <- rnorm(N)
  x2 <- rnorm(N)
  x2b <- ifelse(x2 <= 0.5, 0, 1)

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
  return(multi_data)
}


sigma2_u0 <- 1
sigma2_u1 <- 2
sigma2_01 <- 0 # 0.75
sigma_mat <- matrix(c(sigma2_u0, sigma2_01, sigma2_01, sigma2_u1), 
                    byrow = TRUE, nrow = 2, ncol = 2)

m_data <- gen_two_level_slope_data(5, 10, sigma_mat, 1083)


multi_model <- mixed_model(fixed = Yij ~ X1c + X2b, random = ~ X1c | cluster, 
                           data = m_data, nAGQ = 11, family = binomial())

summary(multi_model)
multi_model$D

sigma_u_hat <- sym_mat_to_vec(multi_model$D) # var-cov matrix (estimate of random effects)
var_ran_effect <- var_slope_rand_effect(multi_model)
var_ran_effect

x1_val <- mean(m_data$X1c)
true_mor <- exp(sqrt(2*sigma2_u0 + 2*x1_val^2*sigma2_u1 + 4*x1_val*sigma2_01) * qnorm(0.75))
true_mor

mor_hat <- exp(
  sqrt(2*sigma_u_hat[1] + 2*x1_val^2*sigma_u_hat[3] + 4*x1_val*sigma_u_hat[2]) 
  * qnorm(0.75)
)
mor_hat

log_mor_expr <- function(x, x1_val) {
  sqrt((2*x[1]) + (2*x1_val^2*x[2]) + (4*x1_val*x[3])) * qnorm(0.75)
}

J <- numDeriv::jacobian(log_mor_expr, c(sigma_u_hat[1], sigma_u_hat[3], sigma_u_hat[2]), 
                        x1_val = x1_val) # success

var_log_mor <- as.numeric(J %*% var_ran_effect %*% t(J))
sqrt(var_log_mor)





