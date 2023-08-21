source(here::here("R/sim_utils.R"))

gen_two_level_slope_data <- function(m, n, fixed_coeff, sigma_mat, data_seed) {
  # m => number of cluster
  # n => size of each cluster
  N = m*n 
  
  set.seed(data_seed)
  
  x1c <- rnorm(N)
  x2 <- rnorm(N)
  x2b <- ifelse(x2 <= 0.5, 0, 1)

  u <- MASS::mvrnorm(m, mu = c(0, 0), Sigma = sigma_mat)
  u0j <- rep(u[, 1], each = n)
  u1j <- rep(u[, 2], each = n)
  
  beta0 <- fixed_coeff[1]
  beta1 <- fixed_coeff[2]
  beta2 <- fixed_coeff[3]
  
  eta_ij <- beta0 + u0j + beta1*x1c + u1j*x1c + beta2*x2b 
  pi_ij <- exp(eta_ij) / (1 + exp(eta_ij))
  
  yij <- rbinom(N, size = 1, prob = pi_ij)
  
  multi_data <- data.frame(
    cluster = rep(seq(m), each = n),
    id = rep(seq(n), times = m),
    u0j = u0j,
    u1j = u1j,
    X1c = x1c,
    X2b = x2b,
    Yij = yij
  )
  return(multi_data)
}


gen_slope_mor_estimate <- function(m, n, fixed_coeff, sigma_mat, data_seed) {

  # data generation ----------------
  multi_data_slope <- gen_two_level_slope_data(m, n, fixed_coeff,
                                           sigma_mat, data_seed)
  
  # sigma mat components -----------
  sigma_u1_sq <- sigma_mat[1, 1]
  sigma_u12_sq <- sigma_mat[1, 2]
  sigma_u2_sq <- sigma_mat[2, 2]
    
  # prevalence ---------------------
  prevalence <- mean(multi_data_slope$Yij, na.rm = TRUE)
  
  # model fitting ------------------
  multi_model_slope <- mixed_model(fixed = Yij ~ X1c + X2b, 
                                   random = ~ X1c | cluster, 
                                   data = multi_data_slope, 
                                   nAGQ = 11, family = binomial("logit"))
  
  # extracting value from model ----
  fixed_effect <- unname(multi_model_slope$coefficients)
  beta0_hat <- fixed_effect[1]
  beta1_hat <- fixed_effect[2]
  beta2_hat <- fixed_effect[3]
  sigma_u_sq_hat <- sym_mat_to_vec(multi_model_slope$D) # D_11 D_12 D_22
  sigma_u1_sq_hat <- sigma_u_sq_hat[1]
  sigma_u12_sq_hat <- sigma_u_sq_hat[2]
  sigma_u2_sq_hat <- sigma_u_sq_hat[3]
    
  # gets variance of rand effect
  # D_11 D_12 D_22 
  var_sigma_u_sq_hat <- var_slope_rand_effect(multi_model_slope) 
  
  
  # mor calc. ----------------------
  x1_val <- mean(multi_data_slope$X1c)
  mor_hat <- exp(
    sqrt(2*sigma_u_sq_hat[1] + 2*x1_val^2*sigma_u_sq_hat[3] + 4*x1_val*sigma_u_sq_hat[2]) 
    * qnorm(0.75)
  )
  log_mor_hat <- log(mor_hat)
  
  
  # delta method -------------------
  log_mor_expr <- function(x, x1_val) {
    sqrt((2*x[1]) + (2*x1_val^2*x[3]) + (4*x1_val*x[2])) * qnorm(0.75)
  }
  
  
  J <- numDeriv::jacobian(log_mor_expr, sigma_u_sq_hat, x1_val = x1_val)
  log_se_mor_hat <- sqrt(as.numeric(J %*% var_sigma_u_sq_hat %*% t(J)))
  se_mor_hat <- exp(log_se_mor_hat)
  
  # coverage calc ------------------
  ci <- log_mor_hat + c(-1, 1) * 1.96 * log_se_mor_hat
  ci_exp <- exp(ci)
  
  true_mor <- exp(
    sqrt(2*sigma_u1_sq + 2*x1_val^2*sigma_u2_sq + 4*x1_val*sigma_u12_sq) * qnorm(0.75)
    )
  coverage <- as.numeric(ci_exp[1] <= true_mor && ci_exp[2] >= true_mor)
  
  is_model_converged <- multi_model_slope$converged
  
  # creating output vector ---------
  out_vec <- c(mor_hat = mor_hat, se_mor_hat = se_mor_hat,
               sigma_u1_sq_hat = sigma_u1_sq_hat,
               sigma_u12_sq_hat = sigma_u12_sq_hat,
               sigma_u2_sq_hat = sigma_u2_sq_hat,
               beta0_hat = beta0_hat, beta1_hat = beta1_hat, 
               beta2_hat = beta2_hat, coverage = coverage,
               prevalence = prevalence,
               converged = is_model_converged)
  # ***remainder `converged` needs to be the last column****
  
  if(is.na(is_model_converged) || !is_model_converged || mor_hat > 40 || se_mor_hat > 30) {
    out_vec_names <- names(out_vec)
    out_vec <- c(rep(NA, length(out_vec) - 1), FALSE)
    names(out_vec) <- out_vec_names
  }
  
  return(out_vec)
}

# -------------------------------------


fixed_coeff <- c(2, 1.75, 0.67)

sigma2_u0 <- 1
sigma2_u1 <- 2
sigma2_01 <- 0 # 0.75
sigma_mat <- matrix(c(sigma2_u0, sigma2_01, sigma2_01, sigma2_u1), 
                    byrow = TRUE, nrow = 2, ncol = 2)




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

gen_slope_mor_estimate(10, 10, fixed_coeff = fixed_coeff, sigma_mat, 1083)



