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



sigma2_u0 <- 2.5
sigma2_u1 <- 2.1
sigma2_01 <- 0.75
sigma_mat <- matrix(c(sigma2_u0, sigma2_01, sigma2_01, sigma2_u1), 
                    byrow = TRUE, nrow = 2, ncol = 2)


m_data <- gen_two_level_slope_data(10, 10, sigma_mat, 1083)



