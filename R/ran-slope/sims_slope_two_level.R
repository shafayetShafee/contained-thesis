fixed_coeff <- c(-3, 1.75, 0.67)

sigma2_u0 <- 1
sigma2_u1 <- 2
sigma2_01 <- 0 # 0.75
sigma_mat <- matrix(c(sigma2_u0, sigma2_01, sigma2_01, sigma2_u1), 
                    byrow = TRUE, nrow = 2, ncol = 2)