source(here::here("R/run_simulations.R"))

fixed_coeff <- c(-4.1, 1.75, 0.67)

sigma2_u1 <- 1 
sigma2_u2 <- 2
sigma2_u12 <- 0 # 0.75
sigma_mat <- matrix(c(sigma2_u1, sigma2_u12, sigma2_u12, sigma2_u2), 
                    byrow = TRUE, nrow = 2, ncol = 2)

# m = 10, 30, 50, 100
cluster_numbers <- c(50)

# n = 5, 10, 30, 50
cluster_size <- c(10)

cluster_params <- tidyr::expand_grid(cluster_size = cluster_size, 
                                     cluster_numbers = cluster_numbers)

log_file <- here::here("log/log_two_lvl_slp_low_prev.txt")
plot_path <- here::here("plots/two-lvl-ran-slope/low-prev/")
plot_name_prefix <- "two_lvl_slp_low_prev"

 