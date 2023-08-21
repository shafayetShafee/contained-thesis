source(here::here("R/run_simulations.R"))

fixed_coeff <- c(2, 1.75, 0.67)

sigma2_u0 <- 1
sigma2_u1 <- 2
sigma2_01 <- 0 # 0.75
sigma_mat <- matrix(c(sigma2_u0, sigma2_01, sigma2_01, sigma2_u1), 
                    byrow = TRUE, nrow = 2, ncol = 2)


cluster_numbers <- c(10)

# n = 5, 10, 30, 50
cluster_size <- c(10)

cluster_params <- tidyr::expand_grid(cluster_size = cluster_size, 
                                     cluster_numbers = cluster_numbers)

log_file <- here::here("test.txt")

tictoc::tic()
res <- purrr::map2_dfr(.x = cluster_params$cluster_numbers, 
                       .y = cluster_params$cluster_size, 
                       .f = ~ run_simulations(m = .x, n = .y, 
                                              fixed_coeff = fixed_coeff,
                                              sigma_u_sq = sigma_mat, 
                                              nsims = 10, seed = 1083,
                                              simulation_type = "slope",
                                              log_file = log_file, append = TRUE)
)

tictoc::toc()

final_res <- res
# final_res <- dplyr::bind_rows(final_res, res)

save(final_res, file=here::here("sim-results/rdata/ran-slope/sim_res_21_aug.RData"))
saveRDS(final_res, file=here::here("sim-results/rds/ran-slope/sim_res_21_aug.rds"))
