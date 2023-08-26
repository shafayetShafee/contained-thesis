source(here::here("R/run_simulations.R"))

fixed_coeff <- c(-2, 1.75, 0.67)

sigma2_u1 <- 1
sigma2_u2 <- 2
sigma2_u12 <- 0 # 0.75
sigma_mat <- matrix(c(sigma2_u1, sigma2_u12, sigma2_u12, sigma2_u2), 
                    byrow = TRUE, nrow = 2, ncol = 2)

# m = 10, 30, 50, 100
cluster_numbers <- c(100)

# n = 5, 10, 30, 50
cluster_size <- c(10)

cluster_params <- tidyr::expand_grid(cluster_size = cluster_size, 
                                     cluster_numbers = cluster_numbers)

log_file <- here::here("log/log_two_lvl_slp_low_prev.txt")
plot_path <- here::here("plots/two-lvl-ran-slope/low-prev/")
plot_name_prefix <- "two_lvl_slp_low_prev"

tictoc::tic()
res <- purrr::map2_dfr(.x = cluster_params$cluster_numbers, 
                       .y = cluster_params$cluster_size, 
                       .f = ~ run_simulations(m = .x, n = .y, 
                                              fixed_coeff = fixed_coeff,
                                              sigma_u_sq = sigma_mat, 
                                              simulation_type = "slope",
                                              nsims = 1000, seed = 1083,
                                              log_file = log_file, append = TRUE,
                                              plot_path = plot_path,
                                              plot_name_suffix = plot_name_prefix)
)

tictoc::toc()

# 12 hour 32 mins (high prev)
# c(2208.19, 2321.421, 2586.075, 3426.504, 1486.212, 1392.129, 1845.036, 
# 2961.062, )

# final_res_slp_low_prev <- res
final_res_slp_low_prev <- dplyr::bind_rows(final_res_slp_low_prev, res)

save(final_res_slp_low_prev, 
     file=here::here("sim-results/rdata/sim_res_two_lvl_slp_low_prev.RData"))
saveRDS(final_res_slp_low_prev, 
        file=here::here("sim-results/rds/sim_res_two_lvl_slp_low_prev.rds"))




# for high prev  ----------------------------------------------------------

# log_file <- here::here("log/log_two_lvl_slp_high_prev.txt")
# plot_path <- here::here("plots/two-lvl-ran-slope/high-prev/")
# plot_name_prefix <- "two_lvl_slp_high_prev"
# 
# tictoc::tic()
# res <- purrr::map2_dfr(.x = cluster_params$cluster_numbers, 
#                        .y = cluster_params$cluster_size, 
#                        .f = ~ run_simulations(m = .x, n = .y, 
#                                               fixed_coeff = fixed_coeff,
#                                               sigma_u_sq = sigma_mat, 
#                                               simulation_type = "slope",
#                                               nsims = 1000, seed = 1083,
#                                               log_file = log_file, append = TRUE,
#                                               plot_path = plot_path,
#                                               plot_name_suffix = plot_name_prefix)
# )
# 
# tictoc::toc()
# 
# # 12 hour 32 mins (high prev)
# 
# # final_res_slp_high_prev <- res
# final_res_slp_high_prev <- dplyr::bind_rows(final_res_slp_high_prev, res)
# 
# save(final_res_slp_high_prev, 
#      file=here::here("sim-results/rdata/sim_res_two_lvl_slp_high_prev.RData"))
# saveRDS(final_res_slp_high_prev, 
#         file=here::here("sim-results/rds/sim_res_two_lvl_slp_high_prev.rds"))

