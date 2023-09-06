source(here::here("R/run_simulations.R"))

# 20, 40
ea_number <- c(20)

# m = 10, 30
hh_number <- c(10)

# n = 5, 15, 30
hh_size <- c(5)

cluster_params <- tidyr::expand_grid(hh_number = hh_number, 
                                     hh_size = hh_size)

fixed_coeff <- c(-1.85, 1.75, 0.67)
sigma_sq <- c(2, 2.5)

# moderate prevalence => -1.85 (30%)
# low prev => -4.1 (10%)

log_file <- here::here("log/log_three_lvl_int_high_prev.txt")
plot_path <- here::here("plots/three-lvl-ran-int/high-prev/")
plot_name_prefix <- "three_lvl_high_prev"

tictoc::tic()
res <- purrr::map2_dfr(.x = cluster_params$hh_number, 
                       .y = cluster_params$hh_size, 
                       .f = ~ run_simulations(m = .x, n = .y, 
                                              fixed_coeff = fixed_coeff,
                                              sigma_u_sq = sigma_sq, 
                                              simulation_type = "three_lvl_int",
                                              nsims = 1000, seed = 1083,
                                              log_file = log_file, append = TRUE,
                                              plot_path = plot_path,
                                              plot_name_suffix = plot_name_prefix,
                                              more_iter = 500, 
                                              l = ea_number)
)

tictoc::toc()
beepr::beep(3)

final_res_int_high_prev <- res
# final_res_int_high_prev <- dplyr::bind_rows(final_res_int_high_prev, res)

save(final_res_int_high_prev, 
     file = here::here("sim-results/rdata/sim_res_three_lvl_int_high_prev.RData"))

saveRDS(final_res_int_high_prev, 
        file=here::here("sim-results/rds/sim_res_three_lvl_int_high_prev.rds"))
