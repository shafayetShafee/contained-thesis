source(here::here("R/run_simulations.R"))

# sigma_u_sq <- 2.5
# MOR <- exp(sqrt(2 * sigma_u_sq) * qnorm(0.75))
# MOR

# m = 10, 30, 50, 100
cluster_numbers <- c(10)

# n = 5, 10, 30, 50
cluster_size <- c(5)

cluster_params <- tidyr::expand_grid(cluster_size = cluster_size, 
                              cluster_numbers = cluster_numbers)

fixed_coeff <- c(2, 1.75, 0.67)
sigma_u_sq <- 2.5

log_file <- here::here("test.txt")

tictoc::tic()
res <- purrr::map2_dfr(.x = cluster_params$cluster_numbers, 
            .y = cluster_params$cluster_size, 
            .f = ~ run_simulations(m = .x, n = .y, 
                              fixed_coeff = fixed_coeff,
                              sigma_u_sq = sigma_u_sq, 
                              simulation_type = "int",
                              nsims = 10, seed = 1083,
                              log_file = log_file, append = TRUE)
            )

tictoc::toc()

final_res <- res
# final_res <- dplyr::bind_rows(final_res, res)

save(final_res, file=here::here("sim-results/rdata/ran-int/sim_res_20_aug.RData"))
saveRDS(final_res, file=here::here("sim-results/rds/ran-int/sim_res_20_aug.rds"))
