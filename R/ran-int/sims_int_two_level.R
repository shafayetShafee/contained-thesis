source(here::here("R/sim_funcs.R"))

# sigma_u_sq <- 2.5
# MOR <- exp(sqrt(2 * sigma_u_sq) * qnorm(0.75))
# MOR

# m = 10, 30, 50, 100
cluster_numbers <- c(50)

# n = 5, 10, 30, 50
cluster_size <- c(30)

cluster_params <- expand_grid(cluster_size = cluster_size, 
                              cluster_numbers = cluster_numbers)
log_file <- here::here("log/log_aug_18.txt")

res <- map2_dfr(.x = cluster_params$cluster_numbers, 
            .y = cluster_params$cluster_size, 
            .f = ~ run_simulations(m = .x, n = .y, sigma_u_sq = 2.5, 
                              nsims = 1000, seed = 1083,
                              log_file = log_file, append = TRUE)
            )

# final_res <- res
final_res <- dplyr::bind_rows(final_res, res)

save(final_res, file=here::here("sim-results/rdata/ran-int/sim_res_18_aug.RData"))
saveRDS(final_res, file=here::here("sim-results/rds/ran-int/sim_res18_aug.rds"))
