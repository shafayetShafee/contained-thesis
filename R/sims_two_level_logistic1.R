source(here::here("R/sim_funcs.R"))

# sigma_b2 <- 2.5
# MOR <- exp(sqrt(2 * sigma_b2) * qnorm(0.75))
# MOR

# m = 10, 30, 50
cluster_numbers <- c(10, 30, 50)

# n = 10, 15, 30, 50
cluster_size <- c(30)

cluster_params <- expand_grid(cluster_size = cluster_size, 
                              cluster_numbers = cluster_numbers)
log_file <- here::here("log/log_aug_05.txt")

res1 <- map2_dfr(.x = cluster_params$cluster_numbers, 
            .y = cluster_params$cluster_size, 
            .f = ~ run_simulations(m = .x, n = .y, sigma_b2 = 2.5, nsims = 1000,
                                   log_file = log_file, append = TRUE)
            )

final_res <- dplyr::bind_rows(final_res, res1)


save(final_res, file="sim_res_05_aug.RData")
saveRDS(final_res, file="sim_res05_aug.rds")
