source(here::here("R/run_simulations.R"))

cluster_numbers <- c(10, 30, 50, 100)
cluster_size <- c(5, 10, 30, 50)

cluster_params <- tidyr::expand_grid(
  cluster_size = cluster_size,
  cluster_numbers = cluster_numbers
)

fixed_coeff <- c(-4.1, 1.75, 0.67)
sigma_u_sq <- 2.5

# moderate prevalence => -1.85 (30%)
# low prev => -4.1 (10%)

log_file <- here::here("log/log_two_lvl_int_low_prev.txt")
plot_path <- here::here("plots/two-lvl-ran-int/low-prev/")
plot_name_prefix <- "two_lvl_low_prev"

res <- purrr::map2_dfr(
  .x = cluster_params$cluster_numbers,
  .y = cluster_params$cluster_size,
  .f = ~ run_simulations(
    m = .x, n = .y,
    fixed_coeff = fixed_coeff,
    sigma_u_sq = sigma_u_sq,
    simulation_type = "two_lvl_int",
    nsims = 1000, seed = 1083,
    log_file = log_file, append = TRUE,
    plot_path = plot_path,
    plot_name_suffix = plot_name_prefix,
    more_iter = 1500
  )
)


# final_res_int_low_prev <- res
final_res_int_low_prev <- dplyr::bind_rows(final_res_int_low_prev, res)

save(final_res_int_low_prev,
  file = here::here("sim-results/rdata/sim_res_two_lvl_int_low_prev.RData")
)

saveRDS(final_res_int_low_prev,
  file = here::here("sim-results/rds/sim_res_two_lvl_int_low_prev.rds")
)
