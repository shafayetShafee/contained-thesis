library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(numDeriv)
library(GLMMadaptive)

source(here::here("R/sim_utils.R"))
source(here::here("R/two_level_int_sim_funcs.R"))
source(here::here("R/two_level_slope_sim_funcs.R"))


run_simulations <- function(m, n, fixed_coeff, sigma_u_sq, nsims = 1000, 
                            simulation_type = c("int", "slope"),
                            seed = 1083, log_file, append = FALSE,
                            plot_path, plot_name_suffix, 
                            more_iter) {
  # sigma_u_sq => a single element for (int)
  # sigma_u_sq => a 2x2 matrix for (slope)

  cat(
    paste("Simulation Log for", Sys.time(), "\n"), 
    file = log_file, append = append
  )
  
  # logging simulation param info ---------------
  cluster_info <- paste0("Simulations for cluster size: ", n, 
                         " and cluster number: ", m, "\n")
  log_output(cluster_info, type = "Info", file = log_file)
  
  # getting simulation matrix -------------------
  sim_type <- match.arg(simulation_type)
  sim_output <- switch(sim_type, 
                       int = simulate_two_lvl_int(m = m, n = n, 
                                          fixed_coeff = fixed_coeff,
                                          sigma_u_sq = sigma_u_sq,
                                          nsims = nsims, 
                                          log_file = log_file,
                                          seed = seed,
                                          more_iter = more_iter), 
                       slope = simulate_two_lvl_slope(m = m, n = n,
                                              fixed_coeff = fixed_coeff,
                                              sigma_mat = sigma_u_sq,
                                              nsims = nsims,
                                              log_file = log_file,
                                              seed = seed, 
                                              more_iter = more_iter)
  )
  out_mat <- sim_output$out_mat
  out_mat_means <- colMeans(out_mat, na.rm = TRUE)
  log_mor_hat <- log(out_mat[, "mor_hat"])
  sim_se_mor_hat <- exp(sd(log_mor_hat, na.rm = TRUE))
  # sim_se_mor_hat <- sd(out_mat[, "mor_hat"], na.rm = TRUE)
  runs_used = unname(sum(out_mat[, "converged"], na.rm = TRUE))
  runs_required = sim_output$runs_required
  
  # relative bias clac --------------------------
  true_mor <- out_mat_means["true_mor"]
  relative_bias <- as.numeric(((out_mat_means["mor_hat"] - true_mor) / true_mor) * 100)
  
  # generating histograms for checking ----------
  hist_plot <- ggplot(drop_na(tibble(log_mor_hat)), aes(x = log_mor_hat)) +
    geom_histogram(bins = 30) +
    labs(x = "log(MOR)",
         y = "frequency") +
    theme_classic()

  ggsave(paste0("hist_", m, "_", n, "_", plot_name_suffix, ".png"), path = plot_path)
  
  return(
    c(cluster_number = m, cluster_size = n, out_mat_means, 
      sim_se_mor_hat = sim_se_mor_hat, 
      relative_bias = relative_bias,
      runs_used = runs_used,
      runs_required = runs_required)
  )
}


# catch_msg <- "boundary (singular) fit"
# msg_prob_detected <- is_valid(str_detect(model_messages, pattern = fixed(catch_msg)))