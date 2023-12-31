library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(numDeriv)
library(GLMMadaptive)

source(here::here("R/sim_utils.R"))
source(here::here("R/two_level_int_sim_funcs.R"))
source(here::here("R/two_level_slope_sim_funcs.R"))
source(here::here("R/three_level_int_sim_funcs.R"))


run_simulations <- function(m, n, fixed_coeff, sigma_u_sq, nsims = 1000,
                            simulation_type = c("two_lvl_int", "two_lvl_slope", "three_lvl_int"),
                            seed = 1083, log_file, append = FALSE,
                            plot_path, plot_name_suffix,
                            more_iter, l = NULL) {
  # sigma_u_sq => a single element for (two_lvl_int)
  # sigma_u_sq => a 2x2 matrix for (two_lvl_slope)
  # sigma_u_sq => a vector of 2 element for three lvl int

  cat(
    paste("Simulation Log for", Sys.time(), "\n"),
    file = log_file, append = append
  )

  sim_type <- match.arg(simulation_type)

  if (sim_type == "three_lvl_int") {
    # logging simulation param info ---------------
    cluster_info <- paste0(
      "Simulations for ea no: ", l,
      ", hh no: ", m, " and hh size: ", n, "\n"
    )
  } else {
    cluster_info <- paste0(
      "Simulations for cluster size: ", n,
      " and cluster number: ", m, "\n"
    )
  }

  log_output(cluster_info, type = "Info", file = log_file)

  # getting simulation matrix -------------------

  sim_output <- switch(sim_type,
    two_lvl_int = simulate_two_lvl_int(
      m = m, n = n,
      fixed_coeff = fixed_coeff,
      sigma_u_sq = sigma_u_sq,
      nsims = nsims,
      log_file = log_file,
      seed = seed,
      more_iter = more_iter
    ),
    two_lvl_slope = simulate_two_lvl_slope(
      m = m, n = n,
      fixed_coeff = fixed_coeff,
      sigma_mat = sigma_u_sq,
      nsims = nsims,
      log_file = log_file,
      seed = seed,
      more_iter = more_iter
    ),
    three_lvl_int = simulate_three_lvl_int(
      l = l, m = m,
      n = n, fixed_coeff = fixed_coeff,
      sigma_sq = sigma_u_sq,
      nsims = nsims,
      log_file = log_file,
      seed = seed,
      more_iter = more_iter
    )
  )

  out_mat <- sim_output$out_mat
  out_mat_means <- colMeans(out_mat, na.rm = TRUE)
  runs_used <- unname(sum(out_mat[, "converged"], na.rm = TRUE))
  runs_required <- sim_output$runs_required

  if (sim_type == "two_lvl_int") {
    log_mor_hat <- log(out_mat[, "mor_hat"])
    sim_se_mor_hat <- exp(sd(log_mor_hat, na.rm = TRUE))

    # relative bias clac --------------------------
    true_mor <- out_mat_means["true_mor"]
    relative_bias <- as.numeric(((out_mat_means["mor_hat"] - true_mor) / true_mor) * 100)

    # generating histograms for checking ----------
    save_plot(m, n,
      log_mor_hat = log_mor_hat,
      plot_name_suffix = plot_name_suffix,
      plot_path = plot_path
    )

    return(c(
      cluster_number = m, cluster_size = n,
      out_mat_means,
      sim_se_mor_hat = sim_se_mor_hat,
      relative_bias = relative_bias,
      runs_used = runs_used,
      runs_required = runs_required
    ))
  } else if (sim_type == "two_lvl_slope") {
    log_mor_hat_q1 <- log(out_mat[, "mor_hat_q1"])
    log_mor_hat_q2 <- log(out_mat[, "mor_hat_q2"])
    log_mor_hat_q3 <- log(out_mat[, "mor_hat_q2"])
    sim_se_mor_hat_q1 <- exp(sd(log_mor_hat_q1, na.rm = TRUE))
    sim_se_mor_hat_q2 <- exp(sd(log_mor_hat_q2, na.rm = TRUE))
    sim_se_mor_hat_q3 <- exp(sd(log_mor_hat_q3, na.rm = TRUE))


    # generating histograms for checking ----------
    plot_name_suffix <- paste0(plot_name_suffix, "_q", 1:3)

    save_plot(m, n,
      log_mor_hat = log_mor_hat_q1,
      plot_name_suffix = plot_name_suffix[1],
      plot_path = plot_path
    )

    save_plot(m, n,
      log_mor_hat = log_mor_hat_q2,
      plot_name_suffix = plot_name_suffix[2],
      plot_path = plot_path
    )

    save_plot(m, n,
      log_mor_hat = log_mor_hat_q3,
      plot_name_suffix = plot_name_suffix[3],
      plot_path = plot_path
    )

    return(c(
      cluster_number = m, cluster_size = n,
      out_mat_means,
      sim_se_mor_hat_q1 = sim_se_mor_hat_q1,
      sim_se_mor_hat_q2 = sim_se_mor_hat_q2,
      sim_se_mor_hat_q3 = sim_se_mor_hat_q3,
      runs_used = runs_used,
      runs_required = runs_required
    ))
  } else if (sim_type == "three_lvl_int") {
    log_mor1_hat <- log(out_mat[, "mor1_hat"])
    log_mor2_hat <- log(out_mat[, "mor2_hat"])
    sim_se_mor1_hat <- exp(sd(log_mor1_hat, na.rm = TRUE))
    sim_se_mor2_hat <- exp(sd(log_mor2_hat, na.rm = TRUE))

    # relative bias clac --------------------------
    true_mor1 <- out_mat_means["true_mor1"]
    relative_bias_1 <- as.numeric(((out_mat_means["mor1_hat"] - true_mor1) / true_mor1) * 100)

    true_mor2 <- out_mat_means["true_mor2"]
    relative_bias_2 <- as.numeric(((out_mat_means["mor2_hat"] - true_mor2) / true_mor2) * 100)

    # generating histograms for checking ----------
    plot_name_suffix <- paste0(plot_name_suffix, "_mor", 1:2)

    save_plot(m, n,
      log_mor_hat = log_mor1_hat,
      plot_name_suffix = plot_name_suffix[1],
      plot_path = plot_path, l = l
    )

    save_plot(m, n,
      log_mor_hat = log_mor2_hat,
      plot_name_suffix = plot_name_suffix[2],
      plot_path = plot_path, l = l
    )

    return(c(
      ea_number = l,
      hh_number = m,
      hh_size = n,
      out_mat_means,
      sim_se_mor1_hat = sim_se_mor1_hat,
      sim_se_mor2_hat = sim_se_mor2_hat,
      relative_bias_1 = relative_bias_1,
      relative_bias_2 = relative_bias_2,
      runs_used = runs_used,
      runs_required = runs_required
    ))
  }
}
