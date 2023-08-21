library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(numDeriv)
library(GLMMadaptive)

source(here::here("R/sim_utils.R"))

gen_two_level_int_data <- function(m, n, fixed_coeff, sigma_u_sq, data_seed) {
  # m => number of cluster
  # n => size of each cluster
  # sigma_u_sq => variance of random effects
  N = m*n 
  
  set.seed(data_seed)
  
  x1c <- rnorm(N)
  x2 <- rnorm(N)
  x2b <- ifelse(x2 <= 0.5, 0, 1)
  
  if(length(fixed_coeff) != 3) {
    stop("fixed_coeff needs to be a vector of three element")
  }
  
  beta0 <- fixed_coeff[1]
  beta1 <- fixed_coeff[2]
  beta2 <- fixed_coeff[3]
  # 2 1.75 0.67
  
  uj <- rep(rnorm(m, 0, sqrt(sigma_u_sq)), each = n)
  eta_ij <- beta0 + beta1*x1c + beta2*x2b + uj
  pi_ij <- exp(eta_ij) / (1 + exp(eta_ij))
  
  yij <- rbinom(N, size = 1, prob = pi_ij)
  
  multi_data <- data.frame(
    cluster = rep(seq(m), each = n),
    id = rep(seq(n), times = m),
    uj = uj,
    X1c = x1c,
    X2b= x2b,
    Yij = yij
  )
  return(multi_data)
}


gen_two_level_slope_data <- function(m, n, fixed_coeff, sigma_mat, data_seed) {
  m = 10 # number of cluster
  n = 10 # size of each cluster
  N = m*n 
  
  set.seed(data_seed)
  
  x1c <- rnorm(N)
  x2 <- rnorm(N)
  x2b <- ifelse(x2 <= 0.5, 0, 1)
  
  u <- MASS::mvrnorm(m, mu = c(0, 0), Sigma = sigma_mat)
  u0j <- rep(u[, 1], each = n)
  u1j <- rep(u[, 2], each = n)
  
  beta0 <- fixed_coeff[1]
  beta1 <- fixed_coeff[2]
  beta2 <- fixed_coeff[3]
  
  eta_ij <- beta0 + u0j + beta1*x1c + u1j*x1c + beta2*x2b 
  pi_ij <- exp(eta_ij) / (1 + exp(eta_ij))
  
  yij <- rbinom(N, size = 1, prob = pi_ij)
  
  multi_data <- data.frame(
    cluster = rep(seq(m), each = n),
    id = rep(seq(n), times = m),
    u0j = u0j,
    u1j = u1j,
    X1c = x1c,
    X2b= x2b,
    Yij = yij
  )
  return(multi_data)
}



gen_int_mor_estimate <- function(m, n, fixed_coeff, sigma_u_sq, data_seed) {
  # data generation ----------------
  multi_data_int <- gen_two_level_int_data(m, n, fixed_coeff,
                                           sigma_u_sq, data_seed)
  
  # prevalence ---------------------
  prevalence <- mean(multi_data_int$Yij, na.rm = TRUE)
  
  # model fitting ------------------
  multi_model_int <- mixed_model(fixed = Yij ~ X1c + X2b, 
                                 random = ~ 1 | cluster,
                                 family = binomial("logit"), 
                                 nAGQ = 20, data = multi_data_int)
  
  # extracting value from model ----
  fixed_effect <- unname(multi_model_int$coefficients)
  beta0_hat <- fixed_effect[1]
  beta1_hat <- fixed_effect[2]
  beta2_hat <- fixed_effect[3]
  sigma_u_sq_hat <- multi_model_int$D[[1]] # rand effect
  var_sigma_u_sq_hat <- vcov_orig_scale(multi_model_int) # gets variance of rand effect
  
  # mor calc. ----------------------
  mor_hat <- exp(sqrt(2 * sigma_u_sq_hat) * qnorm(0.75))
  log_mor_hat <- log(mor_hat)
  

  # delta method -------------------
  log_mor_int_expr <- function(x) {
    # x => variance term
    sqrt(2 * x) * qnorm(0.75)
  }
  
  J <- numDeriv::jacobian(log_mor_int_expr, x = sigma_u_sq_hat)
  log_se_mor_hat <- as.numeric(sqrt(t(J) %*% var_sigma_u_sq_hat %*% J))
  se_mor_hat <- exp(log_se_mor_hat)
  
  # coverage calc ------------------
  ci <- log_mor_hat + c(-1, 1) * 1.96 * log_se_mor_hat
  ci_exp <- exp(ci)
  true_mor <- exp(sqrt(2 * sigma_u_sq) * qnorm(0.75))
  coverage <- as.numeric(ci_exp[1] <= true_mor && ci_exp[2] >= true_mor)
  
  is_model_converged <- multi_model_int$converged
  
  # creating output vector ---------
  out_vec <- c(true_mor = true_mor, mor_hat = mor_hat, 
               se_mor_hat = se_mor_hat,
               sigma_u_sq_hat = sigma_u_sq_hat,
               beta0_hat = beta0_hat, beta1_hat = beta1_hat, 
               beta2_hat = beta2_hat, coverage = coverage,
               prevalence = prevalence,
               converged = is_model_converged)
  
  if(is.na(is_model_converged) || !is_model_converged || mor_hat > 40 || se_mor_hat > 30) {
    out_vec_names <- names(out_vec)
    out_vec <- c(rep(NA, length(out_vec) - 1), FALSE)
    names(out_vec) <- out_vec_names
  }
  
  return(out_vec)
}


gen_slope_mor_estimate <- function(m, n, fixed_coeff, sigma_mat, data_seed) {
  
  # data generation ----------------
  multi_data_slope <- gen_two_level_slope_data(m, n, fixed_coeff,
                                               sigma_mat, data_seed)
  
  # sigma mat components -----------
  sigma_u1_sq <- sigma_mat[1, 1]
  sigma_u12_sq <- sigma_mat[1, 2]
  sigma_u2_sq <- sigma_mat[2, 2]
  
  # prevalence ---------------------
  prevalence <- mean(multi_data_slope$Yij, na.rm = TRUE)
  
  # model fitting ------------------
  multi_model_slope <- mixed_model(fixed = Yij ~ X1c + X2b, 
                                   random = ~ X1c | cluster, 
                                   data = multi_data_slope, 
                                   nAGQ = 11, family = binomial("logit"))
  
  # extracting value from model ----
  fixed_effect <- unname(multi_model_slope$coefficients)
  beta0_hat <- fixed_effect[1]
  beta1_hat <- fixed_effect[2]
  beta2_hat <- fixed_effect[3]
  sigma_u_sq_hat <- sym_mat_to_vec(multi_model_slope$D) # D_11 D_12 D_22
  sigma_u1_sq_hat <- sigma_u_sq_hat[1]
  sigma_u12_sq_hat <- sigma_u_sq_hat[2]
  sigma_u2_sq_hat <- sigma_u_sq_hat[3]
  
  # gets variance of rand effect
  # D_11 D_12 D_22 
  var_sigma_u_sq_hat <- var_slope_rand_effect(multi_model_slope) 
  
  
  # mor calc. ----------------------
  x1_val <- mean(multi_data_slope$X1c)
  mor_hat <- exp(
    sqrt(2*sigma_u_sq_hat[1] + 2*x1_val^2*sigma_u_sq_hat[3] + 4*x1_val*sigma_u_sq_hat[2]) 
    * qnorm(0.75)
  )
  log_mor_hat <- log(mor_hat)
  
  
  # delta method -------------------
  log_mor_expr <- function(x, x1_val) {
    sqrt((2*x[1]) + (2*x1_val^2*x[3]) + (4*x1_val*x[2])) * qnorm(0.75)
  }
  
  
  J <- numDeriv::jacobian(log_mor_expr, sigma_u_sq_hat, x1_val = x1_val)
  log_se_mor_hat <- sqrt(as.numeric(J %*% var_sigma_u_sq_hat %*% t(J)))
  se_mor_hat <- exp(log_se_mor_hat)
  
  # coverage calc ------------------
  ci <- log_mor_hat + c(-1, 1) * 1.96 * log_se_mor_hat
  ci_exp <- exp(ci)
  
  true_mor <- exp(
    sqrt(2*sigma_u1_sq + 2*x1_val^2*sigma_u2_sq + 4*x1_val*sigma_u12_sq) * qnorm(0.75)
  )
  coverage <- as.numeric(ci_exp[1] <= true_mor && ci_exp[2] >= true_mor)
  
  is_model_converged <- multi_model_slope$converged
  
  # creating output vector ---------
  out_vec <- c(true_mor = true_mor, mor_hat = mor_hat, 
               se_mor_hat = se_mor_hat,
               sigma_u1_sq_hat = sigma_u1_sq_hat,
               sigma_u12_sq_hat = sigma_u12_sq_hat,
               sigma_u2_sq_hat = sigma_u2_sq_hat,
               beta0_hat = beta0_hat, beta1_hat = beta1_hat, 
               beta2_hat = beta2_hat, coverage = coverage,
               prevalence = prevalence,
               converged = is_model_converged)
  # ***remainder `converged` needs to be the last column****
  
  if(is.na(is_model_converged) || !is_model_converged || mor_hat > 40 || se_mor_hat > 30) {
    out_vec_names <- names(out_vec)
    out_vec <- c(rep(NA, length(out_vec) - 1), FALSE)
    names(out_vec) <- out_vec_names
  }
  
  return(out_vec)
}


# gen_int_mor_estimate(50, 50, 2.5, 1083)


simulate_int <- function(m, n, fixed_coeff, sigma_u_sq, nsims = 1000, 
                         log_file, seed, ...) {
  
  # creating extra sims to get nsims after accounting 
  # for non-converged cases ---------------------
  total_sims = nsims + min(nsims, 500)
  conv_case_num = 0
  runs_required = 0

  # creating placeholder matrix for result --------
  out_mat <- matrix(NA, nrow = total_sims, ncol = 10)
  colnames(out_mat) <- c("true_mor", "mor_hat", "se_mor_hat", "sigma_u_sq_hat", "beta0_hat", 
                         "beta1_hat", "beta2_hat", "coverage", "prevalence",
                         "converged")
  
  # set the random state initiating seed
  set.seed(seed)
  iteration_seeds <- sample(x = 10:1000000, size = total_sims, replace = FALSE)
  # starting simulation ---------------------------
  for (i in 1:total_sims) {
    
    # checking for converged nsims ----------------
    if (conv_case_num == nsims) { 
      break
    } else {
      runs_required = runs_required + 1
    }
    
    data_seed <- iteration_seeds[i]
    model_output <- safe_and_quietly(fun = gen_int_mor_estimate,
                                     m = m, n = n, 
                                     fixed_coeff = fixed_coeff,
                                     sigma_u_sq = sigma_u_sq, 
                                     data_seed = data_seed)
    
    # extraction ---------------------------------
    model_result <- model_output$result
    model_convergence <- dplyr::coalesce(model_result["converged"], FALSE)
    model_messages <- model_output$messages
    model_warnings <- model_output$warnings
    model_error <- model_output$error
    
    # checking total number of converged cases ----
    if(model_convergence) {
      conv_case_num = conv_case_num + 1
    }
    
    # logging simulation run ----------------------
    log_sim_run(convergence = model_convergence, message = model_messages, 
                warning = model_warnings, error = model_error,
                data_seed = data_seed, log_file = log_file, iter_no = i)
    
    # account for error -------------------------
    msg_prob_detected <- is_valid_str(paste0(model_messages, collapse = ""))
    warning_detected <- is_valid_str(paste0(model_warnings, collapse = ""))
    error_detected <- is_valid_str(paste0(model_error, collapse = ""))

    if(msg_prob_detected || warning_detected || error_detected) {
      NA_result <- rep(NA, times = length(model_result))
      names(NA_result) <- names(model_result)
      model_result <- NA_result
    }

    # storing sim-run ---------------------------
    out_mat[i, ] <- model_result
  }
  
  return(
    list(
      out_mat = out_mat,
      runs_required = runs_required
    )
  )
}


simulate_slope <- function(m, n, fixed_coeff, sigma_mat, nsims = 1000, 
                           log_file, seed, ...) {
  
  # creating extra sims to get nsims after accounting 
  # for non-converged cases ---------------------
  total_sims = nsims + min(nsims, 500)
  conv_case_num = 0
  runs_required = 0
  
  # creating placeholder matrix for result --------
  out_mat <- matrix(NA, nrow = total_sims, ncol = 12)
  colnames(out_mat) <- c("true_mor", "mor_hat", "se_mor_hat", 
                         "sigma_u1_sq_hat", "sigma_u12_sq_hat", "sigma_u2_sq_hat",
                         "beta0_hat", "beta1_hat", "beta2_hat", 
                         "coverage", "prevalence", "converged")
  
  # set the random state initiating seed
  set.seed(seed)
  iteration_seeds <- sample(x = 10:1000000, size = total_sims, replace = FALSE)
  # starting simulation ---------------------------
  for (i in 1:total_sims) {
    
    # checking for converged nsims ----------------
    if (conv_case_num == nsims) { 
      break
    } else {
      runs_required = runs_required + 1
    }
    
    data_seed <- iteration_seeds[i]
    model_output <- safe_and_quietly(fun = gen_slope_mor_estimate,
                                     m = m, n = n, 
                                     fixed_coeff = fixed_coeff,
                                     sigma_mat = sigma_mat, 
                                     data_seed = data_seed)
    
    # extraction ---------------------------------
    model_result <- model_output$result
    model_convergence <- dplyr::coalesce(model_result["converged"], FALSE)
    model_messages <- model_output$messages
    model_warnings <- model_output$warnings
    model_error <- model_output$error
    
    # checking total number of converged cases ----
    if(model_convergence) {
      conv_case_num = conv_case_num + 1
    }
    
    # logging simulation run ----------------------
    log_sim_run(convergence = model_convergence, message = model_messages, 
                warning = model_warnings, error = model_error,
                data_seed = data_seed, log_file = log_file, iter_no = i)
    
    # account for error -------------------------
    msg_prob_detected <- is_valid_str(paste0(model_messages, collapse = ""))
    warning_detected <- is_valid_str(paste0(model_warnings, collapse = ""))
    error_detected <- is_valid_str(paste0(model_error, collapse = ""))
    
    if(msg_prob_detected || warning_detected || error_detected) {
      NA_result <- rep(NA, times = length(model_result))
      names(NA_result) <- names(model_result)
      model_result <- NA_result
    }
    
    # storing sim-run ---------------------------
    out_mat[i, ] <- model_result
  }
  
  return(
    list(
      out_mat = out_mat,
      runs_required = runs_required
    )
  )
}


run_simulations <- function(m, n, fixed_coeff, sigma_u_sq, nsims = 1000, 
                            simulation_type = c("int", "slope"),
                            seed = 1083, log_file, append = FALSE) {
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
                  int = simulate_int(m = m, n = n, 
                                     fixed_coeff = fixed_coeff,
                                     sigma_u_sq = sigma_u_sq,
                                     nsims = nsims, 
                                     log_file = log_file,
                                     seed = seed), 
                  slope = simulate_slope(m = m, n = n,
                                         fixed_coeff = fixed_coeff,
                                         sigma_mat = sigma_u_sq,
                                         nsims = nsims,
                                         log_file = log_file,
                                         seed = seed)
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
  hist_plot <- ggplot(tibble(log_mor_hat), aes(x = log_mor_hat)) +
    geom_histogram(bins = 30) +
    labs(x = "log(MOR)",
         title = cluster_info) +
    theme_classic()
  plot_path <- switch(sim_type,
    int = here::here("plots/ran-int"),
    slope = here::here("plots/ran-slope")
  )
  ggsave(paste0("hist_", m, "_", n, ".png"), path = plot_path)

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