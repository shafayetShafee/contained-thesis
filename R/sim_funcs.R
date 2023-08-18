library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(merDeriv)
library(GLMMadaptive)

source(here::here("R/sim_utils.R"))

gen_two_level_int_data <- function(m, n, sigma_u_sq, seed) {
  # m => number of cluster
  # n => size of each cluster
  N = m*n 
  
  set.seed(seed)
  
  x1c <- rnorm(N)
  x2 <- rnorm(N)
  x2b <- ifelse(x2 <= 0.5, 0, 1)
  
  # sigma_u_sq => variance of random effects
  beta0 <- 2
  beta1 <- 1.75
  beta2 <- 0.67
  
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

gen_int_mor_estimate <- function(m, n, sigma_u_sq, seed) {
  # data generation ----------------
  multi_data_int <- gen_two_level_int_data(m, n, sigma_u_sq, seed)
  
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
  out_vec <- c(mor_hat = mor_hat, se_mor_hat = se_mor_hat,
               sigma_u_sq_hat = sigma_u_sq_hat,
               beta0_hat = beta0_hat, beta1_hat = beta1_hat, 
               beta2_hat = beta2_hat, coverage = coverage,
               converged = is_model_converged)
  
  if(!is_model_converged) {
    out_vec_names <- names(out_vec)
    out_vec <- c(rep(NA, length(out_vec) - 1), is_model_converged)
    names(out_vec) <- out_vec_names
  }
  
  return(out_vec)
}

# gen_int_mor_estimate(50, 50, 2.5, 1083)


run_simulations <- function(m, n, sigma_u_sq, nsims = 1000, 
                            log_file, append = FALSE) {
  # browser()
  cat(
    paste("Simulation Log for", Sys.time(), "\n"), 
    file = log_file, append = append
  )
  
  out_mat <- matrix(NA, nrow = nsims, ncol = 8)
  colnames(out_mat) <- c("mor_hat", "se_mor_hat", "sigma_u_sq_hat", "beta0_hat", 
                         "beta1_hat", "beta2_hat", "coverage", "converged")
  cluster_info <- paste0("Simulations for cluster size: ", n, " and #cluster: ", m, "\n")
  log_output(cluster_info, type = "Info", file = log_file)
  # has_problem <- 0
  
  for (i in 1:nsims) {
    # browser() 
    seed <- floor(runif(1, 100, 1000000))
    model_output <- safe_and_quietly(fun = gen_int_mor_estimate,
                                     m = m, n = n, sigma_u_sq = sigma_u_sq, 
                                     seed = seed)
    model_result <- model_output$result
    model_convergence <- model_result["converged"]
    model_messages <- model_output$messages
    model_warnings <- model_output$warnings
    model_error <- model_output$error
    
    out_mat[i, ] <- model_result
    
    log_output(paste0(rep("-", 50), collapse = ""), type = "", file = log_file)
    log_output(seed, type = "Using seed", file = log_file)
    log_output(model_convergence, type = "converged: ", file = log_file)
    log_output(model_messages, type = "message", file = log_file)
    log_output(model_warnings, type = "warning", file = log_file)
    log_output(model_error, type = "error", file = log_file)
    log_output(
      paste0("Stored output for iteration ", i, "\n"), 
      type="Info", file = log_file
    )
    
    # catch_msg <- "boundary (singular) fit"
    # msg_prob_detected <- is_valid(str_detect(model_messages, pattern = fixed(catch_msg)))
    # 
    # warning_detected <- is_valid_str(paste0(model_warnings, collapse = ""))
    # error_detected <- is_valid_str(paste0(model_error, collapse = ""))
    # 
    # if(msg_prob_detected || warning_detected || error_detected) {
    #   NA_result <- rep(NA, times = length(model_result))
    #   names(NA_result) <- names(model_result)
    #   model_result <- NA_result
    #   has_problem <- has_problem + 1
    # }
  }
  # print(out_mat)
  out_mat_means <- colMeans(out_mat, na.rm = TRUE)
  sim_se_mor_hat <- sd(out_mat[, "mor_hat"], na.rm = TRUE)
  runs_used = out_mat_means["converged"] * 1000
  
  log_mor_hat <- log(out_mat[, "mor_hat"])
  hist_plot <- ggplot(tibble(log_mor_hat), aes(x = log_mor_hat)) +
    geom_histogram(bins = 30) +
    labs(x = "log(MOR)",
         title = cluster_info) +
    theme_classic()
  ggsave(paste0("hist_", m, "_", n, ".png"), path = here::here("plots/ran-int"))
  
  true_mor <- exp(sqrt(2 * sigma_u_sq) * qnorm(0.75))
  relative_bias <- as.numeric(((out_mat_means["mor_hat"] - true_mor) / true_mor) * 100)
  
  return(c(cluster_number = m, cluster_size = n, out_mat_means, 
           sim_se_mor_hat = sim_se_mor_hat, 
           relative_bias = relative_bias,
           runs_used = runs_used))
}
