library(lme4)
library(merDeriv)
library(broom.mixed)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)


source(here::here("R/sim_utils.R"))


gen_two_level_data <- function(m, n, sigma_b2, seed) {
  # m => number of cluster
  # n => size of each cluster
  N = m*n 
  
  set.seed(seed)
  
  x1c <- rnorm(N)
  x2 <- rnorm(N)
  x2b <- ifelse(x2 <= 0.5, 0, 1)
  
  # sigma_b2 => variance of random effects
  beta0 <- 2
  beta1 <- 1.75
  beta2 <- 0.67
  
  uj <- rep(rnorm(m, 0, sqrt(sigma_b2)), each = n)
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

gen_mor_estimate <- function(m, n, sigma_b2, seed) {
  multi_data <- gen_two_level_data(m, n, sigma_b2, seed)
  multi_model <- glmer(Yij ~ X1c + X2b + (1 | cluster), 
                       family = binomial("logit"),
                       nAGQ = 20,
                       data = multi_data)
  
  output_df <- tidy(multi_model)
  fixed_effect <- output_df[output_df$effect == "fixed", ]$estimate
  beta0_hat <- fixed_effect[1]
  beta1_hat <- fixed_effect[2]
  beta2_hat <- fixed_effect[3]
  
  sigma_b2_hat <- as.numeric(output_df[output_df$effect == "ran_pars", "estimate"]) ^ 2
  
  MOR_hat <- exp(sqrt(2 * sigma_b2_hat) * qnorm(0.75))
  
  vv <- merDeriv::vcov.glmerMod(multi_model, full = TRUE, ranpar = "sd")
  se_sigma_b <- sqrt(diag(vv)[4])
  # se_mor_hat <- MOR_hat * sqrt(2) * qnorm(0.75) * se_sigma_b
  
  # log_mor_hat <- sqrt(2 * sigma_b2_hat) * qnorm(0.75)
  log_mor_hat <- log(MOR_hat)
  log_se_mor_hat <- sqrt(2) * qnorm(0.75) * se_sigma_b
  
  ci <- log_mor_hat + c(-1, 1) * 1.96 * log_se_mor_hat
  ci_exp <- exp(ci)
  
  true_MOR <- exp(sqrt(2 * sigma_b2) * qnorm(0.75))
  
  coverage <- as.numeric(ci_exp[1] <= true_MOR && ci_exp[2] >= true_MOR)
  
  se_mor_hat <- exp(sqrt(2) * qnorm(0.75) * se_sigma_b)
  
  out_vec <- c(mor_hat = MOR_hat, se_mor_hat = se_mor_hat,
               sigma_b2_hat = sigma_b2_hat,
               beta0_hat = beta0_hat, beta1_hat = beta1_hat, 
               beta2_hat = beta2_hat, coverage = coverage)
  
  if(MOR_hat > 50) {
    out_vec_names <- names(out_vec)
    out_vec <- rep(NA, length(out_vec))
    names(out_vec) <- out_vec_names
  }
  
  return(out_vec)
}


run_simulations <- function(m, n, sigma_b2, nsims = 1000, 
                            log_file, append = FALSE) {
  # browser()
  cat(paste("Simulation Log for", Sys.time(), "\n"),
      file = log_file, append = append)
  
  out_mat <- matrix(NA, nrow = nsims, ncol = 7)
  colnames(out_mat) <- c("mor_hat", "se_mor_hat", "sigma_b2_hat", "beta0_hat", 
                         "beta1_hat", "beta2_hat", "coverage")
  cluster_info <- paste0("Simulations for cluster size: ", n, " and #cluster: ", m, "\n")
  log_output(cluster_info, type = "Info", file = log_file)
  has_problem <- 0
  
  for (i in 1:nsims) {
    # browser() 
    seed <- floor(runif(1, 100, 1000000))
    model_output <- safe_and_quietly(fun = gen_mor_estimate,
                                     m = m, n = n, sigma_b2 = sigma_b2, 
                                     seed = seed)
    model_result <- model_output$result
    model_messages <- model_output$messages
    model_warnings <- model_output$warnings
    model_error <- model_output$error
    log_output(paste0(rep("-", 50), collapse = ""), type = "", file = log_file)
    log_output(seed, type = "Using seed", file = log_file)
    log_output(model_messages, type = "message", file = log_file)
    log_output(model_warnings, type = "warning", file = log_file)
    log_output(model_error, type = "error", file = log_file)
    
    catch_msg <- "boundary (singular) fit"
    msg_prob_detected <- is_valid(str_detect(model_messages, pattern = fixed(catch_msg)))
    
    warning_detected <- is_valid_str(paste0(model_warnings, collapse = ""))
    error_detected <- is_valid_str(paste0(model_error, collapse = ""))
    
    if(msg_prob_detected || warning_detected || error_detected) {
      NA_result <- rep(NA, times = length(model_result))
      names(NA_result) <- names(model_result)
      model_result <- NA_result
      has_problem <- has_problem + 1
    }
    
    out_mat[i, ] <- model_result
    log_output(paste0("Stored output for iteration ", i, "\n"), type="Info",
               file = log_file)
  }
  # print(out_mat)
  out_mat_means <- colMeans(out_mat, na.rm = TRUE)
  sim_se_mor_hat <- sd(out_mat[, "mor_hat"], na.rm = TRUE)
  runs_used = nrow(out_mat) - sum(apply(is.na(out_mat), 1, all))
  
  log_mor_hat <- log(out_mat[, "mor_hat"])
  hist_plot <- ggplot(tibble(log_mor_hat), aes(x = log_mor_hat)) +
    geom_histogram(bins = 30) +
    labs(x = "log(MOR)",
         title = cluster_info) +
    theme_classic()
  ggsave(paste0("hist_", m, "_", n, ".png"), path = here::here("plots/ran-int"))
  
  true_MOR <- exp(sqrt(2 * sigma_b2) * qnorm(0.75))
  relative_bias <- as.numeric(((out_mat_means["mor_hat"] - true_MOR) / true_MOR) * 100)
  
  return(c(cluster_number = m, cluster_size = n, out_mat_means, 
           sim_se_mor_hat = sim_se_mor_hat, 
           relative_bias = relative_bias,
           problem_perc = has_problem / nsims,
           runs_used = runs_used))
}
