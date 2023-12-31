gen_two_level_int_data <- function(m, n, fixed_coeff, sigma_u_sq, data_seed) {
  # m => number of cluster
  # n => size of each cluster
  # sigma_u_sq => variance of random effects
  N <- m * n

  set.seed(data_seed)

  x1c <- rnorm(N)
  x2 <- rnorm(N)
  x2b <- ifelse(x2 <= 0.5, 0, 1)

  if (length(fixed_coeff) != 3) {
    stop("fixed_coeff needs to be a vector of three element")
  }

  beta0 <- fixed_coeff[1]
  beta1 <- fixed_coeff[2]
  beta2 <- fixed_coeff[3]
  # 2 1.75 0.67

  uj <- rep(rnorm(m, 0, sqrt(sigma_u_sq)), each = n)
  eta_ij <- beta0 + beta1 * x1c + beta2 * x2b + uj
  pi_ij <- exp(eta_ij) / (1 + exp(eta_ij))

  yij <- rbinom(N, size = 1, prob = pi_ij)

  multi_data <- data.frame(
    cluster = rep(seq(m), each = n),
    id = rep(seq(n), times = m),
    uj = uj,
    X1c = x1c,
    X2b = x2b,
    Yij = yij
  )
  return(multi_data)
}


est_two_lvl_int_mor <- function(m, n, fixed_coeff, sigma_u_sq, data_seed) {
  # data generation ----------------
  multi_data_int <- gen_two_level_int_data(
    m, n, fixed_coeff,
    sigma_u_sq, data_seed
  )

  # prevalence ---------------------
  prevalence <- mean(multi_data_int$Yij, na.rm = TRUE)

  # model fitting ------------------
  multi_model_int <- mixed_model(
    fixed = Yij ~ X1c + X2b,
    random = ~ 1 | cluster,
    family = binomial("logit"),
    nAGQ = 20, data = multi_data_int
  )

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
  out_vec <- c(
    true_mor = true_mor, mor_hat = mor_hat,
    se_mor_hat = se_mor_hat,
    sigma_u_sq_hat = sigma_u_sq_hat,
    beta0_hat = beta0_hat, beta1_hat = beta1_hat,
    beta2_hat = beta2_hat, coverage = coverage,
    prevalence = prevalence,
    converged = is_model_converged
  )

  if (is.na(is_model_converged) || !is_model_converged || mor_hat > 40 || se_mor_hat > 30) {
    out_vec_names <- names(out_vec)
    out_vec <- c(rep(NA, length(out_vec) - 1), FALSE)
    names(out_vec) <- out_vec_names
  }

  return(out_vec)
}



simulate_two_lvl_int <- function(m, n, fixed_coeff, sigma_u_sq, nsims = 1000,
                                 log_file, seed, more_iter, ...) {
  # creating extra sims to get nsims after accounting
  # for non-converged cases ---------------------
  total_sims <- nsims + more_iter
  conv_case_num <- 0
  runs_required <- 0

  # creating placeholder matrix for result --------
  out_mat <- matrix(NA, nrow = total_sims, ncol = 10)
  out_colnames <- c(
    "true_mor", "mor_hat", "se_mor_hat", "sigma_u_sq_hat",
    "beta0_hat", "beta1_hat", "beta2_hat", "coverage",
    "prevalence", "converged"
  )
  colnames(out_mat) <- out_colnames

  # set the random state initiating seed
  set.seed(seed)
  iteration_seeds <- sample(x = 10:1000000, size = total_sims, replace = FALSE)
  # starting simulation ---------------------------
  for (i in 1:total_sims) {
    # checking for converged nsims ----------------
    if (conv_case_num == nsims) {
      break
    } else {
      runs_required <- runs_required + 1
    }

    data_seed <- iteration_seeds[i]
    model_output <- safe_and_quietly(
      fun = est_two_lvl_int_mor,
      m = m, n = n,
      fixed_coeff = fixed_coeff,
      sigma_u_sq = sigma_u_sq,
      data_seed = data_seed
    )

    # extraction ---------------------------------
    model_result <- model_output$result
    model_messages <- model_output$messages
    model_warnings <- model_output$warnings
    model_error <- model_output$error

    # account for error -------------------------
    msg_prob_detected <- is_valid_str(paste0(model_messages, collapse = ""))
    warning_detected <- is_valid_str(paste0(model_warnings, collapse = ""))
    error_detected <- is_valid_str(paste0(model_error, collapse = ""))

    if (is_na_result(model_result) ||
      msg_prob_detected || warning_detected || error_detected) {
      # setting model_result to be vec of NA in case of
      # non-convergence
      model_result <- c(rep(NA, ncol(out_mat) - 1), 0)
      names(model_result) <- out_colnames
    }

    model_convergence <- as.logical(model_result["converged"])

    # checking total number of converged cases ----
    if (model_convergence) {
      conv_case_num <- conv_case_num + 1
    }

    # logging simulation run ----------------------
    log_sim_run(
      convergence = model_convergence, message = model_messages,
      warning = model_warnings, error = model_error,
      data_seed = data_seed, log_file = log_file, iter_no = i
    )

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
