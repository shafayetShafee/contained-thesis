library(glmmTMB)

gen_three_level_int_data <- function(l, m, n, fixed_coeff, sigma_sq, data_seed) {
  # l => Number of EA
  # m => number of HH
  # n => size of each HH
  # sigma_sq => variance of (ujk and vk)
  N = l*m*n 
  
  set.seed(data_seed)
  
  x1c <- rnorm(N)
  x2 <- rnorm(N)
  x2b <- ifelse(x2 <= 0.5, 0, 1)
  
  if(length(fixed_coeff) != 3) {
    stop("`fixed_coeff` needs to be a vector of three element")
  }
  
  if(length(sigma_sq) != 2) {
    stop("`sigma_sq` needs to be a vector of two element")
  }
  
  beta0 <- fixed_coeff[1]
  beta1 <- fixed_coeff[2]
  beta2 <- fixed_coeff[3]
  
  sigma_u_sq <- sigma_sq[1]
  sigma_v_sq <- sigma_sq[2]
  
  ujk <- rep(rnorm(l*m, 0, sqrt(sigma_u_sq)), each = n)
  vk <- rep(rnorm(l, 0, sqrt(sigma_v_sq)), each = m*n)
  
  eta_ijk <- beta0 + beta1*x1c + beta2*x2b + ujk + vk
  pi_ijk <- exp(eta_ijk) / (1 + exp(eta_ijk))
  yijk <- rbinom(N, size = 1, prob = pi_ijk)
  
  multi_data <- data.frame(
    ea = rep(seq(l), each = m*n),
    hh = rep(rep(seq(m), each = n), times = l),
    id = rep(seq(n), times = l*m),
    ujk = ujk,
    vk = vk,
    X1c = x1c,
    X2b = x2b,
    Yijk = yijk
  )
  
  return(multi_data)
}


est_three_lvl_int_mor <- function(l, m, n, fixed_coeff, sigma_sq, data_seed) {
  
  # data generation ----------------
  multi_data_int <- gen_three_level_int_data(l, m, n, fixed_coeff, 
                                             sigma_sq, data_seed)
  
  # prevalence ---------------------
  prevalence <- mean(multi_data_int$Yijk, na.rm = TRUE)
  
  # true sigma sq ------------------
  sigma_u_sq <- sigma_sq[1]
  sigma_v_sq <- sigma_sq[2]
  
  # model fitting ------------------
  multi_model_int <- glmmTMB(Yijk ~ X1c + X2b + (1 | ea) + (1 | ea:hh),
                             data = multi_data_int, 
                             family = "binomial")
  
  model_sum <- summary(multi_model_int)
  
  # extracting value from model ----
  fixed_effect <- unname(model_sum$coefficients$cond[, "Estimate"])
  beta0_hat <- fixed_effect[1]
  beta1_hat <- fixed_effect[2]
  beta2_hat <- fixed_effect[3]
  
  # glmmTMB gives the var comp in log-sd scale
  
  # se(log(sd))
  se_logsd <- sqrt(diag(vcov(multi_model_int,full=TRUE)))[4:5] 
  # log(sd)
  logsd <- multi_model_int$sdr$par.fixed[4:5]
  
  # by delta method
  # se(sigma_sq_hat) = sqrt((2 * exp(2*logsd))^2) * se(logsd)
  se_sigma_sq_hat <- diag(unname(se_logsd*2*exp(2*logsd))) 
  # se of (sigma_vk, sigma_ujk)
  
  sigma_sq_hat <- unname(exp(2*logsd)) 
  # (sigma_vk, sigma_ujk)
  var_sigma_sq_hat <- se_sigma_sq_hat^2
  
  if(any(diag(var_sigma_sq_hat) <= 0)) {
    stop("Estimated var of ran effect is negative")
  }
  
  # mor1 calc. ----------------------
  true_mor1 <- exp(sqrt(2 * sigma_u_sq) * qnorm(0.75))
  mor1_hat <- exp(sqrt(2 * sigma_sq_hat[2]) * qnorm(0.75))
  log_mor1_hat <- log(mor1_hat)
  
  log_mor1_expr <- function(x) {
    # x is ran-effect parameterized as sd
    sqrt(2 * x) * qnorm(0.75)
  }
  
  J1 <- numDeriv::jacobian(log_mor1_expr, x = sigma_sq_hat[2])
  log_se_mor1_hat <- as.numeric(sqrt(t(J1) %*% var_sigma_sq_hat[2, 2] %*% J1))
  se_mor1_hat <- exp(log_se_mor1_hat)
  
  # mor2 calc. ----------------------
  true_mor2 <- exp(sqrt(2 * (sigma_u_sq + sigma_v_sq)) * qnorm(0.75))
  mor2_hat <- exp(sqrt(2 * sum(sigma_sq_hat)) * qnorm(0.75))
  log_mor2_hat <- log(mor2_hat)
  
  log_mor2_expr <- function(x) {
    # x is ran-effect parameterized as variance
    sqrt(2 * sum(x)) * qnorm(0.75)
  }
  
  J2 <- numDeriv::jacobian(log_mor2_expr, x = sigma_sq_hat)
  log_se_mor2_hat <- as.numeric(sqrt(J2 %*% var_sigma_sq_hat %*% t(J2)))
  se_mor2_hat <- exp(log_se_mor2_hat)
  
  # coverage calc ------------------
  ci_1 <- log_mor1_hat + c(-1, 1) * 1.96 * log_se_mor1_hat
  ci_1_exp <- exp(ci_1)
  coverage_1 <- as.numeric(ci_1_exp[1] <= true_mor1 && ci_1_exp[2] >= true_mor1)
  
  ci_2 <- log_mor2_hat + c(-1, 1) * 1.96 * log_se_mor2_hat
  ci_2_exp <- exp(ci_2)
  coverage_2 <- as.numeric(ci_2_exp[1] <= true_mor2 && ci_2_exp[2] >= true_mor2)
  
  # model conv check ---------------
  is_model_converged <- performance::check_convergence(multi_model_int)
  
  # creating output vector ---------
  out_vec <- c(true_mor1 = true_mor1, 
               mor1_hat = mor1_hat, 
               se_mor1_hat = se_mor1_hat, 
               coverage_1 = coverage_1,
               true_mor2 = true_mor2, 
               mor2_hat = mor2_hat, 
               se_mor2_hat = se_mor2_hat, 
               coverage_2 = coverage_2,
               sigma_u_sq_hat = sigma_sq_hat[2],
               sigma_v_sq_hat = sigma_sq_hat[1],
               beta0_hat = beta0_hat, 
               beta1_hat = beta1_hat, 
               beta2_hat = beta2_hat,
               prevalence = prevalence,
               converged = is_model_converged)
  
  if(is.na(is_model_converged) || !is_model_converged || mor1_hat > 20 ||
     mor2_hat > 20 || se_mor1_hat > 20 || se_mor2_hat > 20) {
    out_vec_names <- names(out_vec)
    out_vec <- c(rep(NA, length(out_vec) - 1), FALSE)
    names(out_vec) <- out_vec_names
  }
  
  return(out_vec)
}


simulate_three_lvl_int <- function(l, m, n, fixed_coeff, sigma_sq, nsims = 1000, 
                                 log_file, seed, more_iter, ...) {
  
  # creating extra sims to get nsims after accounting 
  # for non-converged cases ---------------------
  total_sims = nsims + more_iter
  conv_case_num = 0
  runs_required = 0
  
  # creating placeholder matrix for result --------
  out_mat <- matrix(NA, nrow = total_sims, ncol = 15)
  out_colnames <- c("true_mor1", "mor1_hat", "se_mor1_hat", "coverage_1", 
                    "true_mor2", "mor2_hat", "se_mor2_hat", "coverage_2", 
                    "sigma_u_sq_hat", "sigma_v_sq_hat", 
                    "beta0_hat", "beta1_hat", "beta2_hat", "prevalence", 
                    "converged")
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
      runs_required = runs_required + 1
    }
    
    data_seed <- iteration_seeds[i]
    model_output <- safe_and_quietly(fun = est_three_lvl_int_mor,
                                     l = l, m = m, n = n, 
                                     fixed_coeff = fixed_coeff,
                                     sigma_sq = sigma_sq, 
                                     data_seed = data_seed)
    
    # extraction ---------------------------------
    model_result <- model_output$result
    model_messages <- model_output$messages
    model_warnings <- model_output$warnings
    model_error <- model_output$error
    
    # account for error -------------------------
    msg_prob_detected <- is_valid_str(paste0(model_messages, collapse = ""))
    warning_detected <- is_valid_str(paste0(model_warnings, collapse = ""))
    error_detected <- is_valid_str(paste0(model_error, collapse = ""))
    
    if(is_na_result(model_result) || 
       msg_prob_detected || warning_detected || error_detected) {
      # setting model_result to be vec of NA in case of 
      # non-convergence
      model_result <- c(rep(NA, ncol(out_mat) - 1), 0)
      names(model_result) <- out_colnames
    }
    
    model_convergence <- as.logical(model_result["converged"])
    
    # checking total number of converged cases ----
    if(model_convergence) {
      conv_case_num = conv_case_num + 1
    }
    
    # logging simulation run ----------------------
    log_sim_run(convergence = model_convergence, message = model_messages, 
                warning = model_warnings, error = model_error,
                data_seed = data_seed, log_file = log_file, iter_no = i)
    
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
