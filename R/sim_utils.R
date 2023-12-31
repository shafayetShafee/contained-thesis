# combined function of purrr::safely() and purrr::quietly() from
# https://github.com/tidyverse/purrr/issues/843 by @Maximilian-Stefan-Ernst

safe_and_quietly <- function(fun, ...) {
  safe_fun <- purrr::quietly(purrr::safely(fun))
  out_safe <- safe_fun(...)

  length_zero_to_na <- function(obj) {
    if (length(obj) == 0) {
      return(NA)
    } else {
      return(obj)
    }
  }

  out <-
    list(
      result = out_safe$result$result,
      error = out_safe$result$error,
      output = out_safe$output,
      warnings = out_safe$warnings,
      messages = out_safe$messages
    )
  if (!is.null(out$error)) {
    out$error <- conditionMessage(out$error)
  }
  out <- map(out, length_zero_to_na)
  return(out)
}


log_output <- function(x, type, file) {
  message(paste0(type, ": ", x, collapse = "\n"))
  cat(paste0(type, ": ", x, collapse = "\n"), "\n", file = file, append = TRUE)
}


log_sim_run <- function(convergence, message, warning, error,
                        data_seed, log_file, iter_no) {
  # started logging ----------------
  log_output(paste0(rep("-", 50), collapse = ""), type = "", file = log_file)
  log_output(data_seed, type = "Using data_seed", file = log_file)
  log_output(as.logical(convergence), type = "converged", file = log_file)
  log_output(message, type = "message", file = log_file)
  log_output(warning, type = "warning", file = log_file)
  log_output(error, type = "error", file = log_file)
  log_output(
    paste0("Stored output for iteration ", iter_no, "\n"),
    type = "Info", file = log_file
  )
}


is_valid_str <- function(x) {
  # `> 5` to prevent "NA"
  !is.null(x) && !is.na(x) && is.character(x) && nchar(trimws(x)) > 5
}

# not used any more
is_valid <- function(x) {
  !is.null(x) && !is.na(x) && x
}

is_na_result <- function(x) {
  length(x) == 1 && is.na(x)
}

sym_mat_to_vec <- function(x) x[upper.tri(x, diag = TRUE)]


D_chol_to_D <- function(x) {
  # transform the log-cholesky parameterized value back to original scale
  D <- GLMMadaptive:::chol_transf(x)
  D[upper.tri(D, diag = TRUE)]
}


vcov_orig_scale <- function(model) {
  D <- model$D

  # transform from covariance matrix to entries of cholesky factor with
  # log-transformed main diagonal
  D_chol_entries <- GLMMadaptive:::chol_transf(D)
  V_chol <- vcov(model, parm = "var-cov")

  J <- numDeriv::jacobian(D_chol_to_D, D_chol_entries)
  # delta method (estimated covariance matrix of entries of D)
  V <- J %*% V_chol %*% t(J)

  colnames(V) <- colnames(V_chol)
  rownames(V) <- rownames(V_chol)
  # se <- sqrt(diag(V))
  # names(se) <- colnames(V_chol)
  # return(se)

  return(V)
}


var_slope_rand_effect <- function(model) {
  rand_var <- vcov_orig_scale(model)
  v <- diag(rand_var)
  # v[1] => D_11
  # v[2] => D_12
  # v[3] => D_22
  v_mat <- diag(v)
  rownames(v_mat) <- c("D_11", "D_12", "D_22")
  colnames(v_mat) <- c("D_11", "D_12", "D_22")
  return(v_mat)
}

# not used so far
calc_slope_se_mor <- function(sigma_u_hat, var_ran_effect, x1_val) {
  log_mor_expr <- function(x, x1_val) {
    sqrt((2 * x[1]) + (2 * x1_val^2 * x[2]) + (4 * x1_val * x[3])) * qnorm(0.75)
  }

  J <- numDeriv::jacobian(log_mor_expr,
    c(sigma_u_hat[1], sigma_u_hat[3], sigma_u_hat[2]),
    x1_val = x1_val
  )

  var_log_mor <- as.numeric(J %*% var_ran_effect %*% t(J))
  return(sqrt(var_log_mor))
}


unnamed_quantile <- function(...) unname(quantile(...))


save_plot <- function(m, n, log_mor_hat, plot_name_suffix, plot_path, l = NULL) {
  hist_plot <- ggplot(drop_na(tibble(log_mor_hat)), aes(x = log_mor_hat)) +
    geom_histogram(bins = 30) +
    labs(
      x = "log(MOR)",
      y = "frequency"
    ) +
    theme_classic()

  # whether plotting is used in three lvl settingss
  sep <- "_"
  if (is.null(l)) sep <- ""

  ggsave(
    filename = paste0("hist_", l, sep, m, "_", n, "_", plot_name_suffix, ".png"),
    plot = hist_plot,
    path = plot_path
  )
}


# copied from
# http://rstudio-pubs-static.s3.amazonaws.com/28864_dd1f084207d54f5ea67c9d1a9c845d01.html
vcov.VarCorr.merMod <- function(object, fit, ...) {
  if (isREML(fit)) {
    warning("refitting model with ML")
    fit <- refitML(fit)
  }
  if (!require("numDeriv")) stop("numDeriv package required")
  useSc <- attr(object, "useSc")
  dd <- lme4:::devfun2(fit, useSc = useSc, signames = FALSE)
  vdd <- as.data.frame(object, order = "lower.tri")
  pars <- vdd[, "sdcor"]
  npar0 <- length(pars)
  if (isGLMM(fit)) {
    pars <- c(pars, fixef(fit))
  }
  hh1 <- optimHess(pars, dd)
  vv2 <- 2 * solve(hh1)
  if (isGLMM(fit)) {
    vv2 <- vv2[1:npar0, 1:npar0, drop = FALSE]
  }
  nms <- apply(
    vdd[, 1:3], 1,
    function(x) paste(na.omit(x), collapse = ".")
  )
  dimnames(vv2) <- list(nms, nms)
  return(vv2)
}


sd_mat_ran_eff <- function(model_obj) {
  sd_mat <- diag(diag(vcov(VarCorr(model_obj), model_obj)))
  return(sd_mat)
}
