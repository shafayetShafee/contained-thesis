# combined function of purrr::safely() and purrr::quietly() from
# https://github.com/tidyverse/purrr/issues/843 by @Maximilian-Stefan-Ernst

safe_and_quietly <- function(fun, ...){
  safe_fun <- purrr::quietly(purrr::safely(fun))
  out_safe <- safe_fun(...)
  
  length_zero_to_na <- function(obj){
    if(length(obj)==0){
      return(NA)
    }else{
      return(obj)
    }
  }
  
  out <- 
    list(
      result = out_safe$result$result,
      error = out_safe$result$error,
      output = out_safe$output,
      warnings = out_safe$warnings,
      messages = out_safe$messages)
  if(!is.null(out$error)){
    out$error <- conditionMessage(out$error)
  }
  out <- map(out, length_zero_to_na)
  return(out)
}


log_output <- function(x, type, file) {
  message(paste0(type, ": ", x, collapse = "\n"))
  cat(paste0(type, ": ", x, collapse = "\n"), "\n", file = file, append = TRUE)
}


is_valid_str <- function(x) {
  # `> 5` to prevent "NA"
  !is.null(x) && !is.na(x) && is.character(x) && nchar(trimws(x)) > 5
}


is_valid <- function(x) {
  !is.null(x) && !is.na(x) && x
}


sym_mat_to_vec <- function(x) x[upper.tri(x, diag = TRUE)]


D_chol_to_D <- function(x) {
  # transform the log-cholesky parameterized value back to original scale
  D <- GLMMadaptive:::chol_transf(x)
  D[upper.tri(D, diag = TRUE)]
}


vcov_orig_scale <- function(model, ...) {
  D <- model$D
  
  # transform from covariance matrix to entries of cholesky factor with 
  # log-transformed main diagonal
  D_chol_entries <- GLMMadaptive:::chol_transf(D)
  V_chol <- vcov(model, parm = "var-cov")
  
  J <- numDeriv::jacobian(D_chol_to_D, D_chol_entries)
  # delta method (estimated covariance matrix of entries of D)
  V <- J %*% V_chol %*% t(J)
  
  # se <- sqrt(diag(V))
  # names(se) <- colnames(V_chol)
  # return(se)
  
  return(V)
}
