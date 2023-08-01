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
