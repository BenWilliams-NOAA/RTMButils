#' Construct a Path Relative to the Project Root, Ignoring NULL Components
#'
#' @description
#' This function is a wrapper around [here::here()] that provides more flexible
#' behavior by automatically removing `NULL` elements from the path components
#' before constructing the final path. This is useful for programmatically
#' building paths where some subdirectories may be optional.
#'
#' @param ... Character vectors providing the path components. Any `NULL`
#'   values passed in `...` will be silently ignored.
#'
#' @return A character vector containing the full path.
#' @export
#'
#' @seealso [here::here()]
#'
#' @examples
#' \dontrun{
#' # The original here::here() returns an empty vector if NULL is included
#' # here::here("data", NULL, "file.csv")
#'
#' # herein() correctly ignores the NULL and builds the path
#' herein("data", NULL, "file.csv")
#'
#' # It's most useful when a path component is an optional variable
#' optional_folder <- NULL
#' herein("output", optional_folder, "plot.png")
#'
#' optional_folder <- "figures"
#' herein("output", optional_folder, "plot.png")
#' }


herein <- function(...) {
  # all arguments into a list
  path = list(...)

  # compact() to remove any NULL elements
  cp = purrr::compact(path)

  # call here::here() with the cleaned list of parts
  base::do.call(here::here, cp)
}

#' Jitter parameters
#'
#' @param pars list of named parameters
#' @param jitter_amt how much to jitter
#' @param map vector of parameters that will be mapped
#' @param seed for reproducibility, default: NULL
#'
#' @returns list of named, jittered parameters
#' @export
#'
jitter_pars <- function(pars, map = NULL, jitter_amt = 0.1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  jitter_value = function(x, jitter_amt) {
    if (length(x) > 1) {
      return(x + rnorm(length(x), mean = 0, sd = jitter_amt))
    } else {
      return(x + rnorm(1, mean = 0, sd = jitter_amt))
    }
  }
  # apply jitter to each parameter
  jits = lapply(names(pars), function(p_name) {
    # if the parameter is in the map, return it unchanged
    if (!is.null(map) && (p_name %in% map || p_name %in% names(map))) {
      return(pars[[p_name]])
    } else {
      return(jitter_value(pars[[p_name]], jitter_amt))
    }
  })

  # preserve names
  names(jits) <- names(pars)
  return(jits)

}


#' Extract and identify high parameter correlations
#'
#' This function converts the fixed-effect covariance matrix from an RTMB/TMB
#' model into a correlation matrix and identifies pairs of parameters that
#' exceed a specified correlation threshold.
#'
#' @param output Output from RTMButils::run_model. A list containing \code{sd}, which is an \code{sdreport} object.
#' @param threshold A numeric value between 0 and 1 indicating the absolute
#'   correlation cutoff for reporting. Defaults to 0.9.
#'
#' @return A square correlation matrix with parameter names as row and column names.
#' @export
#'
#' @examples
#' \dontrun{
#' get_correlations(model_output, threshold = 0.95)
#' }
get_correlations <- function(output, threshold = 0.9) {
  cov_fixed = output$sd$cov.fixed
  cor_matrix = stats::cov2cor(cov_fixed)
  high_cor = which(abs(cor_matrix) > threshold & row(cor_matrix) < col(cor_matrix), arr.ind = TRUE)

  if(nrow(high_cor) > 0) {
    cat("--- High Correlations Detected (>|", threshold, "|) ---\n", sep="")
    for(i in 1:nrow(high_cor)) {
      row_idx = high_cor[i, 1]
      col_idx = high_cor[i, 2]
      cat(sprintf("%s <--> %s: %.3f\n",
                  rownames(cor_matrix)[row_idx],
                  colnames(cor_matrix)[col_idx],
                  cor_matrix[row_idx, col_idx]))
    }
  } else {
    cat("No correlations found above threshold.\n")
  }
  return(cor_matrix)
}


#' Check if parameters are at estimation bounds
#'
#' Checks the estimated parameter values against the lower and upper bounds
#' specified during the optimization (e.g., in \code{nlminb}). Reports
#' parameters that are within a small tolerance of their limits.
#'
#' @param output Output from RTMButils::run_model. A list containing \code{fit}, which is the optimization
#'   result (e.g., from \code{nlminb}) containing \code{par}, \code{lower},
#'   and \code{upper}.
#' @param tol A numeric tolerance. Parameters within this distance of a
#'   bound will be reported. Defaults to 0.001.
#'
#' @return No return value. Prints a summary of parameters at their bounds to the console.
#' @export
#'
#' @examples
#' \dontrun{
#' check_bounds(mod_output)
#' }
check_bounds <- function(output, tol = 0.001) {

  if(is.null(output$fit$lower) & is.null(output$fit$upper)) {
    stop("No bounds found in the model fit object.")
  }

  est = output$fit$par
  low = output$fit$lower
  upp = output$fit$upper

  # identify which are at or near bounds
  at_low = which(est <= (low + tol))
  at_upp = which(est >= (upp - tol))

  cat("--- Boundary Check ---\n")
  if(length(at_low) > 0) {
    cat("Parameters at LOWER bound:\n")
    print(data.frame(Param = names(est)[at_low], Value = est[at_low], Bound = low[at_low]))
  } else { cat("None at lower bound.\n") }

  cat("\n")

  if(length(at_upp) > 0) {
    cat("Parameters at UPPER bound:\n")
    print(data.frame(Param = names(est)[at_upp], Value = est[at_upp], Bound = upp[at_upp]))
  } else { cat("None at upper bound.\n") }
}
