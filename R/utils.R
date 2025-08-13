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
#' @param seed for reproducibility, default: NULL
#'
#' @returns list of named, jittered parameters
#' @export
#'
jitter_pars <- function(pars, jitter_amt = 0.1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  jitter_value = function(x, jitter_amt) {
    if (length(x) > 1) {
      return(x + rnorm(length(x), mean = 0, sd = jitter_amt))
    } else {
      return(x + rnorm(1, mean = 0, sd = jitter_amt))
    }
  }
  # Apply jitter to each parameter
  jits = lapply(pars, function(p) jitter_value(p, jitter_amt))

  # Preserve names
  names(jits) <- names(pars)
  return(jits)

}
