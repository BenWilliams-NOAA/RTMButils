# note that the 2nd time period in the selectivity blocks is not scaled to one for POP
# so have this add on to account for that

#' @export
to_one <- function(x) {
  x / max(x)
}


#' standard MESA logistic selectivity
#'
#' @param age age to examine
#' @param log_a50 inflection point
#' @param delta steepness of slope
#' @param adj based upon index that starts at 1 - unless adjusted to start at recruitment age
#' doesn't really do anything other than scale
#' @export
sel_logistic <- function(age, log_a50, delta, adj) {
  a50 = exp(log_a50)
  sel = 1 / (1 + exp(-log(19) * ((age + adj) - a50) / delta))
  # sel = sel / max(sel)
  sel
}

#' MESA gamma selectivity
#'
#' @param age age to examine
#' @param log_b50 peak point
#' @param delta steepness of slope
#' @param adj based upon index that starts at 1 - unless adjusted to start at recruitment age
#' doesn't really do anything other than scale
#' @export
sel_gamma <- function(age, b50, delta, adj) {
  b50 = exp(log_b50)
  denom = 0.5 * (sqrt(b50^2 + 4 * delta^2) - b50)
  sel = (((age + adj) / b50)^(b50 / denom)) * exp((b50 - (age + adj)) / denom)
  # sel = sel / max(sel)
  sel
}

#' @export
sel_double_normal <- function(age, log_a50, delta, delta2, adj = 0) {
  a50 = exp(log_a50)
  x_val = age + adj - a50
  # delta for the ascending limb (left) and delta2 for the descending limb (right)
  sel <- ifelse(x_val <= 0,
                exp(-(x_val^2) / (2 * delta^2)),
                exp(-(x_val^2) / (2 * delta2^2)))
  # sel = sel / max(sel)
  sel
}

#' @export
sel_double_logistic <- function(age, log_a50_a, k_asc, log_a50_d, k_desc, adj = 0) {
  a50_a = exp(log_a50_a)
  a50_d = exp(log_a50_d)
  # ascending limb
  asc = 1 / (1 + exp(-k_asc * (age + adj - a50_a)))
  # escending limb
  desc = 1 / (1 + exp(-k_desc * (age + adj - a50_d)))
  sel = asc * desc
  # sel = sel / max(sel)
  sel
}

#' @export
get_slx <- function(ages, type, pars, adj) {
  switch(as.character(type),
    "1" = {
      # logistic: pars[1]=a50, pars[2] = delta
      sel_logistic(ages, pars[1], pars[2], adj)
    },
    "2" = {
      # Gamma: pars[1] = b50, pars[2] = delta
      sel_gamma(ages, pars[1], pars[2], adj)
    },
    "3" = {
      # double normal: pars[1] = a50, pars[2] = delta, pars[3] = delta2
      sel_double_normal(ages, pars[1], pars[2], pars[3], adj)
    },
    "4" = {
      # double logistic: pars[1] = a50_asc, pars[2] = k_asc, pars[3] = a50_desc, pars[4] = k_desc
      sel_double_logistic(ages, pars[1], pars[2], pars[3], pars[4], adj)
    }
  )
}