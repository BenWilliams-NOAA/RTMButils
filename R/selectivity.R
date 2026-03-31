# note that the 2nd time period in the selectivity blocks is not scaled to one for POP
# so have this add on to account for that

#' @export
to_one <- function(x) {
  x / max(x)
}


#' standard MESA logistic selectivity
#'
#' @param age age to examine
#' @param a50 inflection point
#' @param delta steepness of slope
#' @param adj based upon index that starts at 1 - unless adjusted to start at recruitment age
#' doesn't really do anything other than scale
#' @export
sel_logistic <- function(age, a50, delta, adj=0) {
  x = age + adj
  sel = 1 / (1 + exp(-log(19) * (x - a50) / delta))
  # sel = sel / max(sel)
  sel
}

#' MESA gamma selectivity
#'
#' @param age age to examine
#' @param b50 peak point
#' @param delta steepness of slope
#' @param adj based upon index that starts at 1 - unless adjusted to start at recruitment age
#' doesn't really do anything other than scale
#' @export
sel_gamma <- function(age, b50, delta, adj=0) {
  x = age + adj
  denom = 0.5 * (sqrt(b50^2 + 4 * delta^2) - b50)
  sel = ((x / b50)^(b50 / denom)) * exp((b50 - x) / denom)
  # sel = sel / max(sel)
  sel
}

#' @export
sel_double_normal <- function(age, a50_a, a50_d, delta, delta2, adj = 0) {

  a50d = a50a + a50_d # plateu width
  x = age + adj
  # delta for the ascending limb (left) and delta2 for the descending limb (right)
  sel = ifelse(x <= a50a,
                exp(-((x - a50a)^2 / (2 * delta^2))),
                ifelse(x <= a50d,
                        1.0,
                        exp(-((x - a50d)^2 / (2 * delta2^2)))))                  
  # sel = sel / max(sel)
  sel
}

#' @export
sel_double_logistic <- function(age, a50_a, delta, a50_d, delta2, adj = 0) {
  a50_d = a50_a + a50_d
  x = age + adj
  # ascending limb
  asc = 1 / (1 + exp(-delta * (x - a50_a)))
  # escending limb
  desc = 1 / (1 + exp(delta2 * (x - a50_d)))
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
      # double logistic: pars[1] = a50_a, pars[2] = delta, pars[3] = a50_d, pars[4] = delta2
      sel_double_logistic(ages, pars[1], pars[2], pars[3], pars[4], adj)
    }
  )
}