#' Run RTMB model
#'
#' @param model the model function
#' @param data named data list
#' @param pars named parameter list
#' @param map named mapping list, default: NULL
#' @param lower unnamed vector of lower parameter limits, default: NULL
#' @param upper unnamed vector of upper parameter limits, default: NULL
#' @param random vector of parameter(s) to be random effects, default: NULL
#'
#' @export
run_model <- function(model, data, pars, map=NULL, lower=NULL, upper=NULL, random = NULL) {
  cmb = function(f, d) function(p) f(p, d)
  obj =  RTMB::MakeADFun(cmb(model, data),
                         pars,
                         map = map,
                         random = random)

  if(!is.null(lower) & !is.null(upper)) {
    fit = nlminb(start = obj$par,
                 objective = obj$fn,
                 gradient = obj$gr,
                 control = list(iter.max=100000,
                                eval.max=20000),
                 lower = lower,
                 upper=upper)
  }
  fit = nlminb(start = obj$par,
               objective = obj$fn,
               gradient = obj$gr,
               control = list(iter.max=100000,
                              eval.max=20000))
  rpt = obj$report(obj$env$last.par.best)
  proj = proj_bio(rpt) # function to project the next 2 years
  sd = sdreport(obj)
  list(obj=obj, fit=fit, rpt=rpt, proj=proj, sd = sd)
}


# use Grant's test for model letter or number
#' @export
model_test <- function(m1, m2) {
  if(sqrt(sum(((m1$rpt$spawn_bio / m2$rpt$spawn_bio - 1)^2) / length(m2$rpt$years))) < 0.1) {
    return('letter')
  } else {
    return ("number")
  }
}

#' check Hessian is positive definite
#' @export
fit_check <- function(fit) {
  sd_fit = fit$sd
  cat("Is the Hessian positive definite:", sd_fit$pdHess,
      "\nThe maximum gradiant is:", max(abs(fit$obj$gr(fit$fit$par))),
      "\nThe gradiant is < 1e-5:", max(abs(fit$obj$gr(fit$fit$par))) < 1e-5)
}
