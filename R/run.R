#' Solve for target F
#' differentiable bisection solver
#' @param target the target F ratio
#' @param M natural mortality
#' @param slx fishery selectivity (dims: A,T)
#' @param wt_mat mature weight at age
#' @param sp_fract when spawning occurs
#' @param A number of age classes
#' @param SB0_val spawning biomass per recruit
#'
#' @export
solve_Fx <- function(target, M, slx, wt_mat, sp_fract, A, SB0_val) {
  low = 0.0
  high = 5.0
  for(i in 1:20) {
    F_test = (low + high) / 2.0
    N_test = rep(1,A)
    for(a in 2:A) {
      N_test[a] = N_test[a-1] * exp(-(M + F * slx[a-1,T]))
    }
    N_test[A] = N_test[A] * exp(-(M + F * slx[A,T]))
    sbpr_test = sum(N_test * wt_mat * exp(-sp_fract * (M + F_test * slx[,T])))
    ratio = sbpr_test / SB0_val
    low = ifelse(ratio > target, F_test, low)
    high = ifelse(ratio < target, F_test, high)
  }
  return((low + high) / 2.0)
}

#' Run RTMB model
#'
#' @param model the model function
#' @param data named data list
#' @param pars named parameter list
#' @param map named mapping list, default: NULL
#' @param proj whether to run the projection module: TRUE
#' @param lower unnamed vector of lower parameter limits, default: NULL
#' @param upper unnamed vector of upper parameter limits, default: NULL
#' @param random vector of parameter(s) to be random effects, default: NULL
#' @param newton_loops number of newton loops to run to reduce gradient: 3 - note only works for unconstrained models
#'
#' @export
run_model <- function(model, data, pars, map=NULL, proj = TRUE, lower=NULL, upper=NULL, random = NULL, newton_loops = 3) {
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

      try_improve <- tryCatch({
        for (i in 1:newton_loops) {
          g <- as.numeric(obj$gr(fit$par))
          h <- optimHess(fit$par, fn = obj$fn, gr = obj$gr)
          fit$par <- fit$par - solve(h, g)
          fit$objective <- obj$fn(fit$par)
        }
      }, error = function(e) {
        # If it fails, print a warning and continue with the original fit
        warning("Newton improvement step failed: ", e$message)
      })

  } else {
  fit = nlminb(start = obj$par,
               objective = obj$fn,
               gradient = obj$gr,
               control = list(iter.max=100000,
                              eval.max=20000))
  }

  rpt = obj$report(obj$env$last.par.best)
  if(isTRUE(proj)) {
    proj = proj_bio(rpt) # function to project the next 2 years
  } else {
    proj = NULL
  }
  sd = RTMB::sdreport(obj)
  list(obj=obj, fit=fit, rpt=rpt, proj=proj, sd=sd, dat=data, model=model)
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
