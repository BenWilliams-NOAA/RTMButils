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
#' @param control optimization settings list(iter.max = 1e5, eval.max = 2e4, rel.tol = 1e-12)
#'
#' @export
run_model <- function(model, data, pars, map=NULL, proj = TRUE, lower=NULL, upper=NULL, random = NULL, newton_loops = 3, control = list(iter.max = 1e5, eval.max = 2e4, rel.tol = 1e-12)) {

  # build AD function
  obj <- RTMB::MakeADFun(cmb(model, data), pars, map = map, random = random)
  
  # default bounds if not supplied
  if (is.null(lower)) lower <- rep(-Inf, length(obj$par))
  if (is.null(upper)) upper <- rep(Inf, length(obj$par))
  
  # optimization 
  fit <- nlminb(
    start = obj$par, 
    objective = obj$fn, 
    gradient = obj$gr, 
    control = control, 
    lower = lower, 
    upper = upper
  )
  
  #  bounded Newton-Raphson improvement Loops
  if (newton_loops > 0) {
    for (i in seq_len(newton_loops)) {
      g = as.numeric(obj$gr(fit$par))
      
      # AD Hessian from RTMB 
      h = tryCatch(obj$he(fit$par), error = function(e) NULL)
      if (is.null(h)) break
      
      # solve for newton step 
      delta = tryCatch(solve(h, g), error = function(e) NULL)
      if (is.null(delta)) break
      
      # new parameter vector & bounds
      prop_par = fit$par - delta
      prop_par = pmax(pmin(prop_par, upper), lower)
      
      # only accept step if NLL improves
      step_size = 1.0
      curr_obj = obj$fn(fit$par)
      improved = FALSE
      
      for (step_iter in 1:5) {
        test_par = fit$par - step_size * delta
        test_par = pmax(pmin(test_par, upper), lower) # bound enforcement
        test_obj = obj$fn(test_par)
        
        if (!is.na(test_obj) && test_obj < curr_obj) {
          fit$par = test_par
          fit$objective = test_obj
          improved= TRUE
          break
        }
        step_size = step_size * 0.5 # halve step size if NLL didn't improve
      }
      
      if (!improved) break # stop newton loops if no further improvement
    }
  }

  # evaluate obj$fn on final par to update obj$env internal state
  final_nll = obj$fn(fit$par)
  
  # reports and standard errors
  rpt = obj$report(fit$par)
  
  if (isTRUE(proj)) {
    proj_out = proj_bio(rpt)
  } else {
    proj_out = NULL
  }
  
  sd <- tryCatch(RTMB::sdreport(obj), error = function(e) {
    warning("sdreport failed: ", e$message)
    return(NULL)
  })
  
  return(list(
    obj   = obj, 
    fit   = fit, 
    rpt   = rpt, 
    proj  = proj_out, 
    sd    = sd, 
    dat   = data, 
    model = model, 
    lower = lower, 
    upper = upper
  ))
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
fit_check <- function(fit, tol = 1e-4, gr_tol = 1e-4) {
  # components
  sd_fit = fit$sd
  par = fit$fit$par
  lower = fit$lower
  upper = fit$upper
  
  # parameter names (handles single or vector parameters)
  par_names = names(par)
  if (is.null(par_names)) par_names = names(fit$obj$par)
  if (!is.null(par_names)) par_names = make.unique(par_names)
  
  # evaluate gradients
  gr = fit$obj$gr(par)
  max_gr = max(abs(gr), na.rm = TRUE)
  max_gr_param = par_names[which.max(abs(gr))]
  
  cat("----------------------------------------------------\n")
  cat("                  MODEL DIAGNOSTICS                 \n")
  cat("----------------------------------------------------\n")
  cat("Hessian Positive Definite :", ifelse(isTRUE(sd_fit$pdHess), "TRUE (Pass)", "FALSE (FAIL)"), "\n")
  cat("Maximum Gradient Value    :", sprintf("%.2e", max_gr), "\n")
  cat("Gradient <", sprintf("%.1e", gr_tol), "          :", ifelse(max_gr < gr_tol, "TRUE (Pass)", "FALSE (WARNING)"), "\n")
  cat("Highest Gradient Parameter:", max_gr_param, "\n")
  cat("----------------------------------------------------\n")
  
  cat("                  BOUNDS CHECK                      \n")
  cat("----------------------------------------------------\n")
  
  if (!is.null(lower) && !is.null(upper)) {
    # identify parameters close to lower or upper limits
    at_lower <- (par - lower) <= tol & !is.infinite(lower)
    at_upper <- (upper - par) <= tol & !is.infinite(upper)
    at_bound <- at_lower | at_upper
    
    if (any(at_bound, na.rm = TRUE)) {
      bound_df <- data.frame(
        Parameter = par_names[at_bound],
        Estimate  = round(par[at_bound], 5),
        Lower     = lower[at_bound],
        Upper     = upper[at_bound],
        Status    = ifelse(at_lower[at_bound], "at lower bound", "at upper bound")
      )
      cat("WARNING: The following parameter(s) are at/near bounds (tol =", tol, "):\n\n")
      print(bound_df, row.names = FALSE)
    } else {
      cat("Pass: All parameters are safely within their upper/lower bounds.\n")
    }
  } else {
    cat("Notice: 'lower' or 'upper' vectors not found in fit object.\n")
  }
  
  cat("----------------------------------------------------\n")
}