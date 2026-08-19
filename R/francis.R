#' Calculate Francis Reweighting Multiplier (Method TA1.8)
#'
#' @param rpt The report object from an RTMB model run.
#' @param data The data used in the RTMB model run
#' @param comp_type A string identifying the composition type, e.g., "fish_age", "srv_age", "fish_size".
#' @return The model report with recommended weight
calculate_francis_weight <- function(rpt, data, comp_type = "fish_size") {

  # Extract relevant data from the report object
  obs = data[[paste0(comp_type, "_obs")]]
  pred = rpt[[paste0(comp_type, "_pred")]]
  iss = data[[paste0(comp_type, "_iss")]]

  # Use ages for age comps, length bins for size comps
  if (grepl("age", comp_type)) {
    bins = data$ages
  } else {
    bins = data$length_bins
  }

  # standardization
  epsilon = 1e-6
  # mean age/length
  obs_mean = apply(obs, 2, function(x) sum(x * bins) / (sum(x) + epsilon))
  pred_mean = apply(pred, 2, function(x) sum(x * bins) / (sum(x) + epsilon))

  # variance of the predicted distribution for each year
  pred_var_dist = colSums(pred * bins^2) / (colSums(pred) + epsilon) - pred_mean^2

  # standardized residuals ---
  # expected standard error of the mean for each year
  sem = sqrt(pmax(pred_var_dist, 0) / iss)
  std_residuals = (obs_mean - pred_mean) / (sem + epsilon)

  # final multiplier (Francis TA1.8)---
  # the new weight is 1 / variance of the standardized residuals
  multiplier = 1 / var(std_residuals, na.rm = TRUE)

  return(multiplier)
}

#' Run model with iterative Francis reweighting for composition data
#'
#' @param output The output from RTMButils::run_model() that contains the model, named data list, and named parameter list.
#' @param iters number of reweighting iterations.
#' @param max_one cap the weight at 1, default: TRUE
#' @export
run_model_reweight <- function(output, iters = 10, max_one = TRUE) {
  data = output$dat
  map = output$obj$env$map
  fit = output$fit
  model = output$model
  pars = output$obj$env$parList()
  lower = fit$lower
  upper = fit$upper
  
  if (is.null(data$fish_age_wt))  data$fish_age_wt  <- 1
  if (is.null(data$srv_age_wt))   data$srv_age_wt   <- 1
  if (is.null(data$fish_size_wt)) data$fish_size_wt <- 1
  
  weight_history <- data.frame(
    iteration = 0, 
    fish_age_wt = data$fish_age_wt, 
    srv_age_wt = data$srv_age_wt, 
    fish_size_wt = data$fish_size_wt
  )
  
  cat("--- Starting Francis Reweighting ---\n")
  if (max_one) cat("--- Capping maximum weights at 1.0 ---\n")
  
  for (i in 1:iters) {
    cat(paste("\n>> Iteration:", i, "/", iters, "\n"))
    cat("Current Weights:\n")
    cat("  Fishery Age :", data$fish_age_wt, "\n")
    cat("  Survey Age  :", data$srv_age_wt, "\n")
    cat("  Fishery Size:", data$fish_size_wt, "\n")
    
    # run model with current weights
    new_run <- run_model(
      model = model, 
      data = data, 
      pars = pars, 
      map = map, 
      lower = lower, 
      upper = upper
    )
    
    # calculate new Francis multipliers
    mult_fa = RTMButils:::calculate_francis_weight(new_run$rpt, data, "fish_age")
    mult_sa = RTMButils:::calculate_francis_weight(new_run$rpt, data, "srv_age")
    mult_fs = RTMButils:::calculate_francis_weight(new_run$rpt, data, "fish_size")
    
    if (!is.finite(mult_fa)) {
      cat("  -- Warning: Fishery Age multiplier is Inf/NaN, resetting to 1.0\n")
      mult_fa <- 1
    }
    if (!is.finite(mult_sa)) {
      cat("  -- Warning: Survey Age multiplier is Inf/NaN, resetting to 1.0\n")
      mult_sa <- 1
    }
    if (!is.finite(mult_fs)) {
      cat("  -- Warning: Fishery Size multiplier is Inf/NaN, resetting to 1.0\n")
      mult_fs <- 1
    }
    
    # apply cap if is TRUE
	if(max_one) {
      mult_fa = min(mult_fa, 1.0)
      mult_sa = min(mult_sa, 1.0)
      mult_fs = min(mult_fs, 1.0)    
	}

    cat("Calculated Multipliers:\n")
    cat("  Fishery Age :", round(mult_fa, 3), "\n")
    cat("  Survey Age  :", round(mult_sa, 3), "\n")
    cat("  Fishery Size:", round(mult_fs, 3), "\n")
    
    # 3. Update data weights for next iteration
    data$fish_age_wt  <- mult_fa
    data$srv_age_wt   <- mult_sa
    data$fish_size_wt <- mult_fs
    

    weight_history[i + 1, ] <- c(i, data$fish_age_wt, data$srv_age_wt, data$fish_size_wt)
    
    # next iteration using fitted parameters from new_run
    pars <- new_run$obj$env$parList()
  }
  
  cat("\n--- Reweighting process finished. ---\n")
  
  return(list(
    model= new_run,
    weight_history = weight_history
  ))
}
