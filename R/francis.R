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
  pred_var_dist <- apply(pred, 2, function(p) {
    mean_p = sum(p * bins) / (sum(p) + epsilon)
    var_p = sum(p * bins^2) / (sum(p) + epsilon) - mean_p^2
    return(var_p)
  })

  # standardized residuals ---
  # expected standard error of the mean for each year
  sem = sqrt(pred_var_dist / iss)
  std_residuals = (obs_mean - pred_mean) / (sem + epsilon)

  # final multiplier ---
  # the new weight is 1 / variance of the standardized residuals
  multiplier = 1 / var(std_residuals, na.rm = TRUE)

  return(multiplier)
}

#' Run model with iterative Francis reweighting for composition data
#'
#' @param model The model function.
#' @param data The initial named data list.
#' @param pars The initial named parameter list.
#' @param iters number of reweighting iterations.
#' @param ... Other arguments to pass to run_model (map, lower, upper, random).
#' @export
run_model_reweight <- function(model, data, pars, iters = 10, threshold = 0.05, ...) {

  # initial weights (start with 1.0 if not specified) ---
  if(is.null(data$fish_age_wt)) data$fish_age_wt <- 1.0
  if(is.null(data$srv_age_wt)) data$srv_age_wt <- 1.0
  if(is.null(data$fish_size_wt)) data$fish_size_wt <- 1.0


  # store the history of weights
  weight_history <- data.frame(
    iteration = 0,
    fish_age_wt = data$fish_age_wt,
    srv_age_wt = data$srv_age_wt,
    fish_size_wt = data$fish_size_wt
  )

  cat("--- Starting Francis Reweighting ---\n")

  for (i in 1:iters) {
    cat(paste("\n>> Iteration:", i, "/", iters, "\n"))
    cat("Current Weights:\n")
    cat("  Fishery Age:", data$fish_age_wt, "\n")
    cat("  Survey Age:", data$srv_age_wt, "\n")
    cat("  Fishery Size:", data$fish_size_wt, "\n")

    # run the model with the current weights
    current_run <- run_model(model = model, data = data, pars = pars, ...)

    mult_fa = calculate_francis_weight(current_run$rpt, data, "fish_age")
    mult_sa = calculate_francis_weight(current_run$rpt, data, "srv_age")
    mult_fs = calculate_francis_weight(current_run$rpt, data, "fish_size")

    # check for Inf/NaN and reset to 1.0 (no change) while issuing a warning
    if (!is.finite(mult_fa)) { cat("  -- Warning: Fishery Age multiplier is Inf/NaN, resetting to 1.0\n"); mult_fa <- 1.0 }
    if (!is.finite(mult_sa)) { cat("  -- Warning: Survey Age multiplier is Inf/NaN, resetting to 1.0\n"); mult_sa <- 1.0 }
    if (!is.finite(mult_fs)) { cat("  -- Warning: Fishery Size multiplier is Inf/NaN, resetting to 1.0\n"); mult_fs <- 1.0 }

    cat("Calculated Multipliers:\n")
    cat("  Fishery Age:", round(mult_fa, 3), "\n")
    cat("  Survey Age:", round(mult_sa, 3), "\n")
    cat("  Fishery Size:", round(mult_fs, 3), "\n")


    # update the weights in the data list for the next run
    data$fish_age_wt <- mult_fa
    data$srv_age_wt <- mult_sa
    data$fish_size_wt <- mult_fs

    # update history
    weight_history[i + 1, ] = c(i, data$fish_age_wt, data$srv_age_wt, data$fish_size_wt)

  }

  cat("\n--- Reweighting process finished. ---\n")
  print(weight_history)

  return(current_run)
}
