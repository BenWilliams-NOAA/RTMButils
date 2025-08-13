#' Profile a parameter using RTMB with re-optimization
#'
#' @param par_name parameter to profile (e.g., "log_M")
#' @param par_values values to profile the parameter over
#' @param derived character vector of derived quantity likelihoods to extract from report
#' @param data input data to RTMB model
#' @param map any other parameters that are fixed
#' @param obj RTMB model object function (from your original model_fn)
#' @param fit Original fit from rtmb::fit()
#' @param data The data input used in the model
#' @param model The model to fit the data to
#'
#' @return A data frame with parameter value, log-likelihood, and derived values
#' @export
profiles <- function(par_name = "log_M", par_values = log(seq(0.03, 0.1, 0.005)),
                     derived = c("like_fish_age", "like_srv_age"),
                     data, map = NULL, obj, fit, model) {

  results <- data.frame(value = par_values,
                        log_like = NA_real_,
                        matrix(NA_real_, nrow = length(par_values), ncol = length(derived)))
  colnames(results)[3:ncol(results)] <- derived

  nms = unique(names(fit$par))
  split_list <- split(fit$par, names(fit$par))
  pars = lapply(split_list, unname)
  pars = pars[nms]

  # put any mapped items back into the pars
  if(!is.null(map)) {
    pars[[names(map)]] = obj$report()[[names(map)]]
  }
  # create map if none provided
  if(is.null(map)) {
    map = list()
  }
  map[[par_name]] = factor(NA)

  for (i in seq_along(par_values)) {
    cat("Profiling", par_name, "at", par_values[i], "\n")

    # copy parameter list
    pars_i <- pars
    pars_i[[par_name]] <- par_values[i]

    # rebuild obj with fixed parameter
    obj_i  = RTMB::MakeADFun(cmb(model, data),,
                             parameters = pars_i,
                             map = map)

    # optimize the rest
    fit_i <- nlminb(start = obj_i$par,
                    objective = obj_i$fn,
                    gradient = obj_i$gr,
                    control = list(iter.max = 100000, eval.max = 20000))

    results$log_like[i] <- fit_i$objective

    # store raw derived values
    rep <- obj_i$report(fit_i$par)
    for (j in seq_along(derived)) {
      key <- derived[j]
      val <- rep[[key]]
      results[i, key] <- if (is.atomic(val) && length(val) > 1) tail(val, 1) else val
    }
  }

  results

}



#' plot parameter profile
#'
#' @param data profile data.frame with columns value, log_like, ...
#' @param exp  exponentiate the parameter value?
#' @export
plot_profile <- function(data, exp = FALSE) {
  if(isTRUE(exp)) {
    data %>%
      dplyr::mutate(value = exp(value)) -> data
  }
  data %>%
    dplyr::mutate(dplyr::across(-value, ~ .x - min(.x))) %>%
    tidyr::pivot_longer(-value, values_to = "nll") %>%
    ggplot2::ggplot(ggplot2::aes(value, nll, color = name)) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 0, lty=3) +
    # ggplot2::coord_cartesian(y = c(0,3)) +
    scico::scale_color_scico_d(palette = 'roma') +
    ylab("Change in -log likelihood")
}

check_model_fits <- function(fit, pars) {
  # Basic convergence info
  cat("=== Convergence Information ===\n")
  cat("Convergence code:", fit$convergence, "\n")
  cat("Message:", fit$message, "\n")
  cat("Iterations:", fit$iterations, "\n")
  cat("Objective value:", fit$objective, "\n\n")

  # Parameter bounds check
  cat("=== Parameter Bounds Check ===\n")
  at_bounds = which(abs(fit$par - unlist(lower)) < 1e-5 |
                       abs(fit$par - unlist(upper)) < 1e-5)
  if(length(at_bounds) > 0) {
    cat("Parameters at bounds:", names(fit$par)[at_bounds], "\n\n")
  } else {
    cat("No parameters at bounds\n\n")
  }

  # Gradient check
  cat("=== Gradient Information ===\n")
  final_grad = obj$gr(fit$par)
  max_grad = max(abs(final_grad))
  cat("Maximum absolute gradient:", max_grad, "\n")
  if(max_grad > 0.001) {
    cat("WARNING: Large gradients present\n\n")
  }

  # Compare initial vs final values
  cat("=== Parameter Changes ===\n")
  initial = unlist(pars)
  final = fit$par
  rel_change = abs((final - initial)/initial)
  changed = which(rel_change > 0.5)
  if(length(changed) > 0) {
    changes_df = data.frame(Parameter = names(initial)[changed],
                            Initial = initial[changed],
                            Final = final[changed],
                            Rel_Change = rel_change[changed]
                          )
    print(changes_df[order(-changes_df$Rel_Change), ])
  } else {
    cat("No parameters changed by more than 50%\n")
  }
}
