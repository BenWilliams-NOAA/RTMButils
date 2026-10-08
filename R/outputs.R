#' pull likelihoods from report
#'
#' @param output The output from the `RTMButils::run_model()` function for the
#'   full dataset. Must contain `obj` (from RTMB::MakeADFun) and `rpt` (the report)
#' @param model name for column
#' @param addl any additional parameters to pull
#' @param exclude any parameters to exclude
#' @export
get_likes <- function(output, model = "Model", addl=NULL, exclude=NULL) {

  report = output$rpt
  items = paste(c("like", "nll", "spr", "regularity", "ssqcatch", addl), collapse = "|")
  selected = report[grep(items, names(report))]
  if (!is.null(exclude)) {
    exclude = paste(exclude, collapse = "|")
    selected = selected[!grepl(exclude, names(selected))]
  }

  df <- data.frame(
    item = names(selected),
    value = round(unlist(selected),4),
    row.names = NULL
  )

  # create a new data frame row for the parameter count
  pars_df <- data.frame(
    item = "n_pars",
    value = length(output$fit$par)
  )

  # Combine the two data frames
  df <- rbind(df, pars_df)

  names(df)[names(df) == "value"] <- model
  df
}

#' pull parameters from report and projection
#'
#' @param output The output from the `RTMButils::run_model()` function for the
#'   full dataset. Must contain `obj` (from RTMB::MakeADFun) and `rpt` (the report)
#' @param model model name for column
#' @param addl any additional parameters to pull
#' @param exclude any parameters to exclude
#' @export
get_pars <- function(output, model = "Model", addl=NULL, exclude=NULL) {
  report = output$rpt
  prj = output$proj[1,]
  items = paste0("^", c("M", "q", "log_mean_R", "log_mean_F", "a50C", "deltaC", "a50S", "deltaS", "sigma", addl, "$"), collapse = "|")
  selected = report[grep(items, names(report))]

  items2 = data.frame(item = c("tot_bio", "spawn_bio", "catch_ofl", "F35", "catch_abc", "F40"),
                      value = round(c(prj$tot_bio, prj$spawn_bio, prj$catch_ofl, prj$F35, prj$catch_abc, prj$F40), 4))

  # Flatten values and preserve indices
  flat = lapply(seq_along(selected), function(i) {
    value = selected[[i]]
    base_name = names(selected)[i]

    # If value is a vector, give indexed names
    if (length(value) > 1) {
      data.frame(
        item = paste0(base_name, seq_along(value)),
        value = value
      )
    } else {
      data.frame(
        item = base_name,
        value = round(value, 4)
      )
    }
  })

  df <- do.call(rbind, flat)

  df = dplyr::bind_rows(df, items2)

  if (!is.null(exclude)) {
    exclude = paste(exclude, collapse = "|")
    df = df[!grepl(exclude, df$item),]
  }


  names(df)[names(df) == "value"] <- model
  df

}

#' Zero Out Recent Data Observations and Refit Model
#'
#' Zeroes out the last \code{yrs} observations of a specified data indicator vector
#' in the model data object and refits the assessment model.
#'
#' @param output The output from the `RTMButils::run_model()`
#' @param item Character string specifying the index to to zero out. Default is \code{"srv_ind"}.
#' @param yrs Integer indicating the number of recent observations/years 
#'   to zero out at the end of the time series. Default is \code{2}.
#'
#' @return An updated model run object returned by \code{run_model()}.
#' @export
rmv <- function(output, item = "srv_ind", yrs = 2) {
  dat =output$dat
  f = output$model
  l = output$lower
  u = output$upper
  map = output$obj$env$map
  pars = output$obj$env$parList(output$fit$par)

  idx <- tail(seq_along(dat[[item]]), yrs)
  dat[[item]][idx] <- 0

  new_run <- run_model(
      model = f, 
      data = dat, 
      pars = pars, 
      map = map, 
      lower = l, 
      upper = u
    )
  new_run
}

#' Drop Data Component Weighting and Refit Model
#'
#' Zeroes out likelihood weights for a specified data component in the model dataset
#' and refits the model, effectively removing its influence from the fit.
#'
#' @param output The output from the `RTMButils::run_model()`
#' @param item Character string specifying the wt to to zero out. Default is \code{"srv_wt"}.
#'
#' @return An updated model run object returned by \code{run_model()}.
#' @export
drop_wt <- function(output, item = "srv_wt") {
  dat = output$dat
  f = output$model
  l = output$lower
  u = output$upper
  map = output$obj$env$map
  pars = output$obj$env$parList()

  if(item == "catch_wt" & length(dat$catch_wt>1)) {
    dat$catch_wt = rep(0, length(dat$catch_wt))
  } else {
    dat[[item]] <- 0
  }

  new_run <- run_model(
      model = f, 
      data = dat, 
      pars = pars, 
      map = map, 
      lower = l, 
      upper = u
    )
  new_run
}

