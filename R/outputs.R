#' pull likelihoods from report
#'
#' @param rpt The output list from the `run_model()` function for the
#'   full dataset. Must contain `obj` (from RTMB::MakeADFun) and `rpt` (the report)
#' @param model name for column
#' @param addl any additional parameters to pull
#' @param exclude any parameters to exclude
#' @export
get_likes <- function(rpt, model = "Model", addl=NULL, exclude=NULL) {

  report = rpt$rpt
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

  names(df)[names(df) == "value"] <- model
  df
}

#' pull parameters from report and projection
#'
#' @param rpt model report
#' @param model name for column
#' @param addl any additional parameters to pull
#' @param exclude any parameters to exclude
#' @export
get_pars <- function(rpt, model = "Model", addl=NULL, exclude=NULL) {
  report = rpt$rpt
  prj = rpt$proj[1,]
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
