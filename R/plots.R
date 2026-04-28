# plots

#' Plot annual catch
#'
#' @param year model year
#' @param output RTMButils model run output
#' @param folder the model folder
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_catch(year=2026, output = m24_26, folder = "m24_26")
#' }
plot_catch <- function(year, output, folder, save=TRUE){

  if (!dir.exists(here::here(year, folder, "figs"))){
    dir.create(here::here(year, folder, "figs"))
  }
  # set view
  ggplot2::theme_set(afscassess::theme_report())

  dat = output$dat
  rpt = output$rpt

  data.frame(year = rpt$years,
             obs = data$catch_obs,
             pred = rpt$catch_pred,
             years = "All years") -> df

  tidytable::filter(df, year %in% (max(df$year) - 20):max(df$year)) %>%
    tidytable::mutate(years = "Recent years") %>%
    tidytable::bind_rows(df) %>%
    tidytable::pivot_longer(c(-year, -years)) %>%
    ggplot2::ggplot(ggplot2::aes(year, value, color = name, lty = name)) +
    ggplot2::geom_line() +
    scico::scale_color_scico_d(name = "", palette = "roma") +
    ggplot2::scale_linetype_manual(name = "",
                                   values = c(1,1)) +
    ggplot2::facet_wrap(~years, scales = "free",
                        dir = "v") +
    ggplot2::ylab("Catch (kt)") +
    ggplot2::xlab("Year") +
    ggplot2::expand_limits(y = 0) +
    tickr::scale_x_tickr(data=df, var=year, by = 10, var_min = 1960) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::theme(legend.justification=c(1,0),
                   legend.position=c(0.98,0.8))

  if(isTRUE(save)) {
    ggplot2::ggsave(here::here(year, folder, "figs", "catch.png"),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)
  }

}

#' Plot biomass
#'
#' @param year model year
#' @param output RTMButils model run output
#' @param folder folder name model is in
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @export
#'
plot_biomass <- function(year, output, folder, save=TRUE) {

  if (!dir.exists(here::here(year, folder, "figs"))){
    create.dir(here::here(year, folder, "figs"))
  }
  # set view
  ggplot2::theme_set(afscassess::theme_report())

  summary(output$sd, "report") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("item")  %>%
    mutate(lci = Estimate - 1.96*`Std. Error`,
           uci = Estimate + 1.96*`Std. Error`) %>%
    mutate(item = gsub("\\..*", "", item),
           year = c( data$years, data$srv_yrs, rep(data$years, 3))) %>%
    dplyr::select(year, item, value = Estimate, se = `Std. Error`, lci, uci) -> df

  df %>%
    tidytable::filter(item  %in% c('spawn_bio', 'tot_bio')) %>%
    tidytable::mutate(item = tidytable::case_when(item == 'spawn_bio' ~ "Spawning biomass",
                                                  item == 'tot_bio' ~ "Total biomass"),
                      value = value / 1000,
                      lci = lci / 1000,
                      uci = uci / 1000) %>%
    ggplot2::ggplot(ggplot2::aes(year, value)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.1) +
    ggplot2::facet_wrap(~item, ncol = 1, scales = "free_y") +
    ggplot2::scale_y_continuous(name = "Biomass (kt)", labels = scales::comma) +
    ggplot2::expand_limits(y = 0) +
    tickr::scale_x_tickr(name = "Year", data=df, var=year, by=10, var_min = 1960)

  if(isTRUE(save)) {
    ggplot2::ggsave(here::here(year, folder, "figs", "est_biomass.png"),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)
  }
}




#' Plot age compositions
#'
#' @param year = assessment year
#' @param output RTMButils model run output
#' @param folder = folder the model lives in
#' @param type 'fishery' or 'survey'
#' @param save default is TRUE, saves fig to the folder the model is in
#' @export
#' @examples
#' plot_age_comps(year, output, folder, type = "fishery")
plot_age_comps <- function(year, output, folder, save = TRUE, type) {
  if (!dir.exists(here::here(year, folder, "figs"))){
    create.dir(here::here(year, folder, "figs"))
  }
  # set view
  ggplot2::theme_set(afscassess::theme_report())
  rpt = output$rpt
  dat = output$dat
  ages = dat$ages

  if(type == 'fishery') {
    obs =  as.data.frame(dat$fish_age_obs)
    pred = as.data.frame(rpt$fish_age_pred)
    yrs = dat$fish_age_yrs
  } else if(type == 'survey') {
    obs =  as.data.frame(dat$srv_age_obs)
    pred = as.data.frame(rpt$srv_age_pred)
    yrs = dat$srv_age_yrs
  } else {
    stop("type must be either 'fishery' or 'survey'")
  }

  cleanup <- function(var, ages, yrs) {
    var_name <- deparse(substitute(var))
    var %>%
      tidytable::bind_cols(age = ages) %>%
      tidytable::pivot_longer(-age) %>%
      tidytable::mutate(year = rep(yrs, each = length(ages)),
                        id = var_name)
  }

  obs = cleanup(obs, ages, yrs)
  pred = cleanup(pred, ages, yrs)

  p1 = obs %>%
    tidytable::mutate(Age = factor(age))  %>%
    dplyr::filter(id == "obs") %>%
    ggplot2::ggplot(ggplot2::aes(age, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Age), width = 1, color = "gray") +
    ggplot2::facet_wrap(~year, strip.position="right",
                        dir = "v",
                        ncol = 1) +
    ggplot2::geom_line(data = pred) +
    ggplot2::theme(panel.spacing.y = grid::unit(0, "mm")) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::xlab("Age") +
    ggplot2::ylab(paste(Hmisc::capitalize(type), "age composition")) +
    ggplot2::theme(legend.position = "none")

  if(isTRUE(save)) {
    ggplot2::ggsave(plot = p1, filename = here::here(year, folder, "figs", paste0(type, "_age_comp.png")),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)
  }
  p1
}


#' Plot size compositions
#'
#' @param year = assessment year
#' @param output RTMButils model run output
#' @param folder = folder the model lives in
#' @param type 'fishery' or 'survey'
#' @param save default is TRUE, saves fig to the folder the model is in
#' @export
#' @examples
#' plot_size_comps(year, output, folder, type = "fishery")
plot_size_comps <- function(year, output, folder, save = TRUE, type) {
  if (!dir.exists(here::here(year, folder, "figs"))){
    create.dir(here::here(year, folder, "figs"))
  }
  # set view
  ggplot2::theme_set(afscassess::theme_report())
  rpt = output$rpt
  dat = output$dat
  lengths = dat$length_bins

  if(type == 'fishery') {
    obs =  as.data.frame(dat$fish_size_obs)
    pred = as.data.frame(rpt$fish_size_pred)
    yrs = dat$fish_size_yrs
  } else if(type == 'survey') {
    obs =  as.data.frame(dat$srv_size_obs)
    pred = as.data.frame(rpt$srv_size_pred)
    yrs = dat$srv_age_yrs
  } else {
    stop("type must be either 'fishery' or 'survey'")
  }

  cleanup <- function(var, lengths, yrs) {
    var_name <- deparse(substitute(var))
    var %>%
      tidytable::bind_cols(length = lengths) %>%
      tidytable::pivot_longer(-length) %>%
      tidytable::mutate(year = rep(yrs, each = length(lengths)),
                        id = var_name)
  }

  obs = cleanup(obs, lengths, yrs)
  pred = cleanup(pred, lengths, yrs)

  p1 = obs %>%
    tidytable::mutate(Length = factor(length))  %>%
    dplyr::filter(id == "obs") %>%
    ggplot2::ggplot(ggplot2::aes(length, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Length), width = 1, color = "gray") +
    ggplot2::facet_wrap(~year, strip.position="right",
                        dir = "v",
                        ncol = 1) +
    ggplot2::geom_line(data = pred) +
    ggplot2::theme(panel.spacing.y = grid::unit(0, "mm")) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::xlab("Length (cm)") +
    ggplot2::ylab(paste(Hmisc::capitalize(type), "Length composition")) +
    ggplot2::theme(legend.position = "none")

  if(isTRUE(save)) {
    ggplot2::ggsave(plot = p1, filename = here::here(year, folder, "figs", paste0(type, "_size_comp.png")),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)
  }
  p1
}
