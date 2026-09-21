# plots

#' Phase plane plot
#'
#' @param year model year
#' @param output RTMB::run_model() output
#' @param folder folder model is in
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @export
#'
#' @examples plot_phase_plane(year=2026, output=m24_26, folder="m24-2026")
plot_phase_plane <- function(year, output, folder=NULL, save = TRUE) {
 if(is.null(folder) & isTRUE(save)) stop("need to name a folder to save the figure in.")
  if(!is.null(folder)) dir.create(here::here(year, folder, "figs"), showWarnings = FALSE)

    ggplot2::theme_set(afscassess::theme_report())
  
  Fabc_ratio <- output$proj$F40[1] / output$proj$F35[1]
  B_ratio <- output$rpt$B40 / output$rpt$B35

  segs <- data.frame(
  # BOTH lines now drop to zero at 0.05 * B_ratio and BOTH kink at B_ratio
  x1 = c(0.05 * B_ratio, B_ratio, 0.05 * B_ratio, B_ratio), 
  x2 = c(B_ratio, 2.8, B_ratio, 2.8),             
  y1 = c(0, 1, 0, Fabc_ratio),
  y2 = c(1, 1, Fabc_ratio, Fabc_ratio),
  group = factor(c("ofl", "ofl", "abc", "abc"),
                 levels = c("ofl", "abc"))
)
  
  p1 = data.frame(year = min(output$rpt$years):(max(output$rpt$years)+2),
           x = c(output$rpt$spawn_bio, output$proj$spawn_bio[1],  
                 output$proj$spawn_bio[2]) / output$rpt$B35,
           y = c(output$rpt$Ft, output$proj$F40 * yld$yld) / output$proj$F35[1]) %>% 
  tidytable::mutate(label = stringr::str_sub(year, 3),
                    decade = (floor(year / 10) * 10)) %>% 
  ggplot2::ggplot(ggplot2::aes(x, y)) +
  geom_path(aes(color = decade), show.legend = FALSE) +
  ggplot2::geom_label(ggplot2::aes(label=label, color = decade), linewidth = 0,
                      show.legend = FALSE, size = 3, family="Times", alpha = 0.5) +
  ggplot2::geom_segment(data = segs, ggplot2::aes(x=x1, y=y1, xend=x2, yend=y2, linetype=group)) +
  ggplot2::scale_linetype_manual(values = c(1, 3),
                                 labels = c(expression(italic(F[OFL])),
                                            expression(italic(F[ABC]))),
                                 name = "") +
  scico::scale_color_scico(palette = "roma") +
  ggplot2::ylab(expression(italic(F/F["35%"]))) +
  ggplot2::xlab(expression(italic(SSB/B["35%"]))) +
  ggplot2::theme(legend.justification=c(1,0),
                 legend.position=c(0.9,0.85)) 
  
  if(isTRUE(save)) {
    ggplot2::ggsave(plot = p1, filename = here::here(year, folder, "figs", "phase_plane.png"),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)
  }
  p1
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
    dir.create(here::here(year, folder, "figs"))
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
    dir.create(here::here(year, folder, "figs"))
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
0000
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


#' Recruitment/SSB plot
#'
#' @param year model year
#' @param output RTMB::run_model() output
#' @param folder folder model is in
#' @param save default is TRUE, saves fig to the folder the model is in
#'
#' @export
#'
#' @examples plot_rec_ssb(year=2026, output=m24_26, folder="m24-2026")
plot_rec_ssb <- function(year, output, folder = NULL, save=TRUE){

  if(is.null(folder) & isTRUE(save)) stop("need to name a folder to save the figure in.")
  if(!is.null(folder)) dir.create(here::here(year, folder, "figs"), showWarnings = FALSE)

    ggplot2::theme_set(afscassess::theme_report())
  rec_age = output$rpt$ages[1]

  p1 = data.frame(year = output$rpt$years,
           spawn_bio = output$rpt$spawn_bio,
           recruits = output$rpt$recruits) %>% 
  dplyr::mutate(spawn_bio = spawn_bio / 1000,
                recruits = dplyr::lead(recruits, n = rec_age),
                label = stringr::str_sub(year, 3),
                decade = (floor(year/10) * 10)) %>% 
  tidyr::drop_na() %>% 
   ggplot2::ggplot(ggplot2::aes(spawn_bio, recruits)) + 
  ggplot2::geom_label(ggplot2::aes(label = label, color = decade), 
                      label.size = 0, show.legend = FALSE, 
                      size = 4, family = "Times", alpha = 0.85) + 
  ggplot2::expand_limits(x = 0, y = 0) + 
  scico::scale_color_scico(palette = "roma") + 
  ggplot2::xlab("Female spawning biomass (kt)") + 
  ggplot2::ylab("Recruitment (millions)")
  
  if(isTRUE(save)) {
    ggplot2::ggsave(plot = p1, filename = here::here(year, folder, "figs", "rec_ssb.png"),
                    width = 6.5, height = 6.5, units = "in", dpi = 200)
  }
  p1

}