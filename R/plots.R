# plots

#' Phase plane plot
#'
#' @param year model year
#' @param output RTMB::run_model() output
#' @param folder folder model is in
#' @param save default is TRUE, saves fig to the folder the model is in
#' @param ... add'l ggplot inputs
#'
#' @export
#'
#' @examples plot_phase_plane(year=2026, output=m24_26, folder="m24-2026")
plot_phase_plane <- function(..., year, output, folder=NULL, save = TRUE, shade = FALSE) {
 if(is.null(folder) & isTRUE(save)) stop("need to name a folder to save the figure in.")
  
  args = list(...)
  is_gg_layer = purrr::map_lgl(args, ~ inherits(.x, "gg") || inherits(.x, "theme") || inherits(.x, "labels"))
  extra_layers = args[is_gg_layer]
  if (!is.null(folder) && !is.character(folder)) {
    folder = NULL
  }

   
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
  
  shade_abc <- data.frame(
    x = c(0.05 * B_ratio, B_ratio, B_ratio, 0.05 * B_ratio),
    y = c(0, Fabc_ratio, 0, 0)
  )

  p1 <- data.frame(
    year = min(output$rpt$years):(max(output$rpt$years) + 2),
    x = c(output$rpt$spawn_bio, output$proj$spawn_bio[1], output$proj$spawn_bio[2]) / output$rpt$B35,
    y = c(output$rpt$Ft, output$proj$F40 * yld$yld) / output$proj$F35[1]
  ) %>% 
    tidytable::mutate(
      label = stringr::str_sub(year, 3),
      decade = (floor(year / 10) * 10)
    ) %>% 
    ggplot2::ggplot(ggplot2::aes(x, y)) +
    afscassess::theme_report() 

  if(isTRUE(shade)) {
  p1 <- p1 +
    # Add light gray shaded ribbon under the Fabc slope
    ggplot2::geom_polygon(
      data = shade_abc, 
      ggplot2::aes(x = x, y = y), 
      fill = "grey85", 
      alpha = 0.3, 
      inherit.aes = FALSE
    ) 
     
  }
  p1 = p1 +
    ggplot2::geom_path(ggplot2::aes(color = decade), show.legend = FALSE) +
    ggplot2::geom_label(
      ggplot2::aes(label = label, color = decade), 
      linewidth = 0, show.legend = FALSE, size = 3, family = "Times", alpha = 0.5
    ) +
    ggplot2::geom_segment(data = segs, ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2, linetype = group)) +
    ggplot2::scale_linetype_manual(
      values = c(1, 3),
      labels = c(expression(italic(F[OFL])), expression(italic(F[ABC]))),
      name = ""
    ) +
    scico::scale_color_scico(palette = "roma") +
    ggplot2::ylab(expression(italic(F/F["35%"]))) +
    ggplot2::xlab(expression(italic(SSB/B["35%"]))) +
    ggplot2::theme(
      legend.justification = c(1, 0),
      legend.position = c(0.9, 0.85)
    )
  
  if (length(extra_layers) > 0) {
    for (layer in extra_layers) {
      p1 <- p1 + layer
    }
  }
  if(isTRUE(save)&& !is.null(folder)) {
    dir.create(herein(year, folder, "figs"), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(plot = p1, filename = herein(year, folder, "figs", "phase_plane.png"),
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


#' Plot Standard Deviation Ribbons for Stock Assessment Model Quantities
#'
#' Extracts specified time series variables and their standard errors from one 
#' or more model objects, calculates 95% confidence intervals, and generates 
#' a comparative ribbon plot. Optionally saves the plot to a specified directory.
#'
#' @param ... One or more model objects (or a single named list of model objects). 
#'   Each model object must contain an \code{sd} summary component (e.g., \code{sdreport}) 
#'   and a \code{rpt$years} vector.
#' @param var_pattern A character string containing a regular expression pattern 
#'   matching the target variable name in the \code{sdreport} summary 
#'   (e.g., \code{"^spawn_bio"}). Default is \code{"^spawn_bio"}.
#' @param palette A character string specifying the \code{scico} palette name to use 
#'   for line and ribbon coloring. Default is \code{"roma"}.
#' @param year Integer or character. The assessment year used for building the output file path.
#' @param folder Optional character string specifying a subdirectory path passed to 
#'   \code{RTMButils::herein}. Default is \code{NULL}.
#' @param base_size Numeric. Base font size for \code{afscassess::theme_report}. Default is \code{11}.
#' @param save Logical. If \code{TRUE}, saves the plot as a PNG to the specified path. 
#'   Default is \code{TRUE}.
#' @param log_space Logical. If \code{TRUE}, assumes estimates and standard errors 
#'   are on the log-scale and exponentiates them (\code{exp(value)}) along with 
#'   log-normal confidence intervals. Default is \code{FALSE}.
#' @param is_deviate Logical. If \code{TRUE}, plots raw log-space deviations 
#'   (e.g., \code{log_Ft}) as linear values without adding mean parameters or 
#'   exponentiating. Default is \code{FALSE}.
#' @param ylab default is \code{"Spawning biomass (t)"}
#' 
#' @return A \code{ggplot2} object showing the trajectory of the target variable 
#'   with 95\% confidence interval ribbons across models.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Plotting individual model objects passed as arguments
#' plot_sd_ribbon(mod1, mod2, var_pattern = "^spawn_bio", year = 2024)
#'
#' # Plotting a named list of models with custom palette
#' model_list <- list("Base" = mod1, "Alternative" = mod2)
#' plot_sd_ribbon(model_list, var_pattern = "^rec", palette = "batlow", year = 2024, save = FALSE)
#' }
#' 
plot_sd_ribbon <- function(..., var_pattern = "^spawn_bio", 
                            palette = "roma", year, folder = NULL, 
                            base_size = 11, save = TRUE, log_space = FALSE,
                            is_deviate = FALSE,
                            ylab = "Spawning biomass (t)") {
  
  args <- list(...)
  
  # separate model objects from additional ggplot layers passed in ...
  is_layer = purrr::map_lgl(args, ~ inherits(.x, "gg") || inherits(.x, "theme") || inherits(.x, "labels"))
  extra_layers = args[is_layer]
  mods = args[!is_layer]
  
  # Passing a single list of models (e.g., plot_sd_ribbon(model_list))
  if (length(mods) == 1 && is.list(mods[[1]]) && !("sd" %in% names(mods[[1]]))) {
    mods <- mods[[1]]
  }
  
  # Names for plot legends
  mod_names = names(mods)
  expr_names = as.character(substitute(list(...)))[-1]
  
  if (is.null(mod_names)) {
    mod_names = expr_names
  } else {
    unnamed = mod_names == ""
    if (any(unnamed)) {
      mod_names[unnamed] = expr_names[unnamed]
    }
  }
  
  # Extract and calculate CIs across all models
  df <- purrr::map2_dfr(mods, mod_names, function(mod, id_name) {
    sd_sum <- base::summary(mod$sd, select = "all") %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column(var = "item") %>%
      tibble::as_tibble() %>%
      dplyr::rename(value = Estimate, se = `Std. Error`)
    
    # Case 1: Special handling for log_Ft deviations
    # if (stringr::str_detect(var_pattern, "log_Ft") && !is_deviate) {
    #   log_mean_F <- sd_sum %>% 
    #     dplyr::filter(item == "log_mean_F") %>% 
    #     dplyr::pull(value)
      
    #   res <- sd_sum %>%
    #     dplyr::filter(stringr::str_detect(item, var_pattern)) %>%
    #     dplyr::mutate(
    #       log_F_total = value + log_mean_F,
    #       se_F = exp(log_F_total) * se,
    #       value = exp(log_F_total),
    #       lci = pmax(0, value - 1.96 * se_F),
    #       uci = value + 1.96 * se_F
    #     )
      
    # # Case 2: Standard log-space variables (e.g., log_recruits)
    # } else 
    if (log_space) {
      res <- sd_sum %>%
        dplyr::filter(stringr::str_detect(item, var_pattern)) %>%
        dplyr::mutate(
          lci = exp(value - 1.96 * se),
          uci = exp(value + 1.96 * se),
          value = exp(value)
        )
      
    # Case 3: Standard natural-scale variables (e.g., spawn_bio, recruits)
    } else {
      res <- sd_sum %>%
        dplyr::filter(stringr::str_detect(item, var_pattern)) %>%
        dplyr::mutate(
          lci = value - 1.96 * se,
          uci = value + 1.96 * se
        )
    }
    
    res %>% 
      dplyr::mutate(
        year = mod$rpt$years, 
        id = id_name
      )
  })
  
  # Render Plot
  ggplot2::ggplot(df, ggplot2::aes(x = year, y = value, color = id, fill = id)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
    ggplot2::geom_line() +
    scico::scale_color_scico_d(name = "Model", palette = palette) +
    scico::scale_fill_scico_d(name = "Model", palette = palette) +
    ggplot2::labs(x = "Year", y = ylab) +
    afscassess::theme_report(base_size = base_size) -> fig

  # Add any extra ggplot layers passed to ... before saving
  if (length(extra_layers) > 0) {
    for (layer in extra_layers) {
      fig = fig + layer
    }
  }

  if (isTRUE(save)) {
    ggplot2::ggsave(
      plot = fig, 
      filename = RTMButils::herein(year, folder, "figs", paste0(sub("\\^", "", var_pattern), "_ribbon.png")),
      width = 6.5, height = 6.5, units = "in", dpi = 200
    )
  }
  
  fig
}

plot_rec_ssb <- function(output, folder = NULL, save = TRUE) {
  rpt = output$rpt
  yrs = rpt$years
  rec_age = rpt$ages[1]
 data.frame(ssb = rpt$spawn_bio[1:(length(yrs) - rec_age)]/1000, 
        rec = rpt$recruits[(rec_age + 1):length(yrs)], 
        year = yrs[1:(length(yrs) - rec_age)]) %>% 
  tidytable::mutate(label = stringr::str_sub(year,3), decade = (floor(year/10) * 10)) %>% 
  ggplot2::ggplot(ggplot2::aes(ssb, rec)) + 
  ggplot2::geom_label(ggplot2::aes(label = label, 
        color = decade), linewidth = 0, show.legend = FALSE, 
        size = 3, family = "Times", alpha = 0.85) + 
  ggplot2::expand_limits(x = 0, 
        y = 0) + 
  scico::scale_color_scico(palette = "romaO") + 
  ggplot2::xlab("SSB (kt)") + 
  ggplot2::ylab("Recruitment (millions)") -> fig
  
  if (isTRUE(save)) {
    ggplot2::ggsave(
      plot = fig, 
      filename = RTMButils::herein(year, folder, "figs", "ssb_rec.png"),
      width = 6.5, height = 6.5, units = "in", dpi = 200)
  }
  
  fig
}

