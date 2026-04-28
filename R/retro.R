#' Run a retrospective analysis (peeling) for an RTMB stock assessment model.
#'
#' This function automates the process of retrospective analysis. It iteratively
#' removes recent years of data, refits the model, calculates Mohn's rho,
#' and generates diagnostic plots for specified quantities of interest.
#'
#' @param output Results from RTMButils::run_model
#' @param n_peels The number of retrospective peels to perform (integer)
#' @param year The assessment year, used for constructing file paths (e.g., 2024)
#' @param folder The folder name for the model or type (e.g., "mgmt", "research"), used for file paths (e.g., "mgmt")
#' @param subfolder used for file paths (e.g., "m24"), default: NULL
#' @param quantities A character vector of variable names from the model report
#'   to analyze (e.g., c("spawn_bio", "tot_bio", "recruits"))
#' @param save_outputs A logical flag. If TRUE, saves the retrospective results
#'   (.RDS) and plots (.png) to a 'retro' subfolder
#'
#' @return A list containing three elements:
#'   1. `retro_data`: A tidy data frame with all estimates from all peels.
#'   2. `rho_metrics`: A data frame with Mohn's rho and MAE for each quantity.
#'   3. `plots`: A named list of the generated ggplot objects.
#'
#' @details
#' The function relies on the following packages: `RTMB`, `dplyr`, `purrr`,
#' `ggplot2`, `patchwork`, `scico`, `tibble`, and `here`. Ensure they are installed.
#' The Mohn's rho is calculated as the mean of `(Peel_terminal - Base_terminal) / Base_terminal`.
#'
#' @export
run_retro <- function(output, n_peels = 10, year, folder, subfolder = NULL,
                      quantities = c("spawn_bio", "tot_bio", "recruits"),
                      save_outputs = TRUE) {

  # setup and load base model components
  message("Starting retrospective analysis...")
  retro_path <- herein(year, folder, subfolder, "retro")
  if (save_outputs && !dir.exists(retro_path)) {
    dir.create(retro_path, recursive = TRUE)
  }

  model = output$model
  d0 = output$dat
  years = output$rpt$years
  obj = output$obj
  map = obj$env$map
  fit = output$fit
  lower = fit$lower
  upper = fit$upper
  nms = unique(names(fit$par))
  split_list = split(fit$par, names(fit$par))
  pars = lapply(split_list, unname)
  p0 = pars[nms]

  # put any mapped items back into the pars
  if(!is.null(map)) {
    p0[[names(map)]] = obj$env$parList()[[names(map)]]
  }
  n_years = length(years)

  # run retrospective peels
  message(paste("Running", n_peels, "retrospective peels..."))
  reps = list()
  N_base_sum_catch = sum(d0$catch_ind)

  for (i in 1:n_peels) {
    data = d0
    pars = p0

    peel_indices = (n_years - i + 1):n_years

    # peel data based on the logic in the original script
    data$years = head(d0$years, -i)
    data$catch_ind = head(d0$catch_ind, -i)
    data$catch_obs = head(d0$catch_obs, -i)
    data$catch_wt = head(d0$catch_wt, -i)

    data$srv_ind[peel_indices] <- 0
    if (i + 2 <= n_years) {
      peel_extra_2 = (n_years - (i + 2) + 1):n_years
      data$fish_age_ind[peel_extra_2] <- 0
    } else {
      # if peel removes all years, zero out the whole vector
      data$fish_age_ind[] <- 0
    }

    if (i + 1 <= n_years) {
      peel_extra_1 =(n_years - (i + 1) + 1):n_years
      data$srv_age_ind[peel_extra_1] <- 0
      data$fish_size_ind[peel_extra_1] <- 0
    } else {
      # if peel removes all years, zero out the whole vector
      data$srv_age_ind[] <- 0
      data$fish_size_ind[] <- 0
    }

    # make sure data works
    data = lapply(data, unname)

    # peel parameters
    pars$log_Ft = head(p0$log_Ft, -i)
    pars$log_Rt = head(p0$log_Rt, -i)

    pars = lapply(pars, unname)

    # refit model for the peel
    new_run = run_model(model = model, data = data, pars = pars, map = map, lower=lower, upper=upper)
    # cmb = function(f, d) function(p) f(p, d)
    # obj <- RTMB::MakeADFun(cmb(model, data), pars, map = map, silent = TRUE)
    # fit <- nlminb(start = obj$par,
    #               objective = obj$fn,
    #               gradient = obj$gr,
    #               control = control)
    reps[[paste0('sd', i)]] <- new_run$sd
  }

  if (save_outputs) {
    saveRDS(reps, file.path(retro_path, 'reps.RDS'))
    message(paste("Retrospective fits saved to:", file.path(retro_path, 'reps.RDS')))
  }

  # process outputs
  message("Processing results and calculating Mohn's rho...")
  # get sdreport summary from the base model fit
  base_sds = as.data.frame(summary(output$sd))
  base_sds = tibble::rownames_to_column(base_sds, "item")

  # process the retrospective fits
  purrr::map(reps, summary) %>%
    purrr::map(., as.data.frame) %>%
    purrr::map(., tibble::rownames_to_column, "item") -> retro_sds_list

  # tidy the retrospective data
  dplyr::bind_rows(retro_sds_list, .id = 'id') %>%
    dplyr::mutate(peel = as.numeric(gsub("sd", "", id)),
                  item = gsub("\\..*", "", item)) %>%
    dplyr::filter(item %in% quantities) %>%
    dplyr::select(item, Estimate, Std_Error = `Std. Error`, peel) -> retro_df

  # tidy the base model data
  base_sds %>%
    dplyr::mutate(peel = 0,
                  item = gsub("\\..*", "", item)) %>%
    dplyr::filter(item %in% quantities) %>%
    dplyr::select(item, Estimate, Std_Error = `Std. Error`, peel) -> base_df

  # combine and add years
  dplyr::bind_rows(retro_df, base_df) %>%
    dplyr::group_by(peel, item) %>%
    dplyr::mutate(year = years[1:dplyr::n()]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(lci = Estimate - 1.96 * Std_Error,
                  uci = Estimate + 1.96 * Std_Error,
                  lci = ifelse(lci < 0, 0, lci),
                  peel = factor(peel, levels = as.character(0:n_peels))) -> all_data

  # calculate Mohn's Rho and MAE
  all_data %>%
    dplyr::filter(peel != 0) %>%
    dplyr::group_by(item, peel) %>%
    dplyr::filter(year == max(year)) %>%
    dplyr::summarise(year = max(year),
                     peel_est = Estimate,
                     .groups = 'drop') -> terminal_peels

  all_data %>%
    dplyr::filter(peel == 0) %>%
    dplyr::select(item, year, base_est = Estimate) -> terminal_base

  dplyr::left_join(terminal_peels, terminal_base, by = c("item", "year")) %>%
    dplyr::mutate(pdiff = (peel_est - base_est) / base_est) %>% # Using your formula
    tidyr::drop_na() %>%
    dplyr::group_by(item) %>%
    dplyr::summarise( rho = mean(pdiff, na.rm = TRUE),
                      mae = median(abs(pdiff), na.rm = TRUE)) -> rho_metrics

  # plots
  message("Generating plots...")
  plot_list <- list()

  for (qty in quantities) {
    metrics = dplyr::filter(rho_metrics, item == qty)
    l1 = paste("Mohn's rho =", round(metrics$rho, 2))
    l2 = paste("MAE =", round(metrics$mae, 2))
    plot_title = gsub("_", " ", qty)
    plot_title = paste0(toupper(substring(plot_title, 1, 1)), substring(plot_title, 2))

    plot_data_qty = dplyr::filter(all_data, item == qty)
    annot_x = min(plot_data_qty$year) + 1 # annotation near the start
    annot_y = max(plot_data_qty$uci, na.rm = TRUE) * 0.9

    # plot 1: time series
    p1 <- ggplot2::ggplot(plot_data_qty,
                          ggplot2::aes(x = year, y = Estimate, color = peel, group = peel, fill = peel)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.1, color = NA) +
      ggplot2::geom_line() +
      scico::scale_color_scico_d("Peel", palette = 'roma', direction = -1) +
      scico::scale_fill_scico_d("Peel", palette = 'roma', direction = -1) +
      ggplot2::annotate(geom = 'text',
                        x = annot_x,
                        y = annot_y,
                        label = paste(l1, l2, sep = "\n"),
                        hjust = 0) +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::labs(y = paste(plot_title, "(t)"), x = "Year", title = paste("Retrospective Pattern:", plot_title)) +
      ggplot2::expand_limits(y = 0)

    # plot 2: relative difference
    p2 <- all_data %>%
      dplyr::filter(item == qty, peel != 0) %>%
      dplyr::left_join(
        all_data %>%
          dplyr::filter(item == qty, peel == 0) %>%
          dplyr::select(year, base_est = Estimate),
        by = "year") %>%
      dplyr::mutate(pdiff = (Estimate - base_est) / base_est) %>% # standard relative diff
      tidyr::drop_na() %>%
      ggplot2::ggplot(ggplot2::aes(year,pdiff, color = peel, group = peel)) +
      ggplot2::geom_line(show.legend = FALSE) +
      scico::scale_color_scico_d("Peel", palette = 'roma', direction = -1) +
      ggplot2::geom_hline(yintercept = 0, linetype = 3) +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(y = "Relative Difference from Base", x = "Year")

    combined_plot = patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(3, 1.5)) +
      patchwork::plot_layout(guides = "collect")

    plot_list[[qty]] = combined_plot

    if (save_outputs) {
      file_name = paste0(gsub(" ", "_", tolower(plot_title)), "_retro.png")
      ggplot2::ggsave(
        filename = file.path(retro_path, file_name),
        plot = combined_plot,
        width = 7, height = 7, units = "in", dpi = 300
      )
    }
  }

  message("Done.")
  return(list(
    retro_data = all_data,
    rho_metrics = rho_metrics,
    plots = plot_list
  ))
}

#' Run a prospective analysis (peeling from the start) for an RTMB stock assessment model.
#'
#' @param output The output list from the `run_model()` function.
#' @param n_peels The number of historical years to remove from the start (integer).
#' @param year The assessment year (for file paths).
#' @param folder The folder name (for file paths).
#' @param subfolder Optional subfolder (for file paths).
#' @param quantities Variable names to analyze (e.g., c("spawn_bio", "recruits")).
#' @param save_outputs Logical. If TRUE, saves results and plots.
#'
#' @export
run_prospective <- function(output, n_peels = 5, year, folder, subfolder = NULL,
                            quantities = c("spawn_bio", "tot_bio", "recruits"),
                            save_outputs = TRUE) {

  # setup
  message("Starting prospective analysis...")
  pros_path <- herein(year, folder, subfolder, "prospective")
  if (save_outputs && !dir.exists(pros_path)) {
    dir.create(pros_path, recursive = TRUE)
  }

  # extract components
  model = output$model
  d0 = output$dat
  obj = output$obj
  map = obj$env$map
  fit = output$fit

  # parameter list
  nms = unique(names(fit$par))
  split_list = split(fit$par, names(fit$par))
  p0 = lapply(split_list, unname)[nms]

  # mapped parameters are included in p0
  if(!is.null(map)) {
    p_full = obj$env$parList()
    for(m_name in names(map)) {
      p0[[m_name]] = p_full[[m_name]]
    }
  }

  base_years = output$rpt$years
  n_years = length(base_years)

  # pProspective peels
  message(paste("Running", n_peels, "prospective peels..."))
  reps = list()

  for (i in 1:n_peels) {
    cat("Peeling year:", i, "\n")

    # Copy data and pars
    data_i = d0
    pars_i = p0
    map_i  = map
    # peel data from the START (tail -i)
    # removes the first 'i' years
    data_i$years = tail(d0$years, -i)
    data_i$catch_ind = tail(d0$catch_ind, -i)
    data_i$catch_obs = tail(d0$catch_obs, -i)
    data_i$catch_wt = tail(d0$catch_wt, -i)
    data_i$srv_ind = tail(d0$srv_ind, -i)
    data_i$fish_age_ind = tail(d0$fish_age_ind, -i)
    data_i$srv_age_ind = tail(d0$srv_age_ind, -i)
    data_i$fish_size_ind = tail(d0$fish_size_ind, -i)

    if(!is.null(d0$srv_yrs)) {
      keep_srv <- d0$srv_yrs %in% data_i$years
      data_i$srv_yrs <- d0$srv_yrs[keep_srv]
      data_i$srv_obs <- d0$srv_obs[keep_srv]
      data_i$srv_sd  <- d0$srv_sd[keep_srv]
    }

    if(!is.null(d0$fish_age_yrs)) {
      # Find which age samples are still within the new year range
      keep_fa <- d0$fish_age_yrs %in% data_i$years
      data_i$fish_age_yrs <- d0$fish_age_yrs[keep_fa]
      data_i$fish_age_iss <- d0$fish_age_iss[keep_fa]
      data_i$fish_age_obs <- d0$fish_age_obs[, keep_fa] # Slice columns
    }
    if(!is.null(d0$srv_age_yrs)) {
      keep_sa <- d0$srv_age_yrs %in% data_i$years
      data_i$srv_age_yrs <- d0$srv_age_yrs[keep_sa]
      data_i$srv_age_iss <- d0$srv_age_iss[keep_sa]
      data_i$srv_age_obs <- d0$srv_age_obs[, keep_sa]
    }
    if(!is.null(d0$fish_size_yrs)) {
      keep_fs <- d0$fish_size_yrs %in% data_i$years
      data_i$fish_size_yrs <- d0$fish_size_yrs[keep_fs]
      data_i$fish_size_iss <- d0$fish_size_iss[keep_fs]
      data_i$fish_size_obs <- d0$fish_size_obs[, keep_fs]
    }

    # parameters from the START
    pars_i$log_Ft = tail(p0$log_Ft, -i)
    pars_i$log_Rt = tail(p0$log_Rt, -i)

    # peel map from the START
    if ("log_Ft" %in% names(map_i)) map_i$log_Ft = tail(map_i$log_Ft, -i)
    if ("log_Rt" %in% names(map_i)) map_i$log_Rt = tail(map_i$log_Rt, -i)

    # Refit model
    # Note: bounds (lower/upper) might need slicing if they are vectors!
    new_run <- run_model(model = model, data = data_i, pars = pars_i,
                         map = map_i, lower = fit$lower, upper = fit$upper)

    reps[[paste0('sd', i)]] = new_run$sd
  }

  # tidying
  message("Processing results...")

  # base model data (peel 0)
  as.data.frame(summary(output$sd)) %>%
    tibble::rownames_to_column("item") %>%
    dplyr::mutate(peel = 0, item = gsub("\\..*", "", item)) %>%
    dplyr::filter(item %in% quantities) %>%
    dplyr::select(item, Estimate, Std_Error = `Std. Error`, peel) -> base_sds

  # peel
  retro_sds_list <- purrr::map(reps, function(x) {
    as.data.frame(summary(x)) %>%
      tibble::rownames_to_column("item")
  })

  dplyr::bind_rows(retro_sds_list, .id = 'id') %>%
    dplyr::mutate(
      peel = as.numeric(gsub("sd", "", id)),
      item = gsub("\\..*", "", item)
    ) %>%
    dplyr::filter(item %in% quantities) %>%
    dplyr::select(item, Estimate, Std_Error = `Std. Error`, peel) -> retro_df

  # combine
  dplyr::bind_rows(retro_df, base_df = base_sds) %>%
    dplyr::group_by(peel, item) %>%
    # For prospective, the first year of the peel is 'peel + 1' index of base_years
    dplyr::mutate(year = base_years[(unique(peel) + 1):n_years]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      lci = pmax(0, Estimate - 1.96 * Std_Error),
      uci = Estimate + 1.96 * Std_Error,
      peel = factor(peel)) -> all_data

  # plotting
  message("Generating plots...")
  plot_list = list()

  for (qty in quantities) {
    plot_data_qty = all_data %>% dplyr::filter(item == qty)
    plot_title = gsub("_", " ", qty)
    plot_title = paste0(toupper(substring(plot_title, 1, 1)), substring(plot_title, 2))

    p1 <- ggplot2::ggplot(plot_data_qty, ggplot2::aes(x = year, y = Estimate, color = peel, group = peel)) +
      ggplot2::geom_line() +
      scico::scale_color_scico_d("Years Removed\nfrom Start", palette = 'roma') +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::labs(y = plot_title, x = "Year", title = paste("Prospective Analysis:", plot_title))

    p2_data <- plot_data_qty %>%
      dplyr::filter(peel != 0) %>%
      dplyr::left_join(plot_data_qty %>%
                         dplyr::filter(peel == 0) %>%
                         dplyr::select(year, base_est = Estimate), by = "year") %>%
      dplyr::mutate(pdiff = (Estimate - base_est) / base_est)

    p2 <- ggplot2::ggplot(p2_data, ggplot2::aes(x = year, y = pdiff, color = peel, group = peel)) +
      ggplot2::geom_line(show.legend = FALSE) +
      scico::scale_color_scico_d(palette = 'roma') +
      ggplot2::geom_hline(yintercept = 0, linetype = 3) +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(y = "% Diff from Base", x = "Year")

    combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(3, 1))
    plot_list[[qty]] <- combined_plot

    if(save_outputs) {
      ggplot2::ggsave(file.path(pros_path, paste0(qty, "_prospective.png")), combined_plot, width = 8, height = 8)
    }
  }

  if (save_outputs) saveRDS(all_data, file.path(pros_path, 'prospective_data.RDS'))

  return(list(prospect_data = all_data, plots = plot_list))
}
