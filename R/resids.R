#' one-step ahead residuals for RTMB model
#'
#' @param obs observed comp - input sample size and comp weighting adjusted
#' @param pred predicted comp
#' @param yrs comp years
#' @param ind vector or ages of lengths
#' @param label axis label
#' @param outlier mark point as outlier if greater than. default: 3
#' @param addCI boolean to include 95% confidence intervals on the SDNR label, default:TRUE
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' obs = fish_age_obs
#' pred = report1$fish_age_pred
#' yrs = fish_age_yrs
#' label = "Age"
#' ind = ages
#' osa(obs = fish_age_obs, pred = report1$fish_age_pred, iss = data$fish_age_iss,
#'     yrs = fish_age_yrs, ind = ages) +  ggtitle('osa fishery age comp residuals')
#' }
osa <- function(obs, pred, yrs, ind, label = 'Age', outlier=3, addCI = TRUE) {

  res = compResidual::resMulti(obs, pred)
  mat = matrix(res, nrow=nrow(res), ncol=ncol(res))
  df = as.data.frame(mat)
  names(df) <- yrs

  df %>%
    dplyr::mutate(ind = head(ind, -1)) %>%
    tidyr::pivot_longer( -ind, names_to = 'year') %>%
    dplyr::mutate(year = as.numeric(year),
           Year = factor(year),
           id = ifelse(value<0 , 'a', 'b'),
           Outlier = ifelse(abs(value) >= outlier, "Yes", "No"),
           Outlier = factor(Outlier, levels = c('No', 'Yes'))) -> df

  df %>%
    ggplot2::ggplot(ggplot2::aes(year, ind, color = value, size = value, shape = Outlier) ) +
    geom_point(show.legend=TRUE) +
    ggplot2::scale_size_area(guide="none", max_size = 3) +
    scico::scale_color_scico(limits = c(-4, 4), palette = 'vik') +
    tickr::scale_y_tickr(data = df, var = ind, var_min=0) +
    tickr::scale_x_tickr(data = df, var = year) +
    ggplot2::scale_shape_manual(values = c(19,8), drop = FALSE) +
    ggplot2::ylab(label) +
    ggplot2::xlab('Year') +
    ggplot2::ggtitle('OSA') -> osa

  res_vec = as.numeric(res) 
  res_vec = na.omit(res_vec) 
  sdnr_est = sd(res_vec)

  if(addCI) {
    df_n = length(res_vec) - 1 # degrees of freedom
    lci = sqrt(qchisq(0.025, df_n) / df_n) # lower 95% CI
    hci = sqrt(qchisq(0.975, df_n) / df_n) # upper 95% CI
    sdnr_text = sprintf("SDNR = %.2f\nExpected (H0: %.2f - %.2f)", sdnr_est, lci, hci)
  } else {
    sdnr_text = sprintf("SDNR = %.2f", sdnr_est)
  }

  df %>% 
  # dplyr::mutate(label = label) %>%
    # tidyr::pivot_longer(-label) %>%
    ggplot2::ggplot() +
    ggplot2::stat_qq(ggplot2::aes(sample = value)) +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::labs(x = 'Theoretical quantiles', y = 'Sample quantiles') +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = sdnr_text, 
                      hjust = -0.1, vjust = 1.2, size = 4) -> qq
  list(osa=osa, qq=qq)

}

#' pearson residuals for RTMB model
#'
#' @param obs observed comp - input sample size and comp weighting adjusted
#' @param pred predictedcomp
#' @param iss input sample size
#' @param yrs comp years
#' @param wt comp data weighting
#' @param ind vector or ages of lengths
#' @param outlier mark point as outlier if greater than
#' @param label axis label
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pearson(obs = fish_age_obs * iss * wt, pred = report1$fish_age_pred, iss = data$fish_age_iss,
#'     yrs = fish_age_yrs, ind = ages) +  ggtitle('osa fishery age comp residuals')
#' }
pearson <- function(obs, pred, iss, yrs, wt, ind, outlier, label) {

  N = iss * wt
  se = sqrt(sweep(pred * (1 - pred), 2, N, "/"))
  resids = (obs - pred) / se

  df = as.data.frame(resids)
  names(df) <- yrs
  df$ind = ind

  df %>%
    tidyr::pivot_longer( -ind, names_to = 'year') %>%
    dplyr::mutate(year = as.numeric(year),
           Year = factor(year),
           id = ifelse(value<0 , 'a', 'b'),
           Outlier = factor(ifelse(abs(value) >= outlier, "Yes", "No"))) -> df

  df %>%
    ggplot2::ggplot(ggplot2::aes(year, ind, color = value, size = abs(value), shape = Outlier) ) +
    ggplot2::geom_point(show.legend=TRUE) +
    ggplot2::scale_size_area(guide="none", max_size = 3) +
    scico::scale_color_scico(limits = c(-4, 4), palette = 'vik') +
    tickr::scale_y_tickr(data = df, var = ind) +
    tickr::scale_x_tickr(data = df, var = year) +
    ggplot2::scale_shape_manual(values = c(19,8), drop = FALSE) +
    ggplot2::ylab(label) +
    ggplot2::xlab('Year') +
    ggplot2::ggtitle('Pearson')

}

#' aggregate residual plot for RTMB model
#'
#' @param obs observed comp - input sample size and comp weighting adjusted
#' @param pred predicted comp
#' @param ind vector of ages or lengths
#' @param label axis label
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obs = fish_age_obs
#' pred = report1$fish_age_pred
#' ind = ages
#' agg(obs, pred, ind, label = 'Age')
#' }
agg <- function(obs, pred, ind, label) {

  df = data.frame(obs = rowSums(obs)/sum(obs),
                   pred = rowSums(pred)/sum(pred),
                   ind = ind)

  df %>%
    ggplot2::ggplot(ggplot2::aes(ind, pred)) +
    ggplot2::geom_bar(ggplot2::aes(y=obs), stat = 'identity', alpha=0.4) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    tickr::scale_x_tickr(data=df, var=ind, var_min=0) +
    ggplot2::xlab(label) +
    ggplot2::ylab('')

}

#' annual composition plot for RTMB model
#'
#' @param obs observed comp - input sample size and comp weighting adjusted
#' @param pred predicted comp
#' @param ind vector of ages or lengths
#' @param label axis label
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obs = fish_age_obs
#' pred = report1$fish_age_pred
#' ind = ages
#' annual(obs, pred, ind, label = 'Age')
#' }
annual <- function(obs, pred, ind, yrs, label) {
  obs = as.data.frame(obs)
  pred = as.data.frame(pred)
  names(obs) <- names(pred) <- yrs
  obs %>%
    dplyr::mutate(type = 'obs',
                  ind = ind) %>%
    dplyr::bind_rows(
      pred %>%
        dplyr::mutate(type = 'pred',
                      ind = ind)
    ) %>%
    tidyr::pivot_longer(-c(ind, type)) %>%
    dplyr::mutate(year = as.numeric(name)) -> df

  df %>%
    dplyr::filter(type=='obs') %>%
    ggplot2::ggplot(ggplot2::aes(ind, value)) +
    ggplot2::geom_col(alpha = 0.5, ggplot2::aes(fill = factor(ind)), show.legend = FALSE) +
    ggplot2::geom_line(data = dplyr::filter(df, type=='pred')) +
    ggplot2::facet_wrap(~year, ncol = 2, dir='v') +
    tickr::scale_x_tickr(data=df, var=ind, var_min=0) +
    scico::scale_fill_scico_d(palette = 'managua') +
    ggplot2::xlab(label) +
    ggplot2::ylab("") +
    ggplot2::theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

}

effn <- function(obs, pred, iss, wt) {
  colSums((1 - pred) * pred) / colSums((obs - pred)^2)
}

sdnr <- function(obs, pred, iss, wt) {
  n = ncol(obs)
  sdnr = vector(length = n)
  for(i in 1:n) {
    # expected variance for multinomial: N*p*(1-p)
    N = iss[i] * wt
		expected_val = pred[,i] * N
    variance = expected_val * (1 - pred[,i])
    # standardized residual
    std_res = (obs[,i] - expected_val) / sqrt(variance)
    sdnr[i] = sd(std_res, na.rm = TRUE)
  }
  sdnr
}

#' Composition sample size
#'
#' @param obs observed comp data- input sample size and comp weighting adjusted
#' @param pred predicted comp data
#' @param iss input sample size
#' @param yrs years
#' @param wt comp weighting
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sample_size(obs, pred, iss, yrs)}
sample_size <- function(obs, pred, iss, yrs, wt){
  obs_p = sweep(obs, 2, colSums(obs), "/")
  data.frame(year = yrs,
             ISS = iss * wt,
             effN = effn(obs_p, pred, iss, wt),
             sdnr = sdnr(obs_p, pred, iss, wt)) %>%
    tidyr::pivot_longer(-year) %>%
    dplyr::mutate(grp = ifelse(name=='sdnr', 'SDNR', 'Sample size')) -> df

  df %>%
    ggplot2::ggplot(ggplot2::aes(year, value, color = name)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~grp, scales = 'free_y', dir = 'h') +
    ggplot2::scale_color_manual("", breaks = c('effN', 'ISS'), values = c("#7E1700","#5DC0D2",1)) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::theme(legend.position=c(0.1,0.1)) +
    tickr::scale_x_tickr(data=df, var=year, by=10, var_min = 1960) +
    ggplot2::xlab('Year') +
    ggplot2::ylab('Value')
}

#' get all residual plots for RTMB model
#'
#' @param object RTMButils::run_model() output
#' @param var variable name, default "fish_age"
#' @param outlier outlier size, default: 3
#' @param label axis label
#' @param addCI add CI to qq plot
#' @export
#'
#' @examples
#' \dontrun{
#' fish_age_resids <- resids(output = model_1)
#' fish_age_resids$osa + ggtitle('osa fishery age comp residuals')
#' }
resids <- function(object, var = "fish_age", outlier = 3, label = 'Age', addCI = TRUE) {
  rpt = object$rpt
	data = object$dat
	id = if(grepl("age", var)) "ages" else "size"
	ind = data[[id]]
	iss = data[[paste0(var,"_iss")]]
	wt = data[[paste0(var,"_wt")]]
	yrs = data[[paste0(var, "_yrs")]]

  # observed as 'effective counts' (needed for osaand pearson)
  obs_counts = round(iss * data[[paste0(var,"_obs")]] * wt)
  # predicted as proportions (standard for most functions)
  pred_prop = rpt[[paste0(var, "_pred")]] 
  # ensure pred_prop sums to 1 across bins for each year
  pred_prop = sweep(pred_prop, 2, colSums(pred_prop), "/")
  # observed as proportions (annual/aggregate plots)
  obs_prop = sweep(iss * data[[paste0(var,"_obs")]] * wt, 2, colSums(iss * data[[paste0(var,"_obs")]] * wt), "/")
  
  osa_obj = osa(obs_counts, pred_prop, yrs, ind, label, outlier, addCI)
  list(osa = osa_obj$osa,
       qq = osa_obj$qq,
       pearson = pearson(obs_prop, pred_prop, iss=iss, wt = wt, yrs=yrs, ind=ind, outlier=outlier, label = label),
       agg = agg(obs_prop, pred_prop, ind, label=label),
       annual = annual(obs, pred_prop, ind, yrs, label=label),
       ss = sample_size(obs_prop, pred_prop, iss, yrs, wt) )
}


