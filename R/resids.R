#' one-step ahead residuals for RTMB model
#'
#' @param obs observed comp
#' @param pred predictedcomp
#' @param iss input sample size
#' @param yrs comp years
#' @param ind vector or ages of lengths
#' @param label axis label
#' @param outlier mark point as outlier if greater than. default: 3
#'
#' @export
#'
#' @examples
#' \dontrun{
#' iss = data$fish_age_iss
#' obs = fish_age_obs
#' pred = report1$fish_age_pred
#' yrs = fish_age_yrs
#' label = "Age"
#' ind = ages
#' osa(obs = fish_age_obs, pred = report1$fish_age_pred, iss = data$fish_age_iss,
#'     yrs = fish_age_yrs, ind = ages) +  ggtitle('osa fishery age comp residuals')
#' }
osa <- function(obs, pred, iss, yrs, ind, label = 'Age', outlier) {

  obs = round(iss * obs / colSums(obs))
  pred = pred / colSums(pred)
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
    ggplot2::scale_size_area(guide=F, max_size = 3) +
    scico::scale_color_scico(limits = c(-4, 4), palette = 'vik') +
    afscassess::scale_y_tickr(data = df, var = ind, start=0) +
    afscassess::scale_x_tickr(data = df, var = year) +
    ggplot2::scale_shape_manual(values = c(19,8), drop = FALSE) +
    ggplot2::ylab(label) +
    ggplot2::xlab('Year') +
    ggplot2::ggtitle('OSA') -> osa

  df %>%
    ggplot2::ggplot() +
    ggplot2::stat_qq(ggplot2::aes(sample=value)) +
    ggplot2::geom_abline(slope=1, intercept = 0, lty=3) +
    ggplot2::labs(x = 'Theoretical quantiles', y = 'Sample quantiles') -> qq
  list(osa=osa, qq=qq)

}

#' pearson residuals for RTMB model
#'
#' @param obs observed comp
#' @param pred predictedcomp
#' @param iss input sample size
#' @param yrs comp years
#' @param ind vector or ages of lengths
#' @param label axis label
#' @param outlier mark point as outlier if greater than. eg., 3
#'
#' @export
#'
#' @examples
#' \dontrun{
#' iss = data$fish_age_iss
#' obs = fish_age_obs
#' pred = report1$fish_age_pred
#' yrs = fish_age_yrs
#' label = "Age"
#' ind = ages
#' pearson(obs = fish_age_obs, pred = report1$fish_age_pred, iss = data$fish_age_iss,
#'     yrs = fish_age_yrs, ind = ages) +  ggtitle('osa fishery age comp residuals')
#' }
pearson <- function(obs, pred, iss, yrs, ind, outlier, label) {

  as.data.frame(iss * (obs -pred) / sqrt(iss * pred)) %>%
    dplyr::mutate(ind = ind) -> df
  names(df) <- c(yrs, 'ind')

  df %>%
    tidyr::pivot_longer( -ind, names_to = 'year') %>%
    dplyr::mutate(year = as.numeric(year),
           Year = factor(year),
           id = ifelse(value<0 , 'a', 'b'),
           Outlier = factor(ifelse(abs(value) >= outlier, "Yes", "No"))) -> df

  df %>%
    ggplot2::ggplot(ggplot2::aes(year, ind, color = value, size = value, shape = Outlier) ) +
    ggplot2::geom_point(show.legend=TRUE) +
    ggplot2::scale_size_area(guide=F, max_size = 3) +
    scico::scale_color_scico(limits = c(-4, 4), palette = 'vik') +
    afscassess::scale_y_tickr(data = df, var = ind) +
    afscassess::scale_x_tickr(data = df, var = year) +
    ggplot2::scale_shape_manual(values = c(19,8), drop = FALSE) +
    ggplot2::ylab(label) +
    ggplot2::xlab('Year') +
    ggplot2::ggtitle('Pearson')

}

#' one-step ahead residuals qq plot for RTMB model
#'
#' @param obs observed comp
#' @param pred predicted comp
#' @param iss input sample size
#' @param yrs comp years
#' @param label axis label
#'
#' @export
#'
#' @examples
#' \dontrun{
#' iss = data$fish_age_iss
#' obs = fish_age_obs
#' pred = report1$fish_age_pred
#' yrs = fish_age_yrs
#' label = "Age"
#' osa(obs = fish_age_obs, pred = report1$fish_age_pred, iss = data$fish_age_iss,
#'     yrs = fish_age_yrs, ind = ages) +  ggtitle('osa fishery age comp residuals')
#' }
qq <- function(obs, pred, iss, yrs, label) {
  obs = round(iss * obs / colSums(obs))
  pred = pred / colSums(pred)
  res = compResidual::resMulti((obs), (pred))
  mat = matrix(res, nrow=nrow(res), ncol=ncol(res))
  df = as.data.frame(mat)
  names(df) <- yrs

  dplyr::mutate(df, label = label) %>%
    tidyr::pivot_longer(-label) %>%
    ggplot2::ggplot() +
    ggplot2::stat_qq(ggplot2::aes(sample = value)) +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::labs(x = 'Theoretical quantiles', y = 'Sample quantiles')
}

#' aggregate residual plot for RTMB model
#'
#' @param obs observed comp
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
    afscassess::scale_x_tickr(data=df, var=ind, start=0) +
    ggplot2::xlab(label) +
    ggplot2::ylab('')

}

#' annual composition plot for RTMB model
#'
#' @param obs observed comp
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
    ggplot2::geom_line(data = filter(df, type=='pred')) +
    ggplot2::facet_wrap(~year, ncol = 2, dir='v') +
    afscassess::scale_x_tickr(data=df, var=ind, start=0) +
    scico::scale_fill_scico_d(palette = 'managua') +
    ggplot2::xlab(label) +
    ggplot2::ylab("") +
    ggplot2::theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

}

effn <- function(obs, pred) {
  colSums((1 - pred) * pred) / colSums((obs - pred)^2)
}

sdnr <- function(obs, pred, iss) {
  n = ncol(obs)
  sdnr = vector(length = n)
  for(i in 1:n) {

    sdnr[i] = sd((obs[,i] - pred[,i]) / sqrt(pred[,i] * (1 - pred[,i]) / iss[i]))
  }
  sdnr
}

#' Composition sample size
#'
#' @param obs observed comp data
#' @param pred predicted comp data
#' @param iss input sample size
#' @param yrs years
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sample_size(obs, pred, iss, yrs)}
sample_size <- function(obs, pred, iss, yrs){
  data.frame(year = yrs,
             ISS = iss,
             effN = effn(obs, pred),
             sdnr = sdnr(obs, pred, iss)) %>%
    tidyr::pivot_longer(-year) %>%
    dplyr::mutate(grp = ifelse(name=='sdnr', 'SDNR', 'Sample size')) -> df

  df %>%
    ggplot2::ggplot(ggplot2::aes(year, value, color = name)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~grp, scales = 'free_y', dir = 'h') +
    ggplot2::scale_color_manual("", breaks = c('effN', 'ISS'), values = c("#7E1700","#5DC0D2",1)) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::theme(legend.position=c(0.1,0.8)) +
    afscassess::scale_x_tickr(data=df, var=year, to=10, start = 1960) +
    ggplot2::xlab('Year') +
    ggplot2::ylab('Value')
}

#' get all residual plots for RTMB model
#'
#' @param obs observed comp
#' @param pred predictedcomp
#' @param iss input sample size
#' @param outlier outlier size, default: 3
#' @param yrs comp years
#' @param ind vector or ages of lengths
#' @param label axis label
#'
#' @export
#'
#' @examples
#' \dontrun{
#' iss = data$fish_age_iss
#' obs = fish_age_obs
#' pred = report1$fish_age_pred
#' yrs = fish_age_yrs
#' label = "Age"
#' ind = ages
#' fish_age_resids <- resids(obs = fish_age_obs, pred = report1$fish_age_pred, iss = data$fish_age_iss,
#'                          yrs = fish_age_yrs, ind = ages)
#'  fish_age_resids$osa + ggtitle('osa fishery age comp residuals')
#' }
resids <- function(obs, pred, iss, outlier = 3, yrs, ind, label = 'Age') {
  list(osa = osa(obs, pred, iss, yrs, ind, label, outlier)$osa,
       qq = osa(obs, pred, iss, yrs, ind, label, outlier)$qq,
       pearson = pearson(obs, pred, iss, yrs, ind, label, outlier),
       agg = agg(obs, pred, ind, label),
       annual = annual(obs, pred, ind, yrs, label),
       ss = sample_size(obs, pred, iss, yrs) )
}

