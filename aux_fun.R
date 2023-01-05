require(dplyr)

#' Plot WAA residuals as bubble plots
#' @param fit from wham
#'
plot_waa_resids <- function(fit){
  ind <- fit$input$data$waa_pointer_ssb
  pears <- (log(fit$rep$pred_waa[ind,,])-log(fit$input$data$waa[ind,,]))/sqrt(log(fit$input$data$waa_cv[ind,,]^2+1))
  pears[!is.finite(pears)] <- NA
  pears <-
    data.frame(age=rep(seq_along(fit$ages), each=length(fit$years)),
               year=rep(fit$years, times=length(fit$ages.lab)),
               residual=as.numeric(pears))
  ggplot(pears, aes(year, size=abs(residual), color=residual<0, y=age)) + geom_point()
}


#' Plot weight-at-age matrix estimates compared to data
#' @param fit Fitted WHAM model that does not use empirical WAA
#' @param plot Whether to plot (default) or not
#' @param by.cohort Whether to plot grouped by cohort (default) or
#'   by age
#' @return Invisibly a data frame of WAA data and results
#'
plot_waa_fit <- function(fit, plot=TRUE, by.cohort=TRUE){
  require(ggplot2);require(dplyr)
  ind <- fit$input$data$waa_pointer_ssb
  ## get SE for SSB matrix
  se <- summary(fit$sdrep, select='all')[,2]
  se <- as.numeric(array(se[names(se)=='pred_waa'],
                         dim=dim(fit$input$data$waa))[ind,,])
  if(length(se)==0)
    stop("No ADREPORT variable for pred_waa found, update WHAM")
  o <- as.numeric(fit$input$data$waa[ind,,])
  e <- as.numeric(fit$rep$pred_waa[ind,,])
  CV <- as.numeric(fit$input$data$waa_cv[2,,])
  lwr <- o/exp(1.96*sqrt(log(1+((o*CV)/o))^2))
  upr <- o*exp(1.96*sqrt(log(1+((o*CV)/o))^2))
  waa <- data.frame(age=rep(fit$ages.lab, each=length(fit$years)),
                    year=rep(fit$years, times=length(fit$ages.lab)),
                    exp=e, obs=o, CV=CV, lwr=lwr, upr=upr,
                    ymin=e-1.96*se, ymax=e+1.96*se) %>%
    mutate(obs=ifelse(CV==0, NA,obs),
           lwr=ifelse(CV==0, NA,lwr),
           upr=ifelse(CV==0, NA,upr)) %>%
    mutate(age=factor(age, levels=fit$ages.lab),
           cohort=year-as.numeric(age))
  if(by.cohort){
    waa <- group_by(waa, cohort) %>% filter(n()>2) %>% ungroup
    g0 <- ggplot(waa, aes(as.numeric(age), obs, ymin=lwr, ymax=upr))+
      facet_wrap('cohort')
  } else {
    g0 <- ggplot(waa, aes(year, obs, ymin=lwr, ymax=upr))+
      facet_wrap('age', scales='free_y')
  }
  g <- g0+
    geom_ribbon(aes(ymin=ymin,ymax=ymax), alpha=.3, fill=2)+
    geom_line(aes(y=exp), col=2, lwd=1)+
    geom_pointrange(fatten=2) +
    labs(y='weight (kg)', x=NULL)
  if(plot) print(g)
  return(invisible(waa))
}


get_caal_from_SS = function(caal_SSdata, fleet, model_years, model_lengths, model_ages) {

  # Filter data and remove age0 if present
  this_data = caal_SSdata[caal_SSdata$FltSvy == fleet,]
  if(colnames(this_data)[10] == 'a0') this_data = this_data[-10]

  #Create output objects
  caal_array = array(0, dim = c(length(model_years), length(model_lengths), length(model_ages)))
  Neff_matrix = matrix(0, nrow = length(model_years), ncol = length(model_lengths))
  use_matrix = matrix(-1, nrow = length(model_years), ncol = 1)

  # Some relevant information:
  data_years = unique(this_data$Yr)
  len_bin = model_lengths[2] - model_lengths[1]
  age_names = colnames(this_data[,10:ncol(this_data)])
  data_ages = as.numeric(gsub(pattern = 'a', replacement = '', x = age_names))

  # New length bins:
  this_data$Lbin_hi = as.numeric(as.character(cut(x = this_data$Lbin_lo,
      breaks = seq(from = model_lengths[1] - len_bin*0.5, to = model_lengths[length(model_lengths)] + len_bin*0.5, by = len_bin),
      labels = model_lengths)))

  # Standardize data (sum = 0) and then multiply by Nsamps:
  this_data[,10:ncol(this_data)] = (this_data[,10:ncol(this_data)]/rowSums(this_data[,10:ncol(this_data)]))*this_data$Nsamp

  # Sum information by new length bins
  prop_data = this_data %>%
              dplyr::group_by(Yr, Lbin_hi) %>%
              dplyr::summarise(across(c(age_names[1]:age_names[length(age_names)]), ~ mean(.x)))
  eff_data = this_data %>%
                dplyr::group_by(Yr, Lbin_hi) %>%
                dplyr::summarise(Neff = sum(Nsamp))


  for(i in seq_along(data_years)) {
    tmp_data = prop_data[prop_data$Yr == data_years[i], ]
    tmp_eff_data = eff_data[eff_data$Yr == data_years[i], ]
    caal_array[match(data_years[i], model_years), match(tmp_data$Lbin_hi, model_lengths), match(data_ages, model_ages)] = as.matrix(tmp_data[,3:ncol(tmp_data)]/rowSums(tmp_data[,3:ncol(tmp_data)]))
    Neff_matrix[match(data_years[i], model_years),match(tmp_data$Lbin_hi, model_lengths)] = tmp_eff_data$Neff
    use_matrix[match(data_years[i], model_years),1] = 1
  }

  output = list(caal = caal_array, Neff = Neff_matrix, use = use_matrix)
  return(output)

}

post_input_pollock = function(input, base_input) {


  # NAA information
  input$par$log_NAA = as.matrix(base_input$par$log_NAA)
  input$map$log_N1_pars = base_input$map$log_N1_pars
  input$map$log_NAA_sigma = base_input$map$log_NAA_sigma
  # F information
  input$par$F_devs = base_input$par$F_devs
  input$par$log_F1 = base_input$par$log_F1
  # Q information
  input$par$logit_q = base_input$par$logit_q
  input$par$q_re = base_input$par$q_re
  input$map$q_repars = base_input$map$q_repars
  input$par$q_repars = base_input$par$q_repars
  input$data$use_q_prior = base_input$data$use_q_prior
  input$data$logit_q_prior_sigma = base_input$data$logit_q_prior_sigma
  input$par$q_prior_re = base_input$par$q_prior_re
  # data agg index sigma
  input$data$agg_index_sigma = base_input$data$agg_index_sigma
  # Ecov
  input$par$Ecov_re = base_input$par$Ecov_re
  # Selectivity
  input$data$selpars_lower[,13:16] = base_input$data$selpars_lower[,13:16]
  input$data$selpars_upper[,13:16] = base_input$data$selpars_upper[,13:16]
  input$par$logit_selpars[,1:16] = base_input$par$logit_selpars
  input$par$selpars_re[1:104] = base_input$par$selpars_re[1:104]
  input$map$selpars_re = factor(c(1:104, rep(NA, 104)))
  input$map$sel_repars = base_input$map$sel_repars

  return(input)
}

get_aging_error_matrix = function(obs_age, sd) {

  out_matrix = matrix(NA, ncol = length(obs_age), nrow = length(obs_age))

  for(i in seq_along(obs_age)) {

    for(j in seq_along(obs_age)) {

      # if(i > 1) {
      #   if((j == 1) | (j == length(obs_age))) {
      #     out_matrix[j,i] = 1 - pnorm(q = (j-obs_age[i])/sd[i])
      #   }
      #   if((j > 1) & (j < length(obs_age))) {
      #     out_matrix[j,i] = pnorm(q = (j+1-obs_age[i])/sd[i]) - pnorm(q = (j-obs_age[i])/sd[i])
      #   }
      #} else {
        if(j == length(obs_age)) {
          out_matrix[j,i] = 1 - pnorm(q = (j-obs_age[i])/sd[i])
        } else {
          out_matrix[j,i] = pnorm(q = (j+1-obs_age[i])/sd[i]) - pnorm(q = (j-obs_age[i])/sd[i])
        }
      #}

    }

  }

  return(out_matrix)

}
