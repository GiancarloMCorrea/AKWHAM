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
#' @param minyr,maxyr The min/max year (or cohort) to plot
#' @return Invisibly a data frame of WAA data and results
#'
plot_waa_fit <- function(fit, plot=TRUE, by.cohort=TRUE, minyr=1900, maxyr=2100, sizeX = 10){
  require(ggplot2);require(dplyr)
  ind <- fit$input$data$waa_pointer_ssb
  nages = length(fit$ages.lab)
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
    waa <- group_by(waa, cohort) %>%
      filter(n()>2, cohort>minyr, cohort<maxyr) %>% ungroup
    g0 <- ggplot(waa, aes(as.numeric(age), obs, ymin=lwr, ymax=upr))+
      facet_wrap('cohort')+
      scale_x_continuous(name="Age", breaks=1:nages) 
  } else {
    waa <- filter(waa, year>minyr, year<maxyr)
    g0 <- ggplot(waa, aes(year, obs, ymin=lwr, ymax=upr))+
      facet_wrap('age', scales='free_y') +
      theme(axis.text=element_text(size=sizeX))
  }
  g <- g0+
    geom_ribbon(aes(ymin=ymin,ymax=ymax), alpha=.3, fill=2)+
    geom_line(aes(y=exp), col=2, lwd=1)+
    geom_pointrange(fatten=2) +
    labs(y='Weight (kg)', x=NULL) 
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
  use_matrix = matrix(-1, nrow = length(model_years), ncol = length(model_lengths))
  
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
    use_matrix[match(data_years[i], model_years),match(tmp_data$Lbin_hi, model_lengths)] = 1
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
  #input$par$Ecov_re = base_input$par$Ecov_re # how this impacts the model?
  # Selectivity
  input$data$selpars_lower[,13:16] = base_input$data$selpars_lower[,13:16]
  input$data$selpars_upper[,13:16] = base_input$data$selpars_upper[,13:16]
  input$par$logit_selpars[,1:16] = base_input$par$logit_selpars
  #input$map$logit_selpars = factor(rep(NA, times = length(input$map$logit_selpars))) # fix parameters
  input$par$selpars_re[1:104] = base_input$par$selpars_re[1:104]
  input$map$selpars_re = factor(c(1:104, rep(NA, 104)))
  #input$map$selpars_re = factor(rep(NA, times = length(input$par$selpars_re))) # fix deviates
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

post_input_GOApcod = function(input, SS_report, NAA_SS) {
  
  years = input$years
  n_years = length(years)
  n_ages = input$data$n_ages

  input$par$log_NAA_sigma = log(SS_report$sigma_R_in) # sigma as in SS
  input$map$log_NAA_sigma = factor(NA) # fix sigma
  input$map$log_N1_pars = factor(c(1,NA))
  # log_NAA initial values:
  input$par$log_NAA = as.matrix(log(NAA_SS)[-1,])
  # Fishing mortality values:
  F_matrix = as.matrix(SS_report$timeseries[SS_report$timeseries$Yr %in% years,grep(pattern = 'F:_', x = colnames(SS_report$timeseries))])
  small_F = 0.0001 # small number F1 for fishery 3
  input$par$log_F1 = c(log(F_matrix[1,1]),log(F_matrix[1,2]),log(small_F)) 
  input$map$log_F1 = factor(c(1,2,NA)) # fix last F
  F_devs = log(F_matrix)[-1,] - log(F_matrix)[-n_years,] # only for fishery 1 and 2
  F_devs[which(is.nan(F_devs))] = 0
  F_devs[10,3] = log(F_matrix[11,3]) - log(small_F)
  input$par$F_devs = F_devs # set F_devs
  input$map$F_devs = factor(c(1:((n_years-1)*2), rep(NA, times = 9), 91:126))
  # Add time block for M 2014-2016:
  input$par$M_re = matrix(rep(log(SS_report$Z_at_age$`0`[SS_report$Z_at_age$Yr %in% wham_data$years]) - log(SS_report$Natural_Mortality[1,5]), times = n_ages), ncol = n_ages)
  tmpMmatrix = matrix(NA, ncol= n_ages, nrow = n_years)
  tmpMmatrix[years %in% 2014:2016] = 1
  input$par$M_repars[1] = log(0.5)
  input$map$M_re = factor(as.vector(tmpMmatrix))
  # Deviations in selectivity parameters (initial values): 
  SSSelex = SS_report$SelSizeAdj[SS_report$SelSizeAdj$Yr %in% years,]
  # FISHERY 1:
  fleet = 1
  tmpSelex = SSSelex[SSSelex$Fleet == fleet, ]
  par1 = -log((input$data$selpars_upper[fleet,25]-tmpSelex$Par1)/(tmpSelex$Par1-input$data$selpars_lower[fleet,25]))-input$par$logit_selpars[fleet,25]
  par2 = -log((input$data$selpars_upper[fleet,26]-tmpSelex$Par2)/(tmpSelex$Par2-input$data$selpars_lower[fleet,26]))-input$par$logit_selpars[fleet,26]
  par3 = -log((input$data$selpars_upper[fleet,27]-tmpSelex$Par3)/(tmpSelex$Par3-input$data$selpars_lower[fleet,27]))-input$par$logit_selpars[fleet,27]
  par4 = -log((input$data$selpars_upper[fleet,28]-tmpSelex$Par4)/(tmpSelex$Par4-input$data$selpars_lower[fleet,28]))-input$par$logit_selpars[fleet,28]
  input$par$selpars_re[1:(n_years*4)] = c(par1, par2, par3, par4)
  # FISHERY 2:
  fleet = 2
  tmpSelex = SSSelex[SSSelex$Fleet == fleet, ]
  par1 = -log((input$data$selpars_upper[fleet,25]-tmpSelex$Par1)/(tmpSelex$Par1-input$data$selpars_lower[fleet,25]))-input$par$logit_selpars[fleet,25]
  par2 = -log((input$data$selpars_upper[fleet,26]-tmpSelex$Par2)/(tmpSelex$Par2-input$data$selpars_lower[fleet,26]))-input$par$logit_selpars[fleet,26]
  par3 = -log((input$data$selpars_upper[fleet,27]-tmpSelex$Par3)/(tmpSelex$Par3-input$data$selpars_lower[fleet,27]))-input$par$logit_selpars[fleet,27]
  par4 = -log((input$data$selpars_upper[fleet,28]-tmpSelex$Par4)/(tmpSelex$Par4-input$data$selpars_lower[fleet,28]))-input$par$logit_selpars[fleet,28]
  input$par$selpars_re[(n_years*4+1):(n_years*8)] = c(par1, par2, par3, par4)
  # FISHERY 3:
  fleet = 3
  tmpSelex = SSSelex[SSSelex$Fleet == fleet, ]
  par1 = -log((input$data$selpars_upper[fleet,25]-tmpSelex$Par1)/(tmpSelex$Par1-input$data$selpars_lower[fleet,25]))-input$par$logit_selpars[fleet,25]
  par2 = -log((input$data$selpars_upper[fleet,26]-tmpSelex$Par2)/(tmpSelex$Par2-input$data$selpars_lower[fleet,26]))-input$par$logit_selpars[fleet,26]
  par3 = -log((input$data$selpars_upper[fleet,27]-tmpSelex$Par3)/(tmpSelex$Par3-input$data$selpars_lower[fleet,27]))-input$par$logit_selpars[fleet,27]
  input$par$selpars_re[(n_years*8+1):(n_years*11)] = c(par1, par2, par3)
  # INDEX 1:
  fleet = 4
  tmpSelex = SSSelex[SSSelex$Fleet == fleet, ]
  par1 = -log((input$data$selpars_upper[fleet,25]-tmpSelex$Par1)/(tmpSelex$Par1-input$data$selpars_lower[fleet,25]))-input$par$logit_selpars[fleet,25]
  par2 = -log((input$data$selpars_upper[fleet,26]-tmpSelex$Par2)/(tmpSelex$Par2-input$data$selpars_lower[fleet,26]))-input$par$logit_selpars[fleet,26]
  par3 = -log((input$data$selpars_upper[fleet,27]-tmpSelex$Par3)/(tmpSelex$Par3-input$data$selpars_lower[fleet,27]))-input$par$logit_selpars[fleet,27]
  par4 = -log((input$data$selpars_upper[fleet,28]-tmpSelex$Par4)/(tmpSelex$Par4-input$data$selpars_lower[fleet,28]))-input$par$logit_selpars[fleet,28]
  par5 = -log((input$data$selpars_upper[fleet,29]-tmpSelex$Par5)/(tmpSelex$Par5-input$data$selpars_lower[fleet,29]))-input$par$logit_selpars[fleet,29]
  par6 = -log((input$data$selpars_upper[fleet,30]-tmpSelex$Par6)/(tmpSelex$Par6-input$data$selpars_lower[fleet,30]))-input$par$logit_selpars[fleet,30]
  input$par$selpars_re[(n_years*11+1):(n_years*17)] = c(par1, par2, par3, par4, par5, par6)
  # Selectivity blocks/deviates (mapping):
  # Fishery 1:
  map_f1_par1 = c(1:13, rep(14, times = 15), rep(15, times = 2), rep(16, times = 10), rep(NA, times = 6))
  map_f1_par2 = c(rep(NA, times = 13), rep(17, times = 15), rep(18, times = 2), rep(19, times = 10), rep(20, times = 6))
  map_f1_par3 = c(21:33, rep(34, times = 15), rep(35, times = 2), rep(36, times = 10), rep(37, times = 6))
  map_f1_par4 = c(38:50, rep(51, times = 15), rep(52, times = 2), rep(53, times = 10), rep(54, times = 6))
  # Fishery 2:
  map_f2_par1 = c(NA, 55:66, rep(67, times = 15), rep(68, times = 2), rep(69, times = 10), rep(70, times = 6))
  map_f2_par2 = c(rep(NA, times = 13), rep(71, times = 15), rep(72, times = 2), rep(73, times = 10), rep(74, times = 6))
  map_f2_par3 = c(NA, 75:86, rep(87, times = 15), rep(88, times = 2), rep(89, times = 10), rep(90, times = 6))
  map_f2_par4 = c(rep(NA, times = 13), rep(NA, times = 15), rep(NA, times = 2), rep(NA, times = 10), rep(NA, times = 6))
  # Fishery 3:
  map_f3_par1 = c(rep(NA, times = 40), rep(91, times = 6))
  map_f3_par2 = c(rep(NA, times = 40), rep(92, times = 6))
  map_f3_par3 = c(rep(NA, times = 40), rep(93, times = 6))
  # Index 1:
  map_i1_par1 = c(rep(NA, times = 19), rep(94, times = 10), rep(95, times = 17))
  map_i1_par2 = c(rep(NA, times = 19), rep(96, times = 10), rep(97, times = 17))
  map_i1_par3 = c(rep(NA, times = 19), rep(98, times = 10), rep(99, times = 17))
  map_i1_par4 = c(rep(NA, times = 19), rep(100, times = 10), rep(101, times = 17))
  map_i1_par5 = c(rep(NA, times = 19), rep(102, times = 10), rep(103, times = 17))
  map_i1_par6 = c(rep(NA, times = 19), rep(NA, times = 10), rep(NA, times = 17))
  # Now merge all vectors:
  input$map$selpars_re = factor(c(map_f1_par1,map_f1_par2,map_f1_par3,map_f1_par4,
                           map_f2_par1,map_f2_par2,map_f2_par3,map_f2_par4,
                           map_f3_par1,map_f3_par2,map_f3_par3,
                           map_i1_par1,map_i1_par2,map_i1_par3,map_i1_par4,map_i1_par5,map_i1_par6))
  #input$map$selpars_re = factor(rep(NA, times = length(input$map$selpars_re)))
  # Fix process error for selectivity:
  input$map$sel_repars = factor(rep(NA, times = length(input$map$sel_repars)))
  # Fix selectivity parameters as in SS
  #input$map$logit_selpars = factor(c(rep(NA, times = 120), 1:15, 16, NA, 17:19, NA, NA, NA, 20, NA, NA, NA, 21:23))
  input$map$logit_selpars = factor(rep(NA, times = length(input$map$logit_selpars)))
  # Fix process error for Ecov:
  input$map$Ecov_process_pars = factor(rep(NA, times = length(input$map$Ecov_process_pars)))

  return(input)
}

post_input_EBSpcod = function(input, SS_report, NAA_SS) {
  
  years = input$years
  n_years = length(years)
  n_ages = input$data$n_ages

  # Update some input information:
  input$par$log_NAA_sigma = log(SS_report$sigma_R_in) # sigma as in SS
  input$map$log_NAA_sigma = factor(NA) # fix sigma
  # log_NAA initial values:
  input$par$log_NAA = as.matrix(log(NAA_SS)[-1,])
  # Fishing mortality values:
  input$par$log_F1 = log(0.18) # as in SS
  Fts = SS_report$derived_quants[grep(pattern = 'F_', x = SS_report$derived_quants$Label),]
  Fts = Fts[1:n_years, 'Value']
  F_devs = log(Fts)[-1] - log(Fts)[-n_years]
  input$par$F_devs[,1] = F_devs # set F_devs
  # Deviations in selectivity parameters: 
  SSSelex = SS_report$SelSizeAdj[SS_report$SelSizeAdj$Yr %in% wham_data$years,]
  ncolSelex = ncol(input$data$selpars_upper)
  # FISHERY 1:
  fleet = 1
  tmpSelex = SSSelex[SSSelex$Fleet == fleet, ]
  par3 = -log((input$data$selpars_upper[fleet,ncolSelex-3]-tmpSelex$Par3)/(tmpSelex$Par3-input$data$selpars_lower[fleet,ncolSelex-3]))-input$par$logit_selpars[fleet,ncolSelex-3]
  par6 = -log((input$data$selpars_upper[fleet,ncolSelex]-tmpSelex$Par6)/(tmpSelex$Par6-input$data$selpars_lower[fleet,ncolSelex]))-input$par$logit_selpars[fleet,ncolSelex]
  input$par$selpars_re[1:(n_years*2)] = c(par3, par6)
  # INDEX 1:
  fleet = 2
  tmpSelex = SSSelex[SSSelex$Fleet == fleet, ]
  par1 = -log((input$data$selpars_upper[fleet,ncolSelex-5]-tmpSelex$Par1)/(tmpSelex$Par1-input$data$selpars_lower[fleet,ncolSelex-5]))-input$par$logit_selpars[fleet,ncolSelex-5]
  par3b = -log((input$data$selpars_upper[fleet,ncolSelex-3]-tmpSelex$Par3)/(tmpSelex$Par3-input$data$selpars_lower[fleet,ncolSelex-3]))-input$par$logit_selpars[fleet,ncolSelex-3]
  input$par$selpars_re[(n_years*2+1):(n_years*4)] = c(par1, par3b)
  # Selectivity blocks/deviates (mapping):
  # Fishery 1:
  map_f1_par3 = 1:n_years
  map_f1_par6 = 1:n_years + n_years
  # Index 1:
  map_f2_par1 = c(rep(NA, times = 5), 93:133)
  map_f2_par3 = c(rep(NA, times = 5), 134:174)
  # Now merge all vectors:
  #input$map$selpars_re = factor(c(map_f1_par3, map_f1_par6, map_f2_par1, map_f2_par3))
  input$map$selpars_re = factor(rep(NA, times = length(input$par$selpars_re)))
  # Fix process error for selectivity:
  #input$par$sel_repars[,1] = log(0.5) # increase sigma selex
  input$map$sel_repars = factor(rep(NA, times = length(input$map$sel_repars)))
  # Fix selectivity parameters as in SS
  #input$map$logit_selpars = factor(c(rep(NA, times = 68), 1:3, NA, 4:5, rep(NA, times = 4), 6, NA)) 
  input$map$logit_selpars = factor(rep(NA, times = length(input$map$logit_selpars)))

  return(input)
}

plot_ecov_fit <- function(mod, label = ""){

  require(ggplot2);require(dplyr)

  dat = mod$env$data
  years_full <- mod$years

  ecov.pred = mod$rep$Ecov_x
  ecov.obs = dat$Ecov_obs[1:dat$n_years_Ecov,,drop=F]
  # ecov.obs.sig = mod$rep$Ecov_obs_sigma # Ecov_obs_sigma now a derived quantity in sdrep
  if(class(mod$sdrep)[1] == "sdreport"){
    sdrep = summary(mod$sdrep)
  } else {
    sdrep = mod$sdrep
  }

  ecov.obs.sig = mod$rep$Ecov_obs_sigma 
  ecov.use = dat$Ecov_use_obs[1:dat$n_years_Ecov,,drop=F]
  ecov.obs.sig = ecov.obs.sig[1:dat$n_years_Ecov,,drop=F]
  ecov.obs.sig[ecov.use == 0] <- NA
  ecov.pred.se = matrix(sdrep[rownames(sdrep) %in% "Ecov_x",2], ncol=dat$n_Ecov)

  # default: don't plot the padded entries that weren't used in ecov likelihood
  ecov.res = (ecov.obs - ecov.pred[1:dat$n_years_Ecov,]) / ecov.obs.sig # standard residual (obs - pred)
  ecovs <- 1:dat$n_Ecov

  plot_dat = data.frame(years = years_full, obs = ecov.obs[,1], se = ecov.obs.sig[,1])
  plot_dat$lwr = plot_dat$obs - 1.96 * plot_dat$se
  plot_dat$upr = plot_dat$obs + 1.96 * plot_dat$se
  plot_dat$exp = ecov.pred[,1]
  plot_dat$exp_se = ecov.pred.se[,1]
  plot_dat$ymin = plot_dat$exp - 1.96 * plot_dat$exp_se
  plot_dat$ymax = plot_dat$exp + 1.96 * plot_dat$exp_se

  g0 <- ggplot(plot_dat, aes(as.numeric(years), obs, ymin=lwr, ymax=upr))+
    geom_ribbon(aes(ymin=ymin,ymax=ymax), alpha=.3, fill=2)+
    geom_line(aes(y=exp), col=2, lwd=1)+
    geom_pointrange(fatten=2) +
    theme_bw()+
    annotate("text", label = label, x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
    labs(y='Temperature Bering10K index', x=NULL) 

  return(g0)
}

plot_data_overview <- function(datlist, sectionCex = 1, fleetCex = 0.5){
  dd <- datlist
  x1 <- with(dd, data.frame(year=styr:endyr, size=cattot, survey='fishery', type='catch'))
  x2 <- with(dd, data.frame(year=fshyrs, size=multN_fsh, survey='fishery', type='ages'))
  x3 <- with(dd, data.frame(year=srvyrs1, size=indxsurv_log_sd1, survey='Shelikof', type='age 3+ index'))
  x4 <- with(dd, data.frame(year=srv_acyrs1, size=multN_srv1, survey='Shelikof', type='ages'))
  x5 <- with(dd, data.frame(year=srvyrs2, size=indxsurv_log_sd2, survey='NMFS BT', type='index'))
  x6 <- with(dd, data.frame(year=srv_acyrs2, size=multN_srv2, survey='NMFS BT', type='ages'))
  x7 <- with(dd, data.frame(year=srvyrs3, size=indxsurv_log_sd3, survey='ADF&G', type='index'))
  x8 <- with(dd, data.frame(year=srv_acyrs3, size=multN_srv3, survey='ADF&G', type='ages'))
  x9 <- with(dd, data.frame(year=srvyrs4, size=indxsurv_log_sd4, survey='Shelikof', type='age 1 index'))
  x10 <- with(dd, data.frame(year=srvyrs5, size=indxsurv_log_sd5, survey='Shelikof', type='age 2 index'))
  x11 <- with(dd, data.frame(year=srvyrs6, size=indxsurv_log_sd6, survey='Summer AT', type='index'))
  x12 <- with(dd, data.frame(year=srv_acyrs6, size=multN_srv6, survey='Summer AT', type='ages'))
  ## lengths special cases
  #x13 <- with(dd, data.frame(year=srv_lenyrs2, size=multNlen_srv2, survey='NMFS BT', type='lengths'))
  #x14 <- with(dd, data.frame(year=srv_lenyrs6, size=multNlen_srv6, survey='Summer AT', type='lengths'))
  dat <- rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)
  dat <- group_by(dat, survey, type) %>% mutate(relsize=size/max(size)) %>% ungroup
  dat <- mutate(dat, id=paste(survey, type, sep='_'))
  ## png('data.png', width=7, height=5, units='in', res=400)
  size.cex <- 1
  maxsize <- 2
  ymax <- 19
  xmid <- mean(unique(dat$year))
  par(mar=c(1.75,.75,2,.75), mgp=c(1.5,.35,0), tck=-.01)
  plot(0, xlim = c(min(dat$year), dd$endyr+16.5), ylim = c(-1, ymax+1), axes = FALSE,  yaxs = "i",
       type = "n", xlab = NA, ylab = "" )
  title('GOA Walleye pollock')
  box()
  mycircles <- function(survey, type,y, color, lab=type){
    text(x=dd$endyr+2, y=ymax-y, label=lab, pos=4, cex = 0.85)
    xx <- dat[dat$survey==survey &  dat$type==type & dat$size>0,]
    if(nrow(xx)>0){
      symbols(x=xx$year, y=rep(ymax-y, length(xx$year)), circles=sqrt(xx$relsize),
              bg = adjustcolor(color, alpha.f=.6),
              add = TRUE, inches = .04)
    }
  }
  cols <- c(rgb(127,201,127,max=256),rgb(190,174,212,max=256),rgb(253,192,134,max=256),rgb(255,255,153,max=256),rgb(56,108,176, max=256))
  text(x=xmid, y=ymax-.15, labels='Fishery', font=2, cex=sectionCex)
  mycircles(survey='fishery', type='catch', y=1, color=cols[1], lab='Catch')
  mycircles(survey='fishery', type='ages', y=2, color=cols[1], lab='Age Comps')
  text(x=xmid, y=ymax-3.15, labels='Shelikof', font=2, cex=sectionCex)
  mycircles(survey='Shelikof', type='ages', y=4, color=cols[2], lab='Age comps')
  mycircles(survey='Shelikof', type='age 3+ index', y=5, color=cols[2], lab='Age 3+ index')
  mycircles(survey='Shelikof', type='age 1 index', y=6, color=cols[2], lab='Age 1 index')
  mycircles(survey='Shelikof', type='age 2 index', y=7, color=cols[2], lab='Age 2 index')
  text(x=xmid, y=ymax-8.15, labels='Summer AT', font=2, cex=sectionCex)
  mycircles(survey='Summer AT', type='index', y=9, color=cols[3], lab='Index')
  mycircles(survey='Summer AT', type='ages', y=10, color=cols[3], lab='Age comps')
  #mycircles(survey='Summer AT', type='lengths', y=11, color=cols[3], lab='Length comps')
  text(x=xmid, y=ymax-12.15, labels='NMFS BT', font=2, cex=sectionCex)
  mycircles(survey='NMFS BT', type='index', y=13, color=cols[4], lab='Index')
  mycircles(survey='NMFS BT', type='ages', y=14, color=cols[4], lab='Age comps')
  #mycircles(survey='NMFS BT', type='lengths', y=15, color=cols[4], lab='Length comps')
  text(x=xmid, y=ymax-16.15, labels='ADF&G BT', font=2, cex=sectionCex)
  mycircles(survey='ADF&G', type='index', y=17, color=cols[5], lab='Index')
  mycircles(survey='ADF&G', type='ages', y=18, color=cols[5], lab='Age comps')
                                       #abline(v=dd$endyr, type=3, col=gray(.5))
  axis(1, at=seq(1970,dd$endyr, by=5))
  return(invisible(dat))
}
