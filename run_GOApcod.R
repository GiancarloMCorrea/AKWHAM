# Install WHAM package (growth branch):
# remotes::install_github(repo = 'GiancarloMCorrea/wham', ref='growth', INSTALL_opts = c("--no-docs", "--no-multiarch", "--no-demo"))

# Set you WD:
setwd("~/GitHub/AKWHAM")

rm(list = ls())
require(dplyr)
require(ggplot2)
require(r4ss)
library(wham)
source('aux_fun.R')

# Load some files:

# SS files and outputs:
# SS model can be found here: 
# https://github.com/pete-hulson/goa_pcod/tree/master/2022/Stock_Synthesis_files/Model19.1a%20(22)%20-%20wADFG
# Change Month column from 1 to 7 for surveys (len comps and CAAL) 
# Fix par2 (selex) for fishery 1
data_file = r4ss::SS_readdat_3.30(file = 'SS_models/GOA_pcod/GOAPcod2022Oct25_wADFG.dat')
SS_report = r4ss::SS_output(dir = 'SS_models/GOA_pcod') # from OM

# Some model parameters:
n_ages = 10
length_vector = SS_report$lbins
min_year = SS_report$startyr
max_year = SS_report$endyr
n_years = length(min_year:max_year)
LWpars = c(SS_report$Growth_Parameters$WtLen1, SS_report$Growth_Parameters$WtLen2)
GWpars = c(SS_report$Growth_Parameters$K, SS_report$Growth_Parameters$Linf, 
           SS_report$endgrowth$Len_Beg[SS_report$endgrowth$Real_Age == 1], # This is L1
           SS_report$endgrowth$SD_Beg[2], SS_report$endgrowth$SD_Beg[11])
Q_pars = c(exp(SS_report$parameters$Value[SS_report$parameters$Label == "LnQ_base_Srv(4)"]),
           exp(SS_report$parameters$Value[SS_report$parameters$Label == "LnQ_base_LLSrv(5)"]))
LAA_SS = as.matrix(SS_report$growthseries[1,6:15])
# Selectivity parameters (main):
int_selPos = grep(pattern =  "Size_DblN", x = SS_report$parameters$Label)[1]
SelecParams = SS_report$parameters[int_selPos:nrow(SS_report$parameters), ]
selpars1 = SelecParams$Value[1:6]
selpars1[5] = -19 # modify initial value
selpars2 = SelecParams$Value[7:12]
selpars2[5] = -19 # modify initial value
selpars3 = SelecParams$Value[13:18]
selpars3[5] = -19 # modify initial value
selpars4 = SelecParams$Value[19:24]
selpars5 = SelecParams$Value[25:30]
selpars5[5] = -19 # modify initial value

# Abundance-at-age info from SS:
NAA_SS = SS_report$natage[SS_report$natage$`Beg/Mid` == 'B' & SS_report$natage$Yr >= min_year & SS_report$natage$Yr <= max_year, 14:(14+n_ages-1)]

# -------------------------------------------------------------------------
# Prepare input data for WHAM:
wham_data = list()
wham_data$ages = 1:n_ages 
wham_data$lengths = length_vector
wham_data$years = as.integer(min_year:max_year)
#Catch information:
wham_data$n_fleets = data_file$Nfleet
wham_data$agg_catch = matrix(data_file$catch$catch[data_file$catch$year > 0], nrow = n_years, ncol = data_file$Nfleet)
wham_data$use_agg_catch = matrix(1L, nrow = n_years, ncol = data_file$Nfleet)
wham_data$use_agg_catch[which(wham_data$agg_catch == 0)] = -1
wham_data$catch_cv = matrix(data_file$catch$catch_se[data_file$catch$year > 0], nrow = n_years, ncol = data_file$Nfleet)
# Survey information:
tmp_data <- data_file$CPUE[data_file$CPUE$index > 0,]
tmp_data2 <- full_join(expand.grid(year = wham_data$years, index = unique(tmp_data$index)),
                       tmp_data, by = c('year', 'index'))
wham_data$n_indices = length(unique(tmp_data2$index))
wham_data$agg_indices = matrix(tmp_data2$obs, nrow = n_years, ncol = wham_data$n_indices)
wham_data$agg_indices[is.na(wham_data$agg_indices)] = 0
wham_data$index_cv = matrix(tmp_data2$se_log, nrow = n_years, ncol = wham_data$n_indices)
wham_data$index_cv[is.na(wham_data$index_cv)] = 0
wham_data$units_indices = matrix(0L, nrow = n_years, ncol = wham_data$n_indices)
wham_data$use_indices = matrix(1L, nrow = n_years, ncol = wham_data$n_indices)
wham_data$use_indices[is.na(tmp_data2$obs)] = -1
# Turn off age comps for fishery and surveys (following WHAM philosophy)
wham_data$use_catch_paa = matrix(0L, nrow = n_years, ncol = wham_data$n_fleets)
wham_data$use_index_paa = matrix(0L, nrow = n_years, ncol = wham_data$n_indices)
# Len comps catch:
wham_lencomps = array(0, dim = c(wham_data$n_fleets, n_years, length(length_vector)))
wham_lenNeff = matrix(0, ncol = wham_data$n_fleets, nrow = n_years)
wham_lenuse = matrix(-1, ncol = wham_data$n_fleets, nrow = n_years)
for(j in 1:wham_data$n_fleets) {
  lencomp_fleet = data_file$lencomp[data_file$lencomp$FltSvy == j, ]
  lencomp_fleet2 = as.matrix(lencomp_fleet[,7:ncol(lencomp_fleet)])
  lencomp_temp = as.matrix(lencomp_fleet2)
  wham_lencomps[j,match(lencomp_fleet$Yr, wham_data$years),] = lencomp_temp/rowSums(lencomp_temp)
  wham_lenNeff[match(lencomp_fleet$Yr, wham_data$years),j] = lencomp_fleet$Nsamp
  wham_lenuse[match(lencomp_fleet$Yr, wham_data$years),j] = 1
}
wham_data$catch_pal = wham_lencomps
wham_data$catch_NeffL = wham_lenNeff
wham_data$use_catch_pal = wham_lenuse
# Len comps index:
wham_lencomps = array(0, dim = c(wham_data$n_indices, n_years, length(length_vector)))
wham_lenNeff = matrix(0, ncol = wham_data$n_indices, nrow = n_years)
wham_lenuse = matrix(-1, ncol = wham_data$n_indices, nrow = n_years)
for(j in 1:wham_data$n_indices) {
  lencomp_fleet = data_file$lencomp[data_file$lencomp$FltSvy == j + wham_data$n_fleets, ]
  lencomp_fleet2 = as.matrix(lencomp_fleet[,7:ncol(lencomp_fleet)])
  lencomp_temp = as.matrix(lencomp_fleet2)
  wham_lencomps[j,match(lencomp_fleet$Yr, wham_data$years),] = lencomp_temp/rowSums(lencomp_temp)
  wham_lenNeff[match(lencomp_fleet$Yr, wham_data$years),j] = lencomp_fleet$Nsamp
  wham_lenuse[match(lencomp_fleet$Yr, wham_data$years),j] = 1
}
wham_data$index_pal = wham_lencomps
wham_data$index_NeffL = wham_lenNeff
wham_data$use_index_pal = wham_lenuse
# CAAL catch (catch 1,2,3):
wham_data$catch_caal = array(NA, dim = c(wham_data$n_fleets, length(wham_data$years), 
                                         length(wham_data$lengths), length(wham_data$ages)))
wham_data$catch_caal_Neff = array(0, dim = c(length(wham_data$years), wham_data$n_fleets, 
                                              length(wham_data$lengths)))
wham_data$use_catch_caal =  array(-1, dim = c(length(wham_data$years), wham_data$n_fleets, 
                                              length(wham_data$lengths)))
out_caal = get_caal_from_SS(caal_SSdata = data_file$agecomp, fleet = 1, model_years = wham_data$years, 
                            model_lengths = wham_data$lengths, model_ages = wham_data$ages)
wham_data$catch_caal[1,,,] = out_caal[[1]]
wham_data$catch_caal_Neff[,1,] = out_caal[[2]]
wham_data$use_catch_caal[,1,] = out_caal[[3]]
out_caal = get_caal_from_SS(caal_SSdata = data_file$agecomp, fleet = 2, model_years = wham_data$years, 
                            model_lengths = wham_data$lengths, model_ages = wham_data$ages)
wham_data$catch_caal[2,,,] = out_caal[[1]]
wham_data$catch_caal_Neff[,2,] = out_caal[[2]]
wham_data$use_catch_caal[,2,] = out_caal[[3]]
out_caal = get_caal_from_SS(caal_SSdata = data_file$agecomp, fleet = 3, model_years = wham_data$years, 
                            model_lengths = wham_data$lengths, model_ages = wham_data$ages)
wham_data$catch_caal[3,,,] = out_caal[[1]]
wham_data$catch_caal_Neff[,3,] = out_caal[[2]]
wham_data$use_catch_caal[,3,] = out_caal[[3]]
# CAAL index (index 4):
wham_data$index_caal = array(NA, dim = c(wham_data$n_indices, length(wham_data$years), 
                                         length(wham_data$lengths), length(wham_data$ages)))
wham_data$index_caal_Neff = array(0, dim = c(length(wham_data$years), wham_data$n_indices, 
                                              length(wham_data$lengths)))
wham_data$use_index_caal = array(-1, dim = c(length(wham_data$years), wham_data$n_indices, 
                                             length(wham_data$lengths)))
out_caal = get_caal_from_SS(caal_SSdata = data_file$agecomp, fleet = 4, model_years = wham_data$years, 
                            model_lengths = wham_data$lengths, model_ages = wham_data$ages)
wham_data$index_caal[1,,,] = out_caal[[1]]
wham_data$index_caal_Neff[,1,] = out_caal[[2]]
wham_data$use_index_caal[,1,] = out_caal[[3]]
# Add aging error:
base_aging_error = get_aging_error_matrix(obs_age = SS_report$age_error_mean$type1[2:(n_ages+1)],
                                          sd = SS_report$age_error_sd$type1[2:(n_ages+1)])
wham_data$index_aging_error = array(NA, dim = c(wham_data$n_indices,n_ages, n_ages))
wham_data$index_aging_error[1,,] = base_aging_error
wham_data$index_aging_error[2,,] = base_aging_error
wham_data$use_index_aging_error = rep(1, times = wham_data$n_indices)
wham_data$catch_aging_error = array(NA, dim = c(wham_data$n_fleets,n_ages, n_ages))
wham_data$catch_aging_error[1,,] = base_aging_error
wham_data$catch_aging_error[2,,] = base_aging_error
wham_data$catch_aging_error[3,,] = base_aging_error
wham_data$use_catch_aging_error = rep(1, times = wham_data$n_fleets)
# selectivity and F options
wham_data$selblock_pointer_fleets = matrix(rep(c(1,2,3), each = n_years), ncol = wham_data$n_fleets, nrow = n_years)
wham_data$F = matrix(0.2, ncol = wham_data$n_fleets, nrow = n_years)
wham_data$selblock_pointer_indices = matrix(rep(c(4,5), each = n_years), ncol = wham_data$n_indices, nrow = n_years)
wham_data$fracyr_indices = matrix(0.5, ncol = wham_data$n_indices, nrow = n_years)
wham_data$fracyr_SSB = matrix(0, ncol = 1, nrow = n_years)
# WAA information 
wham_data$waa_pointer_indices = rep(1, times = wham_data$n_indices)
wham_data$waa_pointer_fleets = rep(2, times = wham_data$n_fleets)
wham_data$waa_pointer_totcatch = 2
wham_data$waa_pointer_ssb = 3
wham_data$waa_pointer_jan1 = 3
wham_data$maturity = matrix(rep(SS_report$endgrowth[2:(n_ages+1),18], times = max_year - min_year + 1),
                            ncol = n_ages, nrow = max_year - min_year + 1, byrow = TRUE) 
wham_data$Fbar_ages = 1L:10L
wham_data$percentSPR = 60
wham_data$percentFXSPR = 100
wham_data$percentFMSY = 100
wham_data$XSPR_R_avg_yrs = 1:n_years
wham_data$XSPR_R_opt = 2
wham_data$simulate_period = c(1,0)
wham_data$bias_correct_process = 1
wham_data$bias_correct_observation = 1

# Ecov information (for Q - index 2):
env1 = numeric(n_years)
env1[match(data_file$envdat$Yr, wham_data$years)] = data_file$envdat$Value
env1 = data_file$envdat$Value
n_env_years = length(data_file$envdat$Yr)

ecov <- list(
  label = c("env1"),
  mean = matrix(c(env1), ncol = 1),
  logsigma = matrix(log(0.01), ncol = 1, nrow = n_env_years), # sigma = 0.2
  year = data_file$envdat$Yr,
  use_obs = matrix(1L, ncol=1, nrow=n_env_years),
  lag = list(rep(0, times = 8)),
  ages = list(1:n_ages),
  process_model = c('ar1'),
  where = list('q'),
  indices = list(2),
  how = c(1))


# -------------------------------------------------------------------------
# Prepare input object: use vB
input_a = prepare_wham_input(model_name="goa_cod_1",
                               selectivity=list(model = rep('len-double-normal', times = 5),
                                                re = c(rep('iid', times = 4), 'none'),
                                                initial_pars=list(selpars1, selpars2, selpars3,
                                                                  selpars4, selpars5),
                                                fix_pars = list(5:6,4:6,4:6,NULL,1:6),
                                                n_selblocks = 5),
                               M = list(model = 'constant', re = 'none',
                                        initial_means = SS_report$Natural_Mortality_endyr[1,5],
                                        est_ages = 1),
                               NAA_re = list(sigma="rec", cor = 'iid', N1_model = 1,
                                             recruit_model = 2,
                                             N1_pars = c(NAA_SS[1,1], 0),
                                             #N1_pars = as.vector(as.matrix(NAA_SS[1,])),
                                             recruit_pars = mean(NAA_SS[,1])),
                               growth = list(model = 'vB_classic',
                                             re = c('none', 'none', 'none'),
                                             init_vals = GWpars[1:3],
                                             est_pars = 1:3,
                                             SD_vals = GWpars[4:5],
                                             SD_est = 1:2),
                               LW = list(re = c('none', 'none'),
                                     init_vals = LWpars),
                               catchability = list(re = c('none', 'none'),
                                                   initial_q = Q_pars,
                                                   q_lower = c(0,0),
                                                   q_upper = c(10,10), prior_sd = c(NA,NA)),
                               ecov = ecov,
                               basic_info = wham_data)

# update some inputs as SS model:
input_a = post_input_GOApcod(input_a, SS_report, NAA_SS)
# no random effects:
input_a$random <- c("Ecov_re")

#Run model:
fit_a = fit_wham(input_a, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton = 0)
check_convergence(fit_a)
save(fit_a, file = 'GOA_pcod/fit_a.RData')

# Make plots
dir.create(path = 'GOA_pcod/fit_a')
plot_wham_output(mod = fit_a, dir.main = 'GOA_pcod/fit_a', out.type = 'png')
