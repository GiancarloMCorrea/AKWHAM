# Install WHAM package (growth branch):
# remotes::install_github(repo = 'gmoroncorrea/wham', ref='growth', INSTALL_opts = c("--no-docs", "--no-multiarch", "--no-demo"))

rm(list = ls())
require(dplyr)
require(ggplot2)
require(r4ss)
library(wham)
source('aux_fun.R')

# SS files and outputs:
# SS model: 21.1 (2021)
# Read SS model:
data_file = r4ss::SS_readdat_3.30(file = 'C:/Users/moroncog/Documents/StockAssessmentModels/Pacific_cod_2021/Model_21_1/Base_2.dat')
SS_report = r4ss::SS_output(dir = 'C:/Users/moroncog/Documents/StockAssessmentModels/Pacific_cod_2021/Model_21_1') # from OM

# Some model parameters:
mycols = c("#000000", "#E69F00")
n_ages = 20
length_vector = SS_report$lbins
min_year = SS_report$startyr
max_year = SS_report$endyr
n_years = length(min_year:max_year)
LWpars = c(SS_report$parameters$Value[SS_report$parameters$Label == "Wtlen_1_Fem_GP_1"],
           SS_report$parameters$Value[SS_report$parameters$Label == "Wtlen_2_Fem_GP_1"])
GWpars = c(SS_report$parameters$Value[SS_report$parameters$Label == "VonBert_K_Fem_GP_1"],
           SS_report$parameters$Value[SS_report$parameters$Label == "L_at_Amax_Fem_GP_1"],
           SS_report$parameters$Value[SS_report$parameters$Label == "L_at_Amin_Fem_GP_1"], 
           SS_report$parameters$Value[SS_report$parameters$Label == "Richards_Fem_GP_1"],
           SS_report$endgrowth$SD_Beg[2], SS_report$endgrowth$SD_Beg[21])
Q_pars = c(exp(SS_report$parameters$Value[SS_report$parameters$Label == "LnQ_base_Survey(2)"]))
# Selectivity parameters (main):
int_selPos = grep(pattern =  "Size_DblN", x = SS_report$parameters$Label)[1]
SelecParams = SS_report$parameters[int_selPos:nrow(SS_report$parameters), ]
selpars1 = SelecParams$Value[1:6]
selpars2 = SelecParams$Value[7:12]

# NAA info from SS:
NAA_SS = SS_report$natage[SS_report$natage$`Beg/Mid` == 'B' & SS_report$natage$Yr >= min_year & SS_report$natage$Yr <= max_year, 14:(14+n_ages-1)]

# -------------------------------------------------------------------------
# Prepare input data for WHAM:
wham_data = list()
wham_data$ages = 1:n_ages 
wham_data$lengths = length_vector
wham_data$years = as.integer(min_year:max_year)
#Catch information:
wham_data$n_fleets = data_file$Nfleet
wham_data$agg_catch = matrix(data_file$catch$catch[-1], nrow = n_years, ncol = 1)
wham_data$use_agg_catch = matrix(1L, nrow = n_years, ncol = 1)
wham_data$catch_cv = matrix(data_file$catch$catch_se[-1], nrow = n_years, ncol = 1)
# Survey information:
wham_data$n_indices = data_file$Nsurveys
tmp_data = data.frame(year = wham_data$years)
tmp_data2 = merge(tmp_data, data_file$CPUE, by = 'year', all.x = TRUE)
tmp_data2$obs[is.na(tmp_data2$obs)] = 0
tmp_data2$se_log[is.na(tmp_data2$se_log)] = 0
wham_data$agg_indices = tmp_data2$obs
wham_data$index_cv = tmp_data2$se_log
wham_data$units_indices = matrix(0L, nrow = n_years, ncol = 1)
tmp_data = data.frame(year = wham_data$years, use = 1)
tmp_data2 = merge(tmp_data, data_file$CPUE, by = 'year', all.x = TRUE)
tmp_data2$use[is.na(tmp_data2$obs)] = -1
wham_data$use_indices = matrix(tmp_data2$use, nrow = n_years, ncol = 1)
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
# Age comps index:
agecomp_index = data_file$agecomp
tmp_data = data.frame(Yr = wham_data$years)
tmp_data2 = merge(tmp_data, agecomp_index, by = 'Yr', all.x = TRUE)
tmp_data3 = t(apply(tmp_data2[,11:ncol(tmp_data2)], 1, function(x) { x/sum(x) }))
tmp_data3[is.na(tmp_data3)] = 0
tmp_data4 = cbind(tmp_data3, matrix(0, ncol = 8, nrow = n_years))
wham_data$index_paa = as.matrix(tmp_data4)
wham_data$index_Neff = matrix(ifelse(test = is.na(tmp_data2$Nsamp), yes = 0, no = tmp_data2$Nsamp),
                              nrow = n_years, ncol = 1)
wham_data$use_index_paa = matrix(ifelse(test = is.na(tmp_data2$Nsamp), yes = -1, no = 1),
                                 nrow = n_years, ncol = 1)
# Add aging error:
wham_data$index_aging_error = array(NA, dim = c(1,n_ages, n_ages))
wham_data$index_aging_error[1,,] = get_aging_error_matrix(obs_age = SS_report$age_error_mean$type1[2:21],
                                                       sd = SS_report$age_error_sd$type1[2:21])
wham_data$use_index_aging_error = 1
# selectivity and F options
wham_data$selblock_pointer_fleets = matrix(1L, ncol = 1, nrow = n_years)
wham_data$F = matrix(0.2, ncol = 1, nrow = n_years)
wham_data$selblock_pointer_indices = matrix(2L, ncol = 1, nrow = n_years)
wham_data$fracyr_indices = matrix(0.5, ncol = 1, nrow = n_years)
wham_data$fracyr_SSB = matrix(0, ncol = 1, nrow = n_years)
wham_data$age_L1 = 1.5
# WAA information
wham_data$waa_pointer_indices = 1
wham_data$waa_pointer_fleets = 2
wham_data$waa_pointer_totcatch = 1
wham_data$waa_pointer_ssb = 3
wham_data$waa_pointer_jan1 = 3
wham_data$maturity = matrix(rep(SS_report$endgrowth[2:(n_ages+1),18], times = max_year - min_year + 1),
                            ncol = n_ages, nrow = max_year - min_year + 1, byrow = TRUE) 

wham_data$Fbar_ages = 3L:20L
wham_data$percentSPR = 40
wham_data$percentFXSPR = 100
wham_data$percentFMSY = 100
wham_data$XSPR_R_avg_yrs = 1:n_years
wham_data$XSPR_R_opt = 2
wham_data$simulate_period = c(1,0)
wham_data$bias_correct_process = 1
wham_data$bias_correct_observation = 1

# Ecov information for LW parameters(fixed)
env1 = data_file$envdat$Value[data_file$envdat$Variable == 2]
env2 = data_file$envdat$Value[data_file$envdat$Variable == 3]

# Prepare input object:
input = prepare_wham_input(model_name="cod_1",
                               selectivity=list(model = c('len-double-normal', 'len-double-normal'),
                                                re = c('iid', 'iid'),
                                                initial_pars=list(selpars1,
                                                                  selpars2),
                                                fix_pars = list(c(1,2,4,5), c(2,4:6)),
                                                n_selblocks = 2),
                               M = list(model = 'constant', re = 'none',
                                        initial_means = SS_report$Natural_Mortality[1,5],
                                        est_ages = 1),
                               NAA_re = list(sigma="rec", cor = 'iid', N1_model = 0,
                                             recruit_model = 2,
                                             N1_pars = as.vector(as.matrix(NAA_SS[1,])),
                                             recruit_pars = mean(NAA_SS[,1])),
                               growth = list(model = 'Richards',
                                             re = c('none', 'none', 'iid_y', 'none', 'none', 'none'),
                                             init_vals = GWpars,
                                             est_pars = c(1:6)),
                               LW = list(re = c('none', 'none'),
                                     init_vals = LWpars),
                               catchability = list(re = c('none'), 
                                                   initial_q = Q_pars, q_lower = 0,
                                                   q_upper = 1000, prior_sd = NA),
                               age_comp = 'dirichlet-pool0',
                               len_comp = 'dir-mult',
                               basic_info = wham_data)

# Update some input information:
input$par$log_NAA_sigma = log(SS_report$sigma_R_in) # sigma as in SS
input$map$log_NAA_sigma = factor(NA) # fix sigma
input$map$log_N1_pars = factor(rep(NA, times = length(input$par$log_N1_pars))) # fix init NAA
# log_NAA initial values:
input$par$log_NAA = as.matrix(log(NAA_SS)[-1,])
# Fishing mortality values:
input$par$log_F1 = log(0.097) # as in SS
Fts = SS_report$derived_quants[grep(pattern = 'F_', x = SS_report$derived_quants$Label),]
Fts = Fts[1:n_years, 'Value']
F_devs = log(Fts)[-1] - log(Fts)[-n_years]
input$par$F_devs[,1] = F_devs # set F_devs
# LW temporal variability (fixed)
input$par$LW_re[,,1] = rep(log(1 + env1/LWpars[1]), times = n_ages)
input$par$LW_re[,,2] = rep(log(1 + env2/LWpars[2]), times = n_ages)
# Deviations in selectivity parameters: 
SSSelex = SS_report$SelSizeAdj[SS_report$SelSizeAdj$Yr %in% wham_data$years,]
# FISHERY 1:
fleet = 1
tmpSelex = SSSelex[SSSelex$Fleet == fleet, ]
par3 = -log((input$data$selpars_upper[fleet,37]-tmpSelex$Par3)/(tmpSelex$Par3-input$data$selpars_lower[fleet,37]))-input$par$logit_selpars[fleet,37]
par6 = -log((input$data$selpars_upper[fleet,40]-tmpSelex$Par6)/(tmpSelex$Par6-input$data$selpars_lower[fleet,40]))-input$par$logit_selpars[fleet,40]
input$par$selpars_re[1:(n_years*2)] = c(par3, par6)
# INDEX 1:
fleet = 2
tmpSelex = SSSelex[SSSelex$Fleet == fleet, ]
par1 = -log((input$data$selpars_upper[fleet,35]-tmpSelex$Par1)/(tmpSelex$Par1-input$data$selpars_lower[fleet,35]))-input$par$logit_selpars[fleet,35]
par3 = -log((input$data$selpars_upper[fleet,37]-tmpSelex$Par3)/(tmpSelex$Par3-input$data$selpars_lower[fleet,37]))-input$par$logit_selpars[fleet,37]
input$par$selpars_re[(n_years*2+1):(n_years*4)] = c(par1, par3)
# Fix deviates and selectivity parameters
input$map$selpars_re = rep(factor(NA), times = length(input$par$selpars_re)) # fix sel deviates
input$map$sel_repars = rep(factor(NA), times = length(input$par$sel_repars)) # fix selre parameters
input$map$logit_selpars = rep(factor(NA), times = length(input$par$logit_selpars)) # fix sel parameters
# Random effects for L1:
input$random = 'growth_re'

#Run model:
fit_a = fit_wham(input, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE)
save(fit_a, file = 'EBS_pcod/fit_a.RData')

# Make plots
dir.create(path = 'EBS_pcod/fit_a')
plot_wham_output(mod = fit_a, dir.main = 'EBS_pcod/fit_a', out.type = 'pdf')


# -------------------------------------------------------------------------
# ADD MODEL WITH ECOV

# -------------------------------------------------------------------------
# Plot SSB:

# Get SSB from SS:
SS_SSB = SS_report$derived_quants[grep(pattern = 'SSB_', x = SS_report$derived_quants$Label),]
data1 = cbind(data.frame(model = 'StockSynthesis', year = fit_a$years), 
              SS_SSB[3:47, c('Value', 'StdDev')])
colnames(data1)[3:4] = c('est', 'sd')
data1$ssb_min = data1$est - 1.96*data1$sd
data1$ssb_max = data1$est + 1.96*data1$sd

# WHAM model:
this_model = fit_a
model_name = 'WHAM'
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
data2 = cbind(model = model_name, year=1977:2021,
              filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
data2$ssb_min = exp(data2$est-1.96*data2$sd)
data2$ssb_max = exp(data2$est+1.96*data2$sd)
data2$est = exp(data2$est)

# Make plot:
plot_data =rbind(data1, data2)
ggplot(plot_data, aes(year, est, ymin=ssb_min, ymax=ssb_max,
                      fill=model, color=model)) +
  ylim(0,NA) + labs(y='SSB') +
  geom_ribbon(alpha=.3, color = NA) + geom_line(lwd=1) +
  labs( color=NULL, fill=NULL) +
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols) +
  theme_bw() +
  theme(legend.position='top')
ggsave(filename = 'EBS_pcod/compare_SSB.png', width = 190, height = 120, units = 'mm', dpi = 500)

# Plot L1 temporal variability:
