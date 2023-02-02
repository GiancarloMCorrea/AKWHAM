# Install WHAM package (growth branch):
# remotes::install_github(repo = 'gmoroncorrea/wham', ref='growth', INSTALL_opts = c("--no-docs", "--no-multiarch", "--no-demo"))

# Set you WD:
setwd("~/GitHub/AKWHAM")

rm(list = ls())
require(dplyr)
require(ggplot2)
require(r4ss)
library(wham)
source('aux_fun.R')

# SS files and outputs:
# SS model can be found here: 
# https://github.com/afsc-assessments/EBS_PCOD/tree/main/2022_ASSESSMENT/NOVEMBER_MODELS/NEW_MODELS/Model19_12A
data_file = r4ss::SS_readdat_3.30(file = 'SS_models/EBS_pcod/BSPcod22_OCT.dat')
SS_report = r4ss::SS_output(dir = 'SS_models/EBS_pcod') # from OM

# Some model parameters:
mycols = c("#D16103", "#52854C", "#4E84C4")
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
wham_data$waa_pointer_totcatch = 2
wham_data$waa_pointer_ssb = 3
wham_data$waa_pointer_jan1 = 3
wham_data$maturity = matrix(rep(SS_report$endgrowth[2:(n_ages+1),18], times = max_year - min_year + 1),
                            ncol = n_ages, nrow = max_year - min_year + 1, byrow = TRUE) 
wham_data$Fbar_ages = 1L:20L
wham_data$percentSPR = 40
wham_data$percentFXSPR = 100
wham_data$percentFMSY = 100
wham_data$XSPR_R_avg_yrs = 1:n_years
wham_data$XSPR_R_opt = 2
wham_data$simulate_period = c(1,0)
wham_data$bias_correct_process = 1
wham_data$bias_correct_observation = 1


# -------------------------------------------------------------------------
# Model with iid_y L1
# Ecov information (for L1):
env1 = read.csv('aux_data/sebs_summer_bottom_temp.csv')
env2 = env1[env1$year %in% (min_year):(max_year) & env1$GCM_scen == 'gfdl_ssp126', ] 
# no lag (model does better and also correlation is higher)
env2$stand_index = (env2$mn_val - mean(env2$mn_val))/sd(env2$mn_val)

ecov <- list(
  label = c("Bering10K"),
  mean = matrix(env2$stand_index, ncol = 1),
  logsigma = matrix(log(0.2), ncol = 1, nrow = n_years), # sigma = 0.2
  #logsigma = 'est_1', # estimate sigma. WHAM estimates Ecov_obs_sigma = 0.598 for all years
  year = min_year:max_year,
  use_obs = matrix(1L, ncol=1, nrow=n_years),
  lag = list(rep(0, times = 7)),
  ages = list(1:n_ages),
  process_model = c('ar1'),
  where = list('growth'),
  where_subindex = 3, # on L1
  how = c(0))

# Prepare input object:
input_a = prepare_wham_input(model_name="ebs_cod_1",
                               selectivity=list(model = c('len-double-normal', 'len-double-normal'),
                                                #re = c('none', 'none'),
                                                re = c('iid', 'iid'),
                                                initial_pars=list(selpars1, selpars2),
                                                #fix_pars = list(1:6, 1:6),
                                                fix_pars = list(c(1,2,4,5), c(2,4:6)),
                                                n_selblocks = 2),
                               M = list(model = 'constant', re = 'none',
                                        initial_means = SS_report$Natural_Mortality[1,5],
                                        est_ages = 1),
                               NAA_re = list(sigma="rec", cor = 'iid', N1_model = 1,
                                             recruit_model = 2,
                                             #N1_pars = as.vector(as.matrix(NAA_SS[1,])),
                                             N1_pars = c(NAA_SS[1,1], 0.12),
                                             recruit_pars = mean(NAA_SS[,1])),
                               growth = list(model = 'Richards',
                                             re = c('none', 'none', 'ar1_y', 'none'),
                                             init_vals = GWpars[1:4],
                                             est_pars = c(1:4),
                                             SD_vals = GWpars[5:6],
                                             SD_est = 1:2),
                               LW = list(re = c('none', 'none'),
                                     init_vals = LWpars),
                               catchability = list(re = c('none'), 
                                                   initial_q = Q_pars, q_lower = 0,
                                                   q_upper = 1000, prior_sd = NA),
                               age_comp = 'dirichlet-pool0',
                               ecov = ecov,
                               basic_info = wham_data)

# update some inputs as SS model:
input_a = post_input_EBSpcod(input_a, SS_report, NAA_SS)
input_a$map$Ecov_beta = factor(rep(NA, times = length(input_a$map$Ecov_beta))) # no effect on L1
input_a$random = c("growth_re", "Ecov_re")
# input_a$random = NULL
# input_a$map$growth_repars = factor(rep(NA, times = length(input_a$map$growth_repars)))

#Run model:
fit_a = fit_wham(input_a, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton = 0)
check_convergence(fit_a)
save(fit_a, file = 'EBS_pcod/fit_a.RData')

# Make plots
dir.create(path = 'EBS_pcod/fit_a')
plot_wham_output(mod = fit_a, dir.main = 'EBS_pcod/fit_a', out.type = 'pdf')

# Get SSB estimates:
this_model = fit_a
model_name = 'AR(1)_WHAM'
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBdata1 = cbind(model = model_name, year=1977:2022,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
SSBdata1$ssb_min = exp(SSBdata1$est-1.96*SSBdata1$sd)
SSBdata1$ssb_max = exp(SSBdata1$est+1.96*SSBdata1$sd)
SSBdata1$est = exp(SSBdata1$est)

# -------------------------------------------------------------------------
# Model with Ecov L1

# Prepare input object:
input_b = prepare_wham_input(model_name="ebs_cod_2",
                             selectivity=list(model = c('len-double-normal', 'len-double-normal'),
                                              #re = c('none', 'none'),
                                              re = c('iid', 'iid'),
                                              initial_pars=list(selpars1, selpars2),
                                              #fix_pars = list(1:6, 1:6),
                                              fix_pars = list(c(1,2,4,5), c(2,4:6)),
                                              n_selblocks = 2),
                             M = list(model = 'constant', re = 'none',
                                      initial_means = SS_report$Natural_Mortality[1,5],
                                      est_ages = 1),
                             NAA_re = list(sigma="rec", cor = 'iid', N1_model = 1,
                                           recruit_model = 2,
                                           #N1_pars = as.vector(as.matrix(NAA_SS[1,])),
                                           N1_pars = c(NAA_SS[1,1], 0.12),
                                           recruit_pars = mean(NAA_SS[,1])),
                             growth = list(model = 'Richards',
                                           re = c('none', 'none', 'none', 'none'),
                                           init_vals = GWpars[1:4],
                                           est_pars = c(1:4),
                                           SD_vals = GWpars[5:6],
                                           SD_est = 1:2),
                             LW = list(re = c('none', 'none'),
                                       init_vals = LWpars),
                             catchability = list(re = c('none'), 
                                                 initial_q = Q_pars, q_lower = 0,
                                                 q_upper = 1000, prior_sd = NA),
                             age_comp = 'dirichlet-pool0',
                             ecov = ecov,
                             basic_info = wham_data)


# update some inputs as SS model:
input_b = post_input_EBSpcod(input_b, SS_report, NAA_SS)
input_b$random = "Ecov_re"

#Run model:
fit_b = fit_wham(input_b, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton = 0)
check_convergence(fit_b)
save(fit_b, file = 'EBS_pcod/fit_b.RData')

# Make plots
dir.create(path = 'EBS_pcod/fit_b')
plot_wham_output(mod = fit_b, dir.main = 'EBS_pcod/fit_b', out.type = 'png')

# Get SSB estimates:
this_model = fit_b
model_name = 'ecov_WHAM'
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBdata2 = cbind(model = model_name, year=1977:2022,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
SSBdata2$ssb_min = exp(SSBdata2$est-1.96*SSBdata2$sd)
SSBdata2$ssb_max = exp(SSBdata2$est+1.96*SSBdata2$sd)
SSBdata2$est = exp(SSBdata2$est)

# -------------------------------------------------------------------------
# Plot SSB:

# Get SSB from SS:
SS_SSB = SS_report$derived_quants[grep(pattern = 'SSB_', x = SS_report$derived_quants$Label),]
SSBdata0 = cbind(data.frame(model = 'SS', year = wham_data$years), 
              SS_SSB[3:48, c('Value', 'StdDev')])
colnames(SSBdata0)[3:4] = c('est', 'sd')
SSBdata0$ssb_min = SSBdata0$est - 1.96*SSBdata0$sd
SSBdata0$ssb_max = SSBdata0$est + 1.96*SSBdata0$sd

# Merge data:
plot_data = rbind(SSBdata0, SSBdata1, SSBdata2)
plot_data$est = plot_data$est*1E-06
plot_data$ssb_min = plot_data$ssb_min*1E-06
plot_data$ssb_max = plot_data$ssb_max*1E-06
plot_data$model = factor(plot_data$model, levels = c('SS','AR(1)_WHAM', 'ecov_WHAM'))

# Make plot:
p1 = ggplot(plot_data, aes(year, est, ymin=ssb_min, ymax=ssb_max,
                      fill=model, color=model)) +
  ylim(0,NA) + labs(y='SSB (million mt)', x = NULL) +
  geom_ribbon(alpha=.3, color = NA) + geom_line(lwd=1) +
  labs( color=NULL, fill=NULL) +
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols) +
  theme_bw() +
  annotate("text", label = 'A', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(fill=guide_legend(title=NULL,nrow = 1), 
         color=guide_legend(title=NULL,nrow = 1),
         shape = guide_legend(override.aes = list(linewidth = 0.8))) +
  theme(legend.position=c(0.5,0.1), 
        legend.text = element_text(size = 10),
        legend.background = element_blank())
#ggsave(filename = 'EBS_pcod/compare_SSB.png', width = 190, height = 120, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------
# Plot L1 temporal variability:
LAA_a = data.frame(year = min_year:max_year, L1 = fit_a$rep$LAA[,1], model = 'AR(1)_WHAM')
LAA_b = data.frame(year = min_year:max_year, L1 = fit_b$rep$LAA[,1], model = 'ecov_WHAM')
LAA_SS = data.frame(year = min_year:(max_year-1), L1 = SS_report$growthseries[,6], model = 'SS')

# Merge datasets:
plot_data = rbind(LAA_a, LAA_b, LAA_SS)
plot_data$model = factor(plot_data$model, levels = c('SS','AR(1)_WHAM', 'ecov_WHAM'))

# Make plot:
p2 = ggplot(plot_data, aes(year, L1, color=model)) +
  geom_line(lwd=1) +
  labs(y='Mean length-at-age 1 (cm)', x = NULL) +
  labs( color=NULL, fill=NULL) +
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols) +
  theme_bw() +
  coord_cartesian(ylim=c(6, 18)) +
  annotate("text", label = 'B', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(fill=guide_legend(title=NULL,nrow = 1), 
         color=guide_legend(title=NULL,nrow = 1),
         shape = guide_legend(override.aes = list(linewidth = 0.8))) +
  theme(legend.position=c(0.5,0.1), 
        legend.text = element_text(size = 10),
        legend.background = element_blank())
# ggsave(filename = 'EBS_pcod/compare_L1.png', width = 190, height = 120, units = 'mm', dpi = 500)

# Plot Ecov:
p3 = plot_ecov_fit(fit_b, label = 'C')

# Merge plots:
png(filename = 'EBS_pcod/main_EBSpcod.png', width = 190, height = 160, units = 'mm', res = 500)
gridExtra::grid.arrange(p1, p2, p3, ncol = 2)
dev.off()


# -------------------------------------------------------------------------

