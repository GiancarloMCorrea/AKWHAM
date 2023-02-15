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

# Load some files:

# SS files and outputs:
# SS model can be found here: 
# https://github.com/pete-hulson/goa_pcod/tree/master/2022/Stock_Synthesis_files/Model19.1a%20(22)%20-%20wADFG
# Change Month column from 1 to 7 for surveys (len comps and CAAL) 
data_file = r4ss::SS_readdat_3.30(file = 'SS_models/GOA_pcod/GOAPcod2022Oct25_wADFG.dat')
SS_report = r4ss::SS_output(dir = 'SS_models/GOA_pcod') # from OM

# Some model parameters:
mycols = c("#D16103", "#52854C")
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
selpars2 = SelecParams$Value[7:12]
selpars3 = SelecParams$Value[13:18]
selpars4 = SelecParams$Value[19:24]
selpars5 = SelecParams$Value[25:30]

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
                                                #re = c(rep('none', times = 5)),
                                                initial_pars=list(selpars1, selpars2, selpars3,
                                                                  selpars4, selpars5),
                                                fix_pars = list(5:6,5:6,4:6,NULL,1:6),
                                                #fix_pars = list(1:6,1:6,1:6,1:6,1:6),
                                                n_selblocks = 5),
                               M = list(model = 'constant', re = 'none',
                                        initial_means = SS_report$Natural_Mortality[1,5],
                                        est_ages = 1),
                               NAA_re = list(sigma="rec", cor = 'iid', N1_model = 1,
                                             recruit_model = 2,
                                             N1_pars = c(NAA_SS[1,1], 0),
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
                                                   q_upper = c(1000,1000), prior_sd = c(NA,NA)),
                               ecov = ecov,
                               basic_info = wham_data)

# update some inputs as SS model:
input_a = post_input_GOApcod(input_a, SS_report, NAA_SS)
# no random effects:
input_a$random <- NULL

#Run model:
fit_a = fit_wham(input_a, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton = 0)
check_convergence(fit_a)
save(fit_a, file = 'GOA_pcod/fit_a.RData')

# Make plots
dir.create(path = 'GOA_pcod/fit_a')
plot_wham_output(mod = fit_a, dir.main = 'GOA_pcod/fit_a', out.type = 'pdf')

# -------------------------------------------------------------------------
# Plot SSB:

# Get SSB from SS:
SS_SSB = SS_report$derived_quants[grep(pattern = 'SSB_', x = SS_report$derived_quants$Label),]
data1 = cbind(data.frame(model = 'SS', year = fit_a$years), 
              SS_SSB[3:48, c('Value', 'StdDev')])
colnames(data1)[3:4] = c('est', 'sd')
data1$ssb_min = data1$est - 1.96*data1$sd
data1$ssb_max = data1$est + 1.96*data1$sd

# WHAM model:
this_model = fit_a
model_name = 'WHAM'
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
data2 = cbind(model = model_name, year=1977:2022,
              filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
data2$ssb_min = exp(data2$est-1.96*data2$sd)
data2$ssb_max = exp(data2$est+1.96*data2$sd)
data2$est = exp(data2$est)

# Make plot:
plot_data =rbind(data1, data2)
p1 = ggplot(plot_data, aes(year, est*1e-06, ymin=ssb_min*1e-06, ymax=ssb_max*1e-06,
                      fill=model, color=model)) +
  ylim(0,NA) + labs(y='SSB (million tons)', x = NULL) +
  geom_ribbon(alpha=.3, color = NA) + geom_line(lwd=1) +
  labs( color=NULL, fill=NULL) +
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols) +
  theme_bw() +
  annotate("text", label = 'A', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(fill=guide_legend(title=NULL,nrow = 1), 
         color=guide_legend(title=NULL,nrow = 1),
         shape = guide_legend(override.aes = list(linewidth = 0.8))) +
  theme(legend.position=c(0.5,0.95), 
        legend.text = element_text(size = 10),
        legend.background = element_blank())
# ggsave(filename = 'GOA_pcod/compare_SSB.png', width = 190, height = 120, units = 'mm', dpi = 500)

# -------------------------------------------------------------------------
# Plot LAA:

# Get LAA from SS:
data1 = data.frame(ages = SS_report$endgrowth$int_Age[2:11], 
                   LAA = SS_report$endgrowth$Len_Beg[2:11], 
                   SD = SS_report$endgrowth$SD_Beg[2:11], 
                   model = 'SS')
data1$len_min = data1$LAA - 1.96*data1$SD
data1$len_max = data1$LAA + 1.96*data1$SD

# Get LAA from WHAM:
data2 = data.frame(ages = 1:n_ages, LAA = fit_a$rep$LAA[1,], 
                   SD = fit_a$rep$SDAA[1,], model = 'WHAM')
data2$len_min = data2$LAA - 1.96*data2$SD
data2$len_max = data2$LAA + 1.96*data2$SD

# Make plot:
plot_data =rbind(data1, data2)
p2 = ggplot(plot_data, aes(ages, LAA, ymin=len_min, ymax=len_max,
                      fill=model, color=model)) +
  ylim(0,NA) + labs(y='Mean length (cm)', x = NULL) +
  geom_ribbon(alpha=.3, color = NA) + geom_line(lwd=1) +
  labs( color=NULL, fill=NULL) +
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols) +
  theme_bw() +
  annotate("text", label = 'B', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  scale_x_continuous(breaks = 1:10, labels = 1:10) +
  guides(fill=guide_legend(title=NULL,nrow = 1), 
         color=guide_legend(title=NULL,nrow = 1),
         shape = guide_legend(override.aes = list(linewidth = 0.8))) +
  theme(legend.position=c(0.5,0.95), 
        legend.text = element_text(size = 10),
        legend.background = element_blank())

# Merge plots:
png(filename = 'GOA_pcod/main_GOApcod.png', width = 190, height = 70, units = 'mm', res = 500)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()
