# Install WHAM package (growth branch):
# remotes::install_github(repo = 'GiancarloMCorrea/wham', ref='growth', INSTALL_opts = c("--no-docs", "--no-multiarch", "--no-demo"))

# Set you WD:
setwd("~/GitHub/AKWHAM")

rm(list = ls())
require(dplyr)
require(ggplot2)
library(wham)
source('functions_pollock.R')
source('aux_fun.R')

# Load some files:
load('aux_data/input.RData') # input data as in original WHAM model (Cole's)
load('aux_data/cv_survey.RData') # CV for observed WAA (survey)
load('aux_data/asdrep.RData') # CV for observed WAA (survey)
# see here: https://github.com/afsc-assessments/GOApollock/blob/dev/pkwham/pkwham_test.R

# Some model parameters:
LAA_jan1 = c(8.8, 21.1, 30.4, 37.4, 42.7, 46.8, 49.8, 52.1, 53.9, 55.2) # obtained from survey data
LW_pars = c(2.847e-06, 3.261) # obtained from survey data
WAA_jan1 = round(LW_pars[1]*LAA_jan1^LW_pars[2], digits = 3)

# -------------------------------------------------------------------------
# Prepare input data for WHAM:
wham_data = list()
wham_data$ages = 1:input$data$n_ages
wham_data$lengths = seq(from = 2, to = 80, by = 2)
wham_data$years = as.integer(input$years)

wham_data$n_fleets = input$data$n_fleets
wham_data$agg_catch = input$data$agg_catch
wham_data$use_agg_catch = input$data$use_agg_catch
wham_data$catch_cv = input$data$agg_catch_sigma
# Age comps fleet:
wham_data$catch_paa = input$data$catch_paa
wham_data$catch_Neff = input$data$catch_Neff
wham_data$use_catch_paa = input$data$use_catch_paa
wham_data$selblock_pointer_fleets = input$data$selblock_pointer_fleets
wham_data$F = matrix(0.2, ncol = 1, nrow = input$data$n_years_catch)

wham_data$n_indices = input$data$n_indices
wham_data$agg_indices = input$data$agg_indices
# Age comps survey
wham_data$index_paa = input$data$index_paa
wham_data$index_Neff =  input$data$index_Neff
wham_data$use_index_paa = input$data$use_index_paa

wham_data$units_indices = input$data$units_indices
wham_data$use_indices = input$data$use_indices
wham_data$units_index_paa = input$data$units_index_paa

wham_data$selblock_pointer_indices = input$data$selblock_pointer_indices
wham_data$fracyr_indices = input$data$fracyr_indices
wham_data$waa = input$data$waa
wham_data$waa_pointer_indices = input$data$waa_pointer_indices
wham_data$waa_pointer_fleets = input$data$waa_pointer_fleets
wham_data$waa_pointer_totcatch = input$data$waa_pointer_totcatch
wham_data$waa_pointer_ssb = input$data$waa_pointer_ssb
wham_data$waa_pointer_jan1 = input$data$waa_pointer_jan1
wham_data$maturity = input$data$mature
wham_data$fracyr_SSB = input$data$fracyr_SSB
wham_data$Fbar_ages = input$data$Fbar_ages
wham_data$percentSPR = input$data$percentSPR
wham_data$percentFXSPR = 100
wham_data$percentFMSY = 100
wham_data$XSPR_R_avg_yrs = 1:length(input$years)
wham_data$XSPR_R_opt = input$data$XSPR_R_opt
wham_data$simulate_period = input$data$simulate_period
wham_data$bias_correct_process = 1
wham_data$bias_correct_observation = 1

# Add CV and turnon likelihood for waa:
cvsurvey = as.matrix(out_data)
cvsurvey[is.na(cvsurvey)] = 0
wham_data$waa_cv = array(0, dim = dim(wham_data$waa))
wham_data$waa_cv[2,,] = cvsurvey
wham_data$use_index_waa = matrix(0L, ncol = wham_data$n_indices, nrow = length(wham_data$years))
wham_data$use_index_waa[17:52,1] = 1L

# -------------------------------------------------------------------------
# 2DAR1 WAA
# use same input data created in iid WAA
input_c = prepare_wham_input(model_name = "pollock_c",
                             selectivity = list(model = c('double-logistic', 'age-specific', 'double-logistic',
                                                          'double-logistic', 'age-specific', 'age-specific',
                                                          'double-logistic'),
                                                re = c('iid', 'none', 'none', 'none', 'none', 'none', 'none'),
                                                initial_pars=list(c(4, 1, 20, 0.36), rep(1, 10),
                                                                  c(4, 1, 20, 0.36),
                                                                  c(4, 1, 20, 0.36), rep(1, 10), rep(1, 10),
                                                                  c(4, 1, 20, 0.36)),
                                                fix_pars = list(NULL, 1:2, 3:4, 3:4, 1:10, 1:10, 1:4),
                                                n_selblocks = 7),
                             M = list(model = 'age-specific', re = 'none',
                                      initial_means = exp(input$par$M_a)),
                             NAA_re = list(sigma="rec", cor = 'iid', N1_model = 1,
                                           recruit_model = 2,
                                           N1_pars = exp(input$par$log_N1_pars),
                                           recruit_pars = exp(input$par$mean_rec_pars)),
                             WAA = list(WAA_vals = WAA_jan1,
                                        re = c('2dar1'),
                                        est_pars = 1:input$data$n_ages),
                             catchability = list(re = c('ar1', 'none', 'ar1', 'none', 'none', 'none'),
                                                 initial_q = rep(1, 6), q_lower = rep(0, 6),
                                                 q_upper = rep(1000, 6), prior_sd = rep(NA, 6)),
                             basic_info = wham_data)

# update some inputs:
input_c = post_input_pollock(input_c, input)
# random WAA only
input_c$random <- c('WAA_re')
#Fix survey selex for age 3,4, and 8 to reach convergence:
tmp = matrix(input_c$map$logit_selpars, nrow = 7)
tmp[2,c(3:4,8)] = NA
input_c$map$logit_selpars = factor(tmp)

#Run model:
fit_c = fit_wham(input_c, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE)

# -------------------------------------------------------------------------
# 3DAR1 WAA:
# wham_data$Var3D_ParamW = 0 # conditional
wham_data$Var3D_ParamW = 1 # marginal
input_d = prepare_wham_input(model_name = "pollock_d",
                             selectivity = list(model = c('double-logistic', 'age-specific', 'double-logistic',
                                                          'double-logistic', 'age-specific', 'age-specific',
                                                          'double-logistic'),
                                                re = c('iid', 'none', 'none', 'none', 'none', 'none', 'none'),
                                                initial_pars=list(c(4, 1, 20, 0.36), rep(1, 10),
                                                                  c(4, 1, 20, 0.36),
                                                                  c(4, 1, 20, 0.36), rep(1, 10), rep(1, 10),
                                                                  c(4, 1, 20, 0.36)),
                                                fix_pars = list(NULL, 1:2, 3:4, 3:4, 1:10, 1:10, 1:4),
                                                n_selblocks = 7),
                             M = list(model = 'age-specific', re = 'none',
                                      initial_means = exp(input$par$M_a)),
                             NAA_re = list(sigma="rec", cor = 'iid', N1_model = 1,
                                           recruit_model = 2,
                                           N1_pars = exp(input$par$log_N1_pars),
                                           recruit_pars = exp(input$par$mean_rec_pars)),
                             WAA = list(WAA_vals = WAA_jan1,
                                        re = c('3dgmrf'),
                                        est_pars = 1:input$data$n_ages),
                             catchability = list(re = c('ar1', 'none', 'ar1', 'none', 'none', 'none'),
                                                 initial_q = rep(1, 6), q_lower = rep(0, 6),
                                                 q_upper = rep(1000, 6), prior_sd = rep(NA, 6)),
                             basic_info = wham_data)

# update some inputs:
input_d = post_input_pollock(input_d, input)
# random WAA only
input_d$random <- c('WAA_re')
# Fix survey selex for age 3,4, and 8 to reach convergence:
tmp = matrix(input_d$map$logit_selpars, nrow = 7)
tmp[2,c(3:4,8)] = NA
input_d$map$logit_selpars = factor(tmp)
# input_d$map$WAA_repars = factor(c(1,2,3,NA)) # fix rho_c

#Run model:
fit_d = fit_wham(input_d, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE)

# -------------------------------------------------------------------------
require(ggplot2)
#save(fit_c, file = '3d_compare/fit_c.RData')
#save(fit_d, file = '3d_compare/fit_d.RData')

# MAIN PERIOD:
df_c = fit_c$rep$pred_waa[2,,]
colnames(df_c) = 1:10
rownames(df_c) = 1970:2021
df_c = reshape2::melt(df_c)
df_c$type = '2dar1'

df_d = fit_d$rep$pred_waa[2,,]
colnames(df_d) = 1:10
rownames(df_d) = 1970:2021
df_d = reshape2::melt(df_d)
df_d$type = '3dgmrf'

plot_data = rbind(df_c, df_d)

ggplot(data = plot_data, aes(x = Var1, y = value, color = factor(type))) +
  geom_line() +
  xlab('Age') + ylab('mean weight (kg)') +
  theme_bw() +
  theme(legend.position = 'top') +
  facet_wrap(~factor(Var2), scales = 'free_y')
ggsave(filename = '3d_compare/compare_2d_3d.png', width = 190, height = 160, units = 'mm', dpi = 500)

# Compare SSB:
this_model = fit_c
model_name = '2dar1'
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBdata2 = cbind(model = model_name, year=1970:2021,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
# Model c:
#load('3d_compare/fit_c.RData')
this_model = fit_d
model_name = '3dgmrf'
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBdata3 = cbind(model = model_name, year=1970:2021,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))

# Compare SSB among models:
plot_data = rbind(SSBdata2, SSBdata3)
plot_data$ssb_min = exp(plot_data$est-1.96*plot_data$sd)
plot_data$ssb_max = exp(plot_data$est+1.96*plot_data$sd)
plot_data$est = exp(plot_data$est)

ggplot(plot_data, aes(year, est*1e-06, ymin=ssb_min*1e-06, ymax=ssb_max*1e-06,
                           fill=model, color=model)) +
  labs(y='SSB (million tons)', x = NULL) +
  geom_ribbon(alpha=.3, color = NA) + geom_line(lwd=1) +
  labs( color=NULL, fill=NULL) +
  coord_cartesian(xlim = c(1970,2021), ylim = c(0,1.2)) +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  annotate("text", label = 'A', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(fill=guide_legend(title=NULL,nrow = 2),
         color=guide_legend(title=NULL,nrow = 2),
         shape = guide_legend(override.aes = list(linewidth = 0.8))) +
  theme(legend.position=c(0.5,0.9),
        legend.text = element_text(size = 10),
        legend.background = element_blank())
ggsave(filename = '3d_compare/compare_biomass_2d_3d.png', width = 190, height = 160, units = 'mm', dpi = 500)

# Compare NLL components (nll_WAA is random effects and nll_waa is observed mean w-at-a)
fit_c$rep[grep('nll',names(fit_c$rep))] %>% lapply(sum) %>% unlist
fit_d$rep[grep('nll',names(fit_d$rep))] %>% lapply(sum) %>% unlist

# Plots and AIC (uses marginal likelihood):
compare_wham_models(mods = list(d2ar1 = fit_c, d3gmrf = fit_d), table.opts = list(calc.rho = FALSE))
