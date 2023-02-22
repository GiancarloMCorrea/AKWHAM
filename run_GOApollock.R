# Install WHAM package (growth branch):
# remotes::install_github(repo = 'gmoroncorrea/wham', ref='growth', INSTALL_opts = c("--no-docs", "--no-multiarch", "--no-demo"))

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
# see here: https://github.com/afsc-assessments/GOApollock/blob/dev/pkwham/pkwham_test.R
load('aux_data/cv_survey.RData') # CV for observed WAA (survey) 
load('aux_data/asdrep.RData') # ADMB pollock model output
# see here: https://github.com/afsc-assessments/GOApollock/blob/dev/pkwham/pkwham_test.R

model_names = c('ADMB', 'EWAA_WHAM', 'iid_WHAM', '2DAR1_WHAM')

# Get SSB estimates from ADMB model
SSBdata0 = cbind(model = model_names[1],
            filter(asdrep, name=='log_ssb') %>% select(year, est, sd=se))

# Some model parameters:
mycols = c("#D16103", "#52854C", "#4E84C4")
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
wham_data$bias_correct_process = input$data$bias_correct_pe
wham_data$bias_correct_observation = input$data$bias_correct_oe


# -------------------------------------------------------------------------
# Empirical WAA
input_a = prepare_wham_input(model_name="pollock_a",
                               selectivity=list(model = c('double-logistic', 'age-specific', 'double-logistic', 
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
                               catchability = list(re = c('ar1', 'none', 'ar1', 'none', 'none', 'none'), 
                                                   initial_q = rep(1, 6), q_lower = rep(0, 6),
                                                   q_upper = rep(1000, 6), prior_sd = rep(NA, 6)),
                               basic_info = wham_data)

# update some inputs as base WHAM model:
input_a = post_input_pollock(input_a, input)
# no random effects:
input_a$random <- NULL
tmp = matrix(input_a$map$logit_selpars, nrow = 7)
tmp[2,c(3:4,8)] = NA
input_a$map$logit_selpars = factor(tmp)

#Run model:
fit_a = fit_wham(input_a, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton = 0)
check_convergence(fit_a)
save(fit_a, file = 'GOA_pollock/fit_a.RData')

# Make plots
dir.create(path = 'GOA_pollock/fit_a')
plot_wham_output(mod = fit_a, dir.main = 'GOA_pollock/fit_a', out.type = 'pdf')

# Get SSB estimates:
this_model = fit_a
model_name = model_names[2]
tmp = data.frame(name = names(this_model$sdrep$value),
                   est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBdata1 = cbind(model = model_name, year=1970:2021,
              filter(tmp, name=='log_SSB') %>% dplyr::select(-name))

# -------------------------------------------------------------------------
# iid WAA

# Add CV and turnon likelihood for waa:
cvsurvey = as.matrix(out_data)
cvsurvey[is.na(cvsurvey)] = 0
wham_data$waa_cv = array(0, dim = dim(wham_data$waa))
wham_data$waa_cv[2,,] = cvsurvey
wham_data$use_index_waa = matrix(0L, ncol = wham_data$n_indices, nrow = length(wham_data$years))
wham_data$use_index_waa[17:52,1] = 1L

input_b = prepare_wham_input(model_name = "pollock_b",
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
                                       re = c('iid'),
                                       est_pars = 1:input$data$n_ages),
                            catchability = list(re = c('ar1', 'none', 'ar1', 'none', 'none', 'none'), 
                                                initial_q = rep(1, 6), q_lower = rep(0, 6),
                                                q_upper = rep(1000, 6), prior_sd = rep(NA, 6)),
                            basic_info = wham_data)

# update some inputs:
input_b = post_input_pollock(input_b, input)
# random WAA only
input_b$random <- c('WAA_re')

#Run model:
fit_b = fit_wham(input_b, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton = 0)
check_convergence(fit_b)
save(fit_b, file = 'GOA_pollock/fit_b.RData')

# Make plots
dir.create(path = 'GOA_pollock/fit_b')
plot_wham_output(mod = fit_b, dir.main = 'GOA_pollock/fit_b', out.type = 'pdf')

# Get SSB estimates:
this_model = fit_b
model_name = model_names[3]
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBdata2 = cbind(model = model_name, year=1970:2021,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))

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

#Run model:
fit_c = fit_wham(input_c, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton = 0)
check_convergence(fit_c)
save(fit_c, file = 'GOA_pollock/fit_c.RData')

# Make plots
dir.create(path = 'GOA_pollock/fit_c')
plot_wham_output(mod = fit_c, dir.main = 'GOA_pollock/fit_c', out.type = 'pdf')

# Get SSB estimates:
this_model = fit_c
model_name = model_names[4]
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBdata3 = cbind(model = model_name, year=1970:2021,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))


# -------------------------------------------------------------------------
# Compare SSB among models

plot_data = rbind(SSBdata0, SSBdata1, SSBdata2, SSBdata3)
plot_data$model = factor(plot_data$model, levels = model_names)

# Plot SSB (mean):
p1 = ggplot(plot_data, aes(x = year, y = exp(est)*1e-06, color=model)) +
  ylim(0,0.8) + labs(y='SSB (million tons)', x= NULL) +
  geom_line(linewidth=1) +
  labs( color=NULL, fill=NULL) +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  annotate("text", label = 'A', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(color=guide_legend(title=NULL,nrow = 2),
         shape = guide_legend(override.aes = list(size = 0.1))) +
  theme(legend.position=c(0.5,0.9), 
        legend.text = element_text(size = 7.5),
        legend.background = element_blank())

p2 = ggplot(plot_data, aes(x = year, y = sd/est, color=model)) +
  ylim(0,0.035) + labs(y='CV SSB', x= NULL) +
  geom_line(linewidth=1) +
  labs( color=NULL, fill=NULL) +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  annotate("text", label = 'B', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(color=guide_legend(title=NULL,nrow = 2),
         shape = guide_legend(override.aes = list(size = 0.1))) +
  theme(legend.position=c(0.5,0.9), 
        legend.text = element_text(size = 7.5),
        legend.background = element_blank())

# Plot LAA and WAA over time ----------------------------------------------

# iid WHAM:
WAA_data = fit_b$rep$WAA_re
colnames(WAA_data) = input$ages.lab
rownames(WAA_data) = input$years
WAA_plot = reshape2::melt(data = WAA_data, varnames = c('year', 'age'))

p3 = ggplot(WAA_plot, aes(x = year, y = value, color=age)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  labs(y='WAA random effects (iid)', x= NULL) +
  labs( color=NULL) +
  scale_color_brewer(palette = 'Spectral') +
  ylim(-0.75, 1) +
  theme_bw() +
  annotate("text", label = 'C', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(fill=guide_legend(title=NULL,nrow = 2), 
         color=guide_legend(title=NULL,nrow = 2),
         shape = guide_legend(override.aes = list(size = 0.1))) +
  theme(legend.position=c(0.5,0.9), 
        legend.text = element_text(size = 7.5),
        legend.background = element_blank()) 


# 2dAR1 WHAM:
WAA_data = fit_c$rep$WAA_re
colnames(WAA_data) = input$ages.lab
rownames(WAA_data) = input$years
WAA_plot = reshape2::melt(data = WAA_data, varnames = c('year', 'age'))

p4 = ggplot(WAA_plot, aes(x = year, y = value, color=age)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  labs(y='WAA random effects (2DAR1)', x= NULL) +
  labs( color=NULL) +
  scale_color_brewer(palette = 'Spectral') +
  ylim(-0.75, 1) +
  theme_bw() +
  annotate("text", label = 'D', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(fill=guide_legend(title=NULL,nrow = 2), 
         color=guide_legend(title=NULL,nrow = 2),
         shape = guide_legend(override.aes = list(size = 0.1))) +
  theme(legend.position=c(0.5,0.9), 
        legend.text = element_text(size = 7.5),
        legend.background = element_blank())
#ggsave(filename = 'GOA_pollock/WAA_re.png', width = 190, height = 120, units = 'mm', dpi = 500)

# Merge plots:
png(filename = 'GOA_pollock/main_GOApollock.png', width = 190, height = 160, units = 'mm', res = 500)
gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()


# -------------------------------------------------------------------------
# Plot fits:

plot_waa_fit(fit = fit_b, minyr=1990, maxyr=2009, by.cohort = FALSE)
ggsave(filename = 'GOA_pollock/summary_WAA_fit_b.png', width = 190, height = 140, units = 'mm', dpi = 500)

plot_waa_fit(fit = fit_c, minyr=1990, maxyr=2009, by.cohort = FALSE)
ggsave(filename = 'GOA_pollock/summary_WAA_fit_c.png', width = 190, height = 140, units = 'mm', dpi = 500)
