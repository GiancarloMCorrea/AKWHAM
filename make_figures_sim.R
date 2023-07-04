rm(list = ls())
source("simulation_functions.R")
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

# -------------------------------------------------------------------------
# GOA pollock simulation plots --------------------------------------------

# Variables to plot:
nages = 10
years = 1970:2021
nyears = length(years)
save_folder = 'GOA_pollock'
model_names = c('wham_ewaa', 'wham_iid', 'wham_2dar1')

# Model a: ------------------------------------------
ts_variables = c('log_F', 'log_SSB')
ts_labels = c('Annual F', 'SSB')
sc_variables = c('log_N1_pars', 'mean_rec_pars') # 'WAA_a', 'sigma_WAA'

results <- readRDS(file.path(save_folder, 'sim_results_a.RDS'))
#results %>% filter(parnum==maxpar) %>% pull(maxgrad) %>% summary
ts <- filter(results, !is.na(year) & par %in% ts_variables)
ts = ts %>% dplyr::mutate(re2 = (exp(est)-exp(true))/exp(true))
sc <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables,"$", sep=""), collapse = '|'),par))
sc = sc %>% dplyr::mutate(re2 = (exp(est)-exp(true))/exp(true))

# Make plot SC:
sc$par = factor(sc$par, levels = sc_variables)
p1 = ggplot(sc, aes(x = par, y = re2)) + 
  geom_boxplot(fill = '#606060', alpha = 0.5) + 
  xlab(NULL) + ylab('Relative error') +
  scale_x_discrete(labels = c('log_N1_pars' = expression(N[1*","*1]),
                              'mean_rec_pars'   = expression(bar(R)))) +
  geom_hline(yintercept=0, color=2) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  ggtitle(label = model_names[1])

# Make plot TS:
ts$par = factor(ts$par, levels = ts_variables, labels = ts_labels)
g = ggplot(ts, aes(year,y= re2)) + 
  facet_wrap('par', nrow=2)+ 
  geom_hline(yintercept=0, color=2) +
  xlab(NULL) + ylab('Relative error') +
  coord_cartesian(ylim=c(-0.6,0.6)) +
  theme(strip.background = element_blank())
g1 = add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#606060', showMedian = TRUE) # '#2171b5'

# Stability plot:
stab_df = results[results$par == 'log_SSB', ]
this_rep = sort(unique(stab_df$rep))[1:100]
save_df = list()
for(i in seq_along(this_rep)) {
  tmp = stab_df[stab_df$rep %in% this_rep[1:i], ]
  med_val = median((exp(tmp$est) - exp(tmp$true))/exp(tmp$true))
  sd_val = sd((exp(tmp$est) - exp(tmp$true))/exp(tmp$true))
  save_df[[i]] = data.frame(rep = i, value = c(med_val, sd_val), type = c('Bias', 'Precision'),
                            model = model_names[1])
}
stab_goa_1 = dplyr::bind_rows(save_df)
  

# Model b: -----------------------------------------------------
ts_variables = c('log_F', 'log_SSB')
ts_labels = c('Annual F', 'SSB')
sc_variables = c('log_N1_pars', 'mean_rec_pars', 'sigma_WAA') # 'WAA_a', 'sigma_WAA'

results <- readRDS('GOA_pollock/sim_results_b.RDS')
results %>% filter(parnum==maxpar) %>% pull(maxgrad) %>% summary
ts <- filter(results, !is.na(year) & par %in% ts_variables)
ts = ts %>% dplyr::mutate(re2 = (exp(est)-exp(true))/exp(true))
sc <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables,"$", sep=""), collapse = '|'),par))
sc$re2 = sc$re
id_log = sc$par %in% sc_variables[1:2]
sc$re2[id_log] = (exp(sc$est[id_log])-exp(sc$true[id_log]))/exp(sc$true[id_log])
## need to carefully splice out the second WAA matrix since
## that's used for SSB.
ind <- as.numeric(array(1:(nyears*nages*4), dim=c(4,nyears,nages))[2,,])
waa <- filter(results, par=='pred_waa' & i %in% ind) %>% group_by(rep) %>%
  mutate(age=rep(1:nages, each=nyears), year=rep(years, times=nages))

# Make plot SC:
sc$par = factor(sc$par, levels = sc_variables)
p2 = ggplot(sc, aes(x = par, y = re2)) + 
  geom_boxplot(fill = '#606060', alpha = 0.5) + 
  xlab(NULL) + ylab('Relative error') +
  scale_x_discrete(labels = c('log_N1_pars' = expression(N[1*","*1]),
                              'mean_rec_pars'   = expression(bar(R)),
                              'sigma_WAA' = expression(sigma[WAA]))) +
  geom_hline(yintercept=0, color=2) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  ggtitle(label = model_names[2])

# Make plot TS:
ts$par = factor(ts$par, levels = ts_variables, labels = ts_labels)
g = ggplot(ts, aes(year,y= re2)) + 
  facet_wrap('par', nrow=2)+ 
  geom_hline(yintercept=0, color=2) +
  xlab(NULL) + ylab('Relative error') +
  coord_cartesian(ylim=c(-0.6,0.6)) +
  theme(strip.background = element_blank())
g2 = add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#606060', showMedian = TRUE)

# Plot WAA:
d <- ggplot(waa, aes(year,re)) + 
  facet_wrap('age')+ 
  geom_hline(yintercept=0, color=2)+
  ylab('Relative error') + xlab(NULL) +
  coord_cartesian(ylim= c(-0.2,0.2)) +
  theme(axis.text.x =element_text(size = 8))
add_ci(d, ci=c(.5,.95), alpha=c(.4,.4), fill = '#606060', showMedian = TRUE)
ggsave(filename = file.path(save_folder, 'sim_waa_plot_b.jpg'), width = 170, height = 130, units = 'mm', dpi = 400)

# Stability plot:
stab_df = results[results$par == 'log_SSB', ]
this_rep = sort(unique(stab_df$rep))[1:100]
save_df = list()
for(i in seq_along(this_rep)) {
  tmp = stab_df[stab_df$rep %in% this_rep[1:i], ]
  med_val = median((exp(tmp$est) - exp(tmp$true))/exp(tmp$true))
  sd_val = sd((exp(tmp$est) - exp(tmp$true))/exp(tmp$true))
  save_df[[i]] = data.frame(rep = i, value = c(med_val, sd_val), type = c('Bias', 'Precision'),
                            model = model_names[2])
}
stab_goa_2 = dplyr::bind_rows(save_df)

# Model c: ------------------------------------------------
ts_variables = c('log_F', 'log_SSB')
ts_labels = c('Annual F', 'SSB')
sc_variables = c('log_N1_pars', 'mean_rec_pars', 'sigma_WAA', 'rho_WAA_a', 'rho_WAA_y') # 'WAA_a', 'sigma_WAA'

results <- readRDS('GOA_pollock/sim_results_c.RDS')
results %>% filter(parnum==maxpar) %>% pull(maxgrad) %>% summary
ts <- filter(results, !is.na(year) & par %in% ts_variables)
ts = ts %>% dplyr::mutate(re2 = (exp(est)-exp(true))/exp(true))
sc <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables,"$", sep=""), collapse = '|'),par))
sc$re2 = sc$re
id_log = sc$par %in% sc_variables[1:2]
sc$re2[id_log] = (exp(sc$est[id_log])-exp(sc$true[id_log]))/exp(sc$true[id_log])
## need to carefully splice out the second WAA matrix since
## that's used for SSB.
ind <- as.numeric(array(1:(nyears*nages*4), dim=c(4,nyears,nages))[2,,])
waa <- filter(results, par=='pred_waa' & i %in% ind) %>% group_by(rep) %>%
  mutate(age=rep(1:nages, each=nyears), year=rep(years, times=nages))

# Make plot SC:
sc$par = factor(sc$par, levels = sc_variables)
p3 = ggplot(sc, aes(x = par, y = re2)) + 
  geom_boxplot(fill = '#606060', alpha = 0.5) + 
  xlab(NULL) + ylab('Relative error') +
  scale_x_discrete(labels = c('log_N1_pars' = expression(N[1*","*1]),
                              'mean_rec_pars'   = expression(bar(R)),
                              'sigma_WAA' = expression(sigma[WAA]),
                              'rho_WAA_a' = expression(rho[age]),
                              'rho_WAA_y' = expression(rho[year]))) +
  geom_hline(yintercept=0, color=2) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  ggtitle(label = model_names[3])

# Make plot TS:
ts$par = factor(ts$par, levels = ts_variables, labels = ts_labels)
g = ggplot(ts, aes(year,y= re2)) + 
  facet_wrap('par', nrow=2)+ 
  geom_hline(yintercept=0, color=2) +
  xlab(NULL) + ylab('Relative error') +
  coord_cartesian(ylim=.6*c(-1,1)) +
  theme(strip.background = element_blank())
g3 = add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#606060', showMedian = TRUE)

# Plot WAA:
g <- ggplot(waa, aes(year,re)) + 
  facet_wrap('age')+ 
  geom_hline(yintercept=0, color=2)+
  ylab('Relative error') + xlab(NULL)+
  coord_cartesian(ylim= c(-0.2,0.2)) +
  theme(axis.text.x =element_text(size = 8))
add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#606060', showMedian = TRUE)
ggsave(filename = file.path(save_folder, 'sim_waa_plot_c.jpg'), width = 170, height = 130, units = 'mm', dpi = 400)

# Stability plot:
stab_df = results[results$par == 'log_SSB', ]
this_rep = sort(unique(stab_df$rep))[1:100]
save_df = list()
for(i in seq_along(this_rep)) {
  tmp = stab_df[stab_df$rep %in% this_rep[1:i], ]
  med_val = median((exp(tmp$est) - exp(tmp$true))/exp(tmp$true))
  sd_val = sd((exp(tmp$est) - exp(tmp$true))/exp(tmp$true))
  save_df[[i]] = data.frame(rep = i, value = c(med_val, sd_val), type = c('Bias', 'Precision'),
                            model = model_names[3])
}
stab_goa_3 = dplyr::bind_rows(save_df)

# Save all sim plots: ---------------------------------------------
lay_mat = matrix(c(1,2,2,3,4,4,5,6,6), ncol = 3)
jpeg(filename = 'GOA_pollock/sim_plot_merged.jpg', width = 170, height = 160, units = 'mm', res = 400)
gridExtra::grid.arrange(p1, g1, p2, g2, p3, g3, layout_matrix = lay_mat)
dev.off()

# Save stability plot: --------------------------------------------
plot_df = rbind(stab_goa_1, stab_goa_2, stab_goa_3)
plot_df$model = factor(plot_df$model, levels = model_names)

ggplot(plot_df, aes(x = rep, y = value)) +
  geom_line() +
  xlab('Number of replicates') +
  ylab(NULL) +
  facet_grid(factor(type) ~ model, scales = 'free_y')
ggsave(filename = 'GOA_pollock/sim_stability.jpg', width = 170, height = 120, units = 'mm', dpi = 400)

# -------------------------------------------------------------------------
# GOA pcod simulation plots --------------------------------------------

# Variables to plot:
nages = 10
nyears = 52
nsim = 100
ts_variables = c('log_Fbar', 'log_SSB')
ts_labels = c('Annual F', 'SSB')
sc_variables1 = c('log_N1_pars', 'mean_rec_pars') 
sc_variables2 = c('K', 'Linf', 'L1', 'SD1', 'SDA')
save_folder = 'GOA_pcod'
model_name = 'wham'

# Model a:
results <- readRDS(file.path(save_folder, 'sim_results_a.RDS'))
results$par[results$par == 'growth_a'] = rep(c('K', 'Linf', 'L1'), times = nsim) # put right name to growth parameters
results$par[results$par == 'SDgrowth_par'] = rep(c('SD1', 'SDA'), times = nsim) # put right name to growth parameters
#results %>% filter(parnum==maxpar) %>% pull(maxgrad) %>% summary
ts <- filter(results, par %in% ts_variables)
ts = ts %>% dplyr::mutate(re2 = (exp(est)-exp(true))/exp(true))
ts$year = rep(1977:2022, times = 2*nsim)
sc1 <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables1,"$", sep=""), collapse = '|'),par))
sc1 = sc1 %>% dplyr::mutate(re2 = (exp(est)-exp(true))/exp(true))
sc2 <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables2,"$", sep=""), collapse = '|'),par))
sc2 = sc2 %>% dplyr::mutate(re2 = (exp(est)-exp(true))/exp(true))

# Make plot SC1:
sc1$par = factor(sc1$par, levels = sc_variables1)
p1 = ggplot(sc1, aes(x = par, y = re2)) + 
  geom_boxplot(fill = '#606060', alpha = 0.5) + 
  xlab(NULL) + ylab('Relative error') +
  scale_x_discrete(labels = c('log_N1_pars' = expression(N[1*","*1]),
                              'mean_rec_pars'   = expression(bar(R)))) +
  geom_hline(yintercept=0, color=2) +
  coord_cartesian(ylim=c(-0.5,0.5)) 

# Make plot SC2:
sc2$par = factor(sc2$par, levels = sc_variables2)
p2 = ggplot(sc2, aes(x = par, y = re2)) + 
  geom_boxplot(fill = '#606060', alpha = 0.5) + 
  xlab(NULL) + ylab('Relative error') +
  scale_x_discrete(labels = c('K' = expression(k),
                              'Linf'   = expression(L[inf]),
                              'L1' = expression(L[1]),
                              'SD1' = expression(SD[1]),
                              'SDA' = expression(SD[A]))) +
  geom_hline(yintercept=0, color=2) +
  coord_cartesian(ylim=c(-0.1,0.1)) 

# Make plot TS:
ts$par = factor(ts$par, levels = ts_variables, labels = ts_labels)
g = ggplot(ts, aes(year,y= re2)) + 
  facet_wrap('par', nrow=1, scales = 'free_y')+ 
  geom_hline(yintercept=0, color=2) +
  xlab(NULL) + ylab('Relative error') +
  theme(strip.background = element_blank()) +
  coord_cartesian(ylim= c(-0.7,0.7)) 
g1 = add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#606060', showMedian = TRUE) # '#2171b5'

# Stability plot:
stab_df = results[results$par == 'log_SSB', ]
this_rep = sort(unique(stab_df$rep))[1:100]
save_df = list()
for(i in seq_along(this_rep)) {
  tmp = stab_df[stab_df$rep %in% this_rep[1:i], ]
  med_val = median((exp(tmp$est) - exp(tmp$true))/exp(tmp$true))
  sd_val = sd((exp(tmp$est) - exp(tmp$true))/exp(tmp$true))
  save_df[[i]] = data.frame(rep = i, value = c(med_val, sd_val), type = c('Bias', 'Precision'),
                            model = model_name)
}
stab_goapcod = dplyr::bind_rows(save_df)

# Save all sim plots: ------------------------------
lay_mat = matrix(c(1,2,3,2), ncol = 2)
jpeg(filename = 'GOA_pcod/sim_plot_merged.jpg', width = 170, height = 160, units = 'mm', res = 400)
gridExtra::grid.arrange(p1, g1, p2, layout_matrix = lay_mat)
dev.off()

# Save stability plot: --------------------------------------------
plot_df = stab_goapcod
plot_df$model = factor(plot_df$model, levels = model_name)

ggplot(plot_df, aes(x = rep, y = value)) +
  geom_line() +
  xlab('Number of replicates') +
  ylab(NULL) +
  facet_grid(model ~ factor(type), scales = 'free_y')
ggsave(filename = 'GOA_pcod/sim_stability.jpg', width = 170, height = 85, units = 'mm', dpi = 400)


# -------------------------------------------------------------------------
# EBS pcod simulation plots --------------------------------------------

# Variables to plot:
nages = 20
years = 1977:2022
nyears = length(years)
save_folder = 'EBS_pcod'
model_names = c('wham_ecov', 'wham_ar1')

# Model a:
ts_variables = c('F', 'SSB')
ts_labels = c('Annual F', 'SSB')
sc_variables1 = c('log_N1_pars', 'mean_rec_pars') 
sc_variables2 = c('K', 'Linf', 'L1', 'gamma', 'SD1', 'SDA')

# For Ecov model: ----------------------------------------
nsim1 = 115
results <- readRDS(file.path(save_folder, 'sim_results_b.RDS')) 
results$par[results$par == 'growth_a'] = rep(c('K', 'Linf', 'L1', 'gamma'), times = nsim1) # put right name to growth parameters
results$par[results$par == 'SD_par'] = rep(c('SD1', 'SDA'), times = nsim1) # put right name to growth parameters
#results %>% filter(parnum==maxpar) %>% pull(maxgrad) %>% summary
ts <- filter(results, par %in% ts_variables)
ts = ts %>% dplyr::mutate(re2 = (est-true)/true)
ts$year = rep(1977:2022, times = 2*nsim1)
sc1 <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables1,"$", sep=""), collapse = '|'),par))
sc1 = sc1 %>% dplyr::mutate(re2 = (exp(est)-exp(true))/exp(true))
sc1 = sc1 %>% mutate(par2 = paste0(par, i))
sc2 <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables2,"$", sep=""), collapse = '|'),par))
sc2 = sc2 %>% dplyr::mutate(re2 = (exp(est)-exp(true))/exp(true))

# Make plot SC1:
sc1$par = factor(sc1$par, levels = sc_variables1)
p1 = ggplot(sc1, aes(x = par2, y = re2)) + 
  geom_boxplot(fill = '#606060', alpha = 0.5) + 
  xlab(NULL) + ylab('Relative error') +
  scale_x_discrete(labels = c('log_N1_pars1' = expression(N[1*","*1]),
                              'log_N1_pars2' = expression(F[1]), 
                              'mean_rec_pars1'   = expression(bar(R)))) +
  geom_hline(yintercept=0, color=2) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  ggtitle(label = model_names[1])

# Make plot SC2:
sc2$par = factor(sc2$par, levels = sc_variables2)
p2 = ggplot(sc2, aes(x = par, y = re2)) + 
  geom_boxplot(fill = '#606060', alpha = 0.5) + 
  xlab(NULL) + ylab('Relative error') +
  scale_x_discrete(labels = c('K' = expression(k),
                              'Linf'   = expression(L[inf]),
                              'L1' = expression(L[1]),
                              'gamma' = expression(gamma),
                              'SD1' = expression(SD[1]),
                              'SDA' = expression(SD[A]))) +
  geom_hline(yintercept=0, color=2) +
  coord_cartesian(ylim=c(-0.1,0.1)) 

# Make plot TS:
ts$par = factor(ts$par, levels = ts_variables, labels = ts_labels)
g = ggplot(ts, aes(year,y= re2)) + 
  facet_wrap('par', nrow=2, scales = 'free_y')+ 
  geom_hline(yintercept=0, color=2) +
  xlab(NULL) + ylab('Relative error') +
  theme(strip.background = element_blank()) +
  coord_cartesian(ylim= c(-0.5,0.5)) 
g1 = add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#606060', showMedian = TRUE) # '#2171b5'

# Stability plot:
stab_df = results[results$par == 'SSB', ]
this_rep = sort(unique(stab_df$rep))[1:100]
save_df = list()
for(i in seq_along(this_rep)) {
  tmp = stab_df[stab_df$rep %in% this_rep[1:i], ]
  med_val = median((tmp$est - tmp$true)/tmp$true)
  sd_val = sd((tmp$est - tmp$true)/tmp$true)
  save_df[[i]] = data.frame(rep = i, value = c(med_val, sd_val), type = c('Bias', 'Precision'),
                            model = model_names[1])
}
stab_ebs_1 = dplyr::bind_rows(save_df)

# For AR1 model: ----------------------------------------------
nsim2 = 138
results <- readRDS(file.path(save_folder, 'sim_results_a.RDS')) 
results$par[results$par == 'growth_a'] = rep(c('K', 'Linf', 'L1', 'gamma'), times = nsim2) # put right name to growth parameters
results$par[results$par == 'SD_par'] = rep(c('SD1', 'SDA'), times = nsim2) # put right name to growth parameters
#results %>% filter(parnum==maxpar) %>% pull(maxgrad) %>% summary
ts <- filter(results, par %in% ts_variables)
ts = ts %>% dplyr::mutate(re2 = (est-true)/true)
ts$year = rep(1977:2022, times = 2*nsim2)
sc1 <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables1,"$", sep=""), collapse = '|'),par))
sc1 = sc1 %>% dplyr::mutate(re2 = (exp(est)-exp(true))/exp(true))
sc1 = sc1 %>% mutate(par2 = paste0(par, i))
sc2 <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables2,"$", sep=""), collapse = '|'),par))
sc2 = sc2 %>% dplyr::mutate(re2 = (exp(est)-exp(true))/exp(true))

# Make plot SC1:
sc1$par = factor(sc1$par, levels = sc_variables1)
p3 = ggplot(sc1, aes(x = par2, y = re2)) + 
  geom_boxplot(fill = '#606060', alpha = 0.5) + 
  xlab(NULL) + ylab('Relative error') +
  scale_x_discrete(labels = c('log_N1_pars1' = expression(N[1*","*1]),
                              'log_N1_pars2' = expression(F[1]), 
                              'mean_rec_pars1'   = expression(bar(R)))) +
  geom_hline(yintercept=0, color=2) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  ggtitle(label = model_names[2])

# Make plot SC2:
sc2$par = factor(sc2$par, levels = sc_variables2)
p4 = ggplot(sc2, aes(x = par, y = re2)) + 
  geom_boxplot(fill = '#606060', alpha = 0.5) + 
  xlab(NULL) + ylab('Relative error') +
  scale_x_discrete(labels = c('K' = expression(k),
                              'Linf'   = expression(L[inf]),
                              'L1' = expression(L[1]),
                              'gamma' = expression(gamma),
                              'SD1' = expression(SD[1]),
                              'SDA' = expression(SD[A]))) +
  geom_hline(yintercept=0, color=2) +
  coord_cartesian(ylim=c(-0.1,0.1)) 

# Make plot TS:
ts$par = factor(ts$par, levels = ts_variables, labels = ts_labels)
g = ggplot(ts, aes(year,y= re2)) + 
  facet_wrap('par', nrow=2, scales = 'free_y')+ 
  geom_hline(yintercept=0, color=2) +
  xlab(NULL) + ylab('Relative error') +
  theme(strip.background = element_blank()) +
  coord_cartesian(ylim= c(-0.5,0.5)) 
g2 = add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#606060', showMedian = TRUE) # '#2171b5'

# Stability plot:
stab_df = results[results$par == 'SSB', ]
this_rep = sort(unique(stab_df$rep))[1:100]
save_df = list()
for(i in seq_along(this_rep)) {
  tmp = stab_df[stab_df$rep %in% this_rep[1:i], ]
  med_val = median((tmp$est - tmp$true)/tmp$true)
  sd_val = sd((tmp$est - tmp$true)/tmp$true)
  save_df[[i]] = data.frame(rep = i, value = c(med_val, sd_val), type = c('Bias', 'Precision'),
                            model = model_names[2])
}
stab_ebs_2 = dplyr::bind_rows(save_df)

# Save all sim plots: ----------------------------------------
lay_mat = matrix(c(1,2,3,3,4,5,6,6), ncol = 2)
jpeg(filename = 'EBS_pcod/sim_plot_merged.jpg', width = 170, height = 210, units = 'mm', res = 400)
gridExtra::grid.arrange(p1, p2, g1, p3, p4, g2, layout_matrix = lay_mat)
dev.off()

# Save stability plot: --------------------------------------------
plot_df = rbind(stab_ebs_1, stab_ebs_2)
plot_df$model = factor(plot_df$model, levels = model_names)

ggplot(plot_df, aes(x = rep, y = value)) +
  geom_line() +
  xlab('Number of replicates') +
  ylab(NULL) +
  facet_grid(factor(type) ~ model, scales = 'free_y')
ggsave(filename = 'EBS_pcod/sim_stability.jpg', width = 170, height = 120, units = 'mm', dpi = 400)
