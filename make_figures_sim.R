rm(list = ls())
source("simulation_functions.R")
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

# GOA pollock simulation plots --------------------------------------------

# Variables to plot:
nages = 10
years = 1970:2021
nyears = length(years)
save_folder = 'GOA_pollock'
model_names = c('wham_ewaa', 'wham_iid', 'wham_2dar1')

# Model a:
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

# Model b:
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
  coord_cartesian(ylim= c(-0.2,0.2)) 
add_ci(d, ci=c(.5,.95), alpha=c(.4,.4), fill = '#606060', showMedian = TRUE)
ggsave(filename = file.path(save_folder, 'sim_waa_plot_b.png'), width = 190, height = 140, units = 'mm', dpi = 500)


## Model c:
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
  coord_cartesian(ylim= c(-0.2,0.2)) 
add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#606060', showMedian = TRUE)
ggsave(filename = file.path(save_folder, 'sim_waa_plot_c.png'), width = 190, height = 140, units = 'mm', dpi = 500)

# Save all sim plots:
lay_mat = matrix(c(1,2,2,3,4,4,5,6,6), ncol = 3)
png(filename = 'GOA_pollock/sim_plot_merged.png', width = 190, height = 180, units = 'mm', res = 500)
gridExtra::grid.arrange(p1, g1, p2, g2, p3, g3, layout_matrix = lay_mat)
dev.off()

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
results <- readRDS(file.path(save_folder, 'sim_results_a2.RDS'))
results$par[results$par == 'growth_a'] = rep(c('K', 'Linf', 'L1'), times = nsim) # put right name to growth parameters
results$par[results$par == 'SD_par'] = rep(c('SD1', 'SDA'), times = nsim) # put right name to growth parameters
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
  coord_cartesian(ylim=c(-0.3,0.3)) 

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
  coord_cartesian(ylim= c(-0.4,0.4)) 
g1 = add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#606060', showMedian = TRUE) # '#2171b5'

# Save all sim plots:
lay_mat = matrix(c(1,2,3,2), ncol = 2)
png(filename = 'GOA_pcod/sim_plot_merged.png', width = 190, height = 180, units = 'mm', res = 500)
gridExtra::grid.arrange(p1, g1, p2, layout_matrix = lay_mat)
dev.off()

# BS pcod simulation plots --------------------------------------------

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

# For Ecov model:
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

# For AR1 model:
nsim2 = 80
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

# Save all sim plots:
lay_mat = matrix(c(1,2,3,3,4,5,6,6), ncol = 2)
png(filename = 'EBS_pcod/sim_plot_merged.png', width = 190, height = 220, units = 'mm', res = 500)
gridExtra::grid.arrange(p1, p2, g1, p3, p4, g2, layout_matrix = lay_mat)
dev.off()