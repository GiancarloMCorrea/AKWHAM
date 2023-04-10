rm(list = ls())

# Set you WD:
setwd("C:/Users/moroncog/Documents/GitHub/AKWHAM")

require(dplyr)
require(ggplot2)
theme_set(theme_bw())
library(wham)
library(snowfall)
source("simulation_functions.R")
source("aux_fun.R")

# Number of replicates:
nsim <- 100
ts_variables = c('log_Fbar', 'log_SSB')
ts_labels = c('Fishing mortality', 'Spawning biomass')
#sc_variables = c('log_F1', 'log_N1_pars', 'logit_q', 'mean_rec_pars')
sc_variables = c('log_N1_pars', 'logit_q', 'mean_rec_pars', 'growth_a', 'M_a')
sc_labels = c('Initial abundance', 'Catchability', 'Mean recruitment', 'Growth parameters', 'Natural mortality')
save_folder = 'EBS_pcod'

### ------------------------------------------------------------
### Model a: AR(1) on L1
## Rerun CM to get OM
load("EBS_pcod/fit_a.RData")
# Extract some parameters relevant for the simulation:
year0 = fit_a$years[1] - 1
nyears = length(fit_a$years)
nages = length(fit_a$ages.lab)

# create OM data:
OMinput = fit_a$input
## vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
OMinput$data$simulate_data <- c(1,1,1) 
OMinput$data$simulate_state <- OMinput$data$simulate_state*0 
# Fix selectivity parameters? Cole, decide this:
OMinput$map$logit_selpars <- factor(NA*OMinput$map$logit_selpars)
OMinput$map$selpars_re <- factor(NA*OMinput$map$selpars_re)
# Run OM again: Make sure you have ~2 hours
om <- fit_wham(OMinput, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton=0)
om$par[which.max(abs(om$gr()))] # check the OM worked ok
om$sdrep

#save(om, file = 'om.RData')
# om2 <- fit_wham(om$input, do.osa = FALSE, do.fit = TRUE,
#                do.retro = FALSE, n.newton=0)
# om2$par[which.max(abs(om2$gr()))] # check the OM worked ok
# om2$sdrep

### Run and save the simulations
## test <- run_em(1)
true <- summary(om$sdrep, select='all')[,1]
sfInit(parallel=TRUE, cpus=5)
sfExport('run_em','true','om')
out <- sfLapply(1:nsim, run_em)
sfStop()
results <- bind_rows(out) %>% group_by(rep, par) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se,
         i=1:n(), re=(est-true)/true,
         year=year0+1:n()+ifelse(n()==nyears,0,NA),
         age=1:n()+ifelse(n()==nages,0,NA)) %>%
  ungroup
saveRDS(results, file.path(save_folder, 'sim_results_a.RDS'))

### Process and plot them
results <- readRDS(file.path(save_folder, 'sim_results_a.RDS'))
results %>% filter(parnum==maxpar) %>% pull(maxgrad) %>% summary
## (t)ime (s)eries we care about
ts <- filter(results, !is.na(year) & par %in% ts_variables)
## (sc)alar quantities we care about
sc <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables,"$", sep=""), collapse = '|'),par))

# Make plot SC:
sc$par = factor(sc$par, levels = sc_variables, labels = sc_labels)
ggplot(sc, aes(factor(i), re)) + 
        geom_violin(fill = '#bdd7e7') + 
        facet_wrap('par', scales='free_x')+
        xlab(NULL) + ylab('Relative error') +
        geom_hline(yintercept=0, color=2) 
        #coord_cartesian(ylim=.2*c(-1,1))
ggsave(filename = file.path(save_folder, 'sim_sc_plot_a.png'), width = 190, height = 160, units = 'mm', dpi = 500)

# Make plot TS:
ts$par = factor(ts$par, levels = ts_variables, labels = ts_labels)
g = ggplot(ts, aes(year,y= re)) + 
      facet_wrap('par', nrow=2)+ 
      geom_hline(yintercept=0, color=2) +
      xlab(NULL) + ylab('Relative error')
      #coord_cartesian(ylim=.15*c(-1,1))
add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#2171b5')
ggsave(filename = file.path(save_folder, 'sim_ts_plot_a.png'), width = 90, height = 130, units = 'mm', dpi = 500)

### ------------------------------------------------------------
### Model b: Ecov on L1
## Rerun CM to get OM
load("EBS_pcod/fit_b.RData")
# Extract some parameters relevant for the simulation:
year0 = fit_b$years[1] - 1
nyears = length(fit_b$years)
nages = length(fit_b$ages.lab)

# create OM data:
OMinput = fit_b$input
## vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
OMinput$data$simulate_data <- c(1,1,1) 
OMinput$data$simulate_state <- OMinput$data$simulate_state*0 
# Fix selectivity parameters? Cole, decide this
OMinput$map$logit_selpars <- factor(NA*OMinput$map$logit_selpars)
OMinput$map$selpars_re <- factor(NA*OMinput$map$selpars_re)
# Run OM again: Make sure you have ~2 hours
om <- fit_wham(OMinput, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton=0)
om$par[which.max(abs(om$gr()))] # check the OM worked ok
om$sdrep

#save(om, file = 'om.RData')
# om2 <- fit_wham(om$input, do.osa = FALSE, do.fit = TRUE,
#                do.retro = FALSE, n.newton=0)
# om2$par[which.max(abs(om2$gr()))] # check the OM worked ok
# om2$sdrep

### Run and save the simulations
## test <- run_em(1)
true <- summary(om$sdrep, select='all')[,1]
sfInit(parallel=TRUE, cpus=5)
sfExport('run_em','true','om')
out <- sfLapply(1:nsim, run_em)
sfStop()
results <- bind_rows(out) %>% group_by(rep, par) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se,
         i=1:n(), re=(est-true)/true,
         year=year0+1:n()+ifelse(n()==nyears,0,NA),
         age=1:n()+ifelse(n()==nages,0,NA)) %>%
  ungroup
saveRDS(results, file.path(save_folder, 'sim_results_b.RDS'))

### Process and plot them
results <- readRDS(file.path(save_folder, 'sim_results_b.RDS'))
results %>% filter(parnum==maxpar) %>% pull(maxgrad) %>% summary
## (t)ime (s)eries we care about
ts <- filter(results, !is.na(year) & par %in% ts_variables)
## (sc)alar quantities we care about
sc <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables,"$", sep=""), collapse = '|'),par))

# Make plot SC:
sc$par = factor(sc$par, levels = sc_variables, labels = sc_labels)
ggplot(sc, aes(factor(i), re)) + 
        geom_violin(fill = '#bdd7e7') + 
        facet_wrap('par', scales='free_x')+
        xlab(NULL) + ylab('Relative error') +
        geom_hline(yintercept=0, color=2) 
        #coord_cartesian(ylim=.2*c(-1,1))
ggsave(filename = file.path(save_folder, 'sim_sc_plot_b.png'), width = 190, height = 160, units = 'mm', dpi = 500)

# Make plot TS:
ts$par = factor(ts$par, levels = ts_variables, labels = ts_labels)
g = ggplot(ts, aes(year,y= re)) + 
      facet_wrap('par', nrow=2)+ 
      geom_hline(yintercept=0, color=2) +
      xlab(NULL) + ylab('Relative error')
      #coord_cartesian(ylim=.15*c(-1,1))
add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#2171b5')
ggsave(filename = file.path(save_folder, 'sim_ts_plot_b.png'), width = 90, height = 130, units = 'mm', dpi = 500)

