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
ts_variables = c('F', 'SSB', 'recruits', 'pred_catch')
ts_labels = c('Fishing mortality', 'Spawning biomass', 'Recruits',
              'Catch')
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
OMinput$par <- fit_a$parList #fit_a$env$parList(fit_a$opt$par)
## vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
OMinput$data$simulate_data <- c(1,1,1)
OMinput$data$simulate_state <- OMinput$data$simulate_state*0
# Fix selectivity parameters to reach convergence (especially for model a)
OMinput$map$logit_selpars <- factor(NA*OMinput$par$logit_selpars)
OMinput$map$selpars_re <- factor(NA*OMinput$par$selpars_re)
## OMinput$map$growth_repars <- factor(NA*OMinput$par$growth_repars)
### no Ecov effect so turn off to speed things up
OMinput$map$Ecov_process_pars <- factor(NA*OMinput$par$Ecov_process_pars)
OMinput$map$Ecov_re <- factor(NA*OMinput$par$Ecov_re)
OMinput$map$M_a <- factor(NA*OMinput$par$M_a)
OMinput$par$SDgrowth_par <- OMinput$par$SD_par # temporary due to version issues
OMinput$map$SDgrowth_par <- OMinput$map$SD_par # temporary due to version issues
OMinput$map$SD_par<- NULL
## OMinput$map$log_NAA <- factor(NA*OMinput$par$log_NAA)
## OMinput$map$SDgrowth_par <- factor(NA*OMinput$par$SDgrowth_par)
## OMinput$map$growth_a <- factor(NA*OMinput$par$growth_a)
## OMinput$map$growth_a <- factor(c(1,NA,3,4))
## OMinput$map$growth_re <- factor(NA*OMinput$par$growth_re)
# OMinput$map$growth_repars <- factor(NA*OMinput$par$growth_repars)
## OMinput$map$logit_q <- factor(NA*OMinput$par$logit_q)
## OMinput$map$log_F1 <- factor(NA*OMinput$par$log_F1)
## OMinput$map$log_N1_pars <- factor(NA*OMinput$par$log_N1_pars)
## OMinput$map$mean_rec_pars <- factor(NA*OMinput$par$mean_rec_pars)
OMinput$random <- 'growth_re'
## Switch to multinomial for simulation
OMinput$data$index_Neff <- 1*OMinput$data$index_Neff
OMinput$data$index_NeffL <- 1*OMinput$data$index_NeffL
OMinput$data$age_comp_model_fleets <- 1
OMinput$data$age_comp_model_indices <- 1
## multinomial has no pars to estimate
OMinput$map$index_paa_pars <- factor(NA*OMinput$par$index_paa_pars)
## turn off ageing error matrices
OMinput$data$use_index_aging_error <- 0
OMinput$data$use_catch_aging_error <- 0

## Build but don't refit it
om <- fit_wham(OMinput, do.osa = FALSE, do.fit = FALSE, do.retro = FALSE,
               n.newton=0, do.sdrep=FALSE)
names(om$par) %>% unique

### Run and save the simulations
sfInit(parallel=TRUE, cpus=10)
nsim <- 10
outdir <- 'EBS_pcod/simfits_a'
sfExport('run_em2','om', 'outdir')
out <- sfLapply(1:nsim, run_em2)
sfStop()
results <- bind_rows(out) %>% group_by(rep, par) %>%
  mutate(i=1:n(), re=(est-true)/true,
         year=year0+1:n()+ifelse(n()==nyears,0,NA),
         age=1:n()+ifelse(n()==nages,0,NA)) %>%
  ungroup
saveRDS(results, file.path(save_folder, 'sim_results_a.RDS'))

### Process and plot them
results <- readRDS(file.path(save_folder, 'sim_results_a.RDS'))
results %>% group_by(rep) %>% slice_head(n=1) %>% pull(maxgrad)
## (t)ime (s)eries we care about
ts <- filter(results, !is.na(year) & par %in% ts_variables)
## (sc)alar quantities we care about
sc <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables,"$", sep=""), collapse = '|'),par))
## Make plot SC:
sc$par = factor(sc$par, levels = sc_variables, labels = sc_labels)
ggplot(sc, aes(factor(i), re)) +
        geom_violin(fill = '#bdd7e7') +
        facet_wrap('par', scales='free_x')+
        xlab(NULL) + ylab('Relative error') +
        geom_hline(yintercept=0, color=2)
        #coord_cartesian(ylim=.2*c(-1,1))
ggsave(filename = file.path(save_folder, 'sim_sc_plot_a.png'), width = 190, height = 160, units = 'mm', dpi = 500)
## Make plot TS:
ts$par = factor(ts$par, levels = ts_variables, labels = ts_labels)
g = ggplot(ts, aes(year,y= re)) +
      facet_wrap('par', nrow=2)+
      geom_hline(yintercept=0, color=2) +
      xlab(NULL) + ylab('Relative error')
      #coord_cartesian(ylim=.15*c(-1,1))
add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#2171b5')
ggsave(filename = file.path(save_folder, 'sim_ts_plot_a.png'), width = 90, height = 130, units = 'mm', dpi = 500)



### ------------------------------------------------------------
### Model b: ecov on L1
load("EBS_pcod/fit_b.RData")
# Extract some parameters relevant for the simulation:
year0 = fit_b$years[1] - 1
nyears = length(fit_b$years)
nages = length(fit_b$ages.lab)


# create OM data:
OMinput = fit_b$input
OMinput$par <- fit_b$parList #fit_b$env$parList(fit_b$opt$par)
## vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
OMinput$data$simulate_data <- c(1,1,1)
OMinput$data$simulate_state <- OMinput$data$simulate_state*0
# Fix selectivity parameters to reach convergence (especially for model a)
OMinput$map$logit_selpars <- factor(NA*OMinput$par$logit_selpars)
OMinput$map$selpars_re <- factor(NA*OMinput$par$selpars_re)
#OMinput$map$growth_repars <- factor(NA*OMinput$par$growth_repars)
## OMinput$map$Ecov_process_pars <- factor(NA*OMinput$par$Ecov_process_pars)
## OMinput$map$Ecov_re <- factor(NA*OMinput$par$Ecov_re)
OMinput$map$M_a <- factor(NA*OMinput$par$M_a)
OMinput$par$SDgrowth_par <- OMinput$par$SD_par # temporary due to version issues
OMinput$map$SDgrowth_par <- OMinput$map$SD_par # temporary due to version issues
OMinput$map$SD_par<- NULL
## OMinput$map$log_NAA <- factor(NA*OMinput$par$log_NAA)
## OMinput$map$SDgrowth_par <- factor(NA*OMinput$par$SDgrowth_par)
## OMinput$map$growth_a <- factor(NA*OMinput$par$growth_a)
## OMinput$map$growth_a <- factor(c(1,NA,3,4))
## OMinput$map$growth_re <- factor(NA*OMinput$par$growth_re)
OMinput$map$growth_repars <- factor(NA*OMinput$par$growth_repars)
## OMinput$map$logit_q <- factor(NA*OMinput$par$logit_q)
## OMinput$map$log_F1 <- factor(NA*OMinput$par$log_F1)
## OMinput$map$log_N1_pars <- factor(NA*OMinput$par$log_N1_pars)
## OMinput$map$mean_rec_pars <- factor(NA*OMinput$par$mean_rec_pars)
OMinput$random <- 'Ecov_re'
## Switch to multinomial for simulation
OMinput$data$index_Neff <- 1*OMinput$data$index_Neff
OMinput$data$index_NeffL <- 1*OMinput$data$index_NeffL
OMinput$data$age_comp_model_fleets <- 1
OMinput$data$age_comp_model_indices <- 1
## multinomial has no pars to estimate
OMinput$map$index_paa_pars <- factor(NA*OMinput$par$index_paa_pars)
## turn off ageing error matrices
OMinput$data$use_index_aging_error <- 0
OMinput$data$use_catch_aging_error <- 0

## Build OM but don't run it
om <- fit_wham(OMinput, do.osa = FALSE, do.fit = FALSE, do.retro = FALSE,
               n.newton=0, do.sdrep=FALSE)
names(om$par) %>% unique
### Run and save the simulations
sfInit(parallel=TRUE, cpus=5)
nsim <- 5
outdir <- 'EBS_pcod/simfits_b'
sfExport('run_em2','om', 'outdir')
out <- sfLapply(1:nsim, run_em2)
sfStop()
results <- bind_rows(out) %>% group_by(rep, par) %>%
  mutate(i=1:n(), re=(est-true)/true,
         year=year0+1:n()+ifelse(n()==nyears,0,NA),
         age=1:n()+ifelse(n()==nages,0,NA)) %>%
  ungroup
saveRDS(results, file.path(save_folder, 'sim_results_a.RDS'))

### Process and plot them
results <- readRDS(file.path(save_folder, 'sim_results_a.RDS'))
results %>% group_by(rep) %>% slice_head(n=1) %>% pull(maxgrad)
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
