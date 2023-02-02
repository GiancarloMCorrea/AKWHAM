rm(list = ls())

# Set you WD:
setwd("~/GitHub/AKWHAM")

require(dplyr)
require(ggplot2)
theme_set(theme_bw())
library(wham)
library(snowfall)
source("simulation_functions.R")
source("aux_fun.R")

# Number of replicates:
nsim <- 100
ts_variables = c('log_F', 'log_SSB')
ts_labels = c('Fishing mortality', 'Spawning biomass')
#sc_variables = c('log_F1', 'log_N1_pars', 'logit_q', 'mean_rec_pars')
sc_variables = c('log_N1_pars', 'logit_q', 'mean_rec_pars', 'WAA_a', 'sigma_WAA')
sc_labels = c('Initial abundance', 'Catchability', 'Mean recruitment', 'Mean weight-at-age', 'Sigma_WAA')
save_folder = 'GOA_pollock'

### ------------------------------------------------------------
### First the self test of the original penalized likelihood model
## Rerun CM to get OM
load("GOA_pollock/fit_a.RData")
# Extract some parameters relevant for the simulation:
year0 = fit_a$years[1] - 1
nyears = length(fit_a$years)
nages = length(fit_a$ages.lab)

## vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
fit_a$input$data$simulate_data <- c(1,1,0)
fit_a$input$data$simulate_state <- fit_a$input$data$simulate_state*0 # why turn this off?
## temporarily turning estimation of selectivity off since it was
## causing some weird convergence failures
fit_a$input$map$logit_selpars <- factor(NA*fit_a$input$par$logit_selpars)
## fit_a$input$map$selpars_re <- factor(NA*fit_a$input$par$selpars_re)
om <- fit_wham(fit_a$input, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton=0)
om$par[which.max(abs(om$gr()))] # check the OM worked ok
om$sdrep
om$input$par <- om$env$parList(par=om$opt$par)
om2 <- fit_wham(om$input, do.osa = FALSE, do.fit = TRUE,
               do.retro = FALSE, n.newton=0)
om2$par[which.max(abs(om2$gr()))] # check the OM worked ok
om$sdrep

### Run and save the simulations
## test <- run_em(1)
true <- summary(om$sdrep, select='all')[,1]
sfInit(parallel=TRUE, cpus=10)
sfExport('run_em','true','om')
out <- sfLapply(1:nsim, run_em)
sfStop()
results <- bind_rows(out) %>% group_by(rep, par) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se,
         i=1:n(), re=(est-true)/true,
         year=year0+1:n()+ifelse(n()==nyears,0,NA),
         age=1:n()+ifelse(n()==nages,0,NA)) %>%
  ungroup
saveRDS(results, file.path(save_folder, 'sim_results_a2.RDS'))

### Process and plot them
results <- readRDS(file.path(save_folder, 'sim_results_a2.RDS'))
results %>% filter(parnum==maxpar) %>% pull(maxgrad) %>% summary
## (t)ime (s)eries we care about
ts <- filter(results, !is.na(year) & par %in% ts_variables)
## (sc)alar quantities we care about
sc <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables,"$", sep=""), collapse = '|'),par))

# Make plot SC:
sc$par = factor(sc$par, levels = sc_variables, labels = sc_labels)
p1 = ggplot(sc, aes(factor(i), re)) + 
        geom_violin(fill = '#bdd7e7') + 
        facet_wrap('par', scales='free_x')+
        xlab(NULL) + ylab('Relative error') +
        geom_hline(yintercept=0, color=2) 
        #coord_cartesian(ylim=.2*c(-1,1))
#ggsave(filename = file.path(save_folder, 'sim_sc_plot_a.png'), width = 190, height = 70, units = 'mm', dpi = 500)

# Make plot TS:
ts$par = factor(ts$par, levels = ts_variables, labels = ts_labels)
g = ggplot(ts, aes(year,y= re)) + 
      facet_wrap('par', nrow=1)+ 
      geom_hline(yintercept=0, color=2) +
      xlab(NULL) + ylab('Relative error') +
      coord_cartesian(ylim=.25*c(-1,1))
g1 = add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#2171b5')
#ggsave(filename = file.path(save_folder, 'sim_ts_plot_a.png'), width = 90, height = 130, units = 'mm', dpi = 500)

# Merge plots:
png(filename = 'GOA_pollock/sim_plot_a.png', width = 190, height = 160, units = 'mm', res = 500)
gridExtra::grid.arrange(p1, g1, ncol = 1)
dev.off()

### End of the fit_a simulation
### ------------------------------------------------------------

### ------------------------------------------------------------
### Test IID random effect smoothing of the WAA
## Rerun CM to get OM
load("GOA_pollock/fit_b.RData")
## vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
fit_b$input$data$simulate_data <- c(1,1,0)
fit_b$input$data$simulate_state <- fit_b$input$data$simulate_state*0
## temporarily turning estimation of selectivity off since it was
## causing some weird convergence failures
fit_b$input$map$logit_selpars <- factor(NA*fit_b$input$par$logit_selpars)
fit_b$input$map$selpars_re <- factor(NA*fit_b$input$par$selpars_re)
om <- fit_wham(fit_b$input, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton=0)
om$par[which.max(abs(om$gr(om$opt$par)))]         # check the OM worked ok
om$sdrep

## quick look at WAA estimates compared to data
plot_waa_fit(om, minyr=1990, maxyr=2009)
ggsave(file.path(save_folder, 'sim_WAA_fit_cohort.png'), width = 190, height = 190, units = 'mm', dpi = 500)
plot_waa_fit(om, by.cohort=FALSE,sizeX = 7.5)
ggsave(file.path(save_folder, 'sim_WAA_fit_year.png'), width = 190, height = 120, units = 'mm', dpi = 500)
plot_waa_resids(om)

### Run and save the simulations
true <- summary(om$sdrep, select='all')[,1]
## test <- run_em(1, TRUE)
nsim <- 100
sfInit(parallel=TRUE, cpus=10)
sfExport('run_em', 'true', 'om')
out <- sfLapply(1:nsim, run_em)
sfStop()
results <- bind_rows(out) %>% group_by(rep, par) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se,
         i=1:n(), re=(est-true)/true,
         year=year0+1:n()+ifelse(n()==nyears,0,NA),
         age=1:n()+ifelse(n()==nages,0,NA)) %>%
  ungroup
saveRDS(results, 'GOA_pollock/sim_results_b2.RDS')

### Process and plot them
results <- readRDS('GOA_pollock/sim_results_b2.RDS')
results %>% filter(parnum==maxpar) %>% pull(maxgrad) %>% summary
## (t)ime (s)eries we care about
ts <- filter(results, !is.na(year) & par %in% ts_variables)
## (sc)alar quantities we care about
sc <- filter(results, is.na(year) &  grepl(paste0(paste("^",sc_variables,"$", sep=""), collapse = '|'),par))
## need to carefully splice out the second WAA matrix since
## that's used for SSB.
ind <- as.numeric(array(1:length(fit_b$input$data$waa), dim=dim(fit_b$input$data$waa))[2,,])
waa <- filter(results, par=='pred_waa' & i %in% ind) %>% group_by(rep) %>%
  mutate(age=rep(1:nages, each=nyears), year=rep(fit_b$years, times=nages))

# Make plot SC:
sc$par = factor(sc$par, levels = sc_variables, labels = sc_labels)
p2 = ggplot(sc, aes(factor(i), re)) + 
  geom_violin(fill = '#bdd7e7') + 
  facet_wrap('par', scales='free_x')+
  xlab(NULL) + ylab('Relative error') +
  geom_hline(yintercept=0, color=2) 
#coord_cartesian(ylim=.2*c(-1,1))
ggsave(filename = file.path(save_folder, 'sim_sc_plot_b.png'), width = 190, height = 140, units = 'mm', dpi = 500)

# Make plot TS:
ts$par = factor(ts$par, levels = ts_variables, labels = ts_labels)
g = ggplot(ts, aes(year,y= re)) + 
  facet_wrap('par', nrow=2)+ 
  geom_hline(yintercept=0, color=2) +
  xlab(NULL) + ylab('Relative error') +
  coord_cartesian(ylim=.2*c(-1,1))
add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#2171b5')
ggsave(filename = file.path(save_folder, 'sim_ts_plot_b.png'), width = 90, height = 130, units = 'mm', dpi = 500)

# Plot WAA:
g <- ggplot(waa, aes(year,re)) + 
  facet_wrap('age')+ 
  geom_hline(yintercept=0, color=2)+
  ylab('Relative error') + xlab(NULL)
  #coord_cartesian(ylim=.5*c(-1,1)) + 
add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#2171b5')
ggsave(filename = file.path(save_folder, 'sim_waa_plot_b.png'), width = 190, height = 140, units = 'mm', dpi = 500)

### End of the fit_b simulation
### ------------------------------------------------------------

