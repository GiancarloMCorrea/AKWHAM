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
sc_variables = c('log_N1_pars', 'logit_q', 'mean_rec_pars', 'growth_a', 'M_a', 'M_re')
sc_labels = c('Initial abundance', 'Catchability', 'Mean recruitment', 'Growth parameters', 'Natural mortality', 'Natural mortality (2014-2016)')
save_folder = 'GOA_pcod'

### ------------------------------------------------------------
### First the self test of the original penalized likelihood model
## Rerun CM to get OM
load("GOA_pcod/fit_a.RData")
# Extract some parameters relevant for the simulation:
year0 = fit_a$years[1] - 1
nyears = length(fit_a$years)
nages = length(fit_a$ages.lab)

# create OM data:
OMinput = fit_a$input
## vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
OMinput$data$simulate_data <- c(1,1,1)
OMinput$data$simulate_state <- OMinput$data$simulate_state*0 # why turn this off?
# Integer Neff for CAAL:
# For some reason, in SS this is a number less than 1 in most cases
# This step is particularly important for the simulation part
OMinput$data$index_caal_Neff = OMinput$data$index_caal_Neff*(1/0.14)
OMinput$data$catch_caal_Neff = OMinput$data$catch_caal_Neff*(1/0.14)
OMinput$data$obs$val[OMinput$data$obs$type == 'catchcaal' | OMinput$data$obs$type == 'indexcaal'] = OMinput$data$obs$val[OMinput$data$obs$type == 'catchcaal' | OMinput$data$obs$type == 'indexcaal']*(1/0.14)
#Fix selectivity:
OMinput$map$logit_selpars = factor(rep(NA, times = length(OMinput$map$logit_selpars)))
OMinput$map$selpars_re = factor(rep(NA, times = length(OMinput$map$selpars_re)))
# Run OM again:
om <- fit_wham(OMinput, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton=0)
om$par[which.max(abs(om$gr()))] # check the OM worked ok
om$sdrep
om$input$par <- om$env$parList(par=om$opt$par)
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
ggplot(sc, aes(factor(i), re)) + 
        geom_violin(fill = '#bdd7e7') + 
        facet_wrap('par', scales='free_x')+
        xlab(NULL) + ylab('Relative error') +
        geom_hline(yintercept=0, color=2) 
        #coord_cartesian(ylim=.2*c(-1,1))
ggsave(filename = file.path(save_folder, 'sim_sc_plot.png'), width = 190, height = 160, units = 'mm', dpi = 500)

# Make plot TS:
ts$par = factor(ts$par, levels = ts_variables, labels = ts_labels)
g = ggplot(ts, aes(year,y= re)) + 
      facet_wrap('par', nrow=2)+ 
      geom_hline(yintercept=0, color=2) +
      xlab(NULL) + ylab('Relative error')
      #coord_cartesian(ylim=.15*c(-1,1))
add_ci(g, ci=c(.5,.95), alpha=c(.4,.4), fill = '#2171b5')
ggsave(filename = file.path(save_folder, 'sim_ts_plot.png'), width = 90, height = 130, units = 'mm', dpi = 500)

