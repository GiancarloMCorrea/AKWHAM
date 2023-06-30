rm(list = ls())
require(dplyr)

# Set you WD:
setwd("C:/Users/moroncog/Documents/GitHub/AKWHAM")

require(dplyr)
require(ggplot2)
library(wham)
library(snowfall)
source("simulation_functions.R")
source("aux_fun.R")
nsim = 100
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
OMinput$data$simulate_data <- c(1,1,0)
OMinput$data$simulate_state <- OMinput$data$simulate_state*0 # why turn this off?
# Integer Neff for CAAL:
# For some reason, in SS this is a number less than 1 in most cases
# This step is particularly important for the simulation part
OMinput$data$index_caal_Neff = OMinput$data$index_caal_Neff*(1/0.14)
OMinput$data$catch_caal_Neff = OMinput$data$catch_caal_Neff*(1/0.14)
OMinput$data$obs$val[OMinput$data$obs$type == 'catchcaal' | OMinput$data$obs$type == 'indexcaal'] = OMinput$data$obs$val[OMinput$data$obs$type == 'catchcaal' | OMinput$data$obs$type == 'indexcaal']*(1/0.14)
#Fix selectivity:
#OMinput$map$logit_selpars = factor(rep(NA, times = length(OMinput$map$logit_selpars)))
#OMinput$map$selpars_re = factor(rep(NA, times = length(OMinput$map$selpars_re)))
# Fix F1:
#OMinput$map$log_F1 = factor(rep(NA, times = length(OMinput$map$log_F1)))
# Fix Q:
#OMinput$map$logit_q = factor(rep(NA, times = length(OMinput$map$logit_q)))
# Fix M: (problematic parameter):
OMinput$map$M_a = factor(NA)
OMinput$map$M_re = factor(rep(NA, times = length(OMinput$map$M_re)))
# Initial values:
OMinput$par = fit_a$parList
# Run OM again:
om <- fit_wham(OMinput, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton=0)
om$par[which.max(abs(om$gr()))] # check the OM worked ok
om$sdrep

#save(om, file = 'GOA_pcod/om.RData')
# om2 <- fit_wham(om$input, do.osa = FALSE, do.fit = TRUE,
#                do.retro = FALSE, n.newton=0)
# om2$par[which.max(abs(om2$gr()))] # check the OM worked ok
# om2$sdrep

### Run and save the simulations
#test <- run_em(1)
true <- summary(om$sdrep, select='all')[,1]
sfInit(parallel=TRUE, cpus=18)
sfExport('run_em','true','om')
out <- sfLapply(1:nsim, run_em)
sfStop()
results <- bind_rows(out) %>% group_by(rep, par) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se,
         i=1:n(), re=(est-true)/true,
         year=year0+1:dplyr::n()+ifelse(dplyr::n()==nyears,0,NA),
         age=1:n()+ifelse(n()==nages,0,NA)) %>%
  ungroup
saveRDS(results, file.path(save_folder, 'sim_results_a2.RDS'))
