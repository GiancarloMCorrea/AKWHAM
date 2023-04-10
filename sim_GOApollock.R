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
save_folder = 'GOA_pollock'

### ------------------------------------------------------------
### First the self test of the original penalized likelihood model
## Rerun CM to get OM
load(paste0(save_folder, "/fit_a.RData"))
# Extract some parameters relevant for the simulation:
year0 = fit_a$years[1] - 1
nyears = length(fit_a$years)
nages = length(fit_a$ages.lab)

## vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
fit_a$input$data$simulate_data <- c(1,1,0)
fit_a$input$data$simulate_state <- fit_a$input$data$simulate_state*0 # why turn this off?
## turn off selectivity estimation for all models (not important):
fit_a$input$map$logit_selpars <- factor(NA*fit_a$input$par$logit_selpars)
fit_a$input$map$selpars_re <- factor(NA*fit_a$input$par$selpars_re)
om <- fit_wham(fit_a$input, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton=0)
om$par[which.max(abs(om$gr()))] # check the OM worked ok
om$sdrep
om$input$par <- om$env$parList(par=om$opt$par)

### Run and save the simulations
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
saveRDS(results, file.path(save_folder, 'sim_results_a.RDS'))

### ------------------------------------------------------------
### Test IID random effect smoothing of the WAA
## Rerun CM to get OM
load(paste0(save_folder, "/fit_b.RData"))
## vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
fit_b$input$data$simulate_data <- c(1,1,0)
fit_b$input$data$simulate_state <- fit_b$input$data$simulate_state*0
## turn off selectivity estimation for all models (not important):
fit_b$input$map$logit_selpars <- factor(NA*fit_b$input$par$logit_selpars)
fit_b$input$map$selpars_re <- factor(NA*fit_b$input$par$selpars_re)
om <- fit_wham(fit_b$input, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton=0)
om$par[which.max(abs(om$gr(om$opt$par)))] # check the OM worked ok
om$sdrep

### Run and save the simulations
true <- summary(om$sdrep, select='all')[,1]
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
saveRDS(results, file.path(save_folder, 'sim_results_b.RDS'))

### ------------------------------------------------------------
### Test 2DAR1 random effect smoothing of the WAA
## Rerun CM to get OM
load(paste0(save_folder, "/fit_c.RData"))
## vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
fit_c$input$data$simulate_data <- c(1,1,0)
fit_c$input$data$simulate_state <- fit_c$input$data$simulate_state*0
## turn off selectivity estimation for all models (not important):
fit_c$input$map$logit_selpars <- factor(NA*fit_c$input$par$logit_selpars)
fit_c$input$map$selpars_re <- factor(NA*fit_c$input$par$selpars_re)
om <- fit_wham(fit_c$input, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE, n.newton=0)
om$par[which.max(abs(om$gr(om$opt$par)))] # check the OM worked ok
om$sdrep

### Run and save the simulations
true <- summary(om$sdrep, select='all')[,1]
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
saveRDS(results, file.path(save_folder, 'sim_results_c.RDS'))
