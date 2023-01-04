## rm(list = ls())
require(dplyr)
require(ggplot2)
library(wham)
library(snowfall)
source("simulation_functions.R")


### ------------------------------------------------------------
### First the self test of the original penalized likelihood model
## Rerun CM to get OM
load("GOA_pollock/fit_a.RData")
## vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
fit_a$input$data$simulate_data <- c(1,1,0)
fit_a$input$data$simulate_state <- fit_a$input$data$simulate_state*0
## temporarily turning estimation of selectivity off since it was
## causing some weird convergence failures
fit_a$input$map$logit_selpars <- factor(NA*fit_a$input$par$logit_selpars)
fit_a$input$map$selpars_re <- factor(NA*fit_a$input$par$selpars_re)
om <- fit_wham(fit_a$input, do.osa = FALSE, do.fit = TRUE,
               do.retro = FALSE, n.newton=0)
om$par[which.max(abs(om$gr()))]         # check the OM worked ok

### Run and save the simulations
## test <- run_em(1)
true <- summary(om$sdrep, select='all')[,1]
nsim <- 100
sfInit(parallel=TRUE, cpus=5)
sfExport('run_em','sim_fn','true', 'om')
out <- sfLapply(1:nsim, run_em)
sfStop()
results <- bind_rows(out) %>% group_by(rep, par) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se,
         i=1:n(), re=(est-true)/true,
         year=1969+1:n()+ifelse(n()==52,0,NA),
         age=1:n()+ifelse(n()==10,0,NA)) %>%
  ungroup
saveRDS(results, 'GOA_pollock/sim_results_a.RDS')

### Process and plot them
results <- readRDS('GOA_pollock/sim_results_a.RDS')
results %>% filter(parnum==maxpar) %>% pull(maxgrad) %>% summary
## (t)ime (s)eries we care about
ts <- filter(results, !is.na(year) & par %in% c('log_F', 'log_SSB'))
## (sc)alar quantities we care about
sc <- filter(results, is.na(year) &  !grepl('NAA|FAA|FXS|resid|F_devs|logit_q_mat|selpars_re|q_re',par))
g1 <- ggplot(sc, aes(factor(i), re)) + geom_violin() + facet_wrap('par', scales='free')+
  geom_hline(yintercept=0, color=2) + coord_cartesian(ylim=.15*c(-1,1))
g <- ggplot(ts, aes(year,y= re)) + #geom_line() +
  facet_wrap('par', nrow=2)+ geom_hline(yintercept=0, color=2)+
  coord_cartesian(ylim=.15*c(-1,1))
g2 <- add_ci(g, ci=c(.1,.5,.8), alpha=c(.2,.4,.6))
g1
g2
### End of the fit_a simulation
### ------------------------------------------------------------
