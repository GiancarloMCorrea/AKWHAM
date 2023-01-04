## rm(list = ls())
require(dplyr)
require(ggplot2)
library(wham)
library(snowfall)
add_ci <- function(g, ci, alpha, ...){
  stopifnot(length(ci)==length(alpha))
  for(ii in 1:length(ci))
    g <- g+stat_summary(fun.data = median_hilow, fun.args = list(conf.int = ci[ii]),
                        geom = 'ribbon', alpha = alpha[ii],...)
  g
}
sim_fn <- function(om, self.fit = FALSE){
  input <- om$input
  simdata = om$simulate(complete=TRUE)
  ## do I have to manually ditch missing values for all observations?
  ##  simdata$agg_indices[omdata$agg_indices==-999] <- -999
  input$data <- simdata
  return(input)
}

## Rerun CM to get OM
load("GOA_pollock/fit_a.RData")
## vector (0/1) if 1 then data type (catch, indices, Ecov obs) will be simulated.
fit_a$input$data$simulate_data <- c(1,1,0)
fit_a$input$data$simulate_state <- fit_a$input$data$simulate_state*0
## Looks like these won't matter unless simulate_state is 1 so
## ignore for now
fit_a$input$data$simulate_period

## temporarily turning estimation of selectivity off since it was
## causing some weird convergence failures
fit_a$input$map$logit_selpars <- factor(NA*fit_a$input$par$logit_selpars)
fit_a$input$map$selpars_re <- factor(NA*fit_a$input$par$selpars_re)
om <- fit_wham(fit_a$input, do.osa = FALSE, do.fit = TRUE,
               do.retro = FALSE, n.newton=0)
om$par[which.max(abs(om$gr()))]

## set.seed(123)
## sim <- sim_fn(om, self.fit = FALSE)
## matplot(cbind(om$input$data$agg_catch, sim$data$agg_catch[,1]))
## matplot(cbind(om$input$data$agg_indices[,2], sim$data$agg_indices[,2]))

true <- summary(om$sdrep, select='all')[,1]
run_em <- function(i){
  library(wham)
  set.seed(i)
  om$retape(FALSE)
  siminput <- sim_fn(om)
  tfit <- fit_wham(siminput, do.fit=TRUE, do.osa = F, do.retro = F,
                   MakeADFun.silent = T, n.newton=0)
    ## check for convergence
  ##  if(tfit$opt$convergence!=0) return(NULL)
  ind <- which.max(abs(tfit$gr(tfit$opt$par)))
  stopifnot(ind %in% 1:length(tfit$opt$par))
  maxgrad <- abs(tfit$gr(tfit$opt$par)[ind])
  tmp <- summary(tfit$sdrep, select='all')
  stopifnot( identical(rownames(tmp), names(true)))
  out <- data.frame(rep=i, parnum=seq_along(rownames(tmp)),
                    maxgrad=maxgrad, maxpar=ind,
                    par=rownames(tmp), est=tmp[,1], se=tmp[,2], true=true)
  return(out)
}
## test <- run_em(1)



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


results <- readRDS('GOA_pollock/sim_results_a.RDS')
## ## Check which models failed to converge
## failed <- (1:nsim)[which(! (1:nsim %in% results$rep))]
## if(length(failed>0)) message("failed runs=",failed)
results %>% filter(parnum==maxpar) %>% pull(maxgrad) %>% hist


## (t)ime (s)eries we care about
ts <- filter(results, !is.na(year) & par %in% c('log_F', 'log_SSB'))
## (sc)alar quantities we care about
sc <- filter(results, is.na(year) &  !grepl('NAA|FAA|FXS|resid|F_devs|logit_q_mat|selpars_re|q_re',par))
unique(sc$par)

ggplot(sc, aes(factor(i), re)) + geom_violin() + facet_wrap('par', scales='free')+
  geom_hline(yintercept=0, color=2) + coord_cartesian(ylim=.15*c(-1,1))

g <- ggplot(ts, aes(year,y= re)) + #geom_line() +
  facet_wrap('par', nrow=2)+ geom_hline(yintercept=0, color=2)+
  coord_cartesian(ylim=.15*c(-1,1))
g <- add_ci(g, ci=c(.1,.5,.8), alpha=c(.2,.4,.6))
g
