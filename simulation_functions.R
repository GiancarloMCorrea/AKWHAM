add_ci <- function(g, ci, alpha, ...){
  stopifnot(length(ci)==length(alpha))
  for(ii in 1:length(ci))
    g <- g+stat_summary(fun.data = median_hilow, fun.args = list(conf.int = ci[ii]),
                        geom = 'ribbon', alpha = alpha[ii],...)
  g
}
sim_fn <- function(om, sample.waa=FALSE){
  input <- om$input
  simdata = om$simulate(complete=TRUE)
  ## no simulator for the WAA so do it manually here
  if(sample.waa)
    simdata$waa[simdata$waa_pointer_ssb,,] <-
      matrix(exp(rnorm(520, mean=log(om$rep$pred_waa_ssb), sd=.1)), nrow=52)
  ## do I have to manually ditch missing values for all observations?
  ##  simdata$agg_indices[omdata$agg_indices==-999] <- -999
  input$data <- simdata
  return(input)
}
run_em <- function(i, sample.waa=FALSE){
  library(wham)
  set.seed(i)
  om$retape(FALSE)
  siminput <- sim_fn(om, sample.waa)
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

