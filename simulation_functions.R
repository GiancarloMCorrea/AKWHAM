add_ci <- function(g, ci, alpha, ...){
  stopifnot(length(ci)==length(alpha))
  for(ii in 1:length(ci))
    g <- g+stat_summary(fun.data = median_hilow, fun.args = list(conf.int = ci[ii]),
                        geom = 'ribbon', alpha = alpha[ii],...)
  g
}

run_em <- function(i){
  library(wham)
  set.seed(i)
  om$retape(FALSE)
  input <- om$input
  simdata = om$simulate(complete=TRUE)
  input$data <- simdata
  siminput <- input
  tfit <- fit_wham(siminput, do.fit=TRUE, do.osa = FALSE, do.retro = FALSE,
                   MakeADFun.silent = TRUE, n.newton=0)
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


run_em2 <- function(i){
  library(wham)
  library(dplyr)
  set.seed(i)
  om$retape(FALSE)
  input <- om$input
  simdata = om$simulate(complete=TRUE)
  input$data <- simdata
  siminput <- input
  true <- om$par
  rm(om);
  gc() ## try to svae on memory
  tfit <- fit_wham(siminput, do.fit=TRUE, do.osa = FALSE, do.retro = FALSE,
                   MakeADFun.silent = 0, do.sdrep=FALSE, n.newton=0)
  ## restart optimizer and do newton step to get convergence better
  siminput$par <- tfit$parList
  tfit <- fit_wham(siminput, do.fit=TRUE, do.osa = FALSE, do.retro = FALSE,
                   MakeADFun.silent = 0, do.sdrep=FALSE, n.newton=1)
  rep <- tfit$report(tfit$par)
  ##  if(tfit$opt$convergence!=0) return(NULL)

  ind <- which.max(abs(tfit$gr(tfit$opt$par)))
  stopifnot(ind %in% 1:length(tfit$opt$par))
  maxgrad <- abs(tfit$gr(tfit$opt$par)[ind])
  ssb <- data.frame(par='SSB', true=simdata$SSB, est=tfit$rep$SSB, year=input$years)
  catch <- data.frame(par='pred_catch', true=simdata$pred_catch,
                      est=tfit$rep$pred_catch, year=input$years)
  rec <- data.frame(par='recruit', true=simdata$NAA[,1], est=tfit$rep$NAA[,1], year=input$years)
  f <- data.frame(par='F', true=simdata$F, est=tfit$rep$F, year=input$years)
  sc <- data.frame(par=names(true), true=true,
                   est=tfit$opt$par, year=NA)
  out <- rbind(ssb,catch,f,rec,sc) %>% cbind(rep=i, maxgrad=maxgrad, maxpar=ind)
  return(out)
}

