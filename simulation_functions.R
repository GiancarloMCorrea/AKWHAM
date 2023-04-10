add_ci <- function(g, ci, alpha, showMedian = FALSE, medianCol = '#000000', ...){
  stopifnot(length(ci)==length(alpha))
  for(ii in 1:length(ci))
    g <- g+stat_summary(fun.data = median_hilow, fun.args = list(conf.int = ci[ii]),
                        geom = 'ribbon', alpha = alpha[ii],...)
   if(showMedian) g = g + stat_summary(fun = median, geom = "line", linewidth = 0.5, colour = medianCol)
  return(g)
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

