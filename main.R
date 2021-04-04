rm(list = ls())
library(fitdistrplus)
library(parallel)
library(gridExtra)
library(rlist)
library("ggcorrplot")

source("lambda-functions.R")
source("sim-multivar.R")
source("estimation-2-stage.R")
par_trans <- function(x){
  return((x/(1+abs(x))))
}
################################################################################
K = 2
# list of PLP baseline parameters (alpha, beta)
baselineparams.org <- list(c(1,1.1), 
                           c(1,1.1))
t_max = 50
#impact matrix
rho.org = matrix(c(0.4,0.0,
                   0.1,0.2), K,K, byrow = T)
################################################################################
multi_proc <- sim_multivar(t_max=t_max, rho = rho.org,
                           K = K, baselineparams = baselineparams.org)
plot_multivar_FR(multi_proc, rho = rho.org, K = K, baselineparams = baselineparams.org)
################################################################################
set.seed(1)
B <- 10 # number of runs
N <- 10 # number of sequences in each run
rho_est <- list()
stage1_est <- list()
for(b_iter in c(1:B)){
  print(b_iter)
  multi_proclist <- multi_proclist.org <- replicate(N, 
                                                    sim_multivar(t_max=t_max, rho = rho.org, 
                                                                 K = K, baselineparams = baselineparams.org),
                                                    simplify = F)
  ##### stage 1
  
  # use this if baseline is unknown (large N is preferred for better convergence)
  # stage1est <- stage1k(multi_proclist)
  
  # use this if Baseline is considered known
  stage1est <- baselineparams.org
  
  ##### stage 2
  params.init <- rep(0,times = K)
  
  # parallelize estimation
  detectCores()
  cl <- makeCluster(mc <- getOption("cl.cores", 2))
  clusterExport(cl=cl, varlist=ls())
  results <- NA
  
  tryCatch(
    results <-  parSapply(cl, X = c(1:K), FUN = function(k){
      optim(par = params.init,loglik2stageindiintensityk, 
            multi_proclist = multi_proclist, k=k)$par } )
    , 
    error=function(e) {
      e
      print(paste("Oops! --> Error in Loop ",i,sep = ""))
    })
  
  stopCluster(cl)
  results <- par_trans(t(results))
  rho_est[[b_iter]] <- (round(results,5))
  stage1_est[[b_iter]] <- stage1est
  rm(results, cl)
}

length(rlist::list.clean(rho_est))
stage1_est_mean <- apply(stage1_est, c(1,2), mean)
rho_est_mean <- apply(simplify2array(rlist::list.clean(rho_est)), c(1,2), mean)

rownames(rho.org) <- colnames(rho.org) <- paste0("k",c(1:K))
p1 <- ggcorrplot(rho.org, lab = T,
                 colors = c("#F8696B", "#FFFFFF", "#63BE7B"), 
                 ggtheme = ggplot2::theme_dark(),
                 legend.title = "Impact",
                 title = "Truth",
                 tl.srt = 0,
                 method = "square")
rownames(rho_est_mean) <- colnames(rho_est_mean) <- paste0("k",c(1:K))
p2 <- ggcorrplot(rho_est_mean, lab = T,
                 colors = c("#F8696B", "#FFFFFF", "#63BE7B"), 
                 ggtheme = ggplot2::theme_dark(),
                 legend.title = "Impact",
                 lab_size = 3.5,
                 title = "Estimated",
                 tl.srt = 0,
                 method = "square")

grid.arrange(p1, p2, nrow = 1)