fitweibcenswrapper <- function(t,ind, curr_ind){
  tleft <- t
  tright <- ifelse(ind == curr_ind, t, NA)
  censdf <- data.frame(left = tleft, right = tright)
  
  fitdistcens(as.data.frame(censdf), distr = "weibull")$estimate
}
stage1k <- function(multi_proclist){
  t_i_fun <- function(multi_proc){
    t_i1 <- multi_proc[[1]][1]
  }
  compo_ind_i_fun <- function(multi_proc){
    compo_ind <- multi_proc[[2]][1]
  }
  t <- unlist(lapply(multi_proclist, FUN = t_i_fun))
  ind <- unlist(lapply(multi_proclist, FUN = compo_ind_i_fun))
  
  fit.para.k.tmp <- t( sapply(c(1:K) ,fitweibcenswrapper, t = t, ind = ind))
  fit.para.k.tmp2 <- fit.para.k.tmp[,c(2,1)]
  colnames(fit.para.k.tmp2) <- NULL
  
  fit.para.k <- lapply(apply(fit.para.k.tmp2, 1, function(x) list(c(x[1], x[2]))), function(y) unlist(y))
  return((fit.para.k))
}
loglik2stageindiintensityk <- function(params, multi_proclist, k){
  
  k=k
  rho <- matrix(0, K, K)
  rho[k,] <- par_trans(params)
  
  l_i_fun <- function(multi_proc){
    t_vec <- multi_proc[[1]]
    compo_ind <- multi_proc[[2]]
    maxT <- max(t_vec)
    if(length(unique(compo_ind))  < K) return(0)
    if(var(compo_ind)==0) return(0)
    termk1 <- sum(log(sapply(t_vec[compo_ind==k], lambdak,t_vec=t_vec, 
                             compo_ind = compo_ind, rho=rho, k=k, baselineparams = stage1est)))
    
    termk2 <- Lambdak(t = maxT, t_vec = t_vec, rho = rho, compo_ind = compo_ind, k=k, baselineparams = stage1est)
    return((-(termk1-termk2)))
    
  }
  t <- unlist(lapply(multi_proclist, FUN = l_i_fun))
  nll <- sum(t[t!=0])
  return(nll)
}