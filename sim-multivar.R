sim_multivar <- function(t_max, rho, K, baselineparams){
  
  t <- 0
  X <- numeric()
  compo_ind <- numeric()
  
  tk_fn <- function(k) rweibull(1, shape = baselineparams[[k]][2], scale = baselineparams[[k]][1])
  
  while(t <= t_max){
    u <- runif(1)
    if(length(X) == 0){
      t1s <- sapply(c(1:K), tk_fn)
      t <- min(t1s)
      X <- c(X,t)
      compo_ind <- c(compo_ind, which.min(t1s))
      rm(t1s)
    }else{
      lambdak_star <- sapply(c(1:K), lambdak, t = t, t_vec = X, compo_ind = compo_ind, 
                             rho = rho, baselineparams = baselineparams)
      sum_lambdas <- sum(lambdak_star)
      
      w <- -1*(log(u)/sum_lambdas)
      t <- t + w
      stopifnot(t>0)
      d <- runif(1)
      
      if(d*sum_lambdas <= ( sum(sapply(c(1:K), lambdak, t = t, t_vec = X, compo_ind = compo_ind, 
                                       rho = rho, baselineparams = baselineparams)) )){
        tmplams <- sapply(c(1:K), lambdak, t = t, t_vec = X, compo_ind = compo_ind, 
                         rho = rho, baselineparams = baselineparams)
        tmpCI <- which( cumsum(tmplams) > (d*sum_lambdas) )[1]
        X <- c(X, t)
        compo_ind <- c(compo_ind, tmpCI)
      }
    }
    
  }
  return(list(X[-c(length(X))], compo_ind[-c(length(compo_ind))]))
}
plot_multivar_FR <- function(multi_proc, rho, K, baselineparams){
  
  Vplotlamk <- Vectorize(lambdak, vectorize.args = "t")
  
  t_vec <- multi_proc[[1]]
  compo_ind <- multi_proc[[2]]
  delT <- max(t_vec)/500 # discretize time axis
  t_series <- seq(0, max(t_vec), by = delT)
  par(mfrow=c(K,2))
  for(k_iter in c(1:K)){
    plot(t_series, Vplotlamk(t=t_series, t_vec = t_vec, rho = rho, compo_ind = compo_ind, k = k_iter,
                             baselineparams = baselineparams),type = "l", ylab = paste0("lambda_",k_iter),xlab="t",col = "green")
    
    plot(t_series, delT*cumsum(Vplotlamk(t=t_series, t_vec = t_vec, rho = rho, compo_ind = compo_ind, k = k_iter,
                                         baselineparams = baselineparams)),type = "l", ylab = paste0("Lambda_",k_iter),xlab="t",col = "green")
    
  }
  
}