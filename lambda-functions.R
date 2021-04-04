base_lambdak <- function(t, k, baselineparams) {
  alpha  = baselineparams[[k]][1]
  beta  = baselineparams[[k]][2]
  alpha*beta*(t^(beta-1))
}
base_Lambdak <- function(t, k, baselineparams) {
  alpha  = baselineparams[[k]][1]
  beta  = baselineparams[[k]][2]
  alpha*(t^(beta))
}
lambdak <- function(t, t_vec, rho, compo_ind, k, baselineparams){
  k=k
  tmpsum <- 0
  t_vec.org <- t_vec
  compo_ind.org <- compo_ind
  t_vec <- t_vec[t_vec.org < t]
  compo_ind <- compo_ind[t_vec.org < t]
  Nt <- length(t_vec)
  # Nt <- sum(t_vec < t)
  if(Nt == 0) return(base_lambdak(t, k, baselineparams))  else{
    for(j in c(0:(Nt-1))){
      # j <- j_iter - 1
      if(j == 0){
        tmpprod <- 1
        r <- rho[k,compo_ind[(Nt-j)]]
        T_Nt_minus_j <- t_vec[(Nt-j)]
        tmpsum <- tmpsum + (r)*tmpprod*base_lambdak(T_Nt_minus_j, k, baselineparams)
      }else{
        r <- rho[k,compo_ind[(Nt-j)]]
        tmpprod <- 1
        for(r_iter in c((Nt-j+1):(Nt))){
          tmpr <- rho[k,compo_ind[(r_iter)]]
          tmpprod <- tmpprod*(1-tmpr)
        }
        T_Nt_minus_j <- t_vec[(Nt-j)]
        tmpsum <- tmpsum + (r)*tmpprod*base_lambdak(T_Nt_minus_j, k, baselineparams)
      }
      
    }
    lam <- base_lambdak(t, k, baselineparams) - tmpsum
    return(lam)
  }
}
Lambdak <- function(t,t_vec, t_max = max(t_vec), rho, compo_ind, k, baselineparams){
  if(length(t_vec[t_vec < t]) == 0){
    return(base_Lambdak(t, k, baselineparams))
  }else{
    t_vec_tmp <- t_vec[t_vec <= t]
    t_vec_tmp <- c(t_vec_tmp,t); t_vec_tmp <- unique(t_vec_tmp)
    compo_ind_tmp <- compo_ind[t_vec <= t]
    B <- sum(sapply(t_vec_tmp, FUN = function(x){
      # Nt <- length(t_vec)
      tmpsum <- 0
      Nt <- sum(t_vec_tmp < x)
      if(Nt == 0) return(0)  else{
        for(j in c(0:(Nt-1))){
          # j <- j_iter - 1
          if(j == 0){
            tmpprod <- 1
            r <- rho[k,compo_ind_tmp[(Nt-j)]]
            T_Nt_minus_j <- t_vec_tmp[(Nt-j)]
            tmpsum <- tmpsum + (r)*tmpprod*base_lambdak(T_Nt_minus_j, k, baselineparams)
          }else{
            r <- rho[k,compo_ind_tmp[(Nt-j)]]
            tmpprod <- 1
            for(r_iter in c((Nt-j+1):(Nt))){
              tmpr <- rho[k,compo_ind_tmp[(r_iter)]]
              tmpprod <- tmpprod*(1-tmpr)
            }
            T_Nt_minus_j <- t_vec_tmp[(Nt-j)]
            tmpsum <- tmpsum + (r)*tmpprod*base_lambdak(T_Nt_minus_j, k, baselineparams)
          }
          
        }
        tmpsum*((min(t,t_vec_tmp[Nt+1], na.rm = T) - t_vec_tmp[Nt]))
      }
    } ))
    
    return(base_Lambdak(t, k, baselineparams) - B)
  }
}