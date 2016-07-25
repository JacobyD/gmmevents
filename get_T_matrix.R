get_T_matrix <- function(Y, DATA) {

  # Ensure Y_hard is matrix
  Y = as.matrix(Y)
	
  K <- dim(Y)[2]
  print(paste("K:", K))
  T = matrix(0, K, 2)
  
  for(k in 1:K) {
    
    print(paste("little k:", k))
    first_obs <- which(Y[, k]!=0)[1]
    last_obs <- tail(which(Y[, k]!=0), 1)
    
    T[k, 1] = DATA[first_obs, 1]
    T[k, 2] = DATA[last_obs, 1]
  }
  return(T)
}
