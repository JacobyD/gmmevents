build_B_from_Y <- function(DATA, Y_hard) {

  # Ensure Y_hard is matrix
  Y_hard = as.matrix(Y_hard)
  
  individuals = unique(DATA[, 2])
  N_individuals = length(individuals)
  
  K = dim(Y_hard)[2]
  
  B = matrix(0, K, N_individuals)
  
  for (i in 1:N_individuals) {
    i_indices = DATA[, 2] == individuals[i]
    i_data <- as.matrix(Y_hard[ i_indices,  ])
    if (sum(i_indices) > 1) {
      B[, individuals[i]] = colSums(i_data, dims = 1)
    } else {
      B[, individuals[i]] = i_data
    }
  }  
  
   
  return(B)
}