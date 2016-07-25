get_directional_cooccurrences_in_bipartite_graph <- function(B, B_first) {

    # assumes [K x N] incidence matrix - K clusters N individuals
    N <- dim(B)[2]
    K <- dim(B)[1]
    
    # Adjacency matrix
    W <- matrix(0, N, N)
    
    # For each indivudal
    for(i in 1:(N - 1)) {
      # For each other individual
      for(j in (i + 1):N) {
        # How frequently these individuals are in each location
        visitation_profile_i <- B[, i]
        visitation_profile_j <- B[, j]
  
        # to vectorise
        for(k in 1:K) {
          cooccurrence = min(visitation_profile_i[k], visitation_profile_j[k])
          if (cooccurrence >= 1) {
            if (B_first[k, i] < B_first[k, j]) {
              W[i, j] = W[i, j] + 1
            } else {
              W[j, i] = W[j, i] + 1
            }
          }
        }
      }
    }
    return(W) 
}



