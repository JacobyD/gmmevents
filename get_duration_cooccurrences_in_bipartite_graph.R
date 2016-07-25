get_duration_cooccurrences_in_bipartite_graph <- function(B, B_first, B_last) {
  
    # assumes [K x N] incidence matrix - K clusters N individuals
    N <- dim(B)[2]
    K <- dim(B)[1]
    
    # Adjacency matrix
    W <- matrix(0, N, N)
    
    # For each individual
    for(i in 1:(N - 1)) {
      # For each other individual
      for(j in (i + 1):N) {
        # How frequently these individuals are in each location
        visitation_profile_i <- B[, i]
        visitation_profile_j <- B[, j]
        
        # to vectorise
        for(k in 1:K) {
          cooccurrence = min(visitation_profile_i[k], visitation_profile_j[k])
          #print("cooccurrence")
          #print(cooccurrence)
          if (cooccurrence >= 1) {
            last_start = max(B_first[k, i], B_first[k, j])
            first_end = min(B_last[k, i], B_last[k, j])
            
            duration = first_end - last_start
            if (duration > 0) {
              W[i, j] = W[i, j] + duration
            } 
          }
        }
      }
    }
    W <- W + t(W)
    return(W)
}