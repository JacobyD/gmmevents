gmmevents <- function(DATA, total_individuals, number_of_randomisations=100, prior_on_K = 999) {
  require(VBmix)
  if (is.null(total_individuals)) {
    total_individuals = length(unique(DATA[, 2]))   
  }
  if (number_of_randomisations <= 1) {
    number_of_randomisations = 0
  }
      
    total_locations = length(unique(DATA[, 3]))
    locations = unique(DATA[, 3])
    location_list = list()
    
    A = matrix(0, total_individuals, total_individuals)
    A_DUR = matrix(0, total_individuals, total_individuals)
    A_PRE = matrix(0, total_individuals, total_individuals)
    A_PRIME = matrix(0, total_individuals, total_individuals)
    A_DIR = matrix(0, total_individuals, total_individuals)
    
    if (number_of_randomisations > 0) {
      Anull_mean = matrix(0, total_individuals, total_individuals)
      Anull_std = matrix(0, total_individuals, total_individuals)
    
      A_PRIME_null_mean = matrix(0, total_individuals, total_individuals)
      A_PRIME_null_std = matrix(0, total_individuals, total_individuals)
    } else {
    
      Anull_mean = NaN
      Anull_std = NaN
      A_PRIME_null_mean = NaN
      A_PRIME_null_std = NaN
    }
    
    
    ids <- DATA[, 2]
    X <- hist(ids, breaks=total_individuals, plot=F)
    
    print("histogram of individuals")
    print(X$breaks)
    print(X$counts)
    
    B <- list()
    B_PRIME <- list()
    
    for (i in 1:total_locations) {
      
      location_index = locations[i]
      
      print(sprintf("Location index %d", location_index))
      print("-----------\n")
      
      DATA_local_worker_copy = DATA
      
      location_indices = DATA_local_worker_copy[, 3] == location_index
      
      if (sum(location_indices) == 0) {
        next
      }
      
      # If there is only a single individual in this location, skip it
      if (length(unique(DATA_local_worker_copy[location_indices, 2])) < 2) {
        next
      }
      
      
      location_list[[i]] = location_index
      
      DATA_LOC = DATA_local_worker_copy[location_indices, ]
      
      # If there's less than 5 data points (this could more objective?), don't run for this location
      if (length(DATA_LOC[, 1]) < 5) {
        print("Skipping location as not enough data (<5)")
        next
      }
      
      output = infer_graph_from_datastream_mm(DATA_LOC, prior_on_K)
      
      print_results(output)
      
      
      A_LOC = output[['A_hard_cooccurrences']]
      A_DIR_LOC = output[['A_directional_cooccurrences']]
      A_DUR_LOC = output[['A_duration_cooccurrences']]
      B_LOC = output[['B_hard_incidence_matrix']]
      A_PRIME_LOC = output[['A_PRIME_hard_cooccurrences']]
      B_PRIME_LOC = output[['B_PRIME_hard_incidence_matrix']]
    
      # Hang on to pre-sig As
      A_LOC_PRE = A_LOC
      
      #DO SIGNIFICANCE TEST HERE
      if (sum(A_LOC) > 0) {
        if (number_of_randomisations > 0) {
          sig_result = do_significance_test_of_adjancency_given_incidence_matrix(A_LOC, B_LOC, .05, number_of_randomisations)
          sig_prime_result = do_significance_test_of_adjancency_given_incidence_matrix(A_PRIME_LOC, B_PRIME_LOC, .05, number_of_randomisations)
          
          
          A_LOC = sig_result[['W_significant']]
          A_NULL_LOC_MEAN = sig_result[['W_null_mean']]
          A_NULL_LOC_STD = sig_result[['W_null_std']]
          A_PRIME_LOC = sig_prime_result[['W_significant']]
          A_NULL_PRIME_LOC_MEAN = sig_prime_result[['W_null_mean']]
          A_NULL_PRIME_LOC_STD = sig_prime_result[['W_null_std']]
        }
        
      }
      
      
      # Set NaNs to 0
      A_LOC_PRE[is.nan(A_LOC_PRE)] = 0
      A_LOC[is.nan(A_LOC)] = 0
      A_PRIME_LOC[is.nan(A_PRIME_LOC)] = 0
      
      print("[POST SIG] Adjacency Matrix (including within event detection counts) 'A_LOC'")
      print(A_LOC)
      print("[POST SIG] Adjacency Matrix (count of shared events) 'A_PRIME_LOC'")
      print(A_PRIME_LOC)
      
      # Get original indices index manager
      im_mapping = output[['im_mapping']] 
      
      for (i in 1:length(im_mapping)) {
        for (j in 1:length(im_mapping)) {
          A_PRE[im_mapping[i], im_mapping[j]] = A_PRE[im_mapping[i], im_mapping[j]] + A_LOC_PRE[i, j]
          A_DUR[im_mapping[i], im_mapping[j]] = A_DUR[im_mapping[i], im_mapping[j]] + A_DUR_LOC[i, j]
          A[im_mapping[i], im_mapping[j]] = A[im_mapping[i], im_mapping[j]] + A_LOC[i, j]
          A_PRIME[im_mapping[i], im_mapping[j]] = A_PRIME[im_mapping[i], im_mapping[j]] + A_PRIME_LOC[i, j]
          A_DIR[im_mapping[i], im_mapping[j]] = A_DIR[im_mapping[i], im_mapping[j]] + A_DIR_LOC[i, j]
        }
      }
      
      B_aux = matrix(0, total_individuals, nrow(B_LOC))

      B_aux[im_mapping, ] = t(B_LOC)
      
      B[[paste(location_index)]] = B_aux
      
      B_PRIME_aux = matrix(0, total_individuals, nrow(B_PRIME_LOC))
      B_PRIME_aux[im_mapping, ] = t(B_PRIME_LOC)
      B_PRIME[[paste(location_index)]] = B_PRIME_aux
      
      if (sum(A_LOC) > 0) {
        if (number_of_randomisations > 0) {
          A_NULL_LOC_MEAN[is.nan(A_NULL_LOC_MEAN)] = 0
          A_NULL_LOC_STD[is.nan(A_NULL_LOC_STD)] = 0
      
          A_NULL_PRIME_LOC_MEAN[is.nan(A_NULL_PRIME_LOC_MEAN)] = 0
          A_NULL_PRIME_LOC_STD[is.nan(A_NULL_PRIME_LOC_STD)] = 0
      
          
          for (i in 1:length(im_mapping)) {
            for (j in 1:length(im_mapping)) {
              Anull_mean[im_mapping[i], im_mapping[j]] = Anull_mean[im_mapping[i], im_mapping[j]] + A_NULL_LOC_MEAN[i, j]
              Anull_std[im_mapping[i], im_mapping[j]] =  Anull_std[im_mapping[i], im_mapping[j]] + (A_NULL_LOC_STD[i, j]^2)
          
              A_PRIME_null_mean[im_mapping[i], im_mapping[j]] = A_PRIME_null_mean[im_mapping[i], im_mapping[j]] + A_NULL_PRIME_LOC_MEAN[i, j]
              A_PRIME_null_std[im_mapping[i], im_mapping[j]] = A_PRIME_null_std[im_mapping[i], im_mapping[j]] + (A_NULL_PRIME_LOC_STD[i, j]^2)
            }
          }
        }
      }
    }

    if (number_of_randomisations > 0) {
      Anull_std = sqrt(Anull_std)
      A_PRIME_null_std = sqrt(A_PRIME_null_std)
    }
  
    results = list(AdjMat = A,
      AdjMat_PRIME = A_PRIME,
      BPGraph = B,
      BPGraph_PRIME = B_PRIME,
      num_occ_individuals = X,
      Anull_mean = Anull_mean,
      A_PRIME_null_mean = A_PRIME_null_mean,
      Anull_std = Anull_std,
      A_PRIME_null_std = A_PRIME_null_std,
      locations =  location_list,
      A_DIR =  A_DIR,
      A_PRE = A_PRE,
      A_DUR = A_DUR)
  
    return(results)

}