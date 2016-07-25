
do_significance_test_of_adjancency_given_incidence_matrix = function(W, B, a_threshold, number_of_randomisations) {
  
  B_dims = dim(B)
  K = B_dims[1]
  N = B_dims[2]
  
  P_values_left = matrix(0, N, N)
  P_values_right = matrix(0, N, N)
  
  #Wall = numpy.zeros((N ** 2, number_of_randomisations))
  someData <- rep(NaN, N*N*number_of_randomisations)
  Wall <- array(someData, c(N, N, number_of_randomisations)) 
  
  for (rand_index in 1:number_of_randomisations) {
      # Start with bipartite graph
      Bnull = B
      # For each individual
      for (n in 1:N) {
        # Get the set of events that they are included in
        visitation_profile = Bnull[, n]
      
        # Work out how many occurrences for that individual
        occurences_n = sum(visitation_profile)
        
        # Normalise each events frequency by the total number of occurances (proportion of visits in each event)
        pi_k = visitation_profile / occurences_n
        
        # Randomise this set of proportions
        pi_k = pi_k[sample(K)]  
        
        
        # From the random set of proportions (probability of being seen in an event) - sample occurences_n times to generate a random set of occurences. And put this in Bnull in the nth column.
        # This will build up a new bipartite graph, mainining individual frequency across locations, but with those locations randomised
        Bnull[, n] = rmultinom(1, occurences_n, prob=pi_k)
        
      }
    
    # Build adjacency matrix from this single random bipartite graph
    Wnull = get_cooccurrences_in_bipartite_graph(Bnull)
    
    
    # Put this random adjacency into the big array
    #Wall[:, rand_index] = numpy.ndarray.flatten(Wnull)
    Wall[, , rand_index] = Wnull
  }
  
  for (i in 1:(N-1)) {
    # For each other individual
    for (j in (i+1):N) {
    
      # Extract the set of all random cooccurences for individuals i and j
      # And append the 'true/observed' occurences on the end
      # This is a random sample of 'number_of_randomisations' possible
      # cooccurences with the observed value on the end
      #a_vector = numpy.append(Wall[k, :], W[i, j])
     
      a_vector = c(Wall[i, j, ], W[i, j])
      
      # What's the maximum number of cooccurences in this random sample
      max_interaction = max(a_vector)
      
      # Get bins [0 - max_interaction]
      bins = seq(0, max_interaction + 1)
      # If there's no bins, just make [0, 1]
      if (length(bins) < 1) { 
        bins = c(0, 1)
      }    
      # Put the set of coocurrences into these bins (histogram)
      result = hist(Wall[i, j, ], breaks=bins, plot=F)
      
      # Get the histogram out (result[1] is the bins)
      h = result$counts
      
      # Normalise to max==1
      h = h / sum(h)
      
      # For this pair of individuals
      # What proportion of our random cooccurences lie to the left
      P_values_left[i, j] = sum(h[1:(W[i, j] + 1)])
      # And right or our observed value
      P_values_right[i, j] = sum(h[(W[i, j] + 1):length(h)])
    }
  }
  # Initial matrices are zeros or empties
  P_significant = matrix(0, N, N)
  P_values = matrix(NaN, N, N)
  
  # For each pair of individuals
  for (i in 1:N) {
    for (j in i:N) {
    
      # This would be if the frequency of cooccurrences are more *or less* likely than random (sig. preference and avoidance)
      #?#P_significant(i,j) = (P_values_left(i,j)<a_threshold) || (P_values_right(i,j)<a_threshold);
      
      # If the proportion of random coocurrences bigger (P_values_right) 
      # then the threshold, then this cooccurence is signficant?
      # Are the observed coocurrences more frequent that random
      P_significant[i, j] = (P_values_right[i, j] < a_threshold)
    
      # If pref or avoid, store approporate p-value
      #?#P_values(i,j) = min(P_values_left(i,j),P_values_right(i,j));
      
      # For just pref, store this p-value
      P_values[i, j] = P_values_right[i, j]
    }
  }
    
  # Make matrix symmetrical
  P_significant = P_significant + t(P_significant)
  
  # Extract only significant associasions (ones with P==0 will be lost)
  W_significant = W * P_significant
  
  # For which adjacencies do we have a significant value
  A = W != 0
  # For those, get the p-value (again values where A==0 will lost)
  P_significant = P_significant * A
  
  #?#original_number_of_links = sum(get_triu_vector(A));
  #?#significant_number_of_links = sum(get_triu_vector(P_significant));
  #?#prune_factor = 1 - significant_number_of_links/original_number_of_links;
  
  # Determine the mean frequency of cooccurrence for all individuals
  #W_null_mean = numpy.mean(Wall, 1)
  #W_null_mean = mean(Wall, 2)
  W_null_mean <- apply(Wall, c(1,2), mean)
  #?#print "W_null_mean"
  #?#print W_null_mean.shape
  # Make it and N*N matrix
  #W_null_mean = numpy.reshape(W_null_mean, (N, N))
  
  # Similarly for std
  #?#print numpy.transpose(Wall).shape
  #W_null_std = numpy.std(Wall, 2)
  W_null_std <- apply(Wall, c(1,2), sd)
  #?#print "W_null_std"
  #?#print W_null_std.shape
  #W_null_std = numpy.reshape(W_null_std, (N, N))
  
  # Return everything...
  result = list(W_significant=W_significant, W_null_mean=W_null_mean, W_null_std=W_null_std, P_values=P_values, P_significant=P_significant, Wall=Wall)

  return(result)
}