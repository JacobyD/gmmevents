generate_directional_bipartite <- function(event_times, DATA) {
  
    individuals <- unique(DATA[, 2])
    N <- length(individuals)
    
    #py K = len(event_times)
    K <- dim(event_times)[1]
    
    print(individuals)
    print(event_times)
    
    # Create matrix to hold first occurrence times for each ind
    B_FIRST_TIME = matrix(0, K, N)
    B_LAST_TIME = matrix(0, K, N)
    
    for(k in 1:K) {
      start = event_times[k, 1]
      end = event_times[k, 2]
      
      # Find 0/1 for those data after start and before end
      #py mask = numpy.nonzero((DATA[:, 0] >= start) & (DATA[:, 0] <= end))
      mask = (DATA[,1] >= start) & (DATA[,1] <= end)
      #print mask.shape
      # Get all that data
      DATA_IN_EVENT = DATA[mask, ]
      
      # Get the individuals in this data
      individuals = unique(DATA_IN_EVENT[, 2])
      
      # For each individual
      for (ind in individuals) {
        # Get only thet data for that individual
        mask = DATA_IN_EVENT[, 2] == ind
        THIS_IND_DATA_IN_EVENT = DATA_IN_EVENT[mask, ]
        
        # The first occurance in this event is the minimum time
        first_occurrence = min(THIS_IND_DATA_IN_EVENT[, 1])
        last_occurrence = max(THIS_IND_DATA_IN_EVENT[, 1])
        
        #print("first_occurrence")
        #print(first_occurrence)
      
        # Put in matrix
        B_FIRST_TIME[k, ind] = first_occurrence
        B_LAST_TIME[k, ind] = last_occurrence
      }
    }
    return(list(first_occurrence=B_FIRST_TIME, last_occurrence=B_LAST_TIME))
} 
