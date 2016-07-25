infer_graph_from_datastream_mm <- function(DATA=NULL, prior_on_K=999) {
  
  # Normalise individuals so that they are 1..N - but keep record of original name
  individuals = DATA[, 2]
  relative_index = as.numeric(as.factor(individuals))
  u_relative_index = unique(relative_index)
  u_individual = unique(individuals)
  im_mapping = rep(NA, length(u_relative_index))
  im_mapping[u_relative_index] = u_individual
  DATA[, 2] = relative_index
  
  print("Creating VBGMM...")
  
  x <- DATA[, 1]
  x <- as.matrix(x)
  print("x")
  print(x)
  print("Solving model...")
  temp <- varbayes(x, ncomp = prior_on_K, thres = 0.1, maxit=NULL)
  
  print("Extracting model...")
  model <- extractSimpleModel(model = temp)
  
  
  print("Extracting event probabilities")
  Y <- getVarbayesResp(x, temp)
  #Y <- getVarbayesResp_mod(x, temp)
  
  print("Determining most likely events")
  event <- apply(Y, 1, which.max)
  

  print("Getting unique, occupied, events")
  events <- unique(event)
  
  if (length(events) > 1) {
    print("more than one event!")
  }
  
  
  print("Extracting centroids, variances and weights...")
  centroids = unlist(temp$model$mean[events])
  model_idx = which(unlist(model$mean) %in% centroids)
  variance = unlist(model$cov[model_idx])
  weights = unlist(model$w[model_idx])
  
  
  Y_hard = matrix(0, dim(Y)[1], dim(Y)[2])
  Y_hard[cbind(1:length(event), event)] = 1
  
  # Get rid of empty dimensions
  Y_hard = Y_hard[, events]
  print(dim(Y_hard))
  
  B_hard_incidence_matrix = build_B_from_Y(DATA, Y_hard)
  B_PRIME_hard_incidence_matrix = build_B_PRIME_from_Y(DATA, Y_hard)
  
  A_hard_cooccurrences = get_cooccurrences_in_bipartite_graph(B_hard_incidence_matrix)
  A_PRIME_hard_cooccurrences = get_cooccurrences_in_bipartite_graph(B_PRIME_hard_incidence_matrix)
  
  print("Y_hard")
  print(Y_hard)
  
  EVENT_TIMES = get_T_matrix(Y_hard, DATA)
  print(EVENT_TIMES)
  
  # Create a matrix with first occurrence time for each individual in each event
  first_last = generate_directional_bipartite(EVENT_TIMES, DATA)
  B_first = first_last[['first_occurrence']]
  B_last = first_last[['last_occurrence']]
  
  
  print("B_first")
  print(B_first)
  
  A_directional_cooccurrences = get_directional_cooccurrences_in_bipartite_graph(B_hard_incidence_matrix, B_first)
  A_duration_cooccurrences = get_duration_cooccurrences_in_bipartite_graph(B_hard_incidence_matrix, B_first, B_last)
  
  
  output = list(im_mapping=im_mapping, centroids=centroids, A_hard_cooccurrences=A_hard_cooccurrences, A_PRIME_hard_cooccurrences=A_PRIME_hard_cooccurrences,B_hard_incidence_matrix=B_hard_incidence_matrix, B_PRIME_hard_incidence_matrix=B_PRIME_hard_incidence_matrix, EVENT_TIMES=EVENT_TIMES, A_directional_cooccurrences=A_directional_cooccurrences, A_duration_cooccurrences=A_duration_cooccurrences)
  return(output)
}