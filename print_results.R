print_results <- function(output) {
  
  print("[PRE SIG] Infer Graph Output... 'output'")

  A_LOC = output[['A_hard_cooccurrences']]
  
  print("[PRE SIG] Adjacency Matrix (including within event detection counts) 'A_LOC'")
  print(A_LOC)
  
  A_PRIME_LOC = output[['A_PRIME_hard_cooccurrences']]
  
  print("[PRE SIG] Adjacency Matrix (count of shared events) 'A_PRIME_LOC'")
  print(A_PRIME_LOC)
  
  B_LOC = output[['B_hard_incidence_matrix']]
  
  print("[PRE SIG] BP Graph (including within event detection counts) 'B_LOC'")
  print(B_LOC)

  A_DIR_LOC = output[['A_directional_cooccurrences']]
  print("[PRE SIG] Directional Adjacency Matrix (including within event detection counts) 'A_DIR_LOC'")
  print(A_DIR_LOC)
  
  B_PRIME_LOC = output[['B_PRIME_hard_incidence_matrix']]
  
  print("[PRE SIG] BP Graph (count of occupied events) 'B_PRIME_LOC'")
  print(B_PRIME_LOC)
  
  
}