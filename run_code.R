if (!require(VBmix)) install.packages("VBmix")

# Load data from file
DATA <- read.csv("4000.csv", header=T)
print("Raw_IDs")

# Print out individual IDs (animal IDs)
Raw_IDs <- unique(DATA[,2])
print(Raw_IDs)

# Print out locations (location IDs)
print("Locations..")
Raw_LOCs <- unique(DATA[,3])
print(Raw_LOCs)

# Re-order the IDs appropriately (chaning from unique IDs to numeric IDs from 1 to N)
DATA[, 2] = as.numeric(factor(DATA[, 2]))
# Save the numerics IDs
num_ids <- unique(as.numeric(factor(DATA[, 2])))
# And save the Raw/original IDs so we can put them back together later
Raw_IDs <- Raw_IDs[order(num_ids)]

# Source all code files
source("build_B_from_Y.R")
source("build_B_PRIME_from_Y.R")
source("generate_directional_bipartite.R")
source("get_cooccurrences_in_bipartite_graph.R")
source("get_directional_cooccurrences_in_bipartite_graph.R")
source("get_duration_cooccurrences_in_bipartite_graph.R")
source("get_T_matrix.R")
source("infer_graph_from_datastream_mm.R")
source("gmmevents.R")
source("do_significance_test_of_adjacency_given_incidence_matrix.R")
source("print_results.R")

# Run main function...
# specify number of individuals, number of randomisations and prior on K
# prior_on_K is initial number of events for the GMM in each location (default = 999)
results = gmmevents(DATA, total_individuals = length(Raw_IDs), number_of_randomisations=1000, prior_on_K = 999)

# Extract resulting matrices from 'results' list
AdjMat_count = results[['AdjMat']]
AdjMat_dur = results[['A_DUR']]
AdjMat_pre_sig_count = results[['A_PRE']]
AdjMat_PRIME = results[['AdjMat_PRIME']]
AdjMat_PRIME_dir = results[['A_DIR']]
Summary = results[['num_occ_individuals']]

# Get number of events by location
BPGraphs <- results$BPGraph[lapply(results$BPGraph, length) > 0] 
print("Number of events in locations with >1 individual")
nevents <- unlist(lapply(BPGraphs, ncol))
nevents

# Plot some images of Adjacency matrices
image(as.matrix(AdjMat_count))
image(as.matrix(AdjMat_PRIME))
image(as.matrix(AdjMat_PRIME_dir))
image(as.matrix(AdjMat_dur))

# Write out adjacency matrices to csv files
write.table(AdjMat_PRIME, file = "AdjMat_PRIME_4000.csv", sep=",", row.names = Raw_IDs, col.names = NA)
write.table(AdjMat_count, file = "AdjMat_count_4000.csv", sep=",", row.names = Raw_IDs, col.names = NA)
write.table(AdjMat_dur, file = "AdjMat_dur_4000.csv", sep=",", row.names = Raw_IDs, col.names = NA)
write.table(AdjMat_pre_sig_count, file = "AdjMat_pre_sig_count_4000.csv", sep=",", row.names = Raw_IDs, col.names = NA)
write.table(AdjMat_PRIME_dir, file = "AdjMat_PRIME_dir_4000.csv", sep=",", row.names = Raw_IDs, col.names = NA)

# Crude network plot
library(igraph)

# Constuct network from adjacency matrix
net <- graph_from_adjacency_matrix(AdjMat_PRIME,mode="undirected",weighted=TRUE,diag=FALSE)
# Specify layout
lay <- layout.fruchterman.reingold(net)
# Alternative layout
lay2 <- layout.circle(net)
# Plot network using layout
plot(net, layout=lay, vertex.label = Raw_IDs, edge.color = "grey60",
     edge.width=(E(net)$weight/(max(E(net)$weight))*10)) 
#NOTE Use different weighting argument if comparing multiple graphs with different max edge weights (e.g. PRIME v's count)
