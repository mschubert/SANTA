MarkovCentrality <- function(
  g, 
  edge.attr = NULL
) {
  # calculate the importance of each vertex using the Markov Centrality method

  # calculate the mfpt between every vertex pair in both directions
	M <- GraphMFPT(g, edge.attr=edge.attr, average.distances=FALSE)

  # if there are unconnected clusters, the mfpt between unconnected vertices is Inf. Replace these values with a value twice as high as the maximum finite value measured.
	M[M == Inf] <- 2 * max(M[M != Inf])

  C <- 1 / apply(M, 2, mean)
  C
}
