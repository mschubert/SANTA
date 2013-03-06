GraphMFPT <- function(
  g,
  v=V(g),
  edge.attr=NULL,
  average.distances=TRUE
) {
  # calculates the mean first-passage time mfpt between every vertex pair
  # algorithms for Estimating Relative Importance in Networks, White & Smyth, 2003

  # if v is not an igraph-class object, convert
	if (class(v) != "igraph.vs") v <- AsiGraph(v, g)

  # identify all of the connected clusters of vertices contained within g
  c        <- clusters(g)
  clustN   <- c$no 
  clustMem <- c$membership
  vertN    <- vcount(g)

  # calculate the mean-first-passage-time for each of the connected clusters of vertices
  M        <- matrix(rep(Inf, vcount(g)^2), vcount(g), vcount(g))
  adj	   <- get.adjacency(g, attr=edge.attr, sparse=FALSE)
  max.dist <- max(adj)
  min.dist <- min(adj[adj > 0])
  
  for (i in 1:clustN) {
    locs   <- which(clustMem==(i))
    M[locs, locs]	<- suppressWarnings(MFPTfct(delete.vertices(g, (1:vertN)[-locs]), edge.attr=edge.attr, max.dist=max.dist, min.dist=min.dist))
  }
  
  # since the distance from vertex A to B may not be the same as the distance from vertex B to A, take the mean of the 2 as the final distance between the two vertices if specified
	if (average.distances) {
    lt    <- lower.tri(M)
    tmp   <- apply(cbind(M[lt], t(M)[lt]), 1, mean)
    M[lt] <- tmp
    M     <- t(M)
    M[lt]	<- tmp
	}
  
  # ensure the diagonal of the distance matrix equals 0
	diag(M) <- 0

	M <- as.matrix(M[v, ])
	rownames(M) <- v
  M
}
