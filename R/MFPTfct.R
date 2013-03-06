MFPTfct <- function(
  g, 
  edge.attr=NULL,
  max.dist=NULL,
  min.dist=NULL
) {
  # calculate the mean first-passage time (mfpt) for a fully connected graph
  # note that this function is unable to deal with graphs that are not fully connected
	
  # adj: the adjacency matrix
  adj <- get.adjacency(g, attr=edge.attr, sparse=FALSE) 
  
  # invert the adjacency matrix, so that small distances result in a greater chance of movement along the edge
  if (is.null(max.dist)) max.dist <- max(adj)
  if (is.null(min.dist)) min.dist <- min(adj[adj > 0])
  adj[adj > 0] <- max.dist - adj[adj > 0] + min.dist

  # A: the transition probability matrix (there is always movement)
  A	  <- adj / apply(adj, 1, sum) 
  if (all(is.nan(A))) A[1, 1] <- 1

  # pi: the stationary distribution of the transition matrix
  I	  <- diag(vcount(g))
  U	  <- matrix(1 / nrow(A), nrow(A), nrow(A)) 
  pi  <- U[1,] %*% solve(I - A + U)

  # Z: the fundamental matrix
  e	  <- rep(1, nrow(A))
  Z	  <- solve(I - A - e %*% pi)

  # M: the mean first passage matrix
  Zdg	<- I * Z
  E	  <- matrix(1, vcount(g), vcount(g))
  M	  <- (I - Z + E %*% Zdg) %*% (I * as.vector(1 / pi)) 
  M
}
