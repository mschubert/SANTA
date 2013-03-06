DistGraph <- function(
  g, 
  v=V(g), 
  edge.attr=NULL,
  dist.method=c("shortest.paths", "diffusion", "mfpt"), 
  correct.inf=TRUE,
  correct.factor=1.0
) {
  # calculates a distance matrix for a graph using a specified method

  # if v is not supplied as an igraph object, convert
  if (class(v) != "igraph.vs") v <- AsiGraph(v, g)
  
  # get edge weights
  if (is.null(edge.attr)) {
    edge.weights <- NULL
  } else {
    edge.weights <- get.edge.attribute(g, edge.attr)
  }

  method <- match.arg(dist.method)
  D <- switch(method,
    # the distance matrices contain a column for each hit and a row for every vertex
    shortest.paths = shortest.paths(g, v=v, weights=edge.weights),			
    diffusion      = GraphDiffusion(g, v=v, edge.attr=edge.attr, correct.factor=correct.factor, correct.neg=TRUE)$dist,
    mfpt           = GraphMFPT(g, v=v, edge.attr=edge.attr, average.distances=TRUE)
	)

  # if there are unconnected vertex pairs, the function changes the distance to twice the maximum finite distance.
  if (correct.inf & dist.method %in% c("shortest.paths", "mfpt")) D[D==Inf] <- correct.factor * max(D[D!=Inf])
  
  # change column names and row names to the names of the specified vertices
  if (nrow(D) == ncol(D)) dimnames(D) <- list(v, v) else rownames(D) <- v
  D
}
