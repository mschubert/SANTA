Kfct <- function(
  B,
  vertex.weights, 
  individual=TRUE
) {
	# setup
	nvertices <- length(vertex.weights)
	diam <- max(B)
	if (ncol(B) != nvertices | nrow(B) != nvertices) stop("dimensions of B should be nvertices by nvertices")
  
  # compute the adjusted vertex weights
  # missing weights (equal to NA) are given vertex.weights and adjusted.vertex.weights of 0 in order to ensure that they do not effect the results
	is.missing <- is.na(vertex.weights)
	vertex.weights[is.missing] <- 0
	adjusted.vertex.weights <- vertex.weights - mean(vertex.weights[!is.missing])
	adjusted.vertex.weights[is.missing] <- 0
	
  # compute the sum of the vertex weights of each vertex within each bin for each hit
  X <- ComputeWeightsInBins(B, adjusted.vertex.weights, nvertices)
  
  # calculate the various statistics
  nodeS  <- t(apply(X, 1, cumsum)) 
  netS   <- apply(nodeS * vertex.weights, 2, sum) 
  netK   <- 2 * netS / sum(vertex.weights) ^ 2 
	netAUK <- sum(netK) / length(netK) 

	if (individual) {
		# if specified, compute the Knode function and the area under the curve for each vertex
    nodeK   <- 2*nodeS / sum(vertex.weights)
    nodeAUK <- apply(nodeK, 1, function(x) sum(x) / length(x)) 
  } else {
    nodeK   <- NA
    nodeAUK <- NA
  }    
 
  # output    
  if (round(netK[length(netK)], 2)!=0) warning("K-function doesn't end at 0")
  list(netK=netK, netAUK=netAUK, nodeK=nodeK, nodeAUK=nodeAUK)	
}
