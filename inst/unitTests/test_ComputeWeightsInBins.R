test_ComputeWeightsInBins <- function() {
  # test that ComputeWeightsInBins can handle SMALL matrix B
  # setup 
  # B matrix used: 
  #   1 2 2 2 3
  #   2 1 3 3 4 
  #   2 3 1 3 4
  #   2 3 3 1 2
  #   3 4 4 2 1
  
  B <- matrix(c(1,2,2,2,3,2,1,3,3,4,2,3,1,3,4,2,3,3,1,2,3,4,4,2,1), 5,5)
  vertex.weights <- c(1,0,0,1,0)
  adjusted.vertex.weights <- vertex.weights - mean(vertex.weights)
  nvertices <- length(vertex.weights)
  res <- SANTA:::ComputeWeightsInBins(B, adjusted.vertex.weights, nvertices)
  checkEquals(round(as.vector(res), 10), 
  	c(0.6, -0.4, -0.4, 0.6, -0.4, -0.2, 0.6, 0.6, 0.2, 0.6, -0.4, 0.2, 0.2, -0.8, 0.6, 0.0, -0.4, -0.4, 0.0, -0.8)) 
 
  # test that ComputeWeightsInBins can handle LARGE matrix B
  # setup
  # B matrix used: 
  #  1 1 1 1 ..
  #  2 2 2 2 ..
  #  3 3 3 3 ..
  #  4 4 4 4 ..
  #  : : : :
  nvertices <- 10000
  B <- matrix(c(rep(1:nvertices, nvertices-1), rep(1, nvertices)), nvertices, nvertices)
  vertex.weights <- rep(c(0,1), nvertices/2)
  adjusted.vertex.weights <- vertex.weights - mean(vertex.weights)
  res <- SANTA:::ComputeWeightsInBins(B, adjusted.vertex.weights, nvertices)
  checkEquals(round(as.vector(res[1:3, 1:3]), 10), 
  	c(0, 0.5, 0.5, 0, -0.5, 0, 0, 0, -0.5))
}
