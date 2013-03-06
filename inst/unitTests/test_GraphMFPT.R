test_GraphMFPT <- function() {
  # setup main graph
  g <- graph.empty(7, directed=FALSE)
  g <- add.edges(g, c(1, 2, 1, 3, 2, 4, 3, 4, 1, 5, 5, 6))
  D1 <- GraphMFPT(g)
  
  # test that GraphMFPT returns correct relative distance on a present graph
  checkEquals(all(D1[1, 1] == 0), TRUE)
  checkEquals(all(D1[1, 4] < D1[1, 6]), TRUE)
  checkEquals(all(D1[4, 2] - D1[4, 3] < 1e-10), TRUE)
  checkEquals(all(D1[4, 5] < D1[4, 6]), TRUE)

  # test that GraphMFPT assigns unconnected vertices infinite vertex distances
  checkEquals(all(is.finite(D1[1, 7]) == FALSE), TRUE)
  
  # test that GraphMFPT correctly takes into account edge distances if specified
  g <- set.edge.attribute(g, "edge.attr", value=c(0.1, 0.9, 0.1, 0.9, 0.5, 0.5))
  D2 <- GraphMFPT(g, edge.attr="edge.attr")
  checkEquals(all(D2[1, 2] < D2[1, 3]), TRUE)
  checkEquals(all(D2[4, 2] < D2[4, 3]), TRUE)
  
  # test that GraphMFPT correctly averages reciprical distances when and when not specified
  D3 <- GraphMFPT(g, average.distances=TRUE)
  checkEquals(as.numeric(D3[5, 6]) == as.numeric(D3[6, 5]), TRUE)
  D4 <- GraphMFPT(g, average.distances=FALSE)
  checkEquals(as.numeric(D4[5, 6]) == as.numeric(D4[6, 5]), FALSE)
  
  # test that GraphMFPT returns a 1x1 matrix if the graph contains only 1 vertex
  g <- graph.empty(1, directed=FALSE)
  D5 <- GraphMFPT(g)
  checkEquals(dim(D5), c(1,1))
 
}
