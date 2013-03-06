test_MarkovCentrality <- function() {
  # setup
  g <- graph.empty(7, directed=FALSE)
  g <- add.edges(g, c(1, 2, 1, 3, 1, 4, 1, 5, 5, 6))
  M <- MarkovCentrality(g)
  
  # test that MarkovCentrality returns correct relative centrality measures on a preset graph
  checkEquals(M[1], max(M))
  checkEquals(M[7], min(M))
  checkEquals(M[5] > M[6], TRUE)
  
  # test that MarkovCentrality takes into account edge distances if specified
  g <- set.edge.attribute(g, "edge.attr", value=c(0.1, 0.9, 0.5, 0.5, 0.5))
  M.without.weights <- MarkovCentrality(g)
  M.with.weights <- MarkovCentrality(g, edge.attr="edge.attr")
  checkEquals(M.without.weights[2] - M.without.weights[3] < 1e-10, TRUE)
  checkEquals(M.with.weights[2] > M.with.weights[4], TRUE)
  checkEquals(M.with.weights[4] > M.with.weights[3], TRUE)
}
