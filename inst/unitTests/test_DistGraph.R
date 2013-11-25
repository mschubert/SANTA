test_DistGraph <- function() {
  # setup
  g <- graph.empty(6, directed=F)
  g <- add.edges(g, c(1, 2, 1, 3, 2, 4, 3, 4, 1, 5, 5, 6))
  g <- set.edge.attribute(g, "weights", value=c(0.1, 0.9, 0.1, 0.9, 0.5, 0.5))
  
  # test that DistGraph correctly returns the correct relative distances using the shortest.paths algorithm
  # without edge weights
  D1 <- DistGraph(g, dist.method="shortest.paths")
  checkEquals(D1[1, 1], 0)
  checkEquals(D1[1, 4], D1[1, 6])
  
  # with edge weights
  D2 <- DistGraph(g, dist.method="shortest.paths", edge.attr="weights")
  checkTrue(all(D2[1, 2] < D2[1, 3]))
  checkTrue(all(D2[4, 2] < D2[4, 3]))
  
  # test that DistGraph correctly returns the correct relative distances using the diffusion algorithm with edge weights
  # without edge weights
  D3 <- DistGraph(g, dist.method="diffusion")
  checkEquals(D3[1, 1], 0)
  checkTrue(D3[1, 4] < D3[1, 6])
  
  # with edge weights
  D4 <- DistGraph(g, dist.method="diffusion", edge.attr="weights")
  checkTrue(all(D4[1, 2] < D4[1, 3]))
  checkTrue(all(D4[4, 2] < D4[4, 3]))
  
  # test that DistGraph correctly returns the correct relative distances using the MFPT algorithm with edge weights
  # without edge weights
  D5 <- DistGraph(g, dist.method="mfpt")
  checkEquals(D5[1, 1], 0)
  checkTrue(D5[1, 4] < D5[1, 6])
  
  # with edge weights
  D6 <- DistGraph(g, dist.method="mfpt", edge.attr="weights")
  checkTrue(all(D6[1, 2] < D6[1, 3]))
  checkTrue(all(D6[4, 2] < D6[4, 3]))
  
  # setup
  g <- graph.empty(7, directed=F)
  g <- add.edges(g, c(1, 2, 1, 3, 2, 4, 3, 4, 1, 5, 5, 6))
  
  # test that DistGraph does and does not correct infinite values when specified
  D7 <- DistGraph(g, dist.method="shortest.paths", correct.inf=F)
  checkTrue(!is.finite(D7[1, 7]))
  D7 <- DistGraph(g, dist.method="shortest.paths", correct.inf=T)
  checkTrue(is.finite(D7[1, 7]))
  
  # test that DistGraph correctly changes infinite values to correct.factor * maximum finite distance measured
  correct.factor <- 1.5
  D8 <- DistGraph(g, dist.method="shortest.paths", correct.inf=T, correct.factor=correct.factor)
  checkEquals(D8[1, 7],correct.factor * (max(D8[D8 != D8[1, 7]])))
  correct.factor <- 4.5
  D9 <- DistGraph(g, dist.method="shortest.paths", correct.inf=T, correct.factor=correct.factor)
  checkEquals(D9[1, 7], correct.factor * (max(D9[D9 != D9[1, 7]])))
  
  # test that DistGraph correctly handles the lack of infinite values on fully connected graphs
  g <- graph.empty(6, directed=F)
  g <- add.edges(g, c(1, 2, 1, 3, 2, 4, 3, 4, 1, 5, 5, 6))
  D10 <- DistGraph(g, dist.method="shortest.paths", correct.inf=T)
  checkTrue(all(is.finite(D10)))
}
