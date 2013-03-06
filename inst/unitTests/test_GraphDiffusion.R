test_GraphDiffusion <- function() {  
  # setup main graph
  g <- graph.empty(7, directed=FALSE)
  g <- add.edges(g, c(1, 2, 1, 3, 2, 4, 3, 4, 1, 5, 5, 6))
  D <- GraphDiffusion(g)$dist
  
  # test that GraphDiffusion does not return infinite distances even if vertices are unconnected
  checkEquals(all(is.finite(D)), TRUE)

  # test that GraphDiffusion returns correct relative distances on a preset graph
  checkEquals(all(D[1, 1] == 0), TRUE)
  checkEquals(all(D[1, 4] < D[1, 6]), TRUE)
  checkEquals(all(D[4, 2] - D[4, 3] < 1e-10), TRUE)
  checkEquals(all(D[4, 5] < D[4, 6]), TRUE)
  checkEquals(all(D[6, 7] == max(D)), TRUE)
  
  # test that GraphDiffusion returns both the diffusion kernel and the diffusion distance
  g2 <- erdos.renyi.game(5, 0.5)
  checkEquals(names(GraphDiffusion(g2)), c("kernel", "dist"))
  
  # test that GraphDiffusion returns a 1x1 matrix if the graph contain only 1 vertex
  g3 <- graph.empty(1, directed=FALSE)
  checkEquals(GraphDiffusion(g3)$dist, matrix(0, 1, 1))
  
  # test that GraphDiffusioncorrectly takes into account edge distances if specified
  g <- set.edge.attribute(g, "edge.attr", value=c(0.1, 0.9, 0.1, 0.9, 0.5, 0.5))
  D1 <- GraphDiffusion(g, edge.attr="edge.attr")$dist
  checkEquals(all(D1[1, 2] < D1[1, 3]), TRUE)
  checkEquals(all(D1[4, 2] < D1[4, 3]), TRUE)
}
