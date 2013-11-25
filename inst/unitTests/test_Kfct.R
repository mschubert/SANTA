test_Kfct <- function() { 
  # setup 
  # B matrix used
  #   1 2 2 2 3 
  #   2 1 3 3 4
  #   2 3 1 3 4
  #   2 3 3 1 2
  #   3 4 4 2 1
  B <- matrix(c(1,2,2,2,3,2,1,3,3,4,2,3,1,3,4,2,3,3,1,2,3,4,4,2,1), 5,5)
  vertex.weights1 <- c(1, 0, 0, 1, 0)

  # test that Kfct outputs the correct list of statistics when individual=TRUE
  res1 <- Kfct(B, vertex.weights1, individual=T)  
  checkEquals(names(res1), c("netK", "netAUK", "nodeK", "nodeAUK"))
  checkTrue(all(is.na(res1) == F))
  
  # test that Kfct outputs nodeK and nodeAUK=NA when individual=FALSE
  res2 <- Kfct(B, vertex.weights1, individual=F)  
  checkEquals(names(res2), c("netK", "netAUK", "nodeK", "nodeAUK"))
  checkTrue(is.na(res2$nodeK) & is.na(res2$nodeAUK))
  checkEquals(is.na(res2$netK[[1]]) & is.na(res2$netAUK[[1]]), F)
  
  # test that Kfct outputs the correct netK, netAUK, nodeK and node AUK values for a square matrix B
  res3 <- Kfct(B, vertex.weights1, individual=T)  
  checkEquals(round(res3$netK, 10), c(0.6, 0.6, 0.0, 0.0))
  checkEquals(round(res3$netAUK, 10), 0.3)
  checkEquals(dim(res3$nodeK), c(nrow(B), max(B)))
  checkEquals(round(res3$nodeAUK, 10), c(0.25, 0.05, 0.05, 0.35, 0.15))
  
  # test that Kfct correctly incorperates different vertex weights 
  vertex.weights2 <- c(0,1,0,0,1)
  res4 <- Kfct(B, vertex.weights2, individual=T) 
  checkEquals(round(res4$netK, 10), c(0.6, 0.2, -0.4, 0.0))
  checkEquals(round(res4$netAUK, 10), 0.1)
  checkEquals(dim(res4$nodeK), c(nrow(B), max(B)))
  checkEquals(round(res4$nodeAUK, 10), c(-0.25, 0.05, -0.45, -0.15, 0.15))
  
  # test that Kfct produces an error if the dimensions of B is not length(vertex.weights1) x length(vertex.weights1)
  checkException(Kfct(B[1:4,], vertex.weights1, individual=T), silent=T)
  checkException(Kfct(B, vertex.weights1[1:4], individual=T), silent=T)
  
  # test that Kfct correctly handles nodes with weight NA
  # these nodes should be included in the network structure, but their weight should effectively be ignored
  # setup 2 networks, 1 with a node with weight NA (inc) and with the same node deleted
  nvertices <- 8
  vertex.to.remove <- 8
  edges <- c(1,2, 2,3, 3,4, 4,1, 1,5, 2,6, 3,7, 4,8) 
  weights <- runif(nvertices)
  weights[vertex.to.remove] <- NA
  vertex.attr <- "weights"
  edge.attr <- "distance"
  network <- graph.empty(nvertices, directed=F)
  network <- add.edges(network, edges)
  network <- set.vertex.attribute(network, name=vertex.attr, value=weights)
  network <- set.edge.attribute(network, name=edge.attr, value=rep(1, ecount(network)))
  network.inc <- network
  network.exc <- delete.vertices(network, v=vertex.to.remove)
  
  D.inc <- DistGraph(g=network.inc, edge.attr=edge.attr, dist.method="shortest.paths") 
  B.inc <- BinGraph(D=D.inc, dist.method="shortest.paths", nsteps=100) 
  D.exc <- DistGraph(g=network.exc, edge.attr=edge.attr, dist.method="shortest.paths") 
  B.exc <- BinGraph(D=D.exc, dist.method="shortest.paths", nsteps=100) 
  
  res.inc <- Kfct(B.inc, get.vertex.attribute(network.inc, name=vertex.attr))
  res.exc <- Kfct(B.exc, get.vertex.attribute(network.exc, name=vertex.attr))
  
  checkEquals(res.inc$netK, res.exc$netK) 
  checkEquals(res.inc$netAUK, res.exc$netAUK)
  checkEquals(res.inc$nodeK[-vertex.to.remove, ], res.exc$nodeK)
  checkEquals(res.inc$nodeAUK[-vertex.to.remove], res.exc$nodeAUK)
}
