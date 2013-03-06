test_Knet <- function() {
  # setup 
  # graph used:
  # 2 - 1 - 4 - 5
  #     | 
  #     3
  
  g <- graph.empty(5, directed=FALSE)
  g <- add.edges(g, c(1,2,1,3,1,4,4,5))
  g <- set.vertex.attribute(g, "vweights1", value=c(1,0,0,1,0))
  g <- set.vertex.attribute(g, "vweights2", value=c(0,1,0,0,1))
  g <- set.edge.attribute(g, "eweights1", value=c(1, 1, 0.01, 1))
  g <- set.edge.attribute(g, "eweights2", value=c(0.01, 0.01, 1, 0.01))
  
  # test that Knet correctly incorperates different vertex weights
  # all edge distances equal 1 as no attribute specified
  suppressMessages(res1 <- Knet(g, nperm=0, vertex.attr="vweights1"))
  suppressMessages(res2 <- Knet(g, nperm=0, vertex.attr="vweights2"))
  checkEquals(res1$AUK.obs, 0.3)
  checkEquals(res2$AUK.obs, 0.1)
  checkEquals(round(res1$nodeAUK, 10), c(0.25, 0.05, 0.05, 0.35, 0.15))
  checkEquals(round(res2$nodeAUK, 10), c(-0.25, 0.05, -0.45, -0.15, 0.15))
  
  # test that Knet correctly incorperates different edge distances when using the different distance measures
  dist.methods <- c("shortest.paths", "diffusion", "mfpt")
  for (dist.method in dist.methods) {
    suppressMessages(res3 <- Knet(g, nperm=0, dist.method=dist.method, vertex.attr="vweights1", edge.attr="eweights1"))
    suppressMessages(res4 <- Knet(g, nperm=0, dist.method=dist.method, vertex.attr="vweights1", edge.attr="eweights2"))
    checkEquals(res3$AUK.obs > res4$AUK.obs, TRUE)
  }

  # test that Knet outputs a list of results containing the correct elements when one vertex attribute is input
  suppressMessages(res5 <- Knet(g, nperm=10, vertex.attr="vweights1", edge.attr="eweights1"))
  checkEquals(names(res5), c("K.obs", "AUK.obs", "K.perm", "AUK.perm", "K.quan", "nodeK", "nodeAUK", "pval"))
  
  # test that Knet outputs a list of lists of results containing the correct elements when multiple vertex attributes are input
  vertex.attributes <- c("vweights1", "vweights2")
  suppressMessages(res6 <- Knet(g, nperm=10, vertex.attr=vertex.attributes, edge.attr="eweights1"))
  checkEquals(length(res6), length(vertex.attributes))
  checkEquals(names(res6[[1]]), c("K.obs", "AUK.obs", "K.perm", "AUK.perm", "K.quan", "nodeK", "nodeAUK", "pval"))

  # test that Knet outputs a list of results, some of which equal NA, when one attribute and no permutations are completed
  suppressMessages(res7 <- Knet(g, nperm=0, vertex.attr="vweights1", edge.attr="eweights1"))
  checkEquals(all(is.na(res7[c("pval", "K.perm", "AUK.perm", "K.quan")])), TRUE)
}
