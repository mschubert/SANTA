test_SpreadHits <- function() {
  # setup
  nvertices <- 100
  g <- CreateGraph(n=nvertices, gen.vertex.weights=FALSE)
  
  # test that SpreadHits adds the start vertex if specified
  for (i in 1:20) {
    start.vertex <- ceiling(runif(1)*nvertices) 
    g1 <- SpreadHits(g, start.vertex=start.vertex)
    checkEquals(get.vertex.attribute(g1, "pheno")[start.vertex], 1)
  }
  
  # test that SpreadHits add vertex attributes "hits", "pheno" and "color"
  checkEquals(list.vertex.attributes(g1), c("hits", "pheno", "color"))
  
  # test that SpreadHits adds the correct number of hits, with and without start vertex specified
  nhits <- 34
  start.vertex <- 20
  g2 <- SpreadHits(g, h=nhits)
  g3 <- SpreadHits(g, h=nhits, start.vertex=start.vertex)
  checkEquals(sum(get.vertex.attribute(g2, "pheno")), nhits)
  checkEquals(sum(get.vertex.attribute(g3, "pheno")), nhits)
  
  # test that SpreadHits produces an error when g is not an igraph object
  checkException(SpreadHits(6), silent=TRUE)
  
  # test that SpreadHits produces an error when lambda is negative
  checkException(SpreadHits(g, lambda=-1), silent=TRUE)
}
