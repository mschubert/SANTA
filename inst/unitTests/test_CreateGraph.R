test_CreateGraph <- function() {
  # test that CreateGraph creates a graph with the correct number of vertices when the erdos renyi algorithm is used
  nvertices <- 1000
  g <- CreateGraph(n=nvertices, type="erdos.renyi")
  checkEquals(vcount(g), nvertices)
  
  # test that CreateGraph creates a graph with the correct number of vertices and edges when the barabasi algorithm is used
  m <- 2
  g <- CreateGraph(n=nvertices, type="barabasi", m=m)
  checkEquals(vcount(g), nvertices)
  checkEquals(ecount(g), m * nvertices - sum(1:m))
  
  # test that CreateGraph creates a graph with the correct number of vertices and edges when the grid algorithm is used
  g <- CreateGraph(n=nvertices, type="grid")
  checkEquals(vcount(g), nvertices)
  if (nvertices==1000) checkEquals(ecount(g), 1936)
  
  # test that CreateGraph does not add vertex or edge attributes if not required 
  g <- CreateGraph(n=nvertices, gen.vertex.weights=FALSE, binary.pheno=FALSE)
  checkEquals(length(list.vertex.attributes(g)), 0)
  checkEquals(length(list.edge.attributes(g)), 0)
  
  # test that CreateGraph generates binary vertex weights when required
  g <- CreateGraph(n=nvertices, gen.vertex.weights=TRUE, binary.pheno=TRUE)
  values <- get.vertex.attribute(g, "pheno")
  checkEquals(all(values == 0 | values == 1), TRUE)
  
  # test that CreateGraph generates continuous vertex weights within the range (0,1) when required
  g <- CreateGraph(n=nvertices, gen.vertex.weights=TRUE, binary.pheno=FALSE)
  values <- get.vertex.attribute(g, "pheno")
  checkEquals(max(values) <= 1, TRUE)
  checkEquals(min(values) >= 0, TRUE)
  
  # test that CreateGraph successfully adds the correct number of links when multiple clusters are added
  nlinks <- 30
  g <- CreateGraph(n=nvertices, type="barabasi", m=m, gen.vertex.weights=TRUE, nclusters=2, nlinks=nlinks)
  checkEquals(ecount(g), 2 * (m * nvertices / 2 - sum(1:m)) + nlinks)
}
