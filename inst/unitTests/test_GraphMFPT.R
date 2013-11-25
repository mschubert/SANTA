test_GraphMFPT <- function() {
    # setup main graph
    g <- graph.empty(5, directed=F)
    g <- add.edges(g, c(2,1, 2,3, 2,4))
    
    # test that function results correct distances
    D1 <- GraphMFPT(g, average.distances=F)
    checkTrue(all(diag(D1) == 0))
    checkEquals(D1[1, 2], 1, checkNames=F)
    checkEquals(D1[3, 2], 1, checkNames=F)
    checkEquals(D1[4, 2], 1, checkNames=F)
    
    # test that function assigns unconnected vertices infinite vertex distances
    checkTrue(all(!is.finite(D1[1:4, 5])))
    checkTrue(all(!is.finite(D1[5, 1:4])))
    
    # test that function correctly averages reciprical distances only when specified
    D2 <- GraphMFPT(g, average.distances=T)
    checkIdentical(D2 + D2, D1 + t(D1))
    checkTrue(!identical(D1, D2))
    
    # test that function correctly takes into account edge distances if specified
    edge.distances <- c(1, 2, 3)
    edge.attr <- "distances"
    ge <- set.edge.attribute(g, name=edge.attr, value=edge.distances)
    D3 <- GraphMFPT(ge, edge.attr=edge.attr, average.distances=T)
    checkTrue(D3[2, 1] < D3[2, 3])
    checkTrue(D3[2, 3] < D3[2, 4])
    
    # test that function returns a 1x1 matrix if the graph contains only 1 vertex
    g1 <- graph.empty(1, directed=F)
    D4 <- GraphMFPT(g1)
    checkIdentical(as.integer(dim(D4)), as.integer(c(1, 1)))
    checkEquals(D4[1, 1], 0, checkNames=F)
    
    # test that the function works on a graph with vertex names
    vertex.names <- paste("vertex", 1:vcount(g), sep="")
    gn <- set.vertex.attribute(g, "name", value=vertex.names)
    D5 <- GraphMFPT(gn)
}

