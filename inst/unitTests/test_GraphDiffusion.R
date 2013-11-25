test_GraphDiffusion <- function() {  
    # setup main graph
    g <- graph.empty(5, directed=F)
    g <- add.edges(g, c(2,1, 2,3, 2,4))
    
    # check that both distance and kernal returned
    res <- GraphDiffusion(g)
    checkIdentical(names(res), c("kernel", "dist"))
    D1 <- res$dist
    
    # test that function does not return infinite distances even if vertices are unconnected
    checkTrue(all(is.finite(D1)))
    
    # test that function returns correct relative distances
    checkTrue(all(diag(D1) == 0))
    checkTrue(D1[1,2] < D1[1,3])
    checkTrue(D1[1,3] - D1[1,4] < 10e-10)
    checkTrue(D1[1,4] < D1[1,5])
    
    # test that function correctly takes into account edge distances if specified
    edge.distances <- c(1, 2, 3)
    edge.attr <- "distances"
    ge <- set.edge.attribute(g, name=edge.attr, value=edge.distances)
    D2 <- GraphDiffusion(ge, edge.attr=edge.attr)$dist
    checkTrue(D2[2, 1] < D2[2, 3])
    checkTrue(D2[2, 3] < D2[2, 4])
    
    # test that the distance between two unconnected nodes is smaller when one of the nodes itself is well connected
    checkTrue(D2[5,2] < D2[5,4])
    
    # test that function returns a 1x1 matrix if the graph contains only 1 vertex
    g1 <- graph.empty(1, directed=F)
    D3 <- GraphDiffusion(g1)$dist
    checkIdentical(as.integer(dim(D3)), as.integer(c(1, 1)))
    checkEquals(D3[1, 1], 0, checkNames=F)
}
