test_MarkovCentrality <- function() {
    # setup main graph
    g <- graph.empty(5, directed=F)
    g <- add.edges(g, c(2,1, 2,3, 2,4))
    
    # test that the function resurts the correct relative centrality measures 
    M <- MarkovCentrality(g)
    checkTrue(M[2] > M[1])
    checkTrue(M[1] - M[3] < 10e-10)
    checkTrue(M[3] - M[4] < 10e-10)
    checkTrue(M[4] > M[5])
    
    # test that MarkovCentrality takes into account edge distances if specified
    edge.distances <- c(1, 2, 3)
    edge.attr <- "distances"
    ge <- set.edge.attribute(g, name=edge.attr, value=edge.distances)
    M <- MarkovCentrality(ge, edge.attr=edge.attr)
    checkTrue(M[2] > M[1])
    checkTrue(M[1] > M[3])
    checkTrue(M[3] > M[4])
    checkTrue(M[4] > M[5])
}
