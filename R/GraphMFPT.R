GraphMFPT <- function(
    g,
    v=V(g),
    edge.attr=NULL,
    average.distances=TRUE
) {
    # calculate the mean first-passage time (MFPT) between every vertex pair
    # Algorithms for Estimating Relative Importance in Networks, White & Smyth, 2003
    # this function will return Inf distances if the graph is not connected
    
    if (class(v) != "igraph.vs") v <- AsiGraph(v, g) # if v is not an igraph-class object, convert
    
    if (!is.null(edge.attr)) {
        # if an attribute name containing the edge distances is supplied, convert to weights through inversion
        edge.distances <- get.edge.attribute(g, edge.attr)
        edge.weights <- max(edge.distances) - edge.distances + min(edge.distances)
        g <- set.edge.attribute(g, name=edge.attr, value=edge.weights)
    }
    
    # obtain the sparse unamed adjacency matrix
    adj <- get.adjacency(g, attr=edge.attr, names=F, sparse=T, type="both")
    
    # compute the mean-first-passage-time for each of the connected clusters of vertices
    D <- matrix(Inf, vcount(g), vcount(g))
    c <- clusters(g) 
    for (i in 1:c$no) {
        indices <- which(c$membership == i)
        D[indices, indices] <- MFPTfct(adj[indices, indices])
    }
     
    if (average.distances) D <- (D + t(D)) / 2 # if required, rake the average of the reciprical distances
    diag(D) <- 0 # the distance between vertex A and A should always be 0
    
    # return rows for each vertex in v
    D <- as.matrix(D[v, ])
    dimnames(D) <- list(v$name, V(g)$name) 
    
    D
}
