GraphMFPT <- function(
    g,
    v=V(g),
    edge.attr=NULL,
    average.distances=TRUE
) {
    # calculate the mean first-passage time (MFPT) between every vertex pair
    # Algorithms for Estimating Relative Importance in Networks, White & Smyth, 2003
    
    if (class(v) != "igraph.vs") v <- AsiGraph(v, g) # if v is not an igraph-class object, convert
    
    # identify all of the connected clusters of vertices contained within g
    ngenes <- vcount(g)
    c <- clusters(g)
    clustN <- c$no 
    clustMem <- c$membership

    # compute the adjacency matrix
    if (is.null(edge.attr)) {
        # without edge weights
        adj  <- get.adjacency(g, attr=NULL, sparse=T, type="both")
    } else {
        # with edge weights
        # modify the network distances, so that smaller distances produce greater chances of movement along the edge
        edge.weights <- get.edge.attribute(g, name=edge.attr)
        g <- set.edge.attribute(g, name="edge.distances.inverted", value=max(edge.weights) - edge.weights + min(edge.weights))
        adj  <- get.adjacency(g, attr="edge.distances.inverted", sparse=T, type="both")
    }
     
    # calculate the mean-first-passage-time for each of the connected clusters of vertices
    M <- matrix(rep(Inf, ngenes^2), ngenes, ngenes)
    for (i in 1:clustN) {
        locs   <- which(clustMem==(i))
        M[locs, locs] <- MFPTfct(adj[locs, locs])
    }

    # since the distance from vertex A to B may not be the same as the distance from vertex B to A, take the mean of the 2 as the final distance between the two vertices if specified
    if (average.distances) M <- (M + t(M)) / 2

    # ensure the diagonal of the distance matrix equals 0
    diag(M) <- 0
    
    M <- as.matrix(M[v, ])
    rownames(M) <- v
    M
}
