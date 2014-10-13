GraphMFPT <- function(
    g,
    v=V(g),
    edge.attr.weight=NULL,
    average.distances=TRUE
) {
    # calculate the mean first-passage time (MFPT) between every vertex pair
    # Algorithms for Estimating Relative Importance in Networks, White & Smyth, 2003
    # this function will return Inf distances if the graph is not connected
    # any edge attribute specified should be the weight of the edge (higher weights -> more significant) not the distance
    
    if (class(v) != "igraph.vs") v <- AsiGraph(v, g) # if v is not an igraph-class object, convert
        
    # obtain the sparse unamed adjacency matrix
    adj <- get.adjacency(g, attr=edge.attr.weight, names=F, sparse=T, type="both")
    
    # compute the distance matrix for each of the connected clusters of vertices
    D <- array(Inf, dim=rep(vcount(g), 2))
    c <- clusters(g) 
    for (i in 1:c$no) {
        indices <- which(c$membership == i)
        D[indices, indices] <- MFPTfct(adj[indices, indices])
    }
     
    if (average.distances) D <- (D + t(D)) / 2 # if required, take the average of the reciprical distances
    diag(D) <- 0 # the distance between vertex A and A should always be 0
    
    # return rows for each vertex in v
    D <- as.matrix(D[v, ])
    dimnames(D) <- list(v$name, V(g)$name) 
    
    D
}



MFPTfct <- function(
    adj
) {
    # calculate the mean first-passage time (MFPT) for a fully connected graph from the adjacency matrix
    # note: this function is unable to deal with graphs that are not fully connected
    
    # if the adjacency matrix contains only a single gene, return a 1x1 matrix containing 0
    if (is.null(dim(adj))) return(matrix(0, 1, 1))
    ngenes <- nrow(adj)
    
    A <- adj / apply(adj, 1, sum)  # A: the transition probability matrix (there is always movement)
    I <- Diagonal(x=rep(as.integer(1), ngenes))
    pi <- as.numeric(rep(1/ngenes, ngenes) %*% solve(I - A + 1 / ngenes)) # pi: the stationary distribution of the transition matrix
    Z <- solve(t(t(I - A) - pi))
    as.matrix(t(t(I - Z) + Z[cbind(1:nrow(Z), 1:nrow(Z))]) %*% (I * (1 / pi))) # M: the mean first passage matrix
}
