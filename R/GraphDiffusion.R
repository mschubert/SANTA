GraphDiffusion <- function(
    g,
    v=V(g), 
    edge.attr.weight=NULL,
    beta=1, 
    correct.neg=TRUE
) { 
    # calculates the diffusion distance using the diffusion kernel method
    # Diffusion Kernels on Graphs and Other Discrete Structures, Kondor and Lafferty, 2002
    # the named numeric distance matrix is returned
    # any edge attribute specified should be the weight of the edge (higher weights -> more significant) not the distance
    
    if (class(v) != "igraph.vs") v <- AsiGraph(v, g) # if v is not an igraph-class object, convert
    if (vcount(g) == 1) return(matrix(0, 1, 1, dimnames=list(v$name, V(g)$name))) # if only a single vertex is contained within the graph, return a single cell matrix
    
    # obtain the sparse unamed adjacency matrix
    adj <- get.adjacency(g, attr=edge.attr.weight, names=F, sparse=T, type="both")
    
    # compute the distance matrix for each of the connected clusters of vertices
    D <- array(Inf, dim=rep(vcount(g), 2))
    c <- clusters(g) 
    for (i in 1:c$no) {
        indices <- which(c$membership == i)
        D[indices, indices] <- Diffusionfct(adj[indices, indices], beta, correct.neg)
    }
    diag(D) <- 0 # ensure that the diagonal is 0 
    
    # return rows for each vertex in v
    D <- as.matrix(D[v, ])
    dimnames(D) <- list(V(g)$name[as.numeric(v)], V(g)$name) 
    
    D
}



Diffusionfct <- function(
    adj,
    beta,
    correct.neg
) {
    # calculates the diffusion-kernel based distance matrix for an adjacency matrix
    
    # if the adjacency matrix contains only a single gene, return a 1x1 matrix containing 0
    if (is.null(dim(adj))) return(matrix(0, 1, 1))
    
    # compute the diffusion kernal
    H <- adj - Diagonal(x=apply(adj, 1, sum))
    x <- eigen(H, symmetric=T)
    K  <- x$vectors %*% diag(exp(beta * x$values)) %*% t(x$vectors)
    Dsub <- outer(diag(K), diag(K), "+") - 2*K
    
    # correct negative distances
    if (any(Dsub < 0) && correct.neg){
        warning("negative distances set to zero")
        Dsub[Dsub < 0] <- 0
    }
    
    sqrt(Dsub)
}
