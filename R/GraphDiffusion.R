GraphDiffusion <- function(
    g,
    v=V(g), 
    edge.attr=NULL,
    beta=1, 
    correct.neg=TRUE
) { 
    # calculates the diffusion distance using the diffusion kernel method
    # Diffusion Kernels on Graphs and Other Discrete Structures, Kondor and Lafferty, 2002
    # a list containing the numeric graph kernel matrix and the named numeric distance matrix is returned
    
    if (class(v) != "igraph.vs") v <- AsiGraph(v, g) # if v is not an igraph-class object, convert
    res.dimnames <- list(v$name, V(g)$name) 
    if (vcount(g) == 1) return(list(kernel=matrix(0, 1, 1, dimnames=res.dimnames), dist=matrix(0, 1, 1, dimnames=res.dimnames))) # if only a single vertex is contained within the graph, return a single cell matrix

    if (!is.null(edge.attr)) {
        # if an attribute name containing the edge distances is supplied, convert to weights through inversion
        edge.distances <- get.edge.attribute(g, edge.attr)
        edge.weights <- max(edge.distances) - edge.distances + min(edge.distances)
        g <- set.edge.attribute(g, name=edge.attr, value=edge.weights)
    }
    
    # obtain the sparse unamed adjacency matrix
    adj <- get.adjacency(g, attr=edge.attr, names=F, sparse=T, type="both")
    
    # compute the diffusion kernal
    H <- adj - Diagonal(x=apply(adj, 1, sum))
    x <- eigen(H, symmetric=T)
    K  <- x$vectors %*% diag(exp(beta * x$values)) %*% t(x$vectors)
    D2 <- outer(diag(K), diag(K), "+") - 2*K
    
    # correct negative distances
    if (any(D2 < 0) && correct.neg){
        warning("negative distances set to zero")
        D2[D2 < 0] <- 0
    }
    D <- sqrt(D2)
    
    # select only the rows corresponding to v
    D <- as.matrix(D[v, ])
    K <- as.matrix(K[v, ])
    dimnames(D) <- res.dimnames
    
    list(kernel=K, dist=D)
}

