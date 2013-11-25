GraphDiffusion <- function(
    g,
    v=V(g), 
    edge.attr=NULL,
    beta=1, 
    correct.neg=TRUE
) { 
    # calculates the diffusion distance using the diffusion kernel method
    # Diffusion Kernels on Graphs and Other Discrete Structures, Kondor and Lafferty, 2002
    
    if (class(v) != "igraph.vs") v <- AsiGraph(v, g) # if v is not an igraph-class object, convert
    if (vcount(g) == 1) return(list(kernel=matrix(0, 1, 1), dist=matrix(0, 1, 1))) # if only a single vertex is contained within the graph, return a single cell matrix
    
    # modify the network distances, so that smaller distances produce greater chances of movement along the edge
    if (!is.null(edge.attr)) {
        edge.weights <- get.edge.attribute(g, name=edge.attr)
        g <- set.edge.attribute(g, name="edge.distances.inverted", value=max(edge.weights) - edge.weights + min(edge.weights))
        adj  <- get.adjacency(g, attr="edge.distances.inverted", sparse=T, type="both")
    } else {
        adj  <- get.adjacency(g, attr=NULL, sparse=T, type="both")
    }
    
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
    D <- as.matrix(D[v, ])
    K <- as.matrix(K[v, ])
    rownames(D) <- v
    
    list(kernel=K, dist=D)
}
