DistGraph <- function(
    g, 
    v = V(g), 
    edge.attr = NULL,
    dist.method = c("shortest.paths", "diffusion", "mfpt"), 
    correct.inf = TRUE,
    correct.factor = 1,
    verbose = TRUE
) {
    # compute the distance matrix for a graph using a specified method
    
    if (verbose) message("computing graph distance matrix... ", appendLF=F)
    if (class(v) != "igraph.vs") v <- AsiGraph(v, g) # if v is not an igraph object, convect
    distances <- if (is.null(edge.attr)) NULL else get.edge.attribute(g, edge.attr)
    method <- match.arg(dist.method)
    
    # the distance matrix D contains a row for each vertex in v and a column for each vertex in g
    D <- switch(method,
                shortest.paths = shortest.paths(g, v=v, weights=distances),			
                diffusion      = GraphDiffusion(g, v=v, edge.attr=edge.attr, correct.neg=T)$dist,
                mfpt           = GraphMFPT(g, v=v, edge.attr=edge.attr, average.distances=T)
    )
    
    # if there are unconnected vertex pairs, the function changes the distance to twice the maximum finite distance.
    if (correct.inf) D[!is.finite(D)] <- correct.factor * max(D[is.finite(D)])
    
    # change column names and row names to the names of the specified vertices
    if (is.null(rownames(D))) rownames(D) <- v$name
    if (is.null(colnames(D))) colnames(D) <- V(g)$name
    if (verbose) message("done")
    D
}
