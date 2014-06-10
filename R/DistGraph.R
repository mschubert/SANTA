DistGraph <- function(
    g, 
    v=V(g), 
    edge.attr=NULL,
    dist.method=c("shortest.paths", "diffusion", "mfpt"), 
    correct.inf=TRUE,
    correct.factor=1,
    verbose=TRUE
) {
    # compute the distance matrix for a graph using a specified method
    # the attribute specified by edge.attr should represent edge distances
    
    method <- match.arg(dist.method)
    if (verbose) message("computing graph distance matrix... ", appendLF=F)
    if (class(v) != "igraph.vs") v <- AsiGraph(v, g) # if v is not an igraph object, convect
    
    # convert the distances to weights
    if (is.null(edge.attr)) {
        distances <- rep(1, ecount(g)) # if NULL, then sometimes and error is produced by the shortest.paths function
        edge.attr.weight <- NULL
    } else {
        distances <- get.edge.attribute(g, edge.attr)
        weights <- max(distances) - distances + min(distances)
        g <- set.edge.attribute(g, name="weights_converted_from_distances", value=weights)
        edge.attr.weight <- "weights_converted_from_distances"
    }
    
    # the distance matrix D contains a row for each vertex in v and a column for each vertex in g
    D <- switch(method,
                shortest.paths = shortest.paths(g, v=v, weights=distances),			
                diffusion      = GraphDiffusion(g, v=v, edge.attr.weight=edge.attr.weight, correct.neg=T),
                mfpt           = GraphMFPT(g, v=v, edge.attr.weight=edge.attr.weight, average.distances=T)
    )
    if (sum(is.na(D)) > 0) warning("NA values returned in D")
    
    # if there are unconnected vertex pairs, the function changes the distance to twice the maximum finite distance.
    if (correct.inf) D[!is.finite(D)] <- correct.factor * max(D[is.finite(D)], na.rm=T)
    
    # change column names and row names to the names of the specified vertices
    if (is.null(rownames(D))) rownames(D) <- v$name
    if (is.null(colnames(D))) colnames(D) <- V(g)$name
    if (verbose) message("done")
    D
}
