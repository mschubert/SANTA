Knode <- function(
    g, 
    dist.method=c("shortest.paths", "diffusion", "mfpt"), 
    vertex.attr="pheno", 
    edge.attr=NULL,
    correct.factor=1,
    nsteps=1000,
    B=NULL, 
    verbose=TRUE
) {  
    # rank the genes exhibiting the phenotype by Knode AUK score. Also produce and list a number of other score for each vertex
    
    # setup
    dist.method <- match.arg(dist.method)
    g <- CheckAttributes(g, vertex.attr, edge.attr, verbose) # check that vertex weights and edge distances are present and suitable
    nvertices <- as.integer(vcount(g))
    if (is.null(get.vertex.attribute(g, "name"))) g <- set.vertex.attribute(g, "name", value=as.character(1:nvertices)) # if vertices do not have a name, add some
    
    # compute D and B if required
    if (is.null(B)) {
        D <- DistGraph(g, edge.attr=edge.attr, dist.method=dist.method, correct.inf=T, correct.factor=correct.factor, verbose=verbose) # compute the vertex pair distances (D)
        B <- BinGraph(D, nsteps=nsteps, verbose=verbose) # compute which bin each vertex pair distance falls into (B)
        rm(D)
    } else {
        if (!identical(dim(B), rep(nvertices, 2))) stop("B is not of the correct dimensions")
    }
    
    # convert to B a vector and compute the maximum of B
    Bv <- as.integer(as.vector(B))
    maxB <- as.integer(max(B))
    rm(B)
    
    # the results for each vertex attrbitute are saved as a data frame. Theses data frames are stored in a list
    res <- vector("list", length(vertex.attr))
    names(res)	<- vertex.attr
    
    # compute the Knode scores for each of the vertex attributes in vertex.attr
    for (attr in vertex.attr) {
        if (verbose) message("running on vertex attribute ", which(vertex.attr == attr), "/", length(vertex.attr), "... ", appendLF=F)
        
        # extract the vertex weights, set any missing weights to the mean
        vertex.weights <- get.vertex.attribute(g, attr)
        vertex.weights[is.na(vertex.weights)] <- mean(vertex.weights, na.rm=T)
        vertex.weights <- as.double(vertex.weights)
        
        # compute nodeAUKs
        nodeAUK <- .Call("computenodeAUK", Bv, vertex.weights, nvertices, maxB)
        names(nodeAUK) <- get.vertex.attribute(g, "name")
        res[[attr]] <- nodeAUK[order(nodeAUK, runif(length(nodeAUK)), decreasing=T)] # sort the Knode scores
        
        if (verbose) message("done")
    }
    
    # if only one vertex attribute is input, don't return a list of lists of results
    if (length(vertex.attr) == 1) res <- res[[1]]
    res
} 
