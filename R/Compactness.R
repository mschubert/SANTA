Compactness <- function(
    g, 
    nperm=100, 
    dist.method=c("shortest.paths", "diffusion", "mfpt"), 
    vertex.attr="pheno",
    edge.attr="distance", # an attribute containing the edge distance
    correct.factor=1,
    D=NULL,
    verbose=T 
) {
    # use the compactness method to measure the strength of association between a gene set and a network
    # this method is based upon the PathExpand method by Glaab et al. 
    
    # check that vertex weights and edge distances are present and suitable. Convert if neccessary
    dist.method <- match.arg(dist.method)
    g <- CheckAttributes(g, vertex.attr, edge.attr, verbose) 
   
    # produce list to hold results
    element.names <- c("score.obs", "score.perm", "pval") 
    res <- sapply(vertex.attr, function(attr) sapply(element.names, function(i) NA, simplify=F), simplify=F)
    
    # compute D if not input
    if (is.null(D)) {
        D <- DistGraph(g, edge.attr=edge.attr, dist.method=dist.method, correct.inf=T, correct.factor=correct.factor, verbose=verbose) 
    } else {
        if (nrow(D) != vcount(g) | ncol(D) != vcount(g)) stop("D does not have the correct dimensions")
    }
    
    # compute the compactness scores for each vertex attribute
    for (attr in vertex.attr) {
        # compute the observed compactness score
        
        # check the weights
        vertex.weights <- get.vertex.attribute(g, attr) 
        if (!all(vertex.weights %in% c(0, 1))) stop("vertex weights not binary")
        hits <- which(vertex.weights == 1)
        res[[attr]][["score.obs"]] <-  mean(D[hits, hits][lower.tri(D[hits, hits], diag=F)])
        
        if (nperm > 0) {
            # permute the hits and compute the permuted compactness scores
            res[[attr]][["score.perm"]] <- rep(0, nperm)
            for (i in 1:nperm) {
                hits.perm <- sample(vcount(g), length(hits))
                res[[attr]][["score.perm"]][i] <- mean(D[hits.perm, hits.perm][lower.tri(D[hits.perm, hits.perm], diag=F)])
            }
            
            # use the z-test to produce a p-value
            res[[attr]][["pval"]] <-  pnorm((res[[attr]][["score.obs"]] - mean(res[[attr]][["score.perm"]])) / sd(res[[attr]][["score.perm"]]), lower.tail=T) 
        }
    }
    
    # output
    if (length(vertex.attr) == 1) res <- res[[1]] # if only a single vertex attr input, don't return list
    res   
}
