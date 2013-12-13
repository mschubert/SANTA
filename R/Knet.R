Knet <- function(
    g, 
    nperm = 100, 
    dist.method = "shortest.paths", 
    vertex.attr = "pheno",
    edge.attr = "distance",
    correct.factor = 1,
    nsteps = 1000, 
    prob = c(0, 0.05, 0.5, 0.95, 1),
    only.pval = F,
    parallel = NULL,
    B = NULL,
    verbose = TRUE
){
    # calculate the Knet function for a graph, along with permutations if required
    
    # check that vertex weights and edge distances are present and suitable. Convert if neccessary
    g <- CheckAttributes(g, vertex.attr, edge.attr)
    
    # the results are saved in a list with an entry for each vertex.attr. Done even if only 1 vertex.attr supplied
    tmp        <- vector("list", 6)
    names(tmp) <- c("K.obs", "AUK.obs", "K.perm", "AUK.perm", "K.quan", "pval")
    class(tmp) <- "Knet"
    res        <- rep(list(tmp), length(vertex.attr)) 
    names(res) <- vertex.attr
    nvertices  <- as.integer(vcount(g))
    
    if (is.null(B)) {
        # compute B and D
        D <- DistGraph(g=g, edge.attr=edge.attr, dist.method=dist.method, correct.inf=T, correct.factor=correct.factor, verbose=verbose) 
        B <- BinGraph(D=D, dist.method=dist.method, nsteps=nsteps, verbose=verbose) 
        rm(D)
    } else {
        # check that B is of the correct dimensions
        if (!identical(dim(B), rep(nvertices, 2))) stop("B is not of the correct dimensions")
    }
    
    # convert B to a vector and compute the maximum of B
    Bv <- as.integer(as.vector(B))
    maxB <- as.integer(max(B))
    rm(B)
    
    # if parallel computing is to be used, set up a cluster
    if (!is.null(parallel)) {
        if (!IsWholeNumber(parallel)) stop("Argument 'parallel' is not a whole number")
        if (verbose) message("setting up cluster of size ", parallel, "...", appendLF=F)
        cl <- makeCluster(parallel, type="SOCK", verbose=F)		
        clusterExport(cl, list("Bv", "nvertices", "maxB"), envir=environment())
        clusterEvalQ(cl, library(SANTA))
        message(" done")
    }
    
    # run the function for each vertex attribute in vertex.attr
    for (attr in vertex.attr) {
        attr.message <- if (length(vertex.attr) == 1) NULL else paste("(", which(vertex.attr == attr), "/", length(vertex.attr), ") ", sep="")
        if (verbose) message("computing the clustering of the '", attr, "' weights ", attr.message, "using ", nperm, " permutations... ", appendLF=F)
        
        # extract vertex weights 
        vertex.weights <- get.vertex.attribute(g, attr)
        vertex.weights[is.na(vertex.weights)] <- 0
        vertex.weights <- as.double(vertex.weights)
        
        # determine whether it is more efficient to use the .c code that removes rows (manyzeros) or not (fewzeros)
        c.function.to.use <- if (sum(vertex.weights == 0) < length(vertex.weights) / 2) "computenetK_fewzeros" else "computenetK_manyzeros"
        
        # calculate the observed netK and netAUK
        res[[attr]]$K.obs <- .Call(c.function.to.use, Bv, vertex.weights, nvertices, maxB)
        res[[attr]]$AUK.obs <- sum(res[[attr]]$K.obs) / length(res[[attr]]$K.obs)
       
        # if specified, run the Knet function on permutations of the graph
        if (!is.null(nperm) & nperm > 0) {
            if (is.null(parallel)) {
                # run permutations without parallel computing
                res[[attr]]$K.perm <- sapply(1:nperm, function(i) .Call(c.function.to.use, Bv, sample(vertex.weights), nvertices, maxB)) 
            } else {
                # run permutations with parallel computing
                clusterExport(cl, "vertex.weights", envir=environment())
                res[[attr]]$K.perm <- parSapply(cl, 1:nperm, function(i) .Call(c.function.to.use, Bv, sample(vertex.weights), nvertices, maxB))
            }
            
            # calculate the quantiles, AUK and p-values (through the z-score) for the permutations
            res[[attr]]$K.quan <- apply(res[[attr]]$K.perm, 1, function(x) quantile(x, prob=prob))
            res[[attr]]$AUK.perm <- apply(res[[attr]]$K.perm, 2, function(x) sum(x) / length(x))
            res[[attr]]$pval <- pnorm((res[[attr]]$AUK.obs - mean(res[[attr]]$AUK.perm)) / sd(res[[attr]]$AUK.perm), lower.tail=F) 
        } else {
            # if no permutations are run, permutation-related statistics are returns equal to NA
            res[[attr]][c("K.perm", "AUK.perm", "K.quan", "pval")] <- NA
        }
        
        # if only the p-value is to be returned, replace the list with the pvalue 
        if (only.pval) res[[attr]] <- res[[attr]]$pval
        
        if (verbose) message(" done")
    }
    
    # cleanup and output
    if (!is.null(parallel)) stopCluster(cl)	# close the cluster
    if (length(vertex.attr) == 1) res <- res[[1]] # if only one vertex attribute is input, don't return a list of lists of results
    res
}

