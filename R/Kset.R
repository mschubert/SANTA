Kset <- function(
    g, 
    nperm=100, 
    dist.method=c("shortest.paths", "diffusion", "mfpt"), 
    seed.vertex.attr="seed",
    query.vertex.attr="query",
    version=c("b", "a"),
    fix.vertex.degree=TRUE,
    edge.attr=NULL,
    correct.factor=1,
    nsteps=1000, 
    prob=c(0, 0.05, 0.5, 0.95, 1),
    parallel=NULL,
    B=NULL,
    verbose=TRUE
){
    # calculate the Kset function for a graph, along with permutations if required

    # setup
    dist.method <- match.arg(dist.method)
    version <- match.arg(version)
    g <- CheckAttributes(g, seed.vertex.attr=seed.vertex.attr, query.vertex.attr=query.vertex.attr, edge.attr=edge.attr, verbose=verbose)
    degrees <- if (fix.vertex.degree) degree(g) else rep(1, vcount(g))
    
    # the results are saved in a list with an entry for each vertex.attr. Done even if only 1 vertex.attr supplied
    tmp        <- vector("list", 6)
    names(tmp) <- c("K.obs", "AUK.obs", "K.perm", "AUK.perm", "K.quan", "pval")
    class(tmp) <- "Kset"
    res        <- rep(list(tmp), length(query.vertex.attr)) 
    names(res) <- query.vertex.attr
    nvertices  <- as.integer(vcount(g))
    
    if (is.null(B)) {
        # compute B and D
        D <- DistGraph(g, edge.attr=edge.attr, dist.method=dist.method, correct.inf=T, correct.factor=correct.factor, verbose=verbose) 
        B <- BinGraph(D, nsteps=nsteps, verbose=verbose) 
        rm(D)
    } else {
        # check that B is of the correct dimensions
        if (!identical(dim(B), rep(nvertices, 2))) stop("B is not of the correct dimensions")
    }
    
    # convert B to a vector and compute the maximum of B
    Bv <- as.integer(as.vector(B))
    maxB <- as.integer(max(Bv))
    if (version == "a") rm(B)
    
    # if parallel computing is to be used, set up a cluster
    if (!is.null(parallel)) {
        if (verbose) message("setting up cluster of size ", parallel, "...", appendLF=F)
        cl <- makeCluster(parallel, type="SOCK", verbose=F)		
        clusterExport(cl, list("Bv", "nvertices", "maxB"), envir=environment())
        clusterEvalQ(cl, library(SANTA))
        message(" done")
    }
    
    # extract seed vertex weights 
    seed.weights <- get.vertex.attribute(g, seed.vertex.attr)
    seed.weights[is.na(seed.weights)] <- 0
    seed.weights <- as.double(seed.weights)
    
    # run the function for each vertex attribute in query.vertex.attr
    for (query.attr in query.vertex.attr) {
        attr.message <- if (length(query.vertex.attr) == 1) NULL else paste("(", which(query.vertex.attr == query.attr), "/", length(query.vertex.attr), ") ", sep="")
        if (verbose) message("measuring the closeness of ", seed.vertex.attr, " and ", query.attr, " weights ", attr.message, "using ", nperm, " permutations... ", appendLF=F)
        
        # extract query vertex weights 
        query.weights <- get.vertex.attribute(g, query.attr)
        query.weights.is.na <- is.na(query.weights) # don't permute the NA values
        query.weights[query.weights.is.na] <- 0
        query.weights <- as.double(query.weights)
        
        # determine whether it is more efficient to use the .c code that removes rows (manyzeros) or not (fewzeros)
        c.function.to.use <- if (sum(query.weights == 0) < length(query.weights) / 2) "computenetK_fewzeros" else "computenetK_manyzeros"

        # calculate the observed netK and netAUK
        res[[query.attr]]$K.obs <- switch(version,
            "a"=.Call(c.function.to.use, Bv, seed.weights, query.weights, nvertices, maxB),
            "b"=compute.Kset.vB(B, seed.weights, query.weights, maxB)
        )
        res[[query.attr]]$AUK.obs <- sum(res[[query.attr]]$K.obs) / length(res[[query.attr]]$K.obs)
       
        # if specified, run the Knet function on permutations of the graph
        if (!is.null(nperm) & nperm > 0) {
            if (is.null(parallel)) {
                # run permutations without parallel computing
                res[[query.attr]]$K.perm <- switch(version,
                    "a"=sapply(1:nperm, function(i) .Call(c.function.to.use, Bv, seed.weights, Shuffle(query.weights, ignore=query.weights.is.na), nvertices, maxB)),               
                    "b"=sapply(1:nperm, function(i) compute.Kset.vB(B, seed.weights, Shuffle(query.weights, degrees=degrees, ignore=query.weights.is.na), maxB))     
                )
            } else {
                # run permutations with parallel computing
                clusterExport(cl, "vertex.weights", envir=environment())
                res[[query.attr]]$K.perm <- switch(version,
                    "a"=parSapply(cl, 1:nperm, function(i) .Call(c.function.to.use, Bv, seed.weights, Shuffle(query.weights, ignore=query.weights.is.na), nvertices, maxB)),                             
                    "b"=parSapply(cl, 1:nperm, function(i) compute.Kset.vB(B, seed.weights, Shuffle(query.weights, degrees=degrees, ignore=query.weights.is.na), maxB))                       
                )
            }
            
            # calculate the quantiles, AUK and p-values (through the z-score) for the permutations
            res[[query.attr]]$K.quan <- apply(res[[query.attr]]$K.perm, 1, function(x) quantile(x, prob=prob))
            res[[query.attr]]$AUK.perm <- apply(res[[query.attr]]$K.perm, 2, function(x) sum(x) / length(x))
            res[[query.attr]]$pval <- pnorm((res[[query.attr]]$AUK.obs - mean(res[[query.attr]]$AUK.perm)) / sd(res[[query.attr]]$AUK.perm), lower.tail=F) 
            if (is.na(res[[query.attr]]$pval)) res[[query.attr]]$pval <- 1
        } else {
            # if no permutations are run, permutation-related statistics are returns equal to NA
            res[[query.attr]][c("K.perm", "AUK.perm", "K.quan", "pval")] <- NA
        }
         
        if (verbose) message(" done")
    }
    
    # cleanup and output
    if (!is.null(parallel)) stopCluster(cl)	# close the cluster
    if (length(query.vertex.attr) == 1) res <- res[[1]] # if only one vertex attribute is input, don't return a list of lists of results
    res
}

compute.Kset.vB <- function(
    B,
    seed.weights,
    query.weights,
    maxB
) {
    # compute Kobs version B
    if (!all(seed.weights %in% c(0, 1))) stop("when Kset version b is being used, all seed weights should be 0 or 1")
    if (!all(query.weights %in% c(0, 1))) stop("when Kset version b is being used, all seed weights should be 0 or 1")
    seed.inds <- which(seed.weights != 0)
    query.inds <- which(query.weights != 0)
    sapply(1:maxB, function(s) sum(apply(B[seed.inds, query.inds] <= s, 2, function(col) !all(!col))), simplify=T) / sum(query.weights)
}
