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
    tmp        <- vector("list", 8)
    names(tmp) <- c("K.obs", "AUK.obs", "K.perm", "AUK.perm", "K.quan", "nodeK", "nodeAUK", "pval")
    class(tmp) <- "Knet"
    res        <- rep(list(tmp), length(vertex.attr)) 
    names(res) <- vertex.attr
    
    if (is.null(B)) {
        # compute the vertex pair distances (D)
        if (verbose) message("computing graph distance matrix...", appendLF=F)
        D <- DistGraph(g=g, edge.attr=edge.attr, dist.method=dist.method, correct.inf=T, correct.factor=correct.factor) 
        if (verbose) message(" done")
        
        # compute which bin each vertex pair distance falls into (B)
        if (verbose) message("computing graph distance bins...", appendLF=F)
        B <- BinGraph(D=D, dist.method=dist.method, nsteps=nsteps) 
        if (verbose) message(" done")   
    } else {
        if (!identical(dim(B), rep(vcount(g), 2))) stop("B is not of the correct dimensions")
    }
    
    # if parallel computing is to be used, set up a cluster
    if (!is.null(parallel)) {
        if (!IsWholeNumber(parallel)) stop("Argument 'parallel' is not a whole number")
        if (verbose) message("setting up cluster of size ", parallel, "...", appendLF=F)
        cl <- makeCluster(parallel, type="SOCK", verbose=F)		
        clusterEvalQ(cl, library(SANTA))
        clusterExport(cl, "B", envir=environment())
        message(" done")
    }
    
    # run the function for each vertex attribute in vertex.attr
    for (attr in vertex.attr) {
        attr.message <- if (length(vertex.attr) == 1) NULL else paste("(", which(vertex.attr == attr), "/", length(vertex.attr), ") ", sep="")
        if (verbose) message("computing the clustering of the '", attr, "' weights ", attr.message, "using ", nperm, " permutations...", appendLF=F)
        
        # extract vertex weights 
        vertex.weights <- get.vertex.attribute(g, attr)
        
        # calculate the observed netK and netAUK
        K                   <- Kfct(B, vertex.weights, individual=T) 
        res[[attr]]$K.obs 	<- K$netK
        res[[attr]]$AUK.obs <- K$netAUK
        res[[attr]]$nodeK 	<- K$nodeK
        res[[attr]]$nodeAUK <- K$nodeAUK
        
        # if specified, run the Knet function on permutations of the graph
        if (!is.null(nperm) & nperm>0) {
            if (is.null(parallel)) {
                # run permutations without parallel computing
                res[[attr]]$K.perm 	<- sapply(1:nperm, function(i) Kfct(B, sample(vertex.weights), individual=F)$netK) 
            } else {
                # run permutations with parallel computing
                clusterExport(cl, "vertex.weights", envir=environment())
                res[[attr]]$K.perm 	<- parSapply(cl, 1:nperm, function(i) Kfct(B, vertex.weights=sample(vertex.weights), individual=F)$netK)
            }
            
            # calculate the quantiles, AUK and p-values (through the z-score) for the permutations
            res[[attr]]$K.quan    <- apply(res[[attr]]$K.perm, 1, function(x) quantile(x, prob=prob))
            res[[attr]]$AUK.perm  <- apply(res[[attr]]$K.perm, 2, function(x) sum(x) / length(x))
            res[[attr]]$pval      <- pnorm((res[[attr]]$AUK.obs - mean(res[[attr]]$AUK.perm)) / sd(res[[attr]]$AUK.perm), lower.tail=F) 
        } else {
            # if no permutations are run, permutation-related statistics are returns equal to NA
            res[[attr]]$K.perm    <- NA
            res[[attr]]$AUK.perm  <- NA
            res[[attr]]$K.quan    <- NA
            res[[attr]]$pval      <- NA	
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
