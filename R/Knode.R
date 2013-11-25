Knode <- function(
    g, 
    dist.method = "shortest.paths", 
    vertex.attr = "pheno", 
    edge.attr = "distance",
    correct.factor = 1,
    nsteps = 1000,
    only.Knode = TRUE,
    vertex.weight = TRUE,
    cluster.id = FALSE,
    vertex.degree = TRUE,
    boncich.power = FALSE,
    markov.centr = FALSE,
    B = NULL, 
    verbose = TRUE
) {  
    # rank the genes exhibiting the phenotype by Knode AUK score. Also produce and list a number of other score for each vertex
    
    # check that vertex weights and edge distances are present and suitable. Convert if neccessary
    g <- CheckAttributes(g, vertex.attr, edge.attr)
    
    # if vertices do not have a name, add some
    if (is.null(get.vertex.attribute(g, "name"))) g <- set.vertex.attribute(g, "name", value=as.character(1:vcount(g)))
    
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
    
    if (!only.Knode) {
        # if required, compute each of the additional graph stats
        stats <- list()
        if(cluster.id) stats[["cluster.id"]]          <- get.vertex.attribute(g, "hits.cluster")
        if(vertex.degree) stats[["vertex.degree"]]    <- degree(graph=g)
        stats[["vertex.betweenness"]]                 <- betweenness(graph=g, directed=F)
        if(boncich.power) stats[["bonpow"]]           <- bonpow(graph=g)
        stats[["constraint"]]                         <- constraint(graph=g, weights=get.edge.attribute(g, edge.attr))
        stats[["evcent"]]                             <- evcent(graph=g, weights=get.edge.attribute(g, edge.attr))$vector
        stats[["pagerank"]]                           <- page.rank(graph=g, directed=F, weights=get.edge.attribute(g, edge.attr))$vector
        stats[["authority"]]                          <- authority.score(graph=g)$vector
        stats[["hub"]]                                <- hub.score(graph=g)$vector
        if(markov.centr) stats[["markov.centrality"]] <- MarkovCentrality(g=g, edge.attr=edge.attr)
    } 
    
    # the results for each vertex attrbitute are saved as a data frame. Theses data frames are stored in a list
    res <- vector("list", length(vertex.attr))
    names(res)	<- vertex.attr
    
    # compute the Knode scores for each of the vertex attributes in vertex.attr
    for (attr in vertex.attr) {
        attr.message <- if (length(vertex.attr) == 1) NULL else paste(" (", which(vertex.attr == attr), "/", length(vertex.attr), ")", sep="")
        if (verbose) message("computing the Knode scores of vertices using '", attr, "' as weights", attr.message, "...", appendLF=F)
        
        vertex.weights <- get.vertex.attribute(g, attr)
        nodeAUK        <- Kfct(B, vertex.weights, individual=T)$nodeAUK
        
        if (!only.Knode) {
            res[[attr]] <- if (vertex.weight) cbind(nodeAUK, vertex.weights, data.frame(stats)) else cbind(nodeAUK, data.frame(stats)) # combine the Kfct results with the statistics
            res[[attr]] <-res[[attr]][order(nodeAUK, runif(length(nodeAUK)), decreasing=T), ] # sort the data frame of Knode scores and statistics by Knode score
        } else {
            names(nodeAUK) <- get.vertex.attribute(g, "name")
            res[[attr]] <- data.frame(nodeAUK[order(nodeAUK, runif(length(nodeAUK)), decreasing=T)]) # sort the Knode scores
        }
        
        if (verbose) message(" done")
    }
    
    # if only one vertex attribute is input, don't return a list of lists of results
    if (length(vertex.attr) == 1) res <- res[[1]]
    res
} 
