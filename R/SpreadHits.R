SpreadHits <- function(
    g, 
    h, 
    clusters=1, 
    distance.cutoff=3,
    lambda=1, 
    dist.method=c("shortest.paths", "diffusion", "mfpt"), 
    edge.attr=NULL, 
    hit.color="red", 
    D=NULL, 
    attempts=1000,
    verbose=TRUE
) {
    # spread hits across a graph in multiple clusters
    # this is done by trial by error
    
    # setup
    dist.method <- match.arg(dist.method)
    if (class(g) != "igraph") stop("g is not a graph")
    if (h * clusters > vcount(g)) stop("too few vertices on network to apply this many hits")
    if (lambda < 0) stop("lambda needs to be greater than or equal to 0")
    
    # compute D if not input
    if (is.null(D)) D <- DistGraph(g, edge.attr=edge.attr, dist.method=dist.method, verbose=verbose)
    D.cutoff <- D >= distance.cutoff
      
    # identify the seed vertices
    seed.flag <- F
    attempt.number <- 0
    while (!seed.flag & attempt.number < attempts) {
        # for a certain maximum number of attempt, try to find a suitable number of seed vertices distant appart
        attempt.number <- attempt.number + 1
        potential.seed <- rep(T, vcount(g))
        seeds <- rep(NA, clusters)
        seeds[1] <- sample(vcount(g), 1)
        
        # try to identify seeds that are suitable distant
        while (!all(potential.seed == F) & !all(!is.na(seeds))) {
            seeds.chosen <- sum(!is.na(seeds))
            potential.seed <- if (seeds.chosen == 1) D.cutoff[seeds[seeds.chosen], ] else apply(D.cutoff[seeds[1:seeds.chosen], ], 2, all)
            if (!all(!potential.seed)) seeds[seeds.chosen + 1] <- sample(which(potential.seed), 1)
        }
        
        # if enough hits have been chosen, proceed
        if (all(!is.na(seeds))) seed.flag <- T
    }
    if (!seed.flag) warning("seed vertices not identified")
    
    if (seed.flag) {
        # spread hits over the network, starting from each seed
        hits <- list()
        for (cluster in 1:clusters) {
            # spread the hits across the graph, probability determined by distance from the start vertex
            prob <- if (lambda > 0) lambda ^ - D[seeds[cluster], ] else rep(1, ncol(D))
            if (sum(prob > 1) > 0) prob <- prob / max(prob) # ensure that no probability are greater than 1
            if (sum(prob > 0) < h + 1) prob[prob == 0] <- min(prob[prob != 0]) # ensure that there aren't too few positive probablility
            prob[seeds] <- 0 # ensure that none of the seeds are chosen
            prob[unlist(hits)] <- 0 # ensure that none of the previously chosen hits are chosen
            
            # choose start.vertex and h-1 other vertices
            hits[[cluster]] <- c(seeds[cluster], sample(as.numeric(V(g)), size=h - 1, prob=prob))
        }
        
        # add the hits and colour to the network
        hit.weights <- rep(0, vcount(g))
        color <- rep("grey", vcount(g))
        hit.weights[1:vcount(g) %in% unlist(hits)] <- 1
        color[1:vcount(g) %in% unlist(hits)] <- hit.color
        g <- set.vertex.attribute(g, name="hits", value=hit.weights)
        g <- set.vertex.attribute(g, name="color", value=color)
    } 
     
    # if both the hits have been successfully applied, return the graph
    if (seed.flag) return(g) else return(NULL)
}
