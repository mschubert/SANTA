BinGraph <- function(
    D, 
    nsteps=1000,
    equal.bin.fill=TRUE,
    verbose=TRUE
) {
    # split a distance matrix into bins
    
    # setup
    if (verbose) message("computing graph distance bins... ", appendLF=F)
    if (!isSymmetric(D)) stop("D is not symmetric")
    if (nsteps < 2) stop("nsteps is < 2")
    Dv <- D[lower.tri(D, diag=T)] # using only the lower triangle speeds computation
    
    if (length(unique(Dv)) <= nsteps) {
        # if there are fewer than nsteps + 1 unique D values, use each unique value as a break
        breaks <- unique(c(-1, 0, sort(unique(Dv))))
    } else {
        # if there are more unique distances than steps, compute suitable breaks
        
        if (equal.bin.fill) {
            # identify breaks that create bins with rougly equal numbers of nodes
            Dv <- sort(Dv)
            
            # place breaks along the sorted distances
            # if there are duplicate breaks, additional breaks are added, until there is the desired number of unique breaks
            n.breaks <- nsteps + 1 # the number of breaks introduced (may contain duplicates)
            n.breaks.unique <- 0
            while (n.breaks.unique < nsteps + 1) {
                breaks <- unique(c(-Inf, 0, Dv[ceiling(seq(1, length(Dv), length.out=n.breaks)[-1])])) # always include 0 as a break
                breaks[length(breaks)] <- Inf
                n.breaks.unique <- length(unique(breaks))
                if (n.breaks.unique < nsteps + 1) n.breaks <- 2 * n.breaks - n.breaks.unique 
                if (n.breaks.unique > nsteps + 1) breaks <- sort(c(breaks[c(1,2)], sample(unique(breaks[c(-1, -2)]), nsteps-1)))
            }
        } else {
            # equally split breaks across the range
            breaks <- c(-Inf, seq(0, max(Dv), length.out=nsteps-1))
            breaks[length(breaks)] <- Inf
        }
    }
    
    # use the breaks to create binned distance matrix
    B <- array(as.numeric(cut(D, breaks=breaks, include.lowest=T, right=T)), dim=dim(D))
    
    # format B
    B <- apply(B, 2, as.integer)
    dimnames(B) <- dimnames(D)
    
    if (verbose) message("done")
    B    
}
