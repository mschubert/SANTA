MarkovCentrality <- function(
    g, 
    edge.attr = NULL
) {
    # calculate the importance of each vertex using the Markov Centrality method
    
    # calculate the MFPT between every vertex pair in both directions
    M <- GraphMFPT(g, edge.attr=edge.attr, average.distances=F)
    
    # if there are unconnected clusters, the mfpt between unconnected vertices is Inf
    # replace these values with a value twice as high as the maximum finite value measured
    M[!is.finite(M)] <- 2 * max(M[is.finite(M)])
    
    1 / apply(M, 2, mean)
}
