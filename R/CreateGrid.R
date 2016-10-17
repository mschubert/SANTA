CreateGrid <- function(
    n=100
) {
    # create a grid-like graph containing n vertices
    nr <- floor(sqrt(n))
    g <- graph.empty(n, directed=F) 
    g <- add.edges(g, rbind(1:(n-1), 2:(n))[, c(rep(T, nr-1), F)]) # add edges from i to i+1 (within layer)
    g <- add.edges(g, as.numeric(rbind(1:(n-nr), (nr+1):n))) # add edges between layers
    g
}
