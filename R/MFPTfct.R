MFPTfct <- function(
    adj
) {
    # calculate the mean first-passage time (MFPT) for a fully connected graph from the adjacency matrix
    # note: this function is unable to deal with graphs that are not fully connected
    
    # if the adjacency matrix contains only a single gene, return a matrix of 1x1 containing 0
    if (is.null(dim(adj))) return(matrix(0, 1, 1))
    ngenes <- nrow(adj)
    
    A <- adj / apply(adj, 1, sum)  # A: the transition probability matrix (there is always movement)
    I <- Diagonal(x=rep(as.integer(1), ngenes))
    pi <- as.numeric(rep(1/ngenes, ngenes) %*% solve(I - A + 1/ngenes)) # pi: the stationary distribution of the transition matrix
    Z <- solve(t(t(I - A) - pi))
    as.matrix(t(t(I - Z) + Z[cbind(1:nrow(Z), 1:nrow(Z))]) %*% (I * (1 / pi))) # M: the mean first passage matrix
}
