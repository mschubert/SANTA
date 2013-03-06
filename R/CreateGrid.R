CreateGrid <- function(
  n=100
) {
  # create a grid-like graph containing n vertices
    
  # create empty graph
  g	<- graph.empty(n, directed=FALSE) 
    
  # add edges from i to i+1 (within layer)
  nr	<- floor(sqrt(n))
  g	<- add.edges(g, rbind(1:(n-1), 2:(n))[, c(rep(TRUE, nr-1), FALSE)])
    
  # add edges between layers
  g	<- add.edges(g, as.numeric(rbind(1:(n-nr), (nr+1):n)))
  g
}
