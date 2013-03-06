SpreadHits <- function(
  g,
  h=10,
  lambda=1,
  dist.method="shortest.paths",
  edge.attr=NULL,
  start.vertex=NULL,
  hit.color="red",
  D=NULL
){
  # spread hits across a graph
  
  # setup
	if (class(g) != "igraph") stop("g is not a graph")
	if (lambda < 0) stop("lambda needs to be greater than or equal to 0")
  
  nvertices	<- vcount(g)
  hits      <- rep(0, nvertices)
  pheno     <- rep(0, nvertices)
  color     <- rep("grey", nvertices)
  
	# randomly choose start vertex
	if (is.null(start.vertex)) start.vertex <- sample(as.character(1:nvertices), 1)

  # spread the hits across the graph, probability determined by distance from the start vertex
  if (is.null(D)) D <- DistGraph(g, edge.attr=edge.attr, dist.method=dist.method)
  
	if (lambda > 0) {
    prob <- lambda*exp(-lambda * D[start.vertex, ]) 
  } else {
    prob <- rep(1, ncol(D))	
  }
  if (sum(prob > 0) < h + 1) prob[prob == 0] <- min(prob[prob != 0])
  prob[as.numeric(start.vertex)] <- 0  
	
  # choose start.vertex and h-1 other vertices
  sampled        <- c(start.vertex, sample(as.character(V(g)), size=h-1, prob=prob))
  sampled        <- as.numeric(sampled)
  
  # add hits, pheno and color attributes
  hits[sampled]	 <- pheno[sampled] <- 1
  color[sampled] <- hit.color	
  g              <- set.vertex.attribute(g, name="hits", value=hits)
  g              <- set.vertex.attribute(g, name="pheno", value=pheno)
  g              <- set.vertex.attribute(g, name="color", value=color)
  g
} 
