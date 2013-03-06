GraphDiffusion <- function(
  g,
  v=V(g), 
  edge.attr=NULL,
  beta=1, 
  correct.factor=1,
  correct.neg=TRUE
) { 
  # calculates the diffusion distance using the diffusion kernel method
  # diffusion Kernels on Graphs and Other Discrete Structures, Kondor and Lafferty, 2002

  # if v is not an igraph-class object, convert
  if (class(v) != "igraph.vs") v <- AsiGraph(v, g)

  # if there is only a single vertex contained within the graph, return a single cell matrix
	if (vcount(g) == 1) return(list(kernel=matrix(0,1,1), dist=matrix(0, 1, 1)))

  # invert the adjacency matrix, so that small distances result in a greater chance of movement along the edge
  adj  <- get.adjacency(g, attr=edge.attr, sparse=FALSE)
  adj[adj > 0] <- max(adj) - adj[adj > 0] + min(adj[adj > 0])
  
  D  <- diag(apply(adj, 1, sum))
  H  <- adj - D	
  x  <- eigen(H)
  K  <- x$vectors %*% diag(exp(beta * x$values)) %*% t(x$vectors)
  D2 <- outer(diag(K),diag(K), "+") - 2*K
  
  if (any(D2 < 0) && correct.neg){
    warning("negative distances set to zero")
    D2[D2 < 0] <- 0
  }
  D <- sqrt(D2)
  
  K <- as.matrix(K[v, ])
  D <- as.matrix(D[v, ])
  rownames(D) <- v
  
  # attribute unconnected vertex pairs with correct.factor * the maximum distance across connected vertex pairs
  connected <- matrix(TRUE, nrow(D), ncol(D))
  clus.members <- clusters(g)$membership
  for (i in 1:max(clus.members)) {
    tmp <- which(clus.members == i)
    connected[tmp, tmp] <- FALSE
  }
  D[connected] <- correct.factor * max(D[!connected])
  
  list(kernel=K, dist=D)
}
