Knode <- function(
  g, 
  dist.method    ="shortest.paths", 
  vertex.attr    ="pheno", 
  edge.attr      ="distance",
  correct.factor =1.0,
  nsteps         =1000,
  only.Knode     =TRUE,
  vertex.weight  =TRUE,
  cluster.id     =FALSE,
  vertex.degree  =TRUE,
  boncich.power  =FALSE,
  markov.centr   =FALSE
) {  
  # rank the genes exhibiting the phenotype by Knode AUK score. Also produce and list a number of other score for each vertex
  
  # check that vertex weights and edge distances are present and suitable. Convert if neccessary
  g <- CheckAttributes(g, vertex.attr, edge.attr)
  
  # if vertices do not have a name, then add one. 
  if (is.null(get.vertex.attribute(g, "name"))) g <- set.vertex.attribute(g, "name", value=as.character(1:vcount(g)))
  g.names <- get.vertex.attribute(g, "name") 
  
  # compute the vertex pair distances (D)
  message("computing graph distances...", appendLF=FALSE)
  D <- DistGraph(g=g, edge.attr=edge.attr, dist.method=dist.method, correct.inf=TRUE, correct.factor=correct.factor)  
  message(" done")
  
  # compute which bin each vertex pair distance falls into (B)
  message("computing graph distance bins...", appendLF=FALSE)
  B <- BinGraph(D=D, dist.method=dist.method, nsteps=nsteps)
  message(" done")
  
  if (!only.Knode) {
    # if required, compute each of the additional graph stats
    stats <- list()
    if(cluster.id) stats[["cluster.id"]]          <- get.vertex.attribute(g, "hits.cluster")
    if(vertex.degree) stats[["vertex.degree"]]    <- degree(graph=g)
    stats[["vertex.betweenness"]]                 <- betweenness(graph=g, directed=FALSE)
    if(boncich.power) stats[["bonpow"]]           <- bonpow(graph=g)
    stats[["constraint"]]                         <- constraint(graph=g, weights=get.edge.attribute(g, edge.attr))
    stats[["evcent"]]                             <- evcent(graph=g, weights=get.edge.attribute(g, edge.attr))$vector
    stats[["pagerank"]]                           <- page.rank(graph=g, directed=FALSE, weights=get.edge.attribute(g, edge.attr))$vector
    stats[["authority"]]                          <- authority.score(graph=g)$vector
    stats[["hub"]]                                <- hub.score(graph=g)$vector
    if(markov.centr) stats[["markov.centrality"]] <- MarkovCentrality(g=g, edge.attr=edge.attr)
  } 
  
  # the results for each vertex attrbitute are saved as a data frame. Theses data frames are stored in a list
  res <- vector("list", length(vertex.attr))
  names(res)	<- vertex.attr
  
  # compute the Knode scores for each of the vertex attributes in vertex.attr
  for (attr in vertex.attr) {
    if (length(vertex.attr)==1) {
      message(paste("computing the Knode scores of vertices using '", attr, "' as weights...", sep=""), appendLF=FALSE)
    } else {
      message(paste("computing the Knode scores of vertices using '", attr, "' (", which(vertex.attr==attr), "/", length(vertex.attr) ,") as weights...", sep=""), appendLF=FALSE)
    }
    
    vertex.weights <- get.vertex.attribute(g, attr)
    nodeAUK        <- Kfct(B, vertex.weights, individual=TRUE)$nodeAUK
    
    if (!only.Knode) {
      # combine the Kfct results with the statistics
      if (vertex.weight) {
        res[[attr]] <- cbind(nodeAUK, vertex.weights, data.frame(stats))	
      } else {
        res[[attr]] <- cbind(nodeAUK, data.frame(stats))
      }
      
      # sort the data frame of Knode scores and statistics by Knode score
      res[[attr]] <-res[[attr]][order(nodeAUK, decreasing=TRUE), ]
    } else {
      # sort the Knode scores
      names(nodeAUK) <- g.names 
      res[[attr]] <- data.frame(sort(nodeAUK, decreasing=TRUE))
    }
    message(" done")
  }
  
  # if only one vertex attribute is input, don't return a list of lists of results
  if (length(vertex.attr)==1) res <- res[[1]]
  
  # output
  res
} 
