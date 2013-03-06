CreateGraph <- function(
  n=100, 
  type="barabasi",
  m=2,
  p.or.m=(2*n-2*floor(sqrt(n))),
  vertex.weights=NULL,
  edge.distances=NULL, 
  gen.vertex.weights=FALSE,
  nclusters=1,
  lambda=1,
  nlinks=ceiling(n/5),
  nhits=ceiling(n/10),
  binary.pheno=TRUE,
  dist.method="shortest.paths",
  mean.hit=1,
  sd.hit=0.05,
  mean.miss=0,
  sd.miss=0.05
) {
  # create an undirected and simple graph according to specified graph types and parameters

  GraphFct <- function(n, type, m, p.or.m, edge.distances) {
    g <- switch(type,
      barabasi    = simplify(barabasi.game(n, m=m, directed=FALSE)),
      erdos.renyi = erdos.renyi.game(n, p.or.m=p.or.m, type="gnm", directed=FALSE),
      grid        = CreateGrid(n)	
    )
    
    if (!is.null(edge.distances)) {
      if(length(edge.distances)!=ecount(g)) stop("Number of edge distances provided not equal to number of edges in graph")
      g <- set.edge.attribute(g, name="distance", value=sample(edge.distances))
    }
    g
  }
    
  # setup
  if (!type %in% c("barabasi", "erdos.renyi", "grid")) stop("Input type unknown")

  if (is.null(vertex.weights) & gen.vertex.weights == FALSE) {
    # generate a simple graph without vertex weights
    g <- GraphFct(n, type, m, p.or.m, edge.distances)
  } else {
    # generate a graph with 1 or more clusters, using the input vertex weights or generated vertex weights

    # calculate the number of vertices and hits to be included within each constituent graph
    nvertices.binned <- BinN(n, nclusters)
    nhits.binned     <- BinN(nhits, nclusters)
    id.boundaries    <- cbind(head(c(1, cumsum(nvertices.binned)+1),-1), cumsum(nvertices.binned)) 

    # create ncluster different constituent graphs
    graphs <- vector("list", nclusters)
    for (i in 1:nclusters) {
      # create graphs and add hits - a hit is equal to 1 while a miss is equal to 0
      graphs[[i]]    <- GraphFct(nvertices.binned[i], type, m, p.or.m, edge.distances=NULL)
      graphs[[i]]    <- SpreadHits(graphs[[i]], h=nhits.binned[i], lambda=lambda, dist.method=dist.method)
      
      # record which cluster the hit originated from
      graphs[[i]]    <- set.vertex.attribute(graphs[[i]], name="hits.cluster", value=ifelse(get.vertex.attribute(graphs[[i]], name="hits")==1, i, 0))
    }
    
    # join graphs and create links between the constituent graphs
    g <- graph.disjoint.union(graphs)
    if (nlinks > 0 & nclusters > 1) {
      # calculate the number of different connections to be made and the number of links to be used to connect each of the conjoined graphs
      nconnect <- sum(1:(nclusters-1))
      nlinks.binned <- BinN(nlinks, nconnect)
      
      # add edges connecting each of the constituent graphs
      c <- 1
      for (i in 1:(nclusters-1)) {
        for (j in (i+1):nclusters) {
          g <- simplify(add.edges(g, as.numeric(rbind(sample(id.boundaries[i,1]:id.boundaries[i,2], nlinks.binned[c], replace=TRUE), sample(id.boundaries[j,1]:id.boundaries[j,2], nlinks.binned[c], replace=TRUE)))))
          c <- c+1
        }
      }	
    }
    
    # add vertex attributes to the composite graph
    g     <- set.vertex.attribute(g, name="hits", value=unlist(lapply(graphs, function(x) get.vertex.attribute(x, name="hits"))))
    g     <- set.vertex.attribute(g, name="hits.cluster", value=unlist(lapply(graphs, function(x) get.vertex.attribute(x, name="hits.cluster"))))
    g     <- set.vertex.attribute(g, name="color", value=unlist(lapply(graphs, function(x) get.vertex.attribute(x, name="color"))))
    whits <- which(get.vertex.attribute(g, name="hits")==1)
    
    # add edge attributes to the composite graph if input - otherwise no edge attributes are added
    if (!is.null(edge.distances)) {
      if(length(edge.distances) != ecount(g)) stop("Number of edge distances provided not equal to number of edges in graph")
      g	<- set.edge.attribute(g, name="distance", value=sample(edge.distances))
    }
      
    if (is.null(vertex.weights)) {
      # if vertex weights are not given, generate then
      if (binary.pheno) {
        # pheno equals 1 for a hit and 0 for a miss
        g             <- set.vertex.attribute(g, name="pheno", value=get.vertex.attribute(g, name="hits"))
      } else {
      	# generate phenotypic data for the hits and non-hits using truncated normal distributions
        pheno         <- double(n)
        pheno[whits]  <- rtnorm(nhits, mean=mean.hit, sd=sd.hit, lower=0, upper=1)
        pheno[-whits] <- rtnorm(n-nhits, mean=mean.miss, sd=sd.miss, lower=0, upper=1)
        g             <- set.vertex.attribute(g, name="pheno", value=pheno)
      }
    } else {
      # if vertex weights are given, apply them
      if (length(vertex.weights) != vcount(g)) stop("Number of vertex weights provided not equal to number of vertices in graph")
      vertex.weights  <- sort(vertex.weights, decreasing=TRUE)
      pheno           <- double(n)
      pheno[whits]    <- sample(vertex.weights[1:nhits])
      pheno[-whits]   <- sample(vertex.weights[(nhits+1):n])
      g               <- set.vertex.attribute(g, name="pheno", value=pheno)
    }
  }
  g
}
