test_SpreadHits <- function() {
    # setup 
    # create string of nodes
    n.nodes <- 20
    edges <- rep(1:n.nodes, each=2)[-c(1, n.nodes*2)]
    g <- graph.empty(n.nodes, directed=F)
    g <- add.edges(g, edges)
    g <- set.edge.attribute(g, "distances1", value=rep(10, ecount(g)))
   
   
    # test that SpreadHits adds the correct number of hits across 1, 2 and 3 clusters
    # test that SpreadHits adds vertex attributes "hits" and "color"
    clusters <- 1:3
    h <- 2
    for (cluster in clusters) {
        g1 <- SpreadHits(g, h=h, clusters=cluster, distance.cutoff=2, lambda=10, dist.method="shortest.paths")
        hits <- which(get.vertex.attribute(g1, "hits") == 1)
        
        checkEquals(length(hits), h*cluster)
        checkTrue(identical(sort(list.vertex.attributes(g1)), c("color", "hits")))
    }

    # test that SpreadHits uses specified D when input
    h <- 2
    clusters <- 3
    D <- DistGraph(g, edge.attr="distances1", dist.method="shortest.paths")
    g2 <- SpreadHits(g, h=h, clusters=clusters, distance.cutoff=20, lambda=10, dist.method="shortest.paths", D=D)
    hits <- which(get.vertex.attribute(g2, "hits") == 1)
    checkEquals(length(hits), h*clusters)
    
    # test that SpreadHits produces an error when g is not an igraph object
    checkException(SpreadHits(6), silent=T)
    
    # test that SpreadHits produces an error when lambda is negative
    checkException(SpreadHits(g, h=h, clusters=clusters, distance.cutoff=20, lambda=-10, dist.method="shortest.paths"), silent=T)
    
    # test that SpreadHits produces error if more hits than vertices are specified to be applied
    clusters <- 2
    checkException(SpreadHits(g, h=vcount(g)/clusters + 1, clusters=clusters, distance.cutoff=2, lambda=10, dist.method="shortest.paths", edge.attr=NULL), silent=T)
    
    # test that SpreadHits returns NULL if it is not possible to identify seed nodes the required distance apart
    g4 <- SpreadHits(g, h=2, clusters=3, distance.cutoff=20, lambda=10, dist.method="shortest.paths", edge.attr=NULL)
    checkTrue(is.null(g4))
}
