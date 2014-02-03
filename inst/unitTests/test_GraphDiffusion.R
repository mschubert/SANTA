test_GraphDiffusion <- function() { 
    # TESTS 
    # 1) returns matrix of dim 1,1 if networks with 1 vertex input
    # 2) uses v vertex names (under the attribute name 'name') as row names if available
    # 3) uses g vertex names (under the attribute name 'name') as column names if available
    # 4) returns only the rows corresultsponding to v, in the order v was input
    # 5) correctly incorperates edge weights
    # 6) returns both the distance matrix and the kernel 
    # 7) doesn't return negative distances when correct.neg is TRUE
    # 8) returns the correct results on a small test graph
    
    # setup
    g.single <- graph.empty(1, directed=F)
    
    edge.attr <- "test.distances"
    edges <- c(1,4, 1,8, 1,9, 1,10, 1,11, 2,5, 2,6, 2,9, 3,7, 3,10, 3,12, 4,2, 4,8, 4,11, 5,3, 5,6, 5,9, 6,9, 7,10, 7,11, 7,12, 8,11, 9,10, 10,12) 
    distances <- rep(1, length(edges) / 2)
    distances[c(2, 7, 13, 16, 18, 22)] <- 10 # move 8 and 6 away
    g <- graph.empty(max(edges), directed=F)
    g <- add.edges(g, edges)
    g <- set.edge.attribute(g, edge.attr, value=distances)
    g.unnamed <- g
    g.named <- set.vertex.attribute(g, "name", value=paste("gene", 1:vcount(g), sep=""))
    
    n.vertex.unconnected <- 10
    g.unconnected <- graph.empty(n.vertex.unconnected, directed=F)
    
    v <- c(4,2,12)
    
    
    # run function
    results <- list()
    results[[1]] <- GraphDiffusion(g.single) # 1 vertex
    results[[2]] <- GraphDiffusion(g.unnamed) # >1 vertex, v = V(g), no vertex names
    results[[3]] <- GraphDiffusion(g.named) # >1 vertex, v = V(g), vertex names
    results[[4]] <- GraphDiffusion(g.named, v=v) # >1 vertex, v = subset, vertex names
    results[[5]] <- GraphDiffusion(g.named, edge.attr=edge.attr) # >1 vertex, v = V(g), no vertex name, edge weights
    results[[6]] <- GraphDiffusion(g.unconnected) # unconnected network
    
    
    # conduct tests
    checkTrue(all(sapply(results[[1]], function(element) dim(element) - c(1,1) < 10e-10))) # 1
    
    checkTrue(is.null(rownames(results[[2]]$dist))) # 2
    checkIdentical(rownames(results[[3]]$dist), V(g.named)$name) # 2
    checkIdentical(rownames(results[[4]]$dist), V(g.named)$name[v]) # 2
    
    checkTrue(is.null(colnames(results[[2]]$dist))) # 3
    checkIdentical(colnames(results[[3]]$dist), V(g.named)$name) # 3
    checkIdentical(colnames(results[[4]]$dist), V(g.named)$name) # 3
    
    checkIdentical(results[[4]]$kernel, results[[3]]$kernel[v, ]) # 4
    checkIdentical(results[[4]]$dist, results[[3]]$dist[v, ]) # 4
    
    checkTrue(all(c(8, 6) %in% tail(order(results[[5]]$dist[1, ]), 2))) # 5
    
    for (result in results) {
        checkIdentical(names(result), c("kernel", "dist")) # 6
        checkTrue(all(sapply(result, function(element) is.matrix(element)))) # 6
     
        checkEquals(sum(result$dist < 0), 0) # 7
    }
    
    for (result in results[c(2, 3)]) {
        checkTrue(all(c(1, 4, 8, 11) %in% head(order(result$dist[1, ]), 4))) # 8
        checkTrue(all(c(2, 5, 6, 9) %in% head(order(result$dist[2, ]), 4))) # 8
        checkTrue(all(c(3, 7, 10, 12) %in% head(order(result$dist[3, ]), 4))) # 8
    }
}
