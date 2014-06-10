test_Knode <- function() {
    # setup 
    # graph used:
    # 2 - 1 - 4 - 5
    #     | 
    #     3
    
    g <- graph.empty(5, directed=F)
    g <- add.edges(g, c(1, 2, 1, 3, 1, 4, 4, 5))
    g <- set.vertex.attribute(g, "vweights1", value=c(1, 0, 0, 1, 0))
    g <- set.vertex.attribute(g, "vweights2", value=c(0, 1, 0, 0, 1))
    g <- set.edge.attribute(g, "eweights1", value=c(0.01, 1, 1, 1))
    g <- set.edge.attribute(g, "eweights2", value=c(1, 0.01, 1, 1))
    
    # test that Knode correctly uses different distance methods
    dist.methods <- c("shortest.paths", "diffusion", "mfpt")
    for (dist.method in dist.methods) {
        res1 <- Knode(g, dist.method=dist.method, vertex.attr="vweights1", verbose=F)
        checkTrue(names(res1)[1] == "4")
    }
    
    # test that Knode correctly incorperates different vertex weights correctly when using different distance methods
    dist.methods <- c("shortest.paths", "diffusion", "mfpt")
    for (dist.method in dist.methods) {
        res2 <- Knode(g, dist.method=dist.method, vertex.attr="vweights1", verbose=F)
        res3 <- Knode(g, dist.method=dist.method, vertex.attr="vweights2", verbose=F)
        checkTrue(all(as.character(c(1,4)) %in% names(res2)[1:2]))
        checkTrue(all(as.character(c(2,5)) %in% names(res3)[1:2]))
    }
    
    # test that Knode correctly incorperates different edge distances when using the different distance measures
    dist.methods <- c("shortest.paths", "diffusion", "mfpt")
    for (dist.method in dist.methods) {
        res4 <- Knode(g, dist.method=dist.method, vertex.attr="vweights1", edge.attr="eweights1", verbose=F)
        res5 <- Knode(g, dist.method=dist.method, vertex.attr="vweights1", edge.attr="eweights2", verbose=F)
        order4 <- names(res4)
        order5 <- names(res5)
        checkTrue(which(order4=="2") < which(order4=="3") & which(order5=="3") < which(order5=="2"))
    }
        
    # test that Knode outputs a list of lists when multiple vertex attributes are input
    vertex.attributes <- c("vweights1", "vweights2")
    res8 <- Knode(g, vertex.attr=vertex.attributes, verbose=F)
    checkEquals(length(res8), length(vertex.attributes))
}
