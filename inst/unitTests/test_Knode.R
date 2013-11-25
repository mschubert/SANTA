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
        suppressMessages(res1 <- Knode(g, dist.method=dist.method, vertex.attr="vweights1", only.Knode=T))
        checkTrue(rownames(res1)[1] == "4")
    }
    
    # test that Knode correctly incorperates different vertex weights correctly when using different distance methods
    dist.methods <- c("shortest.paths", "diffusion", "mfpt")
    for (dist.method in dist.methods) {
        suppressMessages(res2 <- Knode(g, dist.method=dist.method, vertex.attr="vweights1", only.Knode=T))
        suppressMessages(res3 <- Knode(g, dist.method=dist.method, vertex.attr="vweights2", only.Knode=T))
        checkTrue(rownames(res2)[1] == "4" & rownames(res3)[1] == "5")
    }
    
    # test that Knode correctly incorperates different edge distances when using the different distance measures
    dist.methods <- c("shortest.paths", "diffusion", "mfpt")
    for (dist.method in dist.methods) {
        suppressMessages(res4 <- Knode(g, dist.method=dist.method, vertex.attr="vweights1", edge.attr="eweights1", only.Knode=T))
        suppressMessages(res5 <- Knode(g, dist.method=dist.method, vertex.attr="vweights1", edge.attr="eweights2", only.Knode=T))
        order4 <- rownames(res4)
        order5 <- rownames(res5)
        checkTrue(which(order4=="2") < which(order4=="3") & which(order5=="3") < which(order5=="2"))
    }
    
    # test that Knode outputs a list of results containing the correct elements when one vertex attribute is input
    suppressMessages(res6 <- Knode(g, vertex.attr="vweights1", only.Knode=T))
    suppressMessages(res7 <- Knode(g, vertex.attr="vweights1", only.Knode=F))
    checkTrue(nrow(res6) == vcount(g) & nrow(res7) == vcount(g))
    checkTrue(ncol(res6) == 1 & ncol(res7) > 1)
    
    # test that Knode outputs a list of lists when multiple vertex attributes are input
    vertex.attributes <- c("vweights1", "vweights2")
    suppressMessages(res8 <- Knode(g, vertex.attr=vertex.attributes, only.Knode=T))
    checkEquals(length(res8), length(vertex.attributes))
}
