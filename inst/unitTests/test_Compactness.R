test_Compactness <- function() {
    # setup 
    # graph used:
    # 2 - 1 - 4 - 5
    #     | 
    #     3
    
    g <- graph.empty(5, directed=F)
    g <- add.edges(g, c(1,2,1,3,1,4,4,5))
    g <- set.vertex.attribute(g, "vw1", value=c(1, 0, 0, 1, 0))
    g <- set.vertex.attribute(g, "vw2", value=c(0, 1, 0, 0, 1))
    g <- set.edge.attribute(g, "ew1", value=c(1, 1, 0.01, 1))
    g <- set.edge.attribute(g, "ew2", value=c(0.01, 0.01, 1, 0.01))
    
    # correctly incorperates different binary vertex weights 
    res1 <- Compactness(g, nperm=1000, vertex.attr="vw1")
    res2 <- Compactness(g, nperm=1000, vertex.attr="vw2")
    checkEquals(res1$score.obs, 1)
    checkEquals(res2$score.obs, 3)
    checkTrue(res1$pval < res2$pval)
    
    # correctly incorperates different edge distances when using the different distance measures
    dist.methods <- c("shortest.paths", "diffusion", "mfpt")
    for (dist.method in dist.methods) {
        res3 <- Compactness(g, nperm=1000, dist.method=dist.method, vertex.attr="vw1", edge.attr="ew1")
        res4 <- Compactness(g, nperm=1000, dist.method=dist.method, vertex.attr="vw1", edge.attr="ew2")
        checkTrue(res3$score.obs < res4$score.obs)
        checkTrue(res3$pval < res4$pval)
    }
    
    # outputs a list of results containing the correct elements when one vertex attribute is input
    res5 <- Compactness(g, nperm=10, vertex.attr="vw1", edge.attr="ew1")
    checkEquals(names(res5), c("score.obs", "score.perm", "pval") )
    
    # outputs a list of lists of results containing the correct elements when multiple vertex attributes are input
    vertex.attributes <- c("vw1", "vw2")
    res6 <- Compactness(g, nperm=10, vertex.attr=vertex.attributes, edge.attr="ew1")
    checkEquals(length(res6), length(vertex.attributes))
    for (r in res6) checkEquals(names(r), c("score.obs", "score.perm", "pval"))
    
    # outputs a list of results, some of which equal NA, when one attribute and no permutations are completed
    res7 <- Compactness(g, nperm=0, vertex.attr="vw1", edge.attr="ew1")
    checkTrue(is.na(res7$score.perm))
    checkTrue(is.na(res7$pval))
}
