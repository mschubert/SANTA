test_Knet <- function() {
    # setup 
    # graph used:
    # 2 - 1 - 4 - 5
    #     | 
    #     3
    
    g <- graph.empty(5, directed=F)
    g <- add.edges(g, c(1,2,1,3,1,4,4,5))
    g <- set.vertex.attribute(g, "vweights1", value=c(1, 0, 0, 1, 0))
    g <- set.vertex.attribute(g, "vweights2", value=c(0, 1, 0, 0, 1))
    g <- set.vertex.attribute(g, "vweights3", value=c(0.9, 0.1, 0.8, 0.05, 0.15)) # clustering
    g <- set.vertex.attribute(g, "vweights4", value=c(0.1, 0.9, 0.05, 0.15, 0.8)) # no clustering
    g <- set.vertex.attribute(g, "vweights5", value=c(1, 0, 1, NA, NA)) 
    g <- set.edge.attribute(g, "eweights1", value=c(1, 1, 0.01, 1))
    g <- set.edge.attribute(g, "eweights2", value=c(0.01, 0.01, 1, 0.01))
    tol <- 10e-10
    
    # Knet correctly incorperates different binary vertex weights 
    # all edge distances equal 1 as no attribute specified
    suppressMessages(res1 <- Knet(g, nperm=0, vertex.attr="vweights1"))
    suppressMessages(res2 <- Knet(g, nperm=0, vertex.attr="vweights2"))
    checkTrue(tail(res1$K.obs, 1) < tol) # check that the last observed K is basically 0
    checkTrue(tail(res2$K.obs, 1) < tol) # check that the last observed K is basically 0
    checkEquals(res1$AUK.obs, 0.3)
    checkEquals(res2$AUK.obs, 0.1)
   
    # Knet correctly incorperates different continuous vertex weights
    suppressMessages(res3 <- Knet(g, nperm=0, vertex.attr="vweights3"))
    suppressMessages(res4 <- Knet(g, nperm=0, vertex.attr="vweights4"))
    checkTrue(tail(res3$K.obs, 1) < tol) # check that the last observed K is basically 0
    checkTrue(tail(res4$K.obs, 1) < tol) # check that the last observed K is basically 0
    checkEquals(res3$AUK.obs, 0.220625)
    checkEquals(res4$AUK.obs, 0.04875)
     
    # Knet correctly incorperates different edge distances when using the different distance measures
    dist.methods <- c("shortest.paths", "diffusion", "mfpt")
    for (dist.method in dist.methods) {
        suppressMessages(res3 <- Knet(g, nperm=0, dist.method=dist.method, vertex.attr="vweights1", edge.attr="eweights1"))
        suppressMessages(res4 <- Knet(g, nperm=0, dist.method=dist.method, vertex.attr="vweights1", edge.attr="eweights2"))
        checkTrue(res3$AUK.obs > res4$AUK.obs)
    }
    
    # Knet does not shuffle the weights of vertices with weight NA in permutations
    # if this is true, then only AUKs of 0.3 and 0.4 should be seen in the permutations
    # if the weights af these vertices are permuted, then values such as 0.1 would appear
    suppressMessages(res5 <- Knet(g, nperm=100, vertex.attr="vweights5"))
    checkTrue(all(res5$AUK.perm - 0.4 < tol | res5$AUK.perm - 0.3 < tol))
    
    # Knet outputs a list of results containing the correct elements when one vertex attribute is input
    suppressMessages(res6 <- Knet(g, nperm=10, vertex.attr="vweights1", edge.attr="eweights1"))
    checkEquals(names(res6), c("K.obs", "AUK.obs", "K.perm", "AUK.perm", "K.quan", "pval"))
    
    # Knet outputs a list of lists of results containing the correct elements when multiple vertex attributes are input
    vertex.attributes <- c("vweights1", "vweights2")
    suppressMessages(res7 <- Knet(g, nperm=10, vertex.attr=vertex.attributes, edge.attr="eweights1"))
    checkEquals(length(res7), length(vertex.attributes))
    checkEquals(names(res7[[1]]), c("K.obs", "AUK.obs", "K.perm", "AUK.perm", "K.quan", "pval"))
    
    # Knet outputs a list of results, some of which equal NA, when one attribute and no permutations are completed
    suppressMessages(res8 <- Knet(g, nperm=0, vertex.attr="vweights1", edge.attr="eweights1"))
    checkTrue(all(is.na(res8[c("pval", "K.perm", "AUK.perm", "K.quan")])))
    
#     # test that parallel computing works
#     # commented out as it needs to be run on a computer where parallel computing is available
#     suppressMessages(res9 <- Knet(g, nperm=100, vertex.attr="vweights1", parallel=4))
#     checkEquals(res9$AUK.obs, 0.3)
#     checkEquals(round(res9$nodeAUK, 10), c(0.25, 0.05, 0.05, 0.35, 0.15))
}
