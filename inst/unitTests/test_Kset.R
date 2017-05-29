test_Kset <- function() {
    # TESTS
    # 1) Kset correctly uses different binary vertex weights
    # 2) Kset correctly uses different continuous vertex weights
    # 3) Kset correctly incoperates different edge distances when using different distance measures
    # 4) Kset does not shuffle the weights of vertices with weight NA in permutations
    # 5) Kset outputs a list of results containing the correct elements when one query vertex attribute is input
    # 6) Kset outputs a list of lists of results containing the correct elements when multiple quwey vertex attributes are input
    # 7) Kset outputs a list of results, some of which equal NA, when one attribute and no permutations are completed
    
    # setup 
    # graph used:
    # 2 - 1 - 4 - 5
    #     | 
    #     3
    
    g <- graph.empty(5, directed=F)
    g <- add.edges(g, c(1,2,1,3,1,4,4,5))
    g <- set.vertex.attribute(g, "seed", value=c(1, 0, 0, 1, 0))
    g <- set.vertex.attribute(g, "query1", value=c(1, 0, 0, 1, 0)) # clustering 
    g <- set.vertex.attribute(g, "query2", value=c(1, 0, 0, 0, 1)) # no clustering
    g <- set.vertex.attribute(g, "query3", value=c(0.90, 0.10, 0.05, 0.80, 0.15)) # clustering
    g <- set.vertex.attribute(g, "query4", value=c(0.10, 0.90, 0.05, 0.15, 0.80)) # no clustering
    g <- set.vertex.attribute(g, "query5", value=c(1, 0, 1, NA, NA)) 
    g <- set.edge.attribute(g, "weights1", value=c(1.00, 0.01, 0.01, 1.00))
    g <- set.edge.attribute(g, "weights2", value=c(0.01, 1.00, 1.00, 0.01))
    tol <- 10e-10
    
    # Kset (version A) correctly incorperates different binary vertex weights 
    # all edge distances equal 1 as no attribute specified
    suppressMessages(res1 <- Kset(g, nperm=100, seed.vertex.attr="seed", query.vertex.attr="query1", version="a"))
    suppressMessages(res2 <- Kset(g, nperm=100, seed.vertex.attr="seed", query.vertex.attr="query2", version="a"))
    checkTrue(abs(tail(res1$K.obs, 1)) < tol) # check that the last observed K is basically 0
    checkTrue(abs(tail(res2$K.obs, 1)) < tol) # check that the last observed K is basically 0
    checkEquals(res1$AUK.obs, 0.30) # calculated by hand
    checkEquals(res2$AUK.obs, 0.05) # calculated by hand
   
    # Kset (version B) correctly incorperates different binary vertex weights 
    # all edge distances equal 1 as no attribute specified
    suppressMessages(res1 <- Kset(g, nperm=100, seed.vertex.attr="seed", query.vertex.attr="query1", version="b"))
    suppressMessages(res2 <- Kset(g, nperm=100, seed.vertex.attr="seed", query.vertex.attr="query2", version="b"))
    checkTrue(abs(tail(res1$K.obs, 1) - 1) < tol) # check that the last observed K is basically 1
    checkTrue(abs(tail(res2$K.obs, 1) - 1) < tol) # check that the last observed K is basically 1
    checkEquals(res1$AUK.obs, 1.00) # calculated by hand
    checkEquals(res2$AUK.obs, 0.875) # calculated by hand
    
    # Kset (version A) correctly incorperates different continuous vertex weights
    suppressMessages(res3 <- Kset(g, nperm=100, seed.vertex.attr="seed", query.vertex.attr="query3", version="a"))
    suppressMessages(res4 <- Kset(g, nperm=100, seed.vertex.attr="seed", query.vertex.attr="query4", version="a"))
    checkTrue(abs(tail(res3$K.obs, 1)) < tol) # check that the last observed K is basically 0
    checkTrue(abs(tail(res4$K.obs, 1)) < tol) # check that the last observed K is basically 0
    checkEquals(res3$AUK.obs, 0.225) # calculated by hand 
    checkEquals(res4$AUK.obs, -0.1375) # calculated by hand

    # Kset correctly incorperates different edge distances when using the different distance measures
    dist.methods <- c("shortest.paths", "diffusion", "mfpt")
    for (dist.method in dist.methods) {
        suppressMessages(res3 <- Kset(g, nperm=0, dist.method=dist.method, seed.vertex.attr="seed", query.vertex.attr="query1", edge.attr="weights1", version="a"))
        suppressMessages(res4 <- Kset(g, nperm=0, dist.method=dist.method, seed.vertex.attr="seed", query.vertex.attr="query1", edge.attr="weights2", version="a"))
        checkTrue(res3$AUK.obs > res4$AUK.obs)
    }
    
    # Kset does not shuffle the weights of vertices with weight NA in permutations
    # if this is true, then only AUKs of 0.3 and 0.4 should be seen in the permutations
    # if the weights af these vertices are permuted, then values such as 0.1 would appear
    suppressMessages(res5 <- Kset(g, nperm=100, seed.vertex.attr="seed", query.vertex.attr="query5", version="a"))
    checkTrue(all(abs(res5$AUK.perm + 0.2) < tol | abs(res5$AUK.perm - 0.05) < tol))
    
    # Kset outputs a list of results containing the correct elements when one vertex attribute is input
    suppressMessages(res6 <- Kset(g, nperm=10, seed.vertex.attr="seed", query.vertex.attr="query1", edge.attr="weights1"))
    checkEquals(names(res6), c("K.obs", "AUK.obs", "K.perm", "AUK.perm", "K.quan", "pval"))
    
    # Kset outputs a list of lists of results containing the correct elements when multiple vertex attributes are input
    query.vertex.attrs <- c("query1", "query2")
    suppressMessages(res7 <- Kset(g, nperm=10, seed.vertex.attr="seed", query.vertex.attr=query.vertex.attrs, edge.attr="weights1"))
    checkEquals(length(res7), length(query.vertex.attrs))
    checkEquals(names(res7[[1]]), c("K.obs", "AUK.obs", "K.perm", "AUK.perm", "K.quan", "pval"))
    
    # Kset outputs a list of results, some of which equal NA, when one attribute and no permutations are completed
    suppressMessages(res8 <- Kset(g, nperm=0, seed.vertex.attr="seed", query.vertex.attr="query1", edge.attr="weights1"))
    checkTrue(all(is.na(res8[c("pval", "K.perm", "AUK.perm", "K.quan")])))
}
