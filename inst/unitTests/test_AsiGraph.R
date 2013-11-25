test_AsiGraph <- function() {
    # setup
    g <- erdos.renyi.game(5, 1)
    v <- c(1, 2, 3)
    res <- SANTA:::AsiGraph(v, g)
    
    checkEquals(class(res), "igraph.vs") # test that AsiGraph correctly returns vector an igraph.vs object
    checkEquals(length(res), length(v)) # test that AsiGraph correctly returns an igraph.vs object of length x if a vector of lenght x is input
    checkEquals(SANTA:::AsiGraph(res), res) # test that AsiGraph returns object without modificaction if already an igraph.vs object
}
