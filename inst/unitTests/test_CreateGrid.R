test_CreateGrid <- function() {
    # test that function creates a square graph with the correct number of edges
    g <- CreateGrid(9)
    checkEquals(vcount(g), 9)
    checkEquals(ecount(g), 12)
    
    # test that function creates a non-square graph with the correct number of edges
    g <- CreateGrid(11)
    checkEquals(vcount(g), 11)
    checkEquals(ecount(g), 15)
}
