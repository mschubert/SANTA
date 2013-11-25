test_BinGraph <- function() {
    # test that function returns a matrix of the same unequal dimensions when a matrix of unequal dimensions is input
    D1 <- matrix(c(1, 2, 3, 6, 7, 8), 3, 2)
    checkEquals(dim(BinGraph(D1, dist.method="shortest.paths", nsteps=100)), dim(D1))
    
    # test that function splits D correctly (ignores nsteps) when D contains only whole numbers and the dist.method is shortest.paths
    D2 <- matrix(c(0, 1, 2, 4, 5, 6, 8, 9, 10), 3, 3)
    checkEquals(BinGraph(D2, dist.method="shortest.paths", nsteps=100), matrix(c(1, 2, 3, 5, 6, 7, 9, 10, 11), 3, 3))
    
    # test that function splits D correctly when D contains only whole numbers and the dist.method is not shortest.paths
    checkEquals(BinGraph(D2, dist.method="diffusion", nsteps=100), matrix(c(1, 10, 20, 40, 50, 60, 80, 90, 100), 3, 3))
}
