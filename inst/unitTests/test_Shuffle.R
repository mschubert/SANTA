test_Shuffle <- function() {
    # setup
    x <- 1:5
    x.na <- c(1, 2, 3, NA, 5)
    n.trials <- 5
    
    for (i in 1:n.trials) {
        # check that function ingores certain values
        ignore1 <- c(F, T, F, T, F)
        res1 <- SANTA:::Shuffle(x, ignore=ignore1)
        checkEquals(x[ignore1], res1[ignore1])
        
        # check that function can ignore all values
        ignore2 <- c(T, T, T, T, T)
        res2 <- SANTA:::Shuffle(x, ignore=ignore2)
        checkEquals(x[ignore2], res2[ignore2])
        
        # check that function handles NAs
        ignore3 <- c(F, T, F, T, F)
        res3 <- SANTA:::Shuffle(x.na, ignore=ignore3)
        checkEquals(x.na[ignore3], res3[ignore3])
    }
}
