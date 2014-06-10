test_BinGraph <- function() {
    # TESTS
    # 1) works when D values are all 0
    # 2) works when D values are all 1
    # 3) works when D values are all either 0 or 1
    # 4) works when the number of unique D values is nsteps - 1
    # 5) works when the number of unique D values is == nsteps
    # 6) works when the number of unique D values is nsteps + 1
    # 7) works when all D values are unique
    # 8) works when all but the start D values are unique
    # 9) works when all but the middle D values are unique
    # 10) works when all but the end D values are unique
    # 11) returned B has same dimnames as D
    # 12) produces error when nsteps < 2
    # 13) produces error when D is not symmetric
    
    # 1, 11) 
    n.rowcol <- 4
    names.rowcol <- LETTERS[1:n.rowcol]
    D <- array(0, dim=c(n.rowcol, n.rowcol), dimnames=list(names.rowcol, names.rowcol))
    B <- array(apply(D + 1, 2, as.integer), dim=dim(D), dimnames=dimnames(D))
    checkIdentical(BinGraph(D, nsteps=4), B)
    checkIdentical(dimnames(BinGraph(D, nsteps=4)), dimnames(D))
    
    # 2) 
    D <- matrix(1, 4,4)
    diag(D) <- 0
    B <- array(apply(D + 1, 2, as.integer), dim=dim(D), dimnames=dimnames(D))
    checkIdentical(BinGraph(D, nsteps=4), B)
    
    # 3)
    D <- matrix(1, 4,4)
    diag(D) <- 0
    D[1,4] <- D[4,1] <- 0
    B <- array(apply(D + 1, 2, as.integer), dim=dim(D), dimnames=dimnames(D))
    checkIdentical(BinGraph(D, nsteps=4), B)
    
    # 4,5,6)
    D <- matrix(c(0,1,2,3,1,0,1,2,2,1,0,1,3,2,1,0), 4,4)
    B <- array(apply(D + 1, 2, as.integer), dim=dim(D), dimnames=dimnames(D))
    checkIdentical(BinGraph(D, nsteps=5), B) # 4)
    checkIdentical(BinGraph(D, nsteps=4), B) # 5)
    B <- matrix(as.integer(c(1,2,3,3,2,1,2,3,3,2,1,2,3,3,2,1)), 4,4)
    checkIdentical(BinGraph(D, nsteps=3), B) # 6)
    
    # 7) 
    D <- matrix(c(0,1,2,3,4,1,0,5,6,7,2,5,0,8,9,3,6,8,0,10,4,7,9,10,0), 5,5)
    B <- matrix(as.integer(c(1,2,2,2,3,2,1,3,3,3,2,3,1,4,4,2,3,4,1,4,3,3,4,4,1)), 5,5)
    checkIdentical(BinGraph(D, nsteps=4), B)
    
    # 8,9,10)
    nrow.tot <- 8 # number of rows/columns
    nrow.na <- 6 # size of NA square to be replaced with value
    D.template <- array(runif(nrow.tot ^ 2), dim=c(nrow.tot, nrow.tot))
    D.template <- (D.template + t(D.template)) / 2 # ensure D is symmetric
    D.template[cbind(rep(1:nrow.na, nrow.na), rep(1:nrow.na, each=nrow.na))] <- NA
    diag(D.template) <- 0
    nsteps <- 5
    values.na <- c(10e-10, 0.5, 1)
    for (value in values.na) {
        # check that each of the steps is represented, no matter where the group of tied vertices comes
        D <- D.template
        D[is.na(D)] <- value
        B <- BinGraph(D, nsteps)
        checkIdentical(1:nsteps, sort(unique(as.vector(B)))) 
    }
    
    # 12)
    D <- matrix(1, 4,4)
    diag(D) <- 0
    checkException(BinGraph(D, 1), silent=T)
    
    # 13)
    D <- matrix(1, 4,4)
    D[1,2] <- 2
    diag(D) <- 0
    checkException(BinGraph(D, 1), silent=T)
}
