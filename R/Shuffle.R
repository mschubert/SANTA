Shuffle <- function(
    x, 
    ignore=NULL
) {
    # shuffle the values in x
    # ignore should be a logical vector equal in length to x
    # for each T value in ignore, the corresponding value in x is not shuffled
    if (is.null(ignore)) ignore <- rep(F, length(x))
    x[!ignore] <- sample(x[!ignore])
    x    
}
