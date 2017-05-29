Shuffle <- function(
    x, 
    degrees=NULL,
    ignore=NULL
) {
    # shuffle the values in x
    # degrees should be NULL or a numerical vector equal in length to x
    # ignore should be NULL or a logical vector equal in length to x
    # if degrees is specified, values are only shuffled with values of with the same degree
    # if ignore is specified, for each T value in ignore, the corresponding value in x is not shuffled
    if (is.null(degrees)) degrees <- rep(1, length(x))
    if (is.null(ignore)) ignore <- rep(F, length(x))
    degrees.uniq <- unique(degrees)
    for (degree in degrees.uniq) {
        to.shuffle <- !ignore & (degrees == degree)
        if (sum(to.shuffle) > 1) x[to.shuffle] <- sample(x[to.shuffle])
    }
    x    
}
