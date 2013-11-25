BinGraph <- function(
    D, 
    dist.method, 
    nsteps 
) {
    # split a distance matrix into nsteps
    if (dist.method == "shortest.paths" & all(IsWholeNumber(D))) D + 1 else (D %/% ((max(D) * (nsteps + 1)) / nsteps ^2)) + 1
    
    #bin.size <- if (dist.method == "shortest.paths") 1 else (max(D) * (nsteps + 1)) / nsteps ^ 2
    #(D %/% bin.size) + 1 
}
