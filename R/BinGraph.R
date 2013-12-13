BinGraph <- function(
    D, 
    dist.method, 
    nsteps,
    verbose = TRUE
) {
    # split a distance matrix into nsteps
    if (verbose) message("computing graph distance bins... ", appendLF=F)
    B <- if (dist.method == "shortest.paths" & all(IsWholeNumber(D))) D + 1 else (D %/% ((max(D) * (nsteps + 1)) / nsteps ^ 2)) + 1
    if (verbose) message("done")
    B
}
