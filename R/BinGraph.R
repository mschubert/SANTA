BinGraph <- function(
  D, 
  dist.method, 
  nsteps 
) {
  # split a distance matrix into nsteps
  lo <- ifelse(dist.method == "shortest.paths" & all(IsWholeNumber(D)), min(nsteps, max(D) + 2), nsteps + 1)	
  breaks <- seq(from=0, to=max(D) + max(D) / nsteps, length.out=lo)	
  B <- matrix(cut(D, breaks, labels=FALSE, right=FALSE), dim(D))
  B
}
