BinN <- function(x, n) {
  # split an integer x into n integer bins
  
  # check that x and n are positive wholenumbers
  if (!IsWholeNumber(x)) stop("x should be a whole number")
  if (!IsWholeNumber(n) | (n <= 0)) stop("n should be a positive whole number")
  
	res <- rep(floor(x / n), n)
	res[0:(x - sum(res))] <- res[0:(x - sum(res))] + 1
	res
}
