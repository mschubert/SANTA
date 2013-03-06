IsWholeNumber <- function(x, tol=.Machine$double.eps ^ 0.5) {  
  # check whether a number is a whole number
  abs(x - round(x)) < tol
}
