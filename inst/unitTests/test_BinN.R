test_BinN <- function() {
  # test that BinN splits x by a factor of x correctly
  checkEquals(SANTA:::BinN(6, 2), c(3, 3))
  
  # test that BinN splits x by a non-factor of x correctly
  checkEquals(SANTA:::BinN(12, 5), c(3, 3, 2, 2, 2))
  
  # test that BinN splits -x by a non-factor of -x correctly
  checkEquals(SANTA:::BinN(-12, 5), c(-2, -2, -2, -3, -3))
  
  # test that BinN split x=0 correctly
  checkEquals(SANTA:::BinN(0, 3), c(0, 0, 0))
  
  # test that BinN produces an error when x is not a whole number
  checkException(SANTA:::BinN(6.5, 2), silent=TRUE)
  
  # test that BinN produces an error when n is not a whole number
  checkException(SANTA:::BinN(6, 2.5), silent=TRUE)
  
  # test that BinN produces an error when n=0
  checkException(SANTA:::BinN(6, 0), silent=TRUE)
}
