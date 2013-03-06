test_IsWholeNumber <- function() {
  # test that IsWholeNumber correctly identifies 2 as being a whole number
  checkEquals(SANTA:::IsWholeNumber(2), TRUE)
  
  # test that IsWholeNumber correctly identifies 2.5 as not being a whole number
  checkEquals(SANTA:::IsWholeNumber(2.5), FALSE)
  
  # test that IsWholeNumber correctly identifies 0 as being a whole number
  checkEquals(SANTA:::IsWholeNumber(0), TRUE)
  
  # test that IsWholeNumber correctly identifies -2 as being a whole number
  checkEquals(SANTA:::IsWholeNumber(-2), TRUE)
}
