test_IsWholeNumber <- function() {
  # test that IsWholeNumber correctly identifies 2 as being a whole number
  checkTrue(SANTA:::IsWholeNumber(2))
  
  # test that IsWholeNumber correctly identifies 2.5 as not being a whole number
  checkTrue(!SANTA:::IsWholeNumber(2.5))
  
  # test that IsWholeNumber correctly identifies 0 as being a whole number
  checkTrue(SANTA:::IsWholeNumber(0))
  
  # test that IsWholeNumber correctly identifies -2 as being a whole number
  checkTrue(SANTA:::IsWholeNumber(-2))
}
