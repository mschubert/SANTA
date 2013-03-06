test_Kfct <- function() { 
  # setup 
  # B matrix used
  #   1 2 2 2 3 
  #   2 1 3 3 4
  #   2 3 1 3 4
  #   2 3 3 1 2
  #   3 4 4 2 1
  B <- matrix(c(1,2,2,2,3,2,1,3,3,4,2,3,1,3,4,2,3,3,1,2,3,4,4,2,1), 5,5)
  vertex.weights1 <- c(1,0,0,1,0)

  # test that Kfct outputs the correct list of statistics when individual=TRUE
  res1 <- SANTA:::Kfct(B, vertex.weights1, individual=TRUE)  
  checkEquals(names(res1), c("netK", "netAUK", "nodeK", "nodeAUK"))
  checkEquals(all(is.na(res1)==FALSE), TRUE)
  
  # test that Kfct outputs nodeK and nodeAUK=NA when individual=FALSE
  res2 <- SANTA:::Kfct(B, vertex.weights1, individual=FALSE)  
  checkEquals(names(res2), c("netK", "netAUK", "nodeK", "nodeAUK"))
  checkEquals(is.na(res2$nodeK) & is.na(res2$nodeAUK), TRUE)
  checkEquals(is.na(res2$netK[[1]]) & is.na(res2$netAUK[[1]]), FALSE)
  
  # test that Kfct outputs the correct netK, netAUK, nodeK and node AUK values for a square matrix B
  res3 <- SANTA:::Kfct(B, vertex.weights1, individual=TRUE)  
  checkEquals(round(res3$netK, 10), c(0.6, 0.6, 0.0, 0.0))
  checkEquals(round(res3$netAUK, 10), 0.3)
  checkEquals(dim(res3$nodeK), c(nrow(B), max(B)))
  checkEquals(round(res3$nodeAUK, 10), c(0.25, 0.05, 0.05, 0.35, 0.15))
  
  # test that Kfct correctly incorperates different vertex weights 
  vertex.weights2 <- c(0,1,0,0,1)
  res4 <- SANTA:::Kfct(B, vertex.weights2, individual=TRUE) 
  checkEquals(round(res4$netK, 10), c(0.6, 0.2, -0.4, 0.0))
  checkEquals(round(res4$netAUK, 10), 0.1)
  checkEquals(dim(res4$nodeK), c(nrow(B), max(B)))
  checkEquals(round(res4$nodeAUK, 10), c(-0.25, 0.05, -0.45, -0.15, 0.15))
  
  # test that Kfct produces an error if the dimensions of B is not length(vertex.weights1) x length(vertex.weights1)
  checkException(SANTA:::Kfct(B[1:4,], vertex.weights1, individual=TRUE), silent=TRUE)
  checkException(SANTA:::Kfct(B, vertex.weights1[1:4], individual=TRUE), silent=TRUE)
}
