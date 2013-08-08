plot.Knet <- function(
  x, 
  sequential=FALSE,
  ...
){
  # plot the Knet score

  if (all(is.na(x$K.perm), is.na(x$AUK.perm), is.na(x$K.quan))){
    # if no permutations present, plot just the K function line
    plot(x$K.obs, main="Knet function", type="l", lwd=3, col="red", xlab="distance on graph", ylab="Knet")
  } else {
    # if permutations present, produce 2 plots comparing the observed and permutation K function 
    if (!sequential) par(mfrow=c(1, 2))

    # extract Scores
    K.obs    <- x$K.obs
    K.perm   <- x$K.perm
    K.quan   <- x$K.quan
    AUK.obs  <- x$AUK.obs
    AUK.perm <- x$AUK.perm
    g.max	 <- ncol(K.quan)-1
    
    # plot 1
    if (nrow(K.quan) == 5) { 
      lty <- c("dotted", "solid", "solid", "solid", "dotted") 
      lwd <- c(1, 1, 2, 1, 1)
    } else {
      lty	<- rep("solid", nrow(K.quan))
      lwd <- rep(1, nrow(K.quan))
    }
    plot(NA, ylim=range(K.obs, K.quan), xlim=c(0, g.max), xlab="Distance in graph", ylab="Knet", main="Knet-function", ...)
    polygon(x=c(0:g.max, g.max:0), y=c(K.quan["0%",], rev(K.quan["100%",])), col="lightyellow", border=NA)
    for (i in 1:nrow(K.quan)) lines(0:(ncol(K.quan)-1), K.quan[i,], col="grey", lty=lty[i], lwd=lwd[i])
    lines(x=0:g.max, y=K.obs, col="red", lwd="3")
	
    # plot 2
    xlim.values <- range(c(AUK.perm, AUK.obs)) + 1 * (range(c(AUK.perm, AUK.obs)) - mean(range(c(AUK.perm, AUK.obs)))) 
    hist(AUK.perm, xlim=xlim.values, col="grey", border="grey", xlab="Area under K-curve", main="AUK: observed v. permuted")
    abline(v=AUK.obs, lwd=2, col="red")
	}
}
