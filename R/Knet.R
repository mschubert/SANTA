Knet <- function(
  g, 
  nperm          =100, 
  dist.method    ="shortest.paths", 
  vertex.attr    ="pheno",
  edge.attr      ="distance",
  correct.factor =1.0,
  nsteps         =1000, 
  prob           =c(0, 0.05, 0.5, 0.95, 1),
  parallel       =NULL
){
  # calculate the Knet function for a graph, along with permutations if required
  
  # check that vertex weights and edge distances are present and suitable. Convert if neccessary
  g <- CheckAttributes(g, vertex.attr, edge.attr)
  
  # the results are saved in a list with an entry for each vertex.attr. Done even if only 1 vertex.attr supplied
  tmp        <- vector("list", 8)
  names(tmp) <- c("K.obs", "AUK.obs", "K.perm", "AUK.perm", "K.quan", "nodeK", "nodeAUK", "pval")
  class(tmp) <- "Knet"
  res        <- rep(list(tmp), length(vertex.attr)) 
  names(res) <- vertex.attr
  
  # compute the vertex pair distances (D)
  message("computing graph distance matrix...", appendLF=FALSE)
  D <- DistGraph(g=g, edge.attr=edge.attr, dist.method=dist.method, correct.inf=TRUE, correct.factor=correct.factor) 
  message(" done")
  
  # compute which bin each vertex pair distance falls into (B)
  message("computing graph distance bins...", appendLF=FALSE)
  B <- BinGraph(D=D, dist.method=dist.method, nsteps=nsteps) 
  message(" done")
  
  # if parallel computing is to be used, set up a cluster
  if (!is.null(parallel)) {
    if (!IsWholeNumber(parallel)) stop("Argument 'parallel' is not a whole number")
    message(paste("setting up cluster of size ", parallel, "...", sep=""), appendLF=FALSE)
    cl <- makeCluster(parallel, type="SOCK", verbose=FALSE)		
    clusterEvalQ(cl, library(SANTA))	
    message(" done")
  }
  
  # run the function for each vertex attribute in vertex.attr
  for (attr in vertex.attr) {
    if (length(vertex.attr)==1) {
      message(paste("computing the clustering of the '", attr, "' weights using ", nperm, " permutations...", sep=""), appendLF=FALSE)	
    } else {
      message(paste("computing the clustering of the '", attr, "' weights (", which(vertex.attr==attr), "/", length(vertex.attr) ,") using ", nperm, " permutations...", sep=""), appendLF=FALSE)	
    }
    
    # extract vertex weights 
    vertex.weights <- get.vertex.attribute(g, attr)
    
    # calculate the observed netK and netAUK
    K                   <- Kfct(B, vertex.weights, individual=TRUE) 
    res[[attr]]$K.obs 	<- K$netK
    res[[attr]]$AUK.obs <- K$netAUK
    res[[attr]]$nodeK 	<- K$nodeK
    res[[attr]]$nodeAUK <- K$nodeAUK
    
    # if specified, run the Knet function on permutations of the graph
    if (!is.null(nperm) & nperm>0) {
      if (is.null(parallel)) {
        # run permutations without parallel computing
        res[[attr]]$K.perm 	<- sapply(1:nperm, function(i) Kfct(B, sample(vertex.weights), individual=FALSE)$netK) 
      } else {
        # run permutations with parallel computing
        res[[attr]]$K.perm 	<- parSapply(cl, 1:nperm, function(i) Kfct(B, sample(vertex.weights), individual=FALSE)$netK)
      }
      
      # calculate the quantiles, AUK and p-values (through the z-score) for the permutations
      res[[attr]]$K.quan    <- apply(res[[attr]]$K.perm, 1, function(x) quantile(x, prob=prob))
      res[[attr]]$AUK.perm  <- apply(res[[attr]]$K.perm, 2, function(x) sum(x) / length(x))
      res[[attr]]$pval      <- pnorm((res[[attr]]$AUK.obs - mean(res[[attr]]$AUK.perm)) / sd(res[[attr]]$AUK.perm), lower.tail=FALSE) 
    } else {
      # if no permutations are run, permutation-related statistics are returns equal to NA
      res[[attr]]$K.perm    <- NA
      res[[attr]]$AUK.perm  <- NA
      res[[attr]]$K.quan    <- NA
      res[[attr]]$pval      <- NA	
    }
    
    message(" done")
  }
  
  # close the cluster
  if (!is.null(parallel)) stopCluster(cl)	
  
  # if only one vertex attribute is input, don't return a list of lists of results
  if (length(vertex.attr)==1) res <- res[[1]]
  
  # output
  res
}
