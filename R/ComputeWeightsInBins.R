ComputeWeightsInBins <- function(
    B, 
    weight, 
    nvertices
){
	# Use C call to compute the sum of the adjusted vertex weights for each vertex and each bin
	diam <- max(B)
	matrix(
		.Call("computeweights",
			as.integer(as.vector(B)),
			as.double(weight),
			as.integer(nvertices),
			as.integer(diam)
		),
	nvertices, diam)
}
