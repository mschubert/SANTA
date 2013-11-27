\name{BinGraph}
\alias{BinGraph}
\title{
Compute a distance bin matrix from a distance matrix.
}
\description{
In order for \code{Knet} or \code{Knode} to be run on a graph, it is necessary to place the raw distances between each node pair into distance bins. The number and size of the bins that \code{BinGraph} places each distance in to depends on the maximum distance measured, the method used to measure the distance and the number of steps desired. 
}
\usage{
BinGraph(D, dist.method, nsteps)
}
\arguments{
  \item{D}{
A numerical matrix. The distance matrix output by \code{DistGraph}.
}
  \item{dist.method}{
A string. The method used to compute the distance between each vertex pair. Can be either \code{shortest.paths}, \code{diffusion} or \code{mfpt}.
}
  \item{nsteps}{
An integer value. Unless the distance method used is \code{shortest.paths} and every value in \code{D} is a whole number, \code{nsteps} determines the number of bins the distances are placed into.
}
}
\details{
If the graph for which \code{D} was computed contains only integer edge weights and if the \code{shortest.paths} distance measure was used to compute vertex pair distances, then it will be the case that D consists only of integers. If this is the case, then the values of \code{D} can be placed into bins of size 1, as smaller bins sizes would provide no additional information.
}
\value{
A numerical matrix with the same dimensions as \code{D}. The values represent the distance bin that each vertex pair is placed into. 
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk}
}
\seealso{
\code{\link[SANTA:DistGraph]{DistGraph}},
}
\examples{
# Create a graph and calculate the distance matrix using the shortest paths measure
# Place the distances into distance bins. 
g1 <- CreateGraph(n=6, type="barabasi")
plot(g1, layout=layout.fruchterman.reingold)
D1 <- DistGraph(g1, dist.method="shortest.paths")
BinGraph(D1, dist.method="shortest.paths", nsteps=100)

# Create a graph, assign edge distances and calculate the distance matrix using the 
# diffusion kernel-based measure. Place the distances into distance bins
g2 <- CreateGraph(n=6, type="erdos.renyi")
g2 <- set.edge.attribute(g2, name="distance", value=runif(ecount(g2)))
plot(g2, layout=layout.fruchterman.reingold)
D2 <- DistGraph(g2, dist.method="diffusion", edge.attr="distance")
BinGraph(D2, dist.method="diffusion", nsteps=100)
}
