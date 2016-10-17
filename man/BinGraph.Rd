\name{BinGraph}
\alias{BinGraph}
\title{
Compute a distance bin matrix from a distance matrix
}
\description{
In order for the \code{Knet} or \code{Knode} functions to be run on a network, it is necessary to place the raw distances between each vertex pair into discreet distance bins, as the functions are unable to handle continuous distributions of distances. This is done automatically within the \code{Knet} and \code{Knode} functions, but can also be done separately using \code{BinGraph}. 
}
\usage{
BinGraph(D, nsteps=1000, equal.bin.fill=TRUE, verbose=TRUE)
}
\arguments{
  \item{D}{
Numeric matrix, a distance matrix output by \code{DistGraph}.
}
  \item{nsteps}{
Integer value, the desired number of bins across which the distances are to be split. If there are too few unique distances to fill each bin, then fewer bins are returned. 
}
  \item{equal.bin.fill}{
Logical, if \code{TRUE} then the function attempts to fill each bin with an equal number of vertex pairs.
}
  \item{verbose}{
Logical, if \code{TRUE} messages about the progress of the function are displayed. 
}
}
\details{
In order for the \code{Knet} or \code{Knode} functions to be run, the vertex pair distances (as computed by \code{DistGraph}) but be split into bins. This is done as part of the \code{Knet} or \code{Knode}. However, this step is often slow for large networks and therefore the \code{BinGraph} function is provided separately, in order to avoid repeat computation. 

Each vertex pair is placed into a bin, either ranging from 1 to \code{nsteps}, or from 1 to the number of unique distances. 

If \code{equal.bin.fill} is \code{FALSE}, then the bin each vertex pair is placed into is directly proportional to the largest vertex pair distance. For example, if the distance between the pair is 25\% of the largest distance and \code{nsteps} equals \code{100}, then the vertex pair will be placed into bin \code{25}. However, this can create problems when there are a small number of edges with especially large distances, as this can result in the majority of vertex pairs being placed into a small number of bins. This can reduce the effectiveness of the \code{Knet} and \code{Knode} functions. Therefore, when \code{equal.bin.fill} is \code{TRUE}, the function attempts to fill each bin with an equal number of vertex pairs. If there are a large number of tied distances, then the bins may not be filled equally. 
}
\value{
Integer matrix with the same dimensions as \code{D}.
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk}
}
\seealso{
\code{\link[SANTA:DistGraph]{DistGraph}}
}
\examples{
# create network and calculate the distance matrix using the shortest paths measure
g1 <- barabasi.game(6, directed=FALSE)
plot(g1, layout=layout.fruchterman.reingold)
D1 <- DistGraph(g1, dist.method="shortest.paths")
# place the distances into distance bins
BinGraph(D1, nsteps=100)

# create network and calculate the distance matrix using diffusion kernel-based measure
g2 <- erdos.renyi.game(6, p.or.m=0.5, directed=FALSE)
g2 <- set.edge.attribute(g2, name="distance", value=runif(ecount(g2)))
plot(g2, layout=layout.fruchterman.reingold)
# place the distances into distance bins
D2 <- DistGraph(g2, dist.method="diffusion", edge.attr="distance")
BinGraph(D2, nsteps=100)
}
