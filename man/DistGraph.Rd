\name{DistGraph}
\alias{DistGraph}
\title{
Compute the vertex pair distance matrix for a graph
}
\description{
Compute the distances between pairs of vertices in a graph, using a shortest path, diffusion kernel, or mean first-passage time-based measure. 
}
\usage{
DistGraph(g, v=V(g), edge.attr=NULL, dist.method=c("shortest.paths", "diffusion", "mfpt"), 
correct.inf=TRUE, correct.factor=1, verbose=TRUE)
}
\arguments{
  \item{g}{
An igraph object. The graph on which to work.
}
  \item{v}{
An igraph object or numeric vector. The vertices from which each distance is calculated. 
}
  \item{edge.attr}{
A string. The name of the edge attribute to be used as distances along the edges. If left equal to \code{NULL}, then each edge is assumed to have a distance of 1. Smaller edge distances denote stronger interactions between vertex pairs.
}
  \item{dist.method}{
A string. The method used to compute the distance between each vertex pair.
}
  \item{correct.inf}{
A logical value. If \code{TRUE}, then infinite vertex pair distances are replaced with distances equal to the maximum distance measured across the network, multiplied by \code{correct.factor}. If \code{FALSE}, then these distances are returned as \code{Inf}.
}
\item{correct.factor}{
A numeric value. If the graph contains unconnected vertices, then the distance between these vertices is set as the maximum distance between the connected vertices multiplied by \code{correct.factor}.
}
  \item{verbose}{
A logical value. If \code{TRUE}, messages about the progress of the function are displayed. 
}
}
\details{
This function computes a distance matrix for a graph. Different methods can be used to calculate the distance between each pair of vertices. If a set of vertices is specified, a smaller distance matrix containing only vertices corresponding to the vertices is returned.

Descriptions of how the shortest paths (\code{shortest.paths}), diffusion kernel-based (\code{GraphDiffusion}) and mean first-passage time (\code{GraphMFPT}) distance measures work are given in their respection function descriptions.
}
\value{
A numeric matrix, containing the distances between vertex pairs. 
}
\references{
Kondor, R.I. and Lafferty, J. (2002). \emph{Diffusion Kernels on Graph and Other Discrete Structures}. Proc. Intl. Conf. Machine Learning. 

White, S. and Smyth, P. (2003). \emph{Algorithms for Estimating Relative Importance in Networks}. Technical Report UCI-ICS 04-25.
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk}
}
\seealso{
\code{\link[igraph]{shortest.paths}}
\code{\link[SANTA:GraphDiffusion]{GraphDiffusion}},
\code{\link[SANTA:GraphMFPT]{GraphMFPT}},
}
\examples{
# create a graph and compute the distance matrix using the shortest paths measure
g1 <- barabasi.game(6, directed=FALSE)
DistGraph(g1, dist.method="shortest.paths")
plot(g1, layout=layout.fruchterman.reingold)

# create a graph, assign edge distances and compute the distance matrix using the 
# diffusion kernel-based measure
g2 <- erdos.renyi.game(6, p.or.m=0.5, directed=FALSE)
g2 <- set.edge.attribute(g2, name="distance", value=runif(ecount(g2)))
DistGraph(g2, dist.method="diffusion", edge.attr="distance")
plot(g2, layout=layout.fruchterman.reingold)
}
