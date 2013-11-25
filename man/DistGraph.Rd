\name{DistGraph}
\alias{DistGraph}
\title{
Compute the vertex pair distance matrix of a graph.
}
\description{
Compute the distances between pairs of vertices in a graph, using a shortest path, diffusion kernel, or mean first-passage time-based measure. 
}
\usage{
DistGraph(g, v = V(g), edge.attr = NULL,
dist.method = c("shortest.paths", "diffusion", "mfpt"), correct.inf = TRUE, correct.factor = 1)
}
\arguments{
  \item{g}{
An igraph object. The graph on which to work.
}
  \item{v}{
An igraph object or numeric vector. The vertices from which each distance is calculated. 
}
  \item{edge.attr}{
A string. The name of the edge attribute to be used as distances along the edges. If left equal to \code{NULL}, then each edge is assumed to have a distance of 1. 
}
  \item{dist.method}{
A string. The method used to calculate the distance between each vertex pair. Can either be \code{shortest.paths}, \code{diffusion} or \code{mfpt}.
}
  \item{correct.inf}{
A logical value. If \code{TRUE}, then infinite vertex pair distances are replaced with distances equal to the maximum distance measured across the network, multiplied by \code{correct.factor}. If \code{FALSE}, then these distances are returned as \code{Inf}.
}
  \item{correct.factor}{
A numerical value. The factor by which the maximum measured distance is multiplied. 
}
}
\details{
This function calculates a distance matrix for a graph. Different methods can be used to calculate the distance between each pair of vertices. By specifying a set of vertices, a smaller distance matrix containing only the rows of the input vertices can be returned.
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
\code{\link[SANTA:GraphDiffusion]{GraphDiffusion}},
\code{\link[SANTA:GraphMFPT]{GraphMFPT}},
\code{\link[igraph]{shortest.paths}}
}
\examples{
# Create a graph and calculate the distance matrix using the shortest paths measure
g1 <- CreateGraph(n=6, type="barabasi")
DistGraph(g1, dist.method="shortest.paths")
plot(g1, layout=layout.fruchterman.reingold)

# Create a graph, assign edge distances and calculate the distance matrix using the 
# diffusion kernel-based measure
g2 <- CreateGraph(n=6, type="erdos.renyi")
g2 <- set.edge.attribute(g2, name="distance", value=runif(ecount(g2)))
DistGraph(g2, dist.method="diffusion", edge.attr="distance")
plot(g2, layout=layout.fruchterman.reingold)
}
