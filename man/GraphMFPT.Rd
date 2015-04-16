\name{GraphMFPT}
\alias{GraphMFPT}
\title{
Compute mean first-passage time-based distance matrix
}
\description{
Using the mean first-passage time method, compute the distances between vertex pairs in an undirected graph, with or without edge weights.
}
\usage{
GraphMFPT(g, v=V(g), edge.attr.weight=NULL, average.distances=TRUE)
}
\arguments{
  \item{g}{
\code{igraph} object, the graph to work on. 
}
  \item{v}{
\code{igraph} object or numeric vector, the vertices from which each distance is calculated.
}
  \item{edge.attr.weight}{
String, the name of the edge attribute to be used as weights along the edges. Greater weights indicate a stronger interaction between the two genes (this is the opposite to edge distances, where smaller distances indicate stronger interactions). If \code{NULL} then each edge is assumed to have a weight of 1.   
}
  \item{average.distances}{
Logical, if \code{TRUE} then the distance from vertex A to B and the distance from vertex B to A are averaged to give a single distance. Otherwise, two different distances may be returned.
}
}
\details{
The mean first-passage time from vertex A to vertex B is defined as the expected number of steps taken on a random walk emanating from vertex A until the first arrival at vertex B. This provides a method of measuring the distance between pairs of vertices that does not simply take into account the distance along the shortest path, but rather incorporates how well the two vertices are connected across multiple paths.

The mean first-passage time from vertex A to vertex B is not necessarily the same as the mean first-passage time from vertex B to vertex A. If a symmetric distance matrix is required, reciprocal distances can be averaged to give a single value for each vertex pair.

If a vertex pair is unconnected, then the distance between the vertices is \code{Inf}.

The distance from vertex A to vertex A is always 0.
}
\value{
Numeric matrix, containing the mean first-passage time-derived vertex pair distance between each vertex in \code{v} and every vertex in \code{g}.
}
\references{
White, S. and Smyth, P. (2003). \emph{Algorithms for Estimating Relative Importance in Networks}. Technical Report UCI-ICS 04-25.
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk}
}
\seealso{
\code{\link[SANTA:GraphDiffusion]{GraphDiffusion}},
\code{\link[igraph]{shortest.paths}}
}
\examples{
# create a  and compute the mean first-passage time-based vertex pair distance matrix
g <- erdos.renyi.game(6, p.or.m=0.5, directed=FALSE)
GraphMFPT(g) 
plot(g, layout=layout.fruchterman.reingold)
}
