\name{GraphMFPT}
\alias{GraphMFPT}
\title{
Compute the distances between vertex pairs using a mean first-passage time-based method.
}
\description{
Using the mean first-passage time algorithm, compute the distance between vertex pairs in an undirected graph, with or without edge distances. This algorithm provides an alternative to the \code{shortest.paths} and \code{diffusion} measures of vertex pair distance.
}
\usage{
GraphMFPT(g, v = V(g), edge.attr = NULL, average.distances = TRUE)
}
\arguments{
  \item{g}{
An igraph object, the graph to work on. 
}
  \item{v}{
An igraph object or numeric vector, the vertices from which each distance is calculated.
}
  \item{edge.attr}{
A string, the name of the edge attribute to be used as distances along the edges. If left equal to NULL, then each edge is assumed to have a distance of 1.  
}
  \item{average.distances}{
A logical value. If \code{TRUE}, then the distance from vertex A to B and the distance from vertex B to A are averaged to give a single distance. Otherwise, two different distances may be returned.
}
}
\details{
The mean first-passage time from vertex a to vertex b is defined as the expected number of steps taken on a random walk from vertex a until the first arrival at vertex b. This provides a method of measuring the distance between pairs of vertices that does not simply take into account the distance along the shortest path, but rather incorporates how well the two vertices are connected across multiple paths.

The mean first-passage time from vertex a to vertex b is not necessarily the same as the mean first-passage time from vertex b to vertex a. If a symmetric distance matrix is required, reciprocal distances can be averaged to give a single value for each vertex pair.

If a vertex pair is unconnected, then the distance between the vertices is \code{Inf}.

The distance from vertex a to vertex a is always 0.
}
\value{
A numeric matrix, containing the mean first-passage time-derived vertex pair distance between each vertex in \code{v} and every vertex in \code{g}.
}
\references{
White, S. and Smyth, P. (2003). \emph{Algorithms for Estimating Relative Importance in Networks}. Technical Report UCI-ICS 04-25.
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk}
}
\seealso{
\code{\link[SANTA:GraphDiffusion]{GraphDiffusion}},
\code{\link[igraph]{shortest.paths}}
}
\examples{
# Create a graph and calculate the mean first-passage time-based vertex pair distance matrix
g <- CreateGraph(n=6, type="erdos.renyi")
GraphMFPT(g) 
plot(g, layout=layout.fruchterman.reingold)
}
