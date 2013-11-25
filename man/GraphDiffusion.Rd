\name{GraphDiffusion}
\alias{GraphDiffusion}
\title{
Compute the distances between vertex pairs using a diffusion kernel-based method.
}
\description{
Using a diffusion kernel-based algorithm, compute the distance between vertex pairs in an undirected graph, with or without edge distances. This algorithm provides an alternative to the \code{shortest.paths} and \code{mfpt} measures of vertex pair distance.
}
\usage{
GraphDiffusion(g, v=V(g), edge.attr=NULL, beta=1, correct.neg=TRUE)
}
\arguments{
  \item{g}{
An igraph object, the graph to work on. 
}
  \item{v}{
An igraph object or numeric vector, the vertices from which each distance is calculated. 
}
  \item{edge.attr}{
A string, the names of the edge attribute to be used as distances along the edges. If left equal to \code{NULL}, then each edge is assumed to have a distance of 1. 
}
  \item{beta}{
A numeric value, the probability that the lazy random walk will take each of the edges emanating from a vertex. 
}
  \item{correct.neg}{
A logical value. If \code{TRUE}, then negative edge distances are set to 0. 
}
}
\details{
Diffusion across a graph follows a process similar to a random walk. This provides a method of measuring the distance between vertex pairs that does not simply take into account a single path (like the \code{shortest.paths} algorithm) but instead incorporates multiple paths. This function uses a diffusion kernel-based approach to calculate distances. The algorithm implemented is detailed in the referenced paper. 

The distance from vertex a to vertex a is always 0.
}
\value{
A numeric matrix, containing the diffusion kernel-derived vertex pair distance between each vertex in \code{v} and every vertex in \code{g}.
}
\references{
Kondor, R.I. and Lafferty, J. (2002). \emph{Diffusion Kernels on Graph and Other Discrete Structures}. Proc. Intl. Conf. Machine Learning. 
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk}
}
\seealso{
\code{\link[SANTA:GraphMFPT]{GraphMFPT}},
\code{\link[igraph]{shortest.paths}}
}
\examples{
# Create a graph and calculate the diffusion kernel-derived vertex pair distance matrix
g <- CreateGraph(n=6, type="barabasi")
GraphDiffusion(g)
plot(g, layout=layout.fruchterman.reingold)
}
