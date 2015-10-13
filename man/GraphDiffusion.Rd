\name{GraphDiffusion}
\alias{GraphDiffusion}
\title{
Compute diffusion kernel-based distance matrix
}
\description{
Using a diffusion kernel-based algorithm, compute the distance between vertex pairs in an undirected network, with or without edge weights. This algorithm provides an alternative to the \code{shortest.paths} and \code{mfpt} measures of vertex pair distance.
}
\usage{
GraphDiffusion(g, v=V(g), edge.attr.weight=NULL, beta=1, correct.neg=TRUE)
}
\arguments{
  \item{g}{
\code{igraph} object, the network to work on. 
}
  \item{v}{
\code{igraph} object or numeric vector, the vertices from which each distance is calculated. 
}
  \item{edge.attr.weight}{
String, the name of the edge attribute to be used as weights along the edges. Greater weights indicate a stronger interaction between the two genes (this is the opposite to edge distances, where smaller distances indicate stronger interactions). If \code{NULL}, then each edge is assumed to have a weight of \code{1}.  
}
  \item{beta}{
Numeric value, the probability that the diffusion process will take an edge emanating from a vertex. 
}
  \item{correct.neg}{
Logical, if \code{TRUE} then negative edge distances are set to \code{0}. 
}
}
\details{
Diffusion across a network follows a process similar to a random walk. This provides a method of measuring the distance between vertex pairs that does not simply take into account a single path (like the \code{shortest.paths} algorithm) but instead incorporates multiple paths. This function uses a diffusion kernel-based approach to compute distances. The algorithm implemented is detailed in the referenced paper. 

The distance from vertex A to vertex A is always \code{0}.
}
\value{
Numeric matrix, containing the diffusion kernel-based vertex pair distances between each vertex in \code{v} and every vertex in \code{g}.
}
\references{
Kondor, R.I. and Lafferty, J. (2002). \emph{Diffusion Kernels on Graph and Other Discrete Structures}. Proc. Intl. Conf. Machine Learning. 
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk}
}
\seealso{
\code{\link[SANTA:GraphMFPT]{GraphMFPT}},
\code{\link[igraph]{shortest.paths}}
}
\examples{
# create a network and computes the diffusion kernel-derived vertex pair distance matrix
g <- barabasi.game(6, directed=FALSE)
GraphDiffusion(g)
plot(g, layout=layout.fruchterman.reingold)
}
