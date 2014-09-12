\name{GraphDiffusion}
\alias{GraphDiffusion}
\title{
Compute the distances between vertex pairs using a diffusion kernel-based method
}
\description{
Using a diffusion kernel-based algorithm, compute the distance between vertex pairs in an undirected network, with or without edge weights. This algorithm provides an alternative to the \code{shortest.paths} and \code{mfpt} measures of vertex pair distance.
}
\usage{
GraphDiffusion(g, v=V(g), edge.attr.weight=NULL, beta=1, correct.neg=TRUE)
}
\arguments{
  \item{g}{
An igraph object. The network to work on. 
}
  \item{v}{
An igraph object or numeric vector. The vertices from which each distance is calculated. 
}
  \item{edge.attr.weight}{
A string. The name of the edge attribute to be used as weights along the edges. Greater weights indicate a stronger interaction between the two genes (this is the opposite to edge distances, where smaller distances indicate stronger interactions). If left equal to NULL, then each edge is assumed to have a weight of 1.  
}
  \item{beta}{
A numeric value. The probability that the diffusion process will take am edges emanatin from a vertex. 
}
  \item{correct.neg}{
A logical value. If \code{TRUE}, then negative edge distances are set to 0. 
}
}
\details{
Diffusion across a network follows a process similar to a random walk. This provides a method of measuring the distance between vertex pairs that does not simply take into account a single path (like the \code{shortest.paths} algorithm) but instead incorporates multiple paths. This function uses a diffusion kernel-based approach to compute distances. The algorithm implemented is detailed in the referenced paper. 

The distance from vertex A to vertex A is always 0.
}
\value{
A numeric matrix, containing the diffusion kernel-based vertex pair distances between each vertex in \code{v} and every vertex in \code{g}.
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
