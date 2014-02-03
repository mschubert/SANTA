\name{SpreadHits}
\alias{SpreadHits}
\title{
Spread hits across a graph.
}
\description{
Spread hits across a graph using an exponential probability distribution related to the distance of each vertex from the start vertex. 
}
\usage{
SpreadHits(g, h = 10, lambda = 1, dist.method = "shortest.paths", 
edge.attr = NULL, start.vertex = NULL, hit.color = "red", D = NULL)
}
\arguments{
  \item{g}{
An igraph object, the graph to work on.
}
  \item{h}{
An integer value, the number of hits to be added to the graph.
}
  \item{lambda}{
A numeric value, the strength of hit clustering. If \code{lambda=0}, then the hits are randomly distributed. If \code{lambda>0}, then the hits are clustered. The greater the value of \code{lambda}, the greater the strength of the hit clustering. 
}
  \item{dist.method}{
A string, the method used to calculate the distance between each vertex and the start vertex. Can either be \code{shortest.paths}, \code{diffusion} or \code{mfpt}. 
}
  \item{edge.attr}{
A string, containing the name of the edge attribute to be used as distances along the edges. If left equal to \code{NULL}, then each edge is assumed to have a distance of 1.
}
  \item{start.vertex}{
An igraph object, containing the start vertex. 
}
  \item{hit.color}{
A string, the colour that hits will be if the graph is plotted using \code{plot()}.
}
  \item{D}{
A pre-calculated numeric distance matrix. 
}
}
\details{
\code{SpreadHits} can be used to add hits to a graph without hits, or replace hits on a graph with hits. The probability of a vertex being labelled as a hit depends on its ranked distance from a randomly chosen or specified starting vertex. The value of \code{lambda} denotes the shape of the probability distribution used to spread the hits. The greater the value of \code{lambda}, the greater the strength of hit clustering. 

Hits are denoted as 1 under the vertex attributes \code{hits} and \code{pheno}, while misses are denoted as 0. A color can also be chosen to highlight the hits when the graph is plotted. Misses are automatically coloured grey. 
}
\value{
A modified version of the input igraph object. Whether a vertex is a hit or miss is given under the vertex attributes \code{hits} and \code{pheno}. 
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk}
}
\seealso{
\code{\link[SANTA:CreateGraph]{CreateGraph}}
}
\examples{
# Create a graph and add 5 clustered hits
g1 <- CreateGraph(n=30, gen.vertex.weights=FALSE)
g1 <- SpreadHits(g1, h=5, lambda=10)
plot(g1, layout=layout.fruchterman.reingold)

# Create a graph and add 5 unclustered hits
g2 <- CreateGraph(n=30, gen.vertex.weights=FALSE)
g2 <- SpreadHits(g2, h=5, lambda=0)
plot(g2, layout=layout.fruchterman.reingold)
}
