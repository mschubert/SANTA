\name{SpreadHits}
\alias{SpreadHits}
\title{
Spread hits across a network in one or more clusters
}
\description{
Spread hits across a network in one or more clusters using an exponential probability distribution related to the distance of each vertex from a seed vertex. 
}
\usage{
SpreadHits(g, h, clusters=1, distance.cutoff=3, lambda=1, 
dist.method=c("shortest.paths", "diffusion", "mfpt"), 
edge.attr=NULL, hit.color="red", D=NULL, attempts=1000, verbose=TRUE)
}
\arguments{
  \item{g}{
An igraph object. The network to work on.
}
  \item{h}{
An integer value. The number of hits in each cluster.
}
  \item{clusters}{
An integer value. The number of clusters to add the network. 
}
  \item{distance.cutoff}{
An integer value. The minimum distance between the seed vertex of each cluster. 
}
  \item{lambda}{
A numeric value. The strength of hit clustering. If \code{lambda=0}, then the hits are randomly distributed. If \code{lambda>0}, then the hits are clustered. The greater the value of \code{lambda}, the greater the strength of the hit clustering. 
}
  \item{dist.method}{
A string. The method used to compute the distance between each vertex and the seed vertex. 
}
  \item{edge.attr}{
A string. The name of the edge attribute to be used as distances along the edges. If left equal to \code{NULL}, then each edge is assumed to have a distance of 1. Smaller edge distances denote stronger interactions between vertex pairs.
}
  \item{hit.color}{
A string. The colour that the hits will be if the network is plotted using \code{plot()}.
}
  \item{D}{
A pre-computed numeric distance matrix. If \code{NULL}, then \code{D} is computed within the function.
}
  \item{attempts}{
An integer value. The number of attempts made at each stage of the placement algorithm before the function terminates. 
}
  \item{verbose}{
A logical value. If \code{TRUE}, messages about the progress of the function are displayed. 
}
}
\details{
\code{SpreadHits} can be used to add hits to a network without hits, or replace hits on a network with hits. The probability of a vertex being labelled as a hit depends on its distance from a randomly chosen starting vertex. The value of \code{lambda} denotes the shape of the probability distribution used to spread the hits. The greater the value of \code{lambda}, the greater the strength of hit clustering. The probability of vertex i being a hit is proportional to

\deqn{P[i] \sim lambda^{-d[start,i]}}

where \eqn{d[start,i]} is the distance between the start vertex and vertex i.  

Hits can be added in one or more clusters. If multiple clusters are being added, then the function first identifies an equal number of seed vertices at least \code{distance.cutoff} distance apart in the network. If the function cannot identify vertices that satisfy this condition within the specified number of attempts, then it returns \code{NULL}. If the function is able to find suitable seed vertices, then it proceeds to label vertices positioned around each seed vertex as hits, using a probability function shaped by \code{lambda}.

Whether a vertex is a hit or a miss is given in the vertex attribute \code{hits}. If the vertex is a hit, then it has the score 1. If it is a miss, then it has the score 0. A color can also be chosen to highlight the hits when the network is plotted. Misses are automatically coloured grey. 
}
\value{
A modified version of the input \code{igraph} object. Whether a vertex is a hit or miss is given under the vertex attribute \code{hits}. 
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk}
}
\examples{
# create a network and add 1 cluster of 3 hits
g1 <- erdos.renyi.game(30, p.or.m=0.1, directed=FALSE)
g1 <- SpreadHits(g1, h=3, clusters=1, lambda=10)
plot(g1, layout=layout.fruchterman.reingold)

# create a network and add 2 clusters of 3 hits
g2 <- erdos.renyi.game(30, p.or.m=0.1, directed=FALSE)
g2 <- SpreadHits(g2, h=3, clusters=3, distance.cutoff=2, lambda=0)
plot(g2, layout=layout.fruchterman.reingold)
}
