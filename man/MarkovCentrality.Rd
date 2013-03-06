\name{MarkovCentrality}
\alias{MarkovCentrality}
\title{
Compute the Markov centrality score for each vertex in a graph.
}
\description{
The Markov centrality score uses the concept of a random walk through the graph to calculate the centrality of each vertex. The method uses the mean first-passage time from every vertex to every other vertex to produce a score for each vertex. These scores can be used as a ranking of centrality within the graph.
}
\usage{
MarkovCentrality(g, edge.attr = NULL)
}
\arguments{
  \item{g}{
An igraph object, the graph to work on.
}
  \item{edge.attr}{
A string, the name of the edge attribute to be used as distances along the edges. If left equal to \code{NULL}, then each edge is assumed to have a distance of 1. 
}
}
\details{
The mean first-passage time can be used as a measure of how closely each vertex is connected to every other vertex in a graph. The mean first-passage time from vertex a to vertex b is the mean number of steps a random walk emanating from vertex a takes to reach vertex b for the first time. Random walks are more likely to reach well-connected vertices quicker and therefore this method can be used to measure distance. 

In order to calculate the Markov centrality of each vertex in a graph, the inverse of the mean distance between each vertex and every other vertex is taken. Vertices with smaller average distances to all other vertices have higher Markov centrality scores, indicating that they occupy a more central position within the graph. These values can be used to rank the vertices. 
}
\value{
A numeric vector, containing the Markov centrality score of each vertex. 
}
\references{
White, S. and Smyth, P. (2003). \emph{Algorithms for Estimating Relative Importance in Networks}. Technical Report UCI-ICS 04-25.
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk}
}
\seealso{
\code{\link[SANTA:GraphMFPT]{GraphMFPT}}
}
\examples{
# Create a graph and compute the Markov centrality score for each vertex
g <- CreateGraph(n=6, type="barabasi")
MarkovCentrality(g)
plot(g, layout=layout.fruchterman.reingold)
}
