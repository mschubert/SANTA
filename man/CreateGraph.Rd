\name{CreateGraph}
\alias{CreateGraph}
\title{
Generate a simple undirected graph, with or without vertex weights and edge distances.
}
\description{
Generate a simple and undirected graph containing a specified number of vertices, with or without vertex weights and edge distances. Vertex weights added to the graph can be binary or continuous, randomly distributed or grouped into one or more clusters of high-weight vertices (hits). 
}
\usage{
CreateGraph(n=100, type="barabasi", m=2, p.or.m=(2*n-2*floor(sqrt(n))), vertex.weights=NULL, edge.distances=NULL, gen.vertex.weights=FALSE, nclusters=1, lambda=1, nlinks=ceiling(n/5), nhits=ceiling(n/10), binary.pheno=TRUE, dist.method="shortest.paths", mean.hit=1, sd.hit=0.05, mean.miss=0, sd.miss=0.05)
}
\arguments{
  \item{n}{
An integer value, the number of vertices to be included.
}
  \item{type}{
A string, the algorithm to be used to build the graph. Can either be \code{barabasi}, \code{erdos.renyi} or \code{grid}. 
}
  \item{m}{
An integer value, the number of edges to be added in each step of the \code{barabasi} algorithm, if used.
}
  \item{p.or.m}{
An integer value, the number of edges to be included within the graph, if the \code{erdos.renyi} algorithm is used. 
}
  \item{vertex.weights}{
A numeric vector, containing vertex weights to be added to the graph. The length of the vector must equal \code{n}. If the graph contains \code{nhits} hits, then the \code{nhits}-highest vertex weights are randomly assigned to the hit vertices. The remaining weights are randomly assigned to the remaining vertices (the misses). The vertex weights are added under a vertex attribute \code{pheno}.
}
  \item{edge.distances}{
A numeric vector, containing edge distances to be added to the graph. The length of the vector must equal the number of edges in the graph. The edge distances are randomly assigned to the edges under an edge attribute \code{distance}.
}
  \item{gen.vertex.weights}{
A logical constant, denoting whether vertex weights should be generated. If \code{TRUE} and \code{vertex.weights=NULL}, then vertex weights for hits and misses are generated and added under the vertex attribute \code{pheno}.
}
  \item{nclusters}{
An integer value, the number of clusters that the hits are grouped into. 
}
  \item{lambda}{
A numeric value, the strength of the hit clustering. If \code{lambda=0}, then the hits are randomly distributed. The greater the value of \code{lambda}, the greater the strength of hit clustering. Hits are only added if \code{gen.vertex.weights=TRUE} and/or \code{vertex.weights!=NULL}.
}
  \item{nlinks}{
An integer value, the number of linking edges to be added between clusters. 
}
  \item{nhits}{
An integer value, the number of hits to be added to the graph. Hits are only added if \code{gen.vertex.weights=TRUE} and/or \code{vertex.weights!=NULL}.
}
  \item{binary.pheno}{
A logical value, denoting whether generated vertex weights should be binary or continuous. If \code{TRUE}, then hit and miss vertices are assigned weights of 1 and 0 under the vertex attribute \code{pheno}. If \code{FALSE}, then weights are generated using two truncated normal distributions with parameters \code{mean.hit}, \code{sd.hit}, \code{mean.miss} and \code{sd.miss}. These continuous weights also range between 0 and 1. 
}
  \item{dist.method}{
A string, the method used to calculate the distance between each vertex and the start vertex. Can either be \code{shortest.paths}, \code{diffusion} or \code{mfpt}.
}
  \item{mean.hit}{
A numeric value, the mean parameter in the truncated normal distribution used to generate the vertex weight of the hits. 
}
  \item{sd.hit}{
A numeric value, the standard deviation parameter in the truncated normal distribution used to generate the vertex weight of the hits. 
}
  \item{mean.miss}{
A numeric value, the mean parameter in the truncated normal distribution used to generate the vertex weight of the misses. 
}
  \item{sd.miss}{
A numeric value, the standard deviation parameter in the truncated normal distribution used to generate the vertex weight of the misses. 
}
}
\details{
This function creates simple, undirected graphs, with or without associated vertex weights and edge distances.

1 of 3 different algorithms can be used to build the graph. The \code{barabasi} algorithm builds a scale-free graph through preferencial attachment of vertices. The \code{erdos.renyi} algorithm builds a random graph with a fixed number of edges. The \code{grid} algorithm builds a graph with a grid-like structure. The graph produced is always simple, meaning that it contains no loops and no multiple edges. 

Vertex weights can be added to the graph as a vertex attribute under the name \code{pheno}. Larger vertex weights indicate that the vertex is more strongly associated with a certain phenotype or function. The number of vertex weights input must equal the number of vertices present within the graph. Larger vertex weights are randomly assigned to the hit vertices and smaller vertex weights are randomly assigned to the misses. 

Edge distances can be added to the graph as an edge attribute under the name \code{distance}. Smaller edge distances indicate that two vertices are more strongly connected. The number of edge distances input must equal the number of edges present in the graph. The edge distances are randomly assigned to the different edges. 


If \code{vertex.weights!=NULL} and/or \code{gen.vertex.weights=TRUE}, then hits are added across the graph. \code{lambda} shapes the probability distribution under which the hits are added. If \code{lambda=0}, then hits are added to each vertex with equal probability. If \code{lambda>0}, then the hits become clustered around one or more random vertices on the graph. The greater the value of \code{lambda}, the greater the strength of the clustering. The probability of vertex i being a hit is proportional to

\deqn{P[i] \sim lambda*exp^{-lambda * d[start,i]}}{P[i] ~ lambda*exp^{-lambda * d[start,i]}}

where \eqn{d[start,i]} is the distance between the start vertex and vertex i, according to the shortest paths distance measure.  

If \code{binary.pheno=TRUE}, then hits and misses are assigned weights of 1 and 0 respectively. Otherwise, vertex weights ranging between 0 and 1 are generated using 2 truncated normal probability distributions - one distribution for the hits and another for the misses. The shape of these distributions are set using the \code{mean.hit}, \code{sd.hit}, \code{mean.miss} and \code{sd.miss} parameters. 

Multiple clusters of hits can be added to the graph. These clusters are created by adding hits to multiple graphs. These graphs are then connected through the addition of \code{nlinks} linking edges. 

}
\value{
A igraph object. If vertex weights are added, then whether the vertex is a hit or a miss is specified under the vertex attribute \code{hits}, the cluster from which the hit is from under the vertex attribute \code{hits.cluster} and the vertex weight under the vertex attribute \code{pheno}. A vertex attribute named \code{color} is also added to allow for the identification of hits when the graph is plotted using \code{plot}.
}
\references{
Barabasi, A.L. and Albert, R. (1999). \emph{Emergence of scaling in random networks}. Science 286: 509-512. 

Erdos, P. and Renyi, A. (1959). \emph{On random graphs}. Publicationes Mathematicae 6: 290-297. 
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk}
}
\seealso{
\code{\link[SANTA:CreateGrid]{CreateGrid}},
\code{\link[igraph]{barabasi.game}},
\code{\link[igraph]{erdos.renyi.game}}
}
\examples{
# Create a 30-vertex graph with a scale-free structure and no vertex weights
g1 <- CreateGraph(n=30, type="barabasi", gen.vertex.weights=FALSE)
plot(g1, layout=layout.fruchterman.reingold)

# Create a 30-vertex graph with a random structure and 2 clusters of high-weight vertices
g2 <- CreateGraph(n=30, type="erdos.renyi", gen.vertex.weights=TRUE, nhits=8, nclusters=2,
lambda=10, binary.pheno=TRUE)
plot(g2, layout=layout.fruchterman.reingold)

# Create a 36-vertex grid-shaped graph with 6 randomly-distributed high-weight vertices
g3 <- CreateGraph(n=36, type="grid", gen.vertex.weights=TRUE, nhits=6, lambda=0, 
binary.pheno=FALSE)
plot(g3, layout=layout.fruchterman.reingold)
}
