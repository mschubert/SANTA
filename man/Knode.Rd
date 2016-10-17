\name{Knode}
\alias{Knode}
\title{
Rank vertices by their strength of association with high-weight vertices 
}
\description{
Rank vertices by their strength of association with high-weight vertices using a modified version of Ripley's K-statistic. Vertex weights can either be binary or positive and continuous.
}
\usage{
Knode(g, dist.method=c("shortest.paths", "diffusion", "mfpt"), vertex.attr="pheno", 
edge.attr=NULL, correct.factor=1, nsteps=1000, B=NULL, verbose=TRUE)
}
\arguments{
  \item{g}{
\code{igraph} object, the network to work on. 
}
  \item{dist.method}{
String, the method used to calculate the distance between each vertex pair.
}
  \item{vertex.attr}{
Character vector, the name of the vertex attributes under which the vertex weights to be tested are stored. The vector can contain one or more elements.
}
  \item{edge.attr}{
String, the name of the edge attribute to be used as distances along the edges.
}
\item{correct.factor}{
Numeric value, if the network contains unconnected vertices, then the distance between these vertices is set as the maximum distance between the connected vertices multiplied by \code{correct.factor}.
}
  \item{nsteps}{
Integer value, the number of bins into which vertex pairs are placed.
}
  \item{B}{
Symmetric numerical matrix. A precomputed distance bin matrix for \code{g} output by the \code{BinGraph} function. If \code{NULL}, then B is computed within the \code{Knode} function. 
}
  \item{verbose}{
Logical, if \code{TRUE} messages about the progress of the function are displayed. 
}
}
\details{
Using the inner sum of the \code{Knet} equation, it becomes possible to prioritise vertices by how well they are connected, or associated, with high-weight vertices. The inner sum of the Knet equation is

\deqn{K^{node}_i[s]=\frac{2}{p} \sum_j (p_j - \bar{p}) I(d^g(i,j)<=s)}{K^node_i[s] = 2/p * sum_j(p_j - bar{p}) (dg(i,j)<=s)}

where \eqn{p_j} is the weight of vertex \eqn{j}, \eqn{\bar{p}} is the mean vertex weight across all vertices, and \eqn{I(dg[i,j]<=s)} is an identity function, equaling 1 if vertex \eqn{i} and vertex \eqn{j} are within distance \eqn{s} and \eqn{0} otherwise. 

If the name of each vertex is stored within a vertex attribute called \code{name}, then the returned scores are labelled with these names.

Vertex weights should be greater than or equal to 0, or equal to NA if the weight is missing. The Knode statistic is still computed for vertices with missing weights. 

If an edge attribute with this name is not found, then each edge is assumed to have a distance of 1. Smaller edge distances denote stronger interactions between vertex pairs
}
\value{
A sorted named numerical vector of \code{Knode} AUKs for each vertex. 

If one vertex attribute is input, then the \code{Knode} AUKs are calculated and a single numerical vector is returned. If more than one vertex attribute is input, then a list of vectors, one for each set of vertex weights, is returned.
}
\references{
Cornish, A.J. and Markowetz, F. (2014) \emph{SANTA: Quantifying the Functional Content of Molecular Networks.}. PLOS Computational Biology. 10:9, e1003808.
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk} and Florian Markowetz
}
\seealso{
\code{\link[SANTA:Knet]{Knet}}
}
\examples{
# create a network with a single cluster of high-weight vertices
# rank all vertices by their strength of association with the high-weight vertices
g1 <- erdos.renyi.game(15, p.or.m=0.3, directed=FALSE)
g1 <- SpreadHits(g1, h=3, lambda=10)
Knode(g1, vertex.attr="hits")
plot(g1)
}
