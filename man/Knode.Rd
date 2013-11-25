\name{Knode}
\alias{Knode}
\title{
Rank vertices by their strength of association with high-weight vertices. 
}
\description{
Rank vertices by their strength of association with high-weight vertices using a modified version of Ripley's K-statistic. Vertex weights can either be binary or positive and continuous.
}
\usage{
Knode(g, dist.method = "shortest.paths", vertex.attr = "pheno", edge.attr = "distance", 
correct.factor = 1, nsteps = 1000, only.Knode = TRUE, vertex.weight = TRUE, 
cluster.id = FALSE, vertex.degree = TRUE, boncich.power = FALSE, markov.centr = FALSE, 
B = NULL, verbose = TRUE)
}
\arguments{
  \item{g}{
An igraph object, the graph to work on. 
}
  \item{dist.method}{
A string, the method used to calculate the distance between each vertex. Can either be \code{shortest.paths}, \code{diffusion} or \code{mfpt}.
}
  \item{vertex.attr}{
A character vector, containing the name of the vertex attributes under which the vertex weights to be tested are stored. The vector can contain one or more elements. If more than one set of vertex weights are tested, then results for each set of weights are returned in a list. Vertex weights should be greater than or equal to 0, or equal to NA if the weight is missing. The Knode statistic is still counted for nodes with missing weights. However, the weight of these nodes is excluded in the calculation of other node's Knode statistics. 
}
  \item{edge.attr}{
A string, containing the name of the edge attribute to be used as distances along the edges. If an edge attribute with this name is not found, then each edge is assumed to have a distance of 1. 
}
\item{correct.factor}{
A numeric value, the value by which the maximum measured vertex pair is multiplied when used to replace infinite distances. Infinite vertex pair distances can arrise when either the \code{shortest.paths} or \code{mfpt} distance measures are used and not all vertices within the graph are connected.
}
  \item{nsteps}{
An integer value, the number of bins into which the vertex pairs are put according to distance before the Knode function is calculated. Greater values of \code{nsteps} result in greater accuracy and greater run times. 
}
  \item{only.Knode}{
A logical value. If \code{TRUE}, only the Knode AUK for each vertex is calculated. If \code{FALSE}, then other centrality scores are also calculated. These scores can include cluster ID, vertex degree, vertex betweenness score, Boncich power centrality, Burt's constraint, eigenvector centrality, Google PageRank, Kleinberg's hub and authority scores and the Markov centrality score. The calculation of some of these scores can be switched on or off using the parameters below. 
}
  \item{vertex.weight}{
A logical value. If \code{vertex.weight==TRUE} and \code{only.Knode==FALSE}, then the weight of each vertex is returned.
}
  \item{cluster.id}{
A logical value. If \code{cluster.id==TRUE} and \code{only.Knode==FALSE} and the graph contains a \code{hit.cluster} vertex attribute, then the cluster ID for each hit (each vertex with a non-zero weight) is returned. The \code{hit.cluster} attribute is added by the \code{CreateGraph} function. 
}
  \item{vertex.degree}{
A logical value. If \code{vertex.degree==TRUE} and \code{only.Knode==FALSE}, then the degree of each vertex is returned. The degree of a vertex is the number of adjacent edges. 
}
  \item{boncich.power}{
A logical value. If \code{boncich.power==TRUE} and \code{only.Knode==FALSE}, then the Boncick power centrality of each vertex is returned.
}
  \item{markov.centr}{
A logical value. If \code{markov.centr==TRUE} and \code{only.Knode==FALSE}, then the Markov centrality score of each vertex is returned.
}
  \item{B}{
 A numerical matrix. A precomputed distance bin matrix for \code{g} output by the \code{BinGraph} function. If \code{NULL}, then B is computed by the \code{Knode} function. 
}
  \item{verbose}{
A logical value. If \code{verbose==TRUE}, messages about the progress of the function are displayed. 
}
}
\details{
Using the inner sum of the Knet equation, it becomes possible to prioritise vertices by how well they are connected, or associated, with high-weight vertices. The inner sum of the Knet equation is

\deqn{K^{node}_i[s]=\frac{2}{p} \sum_j (p_j - \bar{p}) I(d^g(i,j)<=s)}{K^node_i[s] = 2/p * sum_j(p_j - bar{p}) (dg(i,j)<=s)}

where \eqn{p_j} is the weight of vertex \eqn{j}, \eqn{\bar{p}} is the mean vertex weight across all vertices, and \eqn{I(dg[i,j]<=s)} is an identity function, equaling 1 if vertex \eqn{i} and vertex \eqn{j} are within distance \eqn{s} and \eqn{0} otherwise. 

The \code{Knode} function can also be used to return a number of other centrality measures, including cluster ID, vertex degree, vertex betweenness score, Boncich power centrality, Burt's constraint, eigenvector centrality, Google PageRank, Kleinberg's hub and authority scores and the Markov centrality score. If the name of each vertex is stored within a vertex attribute called \code{name}, then these names are applied to the rows of the returned data frame. Otherwise, the rows are named with the vertex number. 
}
\value{
A sorted data frame containing the \code{Knode} AUK for each vertex and any other centrality scores. 

If one vertex attribute is input, then the \code{Knode} AUK and other centrality scores are calculated and a single sorted data frame containing these scores is returned. If more than one vertex attribute is input, then a list of data frames, one for each set of vertex weights, is returned.
}
\references{
Bonacich, P. (1987). \emph{Power and Centrality: A Family of Measures}. American Journal of Sociology. 92: 1170-1182. 

Brin, S. and Page, L. (1998). \emph{The Anatomy of a Large-Scale Hypertextual Web Search Engine}. Proceedings of the 7th World-Wide Web Conference, Bisbane, Australia. 

Burt, R.S. (2004). \emph{Structural holes and good ideas}. American Journal of Sociology. 110: 349-399.

Kleinberg, J. (1997). \emph{Authoritative sources in hyperlinked environment}. Proc. 9th. ACM-SIAM Symposium on Discrete Algorithms. 

Ulrik Brandes. (2001). \emph{A Faster Algorithm for Betweenness Centrality}. Journal of Mathematical Sociology. 25 (2): 163-177. 

White, S. and Smyth, P. (2003). \emph{Algorithms for Estimating Relative Importance in Networks}. Technical Report: UCI-ICS 04-25. 
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk} and Florian Markowetz
}
\seealso{
\code{\link[SANTA:Knet]{Knet}}
}
\examples{
# Create a graph with a single cluster of high-weight vertices. Rank all vertices 
# by their stength of association with the high-weight vertices. 
g1 <- CreateGraph(n=15, gen.vertex.weights=TRUE, lambda=10, nhits=3, binary.pheno=TRUE)
Knode(g1, only.Knode=FALSE)
plot(g1)

# Create a graph with two clusters of high-weight vertices. Rank all vertices by 
# their strength of association with the high-weight vertices. 
g2 <- CreateGraph(n=15, gen.vertex.weights=TRUE, lambda=10, nhits=6, nclusters=2, 
binary.pheno=FALSE)
Knode(g2, only.Knode=FALSE, cluster.id=TRUE)
plot(g2)
}
