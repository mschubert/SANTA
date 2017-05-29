\name{Kset}
\alias{Kset}
\title{
Measure the strength of association between multiple phenotypes on a network.
}
\description{
Compute the strength of association between phenotypes on a network using a modified version of Ripley's K-statistic. Phenotypes are represented by vertex weights.
}
\usage{
Kset(g, nperm=100, dist.method=c("shortest.paths", "diffusion", "mfpt"),
seed.vertex.attr="seed", query.vertex.attr="query", version=c("b", "a"), fix.vertex.degree=TRUE, edge.attr=NULL, correct.factor=1, nsteps=1000, 
prob=c(0, 0.05, 0.5, 0.95, 1), parallel=NULL, B=NULL, verbose=TRUE)
}
\arguments{
  \item{g}{
\code{igraph} object, the network to work on. 
}
  \item{nperm}{
Integer value, the number of permutations to be completed.
}
  \item{dist.method}{
String, the method used to calculate the distance between vertex pairs.
}
  \item{seed.vertex.attr}{
String, the name of the vertex attribute under which the seed vertex weights are stored.
}
  \item{query.vertex.attr}{
Character vector, the name of the vertex attributes under which the query vertex weights stored. The vector can contain one or multiple elements. If this vector contains multiple elements, then the seed vertex weights are compared against each distribution of query vertex weights. 
}
  \item{version}{
Character vector, the version of Kset to use. Can either be a or b. b by default. 
}
  \item{fix.vertex.degree}{
Logical value, whether to maintain the degree of vertices when the vertex values are shuffled in the permutation analysis. TRUE by default.
}
  \item{edge.attr}{
String, the name of the edge attribute to be used as distances along the edges. If an edge attribute with this name is not found, then each edge is assumed to have a distance of 1. 
}
\item{correct.factor}{
Numeric value, if the network contains unconnected vertices, then the distance between these vertices is set as the maximum distance between the connected vertices multiplied by \code{correct.factor}.
}
  \item{nsteps}{
Integer value, the number of bins into which vertex pairs are placed.
}
  \item{prob}{
Numeric vector, the quantiles to be calculated for the \code{Kset} permutations.
}
  \item{parallel}{  
Numeric value or \code{NULL}. If parallel computing is available, this value denotes the number of cores that the computation will be split over. 
}
  \item{B}{
Symmetrical numeric matrix. A precomputed distance bin matrix for \code{g} output by the \code{BinGraph} function. If \code{NULL}, then B is computed within the \code{Kset} function. 
}
  \item{verbose}{
Logical, if \code{TRUE} messages about the progress of the function are displayed. 
}
}
\details{
TO COMPLETE.
}
\value{
If one query vertex attribute is input, \code{Kset} compares the distribution of weights defined by the seed vertex attribute to the distribution of weights defined by the query vertex attribute, and a single list of statistics (as described below) is returned. If multiple query vertex attributes are input, \code{Kset} compares the distribution of weights defined by the seed vertex attribute to the distribution of weights defined by each query vertex attribute, and a list of lists of statistics (as described is below) is output.  

\item{K.obs}{Kset-function curve for the observed vertex weights.}
\item{AUK.obs}{Area under the Kset-function curve (AUK) for the observed vertex weights.}
\item{K.perm}{Kset-function curve for each permutation of vertex weights. Equals \code{NA} if no permutations are completed.}
\item{AUK.perm}{Area under the Kset-function curve (AUK) for each permutation of vertex weights. \code{NA} if no permutations are completed.}
\item{K.quan}{Quantiles for the permuted Kset-function curves. \code{NA} if no permutations are completed.}
\item{pval}{p-value, calculated from a z-score derived from the observed and permuted AUKs. \code{NA} if no permutations are completed.}
}
\references{
Cornish, A.J. and Markowetz, F. (2014) \emph{SANTA: Quantifying the Functional Content of Molecular Networks.} PLOS Computational Biology. 10:9, e1003808.

Okabe, A. and Yamada, I. (2001). \emph{The K-function method on a network and its computational implementation} Geographical Analysis. 33(3): 271-290.
}
\author{
Alex J. Cornish, Michael Schubert and Florian Markowetz
}
\seealso{
\code{\link[SANTA:Knet]{Knet}}
}
\examples{
# use Kset to determine whether two independantly generated phenotypes, represented by vertex weights, cluster together on a network
n.vertices <- 50
n.hits <- 3
seed.weights <- sample(c(rep(1, n.hits), rep(0, n.vertices - n.hits)))
query.weights <- sample(c(rep(1, n.hits), rep(0, n.vertices - n.hits)))
g <- barabasi.game(n.vertices, directed=FALSE)
g <- set.vertex.attribute(g, name="seed", value=seed.weights)
g <- set.vertex.attribute(g, name="query", value=query.weights)
res.unclustered <- Kset(g, seed.vertex.attr="seed", query.vertex.attr="query")
res.unclustered$pval
plot(res.unclustered)

# use Kset to determine whether two similar phenotype sets, represented by vertex weights, cluster together on a network
n.vertices <- 50
n.hits <- 3
seed.weights <- sample(c(rep(1, n.hits), rep(0, n.vertices - n.hits)))
query.weights <- seed.weights
query.weights[1] <- 1
g <- barabasi.game(n.vertices, directed=FALSE)
g <- set.vertex.attribute(g, name="seed", value=seed.weights)
g <- set.vertex.attribute(g, name="query", value=query.weights)
res.clustered <- Kset(g, seed.vertex.attr="seed", query.vertex.attr="query")
res.clustered$pval
plot(res.clustered)
}
