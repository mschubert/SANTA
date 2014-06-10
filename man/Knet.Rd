\name{Knet}
\alias{Knet}
\title{
Measure the strength of association between a phenotype and a network by computing the strength of hit clustering on the network
}
\description{
Compute the strength of clustering of high-weight vertices (hits) on a network using a modified version of Ripley's K-statistic. This method can be used to measure the strength of association between a phenotype or function and a network. 
}
\usage{
Knet(g, nperm=100, dist.method=c("shortest.paths", "diffusion", "mfpt"),
vertex.attr="pheno", edge.attr=NULL, correct.factor=1, nsteps=1000, 
prob=c(0, 0.05, 0.5, 0.95, 1), parallel=NULL, B=NULL, verbose=TRUE)
}
\arguments{
  \item{g}{
An igraph object. The network to work on. 
}
  \item{nperm}{
An integer value. The number of permutations to be completed.
}
  \item{dist.method}{
A string. The method used to calculate the distance between each hit and every other vertex in the network.
}
  \item{vertex.attr}{
A character vector. The name of the vertex attributes under which the vertex weights to be tested are stored. The vector can contain one or more elements. Vertex weights should be greater or equal that zero or equal to NA if the weight is missing.
}
  \item{edge.attr}{
A string. The name of the edge attribute to be used as distances along the edges. If an edge attribute with this name is not found, then each edge is assumed to have a distance of 1. 
}
\item{correct.factor}{
A numeric value. If the network contains unconnected vertices, then the distance between these vertices is set as the maximum distance between the connected vertices multiplied by \code{correct.factor}.
}
  \item{nsteps}{
An integer value. The number of bins into which vertex pairs are placed. Greater values of \code{nsteps} result in greater accuracy and greater run times. 
}
  \item{prob}{
A numeric vector. The quantiles to be calculated for the \code{Knet} permutations.
}
  \item{parallel}{  
A numeric value or \code{NULL}. If parallel computing is available, this value denotes the number of cores that the computation will be split over. 
}
  \item{B}{
 A symmetrical numeric matrix. A precomputed distance bin matrix for \code{g} output by the \code{BinGraph} function. If \code{NULL}, then B is computed within the \code{Knet} function. 
}
  \item{verbose}{
A logical value. If \code{verbose==TRUE}, messages about the progress of the function are displayed. 
}
}
\details{
The SANTA method uses the 'guilt-by-association' principle to measure the strength of association between a network and a phenotype. It does this by measuring the strength of clustering of the phenotype scores across the network. The stronger the clustering, the greater the association between the network and the phenotype. 

The SANTA method applies Ripley's K-function, a well-established approach to spatial statistics that measures the strength of clustering of points on a plane, and extends it in a number of ways. First, a Knet function is defined by adapting the approach for networks using vertex pair distance measures. Second, vertex weights are incorporated into Knet and the importance of vertices made relative to their own associated weight. Third, the mean vertex weight is subtracted from each individual vertex weight when calculating the Knet function. This means that the Knet function measures the degree of vertex weight clustering relative to a random distribution of vertex weights. The Knet function is defined as 

\deqn{K^{net}[s]=\frac{2}{p^2} \sum_i p_i \sum_j(p_j - \bar{p}) I(d^g[i,j]<=s)}{K^{net}[s] = 2/p^2 * sum_i{(p_i) sum_j{(pj - bar{p}) I(dg(i,j)<=s)}}}

where \eqn{p_i} is the weight of vertex \eqn{i}, \eqn{\bar{p}} is the mean vertex weight across all vertices, and \eqn{I(dg[i,j]<=s)} is an identity function, equaling 1 if vertex \eqn{i} and vertex \eqn{j} are within distance \eqn{s} and \eqn{0} otherwise. 

In order to derive a p-value and quantify the significance of the observed distribution of weights, the observed Knet-curve is compared to Knet-curves obtained using the same network but randomly permuted vertex weights. Vertices with missing weights (\code{NA}) are not included within these permutations. The area under the Knet-curve (AUK) is calculated for the observed network and each of the permuted networks and a z-score used to produce a p-value. This p-value indicates the probability an observed AUK at least this high is seen given the null hypothesis that the vertex weights are randomly distributed.

If parallel computing is possible, \code{parallel} can be used to split permutations over multiple cores. The \code{snow} package is used to manage the parallel computing. If \code{parallel=NULL} or parallel computing is not possible, then only one core is used. If a positive integer is input and parallel computing is possible, then the permutations are split over up to this many cores.
}
\value{
If one vertex attribute is input, \code{Knet} is run on the single set of vertex weights and a list containing the statistics below is returned. If more than one vertex attribute is input, then \code{Knet} is run on each set of vertex weights and a list containing an element for each vertex attribute is returned. Each element contains a sub-list containing the statistics below for the relavent vertex attribute.
\item{K.obs}{The Knet-function curve for the observed vertex weights.}
\item{AUK.obs}{The area under the Knet-function curve (AUK) for the observed vertex weights.}
\item{K.perm}{The Knet-function curve for each permutation of vertex weights. Equals \code{NA} if no permutations are completed.}
\item{AUK.perm}{The area under the Knet-function curve (AUK) for each permutation of vertex weights. Equals \code{NA} if no permutations are completed.}
\item{K.quan}{The quantiles for the permuted Knet-function curves. Equals \code{NA} if no permutations are completed.}
\item{pval}{The p-value, calculated from a z-score derived from the observed and permuted AUKs. Equals \code{NA} if no permutations are completed.}
}
\references{
Paper in preparation.

Okabe, A. and Yamada, I. (2001). \emph{The K-function method on a network and its computational implementation} Geographical Analysis. 33(3): 271-290.
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk} and Florian Markowetz
}
\seealso{
\code{\link[SANTA:Knode]{Knode}}
}
\examples{
# apply Knet to a network with hit clustering
g.clustered <- barabasi.game(50, directed=FALSE)
g.clustered <- SpreadHits(g.clustered, h=10, lambda=10)
res.clustered <- Knet(g.clustered, nperm=100, vertex.attr="hits")
res.clustered$pval
plot(res.clustered)

# apply Knet to a network without hit clustering
g.unclustered <- barabasi.game(50, directed=FALSE)
g.unclustered <- SpreadHits(g.unclustered, h=10, lambda=0)
res.unclustered <- Knet(g.unclustered, nperm=100, vertex.attr="hits")
res.unclustered$pval
plot(res.unclustered)
}
