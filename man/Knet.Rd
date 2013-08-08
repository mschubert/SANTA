\name{Knet}
\alias{Knet}
\title{
Measure the strength of association between a phenotype and a network by computing the strength of hit clustering on the network.
}
\description{
Compute the strength of clustering of high-weight vertices (hits) on a graph using a modified version of Ripley's K-statistic. This method can be used to measure the strength of association between a phenotype or function and a network. 
}
\usage{
Knet(g, nperm = 100, dist.method = "shortest.paths", vertex.attr = "pheno", 
edge.attr = "distance", correct.factor = 1, nsteps = 1000, 
prob = c(0, 0.05, 0.5, 0.95, 1), only.pval = F, parallel = NULL)
}
\arguments{
  \item{g}{
An igraph object, the graph to work on. 
}
  \item{nperm}{
An integer value, the number of permutations to be completed. In each permutation, the vertex weights are randomly redistributed across the graph and the Knet function recalculated. The permuted Knet-function results are then compared to the observed Knet function result in order to derive a p-value for the significance of the hit clustering. 
}
  \item{dist.method}{
A string, the method used to calculate the distance between each hit and every other vertex in the graph. Can either be \code{shortest.paths}, \code{diffusion} or \code{mfpt}.
}
  \item{vertex.attr}{
A character vector, containing the name of the vertex attributes under which the vertex weights to be tested are stored. The vector can contain one or more elements. If more than one set of vertex weights are tested, then results for each set of weights are returned as a list. Vertex weights should be greater or equal that zero or equal to NA if the weight is missing. Vertices with missing weights are still included within the graph. However, their weights do not contribute to the final Knet statistics. 
}
  \item{edge.attr}{
A string, containing the name of the edge attribute to be used as distances along the edges. If an edge attribute with this name is not found, then each edge is assumed to have a distance of 1. 
}
\item{correct.factor}{
A numeric value, the value by which the maximum measured vertex pair distance is multiplied when used to replace infinite distances. Infinite vertex pair distances can arrise when either the \code{shortest.paths} or \code{mfpt} distance measures are used and not all vertices within the graph are connected. 
}
  \item{nsteps}{
An integer value, the number of bins into which the vertex pairs are put according to distance before the Knet function is calculated. Greater values of \code{nsteps} result in greater accuracy and greater run times. 
}
  \item{prob}{
A numeric vector, containing the quantiles to be calculated for the \code{Knet} permutations.
}
  \item{only.pval}{
A logical value. If \code{TRUE}, then only the p-values for each vertex attribute are returned. 
}
  \item{parallel}{  
A numeric value or \code{NULL}. If parallel computing is possible, \code{parallel} can be used to split permutations over multiple cores. The \code{snow} package is used to manage the parallel computing. If \code{parallel=NULL} or parallel computing is not possible, then only one core is used. If a positive integer is input and parallel computing is possible, then the permutations are split over up to this many cores.
}
}
\details{
The SANTA method uses the 'guilt-by-association' principle to measure the strength of association between a network and a phenotype. It does this by measuring the strength of clustering of the phenotype scores across the network. The stronger the clustering, the greater the association between the network and the phenotype. 

The SANTA method applies Ripley's K-function, a well-established approach to spatial statistics that measures the strength of clustering of points on a plane, and extends it in a number of ways. First, a Knet function is defined by adapting the approach for graphs using vertex pair distance measures. Second, vertex weights are incorporated into Knet and the importance of vertices made relative to their own associated weight. Third, the mean vertex weight is subtracted from each individual vertex weight when calculating the Knet function. This means that the Knet function measures the degree of vertex weight clustering relative to a random distribution of vertex weights. The Knet function is defined as 

\deqn{K^{net}[s]=\frac{2}{p^2} \sum_i p_i \sum_j(p_j - \bar{p}) I(d^g[i,j]<=s)}{K^{net}[s] = 2/p^2 * sum_i{(p_i) sum_j{(pj - bar{p}) I(dg(i,j)<=s)}}}

where \eqn{p_i} is the weight of vertex \eqn{i}, \eqn{\bar{p}} is the mean vertex weight across all vertices, and \eqn{I(dg[i,j]<=s)} is an identity function, equaling 1 if vertex \eqn{i} and vertex \eqn{j} are within distance \eqn{s} and \eqn{0} otherwise. 

In order to derive a p-value and quantify the significance of the observed distribution of weights, the observed Knet-curve is compared to Knet-curves obtained using the same graph but randomly permuted vertex weights. The area under the Knet-curve (AUK) is calculated for the observed graph and each of the permuted graphs and a z-score derived. From this, a p-value can be produced. This p-value indicates the probability an observed AUK at least this high is seen given the null hypothesis that the vertex weights are randomly distributed.

Ripley's K-function has previously been applied to geographical networks (such as road networks) in order to identify the clustering of objects along these networks (Okabe and Yamada 2001). However, key differences between the previous implementation and the implementation of the K-function used in this package allows for the function to be applied to numerous biological networks. 

}
\value{
If one vertex attribute is input, \code{Knet} is run on the single set of vertex weights and a list containing the statistics below is returned. If more than one vertex attribute is input, then \code{Knet} is run on each set of vertex weights and a list containing an element for each vertex attribute is returned. Each element contains a sub-list containing the statistics below for the relavent vertex attribute.
\item{K.obs}{The Knet-function curve for the observed vertex weights.}
\item{AUK.obs}{The area under the Knet-function curve (AUK) for the observed vertex weights.}
\item{K.perm}{The Knet-function curve for each permutation of vertex weights. Equals \code{NA} if no permutations are completed.}
\item{AUK.perm}{The area under the Knet-function curve (AUK) for each permutation of vertex weights. Equals \code{NA} if no permutations are completed.}
\item{K.quan}{The quantiles for the permuted Knet-function curves. Equals \code{NA} if no permutations are completed.}
\item{nodeK}{The Knode-function curve for each of the vertices using the observed set of vertex weights.}
\item{nodeAUK}{The area under the Knode-function curve (AUK) for each of the vertices using the observed set of vertex weights.}
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
# Apply Knet to a graph with hit clustering
g.clustered <- CreateGraph(n=50, type="barabasi", gen.vertex.weights=TRUE, 
lambda=10, nhits=10, binary.pheno=FALSE)
res.clustered <- Knet(g.clustered, nperm=100)
res.clustered$pval
plot(res.clustered)

# Apply Knet to a graph without hit clustering
g.unclustered <- CreateGraph(n=50, type="barabasi", gen.vertex.weights=TRUE, 
lambda=0, nhits=10, binary.pheno=FALSE)
res.unclustered <- Knet(g.unclustered, nperm=100)
res.unclustered$pval
plot(res.unclustered)
}
