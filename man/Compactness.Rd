\name{Compactness}
\alias{Compactness}
\title{
Measure the strength of association using compactness scores
}
\description{
The compactness score of set of hits on a network is the mean distance between each pair of hits. By comparing the observed compactness score to the scores of permuted hit sets, it is possible to compute the significance of the strength of association between the phenotype and the network. This method is not as effective as the \code{Knet} function and is included only for comparison.
}
\usage{
Compactness(g, nperm=100, dist.method=c("shortest.paths", "diffusion", "mfpt"), 
vertex.attr="pheno", edge.attr="distance", correct.factor=1, D=NULL, 
verbose=T)
}
\arguments{
  \item{g}{
\code{igraph} object, the network to work on. 
}
  \item{nperm}{
Integer value, the number of permutations to be completed.
}
  \item{dist.method}{
String, the method used to compute the distance between each pair of hits on the network.
}
  \item{vertex.attr}{
Character vector, the name of the vertex attributes under which the hits to be tested are stored. The vector can contain one or more vertex attributes.
}
  \item{edge.attr}{
String, the name of the edge attribute to be used as distances along the edges. If an edge attribute with this name is not found, then each edge is assumed to have a distance of \code{1}. Smaller edge distances denote stronger interactions between vertex pairs.
}
\item{correct.factor}{
Numeric value. If the network contains unconnected vertices, then the distance between these vertices is set as the maximum distance between the connected vertices multiplied by \code{correct.factor}.
}
  \item{D}{
Symmetrical numerical matrix. A precomputed distance matrix for \code{g} output by the \code{DistGraph} function. If \code{NULL}, then \code{D} is computed by the \code{Compactness} function. 
}
  \item{verbose}{
Logical, if \code{TRUE} messages about the progress of the function are displayed. 
}
}
\details{
The compactness score is used by the PathExpand tool by Glaab et al. (2010). It is a measure of the mean distance between a set of genes in a network. By comparing the compactness score of an observed set of hits to sets of permuted hits, it is possible to produce a p-value describing the strength of association between the gene set and the network. This is not some done within the original paper by Glaab et al. (2010). The function is much like the \code{Knet} function, albeit not as effective.  

The compactness score \code{C} is defined as the mean shortest path distance between pairs of vertices in a set \code{P} on network \code{g}. 

\deqn{C(P) = \frac{2  \sum_{i,j \in P; i < j} d^g(i,j)}{|P| * (|P| - 1)}}

The compactness score is only included within the \code{SANTA} package to allow for comparisons to be made. Unlike the \code{Knet} function, it cannot be applied to continuous distributions of vertex weights. It can also result in biases if there is large variability in density across the network. 

The weight of a vertex should be 1 if it is a hit, 0 if it is not a hit or \code{NA} if the information is missing. Vertices with missing weights are still included within the network but are excluded from the permuted sets.
}
\value{
If one vertex attribute is input, \code{Compactness} is run on the single set of vertex weights and a list containing the statistics below is returned. If more than one vertex attribute is input, then \code{Compactness} is run on each set of vertex weights and a list containing an element for each vertex attribute is returned. Each element contains a sub-list containing the statistics below for the relevant vertex attribute.
\item{score.obs}{Observed compactness score}
\item{score.perm}{Permuted compactness scores. \code{NA} if no permutations are completed.}
\item{pval}{p-value, computed from a z-score derived from the observed and permuted compactness scores. \code{NA} if no permutations are completed.}
}
\references{
Cornish, A.J. and Markowetz, F. (2014) \emph{SANTA: Quantifying the Functional Content of Molecular Networks.}. PLOS Computational Biology. 10:9, e1003808.

Glaab, E., Baudot A., Krasnogor N. and Valencia A. (2010). \emph{Extending pathways and processes using molecular interaction networks to analyse cancer genome data}. BMC Bioinformatics. 11(1): 597:607.
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk} 
}
