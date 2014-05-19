\name{Compactness}
\alias{Compactness}
\title{
Measure the strength of association between a phenotype and a network using compactness scores.
}
\description{
The compactness score of set of hits on a network is the mean distance between each pair of hits. By comparing this score to permuted set of hits, a measure of the strength of association between a phenotype and a network can be computed. This method is not as effective as the \code{Knet} function and is included only to allow for comparisons to be made.
}
\usage{
Compactness(g, nperm = 100, dist.method = "shortest.paths", vertex.attr = "pheno", 
edge.attr = "distance", correct.factor = 1, only.pval = F, D = NULL, verbose = TRUE)
}
\arguments{
  \item{g}{
An igraph object, the graph to work on. 
}
  \item{nperm}{
An integer value, the number of permutations to be completed. In each permutation, the vertex weights are randomly redistributed across the graph and the compactness scores recalculated. 
}
  \item{dist.method}{
A string, the method used to calculate the distance between each hit and every other vertex in the graph.
}
  \item{vertex.attr}{
A character vector, containing the name of the vertex attributes under which the hits to be tested are stored. The vector can contain one or more elements. If more than one set of vertex weights are tested, then results for each set of weights are returned as a list. The weight of a vertex should be 1 if it is a hit, 0 if it is not a hit and NA if the information is missing. Vertices with missing weights are still included within the graph. However, their scores do not contribute to the final compactness scores. 
}
  \item{edge.attr}{
A string, containing the name of the edge attribute to be used as distances along the edges. If an edge attribute with this name is not found, then each edge is assumed to have a distance of 1. 
}
\item{correct.factor}{
A numeric value, the value by which the maximum measured vertex pair distance is multiplied when used to replace infinite distances. Infinite vertex pair distances can arrise when either the \code{shortest.paths} or \code{mfpt} distance measures are used and the graph is not fully connected.
}
  \item{only.pval}{
A logical value. If \code{TRUE}, then only the p-values for each vertex attribute are returned. 
}
  \item{B}{
 A numerical matrix. A precomputed distance matrix for \code{g} output by the \code{DistGraph} function. If \code{NULL}, then D is computed by the \code{Compactness} function. 
}
  \item{verbose}{
A logical value. If \code{TRUE}, messages about the progress of the function are displayed. 
}
}
\details{
The compactness score is used by the PathExpand tool by Glaab et al. (2010). It is a measure of the mean distance between a set of genes on a network. By comparing the compactness score of an observed set of hits, to sets of permuted hits, it is possible to produce a p-value describing the strength of association between the gene set and the network, much like the \code{Knet} function. This is not some done within the original paper.

The compactness score is only included within the \code{SANTA} package to allow for comparisons to be made. Unlike the \code{Knet} function, it cannot be applied to continuous distributions of vertex weights. It can also result in biases if there is large variability in density across the network it is being applied to. 
}
\value{
If one vertex attribute is input, \code{Compactness} is run on the single set of vertex weights and a list containing the statistics below is returned. If more than one vertex attribute is input, then \code{Compactness} is run on each set of vertex weights and a list containing an element for each vertex attribute is returned. Each element contains a sub-list containing the statistics below for the relavent vertex attribute.
\item{score.obs}{The observed compactness score}
\item{score.perm}{The permuted compactness scores. \code{NA} if no permutations are completed.}
\item{pval}{The p-value, calculated from a z-score derived from the observed and permuted compactness scores. \code{NA} if no permutations are completed.}
}
\references{
Paper in preparation.

Glaad, E., Baudot A., Krasnogor N. and Valencia A. (2010). \emph{Extending pathways and processes using molecular interaction networks to analyse cancer genome data} BMC Bioinformatics. 11(1): 597:607.
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk} 
}
