\name{plot.Knet}
\alias{plot.Knet}
\title{
Plot the results of the \code{Knet} function
}
\description{
Plot the observed Knet curve against the quantiles of the permuted Knet curves and the observed AUK against the permuted AUKs. 
}
\usage{
\S3method{plot}{Knet}(x, sequential=FALSE, \dots)
}
\arguments{
  \item{x}{
Results from the Knet function.
}
  \item{sequential}{
Logical, if \code{TRUE} then the plots are sequential. Otherwise, the two plots are plotted alongside each other.
}
  \item{\dots}{
Additional arguments to be passed to \code{plot}.  	
}
}
\details{
If the high-weight vertices are clustered, then the observed Knet curve and AUK will be high relative to the permuted Knet curves and AUKs. The greater the strength of clustering, the greater the difference between the observed and permuted statistics. If the strength of clustering is low, then the observed and permuted curves and AUKs will likely overlap.

The first plot displays the the observed curve in red and the quantiles of the permuted curves in yellow. The quantile boundaries are displayed as grey lines. These boundaries are specified in the Knet function. The second plot displays the observed AUK as a red line and the distribution of permuted AUKs in grey. 
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
# plot results with hit clustering
g.clustered <- barabasi.game(100, directed=FALSE)
g.clustered <- SpreadHits(g.clustered, h=10, lambda=10)    
res.clustered <- Knet(g.clustered, nperm=10, vertex.attr="hits")
res.clustered$pval
plot(res.clustered)

# plot results without hit clustering
g.unclustered <- barabasi.game(100, directed=FALSE)
g.unclustered <- SpreadHits(g.unclustered, h=10, lambda=0)
res.unclustered <- Knet(g.unclustered, nperm=10, vertex.attr="hits")
res.unclustered$pval
plot(res.unclustered)
}
