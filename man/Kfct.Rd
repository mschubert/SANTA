\name{Kfct}
\alias{Kfct}
\title{
Function central to the \code{Knet} and \code{Knode} functions. 
}
\description{
This function is called by the \code{Knet} and \code{Knode} functions. It should not be called seperately. 
}
\usage{
Kfct(B, vertex.weights, individual=TRUE)
}
\arguments{
  \item{B}{
A numerical matrix. A matrix of distance bins for \code{g}; output by the \code{BinGraph} function.
}
  \item{vertex.weights}{
The weight of each vertex in \code{g}.
}
  \item{individual}{
Logical value. If \code{TRUE}, then the K and AUK is returned for each individual node. 
}
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk} and Florian Markowetz
}
\seealso{
\code{\link[SANTA:Knet]{Knet}}
\code{\link[SANTA:Knode]{Knode}}
}
