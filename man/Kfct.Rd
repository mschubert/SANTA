\name{Kfct}
\alias{Kfct}
\title{
Compute statistics central to \code{Knet} and \code{Knode} functions. 
}
\description{
Deprecated function. Will be removed in the next version of SANTA.
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
