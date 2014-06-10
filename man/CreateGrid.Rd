\name{CreateGrid}
\alias{CreateGrid}
\title{
Generate a grid-like network
}
\description{
Generate a network with a grid-like arrangement of edges.
}
\usage{
CreateGrid(n=100)
}
\arguments{
  \item{n}{
An integer value. The number of vertices to be included. 
}
}
\details{
This is a simple algorithm that creates a grid-like network. Vertices are arranged in the largest square lattice possible. Vertices not included within this square are added as an additional row. Vertices are connected by edges to their closest neighbours.
}
\value{
An igraph object.
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk}
}
\examples{
# generate and plot a grid-like network containing 100 vertices. 
g <- CreateGrid(n = 100)
plot(g, layout=layout.fruchterman.reingold)
}
