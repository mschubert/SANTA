\name{CreateGrid}
\alias{CreateGrid}
\title{
Generate a grid-like graph.
}
\description{
Generate a graph with a grid-like arrangement of edges.
}
\usage{
CreateGrid(n = 100)
}
\arguments{
  \item{n}{
An integer value, the number of vertices to be included. 
}
}
\details{
This is a simple algorithm that creates a grid-like graph. Vertices are arranged in the largest square lattice possible. Vertices not included within this square are added as an additional row. Vertices are connected by edges to their closest neighbours.
}
\value{
An igraph object.
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk}
}
\seealso{
\code{\link[SANTA:CreateGraph]{CreateGraph}}
}
\examples{
# Generate and plot a grid-like graph containing 100 vertices. 
g <- CreateGrid(n = 100)
plot(g, layout=layout.fruchterman.reingold)
}
