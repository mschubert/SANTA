\name{CheckAttributes}
\alias{CheckAttributes} 
\title{
Check that a graph is associated with the correct vertex and edge attributes and can be passed to the \code{Knet} and \code{Knode} function. 
}
\description{
In order for the \code{Knet} and \code{Knode} functions to operate correctly when applied to a graph, the graph must be associated with vertex and edge attributes that represent the weight of the vertices and the distances along the edges. This function ensures that an igraph object contains these attributes and converts them to the correct format if required. 
}
\usage{
CheckAttributes(g, vertex.attr = "pheno", edge.attr = "distance")
}
\arguments{
  \item{g}{
An igraph object, the graph to work on. 
}
  \item{vertex.attr}{
A character vector, containing the one or more names of the vertex attributes under which vertex weights to be tested are stored.
}
  \item{edge.attr}{
A string, containing the name of the edge attribute to be used as distances along the edges. If no edge distances are associated with the graph, then each edge is assigned a distance of 1 and these distances saved under this attribute name. 
}
}
\details{
In order for the \code{Knet} and \code{Knode} functions to work correctly, one or more vertex attributes representing vertex weights and one edge attribute representing the distances along edges need to be present within the igraph object. \code{CheckAttributes} checks that these attributes are present and in a suitable format. If no edge distances are associated with the graph, then each edge is assigned a distance of 1 by the function. 

This function is included within the \code{Knet} and \code{Knode} functions and does not need to be run separately. 
	
If vertex weights are not present under any of the attribute names in \code{vertex.attr}, an error message is returned. 

If edge distances are not present under the attribute name \code{edge.attr}, then each edge is assigned a distance of 1 under this attribute name.

If the vertex weights or edge distances are not numeric, then they are converted to numerals. If they cannot be converted, an error message is returned.

If any of the vertex weights or edge distances are negative, an error message is returned.
}
\value{
An igraph object, the input graph with any required modifications.
}
\author{
Alex Cornish \email{a.cornish12@imperial.ac.uk}
}
\seealso{
\code{\link[SANTA:Knet]{Knet}},
\code{\link[SANTA:Knode]{Knode}}
}
\examples{
# Create a graph, assign vertex weights and use the CheckAttributes function
# to check that the vertex weights are in the correct format and assign edge 
# distances
g1 <- CreateGraph(10, gen.vertex.weights=FALSE)
g1 <- set.vertex.attribute(g1, name="weights", value=runif(vcount(g1)))
g1 <- CheckAttributes(g1, vertex.attr="weights", edge.attr="distances")
get.vertex.attribute(g1, name="weights")
get.edge.attribute(g1, name="distances")
}
