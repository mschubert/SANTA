AsiGraph <- function(v, g) {
  ## convert vertex indices to an igraph object
  if (class(v) == "igraph.vs") return(v)
  class(v) <- "igraph.vs"            
  ne <- new.env()
  assign("graph", g, envir = ne)
  attr(v, "env") <- ne
  v
}
