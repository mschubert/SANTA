test_CheckAttributes <- function() {
  # test that CheckAttributes produces an error if not all vertex attributes present on graph
  g <- erdos.renyi.game(10, 0.5)
  g <- set.vertex.attribute(g, name="attr1", value=runif(vcount(g)))
  checkException(CheckAttributes(g, vertex.attr=c("attr1", "attr2")), silent=TRUE)
  g <- remove.vertex.attribute(g, "attr1")

  # test that CheckAttributes successfully converts non-numeric vertex weights
  values <- runif(vcount(g))
  values[1] <- "1"
  g <- set.vertex.attribute(g, name="vertex.attr", value=values)
  suppressWarnings(suppressMessages((g <- CheckAttributes(g, vertex.attr="vertex.attr"))))
  checkEquals(get.vertex.attribute(g, name="vertex.attr")[1], 1)
  g <- remove.vertex.attribute(g, "vertex.attr")
  
  # test that CheckAttributes produces an error if it is unable to convert non-numeric vertex weights
  values <- runif(vcount(g))
  values[1] <- "a"
  g <- set.vertex.attribute(g, name="vertex.attr", value=values)
  checkException(CheckAttributes(g, vertex.attr="vertex.attr"), silent=TRUE)
  g <- remove.vertex.attribute(g, "vertex.attr")
  
  # test that CheckAttributes produces an error if not all vertex weights are greater or equal to 0
  values <- runif(vcount(g))
  values[1] <- -1
  g <- set.vertex.attribute(g, name="vertex.attr", value=values)
  checkException(CheckAttributes(g, vertex.attr="vertex.attr"), silent=TRUE)
  g <- remove.vertex.attribute(g, "vertex.attr")
  
  # test that CheckAttributes sets missing edge distances to 1
  values <- runif(vcount(g))
  g <- set.vertex.attribute(g, name="vertex.attr", value=values)
  suppressMessages(g <- CheckAttributes(g, vertex.attr="vertex.attr", edge.attr="edge.attr"))
  checkEquals(get.edge.attribute(g, "edge.attr"), rep(1, ecount(g)))
  g <- remove.edge.attribute(g, "edge.attr")
  
  # test that CheckAttributes successfully converts non-numeric edge weights
  values <- runif(ecount(g))
  values[1] <- "1"
  g <- set.edge.attribute(g, name="edge.attr", value=values)
  suppressWarnings(suppressMessages((g <- CheckAttributes(g, vertex.attr="vertex.attr", edge.attr="edge.attr"))))
  checkEquals(is.numeric(get.edge.attribute(g, name="edge.attr")[1]), TRUE)
  g <- remove.edge.attribute(g, "edge.attr")
  
  # test that CheckAttributes produces an error if it is unable to convert non-numeric edge weights
  values <- runif(ecount(g))
  values[1] <- "a"
  g <- set.edge.attribute(g, name="edge.attr", value=values)
  suppressWarnings(checkException(CheckAttributes(g, vertex.attr="vertex.attr", edge.attr="edge.attr"), silent=TRUE))
  g <- remove.edge.attribute(g, "edge.attr")
  
  # test that CheckAttributes produces an error if not all edge distances are greater or equal to 0
  values <- runif(ecount(g))
  values[1] <- -1
  g <- set.edge.attribute(g, name="edge.attr", value=values)
  suppressWarnings(checkException(CheckAttributes(g, vertex.attr="vertex.attr", edge.attr="edge.attr"), silent=TRUE))
  g <- remove.edge.attribute(g, "edge.attr")
  
  # test that CheckAttributes ensures all edge.attributes range between 0 and 1 (min != 0 but max == 1)
  values <- runif(ecount(g))
  g <- set.edge.attribute(g, name="edge.attr", value=values)
  g <- CheckAttributes(g, vertex.attr="vertex.attr", edge.attr="edge.attr")
  checkEquals(min(get.edge.attribute(g, "edge.attr"))>=0, TRUE)
  checkEquals(max(get.edge.attribute(g, "edge.attr"))==1, TRUE)
}
