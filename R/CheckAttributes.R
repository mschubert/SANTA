CheckAttributes <- function(
    g, 
    vertex.attr="pheno", 
    edge.attr="distance"
) {
    # check that vertex and edge weights are present in a graph and convert them if neccessary
    
    # check that all supplied vertex attributes are present on the graph
    if (!all(vertex.attr %in% list.vertex.attributes(g))) stop("not all vertex attributes found on graph")
    
    for (attr in vertex.attr) {
        vertex.weights <- get.vertex.attribute(g, attr)
        
        # identify the vertex weights that are missing
        missing <- is.na(vertex.weights)
        
        # check that none of the vertex weights are infinite. If there are any, return error
        if (sum(vertex.weights[!missing] == Inf) > 0) stop(attr, " vertex weights contain infinite values")
        
        # check that vertex weights are numeric. If not try to convert. If conversion is not possible, return error
        if (!is.numeric(vertex.weights)) {
            if (is.logical(vertex.weights)) {
                missing.char <- is.na(vertex.weights)
                vertex.weights[missing.char] <- F
            } else {
                missing.char <- vertex.weights == "NA"
                vertex.weights[missing.char] <- "0"     	
            }
            
            suppressWarnings(vertex.weights <- as.numeric(vertex.weights))
            if (sum(is.na(vertex.weights)) > 0) stop("unable to convert non-numeric vertex weights to numerals")
            vertex.weights[missing.char] <- NA
            g <- remove.vertex.attribute(g, attr)
            g <- set.vertex.attribute(g, attr, value=vertex.weights)
            warning("non-numeric vertex weights found and converted to numerals")
        }
        
        # check that all vertex weights are greater or equal to 0. If not, return error
        if (!all(vertex.weights[!is.na(vertex.weights)] >= 0)) stop(attr, " vertex weights contain negative values")
    }
    
    # check that edge distances are present, if not then they are all set equal to 1. This is not an error, as some graphs do not have edge distances
    if (is.null(get.edge.attribute(g, edge.attr))) {
        message(paste("no edge attribute with name '", edge.attr, "' found so all edge distances set to 1", sep=""))
        g <- set.edge.attribute(g, edge.attr, value=rep(1, ecount(g)))
    }
    
    # check that edge distances are numeric. If not try to convert. If conversion is not possible, return error. 
    if (!is.numeric(get.edge.attribute(g, edge.attr))) {
        edge.distances <- as.numeric(get.edge.attribute(g, edge.attr))
        if (sum(is.na(edge.distances)) > 0) stop("unable to convert non-numeric edge distances to numerals")
        g <- remove.edge.attribute(g, edge.attr)
        g <- set.edge.attribute(g, edge.attr, value=edge.distances)
        warning("non-numeric edge distances found and converted to numerals")
    }
    
    # if negative edge distances present, return error
    if (!all(get.edge.attribute(g, edge.attr) >= 0)) stop("edge distances contain negative values")
    
    # ensure edge distances range between 0 and 1
    edge.distances <- get.edge.attribute(g, edge.attr)
    g <- set.edge.attribute(g, edge.attr, value=edge.distances / max(edge.distances))
    g
}
