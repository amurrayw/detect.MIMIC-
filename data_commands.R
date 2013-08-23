


# converts the Tetrad .r.txt representation to a graphNEL object.
read.dag <- function(file){
	orig.mat <- read.table(file=file)
	
	orig.mat[orig.mat==1] <-0
	orig.mat[orig.mat==-1] <-1
	
	final.graph <- igraph.to.graphNEL(graph.adjacency(orig.mat))

	return(final.graph)
}


# Takes in a graphNEL object, and generates normally distributed data from it.
generate.data.from.dag <- function(graph, n=100, errDist="normal"){	

	top.sort <- topological.sort(igraph.from.graphNEL(graph))
	var.names <- nodes(graph)[top.sort]
	
	graph<-igraph.to.graphNEL(graph.adjacency(as(graph, "matrix")[var.names,
	 var.names]))
	
	generated.data<-data.frame(rmvDAG(dag=graph, n=n, errDist=errDist))
	names(generated.data) <- var.names
	
	return(generated.data)
}


