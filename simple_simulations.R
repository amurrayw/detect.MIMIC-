





# Runs both FA test, as well as detect.MIMIC
test.both.methods <- function(file="sim.graph.1.r.txt", sample.size=1000,
 alpha=.01, pval=.05, cut.off=.3, scree=FALSE){
	
	dataset <- generate.data.set(file=file, sample.size=sample.size)[[1]]
	
	clean.dataset <- dataset[,-grep(pattern="[L:digit:]",
	 names(dataset))]
	
	# print(head(clean.dataset))
	
	fa.results <- test.fa(clean.dataset,
		 n.latents=ncol(dataset)-ncol(clean.dataset), cut.off=cut.off,
		 scree==scree)
		
	detect.mimic.results <- find.mimic(clean.dataset, alpha=alpha,
			 pval=pval)
		
			
	return(list(fa.results=fa.results,
		 detect.mimic.results=detect.mimic.results))
		
}

# Runs runs factor analysis on datset.
test.fa <- function(dataset, n.latents, cut.off=.3, scree=FALSE){
	
	# Ensures latents removed from dataset.
	# dataset <- dataset[,-grep(pattern="[L:digit:]",
	#  names(dataset))]
	
	var.names <- names(dataset)
	
	if(scree){
		
		
	}
	else{
		fa.model <- factanal(dataset, factors=n.latents)
	}
	
	
	return(fa.model)
}



generate.data.set <- function(file="sim.graph.1.r.txt", sample.size=10000){
	graph.data <- list(generate.data.from.dag(read.dag(file),
	 n=sample.size))
	
	return(graph.data)
}


test.graphs <- function(n.graphs=1:7, seed.set=NULL, sample.size=1000, plot.graphs=TRUE){
	
	if(!is.null(seed.set)){set.seed(seed.set)}

	found.graphs <- list()
	
	files.to.load<-paste("sim.graph.",n.graphs,".r.txt", sep="")
	for(i in files.to.load){
		
		current.graphs <- test.graph(i, sample.size=sample.size)
		
		found.graphs<-list(found.graphs, current.graphs)
		
		if(plot.graphs){plot.test(current.graphs)}
	}
	return(found.graphs)
}

# Reads in graph file, and generates normally distributed data from it.
test.graph <- function(file, sample.size=1000){
	true.graph<-read.dag(file=file)

	generated.data <- generate.data.from.dag(graph=true.graph, n=sample.size)

	generated.data.no.latents <- generated.data[,-grep(pattern="[L:digit:]",
	 names(generated.data))]
	
	result<-find.mimic(generated.data.no.latents)
	return(list(result=result, true.graph=true.graph))
}

# Plots the various stages of the algorithm used in finding the graph.
plot.test<-function(graph.list){
	
	results<-graph.list$result
	true.graph<-graph.list$true.graph
	
	# No results, so just plots true graph
	if(is.null(results)){plot(true.graph); return()}
	else if(class(results)!="list"){plot(results)}
	else{
		print(results)
		
		par(mfrow=c(2,3))
		if(!is.null(true.graph)){plot(true.graph, main="True Graph")}

		if(!is.null(results$pc.depth.0)){plot(results$pc.depth.0,
			 main="PC Depth=0")}
		
		if(!is.null(results$pre.sober.model)){
				plot(results$pre.sober.model, main="Pre-Sober")}

		if(!is.null(results$mimic.model.graph)){
			plot(results$mimic.model.graph, main="Post-Sober")}	
		
		if(!is.null(results$last.pc)){plot(results$last.pc, 
			main="PC Depth>0")}

		if(!is.null(results$final.model)){plot(results$final.model, 
			main="Final Graph")}
	}
	
}



# Too unreliable with so little data (e.g., can end up with no outputs,
# leading to errors in the adj.matrix construction).

# graph1.100 <- replicate(test.both.methods("sim.graph.1.r.txt",
#  sample.size=100), n=500)
# graph2.100 <- replicate(test.both.methods("sim.graph.2.r.txt",
#  sample.size=100), n=500)
# graph3.100 <- replicate(test.both.methods("sim.graph.3.r.txt",
#  sample.size=100), n=500)
# graph4.100 <- replicate(test.both.methods("sim.graph.4.r.txt",
#  sample.size=100), n=500)
# graph5.100 <- replicate(test.both.methods("sim.graph.5.r.txt",
#  sample.size=100), n=500)



