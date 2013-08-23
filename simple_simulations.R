

source("construct_graph.R")
source("data_commands.R")


graph1.10000 <- replicate(test.both.methods("sim.graph.1.r.txt",
 sample.size=10000), n=500)
graph2.10000 <- replicate(test.both.methods("sim.graph.2.r.txt",
 sample.size=10000), n=500)
graph3.10000 <- replicate(test.both.methods("sim.graph.3.r.txt",
 sample.size=10000), n=500)
graph4.10000 <- replicate(test.both.methods("sim.graph.4.r.txt",
 sample.size=10000), n=500)
graph5.10000 <- replicate(test.both.methods("sim.graph.5.r.txt",
 sample.size=10000), n=500)

graph1.1000 <- replicate(test.both.methods("sim.graph.1.r.txt",
 sample.size=1000), n=500)
graph2.1000 <- replicate(test.both.methods("sim.graph.2.r.txt",
 sample.size=1000), n=500)
graph3.1000 <- replicate(test.both.methods("sim.graph.3.r.txt",
 sample.size=1000), n=500)
graph4.1000 <- replicate(test.both.methods("sim.graph.4.r.txt",
 sample.size=1000), n=500)
graph5.1000 <- replicate(test.both.methods("sim.graph.5.r.txt",
 sample.size=1000), n=500)

graph1.500 <- replicate(test.both.methods("sim.graph.1.r.txt",
 sample.size=500), n=500)
graph2.500 <- replicate(test.both.methods("sim.graph.2.r.txt",
 sample.size=500), n=500)
graph3.500 <- replicate(test.both.methods("sim.graph.3.r.txt",
 sample.size=500), n=500)
graph4.500 <- replicate(test.both.methods("sim.graph.4.r.txt",
 sample.size=500), n=500)
graph5.500 <- replicate(test.both.methods("sim.graph.5.r.txt",
 sample.size=500), n=500)


graph1.250 <- replicate(test.both.methods("sim.graph.1.r.txt",
 sample.size=250), n=500)
graph2.250 <- replicate(test.both.methods("sim.graph.2.r.txt",
 sample.size=250), n=500)
graph3.250 <- replicate(test.both.methods("sim.graph.3.r.txt",
 sample.size=250), n=500)
graph4.250 <- replicate(test.both.methods("sim.graph.4.r.txt",
 sample.size=250), n=500)
graph5.250 <- replicate(test.both.methods("sim.graph.5.r.txt",
 sample.size=250), n=500)


# So unreliable, not worth testing...

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


save.image("inital.tests")


# Runs both FA test, as well as detect.MIMIC
test.both.methods <- function(file="sim.graph.1.r.txt", sample.size=1000,
 alpha=.01, pval=.05, cut.off=.3){
	
	dataset <- generate.data.set(file=file, sample.size=sample.size)[[1]]
	
	clean.dataset <- dataset[,-grep(pattern="[L:digit:]",
	 names(dataset))]
	
	# print(head(clean.dataset))
	
	fa.results <- test.fa(clean.dataset,
		 n.latents=ncol(dataset)-ncol(clean.dataset), cut.off=cut.off)
		
	detect.mimic.results <- find.mimic(clean.dataset, alpha=alpha,
			 pval=pval)
		
			
	return(list(fa.results=fa.results,
		 detect.mimic.results=detect.mimic.results))
		
}

# Runs runs factor analysis on datset.
test.fa <- function(dataset, n.latents, cut.off=.3){
	
	# Ensures latents removed from dataset.
	# dataset <- dataset[,-grep(pattern="[L:digit:]",
	#  names(dataset))]
	
	var.names <- names(dataset)
	
	fa.model <- factanal(dataset, factors=n.latents)
	
	adj.mat <- matrix(FALSE, nrow=(ncol(dataset)+n.latents),
	 ncol=(ncol(dataset)+n.latents), dimnames=list("row"=c(var.names,
		 1:n.latents), "col"=c(var.names, 1:n.latents)))
	
	for(i in 1:n.latents){
		adj.mat[fa.model$loadings[,1]>=cut.off, i+length(var.names)] <- TRUE
	}
	
	adj.mat[(length(var.names)+1):(length(var.names)+n.latents),
	 (length(var.names)+1):(length(var.names)+n.latents)] <- FALSE
	
	return(adj.mat)
}



generate.data.set <- function(file="sim.graph.1.r.txt", sample.size=10000){
	graph.data <- list(generate.data.from.dag(read.dag(file),
	 n=sample.size))
	
	return(graph.data)
}


test.graphs <- function(n.graphs=1:5, seed.set=NULL, sample.size=1000, plot.graphs=FALSE){
	
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



