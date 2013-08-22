

source("construct_graph.R")
source("data_commands.R")


test.graphs <- function(n.graphs=1:5, seed.set=NULL, sample.size=1000){
	if(!is.null(seed.set)){set.seed(seed.set)}
	
	
	
	files.to.load<-paste("sim.graph.",n.graphs,".r.txt", sep="")
	for(i in files.to.load){
		
		print("RUNNING ON GRAPH")
		print(i)
		found.graphs<-test.graph(i, sample.size=sample.size)
		# print(found.graphs$pc.depth.0)
		
		plot.test(found.graphs)
	}
	
}





test.graph <- function(file, sample.size=1000){
	true.graph<-read.dag(file=file)

	generated.data <- generate.data.from.dag(graph=true.graph, n=sample.size)

	generated.data.no.latents <- generated.data[,-grep(pattern="[L:digit:]",
	 names(generated.data))]
	
	result<-find.mimic(generated.data.no.latents)
	return(list(result=result, true.graph=true.graph))
}

# TODO: Add in plotting the pre-last step graph.
plot.test<-function(graph.list){
	
	results<-graph.list$result
	true.graph<-graph.list$true.graph
	# No results, so just plots true graph
	if(is.null(results)){plot(true.graph); return()}
	else if(class(results)!="list"){plot(results)}
	else{
		print(results)
		
		par(mfrow=c(2,2))
		if(!is.null(true.graph)){plot(true.graph, main="True Graph")}

		if(!is.null(results$pc.depth.0)){plot(results$pc.depth.0,
			 main="PC Depth=0")}

		if(!is.null(results$last.pc)){plot(results$last.pc, 
			main="PC Depth>0")}

		if(!is.null(results$final.model)){plot(results$final.model, 
			main="Final Graph")}
	}
	
}



