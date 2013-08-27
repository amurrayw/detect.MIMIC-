# Returns average fa.score and mimic.score. Note that scores are
# as follows: 

# 			tpr	
# 			True Positive Rate: Number of correctly found edges (in estimated
# graph) divided by number of true edges (in true graph)
# 
# 			fpr	
# 			False Positive Rate: Number of incorrectly found edges divided by
#  number of true gaps (in true graph)
# 
# 			tdr	
# 			True Discovery Rate: Number of correctly found edges divided by
#  number of found edges (both in estimated graph)
score.both <- function(graph.sim, graph.name, cut.off=.3){
	
	true.graph <- as(read.dag(graph.name), "matrix")
	
	var.names <- names(data.frame(true.graph))

	latent.positions <- grep(pattern="[L:digit:]", var.names)
	
	non.latent.positions <- (1:length(var.names))[-latent.positions]

	sorted.truth <- (1:length(var.names))[c(non.latent.positions,
		 latent.positions)]

	true.graph <- true.graph[sorted.truth,sorted.truth]
	
	# Scores FA model.
	fa.score <- lapply(graph.sim[1,],
		 function(fa.model){score.fa(fa.model=fa.model, cut.off=cut.off,
			 true.graph=true.graph)})

	# Scores MIMC model.
	mimic.score <- lapply(graph.sim[2,],
		 function(mimic.model){
			# Handles case when only PC output was returned.
			if(length(mimic.model)>1){
				result<-score.mimic(mimic.model=mimic.model$final.model,
			 true.graph=true.graph)
			
				return(result)
				}
			return(c(0,0,0))
			})
		
	n.with.incorrect.latents <- 0
	for(i in mimic.score){
		if(is.null(i)){n.with.incorrect.latents <- n.with.incorrect.latents+1}
	}

	fa.score <- (do.call(rbind, fa.score))
	mimic.score <- (do.call(rbind, mimic.score))

	fa.score <- colMeans(fa.score)
	mimic.score <- colMeans(mimic.score)
			
	return(list(fa.score=fa.score, mimic.score=mimic.score,
		 n.with.incorrect.latents))
	
}

score.mimic <- function(mimic.model, true.graph){
	
	adj.matrix.true <- data.frame(true.graph)
	
	adj.matrix.mimic <- data.frame(as(mimic.model, "matrix"))
	
	n.latents.truth <- grep(pattern="[L:digit:]",
	 names(adj.matrix.true))

	n.latents.mimic <- grep(pattern="[L:digit:]",
	 names(adj.matrix.mimic))
	
	true.graph<-igraph.to.graphNEL(graph.adjacency(true.graph))


	# If the mimic.model found has a different numnber of latents than the
	# true graph, then cannot calculate tpr, fpr, tdr. Have therfore treated
	# those as NULL objects. (i.e., as with FA, correct n.latents is assumed)
	if(length(nodes(mimic.model)) == length(nodes(true.graph))){
		return(compareGraphs(mimic.model, true.graph))
	}
	return(NULL)
}

# TODO: Need to add case where FA got n.latents wrong.
# Scores FA model. Note that it assumes fa.model got n.latents correct.
score.fa <- function(fa.model, cut.off=.3, true.graph){
	
	fa.mat <- as(fa.model$loadings, "matrix")
	n.latents <- ncol(fa.mat)
	var.names <- row.names(fa.mat)
	n.vars <- nrow(fa.mat)

	fa.model <- prune.fa.paths(fa.model, cut.off=cut.off)

	fa.model <- igraph.to.graphNEL(graph.adjacency(fa.model))
	true.graph <- igraph.to.graphNEL(graph.adjacency(true.graph))
	
	if(length(nodes(fa.model))==length(nodes(true.graph))){
		
		graph.comparison <- (compareGraphs(fa.model, true.graph))
		return(graph.comparison)
	}
	return(NULL)
}


prune.fa.paths <- function(fa.model, cut.off=.3){
	
	fa.loadings.matrix <- as(fa.model$loadings, "matrix")
	
	n.latents <- ncol(fa.loadings.matrix)
	var.names <- row.names(fa.loadings.matrix)
	n.vars <- length(var.names)
	
	adj.mat <- matrix(FALSE, nrow=(n.vars+n.latents),
	 ncol=(n.vars+n.latents), dimnames=list("row"=c(var.names,
		 1:n.latents), "col"=c(var.names, 1:n.latents)))
	
	for(i in 1:n.latents){
		adj.mat[abs(fa.model$loadings[,1])>cut.off, 
		i+length(var.names)] <- TRUE
	}
	
	adj.mat[(length(var.names)+1):(length(var.names)+n.latents),
	 (length(var.names)+1):(length(var.names)+n.latents)] <- FALSE
	
	return(adj.mat)
	
}

get.latent.cluster <- function(adj.matrix, n.latents, n.vars){
	
	
	latent.vectors.col <- adj.matrix[,(n.vars+1):(n.vars+n.latents)]
	latent.vectors.row <- t(adj.matrix[(n.vars+1):(n.vars+n.latents),])	
	
	latent.clusters <- (latent.vectors.col+latent.vectors.row)[
	-((n.vars+1):(n.vars+n.latents)),]
	return(latent.clusters)	
}


convert.graph.groups <- function(graph.groups = graph.groups){
	graph.final<-list()
	
	
	for(i in 1:length(graph.groups)){
		graph.list <- unlist(graph.groups[[i]])[c(3,4,1,2)]
		n.null <- c()
		graph.fa.score <- c()
		graph.mimic.score <- c()
		var.names <- c()
		for(j in 1:length(graph.list)){

			graph.score<- get(graph.list[j])
			
			graph.fa.score <- rbind(graph.fa.score, unlist(graph.score[[1]]))
			graph.mimic.score <- rbind(graph.mimic.score,
				 unlist(graph.score[[2]]))
				
			n.null <- c(n.null, unlist(graph.score[[3]]))
		}
		
		row.names(graph.fa.score) <- c("250", "500", "1000", "10000")
		row.names(graph.mimic.score) <- c("250", "500", "1000", "10000")
		
		graph.final[[i]] <- list(fa.scores=graph.fa.score,
			 mimic.scores=graph.mimic.score, n.null=n.null)
		
	}
	return(graph.final)
}

plot.rates <-function(score.list){
		for(i in 1:length(score.list)){
			
			par(mfrow=c(2,2))
			barplot(score.list[[i]]$mimic.scores[,1], 
				main=paste("MIMIC - Graph ",i, " True Positive Rate", sep=""),
				 col=1:4, ylab="Rate", xlab="Number of observations")
				
			barplot(score.list[[i]]$mimic.scores[,2], 
				main=paste("MIMIC - Graph ",i, " False Positive Rate",
				 sep=""), col=1:4, ylab="Rate", xlab="Number of observations")
				
			barplot(score.list[[i]]$mimic.scores[,3], 
				main=paste("MIMIC - Graph ",i, " True Discovery Rate",
				 sep=""), col=1:4, ylab="Rate", xlab="Number of observations")
								
			barplot(score.list[[i]]$fa.scores[,1],
				 main=paste("FA - Graph ",i, " True Positive Rate", sep=""),
				 col=1:4, ylab="Rate", xlab="Number of observations")
				
			barplot(score.list[[i]]$fa.scores[,2],
				 main=paste("FA - Graph ",i, " False Positive Rate", sep=""),
				 col=1:4, ylab="Rate", xlab="Number of observations")
				
			barplot(score.list[[i]]$fa.scores[,3],
				 main=paste("FA - Graph ",i, " True Discovery Rate", sep=""),
				 col=1:4, ylab="Rate", xlab="Number of observations")
				
				barplot(score.list[[i]]$n.null, 
					main=paste("MIMIC - Graph ",i, 
					" Number of False\n Latent Cases", sep=""), col=1:4,
					names=c("250", "500", "1000", "10000"),  
					ylab="Number of incorrect cases", 
					xlab="Number of observations")
					
		}
		par(mfrow=c(2,2))
		
		for(i in 1:length(score.list)){
			plot(read.dag(paste("sim.graph.", i, ".r.txt", sep="")),
			 main=paste("True Graph ", i, sep=""))
			
		}
		
}












