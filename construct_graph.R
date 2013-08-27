# TODO: Clean up code.
# TODO: Break up long functions into smaller helper functions
# TODO: Comment this spaghetti code

# require(pcalg)
# require(igraph)
# require(gRim)
# require(plotrix)

find.mimic <- function(data, alpha=.01, indepTest=gaussCItest, pval=.05, print.intermediate=FALSE){
	pc.model <- find.pc.model(data=data, alpha=alpha, indepTest=indepTest)

	# If no edges are directed, then return the undirected pc graph.
	if(length(degree(pc.model)) < 2 ||
	length(degree(pc.model)) >2){return(pc.model)}

	input.outputs <- find.in.out(pc.model)
	
	# If no outputs have been found
	if(is.null(input.outputs$inputs) ||
		is.null(input.outputs$outputs)){return(pc.model)}
	
	latent.structure <- finding.latent.structure(input.outputs, pc.model)
	
	sobers.step <- sobers.criterion(latent.structure, data,
		 input.outputs, pval)

	last.pc <- final.pc.run(data=data, alpha=alpha, indepTest=indepTest)

	mimic.model.list<-convert.list.to.adj.mat(list.obj=sobers.step,
		 inputs.and.outputs=input.outputs, var.names=names(data))

		# n.latents <- ncol(mimic.model.list)-ncol(data)
		# names(mimic.model.list) <- c(names(data),
		#  paste("L", 1:(n.latents), sep=""))

		
	mimic.model.graph <- igraph.to.graphNEL(graph.adjacency(mimic.model.list))
		if(print.intermediate){
			print("inputs and outputs")
			print(input.outputs)
			print("latent.structure - in find.MIMIC")
			print(latent.structure)
			print("sobers.step")
			print(sobers.step)
			print("mimic.model.list")	
			print(mimic.model.list)
		}
		
	final.model <- last.step(inputs.outputs=input.outputs, pc.graph=last.pc,
		 mimic.graph=mimic.model.graph)
		
		
	pre.sober.model <- convert.list.to.adj.mat(list.obj=latent.structure,
			 inputs.and.outputs=input.outputs, var.names=names(data))
	
	# names(pre.sober.model) <- nodes(final.model)
		
	pre.sober.model <- igraph.to.graphNEL(graph.adjacency(pre.sober.model))
		
	return(list("pc.depth.0"=pc.model, inputs.outputs=input.outputs,
	 latent.structure=latent.structure, pre.sober.model=pre.sober.model, 
	sobers.step=sobers.step, last.pc=last.pc,
	 mimic.model.graph=mimic.model.graph, final.model=final.model))
}

# Finds PC model. Determines optimal depth via recursion.
find.pc.model<-function(data, depth=0, prev.graph=0, indepTest=gaussCItest,
	 alpha=0.01, suffStat=0, n=0, p=0){

	n <- nrow(data)
	p <- ncol(data)

	## define sufficient statistics
	suffStat <- list(C = cor(data), n = n)
	pc.model<-pc(suffStat, indepTest, p, alpha,
			 m.max=depth)@graph
	return(pc.model)
}


# Determines inputs/outputs. Takes a GraphNEL object as input. returns names of inputs and outputs
find.in.out <- function(graph){
	indegree.0<- degree(graph)$inDegree==0
	inputs <- c(which(indegree.0))

	candidate.outputs <- which(!indegree.0)

	adj.matrix <- as(graph, "matrix")
	
	adj.matrix <- adj.matrix[,candidate.outputs]
	
	outputs <- apply(adj.matrix, 2,
		 function(pos.output){
			if(sum(which(pos.output==1)%in%inputs)>0){
				return(pos.output)
			}})
	if(class(outputs)=="list"){		
		outputs<-remove.null.from.list(outputs)
		# Extracts the list vector, containg the names of each output.		
		outputs<-names(outputs)
	}
	else{
		# Extracts the second list vector, containg the col names of matrix
		outputs <- dimnames(outputs)[[2]]
	}
	inputs <- names(inputs)
	return(list(inputs=inputs, outputs=outputs))
}

finding.latent.structure <- function(inputs.and.outputs, graph){
	latent.list <- finding.latents(inputs.and.outputs, graph)
	latent.list <- remove.null.from.list(latent.list)

		if(length(latent.list)==0){
			adj.matrix <- as(graph, "matrix")
			var.names <- names(data.frame(adj.matrix))
			
			latent.list[['1']]$inputs <-
			 var.names[as.numeric(inputs.and.outputs$inputs)]
			
			latent.list[['1']]$outputs <-
			 var.names[as.numeric(inputs.and.outputs$outputs)]
			
			
			latent.list[['1']]$latent <- c(1,1)
			return(latent.list)
		}

	for(i in 1:(length(latent.list))){
		
		latent.pair <- smallest.two.subsets(latent.list,
			 n.inputs=length(inputs.and.outputs$inputs))

		smallest <- latent.pair$smallest
		nextSmallest <- latent.pair$nextSmallest

		if(smallest==nextSmallest){next()}
			
			for(j in 1:length(latent.list))	{
				if(j!=smallest){
					latent.list[[j]]$inputs <-
					 remove.subset(latent.list[[smallest]]$inputs,
							latent.list[[j]]$inputs)
					
				}
				
			}

			latent.list[[nextSmallest]]$latent <-
			 c(latent.list[[nextSmallest]]$latent,
				 as.character(smallest), as.character(nextSmallest))
		}
	return(latent.list)
}

finding.latents <- function(inputs.and.outputs, graph){
	adj.matrix <- as(graph, "matrix")
	var.names <- names(data.frame(adj.matrix))
	
	inputs  <- inputs.and.outputs$inputs
	outputs <- inputs.and.outputs$outputs
	
	input.parents <- adj.matrix[inputs,]
	
	if(class(input.parents)=="matrix"){
		input.parents <- data.frame(input.parents[,
			unique(which(as(input.parents, "matrix")==1, arr.ind=T)[,2])])
	}
	else{		
		input.parents <- data.frame(matrix(which(input.parents==1), nrow=1),
		 row.names=inputs)
		names(input.parents) <- paste("X", input.parents, sep="")
		input.parents[(input.parents>0)]<-1
	}
	latent.list <- construct.latent.list(input.parents,	
			var.names=var.names)

	return(latent.list)	
}

construct.latent.list <- function(input.parents, var.names){
	latent.list <- c()
	
	for(i in 1:ncol(input.parents)){

		if(list.exactly.contains(latent.list,
			 names(col.same(input.parents, column=i)))){next()}
		else{
			
			outputs <- names(input.parents)[col.same(input.parents, column=i)]
			latent.list[[i]] <- list(outputs=outputs,
				 inputs=var.names[as.numeric(get.row.names(input.parents,
					 outputs))])
		}
	}
	return(latent.list)		
}

sobers.criterion <- function(latent.structure, data, inputs.and.outputs,
	 pval=.05){
	inputs<-c()
	outputs<-c()
	
		latents <- get.latents(latent.structure)
		
		for(i in 1:length(latents)){
			if(is.null(latents)){break()}
			if(is.null(latent.structure[[i]]$latent)){next()}
				
			inputs <- c(get.inputs.via.latents(latent.structure,
				 latents[[i]]))
			outputs <- c(get.outputs.via.latents(latent.structure,
					 latents[[i]]))
					
			inputs <- unique(inputs)
			outputs <- unique(outputs)
			
			dsep.inputs <- find.dsep(inputs, outputs, data, pval)
			if(is.null(dsep.inputs)){return(latent.structure)}
			min.set <- 1
			
			for(j in 1:length(dsep.inputs)){
				if(length(dsep.inputs[[j]]) < length(dsep.inputs[[min.set]])){
					min.set<-j
				}
			}
				latent.structure[[as.numeric(latent.structure[[i]]$latent[2])]]$inputs <-
			 unique(c(dsep.inputs[min.set],	latent.structure[[as.numeric(
				 latent.structure[[i]]$latent[2])]]$inputs))
								latent.structure[[as.numeric(latent.structure[[i]]$latent[1])]]$inputs <-
			 unique(c(dsep.inputs[min.set],	latent.structure[[as.numeric(
				 latent.structure[[i]]$latent[1])]]$inputs))


				if(length(latent.structure[[i]]$latent)==2){
					latent.structure[[i]]$latent<-NULL
				}
				else{
					latent.structure[[i]]$latent <- c(latent.structure[[i]]$
						latent[3:length(latent.structure[[i]]$latent)])
					
				}
			
		}
		return(latent.structure)
}

# Last step of the PC run
final.pc.run <- function(data, depth=0, prev.graph=NULL,
	 indepTest=gaussCItest, alpha=0.01, suffStat=0, n=0, p=0){
	
	if(depth==0){
		n <- nrow(data)
		p <- ncol(data)

		## define sufficient statistics
		suffStat <- list(C = cor(data), n = n)
	}
	
	new.graph <- pc(suffStat, indepTest, p, alpha,
		 m.max=depth+1)@graph	
		
	if(is.null(prev.graph)){prev.graph<-pc(suffStat,
		 indepTest, p, alpha, m.max=depth)@graph}
	
	prev.graph <- igraph.from.graphNEL(prev.graph)
	new.graph <- igraph.from.graphNEL(new.graph)
	
	if(isTRUE(all.equal(as(as.undirected(prev.graph), "matrix"),
	 as(as.undirected(new.graph), "matrix")))){
		return(igraph.to.graphNEL(prev.graph))
	}
	else{
		new.graph <- igraph.to.graphNEL(new.graph)
		return(final.pc.run(data=data, depth=depth+1, prev.graph=new.graph,
			 indepTest=indepTest, alpha=alpha, suffStat=suffStat, n=n, p=p))
	}
	
}

# Converts Sobers step to adj mat.
# TODO: Break into helper functions
convert.list.to.adj.mat <- function(list.obj, inputs.and.outputs, var.names){
	inputs <- inputs.and.outputs$inputs
	outputs <- inputs.and.outputs$outputs
	n.variables <- length(var.names)
	n.latent <- c()
	n.unconnected.latents <- 0
	
	
	for(i in 1:length(list.obj)){
		if(!is.null(list.obj[[i]]$latent)){
			n.latent<-(c(n.latent, list.obj[[i]]$latent))
		}
	}
	
	for(i in 1:length(list.obj)){
		if(is.null(list.obj[[i]]$latent) &&
		isFALSE(as.character(i)%in%n.latent)){
			n.unconnected.latents <- n.unconnected.latents+1
		}
	}
	
	n.latent<-length(unique(n.latent))+n.unconnected.latents

	adj.mat <- matrix(nrow=length(var.names)+n.latent,
	 ncol=length(var.names)+n.latent, data=rep(FALSE,
		 (length(var.names)+n.latent)*(length(var.names)+n.latent)))
	
	adj.mat <- data.frame(adj.mat)

	names(adj.mat) <- c(var.names, 1:n.latent)
	row.names(adj.mat) <- c(var.names, 1:n.latent)
	
	# assign latents to their positions
	for(i in 1:n.latent){
		if(n.latent>=1){
			adj.mat[n.variables+i,
			 ] <- c(var.names%in%unlist(list.obj[[i]]$outputs), rep(FALSE,
				 n.latent))

			adj.mat[, n.variables+i
			 ] <- c(var.names%in%unlist(list.obj[[i]]$inputs), rep(FALSE,
				 n.latent))


			if(!is.null(list.obj[[i]]$latent)){
				total.latents <- length(list.obj[[i]]$latent)
				
				
				for(j in 1:(ceiling(total.latents/2))){
					# As latents are stored in ordered, pairs, ensures that
					# the odd position=left, even position=right
					left.lat <-list.obj[[i]]$latent[2*(j)-1] 
					right.lat <- list.obj[[i]]$latent[2*(j)]
					
					adj.mat[length(var.names)+as.numeric(left.lat),
					 length(var.names)+as.numeric(right.lat)]<-TRUE
				}
			}
		}
	}
	# Removes self-causing latents
	diag(adj.mat)<-FALSE
	return(adj.mat)
}

last.step <- function(inputs.outputs, pc.graph, mimic.graph){

	inputs <- as.numeric(inputs.outputs$inputs)
	outputs <- as.numeric(inputs.outputs$outputs)
	
	pc.adj.matrix <- data.frame(as(pc.graph, "matrix"))
	mimic.adj.matrix <- data.frame(as(mimic.graph, "matrix"))

	n.vars <- ncol(pc.adj.matrix)
	
	n.latents <- ncol(mimic.adj.matrix)-n.vars
		
	false.outputs <- names(which(apply(pc.adj.matrix[inputs,
	 outputs], 2, sum)==0))
	
	names(mimic.adj.matrix)	<- c(names(pc.adj.matrix), paste("L",
	 1:n.latents, sep=""))
	

	if(!is.null(false.outputs) && length(false.outputs)>0){
			for(j in false.outputs){
				
				false.connected.latents<-(which(mimic.adj.matrix[
				 j]==1))
				
				mimic.adj.matrix[false.connected.latents,
				 which(names(mimic.adj.matrix)%in%j)] <- 0
				
			}

		for(i in false.outputs){
			false.connected.latents<-(which(mimic.adj.matrix[
			 i]==1))
			mimic.adj.matrix[(1:n.vars),
			 which(names(mimic.adj.matrix)%in%i)] <-
			 pc.adj.matrix[,which(names(mimic.adj.matrix)%in%i)]

	mimic.adj.matrix[which(names(mimic.adj.matrix)%in%i),
	-(false.connected.latents)] <-
	 pc.adj.matrix[which(names(mimic.adj.matrix)%in%i),]
			
			}
			
	}
	return(igraph.to.graphNEL(graph.adjacency(mimic.adj.matrix)))
}



# Helper Methods

# Returns the index of columns that are all identical.
col.same <- function(mat, column){	
	return(which(colSums(mat[,column]==mat)==nrow(mat)))
}

list.exactly.contains <- function(list.object, search.term){
	return(isTRUE(sum(unlist(lapply(list.object, 
		function(item){(search.term %in% item$outputs)})))>0))
}

get.row.names <- function(mat, col.names){
	return(row.names(mat)[unique(which(mat[c(col.names)]==1, arr.ind=T)[,1])])
}

remove.null.from.list <- function(list.object){
	purged.list <- list.object[-which(sapply(list.object,	
			is.null),arr.ind=TRUE)]
	
	if(length(purged.list)==0){return(list.object)}
	else{return(purged.list)}
}

smallest.two.subsets <- function(latent.list, n.inputs){
	smallest <- NULL
	nextSmallest <- NULL

	for(i in 1:(length(latent.list))){		
		for(j in 1:length(latent.list)){

			if(i==j){next()}
			
			if(is.subset(latent.list[[i]]$input, latent.list[[j]]$input)){

				# Ensures that subsets are being compared for smallest/n.small				
				if(is.null(smallest)){smallest<-i; nextSmallest<-j}

				if(isTRUE(length(latent.list[[i]]$inputs) <=
				 length(latent.list[[smallest]]$inputs))){smallest<-i}
				else if(isTRUE(length(latent.list[[i]]$inputs) <=
				 length(latent.list[[nextSmallest]]$inputs))){nextSmallest<-i}

			}
			else if (is.subset(latent.list[[j]]$input,
				 latent.list[[i]]$input)){
					
				# Ensures that subsets are being compared for smallest/n.small
				if(is.null(smallest)){smallest<-j; nextSmallest<-i}
				
				if(isTRUE(length(latent.list[[j]]$inputs) <=
				 length(latent.list[[smallest]]$inputs))){smallest<-j}
				else if(isTRUE(length(latent.list[[j]]$inputs) <=
				 length(latent.list[[nextSmallest]]$inputs))){nextSmallest<-j}
			}
		}
	}
	
	if(is.null(smallest)){return(list(smallest=1, nextSmallest=1))}

	return(list(smallest=smallest, nextSmallest=nextSmallest))
}

there.is.a.subset <- function(latent.list, n.inputs){
	size.list <- smallest.two.subsets(latent.list,
		 n.inputs=length(inputs.and.outputs$inputs))
		return(size.list[[1]]==size.list[[2]])
	
}


remove.subset <- function(small.set, larger.set){
	return(larger.set[!larger.set%in%small.set])
}

is.subset <- function(set.1, set.2){
	
	if(length(set.1)>length(set.2)){return(FALSE)}
	
	joint.membership <- c(0)
	for(element in set.1){
		joint.membership <- c(joint.membership + sum(set.2==element))
	}
	
	if(joint.membership==length(set.1)){
		return(TRUE)
	}
	else{
		return(FALSE)
	}
	
}

get.latents <- function(latent.structure){
	latents <- c()
	for(i in 1:length(latent.structure)){
		latents[[i]] <- latent.structure[[i]]$latent
	}
	return(latents)
}


get.inputs.via.latents <- function(latent.structure, latents){
	latents <- as.numeric(latents)
	inputs <- c()
	
	left.latent <- latents[1]
	right.latent <- latents[2]
	
	inputs <- c(latent.structure[[left.latent]]$inputs,
		latent.structure[[right.latent]]$inputs)
	
	return(inputs)
}

get.outputs.via.latents <- function(latent.structure, latents){
	latents <- as.numeric(latents)
	outputs <- c()
	
	left.latent <- latents[1]
	right.latent <- latents[2]
	
	outputs <- c(latent.structure[[left.latent]]$outputs,
		latent.structure[[right.latent]]$outputs)
		
	return(outputs)
}

find.dsep <- function(inputs, outputs, data, pval=.05){
		variations.list <- c()
		n.inputs <- length(inputs)
		
		for(i in 1:n.inputs){
			variations.list[[i]] <- list(combn(inputs, m=i))
		}

		condi.sets<-variations.list

		for(i in 1:length(condi.sets)){
			dsep.set<-get.cond.combo(data, outputs, condi.sets[i], pval)
			if(!is.null(dsep.set)){
				# print(dsep.set)
				return(dsep.set)
			}
		}		
	return()
}

get.cond.combo <- function(data, outputs, input, pval){

	input <- destroy.list(input)

	# Handles single input cases.
	if(is.null(dim(input))){

		for(i in 1:length(input)){
			if(!ciTest(data, set=c(outputs[c(1, length(outputs))],
			 (input[i])))$p.val<=pval){return(input[i])}
		}
	}
	else{	
		for(i in 1:ncol(input)){
			if(!ciTest(data, set=c(outputs[c(1, length(outputs))],
			 (input[,i])))$p.val<=pval){return(input[,i])}
		}
	return(c())	
	}	
}

destroy.list <- function(list.obj){
	if(class(list.obj)=="list"&&
		listDepth(list.obj)>1){
		
		destroy.list(unlist(list.obj, recursive=F))
	}
	else{return(list.obj[[1]])}
}

isFALSE <- function(truth.vector){return(!isTRUE(truth.vector))}












