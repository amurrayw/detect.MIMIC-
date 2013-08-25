

score.mimic <- function(mimic.model, true.graph){
	
	adj.matrix.true <- data.frame(true.graph)
	
	adj.matrix.mimic <- data.frame(as(mimic.model, "matrix"))
	
	n.latents.truth <- grep(pattern="[L:digit:]",
	 names(adj.matrix.true))

	n.latents.mimic <- grep(pattern="[L:digit:]",
	 names(adj.matrix.mimic))
	
	
	
}

score.fa <- function(fa.model, cut.off=.3, true.graph){
	
	fa.model <- prune.fa.paths(fa.model, cut.off=cut.off)
	
	
	
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
		adj.mat[fa.model$loadings[,1]>cut.off, i+length(var.names)] <- TRUE
	}
	
	adj.mat[(length(var.names)+1):(length(var.names)+n.latents),
	 (length(var.names)+1):(length(var.names)+n.latents)] <- FALSE
	
	return(adj.mat)
	
}