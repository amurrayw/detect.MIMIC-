library(pcalg)


# Load predefined data
data(gmG)
n <- nrow(gmG$x)
p <- ncol(gmG$x)

## define independence test (partial correlations)

indepTest <- gaussCItest

## define sufficient statistics
suffStat <- list(C = cor(gmG$x), n = n)
## estimate CPDAG
alpha <- 0.01
pc.fit <- pc(suffStat, indepTest, p, alpha, verbose = TRUE)

temp<-temp@graph

# use when testing dseperation.
dsep()


plot(graphNEL(nodes=c("a", "b", "c"), edgeL=list("a"=list("edges"=c("b","c")), "b"=list("edges"=c("c", "a")), "c"=list("edges"=c("c", "a"))), edgemode="directed"))

# Get PC

# Get.in.out - use degree function

# 

#  Find latents

# Find subset

# Connect inputs with latents

# 

# Exports a graphNEL object to graviz file.
graph2graphviz()


# library(sna)
library(igraph)
# reads in .dot file, returns adjacency matrix, and converts to graphNEL
igraph.to.graphNEL(read.dot("step2.dot"))


# Read in adjacancy matrix created by tetrad, and create an igraph object, and convert to graph object
igraph.to.graphNEL(graph.adjacency(read.table(file="graph1.r.txt")))

# Finding inputs
temp<-(degree(pc.fit@graph)$inDegree==0)
which(temp)
# Assigning possible outputs
which(!temp)

# TODO: Need to write function to remove candidate outputs without connections to inputs.



isTRUE(all.equal(igraph.to.graphNEL(graph.adjacency(read.table(file="graph1.r.txt"))), igraph.to.graphNEL(graph.adjacency(read.table(file="graph1.r.txt"))), tol = .00001))



library(RBGL)
# in pcalg library
# Generates random data from a graph
rmvDAG(dag=igraph.to.graphNEL(graph.adjacency(read.table(file="graph1.r.txt"))), n=100, errDist="normal")



rmvDAG(dag=igraph.to.graphNEL(graph.adjacency(read.table(file="graph2.r.txt")))[tsort(igraph.to.graphNEL(graph.adjacency(read.table(file="graph2.r.txt"))))] )

# Use to claculate FDR, TDR, and TPR. Is in the pcalg library

compareGraphs(gl, gt)





