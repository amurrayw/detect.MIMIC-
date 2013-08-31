require(pcalg)
require(igraph)
require(gRim)
require(plotrix)
require(nFactors)



source("construct_graph.R")
source("data_commands.R")
source("simple_simulations.R")
source("simulation_scoring.R")


# Currently, this portion is not in a function so as to avoid rerunning
# everything in the event of an error (i.e., an error will not interrupt the 
# function call, leading to already completed simulations being discarded).
# Eventually, a function will be written so that this is cleaner.

# Runs the simulations.
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
graph6.10000 <- replicate(test.both.methods("sim.graph.6.r.txt",
 sample.size=10000), n=500)
graph7.10000 <- replicate(test.both.methods("sim.graph.7.r.txt",
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
graph6.1000 <- replicate(test.both.methods("sim.graph.6.r.txt",
 sample.size=1000), n=500)
graph7.1000 <- replicate(test.both.methods("sim.graph.7.r.txt",
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
graph6.500 <- replicate(test.both.methods("sim.graph.6.r.txt",
 sample.size=500), n=500)
graph7.500 <- replicate(test.both.methods("sim.graph.7.r.txt",
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
graph6.250 <- replicate(test.both.methods("sim.graph.6.r.txt",
 sample.size=250), n=500)
graph7.250 <- replicate(test.both.methods("sim.graph.7.r.txt",
 sample.size=250), n=500)


save.image("simulation.results.latents.RData")




score.graph1.10000<- score.both(graph1.10000, graph.name="sim.graph.1.r.txt")
score.graph1.1000<- score.both(graph1.1000, graph.name="sim.graph.1.r.txt")
score.graph1.500<- score.both(graph1.500, graph.name="sim.graph.1.r.txt")
score.graph1.250<- score.both(graph1.250, graph.name="sim.graph.1.r.txt")

score.graph2.10000<- score.both(graph2.10000, graph.name="sim.graph.2.r.txt")
score.graph2.1000<- score.both(graph2.1000, graph.name="sim.graph.2.r.txt")
score.graph2.500<- score.both(graph2.500, graph.name="sim.graph.2.r.txt")
score.graph2.250<- score.both(graph2.250, graph.name="sim.graph.2.r.txt")

score.graph3.10000<- score.both(graph3.10000, graph.name="sim.graph.3.r.txt")
score.graph3.1000<- score.both(graph3.1000, graph.name="sim.graph.3.r.txt")
score.graph3.500<- score.both(graph3.500, graph.name="sim.graph.3.r.txt")
score.graph3.250<- score.both(graph3.250, graph.name="sim.graph.3.r.txt")

score.graph4.10000<- score.both(graph4.10000, graph.name="sim.graph.4.r.txt")
score.graph4.1000<- score.both(graph4.1000, graph.name="sim.graph.4.r.txt")
score.graph4.500<- score.both(graph4.500, graph.name="sim.graph.4.r.txt")
score.graph4.250<- score.both(graph4.250, graph.name="sim.graph.4.r.txt")

score.graph5.10000<- score.both(graph5.10000, graph.name="sim.graph.5.r.txt")
score.graph5.1000<- score.both(graph5.1000, graph.name="sim.graph.5.r.txt")
score.graph5.500<- score.both(graph5.500, graph.name="sim.graph.5.r.txt")
score.graph5.250<- score.both(graph5.250, graph.name="sim.graph.5.r.txt")

score.graph6.10000<- score.both(graph6.10000, graph.name="sim.graph.6.r.txt")
score.graph6.1000<- score.both(graph6.1000, graph.name="sim.graph.6.r.txt")
score.graph6.500<- score.both(graph6.500, graph.name="sim.graph.6.r.txt")
score.graph6.250<- score.both(graph6.250, graph.name="sim.graph.6.r.txt")

score.graph7.10000<- score.both(graph7.10000, graph.name="sim.graph.7.r.txt")
score.graph7.1000<- score.both(graph7.1000, graph.name="sim.graph.7.r.txt")
score.graph7.500<- score.both(graph7.500, graph.name="sim.graph.7.r.txt")
score.graph7.250<- score.both(graph7.250, graph.name="sim.graph.7.r.txt")




# TODO: Rewrite this portion. It is ugly, and should be in a function.
number.graphs<-7
graph.groups <- c()
for(k in 1:number.graphs){
	graph.groups[[k]] <- list(
	 ls()[grep(pattern=paste("score.graph",
	 k, ".", sep=""), ls(), fixed=T)])
}


save.image("scored.simulations.latents.RData")


