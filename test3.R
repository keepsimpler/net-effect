library(igraph)

graph = as.matrix(M_SD_003)  # get the empirical network
graph[graph > 0] = 1  # remove the possible weights of links
row.num = nrow(graph)
col.num = ncol(graph)
degree.mean = sum(graph) / (row.num + col.num)

# generate list of list of graphs which have different degree heterogeneity 
# but same node number and mean degree with the empirical network
graphs = llply(1:5, .parallel = TRUE, function(i) graphs.rewiring(row.num, col.num, degree.mean))
graphs.assort = llply(1:5, .parallel = TRUE, function(i) graphs.rewiring.swapping(row.num, col.num, degree.mean))

gamma0.max = get.gamma0.max(graph = graph, beta0 = 1, beta1 = 0.1, delta = 0, tol = 0)
parms = parms.lv2(graph = graph, alpha.row.mu = 1, alpha.row.sd = 0, alpha.col.mu = 1, alpha.col.sd = 0, beta0.mu = 1, beta0.sd = 0,
                  beta1.mu = 0.1, beta1.sd = 0, gamma.mu = gamma0.max - 0.1, gamma.sd = 0, h.mu = 0, h.sd = 0, delta = 0)
init = init.lv2(parms)
A = sim.ode.one(model = model.lv2, parms = parms, init = init)
