library(igraph)

graph = as.matrix(M_SD_001)  # get the empirical network
graph[graph > 0] = 1  # remove the possible weights of links
row.num = nrow(graph)
col.num = ncol(graph)
degree.mean = sum(graph) / (row.num + col.num)

# generate list of list of graphs which have different degree heterogeneity 
# but same node number and mean degree with the empirical network
graphs = llply(1:2, .parallel = TRUE, function(i) graphs.rewiring(row.num, col.num, degree.mean))
graphs.assort = llply(1:5, .parallel = TRUE, function(i) graphs.rewiring.swapping(row.num, col.num, degree.mean))

gamma0.max = get.gamma0.max(graph = graph.maxnest, beta0 = 1, beta1 = 0., delta = 0., tol = 0)
#gamma0.max.bydelta = ldply(seq(0, 1.5, 0.1), function(i) get.gamma0.max(graph = graph.maxnest, beta0 = 1, beta1 = 0., delta = i, tol = 0))
parms = parms.lv2(graph = graph, alpha.row.mu = 0.2, alpha.row.sd = 0.15, alpha.col.mu = 0.2, alpha.col.sd = 0.15, beta0.mu = 1, beta0.sd = 0,
                  beta1.mu = 0., beta1.sd = 0, gamma.mu = gamma0.max - 0.001, gamma.sd = 0, h.mu = 0.0, h.sd = 0, delta = 0.)
parms$r = rowSums(parms$C - parms$M)
parms$r =  svd(parms$C - parms$M)$u[,28]

init = init.lv2(parms)
A = sim.ode.one(model = model.lv2, parms = parms, init = init)
A = sim.ode(model = model.lv2, parms = parms, init = init, isout = TRUE, iter.steps = 100, steps = 200,
            perturb = perturb, perturb.type = 'lv2.growth.rate.dec')

ode.nstars = laply(A, function(one) {
  one$nstar
})
matplot(ode.nstars, type = 'l', lwd = 1.)
text(0, ode.nstars[1,],1:length(ode.nstars[1,]))
#legend("right", c("F", "S"), lty = 1:2, bty = "n")
