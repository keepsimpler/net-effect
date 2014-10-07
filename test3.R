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

gamma0.max = get.gamma0.max(graph = graph.rand, beta0 = 1, beta1 = 0., delta = 0.5, tol = 0)
#gamma0.max.bydelta = ldply(seq(0, 1., 0.1), function(i) get.gamma0.max(graph = graph, beta0 = 1, beta1 = 0., delta = i, tol = 0))
parms = parms.lv2(graph = graph.rand, alpha.row.mu = 1., alpha.row.sd = 0.2, alpha.col.mu = 1., alpha.col.sd = 0.2, beta0.mu = 1, beta0.sd = 0,
                  beta1.mu = 0., beta1.sd = 0, gamma.mu = gamma0.max - gamma0.max + 0.01, gamma.sd = 0, h.mu = 0., h.sd = 0, delta = 0.)
parms$r = rowSums(parms$C - parms$M)
parms$r =  - svd(parms$C - parms$M)$u[, nrow(parms$C)]

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


deltas = seq(0, 1, 0.1)
beta1s = seq(0, 0.2, 0.02)
hs = seq(0, 0.2, 0.02)
res = ldply(graphs, function(graphs.2) {
  ldply(graphs.2, .parallel = TRUE, function(graph) {
    print(graph$count)
    graph = graph$B
    ldply(deltas, function(delta) {
      ldply(beta1s, function(beta1) {
        ldply(hs, function(h) {
          gamma0.max = get.gamma0.max(graph = graph, beta0 = 1, beta1 = beta1, delta = delta, tol = 0)
          gamma0s = seq(gamma0.max, 0.01, length.out = 20)
          ldply(gamma0s,.parallel = FALSE, function(gamma0) {
            parms = parms.lv2(graph = graph, beta0.mu = 1, beta0.sd = 0, beta1.mu = beta1, beta1.sd = 0,
                              gamma.mu = gamma0, gamma.sd = 0, h.mu = h, h.sd = 0, delta = delta)
            parms$r = abs(svd(parms$C - parms$M)$u[, nrow(parms$C)])
            solve(parms$C - parms$M) %*% parms$r
            #if (any(solve(parms$C - parms$M) %*% parms$r < 0)) print(parms$r)
            #else print(paste(delta, beta1, h, gamma0))
            #init = init.lv2(parms)
            #A = sim.ode(model = model.lv2, parms = parms, init = init, isout = FALSE, iter.steps = 100, steps = 200,
            #            perturb = perturb, perturb.type = 'lv2.growth.rate.dec')
            
          })
        })
      })
    })
  })
})
