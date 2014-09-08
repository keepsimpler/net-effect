
library(igraph)

graph = M_SD_001
graph[graph > 0] = 1  # remove the possible weights of links
row.num = nrow(graph)
col.num = ncol(graph)
degree.mean = sum(graph) / (row.num + col.num)
graphs = llply(1:5, .parallel = TRUE, function(i) graphs.rewiring(row.num, col.num, degree.mean))

heterogeneity.and.robust = llply(1:5, function(i) {
  llply(graphs[[i]], .parallel = TRUE, function(graph){
    graph = graph$B
    #heterogeneity = get.degree.heterogeneity(graph)
    llply(1:10, function(i) {
      ret = NULL
      parms = parms.lv2(graph)
      init = init.lv2(parms)      
      A = sim.ode.one(model = model.lv2, parms, init)
      if (A$extinct == 0) {
        A = sim.ode(model = model.lv2, parms = parms, init = init, isout = FALSE, iter.steps = 100,
                perturb = perturb, perturb.type = 'lv2.growth.rate.dec')
        ret = list(graph = graph, A = A)
      }
      ret
      #c(heterogeneity = heterogeneity, feasible = feasible)
    })
  })
})

