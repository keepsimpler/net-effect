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

heterogeneity.and.tolerances = ldply(1:5, function(i) {
  ldply(1:length(graphs[[i]]), function(j) {
    graph = graphs[[i]][[j]]$B
    heterogeneity = get.degree.heterogeneity(graph)
    ldply(1:10, function(k) {
      ret = NULL
      if (length(heterogeneity.and.robust[[i]][[j]][[k]]) > 0) {
        nstar.init = sum(heterogeneity.and.robust[[i]][[j]][[k]][[2]][[1]]$nstar)
        tolerances.and.fragility = get.tolerance(heterogeneity.and.robust[[i]][[j]][[k]][[2]])
        tolerances = tolerances.and.fragility$tolerance.species
        fragility = tolerances.and.fragility$fragility
        tolerance.total = sum(tolerances)
        ret = c(heterogeneity, tolerances, tolerance.total, nstar.init, fragility)
      }
      ret
    })
  })
})
