#!/usr/bin/Rscript
#########################################################################################
# Generate bipartite graphs with different degree heterogeneity and degree assortativity.
# 1), generate a random bipartite graph.
# 2), rewiring a link to a new species who has more neighbors than the old species.
# 3), repeat swapping two links to improve assortativity until failed for enough times.
# 4), repeat swapping two links to improve disassortativity until failed for enough times.
# 5), return to step 2)
# warning: this function is very time-consuming!

source('nestedness.r')

graphs.rewiring.swapping <- function(s1, s2, k) {
  result.graphs = list()
  #s1 = 40; s2 = 20  # number of nodes in two groups
  #k = 1.5  # average node degree
  ## generate a connected random bipartite graph as the start point
  ## of generating graphs with different degree heterogeneity and degree assortativity
  G = graph.connected(c(s1, s2), k = k, gtype = 'bipartite')  
  A = get.incidence(G)  # get the incidence matrix of the random bipartite graph
  
  B = A
  count = 0
  repeat {
    count = count + 1
    shouldcontinue = FALSE
    ## rewiring one link to a random node which has more neighbors
    ## if tring [ntry]*5 times, and still fail to rewire, then [shouldcontinue] is false.
    for (i in 1:5) {
      B = rewirelinks.richer.onestep(B, ntry = 1000)
      if (B$flag == TRUE) {
        shouldcontinue = TRUE
        break
      }
      else {
        B = B$B
      }
    }
    if (!shouldcontinue) break
    B = B$B  # the new graph which has more degree heterogeneity
    result.graphs[[length(result.graphs) + 1]] =  list(count = count, count.assort = 0, B = B)
    
    ## swapping two links to increase assortativity while maintain the degree heterogeneity
    ## if tring [ntry]*5 times, and still fail to swap, then [shouldcontinue.assort] is false.
    B2 = B
    count.assort = 0
    repeat {
      count.assort = count.assort + 1
      shouldcontinue.assort = FALSE
      for (j in 1:5) {
        B2 = swaplinks.assort.onestep(B2, ntry = 1000)
        if (B2$flag == TRUE) {
          shouldcontinue.assort = TRUE
          break
        }
        else {
          B2 = B2$B
        }
      }
      if (!shouldcontinue.assort) break
      print(paste(count, count.assort))
      B2 = B2$B  # the new graph
      result.graphs[[length(result.graphs) + 1]] =  list(count = count, count.assort = count.assort, B = B2)
    }
    
    ## swapping two links to decrease assortativity while maintain the degree heterogeneity
    ## if tring [ntry]*5 times, and still fail to swap, then [shouldcontinue.assort] is false.
    B2 = B
    count.assort = 0
    repeat {
      count.assort = count.assort - 1
      shouldcontinue.assort = FALSE
      for (j in 1:5) {
        B2 = swaplinks.disassort.onestep(B2, ntry = 5000)
        if (B2$flag == TRUE) {
          shouldcontinue.assort = TRUE
          break
        }
        else {
          B2 = B2$B
        }
      }
      if (!shouldcontinue.assort) break
      print(paste(count, count.assort))
      B2 = B2$B  # the new graph
      result.graphs[[length(result.graphs) + 1]] =  list(count = count, count.assort = count.assort, B = B2)
    }
  }
  result.graphs 
}
# save(result.graphs, file = paste('result.graphs', s1, s2, k, sep = '-'))

