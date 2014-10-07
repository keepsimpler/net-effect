#!/usr/bin/Rscript
#########################################################################################
# Generate bipartite graphs with different degree heterogeneity by rewiring links.
# 1), generate a random bipartite graph.
# 2), repeat rewiring a link to a new species who has more neighbors than the old species 
#     until failed for enough times.
# 3), return to step 2)

source('nestedness.r')

graphs.rewiring <- function(s1, s2, k) {
  graphs.rewiring = list()  # the generated graphs by rewiring links
  #s1 = 25; s2 = 25  # number of nodes in two groups
  #k = 1.5  # average node degree
  ## generate a connected random bipartite graph as the start point
  ## of generating graphs with different degree heterogeneity by rewiring links
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
      B = rewirelinks.richer.onestep(B, ntry = 5000)
      if (B$flag == TRUE) {  # the rewiring is success
        shouldcontinue = TRUE
        break
      }
      else {
        B = B$B
      }
    }
    if (!shouldcontinue) break
    B = B$B  # the new graph which has more degree heterogeneity
    graphs.rewiring[[length(graphs.rewiring) + 1]] =  list(count = count, B = B)
    print(count)
  }
  graphs.rewiring
}

#save(graphs.rewiring, file = paste('graphs.rewiring', s1, s2, k, sep = '-'))
