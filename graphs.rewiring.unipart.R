#!/usr/bin/Rscript
#########################################################################################
# Generate uni-partite undirect graphs with different degree heterogeneity by rewiring links.
# 1), generate a ER random graph.
# 2), repeat rewiring a link to a new node who has more neighbors than the old node 
#     until failed for enough times.
# 3), return to step 2)

source('nestedness.r')

#' @param n, number of nodes
#' @param k, average node degree
graphs.rewiring.unipart <- function(n, k) {
  graphs.rewiring = list()  # initialize the generated graphs by rewiring links
  ## generate a connected random graph as the start point
  G = graph.connected(s = n, k = k, gtype = 'er')  
  A = get.adjacency(G)  # get the adjacency matrix of the random graph
  A = as.matrix(A)
  
  #B = A
  count = 0
  repeat {
    count = count + 1
    shouldcontinue = FALSE
    ## rewiring one link to a random node which has more neighbors
    ## if tring [ntry]*5 times, and still fail to rewire, then [shouldcontinue] is false.
    for (i in 1:5) {
      A = rewirelinks.richer.onestep.unipart(A, ntry = 5000)
      if (A$flag == TRUE) {  # the rewiring is success
        shouldcontinue = TRUE
        break
      }
      else {
        A = A$A
      }
    }
    if (!shouldcontinue) break
    A = A$A  # the new graph which has more degree heterogeneity
    graphs.rewiring[[length(graphs.rewiring) + 1]] =  list(count = count, A = A)
    print(count)
  }
  graphs.rewiring
}

###############################################################################
#' @title Rewiring algrithm by rewiring links to nodes with more links. (richer get richer)
#'
#' @param A, adjacency matrix of the random graph,
#' @param connected, if the new graph should be connected?
#' @param ntry, how many to try?
#' @return the adjacency matrix whose degree heterogeneity has been increased by rewiring links.
rewirelinks.richer.onestep.unipart <- function(A, connected = TRUE, ntry = 100) {
  # sort node degree descending, ensure the chosen species later has more interactions
  A = A[order(rowSums(A), decreasing = TRUE), order(rowSums(A), decreasing = TRUE)]
  n <- dim(A)[1]
  flag = FALSE # is the rewiring succeed, or the max tried times approach but the rewiring still not succeed 
  ## random choose another node (plant or animal with equal probability), 
  ## and rewire the link to the new selectd species
  for (i in 1:ntry) {
    flag1 = FALSE  #  if this rewiring has succeed?
    flag2 = FALSE  #  if the new graph is connected?
    A2 = A  # copy the original graph, in order to try more than one times
    ## pick one interaction between two random species
    repeat {
      row1 <- sample(1:n, 1)
      col1 <- sample(1:n, 1)
      if (A2[row1, col1] != 0 && row1 != col1) break
    }
    if (runif(1) < 0.5) {  # choose another plant  #NumP/(NumP + NumA)
      row2 =  sample(1:row1, 1)  # choose random plant with more interactions which is ensured by [sortweb]
      # Three exceptions: 1. the new chosen species [row2] is same with the old species [row1]
      # 2. the new chosen species [row2] already has interaction with [col1]
      # 3. the old species [row1] has only one interaction.
      # 4. the new node [row2] is same with the old node [col1]
      if (row2 < row1 && A2[row2, col1] == 0 && sum(A2[row1,]) > 1 && row2 != col1) {
        A2[row2, col1] = A2[row1, col1]
        A2[col1, row2] = A2[col1, row1]
        A2[row1, col1] = 0
        A2[col1, row1] = 0
        flag1 = TRUE  # the link has been rewired to a new plant
      }
    }
    else {  # choose another animal
      col2 =  sample(1:col1, 1)
      if (col2 < col1 && A2[row1, col2] == 0 && sum(A2[,col1]) > 1 && col2 != row1) {
        A2[row1, col2] = A2[row1, col1]
        A2[col2, row1] = A2[col1, row1]
        A2[row1, col1] = 0
        A2[col1, row1] = 0
        flag1 = TRUE  # the link has been rewired to a new animal
      }
    }
    ## if the new graph is connected, [flag2] is TRUE
    G = graph.adjacency(A2)
    if (igraph::is.connected(G)) flag2 = TRUE
    
    ## if the rewiring is success, and (the new graph is connected or that is not required)
    if (flag1 & (flag2 | !connected)) {
      flag = TRUE
      break; 
    }
  }
  if (flag == FALSE) {  # if failed, return the original matrix
    res = list(A = A, flag = flag, tried = i)
    warning(paste(ntry, 'times has been tried, but the rewiring still donot succeed.')) 
  }
  else {
    #A2 = A2[order(rowSums(A2), decreasing = TRUE), order(rowSums(A2), decreasing = TRUE)]
    res = list(A = A2, flag = flag, tried = i)
  }
  res
}

get.nstars <- function(res) {
  ldply(res, function(res.1) {
    ldply(res.1, function(one) {
      if (length(one) > 0) one$nstar
    })
  })
}

get.selfvars <- function(res) {
  ldply(res, function(res.1) {
    ldply(res.1, function(one) {
      if (length(one) > 0) diag(one$vars)
    })
  })
}


#save(graphs.rewiring, file = paste('graphs.rewiring', s1, s2, k, sep = '-'))
