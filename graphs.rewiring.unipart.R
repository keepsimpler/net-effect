#!/usr/bin/Rscript
#########################################################################################
# Generate uni-partite undirect graphs with different degree heterogeneity by rewiring links.
# 1), generate a ER random graph.
# 2), repeat rewiring a link to a new node who has more neighbors than the old node 
#     until failed for enough times.
# 3), return to step 2)

source('nestedness.r')

library(doMC)  # 
registerDoMC()  # register Multi Cores
getDoParWorkers()  # get available Cores


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

get.feasible.stable <- function(res) {
  notfeasible.notstable = 0
  feasible.notstable = 0
  notfeasible.stable = 0
  feasible.stable = 0
  for (i in 1:length(res)) {
    one = res[[i]]
    if (! one$feasible & ! one$stable) notfeasible.notstable = notfeasible.notstable + 1
    else if (one$feasible & ! one$stable) feasible.notstable = feasible.notstable + 1
    else if (! one$feasible & one$stable) notfeasible.stable = notfeasible.stable + 1
    else feasible.stable = feasible.stable + 1
  }
  c(notfeasible.notstable = notfeasible.notstable, feasible.notstable = feasible.notstable,
    notfeasible.stable = notfeasible.stable, feasible.stable = feasible.stable)
}

get.feasible.stable.by.variance <- function(res2) {
    ldply(res2, function(res2.1) {
      ldply(res2.1, function(res2.2) {
        ldply(res2.2, function(res2.3) {
          ldply(res2.3, function(res2.4) {
            r.mean = res2.4$r.mean
            r.sd = res2.4$r.sd
            beta1.mean = res2.4$beta1.mean
            beta1.sd = res2.4$beta1.sd
            res = res2.4$res
            feasible.stable = get.feasible.stable(res)
            c(r.mean = r.mean, r.sd = r.sd, beta1.mean = beta1.mean, beta1.sd = beta1.sd, feasible.stable)
          })
        })
      })
    })  
}

get.nstars.vars.sum <- function(res) {
  ldply(res, function(one) {
    if (one$feasible && one$stable)
      c(hetero = one$hetero, vars.self = sum(diag(one$vars)), vars.sum = sum(one$vars), nstar.sum = sum(one$nstar),
        eigs = sum(1 / (2*(eigen(one$phi)$values))), vars.max = sum(sqrt(diag(one$vars)))^2)
    #         vars.max = sum(sqrt(diag(one$vars)))^2, eigs = sum(one$nstar^2 / (2*(eigen(one$phi)$values))), 
    #          sigs = sum(one$nstar^2 / (2*(svd(one$phi)$d))), eigenvector = sum(eigen(one$phi)$vectors[,1]) )      
  })
}

get.nstars.vars.sum.by.variance <- function(res2) {
  ldply(res2, function(res2.1) {
    ldply(res2.1, function(res2.2) {
      ldply(res2.2, function(res2.3) {
        ldply(res2.3, function(res2.4) {
          r.mean = res2.4$r.mean
          r.sd = res2.4$r.sd
          beta1.mean = res2.4$beta1.mean
          beta1.sd = res2.4$beta1.sd
          res = res2.4$res
          ldply(res, function(one) {
            if (one$feasible && one$stable)
              c(r.mean = r.mean, r.sd = r.sd, beta1.mean = beta1.mean, beta1.sd = beta1.sd, hetero = one$hetero, beta1 = one$beta1, vars.self = sum(diag(one$vars)), vars.sum = sum(one$vars), nstar.sum = sum(one$nstar),
                eigs = sum(1 / (2*(eigen(one$phi)$values))), vars.max = sum(sqrt(diag(one$vars)))^2)
          })
        })
      })
    })
  })  
}

get.nstars <- function(res) {
  ldply(res, function(one) {
    if (one$feasible && one$stable) c(hetero = one$hetero, nstar = one$nstar)    
  })
}

get.selfvars <- function(res) {
  ldply(res, function(one) {
    if (one$feasible && one$stable) c(hetero = one$hetero, selfvars = diag(one$vars))    
  })
}

#' @title get degree heterogeneity of a simple undirected graph
#' @param graph, can be a igraph object or an adjacency matrix
get.degree.hetero <- function(graph) {
  if (class(graph) == 'igraph') {  # if is a igraph object of package <igraph> 
    degrees = degree(graph)
  }
  else {
    degrees = rowSums(graph)    
  }
  degrees = sort(degrees, decreasing = TRUE)
  edges = sum(degrees)
  degrees = degrees / edges
  hetero = - sum( degrees * log(degrees) )
  hetero
}

#' @title random rewiring algorithm that change degree distribution
#'        let node with richer neighbors get richer and increase degree heterogeneity
#' @param degrees, the old degree distribution
#' @param ntried, maximum tried times
#' @return the new degree distribution which is more heterogeneous,
#'         or, NULL if maximum tried times approach which new distribution still not be found
degree.richer.richer <- function(degrees, ntried = 1000) {
  n = length(degrees)  # node number
  flag = FALSE  # is link rewiring success?
  # randomly sample two different nodes whose neighbors are larger than 1 and less than n-1
  for (i in 1:ntried) {
    two.rand = sample(1:n, size = 2, replace = FALSE)
    if ( all(degrees[two.rand] > 1) & all(degrees[two.rand] < n - 1) ) {
      flag = TRUE
      break
    }
  }
  if (flag == TRUE) {  # if two different nodes are sampled
    a = two.rand[1]
    b = two.rand[2]
    # add one neighbor to the node which has more neighbors, and
    # substract one neighbor from the node which has less neighbors
    if (degrees[a] >= degrees[b]) { 
      degrees[a] = degrees[a] + 1
      degrees[b] = degrees[b] - 1
    }
    else if (degrees[b] > degrees[a]) {
      degrees[b] = degrees[b] + 1
      degrees[a] = degrees[a] - 1
    }
    return(degrees)
  }
  else {
    warning('Max tried reached, but still can not rewire link to nodes with more neighbors!')
    return(NULL)
  }
}

#' @title generate graphs with different heterogeneity
#' @param n, number of nodes
#' @param k, average degree
#' @param ntried, maximum tried times for link rewiring algorithm, be transfered to function [degree.richer.richer]
#' @param rep1, 
#' @param rep2, number of random graphs generated for the same degree sequence
get.graphs.hetero <- function(n, k, ntried = 1000, rep1 = 10, rep2 = 10) {
  ret = list()
  for (i in 1:rep1) {
    degrees = rep(k, n)  # initialize a regular network
    repeat {
      degrees = degree.richer.richer(degrees, ntried)  # get new degree distribution by link rewiring
      if (is.null(degrees)) break  # if new degree distribution cannot be found, then exit loops
      if (is.graphical.degree.sequence(degrees)) {  # if new degree sequence can be realized in a simple graph
        for (i in 1:rep2) { # generate [rep2] random graphs according to the new degree sequence
          G = degree.sequence.game(out.deg = degrees, method = 'vl')    
          ret[[length(ret) + 1]] = G        
        }
      }
    }    
  }
  ret
}

#' @title get interaction matrix based on the adjacency matrix of graph
#' @param graph, the adjacency matrix of a graph
#' @param beta0, the diagonal elements of interaction matrix
#' @param beta1.min, beta1.max, the minimum and maximum values of off-diagonal elements of interaction matrix
get.interaction.matrix <- function(graph, beta0 = 1, beta1.mean = 0., beta1.sd = 0., beta1.random.type = 'norm') {
  s = dim(graph)[1]  # number of nodes
  edges = sum(graph > 0)  # number of edges
  B = graph
  # generate off-diagonal elements randomly distributed between minimum and maximum values
  if (beta1.random.type == 'norm')
    B[B > 0] = rnorm(edges, mean = beta1.mean, sd = beta1.sd)
  else if (beta1.random.type == 'unif')
    B[B > 0] = runif(edges, min = beta1.mean - beta1.sd, max = beta1.mean + beta1.sd) 
  diag(B) = rep(beta0, s)  # assign the diagonal elements
  B
}
#B = parms.dlv1.rand(graph, beta0 = 1, beta1.min = 0., beta1.max = 0.)  # no interaction
#B = parms.dlv1.rand(graph, beta0 = 1, beta1.min = 0., beta1.max = 0.5)  # competition interaction
#B = parms.dlv1.rand(graph, beta0 = 1, beta1.min = -0.5, beta1.max = 0.)  # cooperation interaction
#B = parms.dlv1.rand(graph, beta0 = 1, beta1.min = -0.5, beta1.max = 0.5)  # arbitrary interaction

#' @title analytically get the output of LV1 model
#' @param r, the vector of intrinsic growth rates
#' @param B, the interaction matrix
#' @return nstar, species abundances at equilibrium
#'         phi, the community matrix (Jacobian matrix) at equilibrium
#'         extinct, survived, the number of extinct species and survived species when system approaches to equilibrium
analysis.lv1 <- function(r, B) {
  nstar = c(solve(B) %*% r)
  phi = - diag(nstar) %*% B
  extinct = sum(nstar <= 0)
  survived = sum(nstar > 0)
  list(nstar = nstar, phi = phi, extinct = extinct, survived = survived)
}

#' @title get covariance matrix of multivariate OU process
mou.vars <- function(phi, C) {
  s = dim(phi)[1]
  I = diag(1, s)
  - matrix(solve(kronecker(I, phi) + kronecker(phi, I)) %*% as.vector(C), nrow = s, ncol = s)
}

get.results.by.variance <- function(graphs, beta1.random.type = 'norm', r.random.type = 'norm', flag = TRUE) {
  r.mean.min = 1
  r.mean.max = 4
  r.mean.step = 1
  r.means = seq(from = r.mean.min, to = r.mean.max, by = r.mean.step)
  beta1.mean.min = -0.1
  beta1.mean.max = 0.2
  beta1.mean.step = 0.1
  beta1.means = seq(from = beta1.mean.min, to = beta1.mean.max, by = beta1.mean.step)
  llply(r.means, function(r.mean) {
    r.sd.min = 0
    r.sd.max = r.mean
    r.sds = seq(from = r.sd.min, to = r.sd.max, length.out = 4)
    llply(r.sds, function(r.sd) {
      llply(beta1.means, function(beta1.mean) {
        beta1.sd.min = 0
        beta1.sd.max = abs(beta1.mean)
        beta1.sds = seq(from = beta1.sd.min, to = beta1.sd.max, length.out = 4)
        llply(beta1.sds, function(beta1.sd) {
          res = get.results.by.graphs(graphs,beta1.mean = 0.1, beta1.sd = 0., beta1.random.type = beta1.random.type, 
                                      r.mean = r.mean, r.sd = r.sd, r.random.type = r.random.type, flag = flag)
          c(r.mean = r.mean, r.sd = r.sd, beta1.mean = beta1.mean, beta1.sd = beta1.sd, res = list(res))          
        })
      })
    })
  })
}

#' @title get simulation results
#' @param graphs, structure of ecological networks with different heteros
#' @param beta1, interaction strength
#' @param flag, does include x^* in the variance-covariance matrix of environmental fluctuations
get.results.by.graphs <- function(graphs, beta1.mean, beta1.sd, beta1.random.type = 'norm', 
                                  r.mean, r.sd, r.random.type = 'norm', flag = TRUE) {
  llply(graphs, .parallel = TRUE, function(graph) {
    A = as.matrix(get.adjacency(graph))  # get the adjacency matrix of graph
    B = get.interaction.matrix(A, beta1.mean = beta1.mean, beta1.sd = beta1.sd, beta1.random.type = 'norm')  # get the interaction matrix from graph and beta1
    s = dim(B)[1]  # number of nodes
    if (r.random.type == 'norm')
      r = rnorm(s, mean = r.mean, sd = r.sd)
    else if (r.random.type == 'unif')
      r = runif(s, min = r.mean - r.sd, max = r.mean + r.sd)  # intrinsic growth rates of species rep(1, s)
    out.lv1 = analysis.lv1(r, B)  # get the output of LV1 model
    ret = list()
    # if the maximum real part of eigenvalues of the community matrix is not negative, then the equilibrium is not local stable
    if (max(Re(eigen(out.lv1$phi, only.values = TRUE)$values)) >= 0) stable = FALSE else stable = TRUE
    hetero = get.degree.hetero(graph)  # 
    if (out.lv1$extinct == 0 && stable) {  # if the equilibrium is feasible (i.e., all species survived) and stable
      nstar = out.lv1$nstar
      phi = out.lv1$phi
      # if consider species abundances in error covariance matrix [C]
      C = diag(1, s)
      if (flag) C = diag(nstar^2) %*% C 
      vars = mou.vars(phi, C)
      ret = list(hetero = hetero, feasible = TRUE, stable = TRUE, nstar = nstar, phi = phi, vars = vars, B = B)
    }
    else if (out.lv1$extinct == 0 && !stable) {  # if the equilibrium is feasible but not stable
      ret = list(hetero = hetero, feasible = TRUE, stable = FALSE)
    }
    else {  # if the equilibrium is neither feasible nor stable
      ret = list(hetero = hetero, feasible = FALSE, stable = FALSE)
    }
    ret         
  })
}
#save(graphs.rewiring, file = paste('graphs.rewiring', s1, s2, k, sep = '-'))

#' @title get [rep] full graphs with [n] nodes
get.graphs.full <- function(n, rep) {
  graph = graph.full(n)
  llply(1:rep, function(i) {
    graph
  })  
}
