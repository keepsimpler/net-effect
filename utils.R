source('simode.r')

#' @title get different structural measures of bipartite networks
#' @param graphs, list of list of graphs which have different structures such as degree heterogeneity
get.structure.measures <- function(graphs) {
  measures = ldply(1:length(graphs), function(i) {
    ldply(graphs[[i]], .parallel = TRUE, function(graph) {
      graph = graph$B
      heterogeneity = get.degree.heterogeneity(graph)
      assortativity = get.assortativity(graph)
      c(heterogeneity = heterogeneity, assortativity = assortativity)
    })
  })
  measures  
}

#' @title check feasibility of ecological networks with different degree heterogeneity
#' @param graphs, list of list of graphs which have different degree heterogeneity
#' @param rep, repeat times for each graph
get.heterogeneity.and.feasible <- function(graphs, rep = 10) {
  heterogeneity.and.feasible = ldply(1:length(graphs), function(i) {
    ldply(graphs[[i]], .parallel = TRUE, function(graph) {
      graph = graph$B
      heterogeneity = get.degree.heterogeneity(graph)
      assortativity = get.assortativity(graph)
      feasible = 0
      for(j in 1:rep) {
        parms = parms.lv2(graph)
        init = init.lv2(parms)      
        A = sim.ode.one(model = model.lv2, parms, init)
        if (A$extinct == 0) feasible = feasible + 1
      }
      c(heterogeneity = heterogeneity, assortativity = assortativity, feasible = feasible)
    })
  })
  heterogeneity.and.feasible
}

get.graphs.out <- function(graphs, rep = 10) {
  graphs.out = llply(1:length(graphs), function(i) {
    llply(graphs[[i]], .parallel = TRUE, function(graph) {
      graph = graph$B
      #heterogeneity = get.degree.heterogeneity(graph)
      #assortativity = get.assortativity(graph)
      llply(1:rep, .parallel = FALSE, function(k) {
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
      })
    })
  })
  graphs.out
}


#' @title get robust measures (tolerance and fragility) for ecological networks with different degree heterogeneity
#' @param graphs, list of list of graphs which have different degree heterogeneity
#' @param rep, repeat times for each graph
get.heterogeneity.and.robust <- function(graphs, rep = 10) {
  heterogeneity.and.tolerances = ldply(1:length(graphs), function(i) {
    ldply(graphs[[i]], .parallel = TRUE, function(graph) {
      graph = graph$B
      heterogeneity = get.degree.heterogeneity(graph)
      assortativity = get.assortativity(graph)
      ldply(1:rep, .parallel = FALSE, function(k) {
        ret = NULL
        parms = parms.lv2(graph)
        init = init.lv2(parms)      
        A = sim.ode.one(model = model.lv2, parms, init)
        if (A$extinct == 0) {
          B = sim.ode(model = model.lv2, parms = parms, init = init, isout = FALSE, iter.steps = 100,
                      perturb = perturb, perturb.type = 'lv2.growth.rate.dec')        
          nstar.init = sum(B[[1]]$nstar)
          tolerances.and.fragility = get.tolerance(B)
          tolerances = tolerances.and.fragility$tolerance.species
          tolerance.abundance = tolerances.and.fragility$tolerance.abundance
          fragility = tolerances.and.fragility$fragility
          tolerance.total = sum(tolerances)
          tolerance.abundance.total = sum(tolerance.abundance)
          ret = c(heterogeneity = heterogeneity, assortativity = assortativity, tolerances = tolerances, tolerance.total = tolerance.total, 
                  tolerance.abundance = tolerance.abundance, tolerance.abundance.total = tolerance.abundance.total, 
                  nstar.init = nstar.init, fragility = fragility)
        }
        ret
      })
    })
  })
  heterogeneity.and.tolerances
} 


#' @title get the tolerance and trend of species and total communities under intrinsic growh rate delining
#' @param A, the list of ODE output
get.tolerance <- function(A) {
  for( i in 1:length(A)) {  # replace NA with 0
    A[[i]]$nstar[is.na(A[[i]]$nstar)] = 0
  }
  steps = length(A)
  if (sum(A[[steps]]$nstar > 0) > 0) {
    warning('Some species still exist.')
  }
  if (sum(A[[1]]$nstar == 0) > 0) {
    warning('System is not at feasible equilibrium at beginning.')
  }
  n = length(A[[1]]$nstar)  # number of species at beginning
  tolerance.species = rep(0, n)  # initialize
  fragility = 0
  nstar.pre = A[[1]]$nstar  # species abundance of previous step
  nstar.pre.num = sum(nstar.pre > 0)  # species number of previous step
  tolerance.abundance = nstar.pre # 
  for (i in 2:steps) {
    nstar = A[[i]]$nstar  # species abundance of this step
    nstar.num = sum(nstar > 0)  # species number of this step
    # get the summed abundance accross steps
    tolerance.abundance = tolerance.abundance + nstar
    # get trend of species loss
    species.loss = nstar.pre.num - nstar.num
    if (species.loss > 0) {
      fragility = fragility + species.loss * log(species.loss)      
    }
    # chose the species extincted in this step but still exist in previous step, i.e., the tolerance of species
    species.extinct = which(nstar.pre - nstar == nstar.pre & nstar.pre >0)
    if (length(species.extinct) > 0) {  #
      tolerance.species[species.extinct] = i
    }
    if (sum(nstar > 0) == 0) {
      break
    }
    else {
      nstar.pre = nstar
      nstar.pre.num = nstar.num
    }
  }
  list(tolerance.species = tolerance.species, fragility = fragility, tolerance.abundance = tolerance.abundance)
}

#' @title determine if a ODE system is feasible
is.feasible.lv2 <- function(graph, ...) {
  parms = parms.lv2(graph)
  init = init.lv2(parms)
  
  A = sim.ode.one(model = model.lv2, parms, init)
  if (A$extinct == 0) return(TRUE)
  else return(FALSE)
}

#' @title get assortativity (degree-degree corelation), <number of neighbors of neighbors k_2> / <k^2>
#' @param graph, the incidence matrix of bipartite network
get.assortativity <- function(graph) {
  graph[graph != 0] = 1
  degrees2 = c(rowSums(graph %*% t(graph)), rowSums(t(graph) %*% graph)) # two-hop degrees (degrees of neighbors of node)
  degrees = c(rowSums(graph), colSums(graph))
  assortativity = sum( degrees2 / degrees^2 )
  assortativity
}

#' @title get degree heterogeneity, <k^2>/(<k>^2)
get.degree.heterogeneity <- function(graph) {
  graph[graph != 0] = 1
  degrees = c(rowSums(graph), colSums(graph))
  heterogeneity = sum(degrees^2) / sum(degrees)^2
  heterogeneity
}

#' @title generate random modular bipartite graph with different levels of degree heterogeneity
#' @param n1, n2, number of nodes in ROW and COL parts respectively
#' @param modules.row, the members of modules for ROW part
#' @param modules.col, the members of modules for COL part
#' @param degrees.row, node degrees for ROW part
#' @param degrees.col, node degrees for COL part
#' @param Q, expected modularity
gen.modularity <- function(n1, n2, modules.row, modules.col, degrees.row, degrees.col, Q) {
  
}
#' @title measure the modularity of a bipartite network
#' @param graph, the incidence matrix of a bipartite network
#' @param modules.row, the members of modules for ROW part
#' @param modules.col, the members of modules for COL part
get.modularity <- function(graph, modules.row, modules.col) {
  stopifnot(length(modules.row) == length(modules.col))  # assert number of modules for ROW and COL parts are equal
  L = sum(graph)  # total number of links of the bipartite network
  degrees.row = rowSums(graph)  # degrees of node for ROW part
  degrees.col = colSums(graph)  # degrees of node for COL part
  modularity = 0
  for (i in 1:length(modules.row)) {
    m1 = modules.row[[i]]$id  # member list of corresponding module for ROW part
    m2 = modules.col[[i]]$id  # member list of corresponding module for COL part
    module = graph[m1, m2 - nrow(graph)]  # the subgraph of the module
    Lm = sum(module)  # the number of link within module
    expected = sum(degrees.row[m1]) * sum(degrees.col[m2 - nrow(graph)])
    modularity = modularity + Lm / L - expected / L^2
    print(paste(Lm / L, expected / L^2))
  }
  modularity
}

#' @title detect the modules of empirical bipartite networks
#' @param graph, the incidence matrix of a bipartite network
#' @param num.row, number of modules for ROW part
#' @param num.col, number of modules for COL part
#' @param iter, number of iterations
get.modules <- function(graph, num.row, num.col, iter = 100) {
  edges = get.edgelist(graph.incidence(graph))  # get the edge list of bipartite networks
  edges = data.frame(edges)
  types = c(rep(1, nrow(graph)), rep(2, ncol(graph)))  # get the type (part) of nodes belong to
  types = data.frame(types)
  modules = biSBM(data = edges, nodeType = types, ka = num.row, kb = num.col, deg.corr = 1, iter = iter)
  modules = data.frame(id = 1:length(modules), modules)
  members = split(modules, f = modules$modules)  # {1, 2, 3, 2, 3} --> {{1, {1}}, {2,{2, 4}}, {3,{3, 5}}}
  members.row = members[1:num.row]
  members.col = members[(num.row+1):(num.row+num.col)]
  list(members.row, members.col)
}




#' @title generate random bipartite graphs with diferent levels modularity and degree heterogeneity
#' 
# degrees.row = rowSums(graph)
# degrees.col = colSums(graph)
# degrees.out = append(degrees.row, rep(0, length(degrees.col)))
# degrees.in = append(degrees.col, after = 0, rep(0, length(degrees.row)))
# g = degree.sequence.game(out.deg = degrees.out, in.deg = degrees.in)
# g = rewire(g, niter = 1000)  # elimite the multiple edges by rewiring links
# is.multiple(g)

