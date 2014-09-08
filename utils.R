source('simode.r')

#' @title determine if a ODE system is feasible
is.feasible.lv2 <- function(graph, ...) {
  parms = parms.lv2(graph)
  init = init.lv2(parms)
  
  A = sim.ode.one(model = model.lv2, parms, init)
  if (A$extinct == 0) return(TRUE)
  else return(FALSE)
}

get.degree.heterogeneity <- function(graph) {
  graph[graph != 0] = 1
  degrees = c(rowSums(graph), colSums(graph))
  heterogeneity = sum(degrees^2) / sum(degrees)^2
  heterogeneity
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
    m1 = modules.row[[i]]$id  # member list of module for ROW part
    m2 = modules.col[[i]]$id  # member list of module for COL part
    module = graph[m1, m2 - nrow(graph)]  # the subgraph of the module
    Lm = sum(module)  # the number of link within module
    expected = sum(degrees.row[m1]) * sum(degrees.col[m2 - nrow(graph)])
    modularity = modularity + Lm / L - expected / L^2
    cat( Lm / L, expected / L^2)
  }
  modularity
}

#' @title detect the modules of empirical networks
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
  members = split(modules, f = modules$modules)
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

