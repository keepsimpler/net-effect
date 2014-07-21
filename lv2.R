#!/usr/bin/Rscript
library(rootSolve)  # for the Jacobian matrix in equilibrium, function [jacobian.full]
library(simecol)  # for the simulation of ODE, which use the <deSolve> package
require(plyr)
require(bipartite)

library(doMC)  # 
registerDoMC()  # register Multi Cores
getDoParWorkers()  # get available Cores

extinct.threshold.default = .Machine$double.eps * 100  # threshold of species extinction is 100 times of machine precision

#' @title Lotka-Volterra (LV) Equations of Holling type II by Bastolla et al.
#' @param time, time step of simulation
#' @param init, the initial state of the LV system, a vector
#' @param parms parameters passed to LV model
#'        r, the intrinsic growth rate of species, a vector
#'        C, the competition matrix in plants and animals
#'        M, the cooperation matrix between plants and animals
#' @return the derivation
#' @details .
#' @import deSolve
model.lv2 <- function(time, init, parms, ...) {
  r = parms[[1]]  # intrinsic growth rate
  C = parms[[2]]  # the competition matrix
  M = parms[[3]]  # the cooperation matrix
  h = parms[[4]]  # handling time
  N = init  # initial state
  dN <- N * ( r - C %*% N + (M %*% N) / (1 + h * M %*% N) )
  #dN = dN - 0.01
  list(c(dN))
}


#' @title the [parms] of mutualistic lv2 model
#' @param graph, the incident matrix of mutualistic networks which are bipartite
#' @param alpha.mu, alpha.sd, the intrinsic growth rate
#' @param beta0.mu, beta0.sd, the intraspecies competition
#' @param beta1.mu, beta1.sd, the interspecies competition
#' @param gamma.mu, gamma.sd, the interspecies cooperation
#' @param h.mu, h.sd, the Handling time, saturate coefficient
#' @return [parms] of [simObj] class
parms.lv2 <- function(graph, alpha.mu = 0.2, alpha.sd = 0.15, beta0.mu = 1, beta0.sd = 0.2, beta1.mu = 0.03, beta1.sd = 0.02,
                      gamma.mu = 1, gamma.sd = 0.2, h.mu = 0.2, h.sd = 0.1) {
  numP = dim(graph)[1]
  numA = dim(graph)[2]
  s = numP + numA
  r = runif(s) * 2 * alpha.sd + (alpha.mu - alpha.sd)
  C = matrix(runif(s * s), ncol = s, nrow = s) * 2 * beta1.sd + (beta1.mu - beta1.sd)
  diag(C) = runif(s) * 2 * beta0.sd + (beta0.mu - beta0.sd)
  
  edges = sum(graph > 0)  # the number of edges
  M = as.one.mode(graph)  # transform to adjacency matrix (function of package [bipartite])
  M[M > 0] = runif(2 * edges) * 2 * gamma.sd + (gamma.mu - gamma.sd)  # endue values of interspecies cooperation
  
  h = runif(s) * 2 * h.sd + (h.mu - h.sd)
  list(r = r, C = C, M = M, h = h)  # the [parms] of [simObj] class in package [simecol]
}

#' @title LV2 simulation function
#' @param graph, the incidence matrix of graph.
#' @param alpha0,
#' @param beta0,
#' @param gamma0,
#' @param h0,
#' @param isout, if output the time serials of simulation
#' @param steps and stepwise of simulation
sim.lv2.graph <- function(graph, alpha.mu = 0.2, alpha.sd = 0.15, beta0.mu = 1, beta0.sd = 0.2, beta1.mu = 0.03, beta1.sd = 0.02,
                          gamma.mu = 1, gamma.sd = 0.2, h.mu = 0.2, h.sd = 0.1, Xinit = NULL, 
                          isout = FALSE, steps = 10000, stepwise = 0.01) {
  LV2 <- odeModel(
    main = model.lv2, 
    times = c(from = 0, to = steps * stepwise, by = stepwise),
    solver = 'lsoda')

  parms(LV2) = parms.lv2(graph, alpha.mu = alpha.mu, alpha.sd = alpha.sd, beta0.mu = beta0.mu, beta0.sd = beta0.sd,
                             beta1.mu = beta1.mu, beta1.sd = beta1.sd,
                             gamma.mu = gamma.mu, gamma.sd = gamma.sd, h.mu = h.mu, h.sd = h.sd)
  if (is.null(Xinit)) {
    Xinit = solve(parms(LV2)$C - parms(LV2)$M) %*% parms(LV2)$r  # the [init] Value, which is close to the steady state.    
  }
  if (any(Xinit < 0)) {
    print('Initial state values is less than 0 !!')
    #stop('Initial state values is less than 0 !!', init(LV2))
    Xinit = r
  }
  init(LV2) = Xinit
  
  LV2 <- sim(LV2)
  LV2.out = out(LV2)
  
  Nstar = as.numeric(LV2.out[nrow(LV2.out), 2:ncol(LV2.out)]) 
  
  Phi = jacobian.full(y = Nstar, func = model.lv2, parms = parms(LV2))
  if (isout) {
    ret = list(out = LV2.out, Nstar = Nstar, Phi = Phi, params = parms(LV2))
  }
  else {
    ret = list(Nstar = Nstar, Phi = Phi, params = parms(LV2))
  }
  ret
}

## Decrease the intrinsic growth rates
sim.lv2.alpha.dec <- function(dec.steps = 10, dec.stepwise = 0.01, graph, alpha.mu = 0.2, alpha.sd = 0.15, beta0.mu = 1, beta0.sd = 0.2, beta1.mu = 0.03, beta1.sd = 0.02,
                              gamma.mu = 1, gamma.sd = 0.2, h.mu = 0.2, h.sd = 0.1, Xinit = NULL, 
                              isout = FALSE, steps = 10000, stepwise = 0.01) {
  lv2.outs = list()
  
  LV2 <- odeModel(
    main = model.lv2, 
    times = c(from = 0, to = steps * stepwise, by = stepwise),
    solver = 'lsoda')
  
  parms(LV2) = parms.lv2(graph, alpha.mu = alpha.mu, alpha.sd = alpha.sd, beta0.mu = beta0.mu, beta0.sd = beta0.sd,
                         beta1.mu = beta1.mu, beta1.sd = beta1.sd,
                         gamma.mu = gamma.mu, gamma.sd = gamma.sd, h.mu = h.mu, h.sd = h.sd)
  if (is.null(Xinit)) {
    Xinit = solve(parms(LV2)$C - parms(LV2)$M) %*% parms(LV2)$r  # the [init] Value, which is close to the steady state.    
  }
  if (any(Xinit < 0)) {
    print('Initial state values is less than 0 !!')
    #stop('Initial state values is less than 0 !!', init(LV2))
    Xinit = parms(LV2)$r 
  }
  
  for ( i in 1:dec.steps) {
    init(LV2) = Xinit

    LV2 <- sim(LV2)
    LV2.out = out(LV2)
    
    Nstar = as.numeric(LV2.out[nrow(LV2.out), 2:ncol(LV2.out)]) 
    Nstar[Nstar < 10^-10] = 0  # species with biomass less than the threshold is considered to be extinct
    
    Phi = jacobian.full(y = Nstar, func = model.lv2, parms = parms(LV2))
    if (isout) {
      ret = list(out = LV2.out, Nstar = Nstar, Phi = Phi, params = parms(LV2))
    }
    else {
      ret = list(Nstar = Nstar, Phi = Phi, params = parms(LV2))
    }
    lv2.outs[[length(lv2.outs) + 1]] = ret
    #parms(LV2)$r[1] = parms(LV2)$r[1] - dec.stepwise
    parms(LV2)$r = parms(LV2)$r - dec.stepwise
    #Nstar = Nstar - 0.1
    #if (Nstar[1] < 0) Nstar[1] = 0
    Xinit = Nstar  # ret$Nstar
  }
  lv2.outs
}


## Decrease the intrinsic growth rates
sim.lv2.press.dec <- function(dec.steps = 10, dec.stepwise = 0.01, graph, alpha.mu = 0.2, alpha.sd = 0.15, beta0.mu = 1, beta0.sd = 0.2, beta1.mu = 0.03, beta1.sd = 0.02,
                              gamma.mu = 1, gamma.sd = 0.2, h.mu = 0.2, h.sd = 0.1, Xinit = NULL, 
                              isout = FALSE, steps = 10000, stepwise = 0.01) {
  lv2.outs = list()
  
  LV2 <- odeModel(
    main = model.lv2, 
    times = c(from = 0, to = steps * stepwise, by = stepwise),
    solver = 'lsoda')
  
  parms(LV2) = parms.lv2(graph, alpha.mu = alpha.mu, alpha.sd = alpha.sd, beta0.mu = beta0.mu, beta0.sd = beta0.sd,
                         beta1.mu = beta1.mu, beta1.sd = beta1.sd,
                         gamma.mu = gamma.mu, gamma.sd = gamma.sd, h.mu = h.mu, h.sd = h.sd)
  if (is.null(Xinit)) {
    Xinit = solve(parms(LV2)$C - parms(LV2)$M) %*% parms(LV2)$r  # the [init] Value, which is close to the steady state.    
  }
  if (any(Xinit < 0)) {
    print('Initial state values is less than 0 !!')
    #stop('Initial state values is less than 0 !!', init(LV2))
    Xinit = parms(LV2)$r 
  }
  
  for ( i in 1:dec.steps) {
    init(LV2) = Xinit
    
    LV2 <- sim(LV2)
    LV2.out = out(LV2)
    
    Nstar = as.numeric(LV2.out[nrow(LV2.out), 2:ncol(LV2.out)]) 
    Nstar[Nstar < 10^-10] = 0  # species with biomass less than the threshold is considered to be extinct
    
    Phi = jacobian.full(y = Nstar, func = model.lv2, parms = parms(LV2))
    if (isout) {
      ret = list(out = LV2.out, Nstar = Nstar, Phi = Phi, params = parms(LV2))
    }
    else {
      ret = list(Nstar = Nstar, Phi = Phi, params = parms(LV2))
    }
    lv2.outs[[length(lv2.outs) + 1]] = ret
    #parms(LV2)$r[1] = parms(LV2)$r[1] - dec.stepwise
    #parms(LV2)$r = parms(LV2)$r - dec.stepwise
    #Nstar = Nstar - 0.1
    #if (Nstar[1] < 0) Nstar[1] = 0
    Xinit = parms(LV2)$r  # ret$Nstar
  }
  lv2.outs
}


## Decrease the interspecies cooperation
sim.lv2.gamma.dec <- function(dec.steps = 10, dec.stepwise = 0.01, graph, alpha.mu = 0.2, alpha.sd = 0.15, beta0.mu = 1, beta0.sd = 0.2, beta1.mu = 0.03, beta1.sd = 0.02,
                              gamma.mu = 1, gamma.sd = 0.2, h.mu = 0.2, h.sd = 0.1, Xinit = NULL, 
                              isout = FALSE, steps = 10000, stepwise = 0.01) {
  lv2.outs = list()
  
  LV2 <- odeModel(
    main = model.lv2, 
    times = c(from = 0, to = steps * stepwise, by = stepwise),
    solver = 'lsoda')
  
  parms(LV2) = parms.lv2(graph, alpha.mu = alpha.mu, alpha.sd = alpha.sd, beta0.mu = beta0.mu, beta0.sd = beta0.sd,
                         beta1.mu = beta1.mu, beta1.sd = beta1.sd,
                         gamma.mu = gamma.mu, gamma.sd = gamma.sd, h.mu = h.mu, h.sd = h.sd)
  if (is.null(Xinit)) {
    Xinit = solve(parms(LV2)$C - parms(LV2)$M) %*% parms(LV2)$r  # the [init] Value, which is close to the steady state.    
  }
  if (any(Xinit < 0)) {
    print('Initial state values is less than 0 !!')
    #stop('Initial state values is less than 0 !!', init(LV2))
    Xinit = parms(LV2)$r 
  }
  
  for ( i in 1:dec.steps) {
    init(LV2) = Xinit
    
    LV2 <- sim(LV2)
    LV2.out = out(LV2)
    
    Nstar = as.numeric(LV2.out[nrow(LV2.out), 2:ncol(LV2.out)]) 
    Nstar[Nstar < 10^-10] = 0  # species with biomass less than the threshold is considered to be extinct
    extinct = length(Nstar > 0) / length(Nstar)
    cat(extinct)
    
    Phi = jacobian.full(y = Nstar, func = model.lv2, parms = parms(LV2))
    if (isout) {
      ret = list(out = LV2.out, Nstar = Nstar, Phi = Phi, params = parms(LV2))
    }
    else {
      ret = list(Nstar = Nstar, Phi = Phi, params = parms(LV2))
    }
    lv2.outs[[length(lv2.outs) + 1]] = ret
    parms(LV2)$M = parms(LV2)$M - dec.stepwise
    Xinit = ret$Nstar
  }
  lv2.outs
}



  #' @title the [parms] and [init] of mutualistic lv2 model in soft mean field case
#' @param graph, the incident matrix of mutualistic networks which are bipartite
#' @param alpha0, the intrinsic growth rate
#' @param beta0, the mean value of intraspecies competition,
#' @param gamma0, the mean value of interspecies cooperation
#'        which is endued according to the condition of feasible equilibrium
#' @param h0, the Handling time, saturate coefficient.
#' @param xinit, the initialized state values (species abundance)
#' @return a list of [parms] and [init] values of [simObj] class
parms.lv2.softmean <- function(graph, alpha0 = 1, beta0 = 1, gamma0 = NULL, h0 = 0, Xinit = NULL) {
  numP = dim(graph)[1]; numA = dim(graph)[2]
  r = rep(alpha0, numP + numA)  # endue the mean value of intrinsic growth rate
  edges = sum(graph > 0)  # the number of edges
  # [gamma0] take value of squared root of edge number, to ensure the positive definitive of M
  # and thus ensure the only feasible equilibrium.
  if (is.null(gamma0))
    gamma0 = 1 / ceiling( sqrt(edges) )
  C = diag( rep( beta0, numP + numA ) )  # endue the mean value of intraspecies competition
  M = as.one.mode(graph)  # transform to adjacency matrix (function of package [bipartite])
  M[M > 0] = gamma0  # endue the mean value of interspecies cooperation
  h = rep(h0, numP + numA)  # endue the mean value of handling time
  parmsV = list(r = r, C = C, M = M, h = h)  # the [parms] Value of [simObj] class in package [simecol]
  if (is.null(Xinit)) {
    Xinit = solve(C - M) %*% r  # the [init] Value, which is close to the steady state.    
  }
  if (any(Xinit < 0)) {
    print('Initial state values is less than 0 !!')
    Xinit = r
  }
  list(parmsV = parmsV, initV = Xinit)
}

#' @title LV2 simulation function
#' @param graphs, a list of matrices which are incidence matrices of graphs.
#' @param alpha0,
#' @param beta0,
#' @param gamma0,
#' @param h0,
#' @param isout, if output the time serials of simulation
#' @param steps and stepwise of simulation
sim.lv2 <- function(graphs, alpha0 = 1, beta0 = 1, gamma0 = NULL, h0 = 0.01, isout = FALSE, steps = 10000, stepwise = 0.01) {
  LV2 <- odeModel(
    main = model.lv2, 
    times = c(from = 0, to = steps * stepwise, by = stepwise),
    solver = 'lsoda')
  
  graphs.num = length(graphs)  # number of graphs
  result.lv2 = llply(1:graphs.num, .parallel = TRUE, function(i) {
    A = graphs[[i]]$B  # get the incidence matrix of some graph
    print(i)
    parms.and.init = parms.lv2.softmean(A, alpha0 = alpha0, beta0 = beta0, gamma0 = gamma0, h0 = h0)
    parms(LV2) = parms.and.init[[1]]
    init(LV2) = as.numeric(parms.and.init[[2]])
    if (any(init(LV2) < 0)) {
      print('Initial state values is less than 0 !!')
      stop('Initial state values is less than 0 !!', init(LV2))
    }
    
    LV2 <- sim(LV2)
    LV2.out = out(LV2)
    
    Nstar = as.numeric(LV2.out[nrow(LV2.out), 2:ncol(LV2.out)]) 
    
    Phi = jacobian.full(y = Nstar, func = model.lv2, parms = parms(LV2))
    if (isout) {
      ret = list(out = LV2.out, Nstar = Nstar, Phi = Phi)
    }
    else {
      ret = list(Nstar = Nstar, Phi = Phi)
    }
    ret
    #ret = as.matrix(LV2.out, ncol = length(LV2.out))
    #out = LV2.out[2:length(LV2.out)]
  })
  result.lv2
}

#' @title get the max cooperation strength [gamma0], 
#'        that allow a feasible steady state exists, for the soft mean field case
#' @param graph, the incidence matrix of mutualistic networks
#' @param beta0, the intraspecies competition strength
#' @param sigma0, the SD of white noise
get.gamma0.max <- function(graph, beta0 = 1, sigma0 = 0.01) {
  edges = sum(graph > 0)
  numP = dim(graph)[1]
  numA = dim(graph)[2]
  D = diag(beta0, numP + numA)
  A = as.one.mode(graph)
  gamma0 = 1 / sqrt(edges) 
  repeat {
    gamma0 = gamma0 + 0.0001
    A[A > 0] = gamma0
    Theta = D - A  # competition interaction matrix
    if (min(eigen(Theta)$values) <= min( sigma0, 0.01) ) {
      gamma0 = gamma0 - 0.0001
      break
    }
  }
  gamma0
}



sim.lv1 <- function(graphs, alpha0 = 1, beta0 = 1, gamma0 = NULL) {
  result.lv1 = llply(graphs, function(graph) {
    graph = graph$B
    edges = sum(graph > 0)
    numP = dim(graph)[1]
    numA = dim(graph)[2]
    if (is.null(gamma0)) {
      gamma0 = 1 / ceiling( sqrt(edges) )
    }
    D = diag(beta0, numP + numA)
    A = as.one.mode(graph)
    A[A > 0 ] = gamma0
    M = D - A  # competition interaction matrix
    #r = runif(numP + numA, min = alpha0, max = alpha0)
    r = rep(alpha0, numP + numA)  # endue the mean value of intrinsic growth rate
    Nstar = solve(M) %*% r  # the feasible fixed point
    Phi = - M * as.vector(Nstar)  # the community matrix
    list(Nstar = Nstar, Phi = Phi)
  })
  result.lv1
}

sim.lv1.2 <- function(graphs, alpha0 = 1, beta0 = 0.1, gamma0 = 1) {
  result.lv1 = llply(graphs, function(graph) {
    graph = graph$B
    edges = sum(graph > 0)
    numP = dim(graph)[1]
    numA = dim(graph)[2]
    D = diag(1, numP + numA)
    D[1:numP, 1:numP] = beta0
    D[(numP+1):(numP+numA), (numP+1):(numP+numA)] = beta0
    diag(D) = 1
    A = as.one.mode(graph)
    A[A > 0 ] = gamma0
    M = D - A  # competition interaction matrix
    #r = runif(numP + numA, min = alpha0, max = alpha0)
    r = rep(alpha0, numP + numA)  # endue the mean value of intrinsic growth rate
    Nstar = solve(M) %*% r  # the feasible fixed point
    Phi = - M * as.vector(Nstar)  # the community matrix
    list(Nstar = Nstar, Phi = Phi)
  })
  result.lv1
}


#sum(sim.lv1.2(graphs=list(graph1), gamma0=0.2, beta0=0.03)[[1]]$Nstar)

# save(result.lv2, file = 'result.lv2.RData')
