#!/usr/bin/Rscript
library(rootSolve)  # for the Jacobian matrix in equilibrium, function [jacobian.full]
library(deSolve)
#library(simecol)  # for the simulation of ODE, which use the <deSolve> package
require(plyr)
require(bipartite)

library(doMC)  # 
registerDoMC()  # register Multi Cores
getDoParWorkers()  # get available Cores


extinct.threshold.default = .Machine$double.eps * 100  # threshold of species extinction is 100 times of machine precision
extinct.threshold = extinct.threshold.default  # 1e-3


#' @title Lotka-Volterra (LV) Equations of Holling type II by Bastolla et al. for mutualistic communities
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
  list(c(dN))
}

#' @title the parmaters of mutualistic LV2 model
#' @param graph, the incident matrix of mutualistic networks which are bipartite
#' @param alpha.row.mu, alpha.row.sd, the intrinsic growth rate of Plants (i.e., the ROW part)
#' @param alpha.col.mu, alpha.col.sd, the intrinsic growth rate of Animals (i.e., the COL part)
#' @param beta0.mu, beta0.sd, the intraspecies competition
#' @param beta1.mu, beta1.sd, the interspecies competition
#' @param gamma.mu, gamma.sd, the interspecies cooperation
#' @param h.mu, h.sd, the Handling time, saturate coefficient
#' @param delta, trade-off between pairwise mutualistic interaction strengths
#' @return a list of parameters
parms.lv2 <- function(graph, alpha.row.mu = 0.2, alpha.row.sd = 0.15, alpha.col.mu = 0.2, alpha.col.sd = 0.15, 
                      beta0.mu = 1, beta0.sd = 0.2, beta1.mu = 0., beta1.sd = 0.0,
                      gamma.mu = 1., gamma.sd = 0.2, h.mu = 0.1, h.sd = 0.05, delta = 1) {
  numP = dim(graph)[1]
  numA = dim(graph)[2]
  s = numP + numA
  r = c(runif(numP) * 2 * alpha.row.sd + (alpha.row.mu - alpha.row.sd), runif(numA) * 2 * alpha.col.sd + (alpha.col.mu - alpha.col.sd))
  
  C = matrix(0, nrow = s, ncol = s)
  C[1:numP, 1:numP] = runif(numP * numP) * 2 * beta1.sd + (beta1.mu - beta1.sd)
  C[(numP+1):s, (numP+1):s] = runif(numA * numA) * 2 * beta1.sd + (beta1.mu - beta1.sd)
  diag(C) = runif(s) * 2 * beta0.sd + (beta0.mu - beta0.sd)
  
  edges = sum(graph > 0)  # the number of edges
  M = as.one.mode(graph)  # transform to adjacency matrix (function of package [bipartite])
  degrees = rowSums(M)
  M[M > 0] = runif(2 * edges) * 2 * gamma.sd + (gamma.mu - gamma.sd)  # endue values of interspecies cooperation
  M = M / degrees^delta  # trade-off of mutualistic strength    
  
  h = runif(s) * 2 * h.sd + (h.mu - h.sd)
  N = rep(1, s)
  r2 = C %*% N - (M %*% N) / (1 + h * M %*% N)
  list(r = r, C = C, M = M, h = h, r2 = r2)  # the [parms] of ode model
}

init.lv2 <- function(parms) {
  init = solve(parms$C - parms$M) %*% parms$r
  if (any(init < 0)) {
    print('Initial state values is less than 0 !!')
    #stop('Initial state values is less than 0 !!', init(LV2))
    init = parms$r
  }
  init  
}


#' @title One simulation of ODE dynamics
sim.ode.one <- function(model, parms, init, steps = 1000, stepwise = 1) {
  times = seq(from = 0, to = steps * stepwise, by = stepwise)
  ode.out = ode(init, times, model, parms)
  nstar = as.numeric(ode.out[nrow(ode.out), 2:ncol(ode.out)])
  nstar[nstar < extinct.threshold] = 0  # species with biomass less than the threshold is considered to be extinct
  extinct = length(nstar) - sum(nstar > 0) 
  survived = sum(nstar > 0)
  Phi = jacobian.full(y = nstar, func = model, parms = parms)
  ret = list(out = ode.out, nstar = nstar, Phi = Phi, model = model, parms = parms, extinct = extinct, survived = survived)
  ret
}

#' @title iterately change parameters of ODE model, such as decrease intrinsic growth rate
#' @param model, parms, init
sim.ode <- function(model, parms, init, steps = 1000, stepwise = 1, isout = TRUE, 
                    iter.steps = 10, perturb, perturb.type = 'lv2.growth.rate.dec') {
  times = seq(from = 0, to = steps * stepwise, by = stepwise)
  ode.outs = list()
  for(i in 1:iter.steps) {
    print(i)
    ode.out = ode(init, times, model, parms) 
    nstar = as.numeric(ode.out[nrow(ode.out), 2:ncol(ode.out)]) # species biomass at equilibrium
    nstar[nstar < extinct.threshold] = 0  # species with biomass less than extinct threshold is considered to be extinct
    extinct.species = which(nstar == 0)  # extinct species

    Phi = jacobian.full(y = nstar, func = model, parms = parms) # community matrix, Jacobian matrix at equilibrium
    if (isout) {
      ret = list(out = ode.out, nstar = nstar, Phi = Phi, params = parms, extinct.species = extinct.species)
    }
    else {
      ret = list(nstar = nstar, Phi = Phi, params = parms, extinct.species = extinct.species)
    }
    ode.outs[[length(ode.outs) + 1]] = ret
    
#     if (length(extinct.species) > 0) {
#       ret = remove.species(parms, nstar, extinct.species)
#       parms = ret$parms
#       nstar = ret$nstar
#     }
#     if (length(nstar) == 0) break  # if all species are extinct, then stop and quit
    
    perturb.res = perturb(parms, nstar, perturb.type)
    parms = perturb.res$parms
    init = perturb.res$nstar
  }
  ode.outs
}

remove.species <- function(parms, nstar, extinct.species) {
  if (length(extinct.species) > 0) {
    nstar = nstar[- extinct.species]
    parms$r = parms$r[- extinct.species]
    parms$C = parms$C[- extinct.species, - extinct.species]
    parms$M = parms$M[- extinct.species, - extinct.species]
    parms$h = parms$h[- extinct.species]
  }
  list(parms = parms, nstar = nstar)  
}

perturb <- function(parms, nstar, perturb.type, numP = NULL, numA = NULL) {
  if (perturb.type == 'lv2.growth.rate.dec') {
    parms$r = parms$r - 0.01  #runif(length(nstar), min =  0.01, max = 0.02)
  }
  else if(perturb.type == 'lv2.growth.rate.dec.onepart') {
    numP = 16; numA = 25;
    #parms$r[1:numP] = parms$r[1:numP] - 0.02
    parms$r[(numP + 1):(numP+numA)] = parms$r[(numP + 1):(numP+numA)] - 0.02
  }
  else if (perturb.type == 'lv2.primary.extinction') {
    nstar[1] = 0
    #nstar[2] = 0
  }
  else if (perturb.type == 'scheffer.nutrient.inc') {
    parms['N'] = parms['N'] + 0.1
    nstar = init.scheffer()
  }
  list(parms = parms, nstar = nstar)
}


#' @title get the sensitivity matrix
#' @param nstar, species abundance in equilibrium
#' @param Phi, the community matrix in equilibrium
get.sensitivity <- function(nstar, Phi) {
  # delete extinct species, i.e., species with abundance zero
  extinct.species = which(nstar == 0)
  if (length(extinct.species) > 0) {
    nstar = nstar[- extinct.species]
    Phi = Phi[- extinct.species, - extinct.species]    
  }
  sensitivity.matrix = diag(1/nstar) %*% - solve(Phi) %*% diag(nstar)  #%*% diag(1/nstar) plot(nstar, diag(sensitivity))
  #sensitivity.matrix.diag = diag(1/diag(sensitivity.matrix))
  #sensitivity.matrix = sensitivity.matrix %*% sensitivity.matrix.diag
  sensitivity = rowSums(sensitivity.matrix)
  
  if(length(extinct.species) > 0) {
    for(i in 1:length(extinct.species)) 
      sensitivity = append(sensitivity, NaN, extinct.species[i] - 1)
  }
  
  lev = NaN # max(Re(eigen(Phi)$values))
  list(extinct.species = extinct.species, sensitivity = sensitivity, lev = lev)
}

#' @title Lotka-Volterra (LV) model with functional response of Holling Type I
model.lv1 <- function(time, init, parms) {
  r = parms[[1]]  # intrinsic growth rate
  C = parms[[2]]  # the competition matrix
  N = init  # initial state
  dN <- N * ( r - C %*% N )
  list(c(dN))
}

#' @title parameters for random interactions 
#' @param s, number of nodes
#' @param k, node degree, connectance = k / s
#' @param delta, reflect interaction strengths
parms.lv1.rand <- function(s, k, delta) {
  g = graph.connected(s = s, k = k, gtype = 'sf')
  A = as.matrix(get.adjacency(g))
  A[A > 0] = rnorm(n = s * k * 2, mean = 0, sd = delta)
  diag(A) = 1
  N = rep(1, s)
  r = A %*% N
  list(r = r, C = A)
}

init.lv1 <- function(parms) {
  init = solve(parms$C) %*% parms$r
  if (any(init < 0)) {
    print('Initial state values is less than 0 !!')
    init = parms$r
  }
  init  
}

analysis.lv1 <- function(r, B) {
  nstar = c(solve(B) %*% r)
  phi = - diag(nstar) %*% B
  extinct = sum(nstar <= 0)
  survived = sum(nstar > 0)
  list(nstar = nstar, phi = phi, extinct = extinct, survived = survived)
}

mou.vars <- function(phi, C) {
  s = dim(phi)[1]
  I = diag(1, s)
  - matrix(solve(kronecker(I, phi) + kronecker(phi, I)) %*% as.vector(C), nrow = s, ncol = s)
}

