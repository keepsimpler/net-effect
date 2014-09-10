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

#' @references <Rescuing ecosystems from extinction cascades through compensatory perturbations. Motter.>
model.foodweb <- function(time, init, parms) {
  
}

#' @param graph, the binary adjacency matrix of predators and preys
#' @param r, the mass-specific growth rate
#' @param K, the carrying capacity
#' @param x, the mass-specific metabolic rate relative to mass-specific growth rate
#' @param y, the mass-specific ingestion rate relative to the mass-specific metabolic rate
#' @param e, e_ji is the assimilation efficiency of species j when feeding on i
#' @param m, the average body mass of the species increases with the trophic level as m_i = Z^{τ−1}, 
#'           where τ denotes the trophic level and Z is a constant.
#' @param Z, an emperience value, 10
#' @param dx.dr, the ratio d_x / d_r
#' @param dy.dx, the ratio d_y / d_x
#' @param w, the intraspecies competition between individuals of the same predator when they prey
#' @param h, the exponent of functional responses, accounting for different types, such as Holling Type I, II, III
#' 
parms.foodweb <- function(graph, r, K, x, y, e, m) {
  
}

#' @param graph, the adjacency matrix of food web
parms.lv1.foodweb <- function(graph, e = 0.1) {
  s = dim(graph)[1]
  species.basal = which(colSums(graph) == 0)
  species.top = which(rowSums(graph) == 0)
  r = runif(s, min = -1, max = 0) # the non-basal species has negitive growth rates
  r[species.basal] = - r[species.basal]  # the basal species has positive growth rates, 
  edges = sum(graph > 0)
  A = graph
  A[A > 0] = runif(edges, min = -1, max = 0)
  B = - e * t(A)
  C = A + B
  diag(C)[species.basal] = - 0.01
  list(r = r, C = -C)
}

init.lv1.foodweb <- function(graph) {
  s = dim(graph)[1]
  init = runif(s, min = 0, max = 1)
}

#' @title random interactions 
parms.lv1 <- function(s, k, delta) {
  g = graph.connected(s = s, k = k, gtype = 'er')
  A = as.matrix(get.adjacency(g))
  A[A > 0] = rnorm(n = s * k * 2, mean = 0, sd = delta)
  diag(A) = 1
  N = rep(1, s)
  r = A %*% N
  list(r = r, C = A)
}

init.lv1 <- function(s) {
  init = rep(1, s)
}

#' @references <Global stability in many-species systems. Goh.>
parms.lv1.goh.2species.predatorprey <- function() {
  r = c(-11, 5.6)
  C = matrix(c(1, -0.6, 1, -0.5), ncol = 2)
  list(r = r, C = - C)
}
init.lv1.goh.2species.predatorprey <- function() {
  init = c(3, 11)
}
parms.lv1.goh <- function() {
  r = c(2, 2.1, 1.5)
  C = matrix(c(0.8, 0.2, 1, 0.7, 0.9, 0.3, 0.5, 1, 0.2), ncol = 3)
  list(r = r, C = C)
}
init.lv1.goh <- function() {
  init = c(1, 2.126005, 1)  # sp1 : 0.668472 ~ 0.668473 ,  sp2 : 2.126005 ~ 2.126006 ,  sp3: 1.66054 ~ 1.66055
}


#' @title Jacobian matrix of LV1 model
phi.lv1 <- function(r, C, N) {
  J = diag(r - C %*% N)  - diag(N) %*% C
  J
}

analysis.lv1 <- function(Theta, alpha) {
  nstar = c(solve(Theta) %*% alpha)
  extinct = sum(nstar <= 0)
  survived = sum(nstar > 0)
}


#' @title Lotka-Volterra (LV) model with functional response of Holling Type I
model.lv1 <- function(time, init, parms) {
  r = parms[[1]]  # intrinsic growth rate
  C = parms[[2]]  # the competition matrix
  N = init  # initial state
  dN <- N * ( r - C %*% N )
  list(c(dN))
}

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

#' @title press perturbation experiment: continually add members of species by rate of <p>
model.lv2.press <- function(time, init, parms, ...) {
  r = parms[[1]]  # intrinsic growth rate
  C = parms[[2]]  # the competition matrix
  M = parms[[3]]  # the cooperation matrix
  h = parms[[4]]  # handling time
  N = init  # initial state
  dN <- N * ( r - C %*% N + (M %*% N) / (1 + h * M %*% N) )
  p = 0.1
  dN[1] = dN[1] - p
  list(c(dN))
}

#' @title the [parms] of mutualistic LV2 model
#' @param graph, the incident matrix of mutualistic networks which are bipartite
#' @param alpha.mu, alpha.sd, the intrinsic growth rate
#' @param beta0.mu, beta0.sd, the intraspecies competition
#' @param beta1.mu, beta1.sd, the interspecies competition
#' @param gamma.mu, gamma.sd, the interspecies cooperation
#' @param h.mu, h.sd, the Handling time, saturate coefficient
#' @return [parms] of [simObj] class
parms.lv2 <- function(graph, alpha.mu = 0.2, alpha.sd = 0.15, beta0.mu = 1, beta0.sd = 0.2, beta1.mu = 0.0, beta1.sd = 0.0,
                      gamma.mu = 1., gamma.sd = 0.2, h.mu = 0.1, h.sd = 0.05, delta = 0.5) {
  numP = dim(graph)[1]
  numA = dim(graph)[2]
  s = numP + numA
  r = runif(s) * 2 * alpha.sd + (alpha.mu - alpha.sd)
  
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
  list(r = r, C = C, M = M, h = h)  # the [parms] of ode model
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

parms.lv2.mutualism.3species <- function() {
  r = c(0.5, -0.1, -0.1)   # r[1] = 0.28365 ~ 0.28366
  C = matrix(c(0.5, 0, 0, 0, 0.5, 0.1, 0, 0.1, 0.5), ncol = 3)
  M = matrix(c(0, 0.5, 0.5, 0.5, 0, 0, 0.5, 0, 0), ncol = 3)
  h = c(0., 0.1, 0.3)
  list(r = r, C = C, M = M, h = h)  # the [parms] of ode model
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

#' @title the second extinction after ONE species is removed
sim.ode.extinct.one <- function(before.extinct) {
  nstar = before.extinct$nstar
  parms = before.extinct$parms
  model = before.extinct$model
  s = length(nstar)
  ode.outs =
  llply(1:s, function(i) {
    init = nstar
    init = init[- i]  # primary extinction
    parms$r = parms$r[- i]
    parms$C = parms$C[- i, - i]
    parms$M = parms$M[- i, - i]
    parms$h = parms$h[- i]
    ret = sim.ode.one(model, parms, init)
  })
}

#' @title the second extinction after TWO species is removed
sim.ode.extinct.two <- function(before.extinct) {
  nstar = before.extinct$nstar
  parms = before.extinct$parms
  model = before.extinct$model
  s = length(nstar)
  extinct.species = combn(s, 2)
  extinct.species.num = dim(extinct.species)[2]
  ode.outs =
    llply(1:extinct.species.num, .parallel = TRUE, function(i) {
      init = nstar
      init[ extinct.species[1, i] ] = 0  # primary extinction
      init[ extinct.species[2, i] ] = 0  # primary extinction
      ret = sim.ode.one(model, parms, init)
    })
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
    parms$r = parms$r - 0.02 #runif(length(nstar), min =  0.01, max = 0.02)
  }
  else if(perturb.type == 'lv2.growth.rate.dec.onepart') {
    numP = 100; numA = 100;
    #parms$r[1:numP] = parms$r[1:numP] - 0.02
    parms$r[(numP + 1):(numP+numA)] = parms$r[(numP + 1):(numP+numA)] - 0.04
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

sim.ode.old <- function(model, parms, init, steps = 10000, stepwise = 0.01, isout = FALSE, 
                        iter.steps = 10, iter.stepwise = 0.01, perturb, perturb.type = 'lv2.growth.rate.dec') {
  odemodel <- odeModel(
    main = model,
    parms = parms,
    init = init,
    times = c(from = 0, to = steps * stepwise, by = stepwise),
    solver = 'lsoda')
  
  ode.outs = list()
  for(i in 1:iter.steps) {
    odemodel <- sim(odemodel)
    ode.out = out(odemodel)
    
    nstar = as.numeric(ode.out[nrow(ode.out), 2:ncol(ode.out)]) 
    nstar[nstar < extinct.threshold] = 0  # species with biomass less than the threshold is considered to be extinct
    
    Phi = jacobian.full(y = nstar, func = model, parms = parms(odemodel))
    if (isout) {
      ret = list(out = ode.out, nstar = nstar, Phi = Phi, params = parms(odemodel))
    }
    else {
      ret = list(nstar = nstar, Phi = Phi, params = parms(odemodel))
    }
    ode.outs[[length(ode.outs) + 1]] = ret
    
    perturb.res = perturb(parms(odemodel), nstar, perturb.type)
    parms(odemodel) = perturb.res$parms
    nstar = perturb.res$nstar
    init(odemodel) = nstar
  }
  ode.outs
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



#' @references <A primer of ecology with R. Charpter 8: Multi basins of attraction>
model.scheffer <- function (time, init, parms, ...) 
{
  F <- init[1]
  S <- init[2]
  with(as.list(parms), {
    n <- N / (1 + qs * S + qf * F)
    dF <- rf * F * ( n / (n + hf) ) * ( 1 / ( 1 + af * F ) ) - lf * F
    dS <- rs * S * ( n / (n + hs) ) * ( 1 / ( 1 + as * S + b * F + W ) ) - ls * S
    return(list(c(dF, dS)))
  })
}

parms.scheffer <- function(N = 1, as = 0.01, af = 0.01, b = 0.02, qs = 0.075,
                           qf = 0.005, hs = 0, hf = 0.2, ls = 0.05, lf = 0.05, rs = 0.5,
                           rf = 0.5, W = 0) {
  c(N = N, as = as, af = af, b = b, qs = qs,
       qf = qf, hs = hs, hf = hf, ls = ls, lf = lf, rs = rs,
       rf = rf, W = W)
}

init.scheffer <- function(F = 10, S = 10) {
  c(F = F, S = S)
}


#' @references <On the structural stability of mutualistic systems. Bascompte.>
# parms.lv1 <- function(r, C) {
#   list(r = r, C = C)
# }
# A = parms.lv1(r = c(1, 1), C = matrix(c(1, 0.5, 0.5, 1), ncol = 2))  # feasible and global stable
# B = parms.lv1(r = c(1, 2), C = matrix(c(1, 0.5, 0.5, 1), ncol = 2))  # not feasible but global stable
# C = parms.lv1(r = c(1, 1), C = matrix(c(0.5, 1, 1, 0.5), ncol = 2))  # feasible but not global stable (only local stable)

#' @title get the critical mutualistic strength
#' @param tol, tolerance
get.gamma0.max <- function(graph, beta0 = 1, beta1 = 0.2, delta = 0, tol = 0) {
  numP = dim(graph)[1]
  numA = dim(graph)[2]
  s = numP + numA
  
  C = matrix(0, nrow = s, ncol = s)
  C[1:numP, 1:numP] = beta1
  C[(numP+1):s, (numP+1):s] = beta1
  diag(C) = beta0
  
  edges = sum(graph > 0)  # the number of edges
  M = as.one.mode(graph)  # transform to adjacency matrix (function of package [bipartite])
  degrees = rowSums(M)

  gamma0 = 1 / sqrt(edges)
  initCorrect = FALSE 
  repeat {
    M[M > 0] = gamma0  # endue values of interspecies cooperation
    M = M / degrees^delta  # trade-off of mutualistic strength    
    Theta = C - M 
    if (min(Re(eigen(Theta)$values)) > tol) {
      #print('increase one stepwise.')
      initCorrect = TRUE
      gamma0 = gamma0 + 0.0001
    }
    else {
      gamma0 = gamma0 - 0.0001
      M[M > 0] = gamma0  # endue values of interspecies cooperation
      M = M / degrees^delta  # trade-off of mutualistic strength    
      break
    }
  }
  if (initCorrect == FALSE) warning('Init value of gamma0 is not correct!')
  gamma0
}

#' @title get the interaction matrix, including the competition and cooperation sub-matrix
get.interaction.matrix <- function(graph, beta0.mu = 1, beta0.sd = 0., beta1.mu = 0.2, beta1.sd = 0, gamma.mu = 0, gamma.sd = 0, delta = 0) {
  numP = dim(graph)[1]
  numA = dim(graph)[2]
  s = numP + numA
  
  C = matrix(0, nrow = s, ncol = s)
  C[1:numP, 1:numP] = runif(numP * numP) * 2 * beta1.sd + (beta1.mu - beta1.sd)
  C[(numP+1):s, (numP+1):s] = runif(numA * numA) * 2 * beta1.sd + (beta1.mu - beta1.sd)
  #diag(C) = beta0
  diag(C) = runif(s) * 2 * beta0.sd + (beta0.mu - beta0.sd)
  
  edges = sum(graph > 0)  # the number of edges
  M = as.one.mode(graph)  # transform to adjacency matrix (function of package [bipartite])
  degrees = rowSums(M)
  #M[M > 0] = gamma0  # endue values of interspecies cooperation
  M[M > 0] = runif(2 * edges) * 2 * gamma.sd + (gamma.mu - gamma.sd)  # endue values of interspecies cooperation
  M = M / degrees^delta  # trade-off of mutualistic strength    
  
  list(C = C, M = M)
}


#' @title get the effective competition matrix, then get the structural vectors
get.structural.vectors <- function(C, M, numP, numA) {
  T = M %*% solve(C)
  diag(T) = 1
  Theta = T %*% (C - M)
  #Theta[Theta < 1e-10] = 0
  BetaP = Theta[1:numP, 1:numP]
  BetaA = Theta[(numP+1):(numP+numA), (numP+1):(numP+numA)]

  U = svd(BetaP)$u
  V = svd(BetaP)$v
  AlphaPL = U[,1]
  AlphaPR = V[,1]
  U = svd(BetaA)$u
  V = svd(BetaA)$v
  AlphaAL = U[,1]
  AlphaAR = V[,1]
  
  AlphaL = solve(T) %*% c(AlphaPL, AlphaAL)
  AlphaPL = AlphaL[1:numP]
  AlphaAL = AlphaL[(numP+1):(numP+numA)]
  AlphaR = solve(T) %*% c(AlphaPR, AlphaAR)
  AlphaPR = AlphaR[1:numP]
  AlphaAR = AlphaR[(numP+1):(numP+numA)]
  list(AlphaPL = -AlphaPL, AlphaAL = -AlphaAL, AlphaPR = -AlphaPR, AlphaAR = -AlphaAR, 
       AlphaL = -c(AlphaPL, AlphaAL), AlphaR = -c(AlphaPR, AlphaAR))
}

#' @references <On the structural stability of mutualistic systems. Bascompte.>
#' @param graph, the incident matrix of empirical or model-generated mutualistic networks
#' @param alpha, the intrinsic growth rates, calculated from the derivation from the structural vector
#' @param beta0, beta1, the intra- and inter- species competition strength
#' @param gamma0, the cooperation (mutualism) strength, which is up to the critical mutualistic strength.
#'        below the critical mutualistic strength assure the global stability
#' @param delta, the trade-off of mutualistic strength between plants and animals
#' @param h, the handling time
parms.lv2.softmean <- function(alpha = NULL, C, M, h.mu = 0, h.sd = 0) {
  r = alpha
  s = length(alpha)
  h = runif(s) * 2 * h.sd + (h.mu - h.sd)
  list(r = r, C = C, M = M, h = h)  # the [parms] of ode model
}

get.alpha.perturbation <- function(AlphaL, AlphaR) {
  s = length(AlphaL)
  alpha = runif(s, min = -0.001, max = 0.001) + AlphaL
  #r = rnorm(s, mean = 0, sd = 10) + alpha
  alpha
}

#' @title get angle between perturbed vector and structural vector
get.angle <- function(perturbed, AlphaL, AlphaR) {
  angleL = lengths.angle(x = perturbed, y = AlphaL)$angle
  angleR = lengths.angle(x = perturbed, y = AlphaR)$angle
  angle = (1 - cos(angleL) * cos(angleR)) /  (cos(angleL) * cos(angleR))
  c(angle = angle, angleL = angleL, angleR = angleR)
}

#' @title get the random perturbations to structural vectors
get.perturbed <- function(AlphaL, AlphaR, n) {
  s = length(AlphaL)  # number of species
  vars = runif(s, 0, 0.2)
  perturbed = laply(1:s, function(i) {
    rlnorm(n, meanlog = 0, sdlog = vars[i])
  })
  perturbed = AlphaL * perturbed
  perturbed
}


#' @references <Instability of a hybrid module of antagonistic and mutualistic interactions> Kondoh
model.kondoh <- function (time, init, parms, ...) 
{
  X <- init[1]
  Y <- init[2]
  Z <- init[3]
  with(as.list(parms), {
    dX <- (rx - ex * X - a * Y + u * Z / (hz + Z)) * X
    dY <- (g * a * X + d) * Y
    dZ <- (rz - ez * Z + v * X / (hx + X)) * Z
    return(list(c(dX, dY, dZ)))
  })
}

parms.kondoh <- function(rx = 0.2, rz = 0.2, ex = 1, ez = 1, u = 1, v = 1,
                         hz = 0.5, hx = 0.5, a = 0.5, g = 0.25, d = - 0.05) {
  c(rx = rx, rz = rz, ex = ex, ez = ez, u = u, v = v, hz = hz, hx = hx, a = a, g = g, d = d)
}

init.kondoh <- function(X = 1, Y = 1, Z = 1) {
  c(X = X, Y = Y, Z = Z)
}

