## Discrete Lotka-Volterra Model simulation
library(plyr)

extinct.threshold.default = .Machine$double.eps * 100  # threshold of species extinction is 100 times machine precision
extinct.threshold = extinct.threshold.default
steady.threshold = 1e-8  # difference between species abundances of two successive steps is less than [steady.threshold]

#' @title model of Discrete Lotka-Volterra Equations of Holling Type I
model.dlv1 <- function(Nt, r, B, e) {
  Ntplus1 = Nt * exp(r - B %*% Nt + e)
  Ntplus1[Ntplus1 < extinct.threshold] = 0
  Ntplus1
}

#' @title simulate DLV1 model for times
sim.dlv1.times <- function(N0, r, B, C, times = 100) {
  s = dim(B)[1]
  es = mvrnorm(times, mu = rep(0, s), Sigma = C)
  output = matrix(0, nrow = times + 1, ncol = s)
  output[1, ] = N0
  for(t in 1:times) {
    output[t + 1, ] = model.dlv1(output[t, ], r, B, es[t, ])
    output
  }
  output
}

#' @title simuate DLV1 model to steady state
sim.dlv1 <- function(N0, r, B, steady.threshold) {
  s = dim(B)[1]
  output = matrix(0, nrow = 1, ncol = s)
  output[1, ] = N0
  t = 1
  repeat {
    Ntplus1 = model.dlv1(output[t, ], r, B)
    if (all(abs(Ntplus1 - output[t,]) < steady.threshold)) break
    output = rbind(output, c(Ntplus1))
    t = t + 1
  }
  output
}

#' @title get the community matrix for DLV1 model, only in feasbile (positive) equilibrium
phi.dlv1 <- function(r, B) {
  Nstar = solve(B) %*% r 
  if (any(Nstar <= 0)) {  # if a feasbile equilibrium doesn't exist
    warning('Not a feasible (positive) equilibrium!')
    return(NULL)
  }
  else {
    phi = diag(c(Nstar)) %*% B
    I = diag(1, nrow = dim(B)[1])
    phi = I - phi
    return(phi)
  }
}

#' @title parameters for DLV1 model based on the adjacency matrix of graph
parms.dlv1.rand <- function(graph, beta0 = 1, beta1.min = 0., beta1.max = 0.) {
  s = dim(graph)[1]
  edges = sum(graph > 0)
  B = graph
  B[B > 0] = runif(edges, min = beta1.min, max = beta1.max)
  diag(B) = rep(beta0, s)
  B
}
#B = parms.dlv1.rand(graph, beta0 = 1, beta1.min = 0., beta1.max = 0.)  # no interaction
#B = parms.dlv1.rand(graph, beta0 = 1, beta1.min = 0., beta1.max = 0.5)  # competition interaction
#B = parms.dlv1.rand(graph, beta0 = 1, beta1.min = -0.5, beta1.max = 0.)  # cooperation interaction
#B = parms.dlv1.rand(graph, beta0 = 1, beta1.min = -0.5, beta1.max = 0.5)  # arbitrary interaction

parms.dlv1.predator.prey <- function(graph, beta0 = 1, beta1.min = 0., beta1.max = 0.5) {  # predator-prey interaction
  s = dim(graph)[1]
  edges = sum(graph > 0)
  B = graph
  Bij = runif(edges / 2, min = beta1.min, max = beta1.max)
  Bji = runif(edges / 2, min = beta1.min, max = beta1.max)
  signs = sample(c(-1, 1), edges / 2, replace = TRUE)
  B[upper.tri(B) > 0] = Bij * signs  #
  B = t(B)
  B[upper.tri(B) > 0] = - Bji * signs
  diag(B) = rep(beta0, s)
  B
}
#B = parms.dlv1.predator.prey(graph, beta0 = 1, beta1.min = 0, beta1.max = 0.5)

parms.dlv1 <- function(graph, interact.type = 'no', interact.strength = 0.5) {
  if (interact.type == 'arbitrary')
    B = parms.dlv1.rand(graph, beta0 = 1, beta1.min = - interact.strength, beta1.max = interact.strength)
  else if (interact.type == 'competition')
    B = parms.dlv1.rand(graph, beta0 = 1, beta1.min = 0., beta1.max = interact.strength)
  else if (interact.type == 'cooperation')
    B = parms.dlv1.rand(graph, beta0 = 1, beta1.min = - interact.strength, beta1.max = 0.)
  else if (interact.type == 'predator.prey')
    B = parms.dlv1.predator.prey(graph, beta0 = 1, beta1.min = 0., beta1.max = interact.strength)
  else  # none interaction
    B = parms.dlv1.rand(graph, beta0 = 1, beta1.min = 0., beta1.max = 0.)
  B  
}

#' @title get probability of feasible equilibrium and local stable
feasible.and.stable <- function(graph, interact.type = 'no', interact.strength = 0.5, ntried = 1000) {
  ret = llply(1:ntried, function(i) {
    B = parms.dlv1(graph, interact.type, interact.strength)
    s = dim(graph)[1]
    r = rep(1, s)
    Nstar = solve(B) %*% r
    if (any(Nstar <= 0)) {  # not feasible
      ret = list(feasible = 0, stable = 0, B = B, Nstar = Nstar)
    }
    else {  # feasible
      phi = phi.dlv1(r, B)
      eigs = eigen(phi, only.values = TRUE)$values
      if ((any(Mod(eigs) > 1)))  # not stable
        ret = list(feasible = 1, stable = 0, B = B, Nstar = Nstar, phi = phi)
      else {  # stable
        C = diag(1, s)
        vars = matrix(solve(diag(1, s^2) - kronecker(phi, phi)) %*% as.vector(C), nrow = s, ncol = s)
        #vars.self = diag(vars)
        #vars.cor = sum(vars) - sum(vars.self)
        ret = list(feasible = 1, stable = 1, B = B, Nstar = Nstar, phi = phi, vars = vars)        
      }
    }
    ret
  })
  ret
}

disply.eigs.distribution <- function(res) {
  eigs = ldply(res, function(one) {
    res = NULL
    if (one$feasible == 1 && one$stable == 1) {
      res = eigen(one$vars, only.values = TRUE)$values
    }
    res
  })
  #plot(unlist(eigs))
  #hist(unlist(eigs), breaks = 100)
  eigs
}