## Discrete Lotka-Volterra Model simulation

extinct.threshold.default = .Machine$double.eps * 100  # threshold of species extinction is 100 times machine precision
extinct.threshold = extinct.threshold.default
steady.threshold = 1e-8  # difference between species abundances of two successive steps is less than [steady.threshold]

#' @title model of Discrete Lotka-Volterra Equations of Holling Type I
model.dlv1 <- function(Nt, r, B) {
  Ntplus1 = Nt * exp(r - B %*% Nt)
  Ntplus1[Ntplus1 < extinct.threshold] = 0
  Ntplus1
}

#' @title parameters for DLV1 model based on some graph
parms.dlv1.rand <- function(graph, beta0 = 1, beta1.min = 0., beta1.max = 0.9) {
  s = dim(graph)[1]
  edges = sum(graph > 0)
  B = graph
  B[B > 0] = runif(edges, min = beta1.min, max = beta1.max)
  diag(B) = rep(beta0, s)
  B
}

#' @title simulate DLV1 model for times
sim.dlv1.times <- function(N0, r, B, times = 100) {
  s = dim(B)[1]
  output = matrix(0, nrow = times + 1, ncol = s)
  output[1, ] = N0
  for(t in 1:times) {
    output[t + 1, ] = model.dlv1(output[t, ], r, B)
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

w = rep(1, s)
A = phi.dlv1(r, B)
C = diag(0.1, s)
times = 10000
out = mar.sim(w, A, C, times)
matplot(out, type = 'l')
apply(out, 2, mean)
solve(I -A) %*% w
var(out)
matrix(solve(diag(1,s^2) - kronecker(A, A)) %*% as.vector(C), nrow = s, ncol = s)
