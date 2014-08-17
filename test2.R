#' @references <On the structural stability of mutualistic systems. Bascompte.>
M_SD_001 = as.matrix(M_SD_001)
graph = M_SD_001
graph[graph > 0] = 1
numP = dim(graph)[1]
numA = dim(graph)[2]
gamma0.max = get.gamma0.max(graph = graph, beta0 = 1, beta1 = 0.2, delta = 0.5, tol = 0.0)
Theta = get.interaction.matrix(graph = graph, beta0.mu = 1, beta0.sd = 0.2, beta1.mu = 0.2, beta1.sd = 0.1, 
                               gamma.mu = gamma0.max - 0.3, gamma.sd = 0.1, delta = 0.5)  # 0.0402
SV = get.structural.vectors(C = Theta$C, M = Theta$M, numP = numP, numA = numA)  # Structural Vectors
alpha = rowSums(Theta$C - Theta$M)
n = 1000
Alphas = get.perturbed(SV$AlphaL, SV$AlphaR, n)
angles = adply(Alphas, .margins = 2, function(onecol) {
  get.angle(onecol, SV$AlphaL, SV$AlphaR)
})
hist(angles$angle, breaks = 200)
angles = angles %.% arrange(angle)  # sort by angles
A = llply(1:500, function(i) {
  print(i)
  angle = angles[i,]
  index = angle$X1
  alpha = Alphas[, index]
  parms = parms.lv2.softmean(alpha = alpha, C = Theta$C, M = Theta$M, h.mu = 0.2, h.sd = 0.1)
  init = rep(1, numP + numA)  # c(solve(Theta$C - Theta$M) %*% alpha)  # rep(1, numP + numA)  #alpha
  out = sim.ode.one(model = model.lv2, parms = parms, init = init)
  out
})

c(solve(Theta$C - Theta$M) %*% alpha)



parms = parms.lv2(graph)
init = init.lv2(parms)

A = sim.ode.one(model = model.lv2.press, parms, init, stepwise = 0.01)

A = sim.ode.one(model.lv2, parms, init)
C = sim.ode.extinct.one(A)
plot(B[[1]]$nstar - A$nstar,  A$invPhi[,1])

A = sim.ode(model = model.lv2, parms = parms, init = init, isout = TRUE, iter.steps = 90,
            perturb = perturb, perturb.type = 'lv2.growth.rate.dec')

A = sim.ode(model = model.lv2, parms = parms, init = init, isout = TRUE, iter.steps = 2, steps = 10000, stepwise = 0.1,
            perturb = perturb, perturb.type = 'lv2.primary.extinction')

parms = parms.lv2.mutualism.3species()
init = c(1,1,1)
A = sim.ode(model = model.lv2, parms = parms, init = init, isout = TRUE, iter.steps = 2, steps = 10000, stepwise = 0.1,
            perturb = perturb, perturb.type = 'lv2.primary.extinction')
matplot(A[[100]]$out[,2:29], type = 'l')

B = llply(1:110, function(i) {
  sensitivity = get.sensitivity(A[[i]]$nstar, A[[i]]$Phi)
  rowSums(sensitivity$sensitivity)
  #sensitivity$lev
})
plot(A[[1]]$nstar - A[[2]]$nstar,  rowSums(-solve(A[[1]]$Phi) %*% diag(A[[1]]$nstar) %*% (A[[2]]$params$r - A[[1]]$params$r)))
plot(A[[2]]$nstar / A[[1]]$nstar,  rowSums(diag(1/A[[1]]$nstar) %*% -solve(A[[1]]$Phi) %*% diag(A[[1]]$nstar) %*% (A[[2]]$params$r - A[[1]]$params$r)))
plot((A[[2]]$nstar - A[[1]]$nstar) / A[[1]]$nstar[1], -solve(A[[1]]$Phi)[,1] / -solve(A[[1]]$Phi)[1,1])

ode.nstars = laply(A, function(one) {
  one$nstar
})
matplot(ode.nstars, type = 'l', lwd = 1.5)
text(0, ode.nstars[1,],1:length(ode.nstars[1,]))
#legend("right", c("F", "S"), lty = 1:2, bty = "n")

ode.invPhi = laply(A, function(one) {
  Phi = one$Phi
  invPhi = -solve(Phi)
  c(invPhi)
})

ode.neteffects = laply(A, function(one) {
   Phi = one$Phi
   lev = max(Re(eigen(Phi)$values))
   sev = min(Re(eigen(Phi)$values))
   invPhi = -solve(Phi)
  c(neteffects = sum(invPhi) - sum(diag(invPhi)), sum = sum(invPhi), diagsum = sum(diag(invPhi)) )
  sum(c(colSums(invPhi) / diag(invPhi)))
  #c(rowSums(invPhi))
  #c(lev - sev)
})
matplot(ode.neteffects, type = 'l')
legend("right", c("neteffect", "sum", 'diagsum'), lty = 1:4, bty = "n")

parms = parms.scheffer(N = 0.5)
init = init.scheffer()
A = sim.ode(model = model.scheffer, parms = parms, init = init, isout = TRUE, iter.steps = 40,steps = 1000, stepwise = 1,
            perturb = perturb, perturb.type = 'scheffer.nutrient.inc')


repeat{
  graph = niche.model(50, 0.2)
  G = graph.adjacency(graph)
  if (is.connected(G)) break  
}
plot.igraph(G)

parms = parms.lv1.foodweb(graph = graph, e = 0.1)
init = init.lv1.foodweb(graph)
A = sim.ode.one(model = model.lv1, parms = parms, init = init)
matplot(A$out[1:10000,2:51], type = 'l')
sum(A$nstar > 0)


C = laply(1:50, function(i) {
  B = jacobian.full(A$out[i, 2:21], func = A$model, parms = A$parms)
  invB = -solve(B)
  c(det = det(invB), sum = sum(invB))
})



## read the web-of-life files of Bascompte et al.
weboflife = numeric()
for (i in 1:59) {
  networkname = paste('M_PL_', sprintf('%03d', i), sep='')
  weboflife[i] = networkname
  assign(networkname,
         read.csv(paste("~/Data/Bascompte/web-of-life/M_PL_", sprintf("%03d", i), ".csv", sep=""), header=T, row.names = 1))
}
for (i in 1:30) {
  networkname = paste('M_SD_', sprintf('%03d', i), sep='')
  weboflife[59 + i] = networkname
  assign(networkname,
         read.csv(paste("~/Data/Bascompte/web-of-life/M_SD_", sprintf("%03d", i), ".csv", sep=""), header=T, row.names = 1))  
}


## plot straight lines  <abline> 
## Setup up coordinate system (with x == y aspect ratio):
plot(c(-2,3), c(-1,5), type = "n", xlab = "x", ylab = "y", asp = 1)
## the x- and y-axis, and an integer grid
abline(h = 0, v = 0, col = "gray60")
text(1,0, "abline( h = 0 )", col = "gray60", adj = c(0, -.1))
abline(h = -1:5, v = -2:3, col = "lightgray", lty = 3)
abline(a = 1, b = 2, col = 2)
text(1,3, "abline( 1, 2 )", col = 2, adj = c(-.1, -.1))