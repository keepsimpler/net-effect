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

A = sim.ode(model = model.lv2, parms = parms, init = init, isout = TRUE, iter.steps = 200, steps = 100,
            perturb = perturb, perturb.type = 'lv2.growth.rate.dec')

A = sim.ode(model = model.lv2, parms = parms, init = init, isout = TRUE, iter.steps = 2, steps = 10000, stepwise = 0.1,
            perturb = perturb, perturb.type = 'lv2.primary.extinction')

parms = parms.lv2.mutualism.3species()
init = c(1,1,1)
A = sim.ode(model = model.lv2, parms = parms, init = init, isout = TRUE, iter.steps = 2, steps = 10000, stepwise = 0.1,
            perturb = perturb, perturb.type = 'lv2.primary.extinction')
matplot(A[[100]]$out[,2:29], type = 'l')

B = laply(1:80, function(i) {
  sensitivity = get.sensitivity(A[[i]]$nstar, A[[i]]$Phi)
  sensitivity$sensitivity
  #sensitivity$lev
})
plot(A[[1]]$nstar - A[[2]]$nstar,  rowSums(-solve(A[[1]]$Phi) %*% diag(A[[1]]$nstar) %*% (A[[2]]$params$r - A[[1]]$params$r)))
plot(A[[2]]$nstar / A[[1]]$nstar,  rowSums(diag(1/A[[1]]$nstar) %*% -solve(A[[1]]$Phi) %*% diag(A[[1]]$nstar) %*% (A[[2]]$params$r - A[[1]]$params$r)))
plot((A[[2]]$nstar - A[[1]]$nstar) / A[[1]]$nstar[1], -solve(A[[1]]$Phi)[,1] / -solve(A[[1]]$Phi)[1,1])

ode.nstars = laply(A, function(one) {
  one$nstar
})
matplot(ode.nstars, type = 'l', lwd = 1.)
text(0, ode.nstars[1,],1:length(ode.nstars[1,]))
#legend("right", c("F", "S"), lty = 1:2, bty = "n")

ode.invPhi = laply(A, function(one) {
  Phi = one$Phi
  invPhi = -solve(Phi)
  c(det(Phi[-1,-1]), -det(Phi[-2,-1]), det(Phi[-3,-1]))  #invPhi[,1],
})

ode.neteffects = laply(A, function(one) {
   Phi = one$Phi
   lev = max(Re(eigen(Phi)$values))
   #sev = min(Re(eigen(Phi)$values))
   #invPhi = -solve(Phi)
  #c(neteffects = sum(invPhi) - sum(diag(invPhi)), sum = sum(invPhi), diagsum = sum(diag(invPhi)) )
  #sum(c(colSums(invPhi) / diag(invPhi)))
  #c(rowSums(invPhi))
  c(lev)
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




## plot straight lines  <abline> 
## Setup up coordinate system (with x == y aspect ratio):
plot(c(-2,3), c(-1,5), type = "n", xlab = "x", ylab = "y", asp = 1)
## the x- and y-axis, and an integer grid
abline(h = 0, v = 0, col = "gray60")
text(1,0, "abline( h = 0 )", col = "gray60", adj = c(0, -.1))
abline(h = -1:5, v = -2:3, col = "lightgray", lty = 3)
abline(a = 1, b = 2, col = 2)
text(1,3, "abline( 1, 2 )", col = 2, adj = c(-.1, -.1))

r1 <-3
r2 <-2
K1 <- 1.5
K2 <-2
alf12 <- 1
alf21 <- 2
Lotka <-function(t, N, pars)
{
  dN1   <- r1*N[1]*(1-(N[1]+alf12* N[2])/K1)
  dN2   <- r2*N[2]*(1-(N[2]+alf21* N[1])/K2)
  list(c(dN1 , dN2 ))
}

Ax   <- c(0,K2/alf21)
Ay   <- K2 - alf21* Ax
By   <- c(0,K1/alf12)
Bx   <- K1 - alf12* By
xlim <- range(c(Ax, Bx))
ylim <- range(c(Ay, By))
plot  (x=Ax,y=Ay, type="l", lwd=2,
       main="Competition phase-plane",
       # 1st isocline
       xlab="N1",ylab="N2",xlim=xlim,ylim=ylim)
lines (Bx,By,lwd=2,lty=2)            # 2nd isocline

library(deSolve)
trajectory <- function(N1,N2)
{
  times <-seq(0,30,0.1)
  state <-c(N1 = N1, N2 = N2)
  out   <-as.data.frame(ode(state,times, Lotka,NULL))
  lines (out$N1,out$N2,type="l")
  arrows(out$N1[10],out$N2[10],out$N1[11],out$N2[11],
         length=0.1,lwd=2)
}
trajectory (N1=0.05,N2=0.3)
trajectory (0.11,0.3)
trajectory (1.5,1.8)
trajectory (1.0,2.0)
trajectory (0.0,0.1)

# 4 equilibrium points
X     <- c(0, 0 , K1, (K1-alf12*K2)/(1-alf12*alf21))
Y     <- c(0, K2, 0 , (K2-alf21*K1)/(1-alf12*alf21))

require(rootSolve)
ei <- matrix(nrow=4,ncol=2)
for (i in 1:4)
{
  N1 <- X[i]
  N2 <- Y[i]
  # the Jacobian
  Jacob <- jacobian.full(y= c(N1,N2),func=Lotka)
  print(Jacob)
  # eigenvalues
  ei[i,] <- eigen(Jacob)$values
  # colors of symbols
  if (sign(ei[i,1])>0 & sign(ei[i,2])>=0) col <- "white"
  if (sign(ei[i,1])<0 & sign(ei[i,2])<=0) col <- "black"
  if (sign(ei[i,1])* sign(ei[i,2])    <0 ) col <- "grey"
  # equilibrium point plotting
  points(N1,N2,pch=22,cex=2.0,bg=col,col="black")
}
cbind(N1=X,N2=Y,ei)
