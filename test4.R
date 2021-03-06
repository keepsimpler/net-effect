s = 5
k = 1.5
r = rep(1, s)
C = diag(1, s)
graph = graph.connected(s, k, gtype = 'er')
graph = as.matrix(get.adjacency(graph))
B = parms.dlv1.rand(graph, beta1.min = 0.1, beta1.max = 0.1)
nstar = solve(B) %*% r
phi = - diag(c(nstar)) %*% B
Q = eigen(phi)$vectors
eigs = eigen(phi)$values
E = solve(Q) %*% C %*% t(solve(Q))
X = E / outer(eigs, eigs, FUN = '+')
Q %*% X %*% t(Q)


source('graphs.rewiring.unipart.r')

# generate list of list of graphs which have different degree heterogeneity 
# but with same node number and mean degree
graphs = llply(1:5, .parallel = TRUE, function(i) graphs.rewiring.unipart(s, k))
heterogeneity = ldply(graphs, function(graphs.1) {
  ldply(graphs.1, function(graph) {
    graph = graph$A
    get.degree.heterogeneity(graph)
    get.degree.hetero(graph)
  })
})

a = ldply(1:100, function(i) {
  graph = graph.connected(50, 2, gtype = 'regular', expower = 3.5)
  get.degree.hetero(graph)
})




beta1s = seq(from = -0.2, to = 0.2, by = 0.02)
res = llply(graphs, .parallel = TRUE, function(graphs.1) {
  llply(graphs.1, .parallel = TRUE, function(graph) {
    graph = graph$A
    llply(beta1s, function(beta1) {
      print(beta1)
      B = parms.dlv1.rand(graph, beta1.min = beta1, beta1.max = beta1)
      s = dim(B)[1]
      r = rep(1, s)
      out.lv1 = analysis.lv1(r, B)
      ret = list()
      if (max(Re(eigen(out.lv1$phi, only.values = TRUE)$values)) >= 0) stable = FALSE else stable = TRUE
      if (out.lv1$extinct == 0 && stable) {
        heterogeneity = get.degree.heterogeneity(graph)
        nstar = out.lv1$nstar
        phi = out.lv1$phi
        C = diag(1, s)
        C = diag(nstar^2) %*% C
        vars = mou.vars(phi, C)
        ret = list(heterogeneity = heterogeneity, beta1 = beta1, feasible = TRUE, stable = TRUE, nstar = nstar, phi = phi, vars = vars, B = B)
      }
      else if (out.lv1$extinct == 0 && !stable) {
        ret = list(heterogeneity = heterogeneity, beta1 = beta1, feasible = TRUE, stable = FALSE)
      }
      else {
        ret = list(heterogeneity = heterogeneity, beta1 = beta1, feasible = FALSE, stable = FALSE)
      }
      ret         
    })
  })
})

beta1 = 0.10
res = get.results.by.graphs(graphs, beta1)

hetero.beta1.vars = get.hetero.beta1.vars(res)
hetero.beta1.vars = hetero.beta1.vars %.% filter(beta1 != -0.20 & beta1 != -0.18)
qplot(data = hetero.beta1.vars, x = beta1, y = heterogeneity, color = vars.sum)

nstars = get.nstars(res)
selfvars = get.selfvars(res)
nstars.long = melt(nstars, id.vars = c('beta1', 'hetero'), variable.name = 'nstars.name', value.name = 'nstars')
selfvars.long = melt(selfvars, id.vars = c('beta1', 'hetero'), variable.name = 'selfvars.name', value.name = 'selfvars')
nstars.long$selfvars = selfvars.long$selfvars
nstars.long$selfvars.name = selfvars.long$selfvars.name
qplot(data = nstars.long, x = nstars, y = selfvars, color = hetero) + theme_bw()
qplot(data = nstars.long, x = nstars, y = selfvars, log = 'xy') + facet_wrap(~hetero, ncol = 5) + theme_bw()

plot(heterogeneity.and.vars$heterogeneity, heterogeneity.and.vars$vars.self / heterogeneity.and.vars$vars.sum)
pairs(heterogeneity.and.vars)

plot(unlist(nstars), unlist(selfvars))

a = ddply(nstars.long, .(hetero), function(one) {
  summary(lm(log(one$nstars) ~ log(one$selfvars), data=faithful, ))$coefficients[, 'Estimate']      #r.squared
})



ntried = 500
res.competition.sf = llply(1:ntried, .parallel = TRUE, function(i) {
  print(i)
  graph = graph.connected(s, k, gtype = 'sf', expower = 3.5)
  graph = get.adjacency(graph)
  graph = as.matrix(graph)
  #eigs = eigen(graph)$values
  B = parms.dlv1.rand(graph, beta1.min = 0.1, beta1.max = 0.1)
  out.lv1 = analysis.lv1(r, B)
  ret = list()
  if (max(Re(eigen(out.lv1$phi, only.values = TRUE)$values)) >= 0) stable = FALSE else stable = TRUE
  if (out.lv1$extinct == 0 && stable) {
    vars = mou.vars(out.lv1$phi, C)
    ret = list(nstar = out.lv1$nstar, phi = out.lv1$phi, vars = vars)
  }
  ret
})

nstar.competition.sf = ldply(res.competition.sf, function(one) {if (length(one) > 0) one$nstar})
vars.self.competition.sf = ldply(res.competition.sf, function(one) {if (length(one) > 0) diag(one$vars)})
a.competition.sf = ldply(res.competition.sf, function(one) {
  if (length(one) > 0)
    c(nstar.sum = sum(one$nstar), vars.sum = sum(one$vars), vars.max = sum(sqrt(diag(one$vars)))^2, vars.self = sum(diag(one$vars)),
    eigs = sum(1 / (2*Mod(eigen(one$phi)$values))))
  })
plot(unlist(nstar.competition.sf), unlist(vars.self.competition.sf))
plot(a.competition.sf$nstar.sum, a.competition.sf$vars.max)
plot(a.competition.sf$nstar.sum, a.competition.sf$vars.self)
plot(a.competition.sf)







res = feasible.and.stable(graph, interact.type = 'competition', interact.strength = 0.3, ntried = 1000)
ret = disply.eigs.distribution(res)
ret2 = ldply(res, function(one) {
  ret = NULL
  if (one$feasible == 1 && one$stable == 1) {
    eigs = eigen(one$vars, only.values = TRUE)$values
    entropy = sum(log(eigs)*eigs)
    vk = max(eigs)^2
    ret = c(mean = mean(eigs), sd = sd(eigs), syn = sum(one$vars) / sum(sqrt(diag(one$vars)))^2, sumvar = sum(one$vars),
            det = det(one$vars), entropy = entropy, vk = vk)
    #ret = c(one$Nstar)
    #ret = c(diag(one$vars))
  }
  ret
})



B = parms.dlv1.rand(graph)
w = rep(1, s)
A = phi.dlv1(r, B)
C = diag(0.1, s)
times = 500000
out = mar.sim(w, A, C, times)
matplot(out, type = 'l')
apply(out, 2, mean)
solve(I -A) %*% w
var(out)
vars = matrix(solve(diag(1,s^2) - kronecker(phi, phi)) %*% as.vector(C), nrow = s, ncol = s)

vars = matrix(solve(kronecker(diag(1, s), phi) + kronecker(phi, diag(1, s))) %*% as.vector(C), nrow = s, ncol = s)
