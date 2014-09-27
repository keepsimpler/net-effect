s = 10
graph = graph.connected(s, gtype = 'complete')
graph = get.adjacency(graph)
graph = as.matrix(graph)

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
vars = matrix(solve(diag(1,s^2) - kronecker(A, A)) %*% as.vector(C), nrow = s, ncol = s)
