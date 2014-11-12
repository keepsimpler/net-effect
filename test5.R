library('ggplot2')

## generate random networks with different degree heterogeneitise
n = 20
k = 3
graphs = get.graphs.hetero(n = n, k = k, rep1 = 10, rep2 = 1)

## get degree heterogeneities of graphs
heteros = ldply(graphs, function(graph) get.degree.hetero(graph))
hist(heteros$V1, breaks = 50)

## get mean and variance of species abundances at equilibrium
beta1 =  0.1
res = get.results.by.graphs(graphs, beta1)

hetero.beta1.vars = get.nstars.vars.sum(res)
pairs(hetero.beta1.vars)
qplot(data = hetero.beta1.vars, x = beta1, y = hetero, color = vars.sum)

nstars = get.nstars(res)
selfvars = get.selfvars(res)
nstars.long = melt(nstars, id.vars = c('beta1', 'hetero'), variable.name = 'nstars.name', value.name = 'nstars')
selfvars.long = melt(selfvars, id.vars = c('beta1', 'hetero'), variable.name = 'selfvars.name', value.name = 'selfvars')
nstars.long$selfvars = selfvars.long$selfvars
nstars.long$selfvars.name = selfvars.long$selfvars.name
nstars.long = nstars.long %.% filter(hetero > 2.6)
qplot(data = nstars.long, x = nstars, y = selfvars, log = 'xy', color = hetero) + theme_bw()
qplot(data = nstars.long, x = nstars, y = selfvars, log = 'xy') + facet_wrap(~hetero, ncol = 5) + theme_bw()

graphs.full = get.graphs.full(n, 500)
res2 = get.results.by.variance(graphs.full, flag = FALSE)

vars.sum.by.variance = ldply(res2, function(res2.1) {
  ldply(res2.1, function(res2.2) {
    r.mean = res2.2$r.mean
    r.sd = res2.2$r.sd
    res = res2.2$res
    #feasible.stable = get.feasible.stable(res)
    ldply(res, function(one) {
      if (one$feasible && one$stable)
        c(r.mean = r.mean, r.sd = r.sd, hetero = one$hetero, beta1 = one$beta1, vars.self = sum(diag(one$vars)), vars.sum = sum(one$vars), nstar.sum = sum(one$nstar),
          eigs = sum(1 / (2*(eigen(one$phi)$values))), vars.max = sum(sqrt(diag(one$vars)))^2)
    })
  })
})

tmp = vars.sum.by.variance %.% filter(r.mean == 4) %.% select(vars.sum, vars.max, r.sd, hetero)

p <- ggplot(data = tmp, aes(x = vars.sum))
p + geom_histogram(binwidth = 0.1) +
  facet_wrap(~r.sd, scales = 'free') + 
  theme_bw()

ggplot(data = tmp, aes(x = hetero, y = vars.sum)) +
  geom_point() +
  facet_wrap(~r.sd, scales = 'free') + 
  theme_bw()

