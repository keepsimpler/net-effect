## generate random networks with different degree heterogeneitise
n = 20
k = 3
graphs = get.graphs.hetero(n = n, k = k, rep1 = 10, rep2 = 1)

## get degree heterogeneities of graphs
heteros = ldply(graphs, function(graph) get.degree.hetero(graph))
hist(heteros$V1, breaks = 50)

## get mean and variance of species abundances at equilibrium
beta1 = - 0.1
res = get.results.by.graphs(graphs, beta1)

hetero.beta1.vars = get.hetero.beta1.vars(res)
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
