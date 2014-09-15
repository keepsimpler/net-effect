#!/usr/bin/Rscript

#nohup Rscript batch.r &

source('utils.r')
library(doMC)  # 
registerDoMC(20)  # register Multi Cores
getDoParWorkers()  # get available Cores

load(file = 'graphs.assort.RData')
heterogeneity.and.robust.assort = get.heterogeneity.and.robust(list(graphs.assort[[5]]), rep = 10)
save(heterogeneity.and.robust.assort, file = 'heterogeneity.and.robust.assort.RData')

# load(file = 'graphs.RData')
# graphs.out = get.graphs.out(graphs, rep = 10)
# save(graphs.out, file = 'graphs.out.RData')