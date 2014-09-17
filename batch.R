#!/usr/bin/Rscript

#nohup Rscript batch.r &

source('utils.r')
library(doMC)  # 
registerDoMC(20)  # register Multi Cores
getDoParWorkers()  # get available Cores

#load(file = 'graphs.assort.RData')
#heterogeneity.and.robust.assort = get.heterogeneity.and.robust(list(graphs.assort[[5]]), rep = 10)
#save(heterogeneity.and.robust.assort, file = 'heterogeneity.and.robust.assort.RData')

load(file = 'graphs.RData')
#graphs.out = get.graphs.out(graphs, rep = 10)
heterogeneity.and.robust.nobeta1.delta0 = get.heterogeneity.and.robust(graphs)
save(heterogeneity.and.robust.nobeta1.delta0, file = 'heterogeneity.and.robust.nobeta1.delta0.RData')