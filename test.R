library(corrplot)
col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white", 
                           "cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))	
col3 <- colorRampPalette(c("red", "white", "blue"))	
col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F", 
                           "cyan", "#007FFF", "blue","#00007F"))	

A = sim.lv2.graph(graph1, isout = TRUE)

Phi = A[[110]]$Phi
invPhi = - solve(A[[110]]$Phi)

corrplot(invPhi, method = 'circle', is.corr = FALSE )
corrplot(Phi, method = 'circle', is.corr = FALSE )


matplot(A$out[,2:21], type = 'l')

A = sim.lv2.alpha.dec(graph = graph1, isout = TRUE, dec.steps = 150)

lv2.Nstars = laply(A, function(one) {
  one$Nstar
})
matplot(lv2.Nstars, type = 'l')
