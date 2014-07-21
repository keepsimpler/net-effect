library(corrplot)
col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white", 
                           "cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))	
col3 <- colorRampPalette(c("red", "white", "blue"))	
col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F", 
                           "cyan", "#007FFF", "blue","#00007F"))	

A = sim.lv2.graph(graph1, isout = TRUE)

Phi = A[[107]]$Phi
invPhi = - solve(Phi)
diagInvPhi = diag(invPhi)
corInvPhi = diag(1/ sqrt(diagInvPhi)) %*% invPhi %*% diag(1 /sqrt(diagInvPhi))

corrplot(Phi, method = 'circle', is.corr = FALSE, col = col4(500) )
corrplot(corInvPhi, method = 'circle', is.corr = FALSE, col = col2(500) )


matplot(A[[2]]$out[,2:21], type = 'l')

A = sim.lv2.alpha.dec(graph = graph1, isout = TRUE, dec.steps = 150, dec.stepwise = 0.01)

A = sim.lv2.gamma.dec(graph = graph1, isout = TRUE, dec.steps = 100)

lv2.Nstars = laply(A, function(one) {
  one$Nstar
})
matplot(lv2.Nstars, type = 'l', lwd = 1.5)
text(0, lv2.Nstars[1,],1:20)

lv2.invPhi = laply(A, function(one) {
  Theta = - diag(one$Nstar) %*% one$Phi
  lev = max(Re(eigen(Theta)$values))
  sev = min(Re(eigen(Theta)$values))
  invPhi = - solve(one$Phi)
  invPhi.sum = sum(invPhi)
  invPhi.diagsum = sum(diag(invPhi))
  invPhi.neteffect = invPhi.sum - invPhi.diagsum
  invPhi.syn = invPhi.neteffect / invPhi.sum
  diagInvPhi = diag(invPhi)
  corInvPhi = diag(1/ sqrt(diagInvPhi)) %*% invPhi %*% diag(1 /sqrt(diagInvPhi))
  c(sum = invPhi.sum, diagsum = invPhi.diagsum, neteffect = invPhi.neteffect, lev = lev, sev = sev, 
    syn = sum(corInvPhi) - sum(diag(corInvPhi)))
})
matplot(lv2.invPhi, type = 'l', col = 1:6)
legend("right", inset=.05, legend=c('sum', 'diagsum', 'neteffect', 'lev', 'sev', 'syn'), pch=1, col=1:6, horiz=FALSE)

lv2.Nstars.neteffect = cbind(lv2.Nstars, lv2.invPhi[,'neteffect'] / 100)
lv2.Nstars.corInvPhi = cbind(lv2.Nstars, lv2.invPhi[,'syn'] / 30)
matplot(lv2.Nstars.corInvPhi, type = 'l')
text(0, lv2.Nstars[1,],1:20)

