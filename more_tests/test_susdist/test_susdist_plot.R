# R code to plot outputs.
dta=read.table('test_susdist_outfile.out')

png(file='Analytical_compare.png',width=5.4,height=8,units='in',res=300)
par(mfrow=c(2,1))

plot(dta[,3]/3.426391,t='l',ylab='Concentration')
points(dta[,4]/0.01486444,t='l',col=2)
title('Analytical steady state (red) \n and numerical (black) concentration')

plot(dta[,3]/3.426391-dta[,4]/0.01486444,t='l', ylab='Concentration difference')
abline(h=0,col=2)
title('Scaled Numerical minus analytical')
dev.off()
