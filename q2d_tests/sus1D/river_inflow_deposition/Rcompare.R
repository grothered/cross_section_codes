####Compare analytical and numerical solution to the case with a constant depth / velocity in space and time, a constant settling velocity, and a given suspended sediment inflow

#The solution should satisfy
# d/dx(Q*C) = -ws*Width*C
# u*A*dC/dx = -ws*Width*C
# dC/dx = -ws/(depth*u)*C
# so C = C0*(exp(-ws/(depth*u)*x))

source('../../../gettide1.R')

u=0.697538250250218
depth=0.7168066
wset=0.001
num_sects=100
delX=50
x0=num_sects*delX

getf(num_sects)

pdf(file='compare.pdf',width=5.4,height=5.4)
#par(mfrow=c(2,2))
#par(mar=c(4,4,2,0.2))
#par(cex=0.6)

#set k to whatever time step you want.
k=dim(s1)[1]

plot(seq(1,num_sects)*delX, s1[k,],t='o',main=paste('time = ', floor(t1[k]),' seconds'),cex=0.2)

# analytical solution
x = x0-( seq(0,num_sects-1) +0.5) *delX
theory=1.0*exp(-wset/(depth*u)*x)

points(seq(1,num_sects)*delX, theory ,t='l',col=2,cex=0.2)
legend('topleft',c('numerical','analytical'),lty=c(1,1),col=c(1,2))
dev.off()



