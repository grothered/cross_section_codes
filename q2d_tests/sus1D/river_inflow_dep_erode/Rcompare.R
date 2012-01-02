####Compare analytical and numerical solution to the case with a constant depth / velocity in space and time, a constant settling velocity, and a given suspended sediment inflow, and resuspension

# The solution should satisfy
# d/dx(Q*C) = -ws*Width*C + Qe*Width
# u*A*dC/dx = -ws*Width*C + Qe*Width
# dC/dx = -ws/(depth*u)*C + Qe/(depth*u)

# The solution is
# C = C0*exp(-ws/(depth*u)*X) +Qe/ws
# [this follows from the standard methods of solving such things with integrating factors]
# The constant C0 must be chosen so that the equation satisfies its own boundary conditions.

source('../../../gettide1.R')

u=0.697538250250218
depth=0.7168066
wset=0.001
rho=1026
f=0.07
# Critical shear
taucrit=0.2
# Shear
tau=rho*f/8*u^2

# Cohesive erosion formula
Qe=0.001*(tau-taucrit)/taucrit^(0.5)

# When X=0, C=1
# Solve for C0
C0 = (1-Qe/wset)/exp(0)

# Numerical parameters
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
theory=C0*exp(-wset/(depth*u)*x) +Qe/wset

points(seq(1,num_sects)*delX, theory ,t='l',col=2,cex=0.2)
legend('bottomleft',c('numerical','analytical'),lty=c(1,1),col=c(1,2))
dev.off()



