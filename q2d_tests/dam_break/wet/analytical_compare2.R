###R code to compare the wet dam break, numerical and analytical cases.

#We have a rectangular channel with an initial depth of 10m below halfway,  and a depth of 1 above this. An analytical solution is presented e.g. in Zoppou and Roberts (2003). It's a bit complex


source('../../../gettide1.R')
l=1000
getf(l)

#t1=t1-t1[1]
delX=2 #Spatial increment
g=9.8
j=311 #The plot
h1=10
h0=1.0
xz=(seq(1,l)-l/2)*delX

u=rep(0,l)
h=rep(h1,l)

#Now, to get started on the analytical solution, (Zoppou and Roberts, 2003), we need to estimate the shock speed. This requires solving a nonlinear equation.
S2fun<-function(S2){
abs(-S2+2*sqrt(g*h1)+g*h0/(4*S2)*(1+sqrt(1+8*S2^2/(g*h0)) ) -(2*g*h0*sqrt(1+8*S2^2/(g*h0)) -2*g*h0 )^0.5)
}
a11=optimize(S2fun,interval=c(1.0E-6,1000))
S2=a11$minimum
##So now we have the shock speed. We need u2, and h2
u2=S2-g*h0/(4*S2)*(1+sqrt(1+8*S2^2/(g*h0)))

#h2 also needs minimizing
h2fun<-function(h2){
abs(-h2 + h0/2*sqrt(1+8*(2*h2/(h2-h0)*(sqrt(h1)-sqrt(h2))/sqrt(h0) )^2  ) -1/2)
}
a11=optimize(h2fun,interval=c(h0+1E-6,h1)) #This interval is important, should be sensible I guess
h2=a11$minimum

#Now find one of the intervals
b=which( (xz> -t1[j]*sqrt(g*h1))&(xz <= t1[j]*(u2-sqrt(g*h2)))  )

u[b]=2/3*(sqrt(g*h1)+xz[b]/t1[j])
h[b]=4/(9*g)*(sqrt(g*h1)-xz[b]/(2*t1[j]) )^2

##Now find the next interval
b=which( ( xz>t1[j]*(u2-sqrt(g*h2)))&(xz<t1[j]*S2  ))
u[b]=u2
h[b]=h2

#And the final interval
b=which(xz>=t1[j]*S2)
u[b]=0
h[b]=h0


v=which(abs(u)>0)
postscript('Dam_break_wet.eps',width=5.4,height=5.4,onefile=T,horizontal=F)
par(mfrow=c(2,1))
par(mar=c(4,4,1,0.1))
plot(xz,u,t='l',col=2,xlim=c(xz[min(v)-20],xz[ max(v)+20]),xlab='x',ylab='Velocity (m/s)',las=1)
points(xz,d1[j,]/A1[j,],t='o',cex=0.3)
points(xz,u,t='l',col=2)

plot(xz,h,t='l',col=2,xlim=c(xz[min(v)-20],xz[ max(v)+20]),xlab='x',ylab='Water elevation (m)',las=1)
points(xz,w1[j,],t='o',cex=0.3)
points(xz,h,t='l',col=2)
legend('topright',c('Analytical','Numerical'),col=c(2,1),lty=c(1,1),pch=c(NA,1),pt.cex=c(NA,0.3))
dev.off()
#Can compare errors with other studies, e.g. Zoppou and Roberts (2003), Sanders (2001), Vincent et al., (2000).
