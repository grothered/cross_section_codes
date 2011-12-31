####Compare analytical and numerical solution to the case with a constant depth / velocity in space and time, constant diffusion coef, and Criver imposed at a given instant of time. Solution in Chanson (2004), where he studies the reduced diffusion model for river flow - same as this.

#The solution should be 
#C(x,t) = 1/2*(1+exp(u(x)/D))*(1-erf((x0 -u(t-t0))/sqrt(4*D*(t-t0)))  )
#where:
#D=20
#x0=190*100
#u=-0.5 (at least correct to within a few mm/s)
#t0=t[200]
#library(gsl)

source('../../../gettide1.R')

u=0.697538250250218
D=20.00000
num_sects=1000
delX=5
x0=num_sects*delX
t0=40*300

getf(num_sects)

erf<-function(x){
   # This code uses the analytical solutions from Chanson (2004:343), and I
   # use this formula based on his description of erf in the appendices (p
   # 384-385). 
   # Note however that the formulae he gives are inconsistent (For example, his
   # tabulated values do not agree with the small/large parameter
   # approximations that he suggests). 
   #
   # This one agrees with his tabulated values, and with pythons scipy.special.erf.
   # By my calculations it does not agree with Chanson's asymptotic results
   2*pnorm(x,mean=0,sd=sqrt(1/2))-1.0
}

#png(file='compare.png',width=5.4,height=5.4,res=1000,units='in')
pdf(file='compare.pdf',width=5.4,height=5.4)
par(mfrow=c(2,2))
par(mar=c(4,4,2,0.2))
par(cex=0.6)

for(k in (c(280))){#,270,280,290))){
#Set k to whatever time step you want.
    plot(seq(1,num_sects)*delX, s1[k,],t='o',main=paste('Time = ', floor(t1[k]),' seconds'),cex=0.2)
    # Analytical solution
    X = x0-seq(0,num_sects-1)*delX
    erfterm1=erf( (X -u*(t1[k]-t0))/(sqrt(4*D*(t1[k]-t0)) ))
    erfterm2=erf( (X +u*(t1[k]-t0))/(sqrt(4*D*(t1[k]-t0)) ))

    theory=0.5*( (1-erfterm1) + exp(u*(X)/(D))*(1-erfterm2))
    points(seq(1,num_sects)*delX, theory ,t='l',col=2,cex=0.2)
}
legend('bottomright',c('Numerical','Analytical'),lty=c(1,1),col=c(1,2))
dev.off()




#png(file='compare.png',width=5.4,height=5.4,res=1000,units='in')
#par(mfrow=c(2,2))
#par(mar=c(4,4,2,0.2))
#par(cex=0.6)
#for(k in (c(280))){
##Set k to whatever time step you want.
#
##plot(seq(1,220)*100, s1[k,],t='o',main=paste('Time = ', floor(t1[k]),' seconds'),cex=0.2)
#erfterm1=erf( (220*100-seq(1,220)*100 -0.5*(t1[k]-t0))/(sqrt(4*20*(t1[k]-t0)) ))
#erfterm2=erf( (220*100-seq(1,220)*100 +0.5*(t1[k]-t0))/(sqrt(4*20*(t1[k]-t0)) ))
#
#theory=0.5*( (1-erfterm1) + exp(-0.5*(220*100-seq(1,220)*100)/(20))*(1-erfterm2))
#plot(seq(1,220)*100, theory ,t='o',col=2,cex=0.2)
#}
#legend('bottomright',c('Numerical','Analytical'),lty=c(1,1),col=c(1,2))
#dev.off()
