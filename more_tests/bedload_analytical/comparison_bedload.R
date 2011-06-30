#source('cros.R')
#get(2000)

###Here we can compare the analytical and numerical model predictions for the bedload model.

rho=1026
rhos=2600
taue=0.4
D=0.0004*0.1/rhos

# a1*(tau-taue) = rate of resuspension
a1=0.000228/(rhos*sqrt(taue))

# k2*a1*(tau-taue)*dh/dy = rate of lateral bedload transport
k2a1= 0.5*0.000062*(0.000062*((rhos/rho-1)*9.8/1e-06^2)^(1/3))^(-0.3)*1/sqrt(rho)*1/taue*sqrt(taue)

##Doing the calculation to get rb, notice that in the code, I scale downslope transport by sqrt(rho*9.8*(rhos/rho-1)*0.000062/tau), whereas in the notation in the paper, it is sqrt(taue/tau). But, since taue is constant, effectively we have that sqrt(taue) = 0.63245, while sqrt(rho*9.8*(rhos/rho -1)*0.000062) = 0.97790, and so their ratio is 1.546

png('Analytical_compare.png', width=5.4,height=4,units='in', res=300)
k2=k2a1/a1

# Index of 'h' that we plot
iii=50 #:7*200
#iii=200

for(indy in iii){

    ##Note - I had been making a small error in the code, using 1.6 instead of (rhos/rho-1) sometimes. Corrected ,but we can still do comparison with the buggy code using the following hack

    #k2a1= 0.5*0.000062*(0.00003*((rhos/rho-1)*9.8/1e-06^2)^(1/3) )^(-0.3)*1/sqrt(rho)*1/sqrt(taue)*sqrt(rho*9.8*1.6*0.000062/taue)

    Sf=max(tau[indy,])/(-min(h[indy,])*rho*9.8)

    ##Predicted centre depth
    Dpred=(3*D/(2*a1) +taue)/(rho*9.8*Sf)

    ##Predicted bank depth
    BDpred=taue/(rho*9.8*Sf)

    ##Predicted Width
    Bpred=2*sqrt(18*k2*D/(2*a1*rho*9.8*Sf))

    print(c(indy, Sf, Dpred, BDpred, Bpred))

    #par(mfrow=c(1,1))
    ##overplot
    #if(indy==200){
    plot(msc[indy,],h[indy,],xlim=c(85,115),t='p', asp=1, xlab='y', ylab='z', cex=.5)
    #}else{
    #
    #    points (msc[indy,],h[indy,],xlim=c(85,115),t='p', asp=1, xlab='y', ylab='z', cex=.5)
        #predp= 1/(6*k2)*(msc[indy,]-msc[indy,which.min(h[indy,])] )**2 - Dpred
    #}
    predp= 1/(6*k2)*(msc[indy,]-100)**2 - Dpred
    #points(msc[indy,], predp,t='l',col=2)
    predp2= 1/(6*k2)*(seq(-Bpred/2, Bpred/2,len=200) )**2 -Dpred
    #points(seq(-Bpred/2, Bpred/2,len=200)+ msc[indy,which.min(h[indy,])], predp2,t='l',col=2)
    #points(c(-Bpred/2, -Bpred/2)+ msc[indy,which.min(h[indy,])], c(predp2[1], 0), t='l',col=2 )
    #points(c(Bpred/2, Bpred/2)+ msc[indy,which.min(h[indy,])], c(predp2[1], 0), t='l',col=2 )
    #points(c(-100, -Bpred/2)+ msc[indy,which.min(h[indy,])], c(0,0), t='l',col=2)
    #points(c(Bpred/2,1000)+ msc[indy,which.min(h[indy,])], c(0,0), t='l',col=2)

    points(seq(-Bpred/2, Bpred/2,len=200)+ 100, predp2,t='l',col=2)
    points(c(-Bpred/2, -Bpred/2)+ 100, c(predp2[1], 0), t='l',col=2 )
    points(c(Bpred/2, Bpred/2)+ 100, c(predp2[1], 0), t='l',col=2 )
    points(c(-100, -Bpred/2)+ 100, c(0,0), t='l',col=2)
    points(c(Bpred/2,1000)+ 100, c(0,0), t='l',col=2)



}
##Plot the errors
#v=which(h[indy,]< -2)
#plot(h[indy,v]-predp[v])
dev.off()

v=which(h[indy,]< -2)
#v=v[2:(length(v)-1)]
plot(h[indy,v]-predp[v], main =  mean(abs(h[indy,v]-predp[v])))

ddd=lm(h[indy,v]~I((msc[indy,v]-100)^2))
print(summary(ddd))

##########Do it again with msc4,h4
#############
#plot(msc4[indy,],h4[indy,],xlim=c(80,120),t='o')
#predp= 1/(6*k2)*(msc4[indy,]-100)**2 - Dpred
#points(msc4[indy,], predp,t='l',col=2)

##Plot the errors
#v=which(h4[indy,]< -2)
#plot(h4[indy,v]-predp[v])



