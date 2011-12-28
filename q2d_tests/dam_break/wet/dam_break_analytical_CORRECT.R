##Now, in using the dam break solution of Zoppou and Roberts (2003), I have sometimes found strange results. Thus, in the following code, I first implement that solution, and then I implement the solution of Wu et al., (1999) attributed to Stoker. This suggests that the variable h2 in the Zoppou and Roberts solution is wrong IN SOME CASES ONLY.

############################
##Zoppou and Roberts solution - SOMETIMES INCORRECT
###########################

l=880 #Number of x points

t1=seq(0,1000,len=1000) #Range of times that we might like to plot. Actually we only will plot t1[j]
j=91 #The plot

delX=10 #Spatial increment
g=9.8
h1=10
h0=0.1
xz=(seq(1,l)-l/2)*delX

u=rep(0,l)
h=rep(h1,l)

#Now, to get started on the analytical solution, (Zoppou and Roberts, 2003), we need to estimate the shock speed. This requires solving a nonlinear equation.
S2fun<-function(S2){
abs(-S2+2*sqrt(g*h1)+g/(4*S2)*(h0+sqrt(h0^2+8*h0*S2^2/(g)) ) -(2*g*sqrt(h0^2+8*h0*S2^2/(g)) -2*g*h0 )^0.5)
}
a11=optimize(S2fun,interval=c(1.0E-9,1000),tol=1E-12)
S2=a11$minimum
##So now we have the shock speed. We need u2, and h2
u2=S2-g/(4*S2)*(h0+sqrt(h0^2+8*h0*S2^2/(g)))

#h2 also needs minimizing --- however, I think this gives incorrect values
h2fun<-function(h2){
#abs(-h2 + h0/2*sqrt(1+8*(2*h2/(h2-h0)*(sqrt(h1)-sqrt(h2))/sqrt(h0) )^2  ) -1/2)
abs(-h2 + 1/2*sqrt(h0^2+8*h0*(2*h2/(h2-h0)*(sqrt(h1)-sqrt(h2)) )^2  ) -1/2)
}
a11=optimize(h2fun,interval=c(h0+1E-9,h1)) #This interval is important, should be sensible I guess
h2=a11$minimum

#Now find one of the intervals
b=which( (xz> -t1[j]*sqrt(g*h1))&(xz <= t1[j]*(u2-sqrt(g*h2)) ) )

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

#Now, I think the above solution is in error -- I suspect there must be an error in the transcription of Zoppou and Roberts. So let's store it here, then correct later.
uerr=u
herr=h

############################################################################################
################Alternative approach, Stoker, for u2, h2, epdot. Good for checking the above, because there is an error
############################################################################################

##To calculate the Stoker solution, we minimize this function, Wu et al., (1999)
alvars<-function(v){
#u2 is velocity in shock zone
u2=v[1]
#h2 is depth
h2=v[2]
#epdot is wave speed
epdot=v[3]

#Preliminary celerities
C0=sqrt(g*h0)
C1=sqrt(g*h1)
C2=sqrt(g*h2)

#See eqns 44 - 47 in Wu et al. (1999) - these should all be zero for the solution
eq1=epdot/C0 - 1/4*C0/epdot*(1+sqrt(1+8*(epdot/C0)^2)) -u2/C0
eq2= 1/sqrt(2)*(sqrt(1+8*(epdot/C0)^2)-1  )^0.5 -C2/C0
eq3= u2+2*(C2-C1)

abs(eq1)+abs(eq2)+abs(eq3)

}

#Find the minimum of alvars- practically, this should set all equations to 0. For initial values, we use the previous estimates of u2, h2 and S2. 
a11=optim(c(u2,h2,S2),alvars,lower=c(0,h0,0),method='L-BFGS-B')
u2=a11$par[1]
h2=a11$par[2]
S2=a11$par[3]

#At present (6/6/10) this suggests that my implementation of the solution in Zoppou and Roberts is wrong, actually h2 is wrong -- indeed I get funny plots sometimes
#Now, using these values of u2, h2, and S2, we can recalculate the solution.

u=rep(0,l)
h=rep(h1,l)
b=which( (xz> -t1[j]*sqrt(g*h1))&(xz <= t1[j]*(u2-sqrt(g*h2)) ) )
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

##So this is the correct solution, clearly -- e.g. try h0=10, h1=0.1
plot(xz,u,t='l')
points(xz,uerr,t='l',col=2)
