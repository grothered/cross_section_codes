#R code to solve the Shiono and Knight analytical equation in a trapezoidal channel -- useful for analytical comparison purposes. 

#Say the channel has top width 4, depth 1, bed width 2, side slope 1.

#Ycoordinate
ll=250 #Number of points will be 4*ll
y=seq(0,4,len=ll*4)
#Depth
H=c(seq(0,1,len=ll), seq(1,1,len=2*ll),seq(1,0,len=ll))
#Predefine velocity
Ud2=H*0

#Set parameters
g=9.8
f=0.02
S0=0.001
lambda=0.24

##To solve, we need to stitch together 3 analytical solutions, but because of symmetry, there will only be two. On the slope, the solution has the form
# Ud2=( A3*H^alph1 + A4*H^(-alph1-1) +w*H ) 
#where 
alph1=-0.5+0.5*(1 + 1*(1+1^2)^0.5/lambda*(8*f)^0.5)^0.5
w=g*S0/( (1+1^2)^0.5/1*f/8 - lambda/1^2*(f/8)^0.5 )
#This means that A4=0, to ensure that Ud2 is finite as H-->0

#Ud2=( A3*H^alph1 + w*H )
#Ud2[27:101]=0 

#We need to match this with the value and velocity gradient in the channel centre, where
#Ud2=(A1*exp(gma*y)+A2*exp(-gma*y)+8*g*S0*H/f)
gma=sqrt(2/lambda)*(f/8)^0.25*1/H
#and also impose a symmetry condition at the channel centre (d/dy Ud2 == 0)
#d/dy(Ud2) = A1*gma*exp(gma*2)-A2*gma*exp(-gma*2)
#So A2 = A1*exp(gma*2)/exp(-gma*2) =A1*v1 where
v1=exp(gma[ll]*2)/exp(-gma[ll]*2)

#To match the gradient condition at the break in slope
# A3*alph1*1 + w = A1*gma*exp(gma*1)-A2*gma*exp(-gma*1)
# A3*alph1*1 + w = A1*(gma[27]*exp(gma[27]*1)-v1*gma[27]*exp(-gma[27]*1) )

#The continuity leads to A3+w = A1*exp(gma)+A2*exp(-gma)+8*g*S0*1/f 
# A3+w=A1*(exp(gma)+ v1*exp(-gma)) +8*g*S0*1/f

#Let's turn these 2 equations into a matrix and solve:
# |a11,a12| * [A1, A3]^T = [ w , w-8*g*S0*1/f ]
# |a21,a22|

a11=gma[ll]*exp(gma[ll]*1) -v1*gma[ll]*exp(-gma[ll])
a12=-alph1
a21=exp(gma[ll])+v1*exp(-gma[ll])
a22=-1

#Solve these equations
Amat=matrix(c(a11,a12,a21,a22),ncol=2,nrow=2,byrow=T)
ANS=solve(Amat,c(w, w-8*g*S0*1/f))

#Define coefficients
A1=ANS[1]
A3=ANS[2]
A2 = A1*v1
#A4=0

#The first side slope value
Ud2=( A3*H^alph1 + w*H )
#The channel centre value
Indz=(ll+1):(3*ll)
Ud2[Indz]=(A1*exp((gma*y)[Indz])+A2*exp((-gma*y)[Indz])+8*g*S0*H[Indz]/f)
Ud2[(max(Indz)+1):(4*ll)]=rev(Ud2[1:ll]) 

geom=approx(y,-H,n=1000)

write(geom$y,file='hes',ncol=1000)
write(geom$y-1000,file='tendeps',ncol=1000)
write(geom$x,file='lnths',ncol=1000)
write(geom$x,file='tndplnths',ncol=1000)

source('../../cros.R')
get(1000)
png(file='Analytical_compare.png',width=5.4,height=5.4,res=1000,units='in')
plot(ys[3,],sqrt(tau[3,]/(1026*f/8)) ,col=2)
points(y,sqrt(Ud2),t='o')
legend('bottom',c('Analytical','Numerical'),pch=c(1,1),col=c(1,2))
dev.off()
