###

Slope=3/(99*50)
Q=5.
Width=10.
f=0.07
#man_n=0.03

# Can get the analytical solution from these equations:
# Sf=man_n^2*(Q^2/A^2)/(D^(4/3))
# or
# Sf = f/(8*9.8)*(Q^2/A^2)/D

# Q = A*V
# A = Width*D



# Manning
# Sf/(man_n^2*Q^2)=1/(Width^2*D^(10/3))
#D = ( (man_n^2*Q^2)/(Slope*Width^2) )^(3/10)


# Darcy
D= (f/(8*9.8)*Q^2/(Width^2*Slope))^(1/3)
V = Q/(D*Width)


# Comparison with the model
source('../../gettide1.R')
getf(100)
depth1=w1[300,90]-bt[300,90]
depth2=w1[300,70]-bt[300,70]

print(paste(' The difference between the computed and analytical depth at point 1 and 2 is ' , depth1-D, depth2-D))
print(paste(' The correct answer is ', D))
print(paste(' THE RELATIVE ERROR IS ', abs(depth1-D)/D, abs(depth2-D)/D ))
vel1=d2[300,90]/A2[300,90]
vel2=d2[300,70]/A2[300,70]

print(paste(' The difference between the computed and analytical velocity at point 1 and 2 is ' , vel1+V, vel2+V))
print(paste(' The correct answer is ', V))
print(paste(' THE RELATIVE ERROR IS ', abs(vel1+V)/V, abs(vel2+V)/V ))

