### R code to figure out coefficients in cubic interpolation of 4 evenly spaced
### (x,y) points
#   (0, y0), (dx, y1), (2*dx, y2), (3*dx, y3)
# And more particularly, to compute the slope at x = 3/2*dx

# Say poly is
#  y =  a0 + a1*x + a2*x^2 + a3*x^3
# This must pass through the 4 points above.
# Then defining:
# b = ( a0, a1*dx, a2*dx^2, a3*dx^3)
# It follows that:
# A*b = (y0, y1, y2, y3)
# where A is the following matrix
A = matrix(0,nrow=4,ncol=4)

A[1,] = c(1,0,0,0)
A[2,] = c(1,1,1,1)
A[3,] = c(1,2,4,8)
A[4,] = c(1,3,9,27)

# If we calculate Ainverse
Ainv = solve(A)

# Then we can get the coefficients in b in terms of (y0, y1, y2, y3)

# We would like the slope (at 3/2 dx). 
# The slope must be:
# Slope = a1 + 2*a2*x + 3*a3*x^2
# Thus,
# Slope(3/2 dx) = a1 + 3*a2*dx + 27/4*a3*dx^2
# Hence, it follows that
# Slope(3/2 dx)*dx = a1*dx + 3*a2*dx^2 + 27/4*a3*dx^3
# and this form is more easy to compute from 'b'.

# So we have that
# Slope*dx = 
# Ainv[2,].Y + 3*Ainv[3,].Y + 27/4*Ainv[4,].Y
# where Y = (y0, y1, y2, y3)
# The coefficients of the ys are:
# y0
Ainv[2,1] + 3*Ainv[3,1] + 27/4*Ainv[4,1]
# y1
Ainv[2,2] + 3*Ainv[3,2] + 27/4*Ainv[4,2]
# y2
Ainv[2,3] + 3*Ainv[3,3] + 27/4*Ainv[4,3]
# y3
Ainv[2,4] + 3*Ainv[3,4] + 27/4*Ainv[4,4]




