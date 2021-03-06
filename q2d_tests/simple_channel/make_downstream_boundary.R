## R code to make the tidal boundary condition

# Time
t = seq(0, 480*3600,by=60)
# Mouth height
h = 0.8*sin(2*pi*t/(3600*12.4))

out=cbind(t,h)

write.table( out, file='mouth', row.names=F, col.names=F)
