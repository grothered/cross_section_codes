# R code to compare the outputs of 2 runs.

source('../gettide1.R')

# Read in key variables
getf(110)
# Rename them
w2=w1
s2=s1
A11=A1
d11=d1
er1=er

# Change directory, and read in reference versions of the variables
setwd('reference_run')
getf(110)

w_change=range(w2-w1)
s_change=range(s2-s1)
A_change=range(A11-A1)
d_change=range(d11-d1)
er_change=range(er1-er)

print(paste(c('w_change is ', w_change, ' compared with scale ',range(w1))))
print(paste(c('s_change is ', s_change, ' compared with scale ',range(s1))))
print(paste(c('A_change is ', A_change, ' compared with scale ',range(A1))))
print(paste(c('d_change is ', d_change, ' compared with scale ',range(d1))))
print(paste(c('er_change is ', er_change, ' compared with scale ',range(er))))

