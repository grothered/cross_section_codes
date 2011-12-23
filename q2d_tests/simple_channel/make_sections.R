## Quick R routine to create the files sectionsold2 and lnthsold2
## 
a=10
b=30

ys=matrix(NA, nrow=10,ncol=30)
bed=ys

for(i in 1:b){
    ys[,i]=seq(0,20,len=10)
    bed[,i]=c(100,rep(1,2), rep(-2,2), rep(-2,2), rep(1,2), 100)
}

write(ys,file='lnthsold2',ncolumns=10)
write(bed,file='sectionsold2',ncolumns=10)

