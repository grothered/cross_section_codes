## Quick R routine to create the files sectionsold2 and lnthsold2
## 
a=5
b=60
width=10.0

ys=matrix(NA, nrow=a,ncol=b)
bed=ys

for(i in 1:20){
    ys[,i]=c(-0.001, seq(0,width,len=a-2), width+0.001)
    #bed[,i]=c(100,rep(1,2), rep(-2,2), rep(-2,2), rep(1,2), 100)
    bed[,i]=c(100,rep(-2,a-2), 100)
}


for(i in 21:b){
    width=width*1.5
    ys[,i]=c(-0.001, seq(0,width,len=a-2), width+0.001)
    #bed[,i]=c(100,rep(1,2), rep(-2,2), rep(-2,2), rep(1,2), 100)
    bed[,i]=c(100,rep(-2,a-2), 100)
}


write(ys,file='lnthsold2',ncolumns=10)
write(bed,file='sectionsold2',ncolumns=10)

