getf<-function(l, ...){
    # Function to get lots of 1D variables
	a= which(dir()== "Area.1DO")
	b= which(dir()== "Discharge.1DO")
	c=which(dir()=="Width.1DO")
	d= which(dir()=="Bottom.1DO")
	w=which(dir()=="Water.1DO")
	s=which(dir()=="Susconc.1DO")
	t= which(dir()=="Times.ODO")
	e= which(dir()=="Resuspension.1DO")	
	#dee=which(dir()=="deposition")
	dee2=which(dir()=="Discharge_halftime_lim.1DO")
	ayy2=which(dir()=="Area_halftime.1DO")
	
	A1<- scan(dir()[a], what= "numeric", ...)
	d1<- scan(dir()[b], what= "numeric", ...)
	b1<- scan(dir()[c], what= "numeric", ...)
	bt<- scan(dir()[d], what="numeric", ...)	
	w1<- scan(dir()[w], what="numeric", ...)
	s1<- scan(dir()[s],what="numeric", ...)
	t1<- scan(dir()[t], what="numeric", ...)
	e1<-scan(dir()[e], what="numeric", ...)
	#de1<- scan(dir()[dee],what="numeric", ...)	
	d2<- scan(dir()[dee2],what="numeric", ...)	
	A2<- scan(dir()[ayy2],what="numeric", ...)	
	
	A1<<- matrix(as.numeric(A1),ncol=l,byrow=TRUE)
	d1<<- matrix(as.numeric(d1),ncol=l,byrow=TRUE)
	b1<<- matrix(as.numeric(b1), ncol=l,byrow=TRUE)
	bt<<- matrix(as.numeric(bt), ncol=l, byrow=TRUE)
	w1<<- matrix(as.numeric(w1),ncol=l,byrow=TRUE)	
	s1<<- matrix(as.numeric(s1),ncol=l,byrow=TRUE)	
	t1<<- as.numeric(t1)
	er<<- matrix(as.numeric(e1), ncol=l,byrow=TRUE)
	#dep<<- matrix(as.numeric(de1), ncol=l, byrow=TRUE)
	d2<<- matrix(as.numeric(d2), ncol=l, byrow=TRUE)
	#d2<<- sqrt(abs(d2[,2:(l+1)])*abs(d2[,1:l]))*0.5*(sign(d2[,2:(l+1)])+sign(d2[,1:l]))
	A2<<-matrix(as.numeric(A2), ncol=l, byrow=TRUE)
} 

anim1d<-function(d,sleep=0, ...){
	f=dim(d)
	b=range(d,na.rm=TRUE)
	for(i in 1:f[1]){
		plot(c(0,f[2]),b, col=0, ...)
		points(d[i,], ...)
		if(sleep!=0){
			Sys.sleep(sleep)
				} 
	               }
}

anim2<-function(d,x, sleep=0,...){
	f=dim(d)

	
	b=range(c(d,x),na.rm=TRUE)
	
	for(i in 1:f[1]){
		plot(c(0,f[2]),b, col=0, ...)
		points(d[i,],col=2, ...)
		points(x[i,],col=3, ...)

		if(sleep!=0){
			Sys.sleep(sleep)
				} 
	               
	               }
}

danim<-function(A,B, sleep=0,...){
    par("mfrow"=c(2,1))
    f=dim(A)
    ff=dim(B)

    a=range(A)
    b=range(B)

    for(i in 1:f[1]){
    plot(c(0,f[2]),a,col=0,...)
    points(seq(0,f[2],len=f[2]),A[i,],col=2, ...)
    abline(v=30)
    plot(c(0,ff[2]),b,col=0,...)
    points(seq(0,ff[2],len=ff[2]),B[i,],col=3, ...)
    abline(v=30)

            if(sleep!=0){
                Sys.sleep(sleep)
                    } 
    }

}


animv<-function(A1,d1,sleep=0,...){
f=min(dim(A1)[1],dim(d1)[1])


	b=range(d1[1:f,]/A1[1:f,],na.rm=TRUE)
	for(i in 1:f){
		plot(c(0,dim(A1)[2]),b, col=0, main=i)
		points(d1[i,]/A1[i,], ...)
		if(sleep!=0){
			Sys.sleep(sleep)
				} 
	               }
		

}


animts<-function(A,ta,B, tb, n=length(ta),sleep=0, ...){
    f= min(max(ta),max(tb))
    r=range(c(A,B))
    d=dim(A)
        for (i in 1:n){
        
    #	i1= max(which(ta<f/n*i))
    #	i2=max(which(tb<f/n*i))	
        i1=which.min(abs(ta-f/n*i))
        i2=which.min(abs(tb-f/n*i))
        plot(c(0,d[2]),r,col=0)
        points(A[i1,], ...)
        points(B[i2,],col=2, ...)
        if(sleep!=0){
            Sys.sleep(sleep)
                } 
    }
}




wire<-function(lnths,cs1,xincrem,t){

    library(lattice)
    y=lnths[,,t]*0
    for(i in 1:dim(y)[2]){
    y[,i]<- xincrem*(i-1)
     }


    wireframe(cs1[,,t]~lnths[,,t]*y)

}

	
get1d<-function(l){
	a=which(dir()=="Areas.out")	
	b=which(dir()=="Discharge.out")
	c=which(dir()=="times.out")
	f= which(dir()=="Areafil.out")
	g=which(dir()=="Qfil.out")
	h=which(dir()=="Width.out")
	j=which(dir()=="Rando.out")
	
	ds=scan(dir()[b],what="numeric")
	d<<-matrix(as.numeric(ds),ncol=l,byrow=TRUE)
	
	As=scan(dir()[a],what="numeric")
	A<<-matrix(as.numeric(As),ncol=l,byrow=TRUE)
	
	t<<-scan(dir()[c],what="numeric")
		

	Rs= scan(dir()[h],what="numeric")	
	B<<-matrix(as.numeric(Rs),ncol=l,byrow=TRUE)
	
	rss= scan(dir()[j],what="numeric")	
	Rando<<-matrix(as.numeric(rss),ncol=l,byrow=TRUE)	

	Afs=scan(dir()[f],what="numeric")
	Afils<<-matrix(as.numeric(Afs),ncol=l,byrow=TRUE)

	Dfs=scan(dir()[g],what="numeric")
	dfils<<-matrix(as.numeric(Dfs),ncol=l,byrow=TRUE)
	
	n=length(t)
	det<<- (-as.numeric(t[1:(n-2)])+as.numeric(t[3:n]))/2
	
}




gethy<-function(l, ...){
	a= which(dir()== "waters1")
	b= which(dir()== "discharge1")
	#c=which(dir()=="widths")
	#d= which(dir()=="bottom")
	w=which(dir()=="water")
	#s=which(dir()=="susconc")
	t= which(dir()=="times")
	#e= which(dir()=="erosion")	
	#dee=which(dir()=="deposition")

	a1<- scan(dir()[a], what= "numeric", ...)
	b1<- scan(dir()[b], what= "numeric", ...)
	#c1<- scan(dir()[c], what= "numeric", ...)
	#d2<- scan(dir()[d], what="numeric", ...)	
	w1<- scan(dir()[w], what="numeric", ...)
	#s1<- scan(dir()[s],what="numeric", ...)
	t1<- scan(dir()[t], what="numeric", ...)
	#e1<-scan(dir()[e], what="numeric", ...)
	#de1<- scan(dir()[dee],what="numeric", ...)	

	A1<<- matrix(as.numeric(a1),ncol=l,byrow=TRUE)
	d1<<- matrix(as.numeric(b1),ncol=l,byrow=TRUE)
	#b1<<- matrix(as.numeric(c1), ncol=l,byrow=TRUE)
	#bt<<- matrix(as.numeric(d2), ncol=l, byrow=TRUE)
	w1<<- matrix(as.numeric(w1),ncol=l,byrow=TRUE)	
	#s1<<- matrix(as.numeric(s1),ncol=l,byrow=TRUE)	
	t1<<- as.numeric(t1)
	#er<<- matrix(as.numeric(e1), ncol=l,byrow=TRUE)
	#dep<<- matrix(as.numeric(de1), ncol=l, byrow=TRUE)

	} 

genread<-function(name, l, ...){
	s=which(dir()==name)
	s2=scan(dir()[s], what="numeric", ...)
	lnth=floor(length(s2)/prod(l))*prod(l)
	array(as.numeric(s2)[1:lnth], dim=c(l, lnth/prod(l)))

	}

getmor<-function(l){
	#a= which(dir()== "waters1")
	#b= which(dir()== "discharge1")
	c=which(dir()=="widths")
	d= which(dir()=="bottom")
	#w=which(dir()=="water")
	s=which(dir()=="susconc")
	t= which(dir()=="times")
	e= which(dir()=="erosion")	
	dee=which(dir()=="deposition")

	#a1<- scan(dir()[a], what= "numeric")
	#b1<- scan(dir()[b], what= "numeric")
	c1<- scan(dir()[c], what= "numeric")
	d2<- scan(dir()[d], what="numeric")	
	#w1<- scan(dir()[w], what="numeric")
	s1<- scan(dir()[s],what="numeric")
	t1<- scan(dir()[t], what="numeric")
	e1<-scan(dir()[e], what="numeric")
	de1<- scan(dir()[dee],what="numeric")	

	#A1<<- matrix(as.numeric(a1),ncol=l,byrow=TRUE)
	#d1<<- matrix(as.numeric(b1),ncol=l,byrow=TRUE)
	b1<<- matrix(as.numeric(c1), ncol=l,byrow=TRUE)
	bt<<- matrix(as.numeric(d2), ncol=l, byrow=TRUE)
	#w1<<- matrix(as.numeric(w1),ncol=l,byrow=TRUE)	
	s1<<- matrix(as.numeric(s1),ncol=l,byrow=TRUE)	
	t1<<- as.numeric(t1)
	er<<- matrix(as.numeric(e1), ncol=l,byrow=TRUE)
	dep<<- matrix(as.numeric(de1), ncol=l, byrow=TRUE)

	} 







getcs<- function(l, x, ...){
	a=which(dir()=="sections")
#	b=which(dir()=="lengths")
	a1<- scan(dir()[a], what= "numeric", ...)
#	b1<- scan(dir()[b], what="numeric", ...)
	cs1<<- array(as.numeric(a1),dim=c(l,x,floor(length(a1)/(l*x))))
#	lnths<<- array(as.numeric(b1),dim=c(l,x,floor(length(a1)/(l*x))))
			}


getvs<- function(l, x, ...){
	a=which(dir()=="vels")
#	b=which(dir()=="lengths")
	a1<- scan(dir()[a], what= "numeric", ...)
#	b1<- scan(dir()[b], what="numeric", ...)
	vels<<- array(as.numeric(a1),dim=c(l,x,floor(length(a1)/(l*x))))
#	lnths<<- array(as.numeric(b1),dim=c(l,x,floor(length(a1)/(l*x))))
			}



gettendep<-function(a,b,layers, ...){
	aa=which(dir()=="tendep")

	a1<-scan(dir()[aa],what="numeric", ...)

	tendep<<-array(as.numeric(a1),dim=c(a,b,layers) )
}


########Some plotting

 cls<- function(M,zlim){
     b=dim(M)
     a1=M[1:(b[1]-1),1:(b[2]-1)]
     a2=M[1:(b[1]-1),2:(b[2])]
     a3=M[2:(b[1]),2:(b[2])]
     a4=M[2:(b[1]),1:(b[2]-1)]
     
     hs= 0.25*(a1+a2+a3+a4)
    #hs2= 50*(hs-min(hs,na.rm=TRUE))/(max(hs,na.rm=TRUE)-min(hs,na.rm=TRUE))+1
    hs2= 20*(hs-min(zlim))/(max(zlim)-min(zlim))+1

    hs2

 }
 

surf<-function(a,b,M,theta,phi, ...){
    persp(seq(1/dim(M)[1],a,len=dim(M)[1]), seq(1/dim(M)[2],b,len=dim(M)[2]),M,col=as.numeric(cls(M, ...)),theta=theta,phi=phi, scale=FALSE, shade=.5)
}

surf22<-function(a,b,lengths,M,theta=0,phi=40,colzlim=range(M),...){

    l2=lengths*0
    M2=M*0
    for(i in 1:dim(M)[2]){
    b2=approx(lengths[,i],M[,i],seq(0,max(lengths[,i]),len=length(lengths[,i]))  )
    l2[,i]=b2$x
    M2[,i]=b2$y

    }


    persp(seq(1/dim(M2)[1],a,len=dim(M2)[1]), seq(1/dim(M2)[2],b,len=dim(M2)[2]),M2,col=as.numeric(cls(M2, colzlim)),theta=theta,phi=phi, scale=FALSE,...)
}



surf23<-function(lnths,hs,delX, a=4, ...){
    library(fields)
    b<-dim(lnths)[2]
    Xs<-matrix(seq(0,b*delX,len=b), nrow=dim(lnths)[1], ncol=dim(lnths)[2],byrow=T)

    quilt.plot(c(lnths),c(Xs),c(hs), nrow=dim(lnths)[1]*a, ncol=dim(lnths)[2], ...) 

}


surf24<-function(lnths,hs,delX, smth=1,...){
    library(fields)
    b<-dim(lnths)[2]
    Xs<-matrix(seq(0,b*delX,len=b), nrow=dim(lnths)[1]*smth, ncol=dim(lnths)[2],byrow=T)

    l2=Xs*0 #lnths*0
    M2=Xs*0 #hs*0
    for(i in 1:dim(hs)[2]){
    b2=approx(lnths[,i],hs[,i],seq(0,max(lnths[,i]),len=length(lnths[,i])*smth)  )
    l2[,i]=b2$x
    M2[,i]=b2$y

    }


    quilt.plot(c(l2),c(Xs),c(M2), nrow=dim(lnths)[1]*smth, ncol=dim(lnths)[2]*smth, ...) 
    #persp(seq(1/dim(M2)[1],a,len=dim(M2)[1]), seq(1/dim(M2)[2],b,len=dim(M2)[2]),M2,col=as.numeric(cls(M2, zlim)),theta=theta,phi=phi, scale=FALSE,...)

}
#######Plotting to see if conservation laws seem to hold


#######Plotting to see if conservation laws seem to hold

#Integral of discharge over time at a single point. If it creeps beyond what can be attributed to sedimentation/erosion, then we are in trouble
conplot<- function(d,t, ...){
    a=max(length(d), length(t))

    plot(cumsum(d[2:(a-1)]*(-t[1:(a-2)]+t[3:a])*0.5), ...)

    }



#Compare the discharge coming in with the upstream change in area.
conplot2<-function(d,t,A,delX,sectno=1, indds=1:length(t)){
     conplot(d[indds,sectno],t[indds],type='l')
    xx=rowSums(A[indds,sectno:(dim(A)[2])])*delX
     points(xx-xx[1],col=2,t='l')
    xx=rowSums(A[indds,(sectno+1):(dim(A)[2])])*delX
     points(xx-xx[1],col=3,t='l')
}

##This is for use with THE OLD d2 rather than d
conplot3<-function(d2,t,A,delX,sectno=1){
    a=max(length(d2[,sectno]), length(t))
    plot(cumsum(d2[2:a,sectno]*(-t[1:(a-1)]+t[2:a])),t='l')
    xx=rowSums(A[,(sectno):(dim(A)[2])])*delX
     points(xx-xx[1],col=2,t='l')
    yy=rowSums(A[,(sectno+2):(dim(A)[2])])*delX
     points(yy-yy[1],col=3,t='l')
		}

##This is for the NEW d2, with a+1 points. It performs a time integration of discharge2, which is the centered space, centered time discharge
conplot4<-function(d2,t,sectno=1, indds=1:length(t)){
    a=length(indds)
    plot(cumsum(d2[indds[2:a],sectno]*diff(t[indds])))
    }
    ##This is for the new d2, with a+1 points. It performs a time integration of discharge 2 as in conplot4, and then compares that with a calculation of the changing volume upstream of that point.
    conplot5<-function(d2,t,A,delX,sectno=1, indds=1:length(t)){
    a=length(indds)
    par(mfrow=c(2,1))
    plot(cumsum(d2[indds[2:a],sectno]*diff(t[indds])), t='l')
    xx=rowSums(A[indds,(sectno):(dim(A)[2])])*delX
     points(xx-xx[1],col=2,t='l')
    plot( cumsum(d2[indds[2:a],sectno]*diff(t[indds])) -(xx[2:a]-xx[1]),t='l')
    par(mfrow=c(1,1))
    #yy=rowSums(A[indds,(sectno+2):(dim(A)[2])])*delX
    # points(yy-yy[1],col=3,t='l')
}

plotconz<-function(rhos=2600,voidf=0.6,...){
    #conz=matrix(scan('dropout3'),ncol=2,byrow=T)
    conz=matrix(scan('dropout3'),ncol=5,byrow=T)
    x=conz[,1]-(conz[,2]-conz[1,2]+conz[1,1]) #Discharge minus volume of water, normalized
    y=-(conz[,5]-conz[1,5])/rhos*1/(1-voidf)
    par(mfrow=c(3,1)) #Sedimentation plus erosion, normalized
    plot(abs(x/conz[,1]), ...)
    plot(x, ...)
    points(y,col=4, ...)
    plot(abs((x-y)/conz[,1]),t='l')
    #plot(conz[,1],...)
    #points(conz[,2]-conz[1,2]+conz[1,1],...)
    par(mfrow=c(1,1))
}

#local mass conservation as in Sobey
mscon<-function(t,Q,h,b,A, delX){

    a=min(length(t), length(h), length(b), dim(Q)[1] )

    dt=diff(t,lag=2)[1:(a-2)]
    dq=((Q[,3]-Q[,1])[2:(a-1)])

    dA=diff(A,lag=2)[1:(a-2)]

    #dh=diff(h)[1:a]
    #width=b[1:a]
    #plot(width*dh/dt + dq/delX)


    par(mfrow=c(2,1))
    plot(dA + dt*dq/(2*delX),t='o',cex=.3, ylim=c(-1,1)*max(abs(dA + dt*dq/(2*delX))))
    abline(h=0)
    plot( 1 + (dq*dt)/(dA*2*delX), ylim=c(-1,1),t='p',cex=.3, col=dA*10+1)
    abline(h=0)
    abline(h=.1,col=2)
    abline(h=-.1,col=2)
}


# a=length(t1)
# for(i in 2:109){
# mscon(t1[100:a],d1[100:a,(i-1):(i+1)],w1[100:a,i],b1[100:a,i],A1[100:a,i],51.64)
# }

#momentum conservation as in Sobey
#mocon<-function{

#}


##A nice contour plot. Notice how we add a plot at the end = I think we could add whatever here. 
#filled.contour(sus[,,i],color.palette = terrain.colors,zlim=c(0,7), plot.axes={axis(1);axis(2);contour(sus[,,i],add=T)})


##########Function to get rid of wetting points

wet<-function(A,b){

x= which(A>b,arr.ind=TRUE)

x2=0*A

x2[x]= 1

x2
}



#########################For getting the lateral unsteady and non-uniform term
#Assume A is the matrix, and Alast is the A from the previous time step.

##This is just useful

gn<-function(n1,n2,n3){

    gettaus(n1,n2,n3)

    getcs(n1,n2)

    getf(n2)

}


compare<-function(cw,n,l){

    out=c()
    for(i in 1:l){
    dts= sum(depths[,i]*A[,i]*cw) #cw is the lateral width of a cell-- i.e. distance between 2 points.
    d=d1[n,i]

    out=c(out,dts/abs(d))
    }
    out
}





setn<-function(ind,delX,Width){


    A<<- -(taus[,,ind]*8/(0.01*1000))^0.5  #velocity estimates at different times, assuming f=0.01
    Alast<<- -(taus[,,ind-1]*8/(0.01*1000))^0.5
    Anext<<- -(taus[,,ind+1]*8/(0.01*1000))^0.5

    ds=c()
    dslast=c()
    dsnext=c()

    for(i in 1:dim(A)[2]){
    ds=cbind(ds, pmax(w1[ind,i]-cs1[,i,1],0))

    dslast= cbind(dslast, pmax(w1[ind-1,i]-cs1[,i,1],0))

    dsnext=cbind(dsnext, pmax(w1[ind+1,i]-cs1[,i,1],0))
    }

    depths<<- ds  #depth estimates at different times
    dslast<<-dslast
    dsnext<<-dsnext

    ww<<-Width/dim(A)[1]

    Area<<- matrix(apply(depths*ww,2,sum),ncol=ncol(A),nrow=nrow(A),byrow=TRUE) #Area estimate
    Arealast<-  matrix(apply(dslast*ww,2,sum),ncol=ncol(A),nrow=nrow(A),byrow=TRUE) #Area estimate
    Areanext<-  matrix(apply(dsnext*ww,2,sum),ncol=ncol(A),nrow=nrow(A),byrow=TRUE) #Area estimate

    #delX<<-50
    delT<<- 0.5*(t1[ind+1]-t1[ind-1])



    NN(A,Alast,Anext,depths,dslast,dsnext,delX,delT, Area,Arealast,Areanext,ww)

    }


    ##For estimating that unsteady term. Note that I ignore dry points even though they may go wet in future time steps 

    NN<-function(A,Alast,Anext,depths,dslast,dsnext,delX,delT, Area,Arealast,Areanext,w){

    wt<- (depths>0) #wet points
    wtlast<-(dslast>0)
    wtnext<-(dsnext>0)

    #lth<-apply(wt,2,sum) #wetted length
    #lthlast<-apply(wtlast,2,sum)
    #lthnext<-apply(wtnext,2,sum)


    baU<<-apply(A*depths*w,2,sum)/pmax(apply(depths*w,2,sum),0.000001) #Area[1,] #pmax(lth,1) #Average velocity
    baUlast<-apply(Alast*dslast*w,2,sum)/pmax(apply(dslast*w,2,sum),0.00001) #Arealast[1,] ##
    baUnext<-apply(Anext*dsnext*w,2,sum)/pmax(apply(dsnext*w,2,sum),0.00001) #Areanext[1,] #

    mbaU<<-matrix(baU,ncol=ncol(A),nrow=nrow(A),byrow=TRUE)#*wt
    mbaUlast<<-matrix(baUlast,ncol=ncol(Alast),nrow=nrow(Alast),byrow=TRUE)#*wtlast
    mbaUnext<<- matrix(baUnext,ncol=ncol(Anext),nrow=nrow(Anext),byrow=TRUE)#*wtnext


    acU<<-  (A-mbaU)*wt ##matrix(baU,ncol=ncol(A),nrow=nrow(A),byrow=TRUE))*wt #cross sectional variation
    acUlast<<-  (Alast-mbaUlast)*wtlast  #-matrix(baUlast,ncol=ncol(Alast),nrow=nrow(Alast),byrow=TRUE))*wtlast #cross sectional variation
    acUnext<<- (Anext-mbaUnext)*wtnext  ##matrix(baUnext,ncol=ncol(Alast),nrow=nrow(Alast),byrow=TRUE))*wtnext #cross sectional variation


    f=dim(A)
    flast<-dim(Alast)


    nu<- apply( (acU^2*depths*w),2,sum)
    nu<- (-nu[1:(f[2]-2)]+ nu[3:(f[2])])/(2*delX)  #non-uniform convective term
    nu<- c(nu[1],nu,nu[f[2]-2])

    nuc<<-matrix(nu,ncol=f[2],nrow=f[1],byrow=TRUE)

    f=dim(A)
    flast<-dim(Alast)

    ##Note-- dffx is a variable that i constantly redefine-- useful for redefining things cleanly. It seems to me that R functions will always interpret a variable as its global value if they can-- so a variable redefined in a function (say just using <- ,without <<-) will not follow over.  

    dffx<- (-(acU^2*depths)[,1:(f[2]-2)]+(acU^2*depths)[,3:f[2]])/(2*delX) #
    dffx<- cbind(dffx[,1],dffx,dffx[,f[2]-2])
    dffx1<<-dffx


    dffx<-  (-mbaU[,1:(f[2]-2)] + mbaU[,3:f[2]])/(2*delX)
    dffx<- cbind(dffx[,1],dffx,dffx[,f[2]-2])
    dffx2<<- 2*depths*acU*dffx


    dffx<- (-(acU*depths)[,1:(f[2]-2)]+(acU*depths)[,3:f[2]])/(2*delX) #

    dffx<- cbind(dffx[,1],dffx,dffx[,f[2]-2])

    dffx3<<-mbaU*dffx

    dfft<<-   (acUnext*dsnext-acUlast*dslast)/(2*delT) #backward time derivative of acU


    N<<-  (nuc/pmax(Area,0.00001))*(-depths) + dfft + dffx1+dffx2+ dffx3


}


##Function to get taus

gettaus<-function(a,b,t){
    taus<-scan('tausss',nmax=a*b*t)

    taus<<-array(taus,dim=c(a,b,t))
}


##############A NICE PLOT

niceplot<-function(width, converge, cs1){

    library(lattice)

    x=cs1[,,1]*0
    y=cs1[,,1]*0
    for(i in 1:dim(cs1)[2]){
    x[,i]= 50*i
    y[,i]= seq(-width/2*exp(-i/converge), width/2*exp(-i/60), len=dim(cs1)[1])
    }

    b=seq(1,dim(cs1)[1],dim(cs1)[1]/4)
    b2=seq(1,dim(cs1)[2],dim(cs1)[2]/2)

    wireframe(cs1[b,b2,1]~x[b,b2]*y[b,b2],screen=list(z=90, x=150,y=30),aspect=c(4,0.4),drape=TRUE)

}



##fourier power
pow<-function(z){
    x= c(2*z[1],z[2:length(z)]+rev(z[2:length(z)]))/2
    Re(x)
}

surf2<-function(lnths,cs1,long){

	library(scatterplot3d)
	y=matrix(seq(0,long,len=dim(lnths)[1]),ncol=ncol(lnths),nrow=nrow(lnths))

	scatterplot3d(as.numeric(lnths),as.numeric(y),as.numeric(cs1))

	
}


#a<-function(d1,A1){
#bb=dim(d1)
#cc=dim(A1)
#l=min(dim(d1)[1],dim(A1)[1])
#anim1d(d1[1:l,]/A1[1:l,])
#}


#This looks okay -- if you can tell what it means?!!
#scatterplot3d(ab[zz],a[zz],ac[zz],color=colors()[1:120][floor(ac[zz])+6])

#for(i in 1:200){
#k=paste(10000000+floor(t1[i*200, 1]/60),".png",sep="") #The 1000 ensures that the orderis good for the animation program
#jpeg(file=eval(k), quality=100)
#surf22(100,500,lnths[,,i],cs1[,,i],1,5,zlim=c(-10,5))
#title(paste(floor(t1[i*200,1]*50/(60*60*24)),"morphological days (morf=50)"))
#dev.off()
# }

##then on command prompt we type convert *.jepg animation.gif

#par("din")= c(10,10)

#for(i in 1:78){
#k=paste(1200+i,".png",sep="")
#png(file=eval(k))
#surf(400,5000,cs1[,,i],1,3,zlim=c(-5,5))
#dev.off()
#}


#could try par("din") to alter the dimensions and make nice plots

