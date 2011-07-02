# R source for getting data

get<-function(l, gz=F, ...){

	if(!gz){
		a=which(dir()=="taus")
		b=which(dir()=="bed")
		c=which(dir()=="ys")

	
	
	hs=scan(dir()[b],what="numeric", ...)


	ts=scan(dir()[a],what="numeric", ...)

	misc= scan(dir()[c],what="numeric", ...)

#	zz=min(, floor(length(ts)/l), floor(length(misc)/l) )

		}else{
		a=which(dir()=="taus.gz")
		b=which(dir()=="bed.gz")
		c=which(dir()=="ys.gz")

	
	#x2=scan(gzfile('bed.gz', open="rb"))
	
	hs=scan(gzfile(dir()[b], open='rb') ,what="numeric", ...)


	ts=scan(gzfile(dir()[a],open='rb'), what="numeric", ...)

	misc= scan(gzfile(dir()[c],open='rb'), what="numeric", ...)

		}

	h<<-matrix(as.numeric(hs)[1:(floor(length(hs)/l)*l)],ncol=l,byrow=TRUE)
	tau<<-matrix(as.numeric(ts)[1:(floor(length(ts)/l)*l)],ncol=l,byrow=TRUE)
	ys<<- matrix(as.numeric(misc)[1:(floor(length(misc)/l)*l)], ncol=l, byrow=TRUE)


		}


#getsed<-function(l, ...){
#	s=which(dir()=="sed")
#	s2=scan(dir()[s],what="numeric", ...)
#	sed<<-matrix(as.numeric(s2)[1:(floor(length(s2)/l)*l)], ncol=l, byrow=T)
#
#}


genread<-function(name, l){
	s=which(dir()==name)
	s2=scan(dir()[s], what="numeric")
	gen=matrix(as.numeric(s2)[1:(floor(length(s2)/l)*l)], ncol=l, byrow=T)
    gen
	}

#animation function
anim<-function(d,j,...){
	f=dim(d)
	b=range(d,na.rm=TRUE)
	for(i in 1:f[1]){
		plot(c(0,j),b, col=0,asp=1)
		points(seq(j/f[2],j,len=f[2]),d[i,], ...)
	               }
		}



anim2<-function(lengths,bed, sleeper=0, ...){
	f=range(bed,na.rm=TRUE)
	b=range(lengths,na.rm=TRUE)
	for(i in 1:dim(bed)[1]){
		plot(b,f, col=0,asp=1, ...)
		points(lengths[i,],bed[i,], ...)
        title(i)
        if(sleeper>0){
            slptime=paste('sleep', sleeper)
            system(slptime)
                }
	               }
		}

danim2<-function(lengths1,bed1,lengths2,bed2,sleeper=0,...){
	f=range(bed1,na.rm=TRUE)
	b=range(lengths1,na.rm=TRUE)
	for(i in 1:dim(bed1)[1]){
		plot(b,f, col=0,asp=1, ...)
		points(lengths1[i,],bed1[i,], col=2,...)
		points(lengths2[i,],bed2[i,], col=3,...)
        if(sleeper>0){
            slptime=paste('sleep', sleeper)
            system(slptime)
                }
		
	               }
		}



animtendep<-function(lengths,bed,tendeps,...){

	f=range(tendeps,na.rm=TRUE)
	b=range(lengths,na.rm=TRUE)
	for(i in 1:(dim(bed)[1])){
		plot(b,f, col=0)
		points(lengths[i,],bed[i,], t="l",lwd=4)
			for(j in 1:dim(tendeps)[3]){
			points(lengths[i,],tendeps[i,,j],t="l",col=j+1)

						}
	
		
	               }
		}
		
gettendep<-function(l, no){

	b=scan("misc2")
	numrows= floor(length(b)/(l*no))
	b1=matrix(b,ncol=l, nrow= numrows*no, byrow=TRUE)

	tendep1=array(0,dim=c(numrows,l,no))
	for(i in 1:no){
	tendep1[,,i]= b1[seq(1,numrows*no,no)+ i-1,]

	}

tendep<<-tendep1

}


wdth<-function(ys,h, level){
    ##A function to calculate the top width at a given water level. 

    a=dim(h)[1]
    a=min(a,dim(ys)[1])
    l=dim(h)[2]

    mwidths=c()
    mdepths=c()

    for(i in 1:a){

        if( ((max(h[i,])>level)&(min(h[i,])<level)) ){

            low=which((h[i,1:l-1]>level)&(h[i,2:l]<level))+1
            low = min(low) # A fix for some rare problematic cases
            high= which((h[i,1:l-1]<level)&(h[i,2:l]>level))
            high = max(high) # Fix for rare cases
        }else{
            low=1
            high=l
        }
        #print(c(low, high))


        mwidths=c(mwidths,ys[i,high]-ys[i,low])

        mdepths=c(mdepths,sum((level-h[i,low:high])*.5*(c(0,diff(ys[i,(low):high])) +c(diff(ys[i,low:high]), 0) ))/(ys[i,high]-ys[i,low]) ) 


    }


    cbind(mwidths, mdepths)
}



dimsum<-function(msc,h,no){

pts=seq(1,dim(h)[1],no)

wdths=c()
dpths=c()
eyes=c()

for(i in pts){

aa=chan_pts(msc[i,],h[i,], 1)

if(length(aa)>0){
eyes=c(eyes,i)
wdths=c(wdths, msc[i,max(aa)]-msc[i,min(aa)])
dpths=c(dpths,min(h[i,]))
}

}

cbind(wdths,dpths)

}

#################
############
########


chan_pts<-function(ys, hs, wdth){
    #Function to identify channel points using the approach of D'Alpos et al.

    cords<- cbind(ys,hs) 

    cu<-curv(cords, wdth+cords[,1]*0, 1+cords[,1]*0, 0+cords[,1]*0)  #note the trick to get things to be the right dimension.

    b=union(which((abs(cu$curvat))>.1), which(hs< mean(hs[1:10])-0.2))

    b
}

##########################
###################
########

curv<-function(coords,widths,windowmult,windowconst){
    #Function to calculate curvature!!


    ##Need to define coords and widths first

    cntrline=coords  ##Note that this variable name is redefined later-- it is a useful way to represent the channel. 

    #calculate distance between points in coords, then resample evenly

    d=((cntrline[1:(length(cntrline[,1])-1),1]-cntrline[2:length(cntrline[,1]),1])^2)
    d2=((cntrline[1:(length(cntrline[,1])-1),2]-cntrline[2:length(cntrline[,1]),2])^2)
    dist= (d+d2)^0.5
    usdist=c(0,cumsum(dist))   #upstream distance.

    #Still need to evenly space them. This is done by fitting a spline through the x and y separately, as a function of usdist. 

    spacingcoeff=1
    spacing= spacingcoeff*mean(dist)  #the proposed distance between new coordinates. 

    newx= (spline(usdist, cntrline[,1],n= round(usdist[length(usdist)]/spacing+1)))

    newy= (spline(usdist, cntrline[,2],n= round(usdist[length(usdist)]/spacing+1)))

    newwidths= (spline(usdist, widths,n= round(usdist[length(usdist)]/spacing+1)))


    x=newx[[2]]
    y=newy[[2]]

    nwidths= newwidths[[2]]

    #evenly spaced coordinates....and upstream distance Iusdist. 

    Iusdist=newx[[1]]

    cntrline=matrix(data=c(x,y),ncol=2) #This doesn't appear to be being used below. 

    ##So now we are evenly spaced



    ###Step 2 and 3 
    error= 4*spacing  ##This is a little thing that we add to the spacing between the points used to evaluate the derivatives. 

    lowerb= pmax(0,Iusdist-windowmult*nwidths-0*error-windowconst*spacing)  #The term 'spacing/2' is just useful for finding the appropriate index. 
    upperb=pmin(max(Iusdist),Iusdist+windowmult*nwidths+0*error+windowconst*spacing)

    lowerb.ind= findInterval(lowerb, Iusdist) #This is the index of the lower bound on the region that is used in the extimation of the curvature
    upperb.ind=findInterval(upperb, Iusdist) #As above, but the upper bound. 


    ###Step 4. ----
    diff1= coords-coords[lowerb.ind,] #The 'backward' difference. By the mean value theorem, for a continuous function, there exists a point with derivative equal to that approximated as (this difference)/(distance between point and lowerb.ind), and this point is located somewhere between the point where the curvature is  being calculated and lowerb.ind
    diff2= coords[upperb.ind,]-coords #The 'forward' difference. 
    diff3= coords[upperb.ind,]- coords[lowerb.ind,] #The 'central' difference. This can be used to calculate a derivative that occurs somewhere between lowerb,ind and upperb.ind...maybe is is better to have this located more closely to the point where curvature is being calculated??

    dxdyds= diff3/(Iusdist[upperb.ind]-Iusdist[lowerb.ind])  #derivatives of x and y, approximated as the change in x and y at the upper and lower points, divided by the along channel distance between them

    d2xnyds= ( (diff2/(Iusdist[upperb.ind]-Iusdist)) - (diff1/(Iusdist-Iusdist[lowerb.ind]))  )/( 0.5*(Iusdist[upperb.ind]-Iusdist[lowerb.ind]))  #Notice how the second distance is the 1/2 the distance between the 2 endpoints- as our derivatives in this one are sort of being calculated at the mid points of the point and lower, and the point and upper. 

    curvat= (dxdyds[,1]*d2xnyds[,2]- dxdyds[,2]*d2xnyds[,1])/( (dxdyds[,1]^2 +dxdyds[,2]^2)^1.5)

    curvat= c(curvat[2],curvat[2:(length(curvat)-1)],curvat[length(curvat)-1])

    list(curvat=curvat,nwidths= nwidths, usdist=Iusdist)
		}



get_extra<-function(l){
        # Function to get a range of output variables
		taug_file=which(dir()=="taug")
		Cbed_file=which(dir()=="Cbed")
		vel_file=which(dir()=="vel")
        Qe_file = which(dir()=='qe')
        Qbed_file = which(dir()=='Qbed')
        qb_G_file = which(dir()=='qby')
        aref_file = which(dir()=='a_ref')
        #time_file = which(dir()=='time')

	
        taug <<- genread(dir()[taug_file],l)	
        Cbed <<- genread(dir()[Cbed_file],l)	
        vel  <<- genread(dir()[vel_file],l)	
        Qe   <<- genread(dir()[Qe_file],l)	
        Qbed <<- genread(dir()[Qbed_file],l)	
        Qby <<- genread(dir()[qb_G_file],l+1)	
        a_ref<<- genread(dir()[aref_file],l)
        #time1 <<- scan(dir()[time_file])
}


X_conserve<-function(ind, wset=0.016){
    # Function to compare downslope bedload with deposition less erosion
    # Only valid for an even grid
    qbl = diff(Qby[ind,])/msc[1,2]
    D_les_E = -Qe[ind,] + Cbed[ind,]*wset/2600
    
    cbind(qbl, D_les_E)
}
