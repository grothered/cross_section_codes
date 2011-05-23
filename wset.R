# This contains various utility functions for testing / checking / starting the single
# cross-section code.



wset<-function(d){
    # Function to estimate settling velocity with 3 methods

    # d=grain size

    ######## Input parameters
    Rs= 1.65 #Submerged specific gravity of sediment = specific grav -1
    g = 9.8 # Gravity
    kvis = 1.0e-06  # Kinematic water viscosity


    ######## Jimenez and Madsen (2003) settling velocity
    Sstar= d/(4*10^-6)*((Rs)*g*d)^.5
    AA= .954 #Table 1 P=3.5 (natural sediment)
    BB= 5.121 #
    Wstar= 1/(AA+BB/Sstar)
    vs1= Wstar*(Rs*9.8*d)^.5

    ####### Settling velocity, using Rubey (1933)
    eta=0.001 # Viscosity of water (not kinematic). I tested this by comparing my plot to theirs
    vs2= ( sqrt( 4/3*g*1000*(2600-1000)*(d/2)^3 +9*eta^2  ) -3*eta )/(1000*(d/2))

    ###### van Rijn (1984) Sediment Transport, Part II: Suspended Load Transport, Journal Hydraulic Eng 110. 
    # FIXME: NOTE THAT VAN RIJN'S FORMULAE ARE DISCONTINUOUS AT THE BREAK POINTS
    if(d<0.0001) vs3 = 1/18*(Rs*g*d^2)/kvis

    if((d>= 0.0001)&(d<0.001)) vs3 = 10*kvis/d*( (1+ (0.01/kvis^2)*Rs*g*d^3)^(0.5) -1)

    if(d>= 0.001) vs3 = 1.1*(Rs*g*d)^(0.5)

    # Print output 
    output = c(vs1,vs2, vs3)
    names(output) = c('Jimenez and Madsen (2003)', 'Rubey (1933)', 'van Rijn (1884) - discontinuities at 0.1 and 1 mm')

    print(output)
}


tauc<-function(d){
    #R function to estimate the critical shear stress for non-cohesive sediments, following van Rijn (2007a), Equation 2 
    # d = grain size

    if(d<0.00007){
    print('d is small, think about using another method')
    }
    # "Standard" values of physical parameters suitable for a first guess
    rho_s=2600.
    rho_w=1000.
    v=1.0e-06
    g = 9.8
    s = rho_s/rho_w

    d_star = d*( (s-1)*g/v^2)**(1/3)

    # According to van Rijn, 1984 bedload paper, Fig 1
    if(d_star<4){
        theta_c = 0.115*(d_star)^(-0.5)
    } else if((d_star>=4.0)&(d_star<10.0)){
            theta_c = 0.14*d_star^(-0.64)
    } else if((d_star>=10.)&(d_star<20.0)){
            theta_c = 0.04*d_star^(-0.1)
    } else if((d_star>=20.)&(d_star<150.0)){
            theta_c = 0.013*d_star^(0.29)
    } else{
            theta_c = 0.055
    }

    tauc = (rho_s-rho_w)*g*d*theta_c

    output = c(tauc, theta_c)
    names(output) = c('Taucrit', 'Theta_crit')
    output
}


f_to_ks<-function(f, depth){
    # Function to convert an average darcy weisbach 'f' value
    # to a ks value. Note that the 'log' terms in the ks equation are log10!
    #ks1 = 12*depth/exp( sqrt(8*9.8/f)/18 )
    ks2 = 12*depth/10^( sqrt(8*9.8/f)/18 )
    print(paste('ks = ',ks2))
        }


f_colebrook<-function(f_g,depth, vel, d50){
    # Function to study the convergence of direct iteration for the colebrook equation

    k_sg = 10.0*d50
    #Re = max((sqrt(f_g/8.0_dp)*abs(vel))*k_sg/1.0e-06_dp, 10.0_dp)
    Re = abs(vel)*max(depth,20.*k_sg)/1.0e-06

    #Turbulent regime
    for(i in 1:10){
        f_g = 0.25/(log( k_sg/(3.71*4.0*max(depth,20*k_sg)) + 2.51/(Re*sqrt(f_g)),10 )  )**2
    }
    f_g

}


vrijn_bed<-function(vel, h, d50, d90){
    # Function to implement van rijn's (2007) bedload equation

    #f_g =4.0*(8.*9.8/(18.*log10(12.*pmax(h, 20*d90)/(1*d90)+0.0)+0.0e+00)^2)
    f_g =(8.*9.8/(18.*log10(12.*pmax(h, 20*d90)/d90))^2)
    tau = 1000*f_g/8*vel^2
    # note -- in the paper, it is written that tau = rho f/2 vel^2 -- typo, the
    # way it is written here is standard, and gives much better predictions.
    taucrit=as.numeric(tauc(d50)[1])
    Qbed= 0.5*2600*max(0.000062/d50,1.)*d50*(d50*(1.6*9.8/(1.0e-06)**2)**(1./3.))**(-0.3)* 
            1000**(-0.5)*(abs(tau))**(0.5)*
            pmax( abs(tau)-taucrit,0.0)/taucrit*sign(tau)
    cbind(tau, taucrit, Qbed)

            
}

test_vrijn_bed<-function(){
    # Compares the predictions of the vanrijn bedload formula with some data
    # presented in van Rijn (2007a)

    #Measurement data from Table 2 of vanrijn 2007a
    d50 = c(1050, 1050, 950, 530, 530, 530, 530, 530, 690, 400, 400, 400, 400, 400)/1e+06
    d90 = c(1750, 1750, 1550, 700, 700, 700, 900, 1300, 1270, 1000, 1000, 1000, 1000, 1000)/1e+06
    h = c(9, 9.8, 9.5, 4.5, 4.5, 4.5, 4.5, 4.5, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0)
    vel = c(1.4, 1.55, 1.15, 0.45, 0.6, 0.82, 1, 1.2, 0.8, 0.33, 0.45, 0.6, 0.7, 0.83)

    qb = c(0.1, 0.18, 0.06, 0.002, 0.005, 0.031, 0.05, 0.065, 0.035, 0.0005, 0.002, 0.006, 0.012, 0.023)

    qb_pred=qb*0
    tauc = qb*0
    # Compare predictions with measurements
    for(i in 1:length(d50)){
        tmp = vrijn_bed(vel[i], h[i], d50[i], d90[i])
        qb_pred[i] = tmp[3]
        tauc[i] = tmp[1]
    }
     
    plot(qb_pred, qb, main='Predicted and measured bedload transport ', cex=d50*1e+03)
    abline(0,1)
    cbind(tauc, qb_pred,qb)
}


einstein_j1<-function(z,E, n=10){
    # Approximate the einstein integral J1 (Guo and Julien, 2004)
    # This arises because the average suspended sediment concentration
    # Cbar = C_nearbed*J1*1/( (1-deltab)/deltab)^z
    
    # Apply only to non-integer values
    if(abs(z-round(z))<0.0005) z = round(z)+0.0005

    # Compute 1/((1-deltab)/deltab)^z
    db_const = ((1-E)/E)^(-z)
    
    # Compute F1, eqn 8 in their paper
    k = seq(1,n)
    E2 = E/(1-E)
    F1 = ((1-E)^z)/E^(z-1) - z*sum( ((-1)^k)/(k-z)*(E2)^(k-z))
    J1 = z*pi/sin(z*pi) -F1
    J1*db_const
}


test_susdist<-function(ys, bed, water, Cbed, Es, wset, qby, aref, ustar, num_z =5000){
    # Function to check the numerical solution of the suspended sediment
    # distribution equation when the channel is at equilibrium:
    #
    # d(F_l)/dy = Es - Ds
    #
    # and equivalently:
    #
    # F_l + qb_l = 0
    #
    # We do this by taking the model output, using it to estimate the
    # appropriate terms, and then seeing if they balance.  An important aspect
    # is that this routine is coded independently of the original model (so we
    # do not copy-and-paste previous errors), and uses slightly different
    # calculation procedures.

    # ys = y value
    # bed = bed elevation
    # water = water surface elevation
    # Cbed = near bed suspended sediment concentration
    # Es = resuspension rate
    # wset = sediment settling velocity
    # qby = lateral bedload transport rate
    # aref = van_rijun reference height for sediment
    # num_z = number of discrete z points used to discretize the vertical
    #         coordinate.


    # a = 2
    # tmp = test_susdist( ys[a,],h[a,],0.0, Cbed[a,],Qe[a,], 0.014,Qby[a,],a_ref[a,], sqrt(tau[a,]/1026) )

    # Deposition rate    
    Ds = wset*Cbed

    # Compute:
    # F_l = int_{z = bed}^{z = water surface} (epsy*dc/dy) dz
    # where:
    # c = Cbed*f(z)
    # where f(z) defines the vertical distribution of suspended load


    # In the original model, we compute dc/dy for every z using the chain rule.
    # This allows us to to write it in terms of derivatives of variables which
    # vary with y only -- thus, we do not need to 'numerically' differentiate
    # for every z. Good to do things differently in this function, to add
    # diversity and robustness to the check.

    # Represent c(z,y) as a matrix of the form c[zind,yind]. Note that zind,
    # yind correspond to particular z and y values, and z =0 at some arbitrary
    # datum (which is the same for every y)
    c = matrix(NA,ncol=length(bed), nrow = num_z)
    zs = seq(min(bed), water, len=num_z) # z coordinate
    f = zs*NA
    for (i in 1:length(ys)){
        print(i)     
        #for(j in 1:num_z){
            #print(c(i,j))
        f = Rouse(zs, bed[i], water, aref[i], wset, ustar[i])
        #}
        c[, i] = Cbed[i]*f

    }

    # Calculate the lateral eddy viscosity, with a Parabolic model
    epsy = c*NA
    for (i in 1:length(ys)){
        parabola = pmax(zs-bed[i],0.0)*
                   pmin(water-zs,(water-bed[i]))
        epsy[,i] = 0.4*ustar[i]*(water-bed[i])*parabola/(.25*(water-bed[i])^2)
    }
    

    # Calculate dc/dy_(i+1/2)
    dcdy_h = matrix(NA,ncol=length(bed)-1,nrow=num_z)
    epsy_h= dcdy_h
    for(i in 1:(length(bed)-1)){
        dcdy_h[,i] = (c[,i+1] -c[,i])/(ys[i+1]-ys[i]) 
        epsy_h[,i] = 0.5*(epsy[,i+1]+epsy[,i])
        
    }

    # Finally, calculate integrated lateral flux
    integrand = dcdy_h*epsy_h
    dz = zs[2]-zs[1]#diff(zs)
    Fl_h = dcdy_h[1,]*NA
    for(i in 1:(length(bed)-1)){
        Fl_h[i] =  # Trapezoidal integration
            -0.5*sum((integrand[1:(num_z-1),i]*dz)+(integrand[2:(num_z),i]*dz), na.rm=T)
    }

    Fl_h
     
}



Rouse<-function(z,bed, water,aref,wset,ustar){
    # Compute the Rouse profile for suspended sediment
    
    f = z*NA # Predefine the profile

    # Find indexes where f has real values (i.e. above the bed)
    inds = which(z >= bed + aref)
    if((length(inds)>0)&(ustar>0)){

        zstar = wset/(ustar*0.4) 

        # f  
        f[inds] =( ( ( (water - bed) - (z[inds] - bed) )/(z[inds]-bed) )/
                 ( ( (water - bed) - aref )/ (aref)  ) )^(zstar)

        if(any(is.na(f[inds]))){

            print(cbind(inds,z[inds], f[inds]))
            print(c(bed, ustar, water, zstar, aref))
            stop('f is nan when it should not be')
            
        }
    }

    f
    
    }
        

