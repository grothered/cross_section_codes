#R function to estimate settling velocity with 2 methods
wset<-function(d){
    ######## Input parameters
    Rs= 1.65 #Submerged specific gravity of sediment = specific grav -1
    g = 9.8 # Gravity
    kvis = 1.0e-06  # Kinematic water viscosity
    #d=0.0002#grain size

    #vs=0.43*(Rs*9.8*d)^.5 # The constant out the front could be something else!

    #####The end result is highly highly sensitive to the settling velocity.

    #########A better settling velocity formula
    ##############
    ##In Ikeda and Izumi, it is not always clear how to calculate the settling velocity. Here I implement Jimenez and Madsen (2003). They quote Rubey (1933).
    Sstar= d/(4*10^-6)*((Rs)*g*d)^.5
    AA= .954 #Table 1 P=3.5 (natural sediment)
    BB= 5.121 #
    Wstar= 1/(AA+BB/Sstar)
    vs1= Wstar*(Rs*9.8*d)^.5

    ##Settling velocity, using Rubey (1933)
    eta=0.001 #Viscosity of water (not kinematic). I tested this by comparing my plot to theirs
    vs2= ( sqrt( 4/3*g*1000*(2600-1000)*(d/2)^3 +9*eta^2  ) -3*eta )/(1000*(d/2))

    # Settling velocity, using van Rijn (1984) Sediment Transport, Part II: Suspended Load Transport, Journal Hydraulic Eng 110. 
    # NOTE THAT VAN RIJN'S FORMULAE ARE DISCONTINUOUS AT THE BREAK POINTS
    if(d<0.0001) vs3 = 1/18*(Rs*g*d^2)/kvis

    if((d>= 0.0001)&(d<0.001)) vs3 = 10*kvis/d*( (1+ (0.01/kvis^2)*Rs*g*d^3)^(0.5) -1)

    if(d>= 0.001) vs3 = 1.1*(Rs*g*d)^(0.5)

    # Print output 
    output = c(vs1,vs2, vs3)
    names(output) = c('Jimenez and Madsen (2003)', 'Rubey (1933)', 'van Rijn (1884) - discontinuities at 0.1 and 1 mm')

    print(output)
}


#R function to estimate the critical shear stress for non-cohesive sediments, following van Rijn (2007a), Equation 2 
tauc<-function(d){
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

    if(d_star<4){
        theta_c = 0.115*(d_star)^(-0.5)
    }else{
        if((d_star>=4.0)&(d_star<10.0)){
            theta_c = 0.14*d_star^(-0.64)
        }else{
            print(paste('d_star outside range of method, =', d_star))
            print('If we set d_star to 10, we get this result')
            d_star=10.
            theta_c = 0.14*d_star^(-0.64)
        }
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
    # Study the convergence of direct iteration for the colebrook equation

    k_sg = 10.0*d50
    #Re = max((sqrt(f_g/8.0_dp)*abs(vel))*k_sg/1.0e-06_dp, 10.0_dp)
    Re = abs(vel)*max(depth,20.*k_sg)/1.0e-06

    #Turbulent regime
    for(i in 1:10){
        f_g = 0.25/(log( k_sg/(3.71*4.0*max(depth,20*k_sg)) + 2.51/(Re*sqrt(f_g)),10 )  )**2
    }
    f_g

}



