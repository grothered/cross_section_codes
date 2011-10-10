PROGRAM driver_xsect
! Program to run the cross-sectional evolution

USE global_defs !Module with some fundamental constants
USE hydro_xsect !cross-sectional hydrodynamics
USE suspended_xsect ! cross-sectional suspended sediment distribution
USE bed_xsect ! cross-sectional sediment transport / bed evolution
USE util_various ! Random utilities that don't fit nearly elsewhere

IMPLICIT NONE
! Initialise variables -- there is no order to this, and some variables may be
! unused
INTEGER:: nos, i,j, l, u, n, layers, ii, ind(1), indlast, Q_loop, jj, &
          writfreq, jmax, iii, iii2
REAL(dp):: wslope, ar, Q, t, &
           ys,bed, water, waterM, dists, tau,ks,tbst, Qe, Qbed,qb_G, Qd, &
            f, aa, bb, cc, multa, bedold, Qelast, &
            qby,v, bedlast, hss, tss, ddd, hss2, handy, dqbeddx, &
            epsl,epsu, slopes, E, D, C, rmult,Area, inuc, NN,&
             Width,C0,waterlast, dt, slpmx,sllength, vel, &
            DT1,DT1_old, tau_g, f_g, Cbar,water_tmp, &
            taucrit_dep, hlim, mor, Arealast, taucrit_dep_ys, &
            ht,vlast, dst, taucrit, mu, & 
            erconst, lifttodrag, vegdrag, sconc, rho, ysold, &
            lincrem, wset, voidf, smax, rough_coef, man_nveg, veg_ht,rhos,&
            dsand, d50, g, kvis, lambdacon, alpha, &
            ysl,ysu,bedl, bedu, wdthx, TR, storer(9), tmp, tmp2, a_ref, &
            failure_slope, x_len_scale, sus_flux, sed_lag_scale, Clast, &
            lat_sus_flux, int_edif_f, int_edif_dfdy, zetamult
INTEGER::  remesh_freq, num_simulations, too_steep, morbl, morbu
REAL(dp):: discharges(1000), susconcs(1000)
LOGICAL::  compute_twice, susdist, sus2d, readin, geo, remesh, norm, vertical, & 
            tbston, normmov, Qbedon, susQbal, talmon,&
             variable_timestep, high_order_shear, high_order_bedload, &
            taucrit_slope_reduction, evolve_bed
CHARACTER(LEN=20):: friction_type, grain_friction_type, resus_type, &
                    bedload_type, sus_vert_prof, edify_model
NAMELIST /inputdata/ nos,writfreq,jmax, layers, hlim, mor, mu, &
                erconst,lifttodrag,sconc,rho,lincrem,wset, voidf, t, dt, &
                susdist,readin, geo, waterM, width, smax, rough_coef, &
                man_nveg, veg_ht, remesh, remesh_freq, rhos, dsand, d50,&
                g, kvis, norm, &
                vertical, lambdacon, alpha, tbston, normmov, sus2d, &
                Qbedon, susQbal, TR, talmon, variable_timestep, & 
                friction_type, num_simulations, &
                discharges, susconcs, high_order_shear, &
                high_order_bedload, grain_friction_type, &
                resus_type, bedload_type, sus_vert_prof, edify_model, &
                failure_slope, x_len_scale, taucrit_slope_reduction, &
                evolve_bed

ALLOCATABLE ys(:), bed(:), dists(:), tau(:),ks(:),tbst(:),& 
            qby(:), bedlast(:), hss(:), tss(:),  hss2(:), Qe(:),& 
            Qbed(:),qb_G(:), Qd(:),f(:),slopes(:), NN(:),C0(:),Qelast(:),&
            taucrit_dep(:,:), C(:),bedold(:), &
            taucrit_dep_ys(:) ,dst(:,:), taucrit(:,:), slpmx(:,:), &
            vegdrag(:), ysold(:), dqbeddx(:), sllength(:), vel(:), &
            tau_g(:),f_g(:), Cbar(:), a_ref(:), Clast(:), &
            lat_sus_flux(:), int_edif_f(:), int_edif_dfdy(:), &
            zetamult(:), too_steep(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ IN DATA

! Predefine these input vectors to help catch errors
discharges = -1.0_dp
susconcs = -1.0_dp

!INCLUDE 'input_data.modin'
!!Read the input parameters
!open(1001,file= 'input_data.modin')
!close(1001)
! Read from standard input . To do this, type  ./model < input_file ( and then
! you can pipe to an output by adding > outfile.log)
READ(*,nml=inputdata)

PRINT inputdata

! Test that the inputs of discharges, susconcs seem okay
IF(size(discharges) < num_simulations) THEN
    print*, 'Error: the discharges vector can be at most of length ', &
             size(discharges)
    stop
END IF 

IF((minval(discharges(1:num_simulations))<0.0_dp).OR. & 
   (minval(susconcs(1:num_simulations))<0.0_dp)) THEN
    print*, 'Error: Discharges or susconcs are specified incorrectly. &
            Possibly there are not enough values'
    stop
ELSE IF((maxval(discharges(num_simulations+1:1000))>0.0_dp).OR. & 
   (minval(susconcs(num_simulations+1:1000))>0.0_dp)) THEN
    print*, 'Error: Discharges or susconcs are specified incorrectly. &
            Possibly there are too many values'
    stop
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open files which will later be written to

OPEN(1,file="taus")
OPEN(2,file="bed")
OPEN(3,file="ys")
OPEN(4,file="taug")
OPEN(5,file="water")
OPEN(7,file="qe")
OPEN(8,file="Qbed")
OPEN(9,file="Cbed")
OPEN(10,file="vel")
OPEN(11,file="qby")
OPEN(12,file="a_ref")
OPEN(13,file="timestepping_stats")
OPEN(14,file='Fl')
! Allocate memory for arrays
ALLOCATE(ys(nos),bed(nos),dists(nos),tau(nos),ks(nos),tbst(nos),& 
         qby(0:nos), bedlast(nos),hss(nos),tss(nos),hss2(nos),Qe(nos),&
         Qd(nos),f(nos),slopes(nos),Qelast(nos),& 
         NN(nos), C0(nos),taucrit_dep(nos, layers), C(nos), &
         taucrit_dep_ys(nos),dst(nos,0:(layers+1)), &
         taucrit(nos, 0:layers), slpmx(nos,0:(layers+1)), vegdrag(nos),&
         ysold(nos) , Qbed(nos),qb_G(0:nos+1), dqbeddx(nos),sllength(nos), &
         vel(nos), tau_g(nos), f_g(nos), Cbar(nos), bedold(nos), a_ref(nos),&
         Clast(nos), lat_sus_flux(nos+1), int_edif_f(nos+1), &
         int_edif_dfdy(nos+1) , zetamult(0:nos+1), too_steep(nos)) 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!LOOP OVER DIFFERENT VALUES OF DISCHARGE

DO Q_loop= 1, num_simulations!15 

    ! Set suspended sediment concentration for this discharge
    sconc= susconcs(Q_loop) 

    !Create the channel geometry
    If(readin) THEN
        !Cross-channel distance
        OPEN(77,file="lnths",status="old")
        READ(77,*) ys
        CLOSE(77)
        
        ! Bed elevation
        OPEN(78,file="hes",status="old")
        READ(78,*) bed
        CLOSE(78)
        
        ! Bed layer cross-channel distance
        open(79,file="taucrit_dep_lnths")
        read(79,*) taucrit_dep_ys
        close(79)
        
        ! Bed layer elevation
        open(80,file="taucrit_deps")
        read(80,*) taucrit_dep(:,1)
        ! Predefine the deeper bed layers
        do i=1, layers
        taucrit_dep(:,i)= taucrit_dep(:, 1) - i*10._dp
        end do
        close(80)
    ELSE
        !Define the channel geometry. Note that for most of my boundary
        !conditions to be valid, we need the edges of the computational domain
        !to not to ever be sumberged. 
        DO i =1,nos
        
            ! Cross-channel distance, 'y' coordinate
            ys(i)=Width/(nos*1._dp-1._dp)*(i*1._dp-1._dp) 
            
            ! Bed elevation
            !bed(i)= 1.11_dp -8.01_dp*(abs(-1._dp+1._dp*((ys(i)-Width/2._dp))**6/(Width/2._dp)**6)) !-3.01 !+0.00002*(ys(i)-ys(1))*(ys(i)-nos*ys(1))
            !bed(i)= -4.1_dp -0.01_dp*(abs(-1._dp+1._dp*((ys(i)-Width/2._dp))**2/(Width/2._dp)**2)) !-3.01 !+0.00002*(ys(i)-ys(1))*(ys(i)-nos*ys(1))
            !tmp = 40.0 !sqrt(discharges(Q_loop)/2.5_dp)
            tmp = sqrt(discharges(Q_loop)/2.5_dp)
            bed(i) =min( max(abs(ys(i)-Width/2._dp)-1.5_dp*tmp,-tmp ),&
                         0._dp)
            
            !!The SERM channel
            !!IF(handy<=0.75) bed(i)=0.
            !!IF((handy>0.75).AND.(handy<0.9)) bed(i)= (handy-0.75)
            !!IF((handy>=0.9).AND.(handy<3.13)) bed(i)= 0.15
            !!IF(handy>=3.13) bed(i)= (handy-3.13) +0.15
            
            ! Bed layers    
            taucrit_dep(i,:)=-1000._dp
            taucrit_dep_ys(i)= ys(i)

        END DO
        
        
        
    END IF !.readin.

    !!!!!!!!!!!!!!!!!!!!!!
    !INITIALISE VARIABLES
    !
    !!!!!!!!!!!!!!!!!!!!!!

    t=0._dp !Time
    qby=0._dp !A useful variable for storing stuff 
    l=1 !variable for "lower" wet point of cross section
    u=nos !variable for "upper" wet point of cross section
    water=waterM ! Water elevation
    ! Redefine l and u
    call wet(l,u,nos,water, bed)
    tau= 0._dp*tau !Set all the shears to zero prior to update
    E=0._dp ! Cross-sectionally integrated resuspension 
    D=0._dp ! Cross-sectionally integrated deposition (from suspension)
    C=0.0_dp !sconc ! Near bed sediment concentration
    Cbar = 0.0_dp !sconc ! Depth averaged suspended sediment concentration 
    f= 0.02_dp !Friction factor
    f_g= 0.02_dp ! Grain friction factor
    vel = 0.0_dp ! Velocity
    ! rmult = Cross-sectionally averaged friction factor: 
    ! Note that in general, rmult ~= f/depth, because when
    ! interfacing with the 1D St-Venant solver, this is more stable.
    rmult=sum(f)/(8._dp*9.8_dp*nos) 
    vegdrag=0._dp  ! Drag coefficient for vegetation
    ! inuc = Multiplicative factor to account for cross-sectional non-uniformity
    ! in the convective intertial terms in the St Venant equation solver 
    inuc=0._dp
    NN=tau*0._dp ! Defunct constant
    vlast=0._dp
    hss=0._dp
    hss2=0._dp
    Area=0.0000_dp ! Cross-sectional area
    dst(:, layers) = 2000._dp ! Distance from the bed surface to different bed layers
    slpmx(:, layers+1) = 99999999._dp
    dqbeddx=0._dp
    bedlast = bed - 0.001_dp !
    bedold = bed - 0.001_dp !
    DT1 = dt ! An adaptive timestep
    a_ref=0.0_dp !Reference level of vanrijn
    Qe = 0.0_dp
    Qelast=0.0_dp
    int_edif_f=0.0_dp
    int_edif_dfdy=0.0_dp
    sed_lag_scale=1.0_dp
    zetamult=1.0_dp

    !!!Calculate area
    IF(l>0) THEN
        DO i= 1,nos-1
            Area= Area+ max(water-0.5_dp*(bed(i)+bed(i+1)), 0._dp)&
                 *(ys(i+1)-ys(i))
        END DO
    ELSE
        print*, "Totally Dry at start"
        stop
    END IF

    ! Calculate slopes -- note how it is a distance weighted average of the
    ! 1eft and right slopes
    slopes(2:nos-1)= ((bed(3:nos)-bed(2:nos-1))/(ys(3:nos)-ys(2:nos-1))*(ys(2:nos-1)- & 
        ys(1:nos-2)) + (bed(2:nos-1)-bed(1:nos-2))/(ys(2:nos-1)-ys(1:nos-2)) & 
        *(ys(3:nos)-ys(2:nos-1))) /(ys(3:nos)-ys(1:nos-2))

    slopes(1)= (bed(2)-bed(1))/(ys(2)-ys(1))
    slopes(nos)= (bed(nos)-bed(nos-1))/(ys(nos)-ys(nos-1))


    !!!!!!!!!!!!!!
    !THE MAIN LOOP
    !!!!!!!!!!!!!!
    DO j=1, jmax 

        ! Print out to console
        22222 IF( mod(j-1,writfreq).eq.0 ) THEN 
                PRINT*, j, l,u, Q/Area, t, ((Q/Area)*abs(Q/Area)*rmult),&
                     DT1, ys(u)-ys(l)+wdthx, maxval(C), C(nos/2),&
                    rmult*(Area)/(ys(u)-ys(l)+wdthx), f(nos/2), sed_lag_scale
              END IF


        ! Update the critical shear layers to account for any slope effects
        DO i = 1, nos
            DO jj= 0, layers

                multa=1._dp
                
                IF(taucrit_slope_reduction.eqv..TRUE.) THEN 
                    ! Compute slope-related reduction in critical shear stress
                    IF((j.eq.1).and.(i.eq.1).and.(jj==0)) print*, 'WARNING: Critical shear on a slope is reduced'
                    
                    aa= mu**2*(mu*lifttodrag-1._dp)/(mu*lifttodrag+1._dp)
                    bb= -2._dp*mu**2*lifttodrag*cos(atan(slopes(i)))*mu/(1._dp+mu*lifttodrag) 
                    cc= mu**2*cos(atan(slopes(i)))**2 - sin(atan(slopes(i)))**2
                    !multa*(critical shear on a flat bed) = critical shear on a slope.
                    multa= (-bb - sqrt(bb**2-4._dp*aa*cc))/ (2._dp*aa)  

                END IF

                taucrit(i,jj) = erconst*(1._dp+ jj*1._dp)*max(multa,1.0e-01_dp)

                !if(bed(i)>-1.0_dp) taucrit(i,jj) = taucrit(i,jj)*2.0_dp

                IF( isnan(taucrit(i,jj))) THEN
                    PRINT*, "taucrit(", i,",", jj, ") is nan"
                    STOP
                END IF

            END DO 
        END DO

        ! taucrit_dep = the value of taucrit in sediment buried at a
        ! certain depth
        DO  i= 1, nos
            DO n= 1, layers
                taucrit_dep(i,n)= max(min(taucrit_dep(i,n), bed(i)),bed(i)-lincrem*n) !So sediments buried at a certain depth will have their critical shear increase. 
            END DO !N 
        END DO


        !dst = the distance of any critical shear layer to the surface.
        dst(:,0)= 0._dp
        DO i= 1, nos
            DO jj= 1, layers
                dst(i, jj)=max((bed(i)-taucrit_dep(i,jj) ), 0._dp)
                ! print*, hs(i), taucrit_dep(i, 1), dst(i, 1)
            END DO
        END DO


        ! Here either the geotech routine is used, or just the slope is calculated
        IF(geo) THEN
            call geotech(ys, bed,slopes , taucrit_dep_ys,taucrit_dep, dst, nos,nos ,mu, layers, slpmx)
        ELSE
            ! Slope term -- note how it is a distance weighted average of the 1eft and right slopes
            slopes(2:nos-1)= ((bed(3:nos)-bed(2:nos-1))/(ys(3:nos)-ys(2:nos-1))*(ys(2:nos-1)- & 
            ys(1:nos-2)) + (bed(2:nos-1)-bed(1:nos-2))/(ys(2:nos-1)-ys(1:nos-2)) & 
            *(ys(3:nos)-ys(2:nos-1))) /(ys(3:nos)-ys(1:nos-2))  
            
            slopes(1)= (bed(2)-bed(1))/(ys(2)-ys(1))
            slopes(nos)= (bed(nos)-bed(nos-1))/(ys(nos)-ys(nos-1))   
        END IF

        ! Note where the slope is overly large, so we can prevent deposition
        ! there
        too_steep=1
        ! Find where the slope > failure_slope, and prevent deposition on the
        ! upper part of that slope.
        !DO i=1,nos
        !    IF((i==nos).or.(bed(i)>bed(max(i-1,1)))) THEN
        !        IF( (bed(i)-bed(i-1)) >= failure_slope*(ys(i)-ys(i-1))) THEN
        !            too_steep(i)=0
        !        END IF
        !    ELSEIF((i==1).or.(bed(i)>bed(min(i+1,nos)))) THEN
        !        IF( (bed(i)-bed(i+1)) >= failure_slope*(ys(i+1)-ys(i))) THEN
        !            too_steep(i)=0
        !        END IF
        !    END IF
        !END DO

        ! Water elevation (free surface)
        waterlast=water
        water= waterM + (TR/2._dp+.0*sin(2._dp*pi*t/(12.4_dp*3600._dp*100._dp))+0.0_dp)&
                        *sin(2._dp*pi*t/(12.4_dp*3600._dp))  !1.+ 1.*sin(2.*pi*j/500.) 

        ! Find wetted part of section
        call wet(l,u,nos,water, bed) 
        
        ! 'Extra' wetted width associated with the corners of the domain -- do
        ! we really need/ want this?
        wdthx = 0.0
        IF(l>0) THEN
            IF (u<nos) wdthx = wdthx+ (water-bed(u))/(bed(u+1)-bed(u))*(ys(u+1)-ys(u))
            IF (l>1) wdthx = wdthx+ (water-bed(l))/(bed(l-1)-bed(l))*(ys(l)-ys(l-1))
        END IF

        ! Cross sectional area
        Arealast=Area
        Area = compute_area(nos, water,bed,ys,l,u)

        IF ( (Area<0.0002).OR.((l>0).and.(maxval(water-bed)<.01)) )THEN 
            !If we let the area go to zero, then the continuity model will allow
            !very very high velocities (since Q/A can be very large). This is
            !unrealistic because in practise the flat water slope assumption
            !would be very very wrong here. So, we limit the area to a box 20cm
            !by 10cm. 
           print*, "Totally Dry, or maybe just effectively dry", j, maxval(C), &
                    water, Area, maxval(bed), minval(bed)!, bed
           goto 22222 !Back to the start
        END IF

        ! Compute discharge
        IF(TR==0._dp) THEN
            Q=discharges(Q_loop)!*(1.0_dp+0.1_dp*sin(2._dp*pi*t/(12.4_dp*3600._dp))) !/50._dp
        ELSE
            Q= (discharges(Q_loop)+j*0._dp/500._dp)*(ys(u)-ys(l)+wdthx)*abs((water-waterlast))/dt*10._dp 
        END IF

        ! Sanity checks
        IF(isnan(Q)) THEN
            print*, "Q is nan", Q, Area, Arealast
            stop
        END IF
        IF((Q/Area)>5.) THEN
            print*, "Q/A>5", Q/Area, Q, Area, Arealast, j,  dt, water, waterlast, width, t, dt!, l, u, water-bed(l:u)
            !stop     
        END IF



        !!!!!!!!!!!!!!!
        !CALCULATE SHEAR AND UPDATE THE CROSS SECTION
        !!!!!!!!!!!!!!
        qby=0._dp
        tau=0._dp
        Qelast=Qe
        Qe=0._dp
        Qbed=0._dp
        !tss=bed
        ! Determine the y and elevation values of the points just on the dry
        ! edge of the cross-section -- needed for the shear and bed solver 
        IF(l>1) THEN
            ysl=ys(l-1)
            bedl=bed(l-1)
        ELSE
            ysl=ys(l)- 0.001_dp 
            bedl=water 
        END IF
        IF(u<nos) THEN
            ysu=ys(u+1)
            bedu=bed(u+1)
        ELSE
            ysu=ys(u)+0.001_dp
            bedu=water 
        END IF

        DO iii=1, 1 !This can be used to try fancy time stepping techniques
            !IF(iii==1) bedlast=bed
            !f(l)=f(l+1)
            !f(u)=f(u-1)
            !slopes(l)=(bed(l+1)-bed(l))/(ys(l+1)-ys(l))
            !slopes(u)=(bed(u-1)-bed(u))/(ys(u-1)-ys(u))

            !!With the slopes at the edge: How about we extrapolate the wetted boundary
            !location, and use this to calculate the slope at the numerical boundaries

            !slopes(l) = (bed(l+1) - water )/&
            !(ys(l+1) - (ys(l) - (water-bed(l))/(bedl-bed(l))*(ys(l)-ysl) ) )

            !slopes(u) = (bed(u-1) - water )/&
            !(ys(u-1) - (ys(u) + (water-bed(u))/(bedu-bed(u))*(ysu-ys(u)) ) )
            !slopes(l)= (slopes(l)*sqrt(1._dp+slopes(l)**2) +&
            !(bedl-bed(l))/(ysl-ys(l))*sqrt(1._dp+( (bedl-bed(l))/(ysl-ys(l)))**2) )/  &
            !(sqrt(1._dp+( (bedl-bed(l))/(ysl-ys(l)))**2) + sqrt(1._dp+slopes(l)**2) )
            !
            !slopes(u)= (slopes(u)*sqrt(1._dp+slopes(u)**2) +&
            !(bedu-bed(u))/(ysu-ys(u))*sqrt(1._dp+( (bedu-bed(u))/(ysu-ys(u)))**2) )/  &
            !(sqrt(1._dp+( (bedu-bed(u))/(ysu-ys(u)))**2) + sqrt(1._dp+slopes(u)**2) )
       
            ! During the first time step, we compute friction and shear twice,
            ! because each depends on the other. Indicate the need to compute twice
            ! with compute_twice=.TRUE. 
            IF(j==1) compute_twice=.TRUE.            
            
            ! CALCULATE FRICTION on bed 'i' with vel 'i-1'
            1987 call calc_friction(friction_type, grain_friction_type, rough_coef, water, u-l+1,&
                                 bed(l:u), vel(l:u), man_nveg,d50,veg_ht, rhos, rho, g,&
                                 f(l:u), vegdrag(l:u),f_g(l:u), dsand, j, a_ref(l:u)) 

            ! CALCULATE BED SHEAR on bed 'i'
            call calc_shear(u-l+1,DT1,water,Q,bed(l:u),ys(l:u),Area, &
                            water-Area/(ys(u)-ys(l)+wdthx),f(l:u),&
                            rmult,inuc, tau(l:u),& 
                            NN(l:u),j,slopes(l:u), hlim, &
                            u-l+1, vegdrag(l:u), rho & 
                            ,rhos, voidf, d50, g, kvis, norm, vertical, lambdacon, tbston &
                            ,ysl,ysu,bedl,bedu, high_order_shear) 

            ! Calculate depth-averaged velocity
            vel = 0._dp 
            vel(l:u)=sqrt(abs(tau(l:u))/rho*8._dp/f(l:u))*sign(1._dp+0._dp*tau(l:u), tau(l:u))
         
            ! Calculate grain shear stress 
            ! NOTE that vanrijn writes that taug = 0.5*f_g*rho*U^2 -- however, his
            ! data in table2 of the paper are better predicted using the
            ! 'normal' formula, tau = rho f/8 U^2 --- I think the paper must
            ! have a typo. See the function test_vrijn_bed() in wset.R
            tau_g = 0._dp
            !tau_g(l:u) = 0.5_dp*rho*vel(l:u)**2*(f_g(l:u))*sign(1._dp+0._dp*tau(l:u), tau(l:u))
            !Following Abdel-Fattah et al 2004
            tau_g(l:u) = rho*vel(l:u)**2*(f_g(l:u)/8._dp)*sign(1._dp+0._dp*tau(l:u), tau(l:u))
          
            IF((j==1).AND.(compute_twice)) THEN
                ! Compute friction and shear again on the very first
                ! time-step. 
                compute_twice=.FALSE.
                GOTO 1987
            END IF
      
            ! CALCULATE RATES OF RESUSPENSION AND BEDLOAD TRANSPORT on bed 'i'
            ! with shear 'i'
            call calc_resus_bedload(u-l+1,DT1,water,Q,bed(l:u),ys(l:u),Area,&
                                    f(l:u),qby((l-1):u),E,&
                                    C(l:u),wset, 2,tau(l:u), tau_g(l:u),& 
                                    vel(l:u), j,slopes(l:u), hlim, mor, taucrit_dep(l:u,1:layers),&
                                     layers, taucrit_dep_ys(l:u) & 
                                    ,u-l+1, taucrit(l:u, 0:layers) , rho, Qe(l:u) & 
                                    , Qbed(l:u),qb_G((l-1):(u+1)), rhos, voidf, dsand, d50, g, kvis, norm, alpha, &
                                    Qbedon,talmon, ysl,ysu,bedl,bedu, resus_type, bedload_type, a_ref(l:u)) 
        
            !! UPDATE TIME
            IF(iii.eq.1) t=t+DT1

            !Qe(l+1:u-1) = 0.5_dp*Qe(l+1:u-1) + 0.25_dp*(Qe(l:u-2) + Qe(l+2:u))
            !Qe(l) = 0.75_dp*Qe(l) + 0.25_dp*Qe(l+1)
            !Qe(u) = 0.75_dp*Qe(u) + 0.25_dp*Qe(u-1)
            ! Calculate total sediment flux at time = t, 
            ! We will use this to compute the x derivative terms
            ! in dynamic_sus_dist and update_bed
            ! e.g. dQbed/dx = (Qbed - sed_lag_scale*Qbed)/x_len_scale
            ! e.g. dCbar/dx = (Cbar - sed_lag_scale*Cbar)/x_len_scale
            !sus_flux = sum( & 
            !        (Qbed(l:u)+(Cbar(l:u)/rhos)*abs(vel(l:u))*max(water-bed(l:u),0._dp) )*& ! Total load
            !        ( ( (/ ys(l+1:u), ysu /) - (/ ysl, ys(l:u-1) /) )*0.5_dp) &  ! dy
            !              )
            !
            !IF(sus_flux > 1.0e-12_dp) THEN
            !    sed_lag_scale = 1.0_dp*((sconc*Q)/sus_flux)
            !    ! Prevent very high or low values
            !    sed_lag_scale = min(max(sed_lag_scale,0.666_dp),1.5_dp) 
            !    !IF(mod(j,1000).eq.1) PRINT*, 'sed_lag_scale = ', sed_lag_scale
            !    !print*, sed_lag_scale                        
            !ELSE
            !    sed_lag_scale = 1.0_dp
            !END IF
            !IF(mod(j,1000).eq.1) print*, 'sus_flux is:', sus_flux, ' desired flux is:', sconc*Q


            !! CALCULATE THE CROSS-SECTIONAL SUSPENDED LOAD DISTRIBUTION
            IF(susdist) THEN
                !Record old value of C for file output
                Clast=C  
                lat_sus_flux = 0.0_dp ! Preset to zero

                ! Evolve the suspended sediment concentration one timestep
                ! (from i-1 to i)
                call dynamic_sus_dist(u-l+1, DT1, ys(l:u), bed(l:u), water, waterlast, Q, tau(l:u), vel(l:u), wset, & 
                                        0.5_dp*(Qe(l:u)+Qelast(l:u)), lambdacon, rho,rhos, g, d50, bedl,bedu, ysl, ysu, C(l:u),&
                                        Cbar(l:u), Qbed(l:u), sed_lag_scale, j, a_ref(l:u), sus_vert_prof,&
                                        edify_model, x_len_scale, sconc, lat_sus_flux(l:u+1), bedlast(l:u), int_edif_f(l:u+1), &
                                        int_edif_dfdy(l:u+1), zetamult((l-1):(u+1)), too_steep(l:u))

                ! Set C in dry parts of the channel to zero
                ! This is not done earlier, because we needed to store the old
                ! values of C for the dynamic_sus_dist routine
                IF(l>1) THEN
                    Cbar(1:(l-1)) = 0._dp
                    C(1:(l-1)) = 0._dp
                    ! Set the integral terms in the lateral diffusive flux of
                    ! suspended sediment to zero in dry regions.
                    int_edif_f(1:max(l-1,1))=0.0_dp
                    int_edif_dfdy(1:max(l-1,1))=0.0_dp
                    
                END IF
                
                IF(u<nos) THEN
                    Cbar((u+1):nos) = 0._dp
                    C((u+1):nos) = 0._dp
                    
                    int_edif_f(min(u+2,nos):nos)=0.0_dp
                    int_edif_dfdy(min(u+2,nos):nos)=0.0_dp
                END IF
        

            ELSE
                ! Constant suspended sediment concentration
                Clast = C
                C(l:u) = sconc
                Cbar(l:u) = sconc
            END IF
         
            !! UPDATE THE BED
 
            ! Calculate dqbed/dx ~= (Qbed - sed_lag_scale*Qbed)/x_len_scale
            dqbeddx(l:u) = Qbed(l:u)*(1.0_dp-sed_lag_scale)/x_len_scale 
            bedlast= bed ! Record the bed prior to updating
         
            IF(evolve_bed) THEN 
                ! Update the bed from i to i+1, using the rates of erosion and
                ! deposition calculated from bed i. 
                call update_bed(u-l+1,DT1,water,Q,bed(l:u),ys(l:u),Area, &
                                qby((l-1):u),E,D, &
                                C(l:u),2,tau(l:u),tau_g(l:u),& 
                                j,slopes(l:u), hlim, mor, taucrit_dep(l:u,1:layers), &
                                layers, taucrit_dep_ys(l:u) & 
                                ,u-l+1, taucrit(l:u, 0:layers) , rho, &
                                Qe(l:u), Qbed(l:u),qb_G((l-1):(u+1)), wset, dqbeddx(l:u),&
                                rhos, voidf, d50, g, norm, &
                                Qbedon, normmov, sus2d, ysl, ysu, bedl,bedu, iii, bedlast(l:u), &
                                talmon, high_order_bedload, too_steep) 

                ! Correct the banks. In the case that we allow bedload at l-1/2 and
                ! u+1/2, this is very important to ensure mass conservation, because
                ! if there is a downslope bedload flux from l-1/2, or from u+1/2,
                ! then it must come from the dry part of the channel 
                !IF(.FALSE.) THEN
                !    ! Use this case when bedload occurs and l-1/2, u+1/2
                !    IF(l>1) THEN
                !        IF(bed(l-1)>bedl) bed(l-1)=bedl
                !    END IF
                !    IF(u<nos) THEN
                !        IF(bed(u+1)>bedu) bed(u+1)=bedu
                !    END IF
                !END IF

                !IF(.TRUE.) THEN
                !    ! A version of the Delft bank erosion model. 
                !    ! First check that there is no leakage of bedl, bedu in the
                !    ! bed solver (possibly could happen due to matrix round off or
                !    ! coding error).
                !    IF((l>1).and.(u<nos)) THEN
                !        IF((abs(bedl - bedlast(l-1))>0.0e-8_dp).OR.& 
                !           (abs(bedu-bedlast(u+1))>0.0e-08_dp)) THEN
                !            print*, 'ERROR -- there is still erosion of dry points in &
                !                the matrix solution of the bed solver'
                !            print*, abs(bedl - bedlast(l-1)), abs(bedu-bedlast(u-1)), l, u
                !            stop
                !        END IF
                !    END IF
                !    ! If erosion is occuring at the channel margins,
                !    ! then assign it to the neighbouring dry bed point
                !    IF((bed(l)<bedlast(l)).AND.(l>1)) THEN
                !        !IF( abs(bed(l) - bed(l-1))/(ys(l)-ys(l-1))>1.0_dp) THEN
                !        !IF( abs(tau_g(l))>taucrit(l,0)) THEN
                !            bed(l-1) = bed(l-1) - (bedlast(l) - bed(l))
                !            bed(l) = bedlast(l)
                !        !END IF
                !    END IF
                !    IF((bed(u)<bedlast(u)).AND.(u<nos)) THEN
                !        !IF( abs(bed(u+1) - bed(u))/(ys(u+1)-ys(u))>1.0_dp) THEN
                !        !IF( abs(tau_g(u))>taucrit(u,0)) THEN
                !            bed(u+1) = bed(u+1) - (bedlast(u) - bed(u))
                !            bed(u) = bedlast(u)
                !        !END IF
                !    END IF
                !END IF
               
     
                ! BASIC LIMITING OF THE CHANNEL SLOPE -- to circumvent the numerically
                ! difficult problem of allowing infinite banks otherwise
                !IF(mod(j,1)==0)
                !call basic_slope_limit(nos,ys,bed,failure_slope, remesh, 1.0_dp)
                !call basic_jump_limit(nos,ys,bed,0.5_dp, remesh, 1.0_dp)
                call critical_slope_wasting(DT1, nos,ys,bed,failure_slope, 1.0e-06_dp)
                !call critical_bedjump_wasting(DT1, nos,ys,bed,2.0_dp, 1.0e-05_dp)
                !Update Cbar to reflect changes in the bed.
                !DO i=1,nos
                !    IF((water>bedlast(i)).and.(water>bed(i))) THEN
                !        Cbar(i) = Cbar(i)*(water-bedlast(i))/(water-bed(i))
                !    END IF
                !END DO

            END IF

            !! WRITE OUTPUTS -- notice that these are all supposed to be at the
            !same time level -- e.g. tau is calculated using bedlast, so is
            !Clast, Qbed, Qe, etc. 
            IF((mod(j-1,writfreq).EQ.0).AND.(iii.eq.1)) THEN!.or.(j>15250)) THEN 
                !Qd = Clast*wset/rhos !0.5_dp*(Clast +C)*wset/rhos
                print*, 'bed change:', maxval(abs(bed-bedlast)), maxloc(abs(bed-bedlast))
                !print*, 'Resus - dep:', maxval(abs(Qe(l:u) - Qd(l:u)))*mor*DT1, &
                !                        maxloc(abs(Qe(l:u) - Qd(l:u)))
                !print*, 'Qe(mid) = ', Qe(nos/2)

                write(1,*) tau 
                write(2,*) bedlast !Same bed as when tau was calculated
                write(3,*) ys !critical shear
                write(4,*) tau_g !taucrit_dep!water, Q/ar
                write(5,*) water
                write(7,*) Qe
                write(8,*) Qbed
                write(9,*) Clast ! Same C as when bed = bedlast and when tau was calculated
                write(10,*) vel
                write(11,*) qby
                write(12,*) a_ref
                write(13,*) j, l,u, Q/Area, t-DT1, ((Q/Area)*abs(Q/Area)*rmult),&
                            DT1, ys(u)-ys(l)+wdthx, maxval(C),C(nos/2), &
                            rmult*(Area)/(ys(u)-ys(l)+wdthx), f(nos/2)
                write(14,*) lat_sus_flux
                !write(12,*) taucrit_dep_ys

                ! Check for convergence by comparing 'bed' with 'bedold' (=
                ! value of 'bed' at 'writfreq' time steps ago). Note that we
                ! also check against bedlast, which can sometimes be more
                ! different from bed, if oscillations are occurring -- which is
                ! good to catch.    
                tmp =max(maxval(abs(bedold-bed)), maxval(abs(bed-bedlast)))
                IF((evolve_bed).and.(tmp/(DT1*writfreq)<1.0e-12_dp)) THEN
                    print*, 'Converged due to small bed changes'
                    goto 373 !Converged: Go to the end of this loop
                END IF
                bedold=bed
            END IF

            ! DETERMINE THE TIMESTEP -- implicit timestepping will only work
            ! if there are no bed layers, because otherwise DT
            ! has already been used this timestep.
            IF(variable_timestep) THEN
                print*, 'ERROR: VARIABLE TIMESTEP NOT SUPPORTED AT PRESENT'
                stop                
                ! These methods change DT1 during the evolution of the
                ! cross-section
                
                !IF(j.eq.1) print*, ' Warning: Variable timestep is ONLY valid for &
                !        STEADY UNIFORM EQUILIBRIUM computations'
                !DT1 = min(DT, minval(wset/(water-bed(l:u))))
                
                !IF(mod(j,10).eq.1) THEN
                !    tmp = min(max(maxval(abs(wset*C/rhos- Qe)), maxval(abs(qby(l-1:u)))), &
                !              maxval(abs(bed(l+1:u-1) - bedlast(l+1:u-1))) )
                !    tmp2 = minval(ys(2:nos) - ys(1:nos-1)) 
                !    DT1 = min(max(5.0e-03_dp*tmp2/max(tmp,1.0e-020_dp), 1.0e-01_dp, 0.9_dp*DT1), 100.0_dp*3600.0_dp, 1.1_dp*DT1)
                !END IF
            ELSE
                DT1 = DT
            END IF 

        END DO



        !! REMESHING. If the points are becoming too spaced apart, we might have to do this
        IF(remesh) THEN
            print*, 'ERROR: REMESHING PRESENTLY NOT SUPPORTED'
            stop
            IF(mod(j, remesh_freq).EQ.0) THEN
                ysold=ys !Predefine
                !call refit(ys, bed, nos) !Remesh so things are more even.

                ! Find where heights has changed
                !call active_zone(nos-2, bed(2:nos-1), bedold(2:nos-1), morbl, morbu, 6)
                !print*, morbl+morbu, morbl, morbu
                ! Redefine ys
                !call reset_ys(nos-2, ys(2:nos-1), morbl, morbu, 0.1_dp, 5.0_dp)
                ! Also remesh old timestep records of the bed
                call interp3(ysold, bedlast,ys,nos)
                call interp3(ysold, bedold, ys, nos)
                ! Also remesh the suspended sediment concentration.
                call interp3(ysold, C, ys, nos)
                call interp3(ysold, Cbar, ys, nos)
                !Now fix up the layers
                DO jj=1, layers
                    call interp(ysold,  taucrit_dep(:,jj),ys, nos) !Beware the risk that this could cause 'leakage'
                    taucrit_dep(:, jj)= min(bed, taucrit_dep(:,jj)) 
                    !So if we interpolate the critical shear layer, it should of course not be above the bed level.
                END DO

            END IF 
        END IF !REMESH

    END DO ! End of timestepping loop 
        
        373 print*, "run finished, Q_loop=", Q_loop, "j =", j
        !stop
END DO ! End of loop for different Q_loops

close(1)
close(2)
close(3)
close(4)
close(5)
close(7)
close(8)
close(9)
close(10)
close(11)
close(12)
close(13)
close(14)

END PROGRAM      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


