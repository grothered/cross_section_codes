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
           ys,bed, water, waterM, dists, tau,ks,tbst, Qe, Qbed, Qd, &
            f, aa, bb, cc, multa, bedold, &
            recrd,v, bedlast, hss, tss, ddd, hss2, handy, dqbeddx, &
            epsl,epsu, slopes, E, D, C, rmult,Area, inuc, NN,&
             Width,C0,waterlast, dt, slpmx,sllength, vel, &
            integrated_load_flux, DT1,DT1_old, tau_g, f_g, Cbar, &
            taucrit_dep, hlim, mor, Arealast, taucrit_dep_ys, &
            ht,vlast, dst, taucrit, mu, & 
            tauinc, erconst, lifttodrag, vegdrag, sconc, rho, ysold, &
            lincrem, wset, voidf, smax, rough_coef, man_nveg, veg_ht,rhos,&
            dsand, d50, g, kvis, lambdacon, alpha, &
            ysl,ysu,bedl, bedu, wdthx, TR, storer(9), tmp, tmp2, a_ref
INTEGER::  remesh_freq, no_discharges
REAL(dp):: discharges(1000), susconcs(1000)
LOGICAL::  flag, susdist, sus2d, readin, geo, remesh, norm, vertical, & 
            tbston, normmov, Qbedon, susQbal, talmon,&
             variable_timestep, high_order_shear, high_order_bedload, &
            high_order_Cflux
CHARACTER(LEN=20):: friction_type, grain_friction_type, resus_type, &
                    bedload_type, sus_vert_prof
NAMELIST /inputdata/ nos,writfreq,jmax, layers, hlim, mor, mu, tauinc,&
                erconst,lifttodrag,sconc,rho,lincrem,wset, voidf, t, dt, &
                susdist,readin, geo, waterM, width, smax, rough_coef, &
                man_nveg, veg_ht, remesh, remesh_freq, rhos, dsand, d50,&
                g, kvis, norm, &
                vertical, lambdacon, alpha, tbston, normmov, sus2d, &
                Qbedon, susQbal, TR, talmon, variable_timestep, & 
                integrated_load_flux, friction_type, no_discharges, &
                discharges, susconcs, high_order_shear, &
                high_order_bedload, high_order_Cflux, grain_friction_type, &
                resus_type, bedload_type, sus_vert_prof

ALLOCATABLE ys(:), bed(:), dists(:), tau(:),ks(:),tbst(:),& 
            recrd(:), bedlast(:), hss(:), tss(:),  hss2(:), Qe(:),& 
            Qbed(:), Qd(:),f(:),slopes(:), NN(:),C0(:),&
            taucrit_dep(:,:), C(:),bedold(:), &
            taucrit_dep_ys(:) ,dst(:,:), taucrit(:,:), slpmx(:,:), &
            vegdrag(:), ysold(:), dqbeddx(:), sllength(:), vel(:), &
            tau_g(:),f_g(:), Cbar(:), a_ref(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ IN DATA

! Predefine these input vectors to help catch errors
discharges = -1.0_dp
susconcs = -1.0_dp

!!Read the input parameters
!open(1001,file= 'input_data.modin')
!close(1001)
! Read from standard input . To do this, type  ./model < input_file ( and then
! you can pipe to an output by adding > outfile.log)
READ(*,nml=inputdata)

PRINT inputdata

! Test that the inputs of discharges, susconcs seem okay
IF(size(discharges) < no_discharges) THEN
    print*, 'Error: the discharges vector can be at most of length ', &
             size(discharges)
    stop
END IF 

IF((minval(discharges(1:no_discharges))<0.0_dp).OR. & 
   (minval(susconcs(1:no_discharges))<0.0_dp)) THEN
    print*, 'Error: Discharges or susconcs are specified incorrectly. &
            Possibly there are not enough values'
    stop
ELSE IF((maxval(discharges(no_discharges+1:1000))>0.0_dp).OR. & 
   (minval(susconcs(no_discharges+1:1000))>0.0_dp)) THEN
    print*, 'Error: Discharges or susconcs are specified incorrectly. &
            Possibly there are too many values'
    stop
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open files which will later be written to

OPEN(1,file="taus")
OPEN(2,file="bed")
OPEN(3,file="misc")
OPEN(4,file="misc2")
OPEN(5,file="usef")
OPEN(7,file="q")
OPEN(8,file="Area")
OPEN(9,file="sed")
OPEN(10,file="tss")
OPEN(11,file="recrd")
OPEN(12,file="another")
OPEN(13,file="oo")

ALLOCATE(ys(nos),bed(nos),dists(nos),tau(nos),ks(nos),tbst(nos),& 
         recrd(0:nos), bedlast(nos),hss(nos),tss(nos),hss2(nos),Qe(nos),&
         Qd(nos),f(nos),slopes(nos),& 
         NN(nos), C0(nos),taucrit_dep(nos, layers), C(nos), &
         taucrit_dep_ys(nos),dst(nos,0:(layers+1)), &
         taucrit(nos, 0:layers), slpmx(nos,0:(layers+1)), vegdrag(nos),&
         ysold(nos) ,  Qbed(nos), dqbeddx(nos),sllength(nos), &
         vel(nos), tau_g(nos), f_g(nos), Cbar(nos), bedold(nos), a_ref(nos)) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Loop over different values of discharge

DO Q_loop= 1, no_discharges!15 

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
        !conditions to be valid, we need the cross section edges not to ever
        !be sumberged. 
        DO i =1,nos
        
            ! Cross-channel distance, 'y' coordinate
            ys(i)=Width/(nos*1._dp-1._dp)*(i*1._dp-1._dp) 
            
            ! Bed elevation
            !bed(i)= 1.11_dp -8.01_dp*(abs(-1._dp+1._dp*((ys(i)-Width/2._dp))**6/(Width/2._dp)**6)) !-3.01 !+0.00002*(ys(i)-ys(1))*(ys(i)-nos*ys(1))
            !bed(i)= -4.1_dp -0.01_dp*(abs(-1._dp+1._dp*((ys(i)-Width/2._dp))**2/(Width/2._dp)**2)) !-3.01 !+0.00002*(ys(i)-ys(1))*(ys(i)-nos*ys(1))
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
    recrd=0._dp !A useful variable for storing stuff 
    l=1 !variable for "lower" wet point of cross section
    u=nos !variable for "upper" wet point of cross section
    water=waterM
    ! Redefine l and u
    call wet(l,u,nos,water, bed)
    tau= 0._dp*tau !Set all the shears to zero prior to update
    E=0._dp ! Cross-sectionally integrated resuspension 
    D=0._dp ! Cross-sectionally integrated deposition (from suspension)
    C=sconc ! Near bed sediment concentration
    Cbar = sconc ! Depth averaged suspended sediment concentration 
    f= 0.02_dp !A starting value, which will be changed
    f_g= 0.02_dp !A starting value, which will be changed
    vel = 0.0_dp
    rmult=sum(f)/(8._dp*9.8_dp*nos) !A starting value, which will be changed
    inuc=0._dp  
    NN=tau*0._dp
    vlast=0._dp
    hss=0._dp
    hss2=0._dp
    Area=0.0000_dp
    vegdrag=0._dp
    dst(:, layers) = 2000._dp
    slpmx(:, layers+1) = 99999999._dp
    dqbeddx=0._dp
    bedlast = bed - 0.001_dp !
    bedold = bed - 0.001_dp !
    DT1 = dt ! An adaptive timestep
    a_ref=0.0_dp !Reference level of vanrijn

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

        ! Print out to console/check if converged
        22222 IF( mod(j-1,writfreq).eq.0 ) THEN 
                PRINT*, j, l,u, Q/Area, t, ((Q/Area)*abs(Q/Area)*rmult),&
                     DT1, ys(u)-ys(l)+wdthx, maxval(C), &
                    rmult*(Area)/(ys(u)-ys(l)+wdthx), f(nos/2)
              END IF

        ! Update time
        t=t+DT1

        ! Update the critical shear layers to account for any slope effects
        DO i = 1, nos
            DO jj= 0, layers
                !if((bed(i)> 0._dp).and.(jj>1)) mu= mu*100.5_dp
                multa=1._dp
                IF(.FALSE.) THEN !A convenient way to make clear that this code doesn't matter
                    aa= mu**2*(mu*lifttodrag-1._dp)/(mu*lifttodrag+1._dp)
                    bb= -2._dp*mu**2*lifttodrag*cos(atan(slopes(i)))*mu/(1._dp+mu*lifttodrag) 
                    cc= mu**2*cos(atan(slopes(i)))**2 - sin(atan(slopes(i)))**2
                    !!!!if((bed(i)> 0._dp).and.(jj>1)) mu= mu/100.5_dp
                    multa= (-bb - sqrt(bb**2-4._dp*aa*cc))/ (2._dp*aa) !We multiply this by the critical shear on a flat bed to get the critical shear on a slope. 
                    !Correct for any problems
                    !IF((multa<0.2_dp).OR.(isnan(multa))) THEN
                    !    multa=0.2_dp !1000._dp*epsilon(multa)
                    !     !  print*,  'multaprob'
                    !      !stop
                    !END IF
                    IF((j==1).and.(i==1).and.(jj==0)) PRINT*, 'WARNING: Slope dependent critical shear stress' 
                END IF !!GET RID OF THE ABOVE CHUNCK OF CODE

                !if(bed(i)>-100.3_dp) multa=1._dp
                !if((jj>-1).AND.(abs(slopes(i))<mu)) multa=1._dp
                taucrit(i,jj) = erconst*(1._dp+ jj*1._dp)*max(multa,1.0e-01_dp)

                IF( isnan(taucrit(i,jj))) THEN
                    PRINT*, "taucrit(", i,",", jj, ") is nan"
                    STOP
                END IF

            END DO 
        END DO

        !FIXME: Outdated 
        !Update the maximum slope that layers are permitted at (with regard to the geotech routine)
        !IF(.FALSE.) THEN
        !    DO i = 1, nos
        !        DO jj= 0, layers
        !            slpmx(i,jj)=smax
        !            IF(bed(i)>-1._dp) slpmx(i,jj)=smax !4._dp
        !            !end if
        !        END DO 
        !    END DO
        !END IF


        !!Update taucrit_dep
        !$OMP PARALLEL DO
        DO  i= 1, nos
            DO n= 1, layers
                taucrit_dep(i,n)= max(min(taucrit_dep(i,n), bed(i)),bed(i)-lincrem*n) !So sediments buried at a certain depth will have their critical shear increase. 
            END DO !N 
        END DO


        !dst - -the distance of any critical shear layer to the surface.
        dst(:,0)= 0._dp
        DO i= 1, nos
            DO jj= 1, layers
                dst(i, jj)=max((bed(i)-taucrit_dep(i,jj) ), 0._dp)
                ! print*, hs(i), taucrit_dep(i, 1), dst(i, 1)
            END DO
        END DO
        !$OMP END PARALLEL DO


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

        ! Water elevation (free surface)
        waterlast=water
        water= waterM + (TR/2._dp+.0*sin(2._dp*pi*t/(12.4_dp*3600._dp*100._dp))+0.0_dp)*sin(2._dp*pi*t/(12.4_dp*3600._dp))  !1.+ 1.*sin(2.*pi*j/500.) 

        ! Find wetted part of section
        call wet(l,u,nos,water, bed) 
        
        ! 'Extra' wetted width associated with the corners of the domain.
        wdthx = 0.0
        IF (u<nos) wdthx = wdthx+ (water-bed(u))/(bed(u+1)-bed(u))*(ys(u+1)-ys(u))
        IF (l>1) wdthx = wdthx+ (water-bed(l))/(bed(l-1)-bed(l))*(ys(l)-ys(l-1))

        ! Cross sectional area
        Arealast=Area
        Area=0.0_dp
        IF(l>0) THEN 
            DO i = l,u-1
                Area=Area+ 0.5_dp*( max( (water-bed(i)),0._dp) +max( (water-bed(i+1)), 0._dp) )*(ys(i+1)-ys(i))
            END DO
            !!Add edge residuals if needed
            IF(u<nos) Area=Area+0.5_dp*(water-bed(u))**2/(bed(u+1)-bed(u))*(ys(u+1)-ys(u))
            IF(l>1) Area = Area+ 0.5_dp*(water-bed(l))**2/(bed(l-1)-bed(l))*(ys(l)-ys(l-1))

        END IF

        IF(isnan(Area)) THEN
            PRINT*, "Area is nan", l, u, water, maxval(bed(l:u)), minval(bed(l:u)), maxval(tau), minval(tau)
            STOP
        END IF

        IF ((Area<0.02).OR.(maxval(water-bed(l:u))<.05)) THEN 
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

        IF(isnan(Q)) THEN
            print*, "Q is nan", Q, Area, Arealast
        END IF

        IF((Q/Area)>5.) THEN
            print*, "Q/A>5", Q/Area, Q, Area, Arealast, j,  dt, water, waterlast, width, t, dt!, l, u, water-bed(l:u)
            !stop     
        END IF



        !!!!!!!!!!!!!!!
        !CALCULATE SHEAR AND UPDATE THE CROSS SECTION
        !!!!!!!!!!!!!!
        recrd=0._dp
        tau=0._dp
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
        
            ! Calculate friction
            call calc_friction(friction_type, grain_friction_type, rough_coef, water, u-l+1,&
                                 bed(l:u), vel(l:u), man_nveg,d50,veg_ht, rhos, rho, g,&
                                 f(l:u), vegdrag(l:u),f_g(l:u), dsand, j, a_ref(l:u)) 

            ! Calculate bed shear
            call calc_shear(u-l+1,DT1,water,Q,bed(l:u),ys(l:u),Area,ys(u)-ys(l)+wdthx, &
                            water-Area/(ys(u)-ys(l)+wdthx),f(l:u),recrd((l-1):u),E,D, C(l:u),&
                            rmult,inuc, tau(l:u),& 
                            NN(l:u),j,slopes(l:u), hlim, mor, taucrit_dep(l:u,1:layers), layers,&
                            taucrit_dep_ys(l:u) & 
                            ,u-l+1, taucrit(l:u, 0:layers) , vegdrag(l:u), susdist, rho, Qe(l:u) & 
                            , Qbed(l:u), rhos, voidf, d50, g, kvis, norm, vertical, lambdacon, tbston &
                            , ysl,ysu,bedl,bedu, high_order_shear) 
            
            ! Calculate depth-averaged velocity
            vel = 0._dp 
            vel(l:u)=sqrt(abs(tau(l:u))/rho*8._dp/f(l:u))*sign(1._dp+0._dp*tau(l:u), tau(l:u))
            
            tau_g = 0._dp
            !tau_g(l:u) = 0.5_dp*rho*vel(l:u)**2*(f_g(l:u))*sign(1._dp+0._dp*tau(l:u), tau(l:u))
            !Following Abdel-Fattah et al 2004
            tau_g(l:u) = rho*vel(l:u)**2*(f_g(l:u)/8._dp)*sign(1._dp+0._dp*tau(l:u), tau(l:u))
            
            ! Calculate rates of resuspension and bedload transport
            call calc_resus_bedload(u-l+1,DT1,water,Q,bed(l:u),ys(l:u),Area,ys(u)-ys(l)+wdthx,&
                                     water-Area/(ys(u)-ys(l)+wdthx),f(l:u),recrd((l-1):u),E,D, &
                                    C(l:u),wset, rmult,2,inuc, tau_g(l:u),& 
                                    vel(l:u), NN(l:u),j,slopes(l:u), hlim, mor, taucrit_dep(l:u,1:layers),&
                                     layers, taucrit_dep_ys(l:u) & 
                                    ,u-l+1, taucrit(l:u, 0:layers) , vegdrag(l:u), susdist, rho, Qe(l:u) & 
                                    , Qbed(l:u), rhos, voidf, dsand, d50, g, kvis, norm, vertical,alpha, &
                                    tbston, Qbedon, ysl,ysu,bedl,bedu, resus_type, bedload_type, a_ref(l:u)) 

            !Calculate the cross-sectional suspended load distribution
            IF(susdist) THEN
                !! Adjust the erosion factor if erosion is normal to the bed
                IF(norm) THEN
                    sllength = sqrt(1._dp+ slopes**2)
                ELSE
                    sllength=1._dp
                END IF
            
                !IF(.FALSE.) THEN
                !    C=0._dp
                !    call sus_dist_sect(u-l+1, ys(l:u), bed(l:u), water, tau(l:u), vel(l:u), sconc, wset,lambdacon, &
                !            rho,rhos, g, d50, susQbal, Qbed(l:u),talmon, sllength, bedl, bedu, ysl, ysu, C(l:u),&
                !            integrated_load_flux)
                !ELSE
                call dynamic_sus_dist(u-l+1, DT1, ys(l:u), bed(l:u), water, waterlast, Q, tau(l:u), vel(l:u), wset, & 
                                        Qe(l:u), lambdacon, rho,rhos, g, d50, bedl,bedu, ysl, ysu, C(l:u),&
                                        Cbar(l:u), Qbed(l:u), sconc, j, high_order_Cflux, a_ref(l:u), sus_vert_prof)

                ! Set C in dry parts of the channel to zero
                IF(l>1) THEN
                    Cbar(1:(l-1)) = 0._dp
                    C(1:(l-1)) = 0._dp
                END IF
                
                IF(u<nos) THEN
                    Cbar((u+1):nos) = 0._dp
                    C((u+1):nos) = 0._dp
                END IF
        

            END IF
           
        
            
            bedlast= bed ! Record the bed prior to updating
           
            call update_bed(u-l+1,DT1,water,Q,bed(l:u),ys(l:u),Area,ys(u)-ys(l)+wdthx, &
                             water- Area/(ys(u)-ys(l)+wdthx),f(l:u),recrd((l-1):u),E,D, &
                            C(l:u),rmult,2,inuc, tau(l:u),& 
                            NN(l:u),j,slopes(l:u), hlim, mor, taucrit_dep(l:u,1:layers), &
                            layers, taucrit_dep_ys(l:u) & 
                            ,u-l+1, taucrit(l:u, 0:layers) , vegdrag(l:u), susdist, rho, &
                            Qe(l:u), Qbed(l:u), wset, dqbeddx(l:u), rhos, voidf, d50, g, &
                            kvis, norm, vertical, lambdacon, tbston,&
                            Qbedon, normmov, sus2d, ysl, ysu, bedl,bedu, iii, bedlast(l:u), &
                            susQbal, talmon, high_order_bedload) 

            ! Correct the banks. In the case that we allow bedload at l-1/2 and
            ! u+1/2, this is very important to ensure mass conservation, because
            ! if there is a downslope bedload flux from l-1/2, or from u+1/2,
            ! then it must come from the dry part of the channel 
            IF(.FALSE.) THEN
                ! Use this case when bedload occurs and l-1/2, u+1/2
                IF(l>1) THEN
                    IF(bed(l-1)>bedl) bed(l-1)=bedl
                END IF
                IF(u<nos) THEN
                    IF(bed(u+1)>bedu) bed(u+1)=bedu
                END IF
            END IF

            IF(.TRUE.) THEN
                ! A version of the Delft bank erosion model. 
                ! First check that there is no leakage of bedl, bedu in the
                ! bed solver (possibly could happen due to matrix round off or
                ! coding error).
                IF((l>1).and.(u<nos)) THEN
                    IF((abs(bedl - bedlast(l-1))>0.0e-8_dp).OR.& 
                       (abs(bedu-bedlast(u+1))>0.0e-08_dp)) THEN
                        print*, 'ERROR -- there is still erosion of dry points in &
                            the matrix solution of the bed solver'
                        print*, abs(bedl - bedlast(l-1)), abs(bedu-bedlast(u-1)), l, u
                        stop
                    END IF
                END IF
                ! If erosion is occuring at the channel margins,
                ! then assign it to the neighbouring dry bed point
                IF((bed(l)<bedlast(l)).AND.(l>1)) THEN
                    !IF( abs(bed(l) - bed(l-1))/(ys(l)-ys(l-1))>1.0_dp) THEN
                    !IF( abs(tau_g(l))>taucrit(l,0)) THEN
                        bed(l-1) = bed(l-1) - (bedlast(l) - bed(l))
                        bed(l) = bedlast(l)
                    !END IF
                END IF
                IF((bed(u)<bedlast(u)).AND.(u<nos)) THEN
                    !IF( abs(bed(u+1) - bed(u))/(ys(u+1)-ys(u))>1.0_dp) THEN
                    !IF( abs(tau_g(u))>taucrit(u,0)) THEN
                        bed(u+1) = bed(u+1) - (bedlast(u) - bed(u))
                        bed(u) = bedlast(u)
                    !END IF
                END IF
            END IF
            
            ! Write outputs
            IF((mod(j-1,writfreq).EQ.0).AND.(j>0).AND.(iii.eq.1)) THEN!.or.(j>15250)) THEN 
                write(1,*) tau 
                write(2,*) bedlast !Same bed as when tau was calculated
                write(3,*) ys !critical shear
                write(4,*) tau_g !taucrit_dep!water, Q/ar
                write(5,*) water
                write(7,*) Qe
                write(8,*) Qbed
                write(9,*) C
                write(10,*) vel
                write(11,*) recrd
                !write(12,*) taucrit_dep_ys

                ! Check for convergence.
                tmp = maxval(abs(bedold-bed))
                IF(tmp<1.0e-04_dp) THEN
                    goto 373 !Converged: Go to the end of this loop
                    !exit
                    !sconc = sconc*0.5_dp
                END IF
                bedold=bed
            END IF
            
            ! Determine the timestep -- implicit timestepping will only work
            ! if there are no bed layers, because otherwise DT
            ! has already been used this timestep.
            IF(variable_timestep) THEN
                ! These methods change DT1 during the evolution of the
                ! cross-section
                
                !FIXME: At the moment, this slows the morphological evolution,
                !but has no effect on the actual time t - so it can only produce
                !the correct solutions for steady state computations

                IF(j.eq.1) print*, ' Warning: Variable timestep is ONLY valid for &
                        STEADY UNIFORM EQUILIBRIUM computations'

                !IF(.FALSE.) THEN
                !    ! Change the timestep according to the rate of deposition
                !    DT1 = min(DT,0.003_dp/max(abs(maxval(C(l:u)*wset/rhos)*mor), 0.0000001_dp)) 
                !ELSEIF(.FALSE.) THEN
                !    ! Change the timestep according to the rate of morphological
                !    ! change
                !    DT1 = min(DT,0.003_dp/max(maxval(abs(bed(l:u) - bedlast(l:u))), 0.0000001_dp)) 
                !ELSEIF(.TRUE.) THEN
                    !tmp = max(maxval(abs(wset*C/rhos- Qe))/1.0_dp, maxval(abs(bed - bedlast)/DT1))
                    tmp = min(max(maxval(abs(wset*C/rhos- Qe)), maxval(abs(recrd(l-1:u)))), &
                              maxval(abs(bed(l+1:u-1) - bedlast(l+1:u-1))) )
                    tmp2 = minval(ys(2:nos) - ys(1:nos-1)) 
                    !tmp = max(maxval(abs(wset*C/rhos- Qe)), maxval(abs(bed - bedlast)))
                    !tmp = max(maxval(abs(recrd(l-1:u)))*0.1_dp, maxval(abs(bed - bedlast)))
                    !tmp = max(maxval(abs(bed - bedlast)/DT1), 0.00001_dp)
                    !DT1_old = DT1
                    !tmp = maxval(abs(wset*C/rhos- Qe))*1.00_dp
                    DT1 = min(max(1.0e-03_dp*tmp2/max(tmp,1.0e-020_dp), 1.0_dp), 100.0_dp*3600.0_dp)
                    ! Get bed to accelerate 
                    !mor = min(max(3.0e-03_dp/(maxval(abs(bed-bedlast))/DT1_old*DT1), 1.0_dp), 5._dp)
                !END IF

            ELSE
                DT1 = DT
            END IF 
            ! Add random sedimentation to the channel
            ! This is a technique to prevent the channel from converging to a
            ! marginally stable state where tau_g < taucrit, which is dependent
            ! on the imposed initial condition. This way we force tau>=taucrit
            ! in a stable channel, which seems to allow the same solution
            ! irrespective of the initial geometry (RIGOROUSLY CHECK THIS)
            IF((.FALSE.).AND.(mod(j,5000).eq.0).AND.(1.0_dp*j<jmax*1.0_dp)) THEN
                
                call random_seed(size = iii2)
                call random_number(tss)
                tss(3:nos-2) = 0.1_dp*(tss(1:nos-4) + tss(2:nos-3)+tss(3:nos-2)+tss(4:nos-1)+tss(5:nos))/5._dp
                tss(2) = tss(3)
                tss(1) = 0.0_dp
                tss(nos-1) = tss(nos-2)
                tss(nos) = 0.0_dp
   
                ! Zero the perturbation where taug > taucrit
                DO i=2,nos-1
                    IF(max(tau_g(i), tau_g(i+1), tau_g(i-1))>taucrit(i,0)) tss(i) = 0._dp
                    IF(bed(i)>=water) tss(i) = 0._dp
                END DO

                PRINT*, 'Smooth random deposition of max ', maxval(tss), ' metres'
                
                DO i=1,nos
                    bed(i) = min(bed(i)+tss(i), water)
                END DO
            END IF

            ! Alternative method to 'wipe out' points which can't erode.
            !DO i=2,nos-1
            !    IF(tau_g(i)<taucrit(i,0)) bed(i) = min(bed(i) + & 
            !                            max((water-bed(i))*0.01_dp, min(0.02_dp, (water-bed(i))*0.1_dp)), water) 
            !END DO
                
            !bed(l:u)= min(bed(l:u), water-0.00001_dp)

            !slopes((l):(u)) = (bed(l+1:u+1)-bed(l-1:u-1))/(ys(l+1:u+1)-ys(l-1:u-1))
            !slopes(l)=(bed(l+1)-bed(l))/(ys(l+1)-ys(l))
            !slopes(u)=(bed(u-1)-bed(u))/(ys(u-1)-ys(u))
        END DO

        ! Basic limiting of the channel slope -- to circumvent the numerically
        ! difficult problem of allowing infinite banks otherwise
        hss = bed
        IF((remesh.eqv..FALSE.).AND.(.TRUE.).AND.(mod(j,100)==0)) THEN
            !FIXME: At present, this is only valid(mass conserving) with a uniform mesh
            ! Move from centre of channel to left bank
            DO i=nos/2,2,-1
                IF(abs(bed(i)-bed(i-1))>(ys(i)-ys(i-1))) THEN
                    IF(bed(i)>bed(i-1)) THEN
                        tmp = bed(i)-(bed(i-1) + ys(i)-ys(i-1))
                        hss(i) = hss(i)-tmp*0.5_dp
                        hss(i-1) = hss(i-1) + tmp*0.5_dp    
                    ELSE
                        tmp = bed(i-1)-(bed(i) + ys(i)-ys(i-1))
                        hss(i) = hss(i)+tmp*0.5_dp
                        hss(i-1) = hss(i-1) - tmp*0.5_dp    

                    END IF
                    
                END IF
            END DO
            ! Move from centre of channel to right bank
            DO i=nos/2,nos-1,1
                IF(abs(bed(i)-bed(i+1))>(ys(i+1)-ys(i))) THEN
                    IF(bed(i)>bed(i+1)) THEN
                        tmp = bed(i)-(bed(i+1) + ys(i+1)-ys(i))
                        hss(i) = hss(i)-tmp*0.5_dp
                        hss(i+1) = hss(i+1) + tmp*0.5_dp    
                    ELSE
                        tmp = bed(i+1)-(bed(i) + ys(i+1)-ys(i))
                        hss(i) = hss(i)+tmp*0.5_dp
                        hss(i+1) = hss(i+1) - tmp*0.5_dp    

                    END IF
                    
                END IF
            END DO

        !Reset the bed
        bed = hss
        END IF

        !DO i=1,floor(nos*0.50_dp)
        !    IF(abs(bed(i)-bed(nos-i+1))>1.0e-12_dp) THEN
        !        print*, 'bed symmetry loss post perturbation'
        !        stop
        !    END IF
        !END DO        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Delft bank erosion model
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!Erosion at edge wet points is assigned to neighbouring dry points.
        !if((l>1).and.(u<nos)) THEN
        !if(bed(l)<bedlast(l)) THEN
        !bed(l-1)=bed(l)-(bedlast(l)-bed(l))
        !bed(l)=bedlast(l)
        !end if
        !if(bed(u)<bedlast(u)) THEN
        ! bed(u+1)=bed(u)-(bedlast(u)-bed(u))
        !bed(u)=bedlast(u)
        !end if
        !end if



        !Remeshing. If the points are becoming too spaced apart, we might have to do this
        IF(remesh) THEN
            IF(mod(j, remesh_freq).EQ.0) THEN
                !Get distance between consecutive points in the cross seciton. 
                !dists(2:nos)=sqrt((ys(2:nos)-ys(1:nos-1))**2+(bed(2:nos)-bed(1:nos-1))**2) 
                !dists(1)=dists(2)
                !handy=minval(dists)
                !DO i =2, nos-1
                    !IF(abs(slopes(i))>.0_dp) THEN
                        
                    !    IF(dists(i)/handy > 0.0_dp) THEN
                            ! Perform remeshing here
                            ysold=ys !Predefine
                            call refit(ys, bed, nos) !Remesh so things are more even.
                            ! Also remesh old timestep records of the bed
                            call interp3(ysold, bedlast,ys,nos)
                            call interp3(ysold, bedold, ys, nos)
                            ! Also remesh the suspended sediment concentration.
                            call interp3(ysold, C, ys, nos)
                            call interp3(ysold, Cbar, ys, nos)
                            !Now fix up the layers
                            DO jj=1, layers
                                call interp(ysold,  taucrit_dep(:,jj),ys, nos) !Beware the risk that this could cause 'leakage'
                                taucrit_dep(:, jj)= min(bed, taucrit_dep(:,jj)) !So if we interpolate the critical shear layer, it should of course not be above the bed level.
                            END DO
                            ! Break out of loop                        
                            !EXIT

                    !    END IF
                    !END IF
                !END DO

            END IF 
        END IF !REMESH

    END DO ! End of timestepping loop 
        
        373 print*, "converged, Q_loop=", Q_loop
        !stop
END DO ! End of loop for different Q_loops

close(1)
close(2)
close(3)
close(4)
close(5)
close(6)
close(7)
close(8)
close(9)
close(10)

END PROGRAM      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

