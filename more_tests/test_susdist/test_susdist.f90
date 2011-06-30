PROGRAM test_susdist
!Program to provide unit tests for the cross sectional suspended sediment
!routines
IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
!USE crosssection, ONLY: dp

! Useful variables
INTEGER:: a = 1000, i, j
REAL(dp):: lengths, heights, water, tau, vel, C, Cbar, Qbed, sllength,&
             analytical
REAL(dp):: lengthsl,lengthsu,heightsl, heightsu, f, Sf, dT, Q

!Namelist constants
REAL(dp)::susconcs(1000), rho, wset, rhos, dsand, d50, g, kvis, lambdacon, alpha, &
     advection_scale_const, width, TR, integrated_load_flux
NAMELIST /inputdata/ dT,susconcs, rho, wset, rhos, dsand, d50, g, kvis, lambdacon,&
         alpha, width, TR, integrated_load_flux


ALLOCATABLE lengths(:), heights(:), tau(:), vel(:), Qbed(:), C(:), Cbar(:), &
                sllength(:), analytical(:)
ALLOCATE(lengths(a), heights(a), tau(a), vel(a), Qbed(a), C(a), Cbar(a), &
            sllength(a), analytical(a))

!!Read from standard input . To do this, type  ./model < input_file ( and then you can pipe to an output by adding > outfile.log)
READ(*,nml=inputdata)
PRINT inputdata

PRINT*, 'Test the diffusive part of the dynamic single cross-section &
         suspended load equation. &
         This will only work if erosion and deposition are switched off &
         in dynamic_sus_dist'

!!!!!!!!!!!!!!!!!!!!!
!
!Initialise geometry
!
!!!!!!!!!!!!!!!!!!!!!

lengths = (/ ((i-1)*width/(1.0_dp*(a-1)), i=1,a) /)
heights = ( (lengths-width/2._dp)/(width/2._dp) )**2 !Parabolic section
heightsl = heights(1)+1.0e-01_dp
heightsu = heights(a) + 1.0e-01_dp
lengthsl = lengths(1)-1.0e-03_dp
lengthsu = lengths(a) +1.0e-03_dp

! Useful constants
water = 1.001_dp ! Water elevation
IF(water< maxval(heights)) print*, 'Water < maxval(heights)'

f=0.01_dp  ! Friction factor (Darcy)
Sf = 0.0004_dp      ! Friction slope

!Make up some useful variables
tau = rho*g*(water - heights)*Sf
vel = sqrt(tau/(rho*f/8._dp)) !Vel = sqrt(tau/rho(f/8))
Q = sum(vel*max(water - heights,0.0_dp))*(lengths(2)-lengths(1))
Qbed =0._dp
sllength = 1._dp
C = 1.0e-07_dp
Cbar = C

print*, 'dT = ', dT
print*, 'Q = ', Q
print*, 'velmax = ', maxval(vel)

open(3,file='test_susdist_outfile.out')

DO j=1,20
    print*, 'j = ', j
    DO i=1,500
   
        ! Note we set
        ! Qe = wset*C (erosion = deposition)
        ! waterlast = water (Removes any cross-channel velocity)
        ! 
        call dynamic_sus_dist2(a, dT, lengths, heights, water, water,&
                 Q, tau, vel, wset, 0.0*Qbed, lambdacon, rho,rhos, g, d50, &
                 heightsl,heightsu, lengthsl, lengthsu, C, Cbar, Qbed, susconcs(1), i+(j-1)*1000,.FALSE.)
            !SUBROUTINE dynamic_sus_dist(a, delT, ys, bed, water, waterlast, Q, tau, vel, wset, Qe,lambdacon, &
            !                                rho,rhos, g, d50, bedl,bedu, ysl, ysu, cb, Cbar, Qbed, &
            !                                sconc, counter, high_order_Cflux)

        !call dynamic_sus_dist(u-l+1, dT, lengths(l:u), heights(l:u), water, waterlast,&
        !         Q, tau(l:u), vel(l:u), wset, Qe(l:u), lambdacon, rho,rhos, g, d50, &
        !        heightsl,heightsu, lengthsl, lengthsu, C(l:u), Cbar(l:u), & 
        !        advection_scale_const, j)
    END DO

    ! Compute analytical expression for Cbar in this special case
    analytical = 1.0_dp/(water-heights)* & 
                (exp(wset/(sqrt((water-heights)*g*Sf)*0.1_dp) ) -1._dp )**(-2.0_dp) 

    DO i=1,a
        write(3,*) tau(i), vel(i), Cbar(i), analytical(i), Cbar(i)/analytical(i)  
    END DO
END DO
close(3)

PRINT*, 'Finished: Now compare Cbar and analytical -- my initial test found &
    that they agreed well' 

END PROGRAM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HERE IS A COPY OF dynamic_sus_dist from shionocros.f90. The key modification
! is that erosion and deposition are cut out. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dynamic_sus_dist2(a, delT, ys, bed, water, waterlast, Q, tau, vel, wset, Qe,lambdacon, &
                                rho,rhos, g, d50, bedl,bedu, ysl, ysu, cb, Cbar, Qbed, &
                                sconc, counter, high_order_Cflux)
    ! Calculate the cross-sectional distribution of suspended sediment using
    ! some ideas from /home/gareth/Documents/H_drive_Gareth/Gareth_and_colab
    ! s/Thesis/Hydraulic_morpho_model/channel_cross_section/paper/idea_for_
    ! simplified_cross_sectional_sus
    
    ! We solve the equation
    ! depth d (Cbar) / dt + U*depth dCbar/dx +
    ! V*depth* d(Cbar)/dy - d/dy( eddif_y d(depth Cbar)/dy + eddif_y*cb*dh/dy) -
    ! (Es - ws*cb) = 0.
    ! With suitable approximations for dCbar/dx and V
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
    !USE crosssection, ONLY: dp
    INTEGER, INTENT(IN)::a, counter
    REAL(dp), INTENT(IN):: delT, ys, bed, water, waterlast, tau, vel, wset, Qe, lambdacon, rho, rhos,g, & 
                                d50, bedl, bedu,ysl,ysu, sconc, Q, Qbed
    REAL(dp), INTENT(IN OUT):: cb, Cbar ! Near bed suspended sediment concentration, Depth averaged suspended sediment concentration
    LOGICAL, INTENT(IN):: high_order_Cflux
    DIMENSION ys(a), bed(a), tau(a),vel(a), Qe(a), cb(a), Cbar(a),Qbed(a) 

    ! LOCAL VARIABLES
    INTEGER:: i, info
    LOGICAL:: halt
    REAL(dp):: depth(0:a+1), eddif_y(0:a+1), eddif_z(a), zetamult(0:a+1), vd(0:a+1), ys_temp(0:a+1)
    REAL(dp):: M1_lower(a), M1_diag(a), M1_upper(a), M1_upper2(a), M1_lower2(a)
    REAL(dp):: RHS(a), dy_all(a)
    REAL(dp):: tmp1, dy, dy_outer, xlen, tmp2, Cref
    REAL(dp):: dQdx, dhdt, Cbar_old(a), dUd_dx(0:a+1)
    REAL(dp):: DLF(a), DF(a), DUF(a), DU2(a),rcond, ferr, berr, work(3*a), XXX(a, 1)
    REAL(dp):: bandmat(5,a), AFB(7,a), RRR(a), CCC(a)
    INTEGER::  IPV(a), iwork(a)   
    LOGICAL:: const_mesh
    CHARACTER(1):: EQUED
    ! This routine calculates C, Cbar in units m^3/m^3 --- however, elsewhere
    ! they are in kg/m^3 --- so we convert here, and convert back at the end of
    ! the routine
    cb = cb/rhos
    Cbar = Cbar/rhos 

    ! Depth, including at zero-depth boundaries ( 0 and a+1 )
    depth(1:a) = max(water -bed, 0._dp)
    depth(0) = 0._dp
    depth(a+1) = 0._dp

    ! y coordinate, including at zero depth boundaries (0 and a+1)
    ys_temp(1:a) = ys
    ys_temp(0) = ysl
    ys_temp(a+1) = ysu

    ! Lateral eddy diffusivity
    IF(lambdacon>0._dp) THEN
        eddif_y(1:a)= lambdacon*sqrt(abs(tau)/rho)*depth(1:a)
    ELSE
        eddif_y(1:a)= 0.2_dp*sqrt(abs(tau)/rho)*depth(1:a) 
    END IF
   
    IF(.FALSE.) THEN
        eddif_y(1:a)=maxval(eddif_y(1:a)) !eddif_y(1:a)+0.01_dp
        IF(counter.eq.1.) PRINT*, 'WARNING: constant eddy diff'
    END IF 
    
    ! Include zero depth boundaries 0 and a+1
    eddif_y(0) = 0._dp
    eddif_y(a+1) = 0._dp

    IF(.FALSE.) THEN
        eddif_y=0._dp
        IF(counter.eq.1) print*, 'WARNING: Zero eddy diffusivity in dynamic_sus_dist'
    END IF
    
    ! Vertical eddy diffusivity
    eddif_z= 0.1_dp*sqrt(abs(tau)/rho)*depth(1:a) 
    
    IF(.FALSE.) THEN
        eddif_z(1:a)=maxval(eddif_z(1:a)) !eddif_y(1:a)+0.01_dp
        IF(counter.eq.1.) PRINT*, 'WARNING: constant eddy diff'
    END IF 
    
    !zetamult*cbed= depth integrated sediment concentration = Cbar*d
    DO i=1, a
        IF((eddif_z(i)>0._dp).and.(depth(i)>0.0e-00_dp)) THEN 
            zetamult(i)= eddif_z(i)/wset*(1._dp-exp(-(wset/eddif_z(i))*max(depth(i),0._dp)) )
        ELSE 
            zetamult(i)=1.0e-012_dp !1.0e-04_dp
        END IF
    END DO
    zetamult = max(zetamult, 1.0e-012)
    ! Include zero depth boundaries 0 and a+1 -- set to 1, because we divide by
    ! it often
    zetamult(0)   = 1.0e-012_dp
    zetamult(a+1) = 1.0e-012_dp

    !zetamult=1.0_dp

    ! Solve initially for the depth - averaged suspended sediment concentration
    ! depth d (Cbar) / dt + U*depth*dCbar/dx +
    ! V*depth* d(Cbar)/dy - d/dy( eddif_y d(depth Cbar)/dy + eddif_y*cb*dh/dy) -
    ! (Es - ws*cb) = 0.
    !
    ! We treat everything implicitly except the dCbar/dx term
    
    
    M1_diag = 0._dp
    M1_lower = 0._dp
    M1_lower2 = 0._dp
    M1_upper = 0._dp
    M1_upper2 = 0._dp
    RHS = 0._dp
   
    !
    ! Fill out matrices term-by-term
    !

    ! depth*d(Cbar)/dt
    DO i = 1, a
        M1_diag(i) = M1_diag(i) + depth(i)/delT
        RHS(i)     = RHS(i)     + Cbar(i)*depth(i)/delT
    END DO
        
    ! Ud Advection term
    !IF(.FALSE.) THEN
    !
    !    ! Model 1 -- dC/dx is positive in the direction of velocity if Es>Ds, and negative otherwise
    !    ! THIS DIDN'T WORK!
    !    DO i = 1, a

    !        tmp1 = depth(i) 
    !        if(vel(i)==0._dp) tmp1 = 0._dp

    !        IF(zetamult(i) > 0._dp) THEN
    !            M1_diag(i) = M1_diag(i) - tmp1*wset*(depth(i)/zetamult(i)) ! Note depth(i)/zetamult(i)*Cbar = cb
    !        END IF

    !        RHS(i)     = RHS(i)     - tmp1*Qe(i)     

    !    END DO
    !ELSE

    ! Model 2 -- dC/dx is = (C - C/mean_C*desired_C)/x_length_scale
    !tmp1 = sum(Cbar*vel*depth(1:a)*(0.5_dp*(ys_temp(2:a+1)-ys_temp(0:a-1)))) !Cbar flux
    !tmp2 = 3.333333e-05_dp*abs(Q) ! Desired sediment flux
    dy_all = (ys_temp(2:a+1)-ys_temp(0:a-1))*0.5_dp
    tmp1 = sum(Cbar*vel*depth(1:a)*dy_all) !Cbar flux
    tmp2 =sconc*abs(Q) - sum(Qbed(1:a)*dy_all) !Desired Cbar flux = 'Measure of total load less bedload'
    tmp2 = max(tmp2, 0._dp)

    !print*, 'tmp2/tmp1 = ', tmp2/tmp1, tmp2, tmp1, maxval(dy_all), maxval(Cbar), maxval(vel), abs(Q)
    !Check if the mesh is constant -- because if it is, we can use high order
    !derivative methods
    IF(maxval(dy_all) - minval(dy_all) < 1.0e-010) THEN
        const_mesh=.TRUE.
    ELSE
        const_mesh=.FALSE.
    END IF 

    IF(mod(counter,1000).eq.0) PRINT*, 'sus flux is =', tmp1, '; desired flux is', tmp2
    xlen=1000._dp ! dx
    
    DO i=1,a
        !Time centred implicit
        !M1_diag(i) = M1_diag(i) + 0.5_dp*vel(i)*depth(i)*(1._dp-tmp2/tmp1)/xlen !*(1._dp - tmp2/max(tmp1,1.0e-06))/xlen        
        !RHS(i) = RHS(i) - 0.5_dp*vel(i)*depth(i)*(Cbar(i) - Cbar(i)*tmp2/tmp1)/xlen
        !Explicit
        Cref = Cbar(i)*tmp2/tmp1
        !Cref = tmp2/Q 
        !Cref = Cbar(i) + (tmp2-tmp1)/Q 
        !Cref = Cbar(i)*tmp2/tmp1 - xlen/max(vel(i),1.0e-02)*(Qe(i) - 0._dp*wset*depth(i)*Cbar(i)/zetamult(i))
        RHS(i) = RHS(i) - vel(i)*depth(i)*(Cbar(i) - Cref)/xlen

        !M1_diag(i) = M1_diag(i) + vel(i)*depth(i)*(1._dp/max(vel(i),1.0e-02)*wset*depth(i)/zetamult(i) )
    END DO 
             
    !END IF


    DO i = 1, a
        !PRINT*, 'SANITY CHECK a'
        IF(isnan(vd(i))) print*, 'vd(', i,') is NAN'
        IF(isnan(M1_diag(i))) print*, 'M1_diag(', i,') is NaN a'    
        IF(isnan(M1_lower(i))) print*, 'M1_upper(', i,') is NaN a'    
        IF(isnan(M1_upper(i))) print*, 'M1_upper(', i,') is NaN a'    

    END DO
    ! Vd Advection term -- requires prior calculation of vd, using an
    ! approximate model
    vd(0) = 0._dp ! lateral velocity*depth at zero depth boundary
    vd(a+1) = 0._dp 

    dUd_dx = 0._dp !Predefine
 
    dhdt = (water - waterlast)/delT
    dQdx = - (ys_temp(a+1) - ys_temp(0))*dhdt  ! dQ/dx = -B* dh/dt, from the cross-sectionally integrated continuity equation


    DO i = 1, a
        ! Estimate d(ud)/dx, using
        ! d(U*d)/dx ~ Ud/Q*dQ/dx --- approximate model for longitudinal momentum derivative
        IF(Q.ne.0._dp) dUd_dx(i) = (vel(i)*depth(i)/Q)*dQdx  
        ! Integrate to get Vd
        vd(i) = vd(i-1) - (ys_temp(i)-ys_temp(i-1))*(dhdt +0.5_dp*(dUd_dx(i)+dUd_dx(i-1)) )  ! Approximate model for Vd(y) = int_{left_bank}^{y} [- dh/dt -dUd/dx ] dy -- so we are integrating with a trapezoidal type rule
        IF(depth(i)==0._dp) vd(i) =0._dp
        ! Upwind discretization of the advective term (FIXME -- consider making
        ! use of a higher order discretization)
        IF((i.ne.1).and.(i.ne.a)) THEN
            IF(vd(i) < 0.0) THEN 
                dy = ys_temp(i+1) - ys_temp(i)
                M1_upper(i) = M1_upper(i) + vd(i)/dy
                M1_diag(i)  = M1_diag(i)  - vd(i)/dy
            ELSE
                dy = ys_temp(i) - ys_temp(i-1)
                M1_diag(i)  = M1_diag(i)  + vd(i)/dy
                M1_lower(i) = M1_lower(i) - vd(i)/dy
            END IF
        END IF

    END DO
    !write(12,*) vd(1:a)     
    DO i = 1, a
        !PRINT*, 'SANITY CHECK b'
        IF(isnan(vd(i))) print*, 'vd(', i,') is NAN'
        IF(isnan(M1_diag(i))) print*, 'M1_diag(', i,') is NaN b'    
        IF(isnan(M1_lower(i))) print*, 'M1_upper(', i,') is NaN b'    
        IF(isnan(M1_upper(i))) print*, 'M1_upper(', i,') is NaN b'    

    END DO

    ! Diffusion term
    DO i = 1, a
        ! Calculate a centred dy
        dy_outer = 0.5*(ys_temp(i+1) - ys_temp(i-1))
        
        ! d/dy ( eddify*d(depth Cbar)/dy)
        
        tmp1 = 0.5_dp*(eddif_y(i+1)+eddif_y(i))/(dy_outer*(ys_temp(i+1) - ys_temp(i)))
        IF(i<a) THEN
            IF((i>1).and.(i<a-1).and.(const_mesh).and.(high_order_Cflux)) THEN
                ! 4 point derivative approx -- estimate of
                ! 1/dy_outer*(eddify*d(depthCbar)/dy) at i+1/2
                M1_upper2(i) = M1_upper2(i) + 1.0_dp/24.0_dp*tmp1*depth(i+2)
                M1_upper(i) = M1_upper(i) - 9.0_dp/8.0_dp*tmp1*depth(i+1)
                M1_diag(i)  = M1_diag(i)  + 9.0_dp/8.0_dp*tmp1*depth(i)
                M1_lower(i) = M1_lower(i) - 1.0_dp/24.0_dp*tmp1*depth(i-1)
            ELSE
                ! 2 point derivative approx
                M1_upper(i) = M1_upper(i) - tmp1*depth(i+1)
                M1_diag(i)  = M1_diag(i)  + tmp1*depth(i)
            END IF
        END IF
 
        tmp1 = 0.5_dp*(eddif_y(i) + eddif_y(i-1))/((ys_temp(i) - ys_temp(i-1))*dy_outer)
        IF(i>1) THEN
            IF((i>2).and.(i<a).and.(const_mesh).and.(high_order_Cflux)) THEN
                ! 4 point derivative approx -- estimate of
                ! 1/dy_outer*(eddify*d(depthCbar)/dy) at i-1/2
                M1_upper(i) = M1_upper(i) - 1.0_dp/24.0_dp*tmp1*depth(i+1)
                M1_diag(i)  = M1_diag(i)  + 9.0_dp/8.0_dp*tmp1*depth(i)
                M1_lower(i) = M1_lower(i) - 9.0_dp/8.0_dp*tmp1*depth(i-1)
                M1_lower2(i) = M1_lower2(i) + 1.0_dp/24.0_dp*tmp1*depth(i-2)

            ELSE
                ! 2 point derivative approx
                M1_diag(i)  = M1_diag(i)  + tmp1*depth(i)
                M1_lower(i) = M1_lower(i) - tmp1*depth(i-1)
            END IF
        END IF

        ! d/dy ( eddify* dbed/dy * cb)

        ! First compute eddify*dbed/dy at i+1/2
        IF((i>1).and.(i<a-1).and.(const_mesh).and.(high_order_Cflux)) THEN
            !Use higher order dbed/dy approximation
            tmp1 = 0.5_dp*(eddif_y(i+1)+eddif_y(i))
            tmp1 = tmp1*(-1._dp)*(9.0_dp/8.0_dp*(depth(i+1)-depth(i)) &
                            -1.0_dp/24.0_dp*(depth(i+2)-depth(i-1)) )&
                                    /(ys_temp(i+1)-ys_temp(i))
            tmp1 = tmp1/dy_outer
        ELSE
            tmp1 = 0.5_dp*(eddif_y(i+1)+eddif_y(i))
            tmp1 = tmp1*( -(depth(i+1)-depth(i) )/(ys_temp(i+1)-ys_temp(i)))
            tmp1 = tmp1/dy_outer

        END IF

        IF(i<a) THEN
            !Estimate of 1/dy_outer*(eddify* dbed/dy *cb) at i+1/2
                M1_upper(i) = M1_upper(i) - 0.5*tmp1*(depth(i+1)/zetamult(i+1))  ! Note that depth(i)/zetamult(i)*Cbar = cb
                M1_diag(i)  = M1_diag(i)  - 0.5*tmp1*(depth(i)/zetamult(i))
        END IF

        IF((i>2).and.(i<a).and.(const_mesh).and.(high_order_Cflux)) THEN
            !Use higher order dbed/dy approximation
            tmp1 = 0.5_dp*(eddif_y(i)+eddif_y(i-1))
            tmp1 = tmp1*(-1._dp)*(9.0_dp/8.0_dp*(depth(i)-depth(i-1)) &
                            -1.0_dp/24.0_dp*(depth(i+1)-depth(i-2)) )&
                                    /(ys_temp(i+1)-ys_temp(i))
            tmp1 = tmp1/dy_outer
        ELSE
            ! Use central dbed/dy approximation
            tmp1 = 0.5_dp*(eddif_y(i)+eddif_y(i-1))
            tmp1 = tmp1*( -(depth(i)-depth(i-1) )/(ys_temp(i)-ys_temp(i-1)))
            tmp1 = tmp1/dy_outer
        END IF

        IF(i>1) THEN
            !Estimate of 1/dy_outer*(eddify*dbed/dy*cb) at i-1/2
            M1_diag(i)   = M1_diag(i)   + 0.5*tmp1*(depth(i)/zetamult(i))  ! Note that depth(i)/zetamult(i)*Cbar = cb
            M1_lower(i)  = M1_lower(i)  + 0.5*tmp1*(depth(i-1)/zetamult(i-1))
        END IF 

    END DO

    !Sanity-check
    DO i = 1, a
    !    PRINT*, 'SANITY CHECK c'
        IF(isnan(vd(i))) print*, 'vd(', i,') is NAN, c'
        IF(isnan(M1_diag(i))) print*, 'M1_diag(', i,') is NaN, c'    
        IF(isnan(M1_lower(i))) print*, 'M1_upper(', i,') is NaN, c'    
        IF(isnan(M1_upper(i))) print*, 'M1_upper(', i,') is NaN, c'    
    
    END DO

    ! Erosion and deposition
    IF(.FALSE.) THEN
        DO i = 1, a

            RHS(i) = RHS(i) +Qe(i)
            M1_diag(i) = M1_diag(i) + wset*(depth(i)/zetamult(i))  ! Note that depth(i)/zetamult(i)*Cbar = cb


        END DO
    END IF

    !Sanity-check
    DO i = 1, a
        !PRINT*, 'SANITY CHECK d'
        IF(isnan(vd(i))) print*, 'vd(', i,') is NAN'
        IF(isnan(M1_diag(i))) print*, 'M1_diag(', i,') is NaN d' , zetamult(i)   
        IF(isnan(M1_lower(i))) print*, 'M1_upper(', i,') is NaN d'    
        IF(isnan(M1_upper(i))) print*, 'M1_upper(', i,') is NaN d'    

    END DO

    ! Fill out matrix
    bandmat = 0._dp
    DO i=1,a
        IF(depth(i) == 0._dp) THEN
        ! Ensure zero depth has zero suspended sediment
            M1_diag(i) = 1._dp
            M1_upper(i) = 0._dp
            M1_lower(i) = 0._dp
            M1_upper2(i) = 0._dp
            M1_lower2(i) = 0._dp
            RHS(i) = 0._dp
        END IF
        ! Fill out matrix for solver
        IF(i<a-1) bandmat(1, i+2) = M1_upper2(i)
        IF(i<a) bandmat(2, i+1) = M1_upper(i)
        bandmat(3, i) = M1_diag(i)
        IF(i>1) bandmat(4, i-1) = M1_lower(i)
        IF(i>2) bandmat(5, i-2) = M1_lower2(i)
    END DO
   
    ! Solve equations
    !call DGTSV(a,1,M1_lower(2:a),M1_diag,M1_upper(1:a-1), RHS,a, info)
    ! New Cbar, converted to kg/m^3
    !Cbar = RHS*rhos
    
    !call DGTSVX('N', 'N', a, 1, M1_lower(2:a),M1_diag, M1_upper(1:a-1), &
    !    DLF(1:a-1), DF(1:a), DUF(1:a-1), DU2(1:a-2), IPV(1:a), RHS(1:a), a, &
    !    XXX(1:a,1), a, rcond, ferr, berr, work(1:(3*a)), iwork(1:a), info)
    
    ! New Cbar, converted to kg/m^3
    !Cbar = XXX(1:a,1)*rhos


    call DGBSVX('E','N', a,2,2,1, bandmat(1:5,1:a),  2+2+1, AFB(1:7,1:a), 4+2+1, IPV(1:a),EQUED, RRR(1:a), & 
             CCC(1:a), RHS(1:a), a, XXX(1:a,1),a, rcond, ferr,berr, work(1:(3*a)),iwork(1:a), info)
    IF(info.ne.0) THEN
        print*, 'ERROR: info = ', info, ' in DGBSVX, dynamic_sus_dist'
        stop
    END IF
    Cbar = XXX(1:a,1)*rhos
    !bed(1:a)=XXX(1:a,1)
    
    ! Sometimes, the near bed diffusive flux can lead to small negative values
    ! of Cbed being predicted right near the channel edge. This is because when
    ! the depth is changing by order of magnitude (e.g. from 1.0e-05 to 1.0e-02)
    ! outgoing flux from the small depth point [say (cbed*eddvis*slope)(i+1/2)
    ! ] can be much greater than the amount of suspended sediment over the small
    ! depth point. In my experience, this often happens at points with a tiny
    ! depth right next to the dry bank. (e.g. depth(i-1)=0.0)
    ! A simple fix is to make Cbed = 0.0 if negative, and remove the required
    ! amount of suspended sediment from the neighbouring point. 
    !tmp1 = maxval(Cbar)
    DO i = 1,a
        IF((Cbar(i)<0._dp).and.(.TRUE.)) THEN
            ! Clip negligably small Cbar values
            IF(abs(Cbar(i))<1.0e-10) THEN
                Cbar(i) = 0._dp
                goto 21121            
            END IF
            ! Look at neighbouring points
            IF((i>1).and.(i<a)) THEN
                IF(depth(i-1)>depth(i+1)) THEN
                    Cbar(i-1) =(Cbar(i-1)*depth(i-1)+Cbar(i)*depth(i))/depth(i-1)
                    Cbar(i) = 0._dp
                ELSE IF(depth(i+1)>depth(i-1)) THEN
                    Cbar(i+1) =(Cbar(i+1)*depth(i+1)+Cbar(i)*depth(i))/depth(i+1)
                    Cbar(i) = 0._dp
                ELSE IF((depth(i-1)==0._dp).and.(depth(i+1)==0._dp)) THEN
                    ! Isolated point
                    Cbar(i) = 0._dp
                END IF

            ELSE
                ! Treat boundary values
                IF(i==1) THEN
                    IF(depth(i+1)>0._dp) THEN
                        Cbar(i+1) =(Cbar(i+1)*depth(i+1)+Cbar(i)*depth(i))/depth(i+1)
                        Cbar(i) = 0._dp
                    ELSE
                        ! Isolated point
                        Cbar(i) = 0._dp
                    END IF 
                ELSE IF(i==a) THEN
                    IF(depth(i-1)>0._dp) THEN
                        Cbar(i-1) =(Cbar(i-1)*depth(i-1)+Cbar(i)*depth(i))/depth(i-1)
                        Cbar(i) = 0._dp
                    ELSE
                        Cbar(i)=0._dp
                    END IF
                END IF

            END IF
    21121 CONTINUE
        END IF
    END DO
    
    IF((.TRUE.).and.(minval(Cbar)<0._dp)) THEN
        print*, 'Cbar clip', minval(Cbar)
        Cbar = max(Cbar, 0._dp)
        IF(counter.eq.1) print*, 'WARNING: Negative Cbar values are being clipped &
                                to zero'
        ! IDEA -- try approaching this with an approach along the lines of the
        ! one in that paper that Steve Roberts pointed you to.
    END IF 

    ! New near bed cb
    DO i = 1, a
        IF(zetamult(i) > 0._dp) THEN
            cb(i) = Cbar(i)*depth(i)/zetamult(i)
        ELSE
            cb(i) = 0._dp
        END IF
    END DO
    

    ! Sanity Checks
    halt = .FALSE.
    DO i = 1,a
        IF((isnan(Cbar(i))).or.(Cbar(i)<0._dp)) THEN
            print*, 'Cbar(', i, ') is ', Cbar(i), cb(i), depth(i),zetamult(i) 
            halt=.TRUE.
        END IF
    END DO
    IF(halt) THEN
        DO i = 1,a 
            print*, Cbar(i), cb(i), depth(i)
        END DO
        stop
    END IF
END SUBROUTINE dynamic_sus_dist2
!!!!!!!!!!!!!!!!!!!!!!!!!
