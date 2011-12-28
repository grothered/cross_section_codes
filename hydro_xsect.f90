MODULE hydro_xsect ! Module for hydrodynamics
!
!
!

!Global constants
use global_defs
! MODULE TO HOLD PIZZUTO SHEAR MODEL
use pizzutotry

IMPLICIT NONE

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dAdP(nn, ys,bed, water, B, dAP,dPer,tbst, f,slopes, rho, lambdacon,tbston)
    !This just calculates some useful things for the shear routine.
    !When using the Pizzuto method, this routine is more complex, hence it is
    !separated from the main shear routine, although it is simple for the modified
    !version of the SKM. 

    ! nn = number of points
    ! ys = yvalues
    ! bed = bed values
    ! water = water elevation
    ! B = variable that stores the output
    ! dAP = defunct parameter (useful with a different shear routine)
    ! dPer = defunct parameter (useful with a different shear routine)
    ! tbst = defunct parameter (useful with a different shear routine)
    ! f = friction factor (in darcy-weisbach form, so tau = rho f/8 U^2)
    ! slopes = d(bed)/dy
    ! rho = density of water
    ! lambdacon = dimensionless eddy viscosity
    ! tbston = sqrt(1+(d(bed)/dy)^2 ), not needed here anymore

    INTEGER, INTENT(IN)::nn
    REAL(dp), INTENT(IN):: ys, bed,water,f,slopes, rho, lambdacon !So ys is a vector of y coordinates for the bed points, bed is their bed, water is the water depth, f is roughness
    REAL(dp),INTENT(IN OUT)::B, dAP, dPer,tbst ! the integral term, and the final derivative term (area centered around each point), and the wetted perimeter increment (b!etween each point, and the transverse bedslope term
    LOGICAL, INTENT(IN)::tbston

    DIMENSION:: B(0:nn),dAP(nn),ys(nn),bed(nn), dPer(nn),tbst(nn), f(nn),slopes(nn)

    REAL(dp):: lambda(0:nn) !, maxval
    REAL(dp):: tptx(nn) ! water surface point of each normal,x value only
    REAL(dp):: btx(nn) ! bottom of trapezium (shifted point on normal), x values only
    REAL(dp)::DA(nn) !Area
    INTEGER::i
    LOGICAL::flag

    !Dimensionless eddy viscosity
    lambda(1:nn) = lambdacon
    lambda(0)= lambda(1) 

    !!The function inside the derivative, pointwise
    B(1:nn)= rho*0.5_dp*lambda(1:nn)*(  ( water-bed(1:nn) )**2 )*sqrt((f(1:nn)/8._dp))
    B(0)= 0._dp


END SUBROUTINE dAdP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE shear(nn,ys,bed,water,wslope,taus,f,NNN,slopes, counter, Q, vegdrag, rho,g, lambdacon, tbston, & 
                ysl,ysu, bedl,bedu, high_order_shear)
    ! Solve a model for the distribution of bed shear and velocity over a
    ! cross-section

    ! Solves the equation:
    ! rho*g*Sf*depth = rho*f/8*U^2*sqrt(1+(dbed/dy)^2) + rho*vegdrag*depth*U^2
    !                 - d/dy(rho*0.5*lambdacon*depth^2*sqrt(f/8)*dU^2/dy) 
    ! for U^2 (and thus tau = rho*f/8*U^2)

    !nn is the number of bed points
    !ys is a vector of y coordinates for the bed points
    !bed is a vector of z coordinates for the bed points
    !water is the free surface elevation
    !wslope is the water slope. Note that even if it is incorrect, the shear
    !distribution will be correct to a constant multiple, because the equation
    !is linear in wslope. 
    !taus is the bed shear 
    !f is the bed roughness vector
    !NNN is a defunct parameter
    !slopes is a vector of the lateral bed slopes
    !counter tells us how many time steps we have taken
    !Q is the discharge
    !vegdrag is a vector of vegetation drag coefs
    !rho is the density of water
    !lambdacon is the dimensionless eddy viscosity. 
    !tbston switches on/off the term sqrt(1+slopes^2) in the SKM equation
    !ysl is the y value of the left-side point neighbouring the first wet point
    !ysu is the y value of the right-side point neighbouring the last wet point 
    !bedl is the z value of the left-side point neighbouring the first wet point
    !bedu is the z value of the right-side point neighbouring the last wet point 
    !high_order_shear (TRUE/FALSE) determines whether or not we try to use a higher
    !           order derivative approximation at interior points. 


    INTEGER, INTENT(IN)::nn, counter
    REAL(dp), INTENT(IN):: ys, bed,water,wslope, f,NNN,slopes, Q, vegdrag, rho, lambdacon, ysl, ysu, bedl, bedu,g 
    REAL(dp), INTENT(IN OUT):: taus !shear, side slope
    LOGICAL, INTENT(IN):: tbston, high_order_shear
    DIMENSION:: ys(nn),bed(nn),taus(nn), f(nn),NNN(nn), slopes(nn), vegdrag(nn)

    !logical, intrinsic:: isnan
    LOGICAL:: flag, stopper, const_mesh
    INTEGER:: info, i , j  !For the lapack routine
    REAL(dp)::Res, Res2=0._dp
    REAL(dp)::x(nn), l(nn),ds(nn),u(nn), tbst(nn)
    REAL(dp):: tol=1.E-8_dp, y0, ymax, power, Bconst, Aconst, b1, upart,constpart, ss, tmp,tmp1, dy_outer
    REAL(dp):: dAP(nn), B(0:(nn)),Dn(nn), dPer(nn), Bf(1:nn) 
    REAL(dp)::alpht(nn),alphb(nn),alpht2(nn), alphb2(nn), s(nn)  !first 2 are used in the matrix to represent the diagonals (top and bottom), and last is the right hand side
    REAL(dp):: diag(nn)  !the main diagonal of the matrix
    REAL(dp):: second(nn), vslopes(nn), K(nn) , Bderiv(nn), vel(nn), velsq(nn), dyf(nn)  
    REAL(dp):: DLF(nn), DF(nn), DUF(nn), DU2(nn), XXX(1:nn,1),rcond,ferr,berr,work(3*nn) 
    INTEGER:: IPV(nn),iwork(nn)
    CHARACTER(1)::EQUED
    REAL(dp):: bandmat(5,nn), AFB(7,nn), RRR(nn), CCC(nn)

    ! We calculate the shear with a guessed water slope
    ! value (we could use anything, though in practise we use a nearly correct value). 
    ! This gives us the correct shape of the tau distribution, but perhaps not the
    ! correct scale. We fix that by scaling later on.

    !Quick exit
    IF(wslope.eq.0._dp) THEN
        taus=0._dp
        RETURN
    END IF

    !Take care of the case in which there are not enough points to use the LDM
    IF (nn<=5) THEN 
        ! rho f/8 U^2 + rho*vegdrag*d *U^2 = rho g Sf d
        ! Solve for U^2, then multiply by rho f/8 to get taus
        taus(1:nn) = rho*f/8._dp*(rho*g*(water-bed(1:nn))*wslope)& 
                        /(rho*f/8._dp + rho*vegdrag(1:nn)*(water-bed(1:nn))  ) 
        tbst= 1._dp  ! A rough assumption. 
        RETURN
    END IF

    !Suppose the second derivative term in the equation is written as d/dy(B d(U^2)/dy).
    !The following subroutine evaluates B
    call dAdP(nn, ys,bed, water, B, dAP,dPer,tbst, f,slopes, rho, lambdacon, tbston)
    
    !Value of B at i+1/2
    Bf(1:nn-1)=0.5_dp*(B(2:nn)+B(1:nn-1))
    ! Forward dy increment
    dyf(1:nn-1)= ys(2:nn)-ys(1:nn-1) 

    ! Calculate length element
    IF(tbston) THEN
        tbst=sqrt( 1._dp + slopes*slopes )
    ELSE
        tbst=1._dp
    END IF
    
    ! The friction slope term, on the right hand side
    s(1:nn)= rho*(g*(water-bed)*wslope ) 

    !! Set up the diagonals. 
    !! Upper diagonal
    !alpht(2:nn-1) = - 1._dp/(.5_dp*(dyf(2:nn-1)+dyf(1:nn-2)))*( (Bf(2:nn-1))*1._dp/dyf(2:nn-1)) 
    !! Lower diagonal
    !alphb(2:nn-1) = - 1._dp/(.5_dp*(dyf(2:nn-1)+dyf(1:nn-2)))*((Bf(1:nn-2))*1._dp/dyf(1:nn-2)) 
    !! Main diagonal
    !diag(2:nn-1)=  rho*(f(2:nn-1)/8._dp)*tbst(2:nn-1)+rho*vegdrag(2:nn-1)*(water-bed(2:nn-1)) - alphb(2:nn-1) -alpht(2:nn-1) 

    ! Determine if the mesh has a constant increment -- this allows for easier
    ! high order methods
    const_mesh=.TRUE.
    tmp = ys(2)-ys(1)
    DO i = 3, nn
        IF( abs((ys(i)-ys(i-1))-tmp) > 1.0e-10_dp) THEN
            const_mesh=.FALSE.
            exit
        END IF
    END DO
    !const_mesh=.FALSE.

    diag =0.0_dp
    alpht=0.0_dp
    alphb=0.0_dp
    alpht2=0.0_dp
    alphb2=0.0_dp
    
    ! Define diagonal terms in the matrix which will approximate the differential equation
    ! diag = main diagonal
    ! alpht = first upper diagonal
    ! alpht2 = second upper diagonal
    ! alphb = first lower diagonal
    ! alphb2 = second lower  diagonal
    DO i = 2, nn-1
        !Non-derivative terms
        diag(i) = rho*(f(i)/8._dp)*tbst(i) + rho*vegdrag(i)*(water-bed(i))
        
        ! Derivative term
        ! d/dy( B dU^2/dy)
        dy_outer= 2._dp/(ys(i+1) - ys(i-1))
        !IF((i>1).and.(i<nn-1).and.(const_mesh).and.(high_order_shear)) THEN
        !    ! Use high order derivative
        !    tmp = dy_outer*1.0_dp/24.0_dp*Bf(i)*1._dp/dyf(i)
        !    tmp1 = dy_outer*9.0_dp/8.0_dp*Bf(i)*1._dp/dyf(i)

        !    alpht2(i) = alpht2(i) + tmp
        !    alpht(i) = alpht(i)   - tmp1           
        !    diag(i) = diag(i)     + tmp1
        !    alphb(i) = alphb(i)   - tmp
        !ELSE
            ! Use central derivative
            tmp = dy_outer*Bf(i)*1._dp/dyf(i)
            alpht(i) =alpht(i) -tmp
            diag(i) = diag(i)  +tmp
        !END IF

        !IF((i>2).and.(i<nn).and.(const_mesh).and.(high_order_shear)) THEN
        !    ! Use high order derivative
        !    tmp = dy_outer*1.0_dp/24.0_dp*Bf(i-1)*1._dp/dyf(i-1)
        !    tmp1 = dy_outer*9.0_dp/8.0_dp*Bf(i-1)*1._dp/dyf(i-1)
        !    alpht(i) = alpht(i)  -tmp
        !    diag(i) =  diag(i)   +tmp1
        !    alphb(i) = alphb(i)  -tmp1
        !    alphb2(i) = alphb2(i)+tmp
        !ELSE
            ! Use central derivative
            tmp = dy_outer*Bf(i-1)*1._dp/dyf(i-1)
            diag(i) = diag(i) +tmp
            alphb(i) = alphb(i) -tmp
        !END IF
     
    END DO

    !Boundary conditions
    !!!Assume that U^2 -> 0 as depth -> 0, and we be careful about the edge
    IF(.FALSE.) THEN
        alpht(1)=  0._dp - 1._dp/(.5_dp*(dyf(1)+(ys(1)-ysl)))*( .5_dp*(B(2)+B(1))*1._dp/dyf(1)) 
        alphb(1)= 0._dp
        diag(1)= rho*(f(1)/8._dp)*tbst(1) +rho*vegdrag(1)*(water-bed(1)) -alpht(1) + & 
                 1._dp/(.5_dp*(dyf(1)+(ys(1)-ysl)))*.5_dp*B(1)*1._dp/(ys(1)-ysl)
        !end if
        !if(.true.) THEN
        !
        !
        alpht(nn)=0._dp
        alphb(nn)= 0._dp - 1._dp/(.5_dp*(dyf(nn-1)+(ysu-ys(nn) ) ))*(.5_dp*(B(nn)+B(nn-1))*1._dp/dyf(nn-1))
        diag(nn) =  rho*(f(nn)/8._dp)*tbst(nn)+rho*vegdrag(nn)*(water-bed(nn))  -alphb(nn) +&
                    1._dp/(.5_dp*(dyf(nn-1)+(ysu-ys(nn) ) ))*.5_dp*B(nn)*1._dp/(ysu-ys(nn))
        !!!!!!!!
    END IF
   
    ! BOUNDARY CONDITIONS 
    ! Assume here that U^2 -> 0 as depth -> 0, and be careful to estimate the
    ! exact location of the waters edge. 
    IF(.TRUE.) THEN
        ! Left boundary
        alpht(1)=  0._dp - 1._dp/(.5_dp*(dyf(1)+(ys(1)-ysl)*(bed(1)-water)/(bed(1)-bedl)))*( (Bf(1))*1._dp/dyf(1)) 
        alphb(1)= 0._dp
        diag(1)= rho*(f(1)/8._dp)*tbst(1) +rho*vegdrag(1)*(water-bed(1)) -alpht(1) + & 
            1._dp/(.5_dp*(dyf(1)+(ys(1)-ysl)*(bed(1)-water)/(bed(1)-bedl)))*.5_dp*B(1)*& 
            1._dp/( (ys(1)-ysl)*(bed(1)-water)/(bed(1)-bedl))
       
        ! Right boundary 
        alpht(nn)=0._dp
        alphb(nn)= 0._dp - 1._dp/(.5_dp*(dyf(nn-1)+(ysu-ys(nn) )*(water-bed(nn))/(bedu-bed(nn)) ))*((Bf(nn-1))*1._dp/dyf(nn-1))
        diag(nn) =  rho*(f(nn)/8._dp)*tbst(nn)+rho*vegdrag(nn)*(water-bed(nn))  -alphb(nn) +&
            1._dp/(.5_dp*(dyf(nn-1)+(ysu-ys(nn) )*(water-bed(nn))/(bedu-bed(nn)) ))*.5_dp*B(nn)*1._dp/& 
            ((ysu-ys(nn))*(water-bed(nn))/(bedu-bed(nn)))
        !!!!!!!!
    END IF
    
    ! BOUNDARY CONDITIONS 
    ! Assume here that the momentum flux at (i=1/2) is 0 (so the system is not
    ! losing any momentum), and same for (i=nn-1/2)
    ! FIXME:(10/7/11)  A problem with this boundary condition has always been
    ! that it prevents the channel reaching  a steady state in the case of pure
    ! suspension of constant concentration -- so it is a bit of a worry. The
    ! extra momentum loss in the above condition is small as long as the depth
    ! at the channel edge is small.
    IF(.FALSE.) THEN
        ! Left boundary
        !alpht(1)=  0._dp - 1._dp/(.5_dp*(dyf(1)+(ys(1)-ysl)*(bed(1)-water)/(bed(1)-bedl)))*( (Bf(1))*1._dp/dyf(1)) 
        alpht(1)=  0._dp - 1._dp/(.5_dp*(dyf(1)+(ys(1)-ysl)))*( (Bf(1))*1._dp/dyf(1)) 
        alphb(1)= 0._dp
        diag(1)= rho*(f(1)/8._dp)*tbst(1) +rho*vegdrag(1)*(water-bed(1)) -alpht(1)  
       
        ! Right boundary 
        alpht(nn)=0._dp
        !alphb(nn)= 0._dp - 1._dp/(.5_dp*(dyf(nn-1)+(ysu-ys(nn) )*(water-bed(nn))/(bedu-bed(nn)) ))*((Bf(nn-1))*1._dp/dyf(nn-1))
        alphb(nn)= 0._dp - 1._dp/(.5_dp*(dyf(nn-1)+(ysu-ys(nn) ) ))*((Bf(nn-1))*1._dp/dyf(nn-1))
        diag(nn) =  rho*(f(nn)/8._dp)*tbst(nn)+rho*vegdrag(nn)*(water-bed(nn))  -alphb(nn) 
        !!!!!!!!
    END IF

    ! Fill out matrix for LAPACK solver
    bandmat=0._dp
    DO i=1,nn
        ! Make sure that zero depth = zero shear
        IF(s(i).eq.0._dp) THEN
            diag(i) = 1._dp
            alpht(i)=0._dp
            alphb(i)=0._dp
            alpht2(i)=0._dp
            alphb2(i)=0._dp
        END IF
        
        IF(i<nn-1) bandmat(1,i+2) = alpht2(i)
        IF(i<nn) bandmat(2,i+1) = alpht(i)
        bandmat(3,i)   = diag(i)
        IF(i>1) bandmat(4,i-1) = alphb(i)
        IF(i>2) bandmat(5,i-2) = alphb2(i)
    END DO

    ! Set up the matrix solution. 

    ! Matrix solver -- use banded matrix solver, or if possible, tridiagonal
    ! solver
    IF((const_mesh).AND.(high_order_shear)) THEN
        ! Banded matrix solver, LAPACK
        call DGBSVX('E','N', nn,2,2,1, bandmat(1:5,1:nn),  2+2+1, AFB(1:7,1:nn), 4+2+1, IPV(1:nn),EQUED, RRR(1:nn), & 
                 CCC(1:nn), s(1:nn), nn, XXX(1:nn,1),nn, rcond, ferr,berr, work(1:(3*nn)),iwork(1:nn), info)
    ELSE
        ! Tridiagonal matrix solver, LAPACK
        XXX(1:nn,1) = s(1:nn)
        call DGTSV(nn, 1, bandmat(4,1:nn-1), bandmat(3,1:nn), bandmat(2,2:nn), XXX(1:nn,1), nn, info)
    END IF
    ! Determine tau
    taus = XXX(1:nn,1)*rho*(f/8._dp)

    !Check
    IF(info.ne.0) THEN
        IF((const_mesh).AND.(high_order_shear)) THEN
             PRINT*, 'ERROR: info = ', info, ' in DGBSVX, in shear'
        ELSE
             PRINT*, 'ERROR: info = ', info, ' in DGTSV, in shear'
        END IF
        stop
    END IF
        

    ! Detailed testing, doesn't work with newer matrix solvers
    IF(.FALSE.) THEN
        !do i=1,nn
        !        if(i>1) bandmat(1+1+i-(i-1),i-1)=alphb(i)
        !        bandmat(1+1+i-(i),i)=diag(i)
        !        if(i<nn) bandmat(1+1+i-(i+1),i+1)=alpht(i)
        !end do
        !call DGBSVX('E','N', nn,1,1,1, bandmat(1:3,1:nn),  1+1+1, AFB(1:4,1:nn), 2+1+1, IPV,EQUED, RRR(1:nn),  CCC(1:nn), velsq, nn, &
        !        XXX(1:nn,1),nn, rcond, ferr,berr, work(1:(3*nn)),iwork(1:nn), info)      
        !if(info.ne.0) THEN
        !print*, 'Tau solver DGBSVX ERROR, info=', info, 'rcond = ', rcond
        !end if
        !print*, 'EQUED = ', EQUED
        !taus= XXX(1:nn,1)*rho*(f/8._dp)

        ! Check for errors
        ! if(taus(nn)<0._dp) print*, 'negtaus ', taus(nn),wslope 
        IF(info.NE.0) PRINT*, "info .ne. 0", info

        DO i = 1, nn
            flag= (isnan(taus(i)))
            IF(flag) THEN
                   ! taus(i)=0._dp
                print*, "Tau is nan " , taus(i),i, maxval(s), minval(s), nn, maxval(slopes), & 
                        minval(slopes), maxval(water-bed), minval(water-bed) 
                        !, minval(water-bed) ! , i, alphb(i:i+3),"---", alpht(i:i+3), "---", s(i:i+3)/tbst(i:i+3),  "---", diag(i:i+3)
                print*, info      
                print*, "----"
                !  print*, y0, ys(1), ys(2), "....", ys(nn-1), ys(nn), ymax
                    !   print*, "caught tau nan just after dgtsv"!alpht 
                !stop
            END IF
        END DO

        !!!!!!!!Error checking

        x=taus/(rho*f/8._dp)

        stopper=.false.
        Res = s(1)-diag(1)*x(1)-alpht(1)*x(2)
        Res=Res/diag(1)
        IF(abs(Res).GT.tol) PRINT*, " thomas: alarm", Res !, s(1)
        DO i=2,nn-1
            Res = s(i)-alphb(i)*x(i-1)-diag(i)*x(i)-alpht(i)*x(i+1)
            IF (abs(Res).gt.(tol*diag(i))) THEN
                !Res2=max(abs(Res)/tbst(i),Res2)
                PRINT*, " Matrix solution : alarm", Res/tbst(i), s(i), diag(i),taus(i), i ,nn, rcond
                stopper=.true.
            END IF
        END DO     
        IF(stopper) STOP
    END IF
    !!After all the routines are finished, we correct the tau distribution so that the discharge is correct. We don't do this here, because it is possible that internal dry points have partitioned the cross section, so we don't know what the discharge is.

END SUBROUTINE shear
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_shear(a, dT, water, Q, bed,ys,Area, bottom, ff,rmu,inuc,tau,NN,counter &
            ,slopes, hlim,nos, vegdrag, rho, & 
            rhos, voidf, d50, g, kvis, vertical, lambdacon, tbston, ysl,ysu,bedl,bedu, & 
            high_order_shear) 
    ! A routine which calculates the shear over the cross-section
    ! It uses the routine 'shear' above.
    ! The routine ensures that the discharge is correct,
    ! allows the suitable treatment of interior dry points in the model domain,
    ! calculates the 1D friction factor,
    ! and calculates a non-uniform convective term for interaction with a 1D
    ! unsteady solver

    ! a = number of points on the cross-section
    ! dT = the time step
    ! water = the water surface elevation
    ! Q = discharge
    ! bed = the bed elevation
    ! ys = the y coordinates
    ! Area = the cross-sectional area
    ! bottom = the average bed elevation
    ! ff = the darcy-weisbach friction factor
    ! rmu = the 'roughness multiplier' -- defined so that the friction slope Sf = rmu * (Q/Area)^2
    !                                  -- so in some sense rmu = (f/8)/mean_depth
    ! inuc = non-uniform convective intertia multiplier term for interaction with a 1D unsteady solver
    ! tau = bed shear
    ! NN = defunct term that I had been using to add secondary flows to the shear model
    ! counter = variable which records the loop iterator in the routine that calls this
    ! slopes = dbed/dy
    ! hlim = do not compute the shear if the average depth is less than hlim
    ! nos = number of points in taucrit_dep_ys
    ! vegdrag = vegetation drag coefficient
    ! rho = density of water
    ! Qe = rate of resuspension
    ! Qbed = bedload flux 
    ! voidf = void fraction of the bed (porosity)
    ! d50 = median sediment size
    ! g = gravity
    ! kvis = kinematic viscosity
    ! vertical = logical -- do we use the 'shear' routine, or one based on Pizzuto (1991)
    ! lambdacon = dimensionless eddy viscosity
    ! tbston = logical -- do we retain the sqrt(1+(dbed/dy)^2) term in the shear routine (true) or treat it as 1 (false)
    ! ysl = y-coordinate just left of the wetted part of the cross-section
    ! ysu = y-coordinate just right of the wetted part of the cross-section
    ! bedl = bed-coordinate just left of the wetted part of the cross-section
    ! bedu = bed-coordinate just right of the wetted part of the cross-section
    ! high_order_shear = logical -- do we try to use a high order approximation for
    ! interior derivatives in the 'shear' routine


    INTEGER, INTENT(IN)::a,counter,nos
    REAL(dp), INTENT(IN)::water,Q, Area, bottom, ff, hlim, vegdrag, dt, rho, rhos, voidf,&
         d50, g, kvis, lambdacon, ysl,ysu,bedl,bedu, bed, ys,slopes
    REAL(dp), INTENT(IN OUT):: rmu,inuc,tau, NN 
    LOGICAL, INTENT(IN):: vertical, tbston, high_order_shear
    DIMENSION bed(a),ys(a), ff(a),tau(a), NN(a),slopes(a), vegdrag(a)! 
    
    INTEGER::i, j, bgwet, up,  jj,  info,ii, n(a)
    REAL(dp)::wslope,  Qelocal, tt, corfact
    REAL(dp):: kkkk(a), tbst(a), f(a)  
    REAL(dp)::vel(a), Qb(a),& 
        bedlast(a), sinsl(a), mu_d, Qtemp, useful(a), Ceq(a) 
    !logical::writ_tau=.TRUE.  !This will write out the cross sectional taus-- will lead to massive files if you're not careful. 
    LOGICAL::  dry(a)

    ! Water slope term. 
    IF(vertical) THEN
        !Water slope, accounting for the (shape factor)/(depth)= rmu.
        wslope= ((Q/Area)*abs(Q/Area)/((1._dp))*rmu) 
    ELSE 
        ! at present we have no algorithm for rmu for the normal depth method.
        wslope=  ((Q/Area)*abs(Q/Area)/(((water-bottom)))*sum(ff)/(8._dp*9.8_dp*a)) 
    END IF

    ! DEBUG
    IF(mod(counter,1000).eq.0) THEN 
        print*, wslope, rmu, ff(a/2), bedl, bedu, ysl, ysu
    END IF

    IF(abs(wslope)>.1) print*, "|wslope| >.1", wslope, water-bottom, rmu, Q, Area, Q/Area
    IF(isnan(wslope)) print*, "wslope is nan", Q, A, water-bottom, rmu

    IF(a<3) THEN
        print*, 'ERROR: less than 3 wet points on a cross-section -- shear &
                 cannot work with this', a
        stop
    END IF


    ! CALCULATE SHEAR
    ! This can deal with multiple wet sections/ interior dry points. 
    ! It also has an option as to how we calculate the shear

    ! Predefine
    tau(1:a)=0._dp
    tbst(1:a)=1._dp
    kkkk= 0.03_dp
    !kkkk=((water-bed)/cos(atan(slopes)))/(exp(1._dp+.4_dp*sqrt(8._dp/ff))) !For a justification of this formula, see my notes on the friction factor.
    DO i= 1, a
        IF(isnan(bed(i))) THEN
            print*, "bed(i) nan before shear", bed(i-1:i+1)
        END IF
    END DO

    ! If the bed is too shallow, skip the shear calculation
    IF(water-bottom<hlim) THEN
        vel=0._dp
        GOTO 131  !This will jump past the shear calculation for very shallow channels
    END IF

    ! Now we do the shear calculations. The way that we do it allows for
    ! a section with mid-channel bars. Note that we calculate a whole bunch of shears
    ! for an arbitary friction slope, and then correct it so that the integrated
    ! discharge matches the total discharge at the end.
    
    dry= (water-bed<0.0_dp) ! I have ensured that bed(1) and bed(a) are wet, but it is still possible that there are interior dry points
    bgwet= 1 ! This will record the location of the leftmost point which is connected to i by water  
    IF(vertical) THEN
            !if(counter>5011) print*, "we get to the shear loop"
        DO i=2, a

            IF ( ( (dry(i)).AND.(.NOT.dry(i-1)) ) )  THEN ! Reached an interior dry point folling a wet point
                up=i-1
                IF(maxval(water-bed(bgwet:up))>.01_dp) THEN
                    call shear(up-bgwet+1,ys(bgwet:up),bed(bgwet:up),water, wslope ,tau(bgwet:up), & 
                            ff(bgwet:up),NN(bgwet:up), slopes(bgwet:up), counter, Q, vegdrag(bgwet:up), rho,g, & 
                            lambdacon, tbston, ysl,ysu, bedl, bedu, high_order_shear)
                END IF
            END IF

            IF ( (.NOT.dry(i)).AND.(dry(i-1))) THEN !So we reached a new wet boundary following an interior dry point 
                bgwet= i
            END IF

            IF((i==a)) THEN 
                IF (maxval(water-bed(bgwet:a))>.01_dp ) THEN 
                    call shear(i-bgwet+1,ys(bgwet:a),bed(bgwet:a),water, wslope ,tau(bgwet:a),ff(bgwet:a) &
                        ,NN(bgwet:a),slopes(bgwet:a),counter, Q, vegdrag(bgwet:a), rho,g, lambdacon, tbston, ysl,ysu, &
                         bedl, bedu, high_order_shear)
                END IF
            END IF

        END DO

        ! Make sure integrated discharge is correct 
        vel= sqrt(8._dp*abs(tau)/(ff*rho))*sign(1._dp+0._dp*tau, tau) 
        ! vel = velocity, if the tau distribution is correct. 
        ! In general there could be a scale error, which we correct below using corfact
        ! to ensure that INTEGRAL ( vel*(water-bed) ) dy = Q 
        IF((abs(Q)>0._dp).AND.(maxval(abs(vel))>0._dp)) THEN
            !if(counter>5011) print*, "we get to the correction of tau", maxval(vel),	
            corfact=(Q/sum(0.5_dp*( vel*(water-bed) )*(& 
                (/1._dp*(ys(2)-ys(1) +1._dp*((water-bed(1))/(bedl-bed(1))*(ys(1)-ysl) )),&
                (ys(3:a)-ys(1:(a-2))), & 
                (1._dp*((water-bed(a))/(bedu-bed(a))*(ysu-ys(a)) )+ ys(a)-ys(a-1))/))))**2
            !corfact=1._dp
            tau = tau*corfact
            vel=vel*sqrt(corfact)
            ! DEBUG
            IF(mod(counter,1000).eq.0) THEN
               corfact= sum(0.5_dp*( vel*(water-bed) )*(& 
                            (/1._dp*(ys(2)-ys(1) +1._dp*((water-bed(1))/(bedl-bed(1))*(ys(1)-ysl) )),&
                            (ys(3:a)-ys(1:(a-2))), & 
                            (1._dp*((water-bed(a))/(bedu-bed(a))*(ysu-ys(a)) )+ ys(a)-ys(a-1))/)))
               print*, 'Corfact:', corfact, water-bed, vel, ys, (water-bed(a))/(bedu-bed(a))*(ysu-ys(a)) 
            END IF
        ELSE
            tau=0._dp
            vel=0._dp
        END IF
    END IF !END OF THE VERTICAL METHOD

    ! Normal shear method of Pizzuto
    IF(vertical.eqv..false.) THEN
        PRINT*, "the normal shear method is not supported in this version of the code - where we allow a signed tau"
        STOP
        useful=0._dp

        DO i=2, a

            IF ( ( (dry(i)).AND.(.NOT.dry(i-1)) ) )  THEN !So we reached an interior dry point folling a wet point, or we got to the end (which we know is wet) 
                up=i-1
                IF(maxval(water-bed(bgwet:up))>.01_dp) THEN
                    call shearP(up-bgwet+1,ys(bgwet:up),bed(bgwet:up),water, wslope ,tau(bgwet:up),kkkk(bgwet:up),&
                        ff(bgwet:up),NN(bgwet:up), slopes(bgwet:up), counter, Q, vegdrag(bgwet:up), rho, useful(bgwet:up))
                END IF
            END IF

            IF ( (.NOT.dry(i)).AND.(dry(i-1))) THEN !So we reached a new wet boundary following an interior dry point 
            bgwet= i
            END IF

            IF((i==a)) THEN 
                IF (maxval(water-bed(bgwet:a))>.01_dp ) THEN 
                    call shearP(i-bgwet+1,ys(bgwet:a),bed(bgwet:a),water, wslope ,tau(bgwet:a),kkkk(bgwet:a),ff(bgwet:a) &
                        ,NN(bgwet:a),slopes(bgwet:a),counter, Q, vegdrag(bgwet:a), rho, useful(bgwet:a))
                END IF
            END IF
        END DO

        ! Now correct the shear distribution
        ! to ensure that the integrated velocity is equal to the discharge 
        IF(maxval(useful)>0._dp) THEN 
            Qtemp= sum( (sqrt(abs(tau)/rho)/.4_dp)*useful )  ! 'useful' will include everything we need
            IF(abs(Q)>0._dp) THEN
                tau = tau*(Q/Qtemp)**2
            ELSE
                tau=0._dp
            END IF
        END IF

        !End of the corrected shear distribution
    END IF !END OF THE NORMAL METHOD


    ! Check that the results are reasonable!
    IF(maxval(abs(tau))>500._dp) THEN
        PRINT*, 'taumax>500', maxval(abs(tau)), rmu, wslope, vel(maxloc(abs(tau))), & 
            rho*ff(maxloc(abs(tau)))/8._dp, maxloc(abs(tau)), a, water, bed(maxloc(abs(tau))) 
            !, counter, a, Q, Qtemp, Area, elev, maxval(tau), minval(tau), maxval(useful), minval(useful)
            !elev-bed !"taumax>10", Q, maxval(vel), minval(ff), maxval(elev-bed), a, maxval(tau), minval(tau)
            !'vel=', vel, 'a=', a, elev-bed!maxval(s), minval(s), maxval(elev-bed)
        !stop
    END IF

    131 DO i= 1, a  !The label allows dry points to jump past the shear calculations
            IF((isnan(tau(i)) )) THEN
                print*, "tau is less than 0. ", a, maxval(tau), minval(tau)
                print*, water
                !tau!minval(tau)
                print*, (water-bed)
               ! write(2,*) bed
                !tau(i)=0.
              !  print*, tau
                print*, "____"
               ! print*, dry
                print*, "----"
                
                print*, "caught tau nan in sectupdate" !elev-bed
                !print*, bgwet, size(tau)
                stop

            END IF
        END DO
    ! End of check


    ! Calculate new value of rmu = correction factor for shear
    IF(tbston) THEN !When integrating the friction factor, include the sqrt(1+slopes^2) term
        call roughmult(a,rmu, vel, Q, Area, sqrt(1._dp+slopes**2),& 
            max(water-bed,0._dp), ys, ff, vegdrag, counter, ysl, ysu, bedl, bedu, bed, water,g)
    ELSE
        call roughmult(a,rmu, vel, Q, Area, 1._dp +0._dp*slopes ,& 
            max(water-bed,0._dp), ys, ff, vegdrag, counter,ysl, ysu, bedl, bedu, bed, water,g)
    END IF

    ! Calculate the integral of the nuc term.
    call intnuc(inuc,water,ys,bed,vel, a, Q, Area)


    IF(isnan(rmu)) THEN
        print*, "rmu is nan"
        print*, water-bed
        print*, "______"
        print*, water
        !stop
    END IF   
           

END SUBROUTINE calc_shear
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE roughmult(aa,rmu, vel, Q, A, tbst,depths, ys, f, vegdrag, &
                     counter, ysl, ysu, bedl, bedu,bed,water,g)
    ! Compute:
    ! integral[ Bed + vegatation shear ] / (1D_velocity^2 *g *width*av_depth)
    ! Note that this is some sort of averaged version of:
    ! total friction / av_depth

    ! Purpose: for calculating the constant that we should multiply the friction
    ! slope in the 1D model by, in order to account for the lateral distribution of
    ! velocity. 
    ! So it basicaly gives a roughness coefficient, although actually,
    ! we get the (roughness coefficient/ depth) for numerical stability reasons

    INTEGER, INTENT(IN)::aa, counter
    REAL(dp), INTENT(IN OUT):: rmu
    REAL(dp), INTENT(IN):: vel, Q, A, tbst, depths, ys, f, vegdrag, ysl, ysu,&
                           bedl, bedu, bed, water,g
    DIMENSION vel(aa),tbst(aa), depths(aa), ys(aa), f(aa), vegdrag(aa), bed(aa)
    
    INTEGER:: a2, i
    REAL(dp):: depths2(aa), rmutop, rmubot, rmulasts

    !Calculation requires |velocity| > 0
    IF((abs(Q)>1.0E-10_dp).AND.(maxval(abs(vel))>1.0E-10_dp)) THEN
        
        rmutop=0._dp
        ! Numerator integrated using trapezoidal rule
        DO i=2,aa-1
            rmutop= rmutop+ & 
                    (f(i)/8._dp*vel(i)**2*tbst(i) + vegdrag(i)*vel(i)**2*depths(i))*0.5_dp*&
                    (ys(i+1)-ys(i-1))
        END DO
        rmutop=rmutop+ &
               (f(1)/8._dp*vel(1)**2*tbst(1) + vegdrag(1)*vel(1)**2*depths(1))*0.5_dp*&
               (ys(2)-ys(1)+2._dp*(water-bed(1))/(bedl-bed(1))*(ys(1)-ysl))
        rmutop=rmutop+ &
               (f(aa)/8._dp*vel(aa)**2*tbst(aa) + vegdrag(aa)*vel(aa)**2*depths(aa))*0.5_dp*&
               (ys(aa)-ys(aa-1) +2._dp*(water-bed(aa))/(bedu-bed(aa))*(ysu-ys(aa)))

        !Denominator 
        rmubot= (Q/A)**2*A*g  
        ! (roughness/depth)
        rmu= rmutop/rmubot
        ! DEBUG
        IF(mod(counter,1000).eq.0) THEN
            print*, 'roughmult', rmu, rmutop, rmubot, sum(vel)/(1.0_dp*aa), Q, A
        END IF
    ELSE
        ! Default behaviour if there is no discharge
        !rmu=sum(f)/(8._dp*g*sum(depths)) !The mean value of (f/8g)/depth
        rmutop=0.0_dp
        DO i=2,aa-1
            rmutop= rmutop + &
                    (f(i)/8.0_dp*tbst(i) + vegdrag(i)*depths(i))*0.5_dp*(ys(i+1)-ys(i-1))
        END DO
        !Denominator
        rmubot=A*g
        ! (roughness/depth)
        rmu=rmutop/rmubot 
    END IF

    ! Sanity checks
    IF(rmu.eq.0._dp) THEN
        print*, " ERROR: rmu is zero", abs(Q), aa, sum(f), maxval(vel), minval(vel)
        stop
    END IF

    IF(isnan(rmu)) THEN
            print*, "rmu is nan"
            print*, vel
            print*, "...."
            print*, tbst
            print*, "..."
            print*, ys
            stop
    END IF

END SUBROUTINE roughmult
!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_friction(friction_type, grain_friction_type, rough_coef, water,&
                            a, bed, vel, man_nveg,d50,veg_ht, rhos, rho, g, f,&
                            vegdrag,f_g, dsand, counter, a_ref) 
    INTEGER, INTENT(IN):: a, counter
    CHARACTER(char_len), INTENT(IN):: friction_type, grain_friction_type
    REAL(dp), INTENT(IN):: rough_coef, water, man_nveg, d50, veg_ht, rhos, rho, g, dsand
    REAL(dp), INTENT(IN):: bed, vel
    REAL(dp), INTENT(IN OUT):: f, vegdrag, f_g, a_ref

    DIMENSION bed(a), vel(a), f(a), vegdrag(a), f_g(a), a_ref(a)
    
    !INPUT: Various friction related parameters
    !PURPOSE: calculate the value of the friction parameters 'f' and 'vegdrag'
    !and 'f_g' (which is the grain related roughness)
    
    !LOCAL VARIABLES
    REAL(dp)::ks, f_tmp(a), Re(a), f_glast(a)
    REAL(dp):: dgravel, dsilt, f_cs, si, k_scr, k_scmr, k_scd, f_fs, k_sc, k_sg
    INTEGER:: i, j
    REAL(dp):: onethird = (1.0_dp/3.0_dp)


    ! Friction factors
    SELECT CASE(friction_type)

        CASE('manning')
            ! Manning style friction
            DO i= 1, a
                !f(i)= rough_coef**2*g*8._dp/(max( water-bed(i),200._dp*d50)**(onethird))!
                f(i)= rough_coef**2*g*8._dp/(max( water-bed(i),1.0e-08_dp)**(onethird))!
            END DO 

            ! Now spatially average the friction
            !DO i = 2, a-1
            !    f_tmp(i) = 0.25_dp*(2.0_dp*f(i) + f(i+1) + f(i-1)) ! Average value of f in a cell centred at i 
            !END DO
            !    f_tmp(1) = 0.25_dp*(3.0_dp*f(1) + f(1+1))
            !    f_tmp(a) = 0.25_dp*(3.0_dp*f(a) + f(a-1))
           
            !    f = f_tmp 

        CASE('darcy')
            !Darcy-weisbach style friction factor
            f= rough_coef 

        CASE('ks') 
            ! Roughness height style friction factor
            ks = rough_coef

            DO i=1,a
            !f(i) = 8._dp*g/(18.0_dp*log10(12.0_dp*max((water-bed(i)),ks)/ks))**2
            !f(i) = 8._dp*g/(18.0_dp*log10(12.0_dp*max((water-bed(i)),ks)/ks))**2
            f(i) = 8.0_dp*(0.4_dp/log(max(water-bed(i), 3.0_dp*ks)/ks - 1.0_dp))**2
            END DO

        CASE('vanrijn') 
            ! van Rijn (2007) friction -- note that I am concerned that there could
            ! easily be typos in this paper (or in my coding of it!) 
            ! First need to calculate roughness height
            dgravel=0.002_dp !According to van Rijn -- yes this is 2mm, not a typo
            dsilt = 0.000032_dp !According to van Rijn -- yes this is 32 micro m, not a typo

            !Useful fudge constants f_cs, f_fs
            f_cs = min( (0.25_dp*dgravel/ d50)**(1.5_dp), 1.0_dp)
            f_fs = min(d50/(1.5_dp*dsand),1.0_dp) 

            DO i = 1, a
                si = (vel(i)**2)/((rhos-rho)/rho*g*d50)

                ! Roughness height due to ripples, eqn 5e
                k_scr = f_cs*d50*(85.0_dp-65.0_dp*tanh(0.015_dp*(si - 150.0_dp)))
                k_scr = max(k_scr, 20.0_dp*dsilt) !Lower limit according to text


                !Roughness height due to mega_ripples - eqn 6e
                !This seems to have a discontinuity, so I cut it for now.
                PRINT*, 'ERROR: Mega-ripples roughness not yet implemented in vanrijn &
                            roughness, in calc_friction '
                stop
                k_scmr=0.0_dp
                !k_scmr = 2.0e-05_dp*f_fs*(water-bed(i))*(1.0_dp - exp(-0.05_dp*si))&
                !            *(550._dp-si)
                !!6c, 6d, and the next one
                !IF((si>550.0_dp).and.(d50>=1.5_dp*dsand)) k_scmr = 0.02_dp
                !IF((si>550.0_dp).and.(d50<1.5_dp*dsand)) k_scmr = 200._dp*d50
                !IF(d50<dsilt) k_scmr = 0._dp
                !! In the text, states that k_scmr <= 0.2
                !k_scmr = min(k_scmr, 0.2_dp)

                !Roughness height due to dunes, eqn 7e -- note van Rijn typo in
                !variable name
                k_scd = 8.0e-05_dp*f_fs*(water-bed(i))*(1.0_dp - exp(-0.02_dp*si))*(600._dp-si)            
                !Eqns 7c, &7d
                IF((si>600._dp).or.(d50<dsilt))THEN
                    k_scd=0._dp
                END IF
                !In the text, states that k_scd<1.0m
                k_scd = min(k_scd, 1.0)

                ! Total physical current roughness
                k_sc = sqrt(k_scr**2 + k_scmr**2 + k_scd**2)        

                ! Friction factor
                f(i) =(8._dp*g/(18._dp*log10(12._dp*max(water-bed(i), k_sc)/(k_sc)))**2)  

                IF((mod(counter,10000).eq.0).and.(i.eq.(a/2))) THEN
                    print*, 'Roughnesses in channel centre:', k_scr, k_scmr, k_scd, k_sc, f(i), & 
                        'depth', water-bed(i), 'si', (vel(i)**2)/((rhos-rho)/rho*g*d50)
                END IF

            END DO
    
        CASE DEFAULT
            PRINT*, 'ERROR -- friction_type does not appear to be correctly specified in calc_shear'
            stop
        
    END SELECT

    ! Vegetation drag
    WHERE (bed > veg_ht)
            vegdrag= man_nveg**2*g*8._dp 
    ELSEWHERE
            vegdrag=0._dp
    END WHERE
        

    ! Grain roughness
    SELECT CASE (grain_friction_type)
        CASE('vanrijn')
            ! Van Rijn, fully turbulent flow
             ! NOTE THAT vanrijn writes that tau = 0.5*f*rho*U^2 -- however, his
             ! data in table2 of the paper are better predicted using the
             ! 'normal' formula, tau = rho f/8 U^2 --- I think the paper must
             ! have a typo
             f_g =(8._dp*g/(18._dp*log10(12._dp*max(water-bed, 20.0_dp*3.0_dp*d50)/(3._dp*d50)+0.0_dp)+0.0e+00_dp)**2)
             
             ! Note that this version (log instead of log10) gave me nice
             ! results in an early version of the code --- but it is not
             ! correct!
             !f_g =(8._dp*g/(18._dp*log(12._dp*max(water-bed, 20.0_dp*3.0_dp*d50)/(3._dp*d50)+0.0_dp)+0.0e+00_dp)**2)

        CASE('colebrook') 
            ! Colebrook and White, following Chanson (2004)
            ! f = 0.25/(log10(ks/(3.71*4*d) + 2.51/(Re*sqrt(f))))^2
            ! where Re = u_star*ks/kinematic_visc 
            k_sg = 10.0_dp*d50
            !!Re = max((sqrt(f_g/8.0_dp)*abs(vel))*k_sg/1.0e-06_dp, 10.0_dp)
            Re = max(abs(vel)*max(water-bed, 20.0_dp*k_sg)/1.0e-06_dp, 10.0_dp)
            !! Solve through iteration
            f_glast=f_g+0.001_dp
            DO i=1,a
                j=0
                DO WHILE(abs(f_glast(i) - f_g(i))>1.0e-08_dp)
                    j=j+1
                    f_glast(i)=f_g(i)
                    f_g(i) = 0.25_dp/(log10( k_sg/(3.71_dp*4.0_dp*max(water-bed(i),20._dp*k_sg)) &
                                  + 2.51_dp/(Re(i)*sqrt(f_glast(i))))  )**2
                    IF(j>1000) THEN
                        PRINT*,'ERROR in CALC_FRICTION: f_g did not converge in ', j,' iterations' 
                        stop
                    END IF
                END DO
            END DO

        CASE('onethird')
            f_g = f/3._dp

        CASE('one')
            f_g = f
        
        CASE DEFAULT
            PRINT*, 'ERROR: grain_friction_type does not appear to be specified correctly in calc_shear'
            stop

    END SELECT

    ! The overall roughness should be >= the grain roughness.
    f = max(f,f_g)

    ! Compute vanrijn's (2007) roughness height 'a'
    dgravel=0.002_dp !According to van Rijn -- yes this is 2mm, not a typo
    IF(d50< 0.25_dp*dgravel) THEN
        f_cs=1.0_dp
    ELSE
        f_cs = (0.25_dp*dgravel/ d50)**(1.5_dp)
    END IF
    
    DO i=1,a 
        si = (vel(i)**2)/((rhos-rho)/rho*g*d50)
        k_scr = f_cs*d50*(85.0_dp-65.0_dp*tanh(0.015_dp*(si - 150.0_dp)))
        !a_ref = 0.5_dp*k_scr !, 0.01_dp*(water-bed(i))) !Reference level (m) 
        !a_ref =  max(0.01_dp, 0.5_dp*k_scr)!, 0.99_dp*(water-bed(i))) !, 0.01_dp*(water-bed(i))) !Reference level (m) 
        !a_ref = min(max(0.01_dp, 0.5_dp*k_scr), 0.99_dp*(water-bed(i))) !, 0.01_dp*(water-bed(i))) !Reference level (m) 
        a_ref(i) = max(0.5_dp*k_scr, 0.01_dp*(water-bed(i))) !Reference level (m), van Rijn 1984
        !a_ref(i) = max(0.5_dp*k_scr, 0.01_dp) !Reference level (m), van Rijn 2007
        !a_ref(i) = max(0.5_dp*k_scr, 0.01_dp*(water-bed(i)), 0.01_dp) !Reference level (m), van Rijn 2007
        !a_ref(i) = 0.05_dp*(water-bed(i)) ! Garcia and Parker (1991), and many others
    END DO
    

END SUBROUTINE calc_friction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE intnuc(inuc,waters,ys,bottom, vels,b, Q, A)
    !! This is for getting the nuc term-- we get the portion inside the derivative,
    ! (called inuc) and then differentiate it within the mcCormack scheme

    INTEGER, INTENT(IN):: b
    REAL(dp), INTENT(IN):: vels, waters,bottom,Q,A, ys
    REAL(dp), INTENT(IN OUT):: inuc
    DIMENSION vels(b), bottom(b), ys(b)

    INTEGER:: i,j
    REAL(dp):: meanvel,depths(b), easy(b)

    DO i= 1, b
        depths(i)= max(waters-bottom(i),0.)
    END DO

    !First integrate the mean velocity.- trapz rule
    !Discharge - we call it meanvel for economy!
    !meanvel= .5*sum(vels(2:b)*(depths(2:b))*(ys(2:b)-ys(1:b-1)) + &
    !vels(1:b-1)*(depths(1:b-1))*(ys(2:b)-ys(1:b-1)))
    easy= vels*depths
    meanvel= .5_dp*sum((easy(2:b)+easy(1:b-1))*(ys(2:b)-ys(1:b-1)))

    !REAL Mean velocity
    meanvel=meanvel/(.5_dp*sum( (depths(2:b)+depths(1:b-1))*(ys(2:b)-ys(1:b-1)) ) )
    !meanvel=meanvel/ (.5*sum(depths(2:b)*(ys(2:b)-ys(1:b-1)) + &
    !depths(1:b-1)*(ys(2:b)-ys(1:b-1)) ) )


    !Now we calculate the actual integral- trapz rule
    !inuc= .5*sum((vels(2:b)-meanvel)**2*(depths(2:b))*(ys(2:b)-ys(1:b-1)) +&
    !(vels(1:b-1)-meanvel)**2.*(depths(1:b-1))*(ys(2:b)-ys(1:b-1))) 
    easy=(vels-meanvel)**2*depths

    !easy=(vels-Q/A)**2*depths

    inuc= .5_dp*sum((easy(2:b)+easy(1:b-1))*(ys(2:b)-ys(1:b-1))) 

    IF(abs(Q)>0._dp) THEN
            inuc=inuc/(Q**2/A)
    ELSE
            inuc=0._dp
    END IF
    !print*, inuc
    !inuc=0._dp

END SUBROUTINE intnuc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE findn(NN,inuc,A,Q,vels,aa,b, delX,DT,waters, bed,acUdlast,NN_old)
    ! This is to estimate 'N', the neglected terms in the longitudinal momentum
    ! balance
    INTEGER, intent(in):: aa,b
    REAL(dp), intent(in out):: NN
    REAL(dp), intent(in)::inuc, A, Q, vels, delX, DT, waters, bed,acUdlast, NN_old

    DIMENSION NN(aa,b), inuc(b), A(b), Q(b), vels(aa,b), waters(b), bed(aa,b), acUdlast(aa,b), NN_old(aa,b)

    integer:: i,j,jj,jjj
    REAL(dp):: barU(b), acU(aa,b), depth(aa,b), NNlast(aa,b),xx(aa,b)

    NNlast=NN

    !print*, "b=", b
    DO i = 1, b

            barU(i)= Q(i)/A(i) !Mean velocity
            acU(:,i)= vels(:,i)-barU(i) !Deviation from mean velocity
            depth(:,i)= max(waters(i)-bed(:,i), 0._dp)
    END DO

    !Calculate everywhere but at the boundaries
    DO i=2,b-1

            IF(A(i)<=0._dp) THEN
                    NN = 0._dp
                    return
            END IF


            DO j= 1, aa
                    IF (vels(j,i)>=0._dp) THEN 
            !	IF (barU(i)>=0._dp) THEN 
                            jj=1
                            jjj=1
                    ELSE
                            jj=0
                            jjj=0
                    END IF

                    !IF(barU(i)>0._dp) THEN
                    !	jjj=1
                    !ELSE
                    !	jjj=0
                    !END IF

                    NN(j,i)= -(depth(j,i))*(inuc(i+1-jj)*Q(i+1-jj)**2/A(i+1-jj)-inuc(i-jj)*Q(i-jj)**2/A(i-jj))/(delX*A(i)) + & 
                    (acU(j,i+1-jj)**2*depth(j,i+1-jj)-acU(j,i-jj)**2*depth(j,i-jj))/(delX) + &
                    2._dp*depth(j,i)*acU(j,i)*(barU(i+1-jj)-barU(i-jj))/(delX) + & 
                    barU(i)*(depth(j,i+1-jjj)*acU(j,i+1-jjj) - depth(j,i-jjj)*acU(j,i-jjj))/(delX)+&
                    1._dp*(depth(j,i)*acU(j,i) - acUdlast(j,i))/DT
                    
                    if(vels(j,i)==0._dp) THEN
                    NN(j,i)= -(depth(j,i))*(inuc(i+1)*Q(i+1)**2/A(i+1)-inuc(i-1)*Q(i-1)**2/A(i-1))/(2._dp*delX*A(i)) + & 
                    (acU(j,i+1)**2*depth(j,i+1)-acU(j,i-1)**2*depth(j,i-1))/(2._dp*delX) + &
                    2._dp*depth(j,i)*acU(j,i)*(barU(i+1)-barU(i-1))/(2._dp*delX) + & 
                    barU(i)*(depth(j,i+1)*acU(j,i+1) - depth(j,i-1)*acU(j,i-1))/(2._dp*delX)+&
                    1._dp*(depth(j,i)*acU(j,i) - acUdlast(j,i))/DT
                    end if
            END DO
    END DO

    !Calculate at boundaries
    do j=1,aa
            if(vels(j,1)>0._dp) THEN
                    NN(j,1)=0._dp
            !	NN(j,1)=(depth(j,1)*acU(j,1) - acUdlast(j,1))/DT !NN(:,2)!2._dp*NN(:,2)-NN(:,3)
            ELSE
                    NN(j,1)=-(depth(j,1))*(inuc(1+1)*Q(1+1)**2/A(1+1)-inuc(1)*Q(1)**2/A(1))/(delX*A(1)) + & 
                    (acU(j,1+1)**2*depth(j,1+1)-acU(j,1)**2*depth(j,1))/(delX) + &
                    2._dp*depth(j,i)*acU(j,i)*(barU(1+1)-barU(1))/(delX) + & 
                    barU(i)*(depth(j,1+1)*acU(j,1+1) - depth(j,1)*acU(j,1))/(delX)+&
                    1._dp*(depth(j,i)*acU(j,i) - acUdlast(j,i))/DT
            END IF

            if(vels(j,b)>0._dp) THEN
                    NN(j,b)=NN(j,b-1)!NN(:,b-1)!2._dp*NN(:,b-1)-NN(:,b-2)
            ELSE
                    NN(j,b)=0._dp
            END IF
    END DO
    !NN(:,1)=NN(:,2)
    !NN(:,b)= NN(:,b-1)

    !NN=.5_dp*NN+.5_dp*NN_old

    !call random_number(xx)
    !write(12, *) xx
    !do i=1,b
    !	do j=1,aa
    !	if(xx(j,i)>.3_dp) NN(j,i)=.9_dp*NNlast(j,i)+.1_dp*NN(j,i)	
    !NN(j,i) = NN(j,i) +.5_dp*NN_old!NN*pi/2._dp- (pi/2._dp-1._dp)*NNlast!1.1_dp*NN-.1_dp*NNlast
    !	END DO
    !END DO
    !print*, "NNrange", maxval(NN), minval(NN)
    !print*, "acUrange", maxval(acU), minval(acU)
    !print*, "inucrange", maxval(inuc), minval(inuc)
    !NN=0.

end subroutine findn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine findn2(NN,inuc,A,Q,vels,aa,b, delX,DT,waters, bed,acUlast,NN_old)
    ! This is to estimate 'N', the neglected terms in the longitudinal momentum
    ! balance
    INTEGER, intent(in):: aa,b
    REAL(dp), intent(in out):: NN
    REAL(dp), intent(in)::inuc, A, Q, vels, delX, DT, waters, bed,acUlast, NN_old

    DIMENSION NN(aa,b), inuc(b), A(b), Q(b), vels(aa,b), waters(b), bed(aa,b), acUlast(aa,b), NN_old(aa,b)

    integer:: i,j,jj,jjj
    REAL(dp):: barU(b), acU(aa,b), depth(aa,b), NNlast(aa,b),xx(aa,b)

    NNlast=NN

    !print*, "b=", b
    DO i = 1, b

            barU(i)= Q(i)/A(i) !Mean velocity
            acU(:,i)= vels(:,i)-barU(i) !Deviation from mean velocity
            depth(:,i)= max(waters(i)-bed(:,i), 0._dp)
    END DO

    !Calculate everywhere but at the boundaries
    DO i=2,b-1

            IF(A(i)<=0._dp) THEN
                    NN = 0._dp
                    return
            END IF


            DO j= 1, aa
                    IF (vels(j,i)>=0._dp) THEN 
            !	IF (barU(i)>=0._dp) THEN 
                            jj=1
                            jjj=1
                    ELSE
                            jj=0
                            jjj=0
                    END IF

                    !IF(barU(i)>0._dp) THEN
                    !	jjj=1
                    !ELSE
                    !	jjj=0
                    !END IF

                    NN(j,i)= -(depth(j,i))*(inuc(i+1-jj)*Q(i+1-jj)**2/A(i+1-jj)-inuc(i-jj)*Q(i-jj)**2/A(i-jj))/(delX*A(i)) + &  !-dnuc/A
                    depth(j,i)*acU(j,i)*(vels(j,i+1-jj)-vels(j,i-jj))/(delX) + & !d*u'*dUbar/dx
                    depth(j,i)*barU(i)*(acU(j,i+1-jjj) - acU(j,i-jjj))/(delX)+& !d*Ubar*du'/dx
                    depth(j,i)*(acU(j,i) - acUlast(j,i))/DT !d*du'/dt
                    
                    if(vels(j,i)==0._dp) THEN
                    NN(j,i)= -(depth(j,i))*(inuc(i+1)*Q(i+1)**2/A(i+1)-inuc(i-1)*Q(i-1)**2/A(i-1))/(2._dp*delX*A(i)) + & 
                    depth(j,i)*acU(j,i)*(vels(j,i+1)-vels(j,i-1))/(2._dp*delX) + & 
                    depth(j,i)*barU(i)*(acU(j,i+1) - acU(j,i-1))/(2._dp*delX)+&
                    depth(j,i)*(acU(j,i) - acUlast(j,i))/DT
                    end if
            END DO
    END DO

    !Calculate at boundaries
    do j=1,aa
            if(vels(j,1)>0._dp) THEN
                    NN(j,1)=0._dp
            !	NN(j,1)=(depth(j,1)*acU(j,1) - acUlast(j,1))/DT !NN(:,2)!2._dp*NN(:,2)-NN(:,3)
            ELSE
                    NN(j,1)=-(depth(j,1))*(inuc(1+1)*Q(1+1)**2/A(1+1)-inuc(1)*Q(1)**2/A(1))/(delX*A(1)) + & 
                    depth(j,i)*acU(j,i)*(vels(j,1+1)-vels(j,1))/(delX) + & 
                    depth(j,i)*barU(i)*(acU(j,1+1) - acU(j,1))/(delX)+&
                    depth(j,i)*(acU(j,i) - acUlast(j,i))/DT
            END IF

            if(vels(j,b)>0._dp) THEN
                    NN(j,b)=NN(j,b-1)!NN(:,b-1)!2._dp*NN(:,b-1)-NN(:,b-2)
            ELSE
                    NN(j,b)=0._dp
            END IF
    END DO
    !NN(:,1)=NN(:,2)
    !NN(:,b)= NN(:,b-1)

    !NN=.5_dp*NN+.5_dp*NN_old

    !call random_number(xx)
    !write(12, *) xx
    !do i=1,b
    !	do j=1,aa
    !	if(xx(j,i)>.3_dp) NN(j,i)=.9_dp*NNlast(j,i)+.1_dp*NN(j,i)	
    !NN(j,i) = NN(j,i) +.5_dp*NN_old!NN*pi/2._dp- (pi/2._dp-1._dp)*NNlast!1.1_dp*NN-.1_dp*NNlast
    !	END DO
    !END DO
    !print*, "NNrange", maxval(NN), minval(NN)
    !print*, "acUrange", maxval(acU), minval(acU)
    !print*, "inucrange", maxval(inuc), minval(inuc)
    !NN=0.

end subroutine findn2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine findn3(NN,inuc,A,Q,vels,aa,b, delX,DT,waters, bed,acUdlast,NN_old, ys, l, u)
    ! This is to estimate 'N', the neglected terms in the longitudinal momentum
    ! balance
    INTEGER, intent(in):: aa,b, l, u
    REAL(dp), intent(in out):: NN
    REAL(dp), intent(in)::inuc, A, Q, vels, delX, DT, waters, bed,acUdlast, NN_old, ys

    DIMENSION NN(aa,b), inuc(b), A(b), Q(b), vels(aa,b), waters(b), bed(aa,b), acUdlast(aa,b), &
                NN_old(aa,b), ys(aa,b), l(b), u(b)

    integer:: i,j,jj,jjj
    REAL(dp):: barU(b), acU(aa,b), depth(aa,b), NNlast(aa,b),xx(aa,b), VDF(aa,b)

    NNlast=NN

    !print*, "b=", b
    DO i = 1, b

            barU(i)= Q(i)/A(i) !Mean velocity
            acU(:,i)= vels(:,i)-barU(i) !Deviation from mean velocity
            depth(:,i)= max(waters(i)-bed(:,i), 0._dp)
    END DO

    !Calculate V*depth, by integrating - dY/dt - dUd/dx
    VDF=0._dp 
    do i=2, b-1
    !Calculate the integral of dUd/dx. Note that we use u(i)+1 as a trick
    if(vels(aa/2,i)>0._dp) THEN
    jj=1
    else
    jj=0
    end if
            do j=max(l(i),2), min(u(i)+1,aa)
                    VDF(j,i) = VDF(j-1,i) + 0.5_dp*(-vels(j,i+1-jj)*max(waters(i+1-jj)-bed(j,i+1-jj),0._dp) - &
                     vels(j-1,i+1-jj)*max(waters(i+1-jj)-bed(j-1,i+1-jj),0._dp) + & 
                    vels(j,i-jj)*max(waters(i-jj)-bed(j,i-jj), 0._dp) + &
                    vels(j-1,i-jj)*max(waters(i-jj)-bed(j-1,i-jj), 0._dp) )& 
                    /(delX)*(ys(j,i)-ys(j-1,i))
            end do
            !Now calculate dY/dT and correct - note the use of VDF(u(i)+1,i), which is later set to zero
            VDF(l(i):u(i),i) = VDF(l(i):u(i),i) - &
            (ys(l(i):u(i),i)-(ys(max(l(i)-1,1),i)) )*&
            (VDF(min(u(i)+1,aa),i)/((ys(min(u(i)+1,aa),i))-(ys(max(l(i)-1,1),i) ) )) !Note that the last (VDN(u(i).. ) is dY/dt 

            VDF(min(u(i)+1,aa),i)=0._dp
    end do
    !Calculate everywhere but at the boundaries
    DO i=2,b-1

            IF(A(i)<=0._dp) THEN
                    NN = 0._dp
                    return
            END IF


            DO j= 1, aa
                    IF (vels(j,i)>=0._dp) THEN 
            !	IF (barU(i)>=0._dp) THEN 
                            jj=1
                            jjj=1
                    ELSE
                            jj=0
                            jjj=0
                    END IF

                    !IF(barU(i)>0._dp) THEN
                    !	jjj=1
                    !ELSE
                    !	jjj=0
                    !END IF

                    NN(j,i)= -(depth(j,i))*(inuc(i+1-jj)*Q(i+1-jj)**2/A(i+1-jj)-inuc(i-jj)*Q(i-jj)**2/A(i-jj))/(delX*A(i)) + & 
                    (acU(j,i+1-jj)**2*depth(j,i+1-jj)-acU(j,i-jj)**2*depth(j,i-jj))/(delX) + &
                    2._dp*depth(j,i)*acU(j,i)*(barU(i+1-jj)-barU(i-jj))/(delX) + & 
                    barU(i)*(depth(j,i+1-jjj)*acU(j,i+1-jjj) - depth(j,i-jjj)*acU(j,i-jjj))/(delX)+&
                    1._dp*(depth(j,i)*acU(j,i) - acUdlast(j,i))/DT+ &!Unsteady
                    (acU(min(j+1,aa),i)*VDF(min(j+1,aa),i) - acU(max(j-1,1),i)*VDF(max(j-1,1),i))/&
                    (ys(min(j+1,aa),i)-ys(max(j-1,1),i)) !Secondary flow term
                    
                    if(vels(j,i)==0._dp) THEN
                    NN(j,i)= -(depth(j,i))*(inuc(i+1)*Q(i+1)**2/A(i+1)-inuc(i-1)*Q(i-1)**2/A(i-1))/(2._dp*delX*A(i)) + & 
                    (acU(j,i+1)**2*depth(j,i+1)-acU(j,i-1)**2*depth(j,i-1))/(2._dp*delX) + &
                    2._dp*depth(j,i)*acU(j,i)*(barU(i+1)-barU(i-1))/(2._dp*delX) + & 
                    barU(i)*(depth(j,i+1)*acU(j,i+1) - depth(j,i-1)*acU(j,i-1))/(2._dp*delX)+&
                    1._dp*(depth(j,i)*acU(j,i) - acUdlast(j,i))/DT + &
                    (acU(min(j+1,aa),i)*VDF(min(j+1,aa),i) - acU(max(j-1,1),i)*VDF(max(j-1,1),i))/&
                    (ys(min(j+1,aa),i)-ys(max(j-1,1),i)) !Secondary flow term
                    end if
            END DO
    END DO

    !Calculate at boundaries
    do j=1,aa
            if(vels(j,1)>0._dp) THEN
                    NN(j,1)=0._dp
            !	NN(j,1)=(depth(j,1)*acU(j,1) - acUdlast(j,1))/DT !NN(:,2)!2._dp*NN(:,2)-NN(:,3)
            ELSE
                    NN(j,1)=-(depth(j,1))*(inuc(1+1)*Q(1+1)**2/A(1+1)-inuc(1)*Q(1)**2/A(1))/(delX*A(1)) + & 
                    (acU(j,1+1)**2*depth(j,1+1)-acU(j,1)**2*depth(j,1))/(delX) + &
                    2._dp*depth(j,i)*acU(j,i)*(barU(1+1)-barU(1))/(delX) + & 
                    barU(i)*(depth(j,1+1)*acU(j,1+1) - depth(j,1)*acU(j,1))/(delX)+&
                    1._dp*(depth(j,i)*acU(j,i) - acUdlast(j,i))/DT + &
                    (acU(min(j+1,aa),1)*VDF(min(j+1,aa),1) - acU(max(j-1,1),1)*VDF(max(j-1,1),1))/&
                    (ys(min(j+1,aa),1)-ys(max(j-1,1),1)) !Secondary flow term
            END IF

            if(vels(j,b)>0._dp) THEN
                    NN(j,b)=NN(j,b-1)!NN(:,b-1)!2._dp*NN(:,b-1)-NN(:,b-2)
            ELSE
                    NN(j,b)=0._dp
            END IF
    END DO
    !NN(:,1)=NN(:,2)
    !NN(:,b)= NN(:,b-1)

    !NN=.5_dp*NN+.5_dp*NN_old

    !call random_number(xx)
    !write(12, *) xx
    !do i=1,b
    !	do j=1,aa
    !	if(xx(j,i)>.3_dp) NN(j,i)=.9_dp*NNlast(j,i)+.1_dp*NN(j,i)	
    !NN(j,i) = NN(j,i) +.5_dp*NN_old!NN*pi/2._dp- (pi/2._dp-1._dp)*NNlast!1.1_dp*NN-.1_dp*NNlast
    !	END DO
    !END DO
    !print*, "NNrange", maxval(NN), minval(NN)
    !print*, "acUrange", maxval(acU), minval(acU)
    !print*, "inucrange", maxval(inuc), minval(inuc)
    !NN=0.

end subroutine findn3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE hydro_xsect
