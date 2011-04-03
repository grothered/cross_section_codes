!!!!!!!!!!!!!!!!!!!!

Module crosssection

use pizzutotry
!use m_mrgrnk  !A sorting routine
!use matrix_solvers !Tridiagonal matrix solvers from Linpack and Lapack
!use omp_lib

IMPLICIT NONE
!INTEGER, PARAMETER, private  :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER, PARAMETER, PUBLIC  :: dp = SELECTED_REAL_KIND(12, 60)
REAL(dp), PARAMETER, private ::  pi=atan(1._dp)*4._dp !3.14159265259_dp
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
SUBROUTINE shear(nn,ys,bed,water,wslope,taus,ks, f,NNN,slopes, counter, Q, vegdrag, rho,g, lambdacon, tbston, & 
                ysl,ysu, bedl,bedu, high_order_shear)
    ! Solve a model for the distribution of bed shear and velocity over a
    ! cross-section

    ! Solves the equation:
    ! rho*g*Sf*depth = rho*f/8*U^2*sqrt(1+(dbed/dy)^2) + rho*vegdrag*depth*U^2
    !                 + d/dy(rho*0.5*lambdacon*depth^2*sqrt(f/8)*dU^2/dy) 
    ! for U^2 (and thus tau = rho*f/8*U^2)

    !nn is the number of bed points
    !ys is a vector of y coordinates for the bed points
    !bed is a vector of z coordinates for the bed points
    !water is the free surface elevation
    !wslope is the water slope. Note that even if it is incorrect, the shear
    !distribution will be correct to a constant multiple, because the equation
    !is linear in wslope. 
    !taus is the bed shear 
    !ks is a defunct parameter
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
    REAL(dp), INTENT(IN):: ys, bed,water,wslope, ks, f,NNN,slopes, Q, vegdrag, rho, lambdacon, ysl, ysu, bedl, bedu,g 
    REAL(dp), INTENT(IN OUT):: taus !shear, side slope
    LOGICAL, INTENT(IN):: tbston, high_order_shear
    DIMENSION:: ys(nn),bed(nn),taus(nn), ks(nn), f(nn),NNN(nn), slopes(nn), vegdrag(nn)

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
        !print*, 'nn = ', nn
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
        IF((i>1).and.(i<nn-1).and.(const_mesh).and.(high_order_shear)) THEN
            ! Use high order derivative
            tmp = dy_outer*1.0_dp/24.0_dp*Bf(i)*1._dp/dyf(i)
            tmp1 = dy_outer*9.0_dp/8.0_dp*Bf(i)*1._dp/dyf(i)

            alpht2(i) = alpht2(i) + tmp
            alpht(i) = alpht(i)   - tmp1           
            diag(i) = diag(i)     + tmp1
            alphb(i) = alphb(i)   - tmp
        ELSE
            ! Use central derivative
            tmp = dy_outer*Bf(i)*1._dp/dyf(i)
            alpht(i) =alpht(i) -tmp
            diag(i) = diag(i)  +tmp
        END IF

        IF((i>2).and.(i<nn).and.(const_mesh).and.(high_order_shear)) THEN
            ! Use high order derivative
            tmp = dy_outer*1.0_dp/24.0_dp*Bf(i-1)*1._dp/dyf(i-1)
            tmp1 = dy_outer*9.0_dp/8.0_dp*Bf(i-1)*1._dp/dyf(i-1)
            alpht(i) = alpht(i)  -tmp
            diag(i) =  diag(i)   +tmp1
            alphb(i) = alphb(i)  -tmp1
            alphb2(i) = alphb2(i)+tmp
        ELSE
            ! Use central derivative
            tmp = dy_outer*Bf(i-1)*1._dp/dyf(i-1)
            diag(i) = diag(i) +tmp
            alphb(i) = alphb(i) -tmp
        END IF
     
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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE wet(l, u ,nos,water, bed)
    !Determine the wetted boundaries of the cross-section

    ! l is the lower (i.e. leftmost) wetted cross-sectional point
    ! u is the upper (i.e. rightmost) wetted cross-sectional point    
    ! nos is the number of points on the cross-section
    ! water is the water surface elevation
    ! bed is the bed

    INTEGER, INTENT(IN)::nos !number of points representing cross section
    INTEGER, INTENT(OUT):: l, u
    REAL(dp), INTENT(IN):: water, bed
    DIMENSION:: bed(nos)

    INTEGER:: wetdry(nos), lenwet(nos), lastdry, width, i, maxloc


    !Predefine  in case none are wet - the code will recognise this. 
    l=0
    u=0

    !Step in from the left until we find a wet point
    DO i =1, nos

        IF(water-bed(i)>0.0_dp) THEN
            l=i
            exit
        END IF


    END DO

    !Now step in from the right until we find a wet point
    DO i =1, nos

        IF(water-bed(nos-i+1)>0.0_dp) THEN
            u=nos-i+1
            exit
        END IF

    END DO

    RETURN 

END SUBROUTINE wet
!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE reader(a,l)

    !Character, Intent(In):: b
    INTEGER, INTENT (IN) :: l
    REAL(dp), INTENT(IN OUT):: a(l)
    
    open(66,file= 'depths', status="old")
    read(66,*) a
    close(66)

    !open(66,file= 'widths', status="old")
    !read(66,*) b
    !close(66)


END SUBROUTINE reader
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE readcs(bed,ys,a,b)
    INTEGER, INTENT(IN)::a, b
    REAL(dp), INTENT (IN OUT):: bed(a,b), ys(a,b)
    REAL(dp):: holder(a*b), holder2(a*b)
    
    open(77,file="sectionsold2",status="old")
    read(77,*) holder

    open(88,file="lnthsold2",status="old")
    read(88,*) holder2

    bed= reshape(holder,(/a,b/))
    ys= reshape(holder2,(/a,b/))
    
    close(77)
    close(88)

END SUBROUTINE readcs
!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE set_geo(bed, ys, waters,recrds,fs,a,b,hlim, readin, water_m, water_mthick) !initial geometry
    ! A routine to set initial geometry conditions for the quasi-2D model

    INTEGER, INTENT(IN):: a, b 
    REAL(dp), INTENT(IN OUT):: bed, waters,recrds,fs,ys
    REAL(dp), INTENT(IN):: hlim, water_m, water_mthick
    LOGICAL, INTENT(IN):: readin
    DIMENSION bed(a,b),waters(b),recrds(a,b), ys(a,b),fs(a,b)
    
    INTEGER:: i, j,m,n


    IF(readin) call readcs(bed,ys,a,b)  

    DO j= 1, b

        DO i= 1, a

            IF(readin.EQV..false.) bed(i,j) = 5.5_dp- 14.5*(((1._dp*(b-j))/(1._dp*b)))**1.2_dp+ &
            (6._dp)*abs((i*1._dp-a*0.5_dp-0.5_dp)/(a*0.5_dp) )**1.2_dp  !4.5_dp- 8.5*(((1._dp*(b-j))/(0.6_dp*b)))**1.2_dp+ &
            !(6._dp)*abs((i/1._dp-a/2._dp-0.5_dp)/(a/2._dp) )**1.2_dp  !- 4.*(((1.*(b-j))/(1.*b)))**1.2  !+ (4.)*abs((i-a/2-0.5)/(a/2) )**1.2  ! -7.*exp(-j*1.0/(4.*b)) !initial valley waterations  !0.000002*abs((i-a/2))**2 - 5.*(100.-j)/100.
            fs(i,j)=0.032_dp

        END DO !i
    END DO !j

    bed(1,:)= 100. !Side wall!
    bed(a, :)= 100. !Side Wall!

    DO j = 1, b
        DO i= 1, a
            IF(readin.EQV..false.) ys(i,j)= 1._dp*(i-1)*1000._dp/(1._dp*a-1._dp)
        END DO
    END DO

    DO j= 1, b
        waters(j)= max(bed(a/2,j)+hlim/1.001_dp+water_mthick, water_m) !0.d0 !hs(a/2,j)+2.!0d0 !initial water surface wateration
    END DO 

END SUBROUTINE set_geo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE meanvars(bed,ys,waters,fs,a,b,u,l,Width,Area,bottom, dbdh, even,hlim) !use this to get suitably averaged variables for input into the 1D code 
    INTEGER, INTENT(IN)::a,b
    INTEGER, INTENT(IN OUT):: l,u
    LOGICAL, INTENT(IN):: even !!To enforce even cross sections. Sometimes the wetting and drying routine can fall over if we require evenness, in situations with say several pools of water 
    REAL(dp), INTENT(IN):: bed,ys,waters,hlim
    REAL(dp), INTENT(INOUT):: Width, Area, bottom, dbdh !Averaged variables
    REAL(dp), INTENT(IN OUT):: fs
    DIMENSION bed(a,b), ys(a,b), waters(b), fs(a,b), Width(b), Area(b), bottom(b),  l(b), u(b),dbdh(b,2)
    INTEGER:: i,ll,uu,j, wetpts(a)

    REAL(dp):: increm(a) , edgel, edgeu, nom, db, Area_old(b)
    LOGICAL:: alldry, allwet
    !print*, "here I am"

    DO i= 1, b !For every cross-profile
        !ll=l(i)
        !uu=u(i)
        call  wet(l(i), u(i), a, waters(i), bed(:,i))  !Find the rightmost and leftmost wetted points
        wetpts(1:a)=0 !!!Find those points that might be dry internally. This is done because it is
        !possible that between l(i) and u(i) there are "mid channel bars" which are dry. 
        alldry=.true. !Are all the points dry?
        allwet=.true. !Are all the points wet?
        !print*, l(i), u(i)
        IF ((u(i).ge.l(i)).AND.(l(i)>0)) THEN
            DO j= l(i), u(i)
                IF(waters(i)-bed(j,i)<0._dp) THEN !This point is dry
                        !wetpts(j)=0
                        IF(allwet) allwet=.false.

                ELSE
                        wetpts(j)=1
                        IF(alldry) alldry=.false.

                END IF
            END DO 
        END IF

        ll=l(i)
        uu=u(i)

        Area_old(i)=Area(i)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Now we calculate a bunch of averaged variables.
        !Note that since we account for the possibility of interior dry points,
        !everything is a little more complex. So we first treat the case where there are
        !no interior dry points, then the case when there are mid channel bars, and then
        !the case where the cross section is basically dry

        !FIXME: Convert this IF to a SELECT CASE statement
        IF((allwet).and.(.not.alldry)) THEN !!First case - no dry interior points
            IF(ll==0) print*, 'll=0', ll, uu

            !Calculate the 'edge bits' that we add to the width	
            IF((uu<a).and.(ll>1)) THEN
                edgel=(ys(ll,i)-ys(ll-1,i))*(waters(i)-bed(ll,i))/(bed(ll-1,i)- bed(ll,i))  !So this is the wet distance to the dry side of ll
                edgeu= (ys(uu+1,i)-ys(uu,i))*(waters(i)-bed(uu,i))/(bed(uu+1,i)- bed(uu,i))
                
                IF(edgel>1.00001_dp*(ys(ll,i)-ys(ll-1,i))) THEN 
                    print*,"edgel prob", bed(ll-1,i), bed(ll,i), waters(i), ys(ll,i), ys(ll-1, i), ll, i
                ENDIF
                IF(edgeu>1.00001_dp*(ys(uu+1,i)-ys(uu,i))) THEN 
                    print*, "edgeu prob", edgeu, (ys(uu+1,i)-ys(uu,i)),& 
                    edgeu-(ys(uu+1,i)-ys(uu,i)), bed(uu+1,i), bed(uu,i),&
                    waters(i), ys(uu+1,i), ys(uu,i), ll, i
                END IF
            ELSE
                edgel=0._dp
                edgeu=0._dp

            END IF

            increm(ll:uu-1)=ys((ll+1):uu,i)-ys(ll:(uu-1),i) !like a dx term-- it helps us calculate averaged variables like mean depth etc. Not stored	
            Width(i)= ys(uu,i)-ys(ll,i) +edgel+edgeu !The water surface width of the section 
            Area(i)= .5_dp*(sum(increm(ll:uu-1)*(max(waters(i)-bed(ll:uu-1,i),0._dp)+max(waters(i)-bed(ll+1:uu,i),0._dp))) + & 
            edgel*(waters(i)-bed(ll,i)) +edgeu*(waters(i)-bed(uu,i)) ) !Cross sectional area
            bottom(i)= waters(i) - (Area(i)/Width(i)) !Mean bottom elevation
            !bottom(i)=sum(increm(ll:(uu-1))*0.5_dp*( bed(ll:(uu-1), i)*wetpts(ll:(uu-1))+bed((ll+1):uu,i)*wetpts((ll+1):uu))) & 
            !+edgel*(bed(ll,i)+(waters(i)-bed(ll,i))/2._dp) +edgeu*(bed(uu,i)+(waters(i)-bed(uu,i))/2.) !!Note that the latter terms are the width of the edge bit, times the average depth over this region -- doing it this way seemed to be a big improvement. 
            !	bottom(i)= bottom(i)/Width(i)  !Mean bottom elev
                         
            IF((uu<a).AND.(ll>1)) THEN !The if statements ensure that we are not on the edge of the domain, and that there are at least 2 points between l and u. Still I have not accounted for mid channel dry points. 
                dbdh(i,1:2)= -(ys(ll,i) -ys(ll-1,i))/(bed(ll,i)-bed(ll-1,i)) &
                +  (ys(uu+1,i) -ys(uu,i))/(bed(uu+1,i)-bed(uu,i))
            ELSE
                dbdh(i,1:2)=0._dp
            END IF
           
            
        ELSE    !!Now we consider the case with mid channel bars
            IF ( (uu>ll).AND.(.not.alldry)) THEN 
                Width(i)=0._dp 
                Area(i)=0._dp
                dbdh(i,1:2)=0._dp 

                DO j=ll,uu+1
                    IF(j<a+1) THEN !if j=a+1 we don't need to do this
                        IF(wetpts(j)>0) THEN
                            IF(j>1) THEN
                                IF(wetpts(j-1)>0) THEN
                                    db=ys(j,i)-ys(j-1,i)	
                                    Width(i)= Width(i)+ db
                                    Area(i)= Area(i)+ db*(waters(i)-.5_dp*(bed(j,i)+bed(j-1,i)))
                                ELSE
                                    db=abs( (ys(j,i)-ys(j-1,i))*(waters(i)-bed(j,i))/(bed(j-1,i)-bed(j,i)))
                                    Width(i)= Width(i)+ db
                                    Area(i)= Area(i)+ db*(waters(i)-bed(j,i))*.5_dp
                                    dbdh(i,1)= dbdh(i,1)+ (abs((ys(j,i)-ys(j-1,i))*&
                                                (waters(i)+1.0E-06_dp-bed(j,i))/(bed(j-1,i)-bed(j,i)))-db)/1.0E-06_dp
                                        !Note that the above expression is = ( db(h+delh) - db(h))/delh
                                    dbdh(i,2)=dbdh(i,1)
                                END IF
                            END IF	
                        ELSE ! so wetpts(j)=0, but point j-1 might be wet and we need to account for that
                            IF (j>1) THEN
                                IF(wetpts(j-1)>0) THEN
                                    db=abs((ys(j,i)-ys(j-1,i))*(waters(i)-bed(j-1,i))/(bed(j,i)-bed(j-1,i)))
                                    Width(i)= Width(i)+ db
                                    Area(i)=Area(i)+db*(waters(i)-bed(j-1,i))*.5_dp
                                    dbdh(i,1)= dbdh(i,1)+ &
                    ( abs((ys(j,i)-ys(j-1,i))*(waters(i)+1.0E-06_dp-bed(j-1,i))/(bed(j,i)-bed(j-1,i)))-db)/1.0E-06_dp
                                                        dbdh(i,2)=dbdh(i,1)
                                END IF
                            END IF	
                        END IF
                    END IF	
                END DO

                bottom(i)= waters(i)-Area(i)/Width(i)

            ELSE  ! So this is where ll= u(i) -- i.e. fairly dry! So we put the min water over things to calculate their depth. 
                    IF(alldry) THEN
                       nom=  minval(bed(:,i))+hlim/1.001_dp  !Nominal water elevatioon
                       call wet(l(i),u(i),a,nom, bed(:,i)) !This lets us tell what the width should be when we get inundated-- assume the constant is of adequate size. Prevents flooding with a 0 width.   
                        ll=l(i)
                        uu=u(i)
                    ELSE
                       nom=waters(i) !Nominal water elevation
                    END IF
                    !if(.not.alldry) THEN
                    !So here ll=uu, i.e. we have only one wet point on the cross
                    !section
            
                    IF((ll>1).and.(uu<a)) THEN
                        edgel=max((ys(ll,i)-ys(ll-1,i))*( (nom-bed(ll,i))/(bed(ll-1,i)- bed(ll,i))),0._dp)  !So this is the wet distance to the dry side of l(i)
                       edgeu=max((ys(uu+1,i)-ys(uu,i))*((nom-bed(uu,i))/(bed(uu+1,i)- bed(uu,i))), 0._dp)
                    ELSE
                        edgel=0._dp
                        edgeu=0._dp
                    END IF

                    Width(i)=edgel+edgeu  
                    Area(i)= 0.5_dp*(nom-bed(ll,i))*Width(i) 
                    bottom(i)= nom-Area(i)/Width(i) !Mean bottom elevation
                    
                    dbdh(i,1)= -(ys(ll,i) -ys(ll-1,i))/(bed(ll,i)-bed(ll-1,i)) &
                    +  (ys(uu+1,i) -ys(uu,i))/(bed(uu+1,i)-bed(uu,i)) 
                    dbdh(i,2)=0._dp
            END IF !END of the if statement where we deal with cross sections that either have internal dry points or are entirely dry



        END IF !End of the main calculation, i.e. end of IF(minval(wetpts(ll:uu))>0) THEN 

        !dbdh=0._dp
         
        IF(dbdh(i,1)<0._dp) THEN 
          PRINT*, "dbdh<0", i, l(i),u(i), bed(l(i)-1,i),bed(l(i),i), bed(u(i)+1,i), bed(u(i),i) , ys(l(i)-1:l(i),i), &
          ys(u(i):u(i)+1,i), dbdh(i,1) 
          STOP
        END IF

        !if(bottom(i)>waters(i)-0.03) print*, 'somethings less is going on'
        IF(Width(i)>1000._dp) THEN
            PRINT*, 'Width(', i, ') > 1000 =', Width(i), ll, uu, edgel, edgeu, alldry, allwet, l(i), u(i)
            STOP
        END IF

    END DO

    !dbdh=0._dp
    !if(b>80) write(23,*) Area_old(80)-Area(80)

END SUBROUTINE meanvars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_shear(a, dT, water, Q, bed,ys,Area, Width,bottom, ff,recrd, E, D,C,rmu,inuc,tau,NN,counter &
            ,slopes, hlim,mor,taucrit_dep,layers, taucrit_dep_ys, nos, taucrit, vegdrag, susdist, rho, Qe, & 
            Qbed, rhos, voidf, d50, g, kvis, norm, vertical, lambdacon, tbston, ysl,ysu,bedl,bedu, & 
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
    ! Width = the wetted cross-sectional width
    ! bottom = the average bed elevation
    ! ff = the darcy-weisbach friction factor
    ! recrd = a useful parameter to record things (might not be used here)
    ! E = the cross-sectionally integrated rate of erosion
    ! D = the cross-sectionally integrated rate of depositiona
    ! C = the near bed suspended sediment concentration
    ! rmu = the 'roughness multiplier' -- defined so that the friction slope Sf = rmu * (Q/Area)^2
    !                                  -- so in some sense rmu = (f/8)/mean_depth
    ! inuc = non-uniform convective intertia multiplier term for interaction with a 1D unsteady solver
    ! tau = bed shear
    ! NN = defunct term that I had been using to add secondary flows to the shear model
    ! counter = variable which records the loop iterator in the routine that calls this
    ! slopes = dbed/dy
    ! hlim = do not compute the shear if the average depth is less than hlim
    ! mor = morphological factor (for timestep accelleration)
    ! taucrit_dep = z coordinate of critical shear stress bed layers (at certain depths from the bed surface)
    ! layers = number of taucrit layers
    ! taucrit_dep_ys = y coordinates of taucrit_dep
    ! nos = number of points in taucrit_dep_ys
    ! taucrit = the critical shear stress at every point in taucrit_dep
    ! vegdrag = vegetation drag coefficient
    ! susdist = logical variable -- is lateral variation in C allowed (true) or not (false)
    ! rho = density of water
    ! Qe = rate of resuspension
    ! Qbed = bedload flux 
    ! voidf = void fraction of the bed (porosity)
    ! d50 = median sediment size
    ! g = gravity
    ! kvis = kinematic viscosity
    ! norm = logical -- is resuspension/bedload treated as occuring normal to the bed (true) or in the vertical (false)
    ! vertical = logical -- do we use the 'shear' routine, or one based on Pizzuto (1991)
    ! lambdacon = dimensionless eddy viscosity
    ! tbston = logical -- do we retain the sqrt(1+(dbed/dy)^2) term in the shear routine (true) or treat it as 1 (false)
    ! ysl = y-coordinate just left of the wetted part of the cross-section
    ! ysu = y-coordinate just right of the wetted part of the cross-section
    ! bedl = bed-coordinate just left of the wetted part of the cross-section
    ! bedu = bed-coordinate just right of the wetted part of the cross-section
    ! high_order_shear = logical -- do we try to use a high order approximation for
    ! interior derivatives in the 'shear' routine


    INTEGER, INTENT(IN)::a,counter,layers, nos
    REAL(dp), INTENT(IN)::water,Q, Width, Area, bottom, ff, hlim,mor, vegdrag, dt, rho, rhos, voidf,&
         d50, g, kvis, lambdacon, ysl,ysu,bedl,bedu 
    REAL(dp), INTENT(IN OUT):: bed, recrd, E, D,rmu,inuc,tau, NN, ys,C,taucrit_dep, taucrit_dep_ys, slopes, & 
        taucrit, Qe, Qbed
    LOGICAL, INTENT(IN):: susdist, norm, vertical, tbston, high_order_shear
    DIMENSION bed(a),ys(a), ff(a),recrd(a),tau(a), NN(a),slopes(a),taucrit_dep(nos,layers),C(a),taucrit_dep_ys(nos), & 
        taucrit(nos,0:layers), vegdrag(a), Qe(a), Qbed(a) ! 
    
    INTEGER::i, j, bgwet, up,  jj,  info,ii, n(a)
    REAL(dp)::wslope,  Qelocal, tt, corfact
    REAL(dp):: kkkk(a), tbst(a), f(a)  
    REAL(dp)::dst(a,0:(layers+1)), vel(a), Qb(a),& 
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

    IF(abs(wslope)>.1) print*, "|wslope| >.1", wslope, water-bottom, rmu, Q, Area, Q/Area
    IF(isnan(wslope)) print*, "wslope is nan", Q, A, water-bottom, rmu


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
                    call shear(up-bgwet+1,ys(bgwet:up),bed(bgwet:up),water, wslope ,tau(bgwet:up),kkkk(bgwet:up), & 
                            ff(bgwet:up),NN(bgwet:up), slopes(bgwet:up), counter, Q, vegdrag(bgwet:up), rho,g, & 
                            lambdacon, tbston, ysl,ysu, bedl, bedu, high_order_shear)
                END IF
            END IF

            IF ( (.NOT.dry(i)).AND.(dry(i-1))) THEN !So we reached a new wet boundary following an interior dry point 
                bgwet= i
            END IF

            IF((i==a)) THEN 
                IF (maxval(water-bed(bgwet:a))>.01_dp ) THEN 
                    call shear(i-bgwet+1,ys(bgwet:a),bed(bgwet:a),water, wslope ,tau(bgwet:a),kkkk(bgwet:a),ff(bgwet:a) &
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
                (1._dp*((water-bed(a))*(bedu-bed(a))*(ysu-ys(a)) )+ ys(a)-ys(a-1))/))))**2
            !corfact=1._dp
            tau = tau*corfact
            vel=vel*sqrt(corfact)
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
        call roughmult(a,rmu, vel, Q, Area, width,sqrt(1._dp+slopes**2),& 
            max(water-bed,0._dp), ys, ff, vegdrag, counter, ysl, ysu, bedl, bedu, bed, water,g)
    ELSE
        call roughmult(a,rmu, vel, Q, Area, width,1._dp +0._dp*slopes ,& 
            max(water-bed,0._dp), ys, ff, vegdrag, counter,ysl, ysu, bedl, bedu, bed, water,g)
    END IF

    ! Calculate the integral of the nuc term.
    call intnuc(inuc,water,ys,bed,vel, a, width/a, Q, Area)


    IF(isnan(rmu)) THEN
        print*, "rmu is nan"
        print*, water-bed
        print*, "______"
        print*, water
        !stop
    END IF   
           

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE calc_shear
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_resus_bedload(a, dT, water, Q, bed,ys,Area, Width,bottom, ff,recrd, E, D,C, wset, rmu,a2, inuc,tau,vel,NN & 
    ,counter,slopes, hlim,mor,taucrit_dep,layers, taucrit_dep_ys, nos, taucrit, vegdrag, susdist, rho, Qe, Qbed, rhos,& 
    voidf, dsand, d50, g, kvis, norm, vertical,alpha, tbston, Qbedon, ysl,ysu,bedl,bedu, resus_type, bedload_type, a_ref) 
    ! Calculate the rate of resuspension and bedload transport over a
    ! cross-section

    INTEGER, INTENT(IN)::a,a2,counter,layers, nos
    REAL(dp), INTENT(IN)::water,Q, Width, Area, bottom, ff, hlim,mor, vegdrag,dt, rho, rhos, voidf, dsand, d50, g, & 
        kvis, alpha, wset, a_ref
    REAL(dp), INTENT(IN OUT):: bed, recrd, E, D,rmu,inuc,tau,vel, NN, ys,C,taucrit_dep, taucrit_dep_ys, slopes, taucrit,& 
         Qe, Qbed
    REAL(dp), INTENT(IN):: ysl,ysu,bedl,bedu 
    LOGICAL, INTENT(IN):: susdist, norm, vertical, tbston, Qbedon
    CHARACTER(LEN=20), INTENT(IN):: resus_type, bedload_type
    DIMENSION bed(a),ys(a), ff(a),recrd(a),tau(a),vel(a), NN(a),slopes(a),taucrit_dep(nos,layers),C(a),& 
                taucrit_dep_ys(nos), taucrit(nos,0:layers), vegdrag(a), Qe(a), Qbed(a), a_ref(a) ! 

    INTEGER::i, j, bgwet, up,  jj,  info,ii, n(a)
    REAL(dp)::wslope,  Qelocal, tt, corfact, useme, dgravel
    REAL(dp):: kkkk(a), tbst(a), f(a)  
    REAL(dp)::sllength(a),    dst(a,0:(layers+1)), Qb(a), & 
        bedlast(a), sinsl(a), mu_d, Qtemp, useful(a), Ceq(a) 
    REAL(dp)::writout(a2), d_star, c_a, k_scr, f_cs, si, tmp
    !logical::writ_tau=.TRUE.  !This will write out the cross sectional taus-- will lead to massive files if you're not careful. 
    LOGICAL::  dry(a)


    ! To direct erosion normal to the bed, adjust the erosion factor accordingly.
    IF(norm) THEN
        sllength = sqrt(1._dp+ slopes(1:a)**2)
    ELSE
        sllength=1._dp
    END IF

    Qb=0._dp !Predefine the bedload vector
    
    !LOOP OVER ALL BED POINTS AND CALCULATE RATE OF EROSION/BEDLOAD
    DO i=1, a
        jj=0 !Predefine this variable, which tells us which shear layer we are in
        
        !NOTE THAT DST IS NOT AN ARGUMENT TO SECTUPDATE AT PRESENT -- this confused me
        !before. It means the same thing, although here we have an extra column which
        !accounts for the layer past the last layer -- presumed a long way away!
        dst(i, 0)=0._dp !The distance to the zeroth layer -- always 0.
        !The distance to the shear transition layers, in the vertical
        DO j= 1, layers
            dst(i, j)=max((bed(i)-taucrit_dep(i,j) ), 0._dp) !The distance from the bed surface to the next shear layer 
            !Check for consistency
            IF(bed(i)<taucrit_dep(i,j)) THEN 
                PRINT*, 'taucrit_dep prob again', i, j, bed(i), taucrit_dep(i,j)        
            END IF
        END DO

        dst(i, layers+1)=9.9E+10_dp !The distance to the layer past the last layer -- a convenience 


        !!!!Now figure out into which layer we will initially erode
        DO ii=1, layers
            IF (dst(i,ii)>0._dp) THEN  !If this is true then we are in the layer ii-1 (where layer zero is the bedload type layer)
                jj=ii-1
                goto 333
            ELSE
                   
            END IF
            IF(ii==layers) THEN 
                jj= ii
            END IF
        END DO      
        333 CONTINUE
        ! So now we should know which layer we are cutting into (it is called jj, and ranges over 0,1,.., layers) , and the critical shear in this layer (taucrit(i,jj)).
         IF((minval(dst(i,:))<0._dp).AND.(jj<layers)) THEN 
             PRINT*, "dst<0", dst, bed(i), taucrit_dep(i,:), abs(slopes(i)), i, a 
             STOP
         END IF
        
        ! Calculate the amount of erosion
        Qe(i)=0._dp !This will be redefined below     
        !Predefine for loop  
        tt=0._dp
        Qelocal=0._dp
        !Check
        IF(taucrit(i,jj)<0._dp) THEN
            PRINT*, 'taucrit <0', counter, i,j, slopes(i) 
            STOP
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !EROSION -- There are a number of cases to consider
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(mor>0._dp) THEN

          1456  IF ((abs(tau(i))>taucrit(i, jj)).and.(.true.)) THEN !At least some erosion through this layer occurs
                    
                    SELECT CASE (resus_type)
                        CASE('cohesive') 
                            !! This is the rate of erosion in this layer, in
                            !meters per second of SOLID material, i.e. not
                            !accounting for porosity, which is accounted for
                            !later. 
                            Qelocal=  (alpha/rhos)*( abs(tau(i))-taucrit(i,jj) )**1._dp*& 
                                      (taucrit(i,jj)**(-.5_dp))*sllength(i) 
                                !min(taucrit(i,jj)**(-.5_dp), 50._dp)*sllength(i) 
                        CASE('vanrijn')
                            ! vanrijn 2007 reference concentration method
                            d_star = (d50*((rhos/rho-1._dp)*g/kvis**2)**(1._dp/3._dp))  !Van Rijn d_star parameter
                            c_a = 0.015_dp*max(dsand/d50,1.0_dp)*d50/(a_ref(i)*d_star**0.3_dp)* & 
                                    (max(0._dp,abs(tau(i))-taucrit(i,jj))/taucrit(i,jj))**1.5_dp ! Van Rijn reference concentration, in m^3/m^3     
                            Qelocal = wset*c_a*sllength(i) !/rhos !Rate of erosion in m/s of SOLID material
                            
                        CASE('smithmac')
                            !Based on description in Garcia and Parker
                            tmp = 2.4e-03_dp*(max(0._dp,abs(tau(i))-taucrit(i,jj))/taucrit(i,jj))
                            c_a = 0.65_dp*tmp/(1.0_dp+tmp)
                            Qelocal = wset*c_a*sllength(i) !/rhos !Rate of erosion in m/s of SOLID material

                        CASE DEFAULT
                            print*, 'ERROR: resus_type does not have the correct value in calc_resus_bedload'

                    END SELECT
                    !Now lets make sure that we don't cut through new layers. The next layer
                    !is jj+1 (unless jj == layers in which case there is no new layer). 
                    IF((jj<layers).AND.(Qelocal*mor*(dt-tt)/(1._dp-voidf)> dst(i,jj+1)-dst(i,jj) )) THEN !We are moving into a new shear resistance level, and need to account for this. 
                        !  print*, "cuttin", jj, Qelocal*mor*dt/.6_dp, dst(i,jj+1), dst(i,jj)
                        jj=jj+1 !So this is the next layer     
                        tt=  tt+ (dst(i,jj)-dst(i,jj-1))/(Qelocal*mor)*(1._dp-voidf) !This time taken to cut to this layer
                        GOTO 1456
                    ELSE !So we have eroded, but not entirely to the next layer
                        Qe(i)= dst(i,jj)/(mor*dT)*(1._dp-voidf) + Qelocal*(dt-tt)/dt  
                    END IF !jj<layers and qelocal        
                ELSE 
                    !So here, we (may) have eroded, and now we hit a layer that we
                    !couldn't cut through
                    Qe(i)= (dst(i,jj))/(mor*dT)*(1._dp-voidf)  
                    !But, if we were already in a higher shear layer, then we
                    !better move the points down to reflect that fact.

                END IF !tau>taucrit 
                    

                IF((Qe(i)<0._dp).OR.(isnan(Qe(i)))) THEN 
                    print*, "erosion < 0 or nan", taucrit(i, jj), tau(i), Qe(i)
                    print*, ".........."
                    !print*, i, jj, ys(i), hs(i), taucrit_dep_ys(indd(i,:)), taucrit_dep(indd(i,:), jj), & 
                    !    taucrit_dep(indd(i,:), jj+1)  
                    STOP
                END IF

                DO j=0, jj
                    IF(isnan(taucrit(i,j))) THEN
                        PRINT*, 'isnan taucrit',i,j,jj, slopes(i)
                        STOP
                    END IF
                END DO

        END IF !MOR>0


       !!Bedload - note it is predefined to be zero, so there is no need to
       !change it if tau<taucrit.
        useme=(rhos-rho)*g*d50
        IF( (abs(tau(i))>taucrit(i,0)).AND.(Qbedon)) THEN
            !Qb(i)= C(i)/rhos*(water-bed(i))*sqrt(tau(i)/(rho*ff(i)/8._dp))*mor
            !Qb(i)=0._dp!5._dp*Qe(i)
            SELECT CASE (bedload_type)
                CASE('vanrijn')
                    !See van Rijn (2007 paper 1 eqn 10). This is the flux in
                    !m^3/m/s of SOLID MATERIAL, i.e. not accounting for
                    !porosity. Here we are assuming that the bedload layer has
                    !critical shear equal to the critical shear for suspension
                    !of the upper layer. Now, this might not be the best way to
                    !go.
                    Qbed(i)= 0.5_dp*max(dsand/d50,1._dp)*d50*&
                            (d50*((rhos/rho-1._dp)*g/kvis**2)**(1._dp/3._dp))**(-0.3_dp)*& 
                            rho**(-.5_dp)*(abs(tau(i)))**(.5_dp)*&
                            ( abs(tau(i))-taucrit(i,0))/taucrit(i,0)*sign(1._dp, tau(i)) &
                            *sllength(i)                 
                CASE('mpm')
                ! Meyer-Peter and Muller
                    Qbed(i) = 8.0_dp*((abs(tau(i))-taucrit(i,0))/useme)**(1.5_dp)*sign(1._dp,tau(i))&
                                *sqrt(useme/rho*d50**2._dp)*sllength(i)
                !Qbed(i) = 8.0_dp*(abs(tau(i))/useme-0.047_dp)**(1.5_dp)*sign(1._dp,tau(i))&
                !            *sqrt(useme/rho*d50**2._dp)*sllength(i)
                CASE DEFAULT
                    PRINT*, 'ERROR: bedload_type was not correctly specified in calc_resus_bedload'
            END SELECT
        END IF

    END DO
    
    !!Total erosion - useful for the 1d sus sed 
    E=.5_dp*(sum( (Qe(2:a-1)*(ys(3:a)-ys(1:a-2)))) & 
        + Qe(1)*(ys(2)-ys(1) +1._dp*((water-bed(1))/(bedl-bed(1))*(ys(1)-ysl) )) & 
        + Qe(a)*(1._dp*((water-bed(a))/(bedu-bed(a))*(ysu-ys(a)) )+ys(a)-ys(a-1)) )


END SUBROUTINE calc_resus_bedload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_bed(a, dT, water, Q, bed,ys,Area, Width,bottom, ff,recrd, E, D,C,rmu,a2, inuc,tau,NN,counter&
    ,slopes, hlim,mor,taucrit_dep,layers, taucrit_dep_ys, nos, taucrit, vegdrag, susdist,rho, Qe, Qbed, & 
    wset,dqbeddx, rhos, voidf, d50, g, kvis, norm, vertical, lambdacon, tbston, Qbedon, normmov,sus2d,ysl, & 
    ysu,bedl,bedu,iii, bedlast, susQbal, talmon, high_order_bedload)

    INTEGER, INTENT(IN)::a,a2,counter,layers, nos, iii
    REAL(dp), INTENT(IN)::water,Q, Width, Area, bottom, ff, hlim,mor, vegdrag, dt, rho, Qbed, Qe, dqbeddx, &
        rhos, voidf, d50, g, kvis,wset, lambdacon, ysl,ysu, bedlast!QbedI, dQbedI
    REAL(dp), INTENT(IN OUT):: bed, recrd, E, D,rmu,inuc,tau, NN, ys,C,taucrit_dep, taucrit_dep_ys, slopes, & 
        taucrit, bedu, bedl
    LOGICAL, INTENT(IN):: susdist, norm, vertical, tbston, Qbedon, normmov, sus2d, susQbal, talmon, high_order_bedload
    DIMENSION bed(a),ys(a), ff(a),recrd(0:a),tau(a), NN(a),slopes(a),taucrit_dep(nos,layers),C(a),taucrit_dep_ys(nos), & 
        taucrit(nos,0:layers), vegdrag(a), Qbed(a), Qe(a), dqbeddx(a), bedlast(a) ! 

    INTEGER::i, j, bgwet, up, bfall, jj,dstspl, jjj, minX,maxX, storindex(a), info,ii, indd(a,layers), n(a), b(1)
    REAL(dp):: val, tmp1,tmp2 
    REAL(dp)::wslope, Qeold, tt, Qelocal, bed_tmp(0:a+1), ys_tmp(0:a+1)!, dqbeddx(a)
    REAL(dp):: kkkk(a), tbst(a), f(a), Qd(a), h_rhs(0:a+1), newxs(a), newys(a), lcslp(a)
    REAL(dp)::vinext,vi, vilast, epsl, epsu , mxslopes(a), erode, newslope(a), slim,p1, w1,w2,eqspc,dsts(a), dsts2(a), taud(a) &
        , sllength(a), hlast1,hlast2, storeys(a), C_diag(a),C_upper(a),C_lower(a), C_out(a),dDdy(a), dst(a,0:(layers+1)), taucinc,& 
        Bdist,p11,p12,taucrit_depnew(nos,layers),edvis(a),dedvisdy(a),edvisd(a),dedvisddy(a),d2hdy2(a),rnd(a),hss2_deriv(a), & 
        h_diag(0:a+1), h_upper(0:a+1), h_lower(0:a+1), h_rhs2(0:a+1), vel(a), qb_G(0:a+1), sinf(a), bednew(a),&
        hss22(a),bedlast_tmp(a), sinsl(a), mu_d, Qtemp, useful(a), Ceq(a), h_lower2(0:a+1),h_upper2(0:a+1), &
        ecu, zetamult(a), upper(a), lower(a), diag(a), rhs(a), dyf(a), upper2(a), lower2(a), diag2(a),rhs2(a), &
        edvisy(a), dum1(a), dum2(a), dum3(a), dum4(a),spec(a),homo(a), sllength2(a)  

    REAL(dp)::writout(a2), impcon, useme
    REAL(dp)::bandmat(7,0:a+1), dbeddyH(4,0:a), AFB(7,0:a+1), RRR(0:a+1), CCC(0:a+1), XXX(0:a+1,1), rcond, ferr, berr,work(3*(a+2))
    CHARACTER(1)::EQUED
    INTEGER:: IPV(a+2), iwork(a+2)
    LOGICAL:: crazy_check=.false., dry(a), cnan

    !! Adjust the erosion factor if erosion is normal to the bed
    IF(norm) THEN
        sllength = sqrt(1._dp+ slopes(1:a)**2)
    ELSE
        sllength=1._dp
    END IF

    !!!!!!!!!!!!Check is we should apply the LDM
    !Don't update if the water is very shallow. This is pragmatic as the hydro
    !solver is not very good in this region 
    IF((water-bottom<0.0_dp).OR.(a<2)) THEN 
    !        tau=0._dp+0._dp*bed
    !        E=0._dp
    !        D=0._dp
    !        rmu=maxval(ff)
            RETURN
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate deposition rate
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SELECT CASE(sus2d)
        CASE(.FALSE.)
            ! When we do not use the fully 2D suspended load routine
            WHERE (bed<water)
                Qd=(wset/rhos)*C
            ELSEWHERE
                Qd = 0.0_dp
            END WHERE

        CASE(.TRUE.)
            ! Where we do use the 2D suspended load routine
            WHERE((abs(tau)>0._dp).and.(wset>0._dp))
                    Qd=(wset/rhos)*C*& 
                        min( &
                        (water-bed)/(.1_dp*sqrt(abs(tau)/rho)*(water-bed)/wset*& 
                        (1._dp-exp(-wset/(.1_dp*sqrt(abs(tau)/rho)*(water-bed))*(water-bed)))) &
                        , 20._dp) 
            ELSEWHERE
                    Qd= (wset/rhos)*C*20._dp 
                    !So in this case the shear is zero, and whatever deposits
                    !should deposit fast 
            END WHERE 
    END SELECT

    ! Reality check
    DO i=1,a
        IF(isnan(Qd(i)).or.(Qd(i)<0._dp)) THEN
            PRINT*, "Qd(",i,") is nan", Qd(i)
            STOP
        END IF 
    END DO                

    !FIXME: Force intertidal flat elevation to be a set value
    IF(.FALSE.) THEN
        IF(counter < 10) print*, ' WARNING: Intertidal flat elevation is artificially limited'
        DO i=1,a
            IF(bed(i)>0.7_dp) THEN
            Qd(i) = 0._dp
            END IF 
        END DO
    END IF


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SOLVE FOR THE UPDATED GEOMETRY
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    bedlast_tmp=bed
    qb_G=0._dp
    
    IF(Qbedon) THEN
        ! Calculate downslope bedload transport coefficients
        ! Qby = -qb_G*d(bed)/dy
        tmp1 = (rho*g*(rhos/rho-1._dp)*d50) ! A constant in the lateral bedload formula
        DO i=1, a
            IF(tau(i).ne.0._dp) THEN
                !Choose downslope bedload relation
                IF(talmon.eqv..FALSE.) THEN
                    ! Simple downslope bedslope relation 
                    qb_G(i)= -abs(Qbed(i))*sqrt(abs(taucrit(i,0))/tau(i))
                ELSE 
                ! Talmon (1995) relation
                    qb_G(i)=-abs(Qbed(i))*sqrt(tmp1/tau(i)) &
                                    *1._dp/9._dp*(max(water-bed(i),0._dp)/d50)**.3_dp 
                END IF
            ELSE !tmp1>0
                qb_G(i)=0._dp
            END IF !tmp1>0
        END DO
        
        qb_G(0)=0._dp
        qb_G(a+1)=0._dp
        ! qbh_approx takes the values of qb_G at y(i), and replaces them with the
        ! values at y(i+1/2) calculated using parabolic interpolation with slope
        ! based upwinding, or standard linear interpolation. 
        ! Note that there is 1 less output point than input
        ! point. Hence, there is a funny organisation of arguments.
        call qbh_approx(a,ys,qb_G(0:a),qb_G(a+1),bed, bedl, bedu, ysl, ysu, 2)                
    END IF !Qbedon

    IF(.TRUE.) THEN
        IF(counter.eq.1) print*, 'WARNING: EDGE BEDLOAD VALUES DROPPED TO ZERO'
        qb_G(0) = 0.0_dp
        qb_G(a) = 0.0_dp
    END IF
    !Next we check if qb_G =0 everywhere. If so, then there is no point solving the
    !matrix equation.
    IF(maxval(abs(qb_G)).EQ.0._dp) THEN 
        IF(normmov.eqv..false.) THEN

            IF(iii==1) bed = bed + (-dqbeddx + Qd - Qe)*mor*dT/(1._dp-voidf)
            IF(iii==2) bed = bedlast + (-dqbeddx + Qd - Qe)*mor*dT/(1._dp-voidf)

        ELSE !!USE POINTS WHICH MOVE IN THE DIRECTION OF THE VECTOR D-E
        !    print*, 'FIXME: This is not consistent with the new definition of Qe (30/12/2010), which already includes sllength'
        !    stop
        !    bed=bed+( Qd - Qe*cos(atan(slopes)))*mor*dT/(1._dp-voidf)
        !    ys=ys+Qe*sin(atan(slopes))*mor*dT/(1._dp-voidf)

        !    IF(.true.) THEN !Turn on/off undercutting
        !        !Now Undercut as needed
        !        !!March back from the centre toward the left bank (i=1)
        !        DO i= floor(a/2._dp), 2, -1 !Found that this could produce spikes. 
        !            !March back from the left edge to the centre
        !            !do i= 2, floor(a/2._dp) !This can also produce spikes
        !            DO j= 1, i-1
        !                !Look to see if i has undercut i-j
        !                IF(ys(i)<ys(i-j)) THEN
        !                    CONTINUE !Do nothing until we find a point that i has not undercut
        !                ELSE 
        !                    IF(j==1) GOTO 1010  !So in this case there was no undercutting
        !                    
        !                    DO jj=1, j-1 !So in this case we have undercut j-1 points
        !                        ys(i-j+jj) = ys(i-j)+ (ys(i)-ys(i-j))/(1._dp*(j))*jj
        !                        !bed(i-j+jj)= minval(bed((i-j+1):(i)))

        !                        !Should we sort the bed here also?
        !                        b=maxloc( bed( (i-j+jj):(i-1) ) ) !The location of the max height between i-(j-jj) and i-1 
        !                        jjj=b(1)
        !                        bed( (/i-(j-jj), i-j+jj+jjj-1/) ) = bed( (/i-j+jj+jjj-1, i-(j-jj)/)) !So we swap the max height and the i-(j-jj) height 
        !                    END DO
        !                    GOTO 1010
        !                END IF

        !            END DO
        !           
        !           1010 CONTINUE
        !        END DO 
        !        
        !        !!March back from the centre toward the right bank (i=a)
        !        DO i= floor(a/2._dp)+1, a !Found that this can produce spikes
        !            !March back from the right edge to the centre
        !            !do i= a, floor(a/2._dp)+1, -1 !This can also produce spikes
        !            !Look to see if we have undercut point i+j
        !            DO j= 1, a-i
        !                IF(ys(i)>ys(i+j)) THEN
        !                    CONTINUE !No need to do anything until we find a point that we have not undercut
        !                ELSE 
        !                    IF(j==1) GOTO 1011  !So in this case there was no undercutting
        !                    !If we got here with j>1, then there is undercutting to do 
        !                    DO jj=1, j-1
        !                        ys(i+jj) = ys(i)+ (ys(i+j)-ys(i))/(1._dp*(j))*jj
        !                        !bed((i+1):(i+jj-1)) = minval(bed(i:i+jj-1))
        !                        !Should we sort the bed here also?
        !                        b=minloc(bed( (i+jj): (i+(j-1)) ) ) !The location of the min height between i+jj and i+(j-1) 
        !                        jjj=b(1)
        !                        bed( (/i+jj, i+jj+jjj-1/)) = bed((/i+jj+jjj-1, i+jj/)) !So we swap the min height and the i+jj height 
        !                    END DO

        !                    GOTO 1011
        !                END IF
        !            END DO
        !            1011 CONTINUE
        !        END DO       
        !    END IF  !Turn off the undercutting

        END IF !If normmov.eqv..false.

    ELSE !IF(maxval(abs(qb_G)).EQ.0._dp)

        
        !
        h_lower2=0._dp
        h_lower=0._dp
        h_diag=0._dp
        h_upper=0._dp
        h_upper2=0._dp

        !Define convenient new ys and bed variables which include dry
        !boundary points - these will let us include those boundaries in the
        !implicit solution
        ys_tmp(1:a)=ys
        bed_tmp(1:a)=bed
        ys_tmp(0)=ysl
        ys_tmp(a+1)=ysu
        bed_tmp(0)=bedl
        bed_tmp(a+1)=bedu
        !print*, maxval(qb_G(1:a-1))
        impcon = 0.5_dp
        tmp2 = mor*dT*1._dp/(1._dp-voidf) ! Time evolution constants
        
        !Calculate coefficients for bed slope terms
        IF(high_order_bedload) THEN
            call dbeddyH_approx(a,ys,bed,dbeddyH, ysl, ysu, bedl, bedu, 3)
        ELSE
            call dbeddyH_approx(a,ys,bed,dbeddyH, ysl, ysu, bedl, bedu, 2)
        END IF
        !!Set diagonals for Matrix -- termed 'h'
        !! h stores the implicit terms for [1 / (1-voidf)] dbed/dt = -d/dy (qb_G*dbed/dy) + Ds - Es
    
        DO i=1,a
            tmp1 = tmp2*impcon*(2._dp/(ys_tmp(i+1)-ys_tmp(i-1)))
            IF(i>1) h_lower2(i)= tmp1*(              - qb_G(i-1)*dbeddyH(1,i-1))
            h_lower(i)=  tmp1*(qb_G(i)*dbeddyH(1,i) - qb_G(i-1)*dbeddyH(1+1,i-1))
            h_diag(i)=   tmp1*(qb_G(i)*dbeddyH(2,i) - qb_G(i-1)*dbeddyH(2+1,i-1))
            h_upper(i)=  tmp1*(qb_G(i)*dbeddyH(3,i) - qb_G(i-1)*dbeddyH(3+1,i-1))
            IF(i<a) h_upper2(i)= tmp1*(qb_G(i)*dbeddyH(4,i))
      
            ! Store qb_G*dbed/dy 
            IF((i.ne.a)) recrd(i) = qb_G(i)*(dbeddyH(4,i)*bed_tmp(i+2) + &
                                                          dbeddyH(3,i)*bed_tmp(i+1) + &
                                                          dbeddyH(2,i)*bed_tmp(i)   + &
                                                          dbeddyH(1,i)*bed_tmp(i-1))

            IF(i==1)  recrd(i-1) = qb_G(i-1)*(dbeddyH(4,i-1)*bed_tmp(i+1) + &
                                              dbeddyH(3,i-1)*bed_tmp(i)   + &
                                              dbeddyH(2,i-1)*bed_tmp(i-1)   )


            IF(i==a)                  recrd(i) = qb_G(i)*(0._dp                     + &
                                                          dbeddyH(3,i)*bed_tmp(i+1) + &
                                                          dbeddyH(2,i)*bed_tmp(i)   + &
                                                          dbeddyH(1,i)*bed_tmp(i-1) )

        END DO
        ! Make right hand side -- note that the extra division by 'impcon'
        ! cancels the multiplication by impcon above, which is what we need
        h_rhs2(1:a)= -(1._dp-impcon)/impcon*(+h_upper2(1:a)*(/bed_tmp(3:a+1),0._dp /) & 
                        + h_upper(1:a)*bed_tmp(2:a+1) +h_diag(1:a)*bed_tmp(1:a) +h_lower(1:a)*bed_tmp(0:a-1) & 
                        +h_lower2(1:a)*(/0._dp,bed_tmp(0:a-2)/) ) 
        h_diag(1:a)=h_diag(1:a)+1._dp

        !!Boundary conditions. 
        h_lower2(0)=0._dp
        h_lower(0)=0._dp !This will be zero - though it is included in the in the above implicit terms,it will be zero there (dbeddyH(1,0)=0
        !h_lower(0)=tmp2*impcon*(2._dp/(ys_tmp(2)-ys_tmp(0)))*(qb_G(0)*dbeddyH(1,0))
        h_diag(0)=tmp2*impcon*(2._dp/(ys_tmp(2)-ys_tmp(0)))*(qb_G(0)*dbeddyH(2,0))
        h_upper(0)=tmp2*impcon*(2._dp/(ys_tmp(2)-ys_tmp(0)))*(qb_G(0)*dbeddyH(3,0))
        h_upper2(0)=tmp2*impcon*(2._dp/(ys_tmp(2)-ys_tmp(0)))*(qb_G(0)*dbeddyH(4,0))
        h_rhs2(0)= - (1._dp-impcon)/impcon*( h_diag(0)*bed_tmp(0)+h_upper(0)*bed_tmp(1) +h_upper2(0)*bed_tmp(2))
        !Correct h_diag
        h_diag(0)=1._dp+h_diag(0)

        h_lower2(a+1)=tmp2*impcon*(2._dp/(ys_tmp(a+1)-ys_tmp(a-1)))*( - qb_G(a)*dbeddyH(1,a))
        h_lower(a+1)=tmp2*impcon*(2._dp/(ys_tmp(a+1)-ys_tmp(a-1)))*( - qb_G(a)*dbeddyH(2,a))
        h_diag(a+1)=tmp2*impcon*(2._dp/(ys_tmp(a+1)-ys_tmp(a-1)))*( - qb_G(a)*dbeddyH(3,a))
        h_upper(a+1)=0._dp !This will be zero, because dhdy(4,a) will be zero.
        h_upper2(a+1)=0._dp
        h_rhs2(a+1)=- (1._dp-impcon)/impcon*( h_diag(a+1)*bed_tmp(a+1)+h_lower(a+1)*bed_tmp(a) +h_lower2(a+1)*bed_tmp(a-1))
        !Correct h_diag
        h_diag(a+1)=1._dp+h_diag(a+1)


        IF(iii==1) THEN
            !Now define the right_hand side
            h_rhs(1:a)= (bed_tmp(1:a) + (-dqbeddx + Qd - Qe)*tmp2) + h_rhs2(1:a)
            h_rhs(0)=bed_tmp(0) + h_rhs2(0)
            h_rhs(a+1)=bed_tmp(a+1)+h_rhs2(a+1)

            !Call the matrix solver
            !call DGTSV(a,1,h_lower(2:a), h_diag,h_upper(1:a-1),h_rhs,a ,info)
            !if(mod(counter,2).eq.0) THEN
            !call DGTSV(a,1,h_lower(2:a),h_diag,h_upper(1:a-1), h_rhs,a, info)
            !else
            !call DGTSV(a,1,h_upper((a-1):1:-1),h_diag(a:1:-1),h_lower(a:2:-1), h_rhs(a:1:-1),a, info)
            !end if
            
            ! Fill out matrix for LAPACK solution
            DO i=0,a+1
                IF(i<a) bandmat(2+2+1+i-(i+2),i+2)=h_upper2(i)
                IF(i<a+1)  bandmat(2+2+1+i-(i+1),i+1)=h_upper(i)
                bandmat(2+2+1+i-(i),i)=h_diag(i)
                IF(i>0) bandmat(2+2+1+i-(i-1),i-1)=h_lower(i)
                IF(i>1) bandmat(2+2+1+i-(i-2),i-2)=h_lower2(i)
            END DO
            !call DGBSV(a+2,2,2,1, bandmat, 4+2+1, IPV,h_rhs,a+2,info)      
            !hs=h_rhs(1:a)
            !hsu=h_rhs(a+1)
            !hsl=h_rhs(0)
            !The following version of dgbsv supposedly improves the solution and checks for
            !badness
            IF(high_order_bedload) THEN
                !Banded matrix solver from LAPACK for the high order case
                call DGBSVX('E','N', a+2,2,2,1, bandmat(3:7,0:a+1),  2+2+1, AFB(1:7,0:a+1), 4+2+1, IPV(1:a+2),EQUED, RRR(0:a+1), & 
                         CCC(0:a+1), h_rhs(0:a+1), a+2, XXX(0:a+1,1),a+2, rcond, ferr,berr, work(1:(3*(a+2))),iwork(1:(a+2)), info)
            ELSE
                !Tridiagonal matrix solver from LAPACK for the lower order
                !(standard) case
                XXX(0:a+1,1) = h_rhs(0:a+1)
                call DGTSV(a+2, 1, bandmat(6,0:a), bandmat(5,0:a+1), bandmat(4,1:a+1), XXX(0:a+1,1), a+2, info)  
            END IF 
            bed(1:a)=XXX(1:a,1)
            bedl=XXX(0,1)
            bedu=XXX(a+1,1)
            !Check
            IF(info.ne.0) THEN
                IF(high_order_bedload) THEN
                     print*, 'Bed solver DGBSVX ERROR, info=', info, 'rcond = ', rcond
                ELSE
                    print*, 'Bed solver DGTSV ERROR, info=', info, 'rcond = ', rcond 
                END IF
                stop
            END IF

 
            !FIXME: THIS WAS NEEDED TO PREVENT CUMULATIVE GROWTH OF ROUND-OFF
            IF(.FALSE.) THEN
                IF(counter.eq.1) print*, 'Symmetry correction in bed routine'
                DO i=1, floor((a*1._dp)/2._dp)
                    bed(i)=0.5_dp*(bed(i)+bed(a-i+1))
                    bed(a-i+1)=bed(i)
                END DO
                bedl=0.5_dp*(bedl+bedu)
                bedu=bedl
            END IF

        ELSE !iii=1
            !FIXME
            PRINT*, 'TWO STEP VERSION NEEDS EDITING BEFORE IT CAN BE SUPPORTED'
            STOP
            h_rhs= bedlast + (-dqbeddx + Qd - Qe)*mor*dT/(1._dp-voidf) !+Qd*mor*dT/(1._dp-voidf) - dT*mor*Qe/(1._dp-voidf) 
        !!!ENFORCE BOUNDARY CONDITIONS
            h_rhs(1)= bedlast(1) + 2._dp/(ys(2)-ysl)*qb_G(0)*(-bedl/(ys(1)-ysl))*tmp2
            h_rhs(a) = bedlast(a) +2._dp/(ysu-ys(a-1))*qb_G(a)*(-bedu/(ysu-ys(a)))*tmp2
            !Call the matrix solver
            call DGTSV(a,1,h_lower(2:a), h_diag,h_upper(1:a-1),h_rhs,a, info)
            bed=h_rhs
        END IF !iii=1
    END IF !IF(maxval(abs(qb_G)).EQ.0._dp)


    DO i=1, a
        IF(isnan(bed(i))) THEN
            PRINT*, "bed nan", i, Qe(i), Qd(i), C(i)!, QbedI, dQbedI, bedlast_tmp(i-1:i+1), Ceq(i)
            PRINT*, info
            STOP
        END IF

        IF(abs(bedlast_tmp(i)-bed(i))>0.1_dp) THEN
            print*, "bedjump", i, bedlast_tmp(i)-bed(i), bedlast_tmp(i), bed(i), Qe(i), C(i),&
            Qd(i), dqbeddx(i)
        END IF
    END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!Erosion due to undercutting. This does not seem to be the main thing
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !IF(.FALSE.) THEN
    !    !if(normmov.eqv..false.) THEN !Switch on/off undercutting - if normmov is true, then undercutting has already happened
    !    h_rhs=bed
    !    DO i=2,a-1
    !        ecu= abs(Qe(i)*sin(atan(slopes(i))))*mor*dT/(1._dp-voidf)
    !        !So this is the lateral projection of the 'normal to the bed' amount of
    !        !erosion. To account for undercutting, if ecu cuts past neighbouring points,
    !        !then such points are replaced with a weighted average of the new bed(i) and the
    !        !'normally eroded' position of bed(i).
    !        IF(slopes(i)<0._dp) THEN
    !            DO j=1, max(i+1,1)
    !                IF( ecu > ys(i)-ys(i-j)) THEN
    !                    !print*, 'undercut left ', ecu
    !                    h_rhs(i-j)=min( ((ys(i)-ys(i-j))*(bedlast_tmp(i)+ ecu/slopes(i)+Qd(i)*mor*dT/(1._dp-voidf)) +& 
    !                    (ys(i-j)-( ys(i)-ecu) )*bed(i))/ecu, h_rhs(i-j)) !Note that the slopes(i) is tan(atan(slopes(i))). This is an appropriately weighted average of the height at the perpendicularly eroded point, and the height at ys(i)
    !                ELSE
    !                    GOTO 1415
    !                END IF	
    !            END DO

    !        ELSE !slopes >0
    !            DO j=1, max(a-i,1)
    !                IF(ecu > ys(i+j)-ys(i)) THEN			
    !                        !print*, 'undercut right ', ecu
    !                    h_rhs(i+j)= min( ( (ys(i+j)-ys(i))*(bedlast_tmp(i)-ecu/slopes(i)+ Qd(i)*mor*dT/(1._dp-voidf)) +&
    !                    (ys(i)+ecu -ys(i+j) )*bed(i))/ecu, h_rhs(i+j))
    !                ELSE
    !                    GOTO 1415
    !                END IF		
    !                !bed(i+1)=bed(i)!max(bedlast_tmp(i), bed(i))
    !            END DO
    !        END IF !slopes > 0 

    !        1415 CONTINUE 
    !    END DO
    !    bed=h_rhs
    !END IF !Switch on/off undercutting


    !IF(.FALSE.) THEN
    !    DO i=1, a
    !        DO j= 1, layers
    !            IF(taucrit_dep(i,j)>bedlast_tmp(i)) THEN
    !                print*, 'taucrit_dep prob', i, j, bedlast_tmp(i), taucrit_dep(i,j), counter
    !            END IF

    !            IF((taucrit_dep(i,j)> bedlast_tmp(i)- dt*mor*Qe(i)/(1._dp-voidf)).AND.(tau(i)<taucrit(i,j))) THEN
    !                PRINT*, 'problem with taucrit_dep sinking', i, tau(i), j, taucrit(i,:), "and more", taucrit_dep(i,j), &
    !                bedlast_tmp(i),dt*mor*Qe(i)/(1._dp-voidf)
    !                STOP
    !            END IF
    !            !!!THIS NEEDS FURTHER THOUGHT IF it will be used for bedload
    !            taucrit_dep(i,j)= min(taucrit_dep(i,j), bedlast_tmp(i)- dt*mor*Qe(i)/(1._dp-voidf))!, bed(i)-Qd(i)*mor*dt/(1._dp-voidf) ) !Note that we only consider old height - total erosion. Deposition (either suspended or through bedload) is assumed to add new sediment. Note that this adjustment to taucrit_dep will only change it if the old height less erosion cut into the taucrit_dep layer.
    !        END DO
    !    END DO
    !END IF       


END SUBROUTINE update_bed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE roughmult(aa,rmu, vel, Q, A, width,tbst,depths, ys, f, vegdrag, counter, ysl, ysu, bedl, bedu,bed,water,g)
    !This is for calculating the constant that we should multiply the friction
    !slope in the 1D model by, in order to account for the lateral distribution of
    !velocity. So it basicaly gives a roundness coefficient, although actually,
    !we get the (roughness coefficient/ depth) for numerical stability reasons
    INTEGER, INTENT(IN)::aa, counter
    REAL(dp), INTENT(IN OUT):: rmu
    REAL(dp), INTENT(IN):: vel, Q, A, width,tbst, depths, ys, f, vegdrag, ysl, ysu, bedl, bedu, bed, water,g
    DIMENSION vel(aa),tbst(aa), depths(aa), ys(aa), f(aa), vegdrag(aa), bed(aa)
    
    INTEGER:: a2, i
    REAL(dp):: depths2(aa), rmutop, rmubot, rmulasts

    !Calculation requires |velocity| > 0
    IF((abs(Q)>1.0E-10_dp).AND.(maxval(abs(vel))>1.0E-10_dp)) THEN
        
        rmutop=0._dp
        ! Numerator integrated using trapezoidal rule
        DO i=2,aa-1
            rmutop= rmutop+(f(i)/8._dp*vel(i)**2*tbst(i) + vegdrag(i)*vel(i)**2*depths(i))*(ys(i+1)-ys(i-1))*0.5_dp
        END DO
        rmutop=rmutop+ (f(1)/8._dp*vel(1)**2*tbst(1) + vegdrag(1)*vel(1)**2*depths(1))*0.5_dp*(ys(2)-ys(1) &
             + 2._dp*(water-bed(1))/(bedl-bed(1))*(ys(1)-ysl))
        rmutop=rmutop+ (f(aa)/8._dp*vel(aa)**2*tbst(aa) + vegdrag(aa)*vel(aa)**2*depths(aa))*0.5_dp*(ys(aa)-ys(aa-1) & 
            + 2._dp*(water-bed(aa))/(bedu-bed(aa))*(ysu-ys(aa)))

        !Denominator 
        rmubot= (Q/A)**2*A*g  
        ! (roughness/depth)
        rmu= rmutop/rmubot
    ELSE
        rmu=sum(f)/(8._dp*g*sum(depths)) !The mean value of (f/8g)/depth
    END IF

    IF(rmu.eq.0._dp) THEN
        print*, "rmu is zero", abs(Q), aa, sum(f), maxval(vel), minval(vel)
    END IF

    IF(isnan(rmu)) THEN
            print*, "rmu is nan"
            print*, vel
            print*, "...."
            print*, tbst
            print*, "..."
            print*, ys
    END IF

END SUBROUTINE roughmult
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE intnuc(inuc,waters,ys,bottom, vels,b, wincrem, Q, A)
    !! This is for getting the nuc term-- we get the portion inside the derivative,
    ! (called inuc) and then differentiate it within the mcCormack scheme

    INTEGER, INTENT(IN):: b
    REAL(dp), INTENT(IN):: vels, waters,bottom,wincrem,Q,A, ys
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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE refit(ys,bed,a)
    ! Given a set of 'a' points (ys, bed), this subroutine tries to respace them
    ! over the curve (ys, bed).
    ! The idea is to have a high density of points where that is needed, and a
    ! low density of points where that is okay

    ! ys = vector of y values
    ! bed = vector of bed elevationa values
    ! a = length(bed)

    INTEGER, INTENT(IN):: a
    REAL(dp), INTENT (IN OUT)::ys(a), bed(a)

    INTEGER:: dstspl,i,j, bankl, bankr, num_bank_pts, tmp(1), n1, n2, num_pts(5)
    REAL(dp):: dsts(a),dsts2(a),eqspc(a),p1,newxs(a),newys(a),w2,w1, wdth,&
                av,mxch(a),b,newxs2(a), tmpR, res_pts(6), high_res_width, bank_old
    SAVE res_pts ! This will record the boundaries between zones of different resolutions
    DATA res_pts /6*0.0_dp/
    LOGICAL:: okay, reint
    
    
    !Check the input data
    DO i= 1, a
    IF ( (isnan(ys(i))).or.(isnan(bed(i))) ) THEN
            print*, "isnan before refit "
            stop
    END IF
    END DO

    ! Bed slope
    dsts(2:a-1) = ((bed(3:a)-bed(2:a-1))/(ys(3:a) - ys(1:a-2))*(ys(2:a-1) - ys(1:a-2)) + &
                   (bed(2:a-1)-bed(1:a-2))/(ys(2:a-1) - ys(1:a-2))*(ys(3:a) - ys(2:a-1)) )/ &
                   (ys(3:a) - ys(1:a-2))
    dsts(1) = 0._dp
    dsts(a) = 0._dp

    ! Find location of the maximum slope on the left half of the channel
    tmp = maxloc(abs(dsts(2:a/2))) + 1
    bankl=tmp(1)

    ! Same on the right half of the channel
    tmp = maxloc(abs(dsts(a/2:a))) + a/2 -1
    bankr=tmp(1)

    ! Check whether we need to remesh
    ! If the distance between 'bankl' and 'the value of bankl last time we remeshed'
    ! is small, and if this is also true for 'bankr', then we don't remesh.
    ! This is good because remeshing can introduce numerical
    ! diffusion/disturbance, which we don't need. 
    bank_old = 0.5_dp*(res_pts(3) + res_pts(2)) !Value of the left bank last time we remeshed
    IF(abs(ys(bankl) - bank_old)< 0.25_dp*(res_pts(3)-res_pts(2)) ) THEN
        ! Check right bank
        bank_old = 0.5_dp*(res_pts(4) + res_pts(5))
        IF(abs(ys(bankr) - bank_old)< 0.25_dp*(res_pts(5)-res_pts(4)) ) THEN
            PRINT*, 'No need to remesh'
            RETURN
        END IF
    END IF 


    ! Define the y value where the resolution (dy) will change The idea is to
    ! have high resolution between [res_pts(2), res_pts(3)], and [res_pts(4),
    ! res_pts(5)], and low resolution elsewhere
    high_res_width=min(10._dp,(ys(bankr)-ys(bankl))/4.0_dp, &
                        (ys(bankl)-ys(1))*0.9_dp,(ys(a) - ys(bankr))*0.9_dp)
    res_pts(1) = ys(1)
    res_pts(6) = ys(a)
    res_pts(2) = ys(bankl) - high_res_width
    res_pts(3) = ys(bankl) + high_res_width
    res_pts(4) = ys(bankr) - high_res_width
    res_pts(5) = ys(bankr) + high_res_width
   
 
    IF((res_pts(2)<0._dp).or.(res_pts(5)>ys(a))& 
        .or.(res_pts(4)<res_pts(3)) ) THEN
        PRINT*, 'ERROR - the bank region points are not ordered correctly. &
                 Need to recode the refit routine to make this more general'
        PRINT*, res_pts
        stop
    END IF 

    n1 = min(2*floor( a*(high_res_width/ys(a))*3.0_dp), floor(0.4_dp*a)) ! Total number of points in the two high res regions
    n2 = a-2*n1 ! Total number of points in the three low res regions.
    ! num_pts holds the number of points in each of the 5 regions
    num_pts(2) = n1
    num_pts(4) = n1
    
    print*, 'REMESHING' 
    print*, 'n1 = ', n1
    print*, 'n2 = ', n2
    print*, 'ys(bankl) = ', ys(bankl)
    print*, 'ys(bankr) = ', ys(bankr)
    print*, 'high_res_width = ', high_res_width
    print*, 'END REMESHING'

    ! Calculate number of points in the first low res region  =
    ! n2*first_region_width/total_lowres_region_width
    tmpR = n2*(res_pts(2)-res_pts(1))/& 
           (res_pts(2)-res_pts(1) + res_pts(4) - res_pts(3) + res_pts(6) - res_pts(5) )
    num_pts(1) = floor(tmpR)
    num_pts(5) = floor(tmpR)
    num_pts(3) = n2 - 2*floor(tmpR)   
    
    ! Define newx values, which will be used to re-interpolate over the
    ! cross-section

    ! Region 1
    newxs(1:num_pts(1)) = res_pts(1) + &
                          (res_pts(2) - res_pts(1))*(/ (i-1,i=1, num_pts(1)) /)/(1.0_dp*(num_pts(1)-1))  
    ! Region 2
    newxs((num_pts(1)+1):(num_pts(2)+num_pts(1))) = res_pts(2) + &
                                                (res_pts(3) - res_pts(2))*(/ (i,i=1, num_pts(2)) /)/(1.0_dp*(num_pts(2)+1))
    ! Region 3
    tmp(1) = num_pts(1) + num_pts(2)
    newxs((tmp(1)+1):(num_pts(3)+ tmp(1))) = res_pts(3) + &
                                         (res_pts(4) - res_pts(3))*(/ (i-1,i=1, num_pts(3)) /)/(1.0_dp*(num_pts(3)-1))
    ! Region 4
    tmp(1) = num_pts(1) + num_pts(2) + num_pts(3)
    newxs((tmp(1)+1):(num_pts(4)+tmp(1))) = res_pts(4) + &
                                        (res_pts(5) - res_pts(4))*(/ (i,i=1, num_pts(4)) /)/(1.0_dp*(num_pts(4)+1))
    ! Region 5
    tmp(1) = num_pts(1) + num_pts(2) + num_pts(3) + num_pts(4)
    newxs((tmp(1)+1):(num_pts(5)+ tmp(1))) = res_pts(5) + &
                                        (res_pts(6) - res_pts(5))*(/ (i-1,i=1, num_pts(5)) /)/(1.0_dp*(num_pts(5)-1))
    
      
    ! Get the bed values with cubic interpolation
    call interp3(ys, bed, newxs, a)

    ys(2:a-1)=newxs(2:a-1)

    !Check that we didn't have problems
    DO i= 1, a
    IF ( (isnan(ys(i))).or.(isnan(bed(i))) ) THEN
            print*, "isnan after refit"
            stop
    END IF
    END DO


end subroutine refit
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE geotech(ys, bed, slopes, taucrit_dep_ys,taucrit_dep, dst, nos,a, mu,layers, slpmx)
    ! A defunct routine for ensuring that the lateral bed slope does not become
    ! too large
    INTEGER, INTENT(IN):: nos, a,layers
    REAL(dp), INTENT(IN):: mu, ys(a),taucrit_dep_ys(nos),taucrit_dep(nos,layers),slpmx(nos,0:layers+1)
    REAL(dp), INTENT(IN OUT):: dst(a,0:layers+1), bed(a), slopes(a)
    !Integer, intent(in out):: indd(a,layers)

    INTEGER::i, lyr, lyr1(1)
    REAL(dp):: lcslp(a), bedtrial, hfirst(nos)


    !call findten(ys(floor(a/2.)+1), bed(floor(a/2.)+1), slopes(floor(a/2.)+1), taucrit_dep_ys, & 
    !taucrit_dep(:,1), dst(floor(a/2.)+1,1), nos, indd(floor(a/2.)+1,1), .false.) 

    !print*, dst(:, 1)
    hfirst=bed

    DO i= floor(a/2.), 2, -1 !March back from the center checking for slope failure.

        IF(dst(i, 1)>0._dp) THEN  !We are in the first layer
            lyr= 0
        ELSE
            lyr1=minloc(dst(i,:), dst(i,:)>0._dp)  !This is the layer we are in.
            lyr=lyr1(1)
            !print*, lyr, lyr1(1)
        END IF

        214 lcslp(i)= (bed(i)-bed(i-1))/(ys(i)-ys(i-1)) ! The local slope. Notice how this is a less conservative slope estimate than is the weighted central difference. It should take larger values than the central difference, thus keeping the central difference below mu. 
        IF(abs(lcslp(i)).GE. slpmx(i,lyr) ) THEN ! The slope limit is exceeded, and so erosion occurs at the higher point. It erodes either until the cohesive layer is reached, or until the local slope is mu
            !print*, "changing 1", i, lcslp(i), bed(i), bed(i-1)
            !Here the first if statement checks whether this point or the one above
            !should be eroded -- one of them must move if the slope is to
            !decrease. Normally for this half of the channel, bed(i)>bed(i+1) 
            IF((bed(i)>bed(i-1))) THEN
                bedtrial=  bed(i-1)+slpmx(i,lyr)*( ys(i)-ys(i-1) ) !This is the height value possible for this layer. It will be accepted if we haven't cut into a new layer in the process
                IF(bedtrial>bed(i)-dst(i,lyr+1)) THEN !Accept
                    bed(i)= bedtrial
                ELSE !We need to move down to the next layer.
                    bed(i)= bed(i)-dst(i,lyr+1)
                    lyr=lyr+1
                    GOTO 214
                END IF

            ELSE
                bedtrial=  bed(i)+slpmx(i,lyr)*(ys(i)-ys(i-1) )
                IF(bedtrial>bed(i-1)-dst(i-1, lyr+1)) THEN !Accept
                    bed(i-1) = bedtrial
                ELSE !We need to move down a layer
                    bed(i-1) = bed(i-1)-dst(i-1, lyr+1)
                    lyr=lyr+1
                    GOTO 214
                END IF

            END IF

             ! print*, hs(i), hs(i-1)
        END IF
    END DO

    !call findten(ys(ceiling(a/2.)), bed(ceiling(a/2.)), slopes(ceiling(a/2.)), taucrit_dep_ys, taucrit_dep(:,1), & 
    !dst(ceiling(a/2.),1), nos, indd(ceiling(a/2.),1), .false.) 

    DO i= ceiling(a/2.), a-1 !March back from the center checking for slope failure.
        !indd(i+1,1)= indd(i,1)-5 !Predefine this at a 'lower likely bound' to shorten the findten search. Could be bad if the defined value is not a lower bound
        !call findten(ys(i+1), bed(i+1), slopes(i+1), taucrit_dep_ys, taucrit_dep(:,1), dst(i+1,1), nos, indd(i+1,1), .true.) 

        IF(dst(i, 1)>0._dp) THEN  !We are in the first layer
            lyr= 0
        ELSE
            lyr1=minloc(dst(i,:), dst(i,:)>0._dp)  !This is the layer we are in.
            lyr=lyr1(1)
        END IF

        314 lcslp(i)= (bed(i+1)-bed(i))/(ys(i+1)-ys(i)) !Notice how this is a less conservative slope estimate than is the weighted central difference. It should take larger values than the central difference, thus keeping the central difference below mu. 
        IF(abs(lcslp(i)).GE. slpmx(i,lyr) ) THEN ! The slope limit is exceeded, and so erosion occurs at the higher point. It erodes either until the cohesive layer is reached, or until the local slope is mu
        !    print*, "changing 2", i, lcslp(i), bed(i), bed(i+1)
            !Here the first if statement checks whether this point or the one above
            !should be eroded -- one of them must move if the slope is to
            !decrease. Normally for this half of the channel, bed(i)<bed(i+1) 
            IF((bed(i)>bed(i+1))) THEN
                bedtrial=  bed(i+1)+slpmx(i,lyr)*(ys(i+1)-ys(i) ) !This is the height value possible for this layer. It will be accepted if we haven't cut into a new layer in the process
                If(bedtrial>bed(i)-dst(i,lyr+1)) THEN !Accept
                    bed(i)= bedtrial
                ELSE !We need to move down to the next layer.
                    bed(i)= bed(i)-dst(i,lyr+1)
                    lyr=lyr+1
                    GOTO 314
                END IF

            ELSE
                bedtrial=  bed(i)+slpmx(i,lyr)*(ys(i+1)-ys(i) )
                IF(bedtrial>bed(i+1)-dst(i+1, lyr+1)) THEN !Accept
                    bed(i+1) = bedtrial
                ELSE !We need to move down a layer
                    bed(i+1) = bed(i+1)-dst(i+1, lyr+1)
                    lyr=lyr+1
                    GOTO 314
                END IF

            END IF
        END IF
    END DO
    !Recalculate the slopes
    slopes(2:a-1)= ((bed(3:a)-bed(2:a-1))/(ys(3:a)-ys(2:a-1))*(ys(2:a-1)-ys(1:a-2)) + &
    (bed(2:a-1)-bed(1:a-2))/(ys(2:a-1)-ys(1:a-2))*(ys(3:a)-ys(2:a-1) ) ) & 
    /(ys(3:a)-ys(1:a-2)) 

    slopes(1)= (bed(2)-bed(1))/(ys(2)-ys(1))
    slopes(a)= (bed(a)-bed(a-1))/(ys(a)-ys(a-1))

    DO i= 1, a
        IF(isnan(bed(i))) THEN
            PRINT*, "bed nan in geotech", bed(i-1:i+1), slopes(i-1: i+1)
            STOP
        END IF
    END DO

END SUBROUTINE geotech
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interp(oldxs, oldys, newxs, a)
    !!!Linear interpolation of the curve oldxs,oldys at the x values newxs - output
    !the new y values in the oldys vector.
    INTEGER, INTENT(IN) :: a
    REAL(dp), INTENT(IN) :: oldxs(a), newxs(a)
    REAL(dp), INTENT(IN OUT):: oldys(a)

    INTEGER:: i, ind, lwer, j
    REAL(dp):: newys(a), tester(a)


    lwer= 1
    DO i=2, a-1
        DO j= lwer, a-1
            IF( (newxs(i)>=oldxs(j)).AND.(newxs(i)<=oldxs(j+1))) THEN
                newys(i)= ( oldys(j+1)*(newxs(i)-oldxs(j))+ oldys(j)*(oldxs(j+1)-newxs(i)))/(oldxs(j+1)-oldxs(j)) !The linearly interpolated value of newys(i) at newxs(i)
                lwer=j
                goto 2232
            END IF
        END DO
        2232 CONTINUE
    END DO

    !FIXME: POTENTIAL BUG: Here we test for symmetry in oldys. If it is symmetric, then we assume that newys should also be symmetric. Now, technically, this might not always be true. However, in my situations, I think it will always be true. And loss of roundoff is a problem. So here we go
    tester=oldys(1:a)-oldys(a:1:-1)
    IF(maxval(abs(tester))<1.0E-9) THEN
            !Enforce symmetry
            DO i=1,floor(a*0.5_dp)
            newys(i)=0.5_dp*(newys(i)+newys(a-i+1))
            newys(a-i+1)=newys(i)
            END DO
    END IF

    oldys(2:a-1)=newys(2:a-1) !Update the y values

END SUBROUTINE interp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE interp3(oldxs, oldys, newxs, a)
    !!!Cubic interpolation of the curve oldxs,oldys at the x values newxs - output
    !the new y values in the oldys vector.
    INTEGER, INTENT(IN) :: a
    REAL(dp), INTENT(IN) :: oldxs(a), newxs(a)
    REAL(dp), INTENT(IN OUT):: oldys(a)

    INTEGER:: i, ind, lwer, j, ierr
    LOGICAL:: skip=.true.
    REAL(dp):: newys(a), tester(a), deriv(a), wk(4*a)



    !Here we call a slatec routine which calculates derivatives as needed for
    !monotonic (locally) hermitite cubic spline interpolation
    call DPCHIC((/0,0/),(/0,0/),0, a,oldxs,oldys,deriv,1,wk, 4*a, ierr)

    !call DPCHIM(a,oldxs,oldys,deriv,1, ierr) !Interestingly, this seems to cause an error

    !Check Error flag
    IF(ierr.NE.0) THEN
        print*, "Error in spline derivs: Flag =", ierr
        stop
    END IF

    !Actually do the interpolation
    call DPCHFE(a,oldxs,oldys,deriv,1,skip,a,newxs,newys,ierr)


    !lwer= 1
    !do i=2, a-1
    !
    !        do j= lwer, a-1
    !        IF( (newxs(i)>=oldxs(j)).and.(newxs(i)<=oldxs(j+1))) THEN
    !                newys(i)= ( oldys(j+1)*(newxs(i)-oldxs(j))+ oldys(j)*(oldxs(j+1)-newxs(i)))/(oldxs(j+1)-oldxs(j)) !The linearly interpolated value of newys(i) at newxs(i)
    !             
    !              lwer=j
    !                goto 2232
    !        end if
    !        end do
    !        2232 continue
    !
    !!end do

    !FIXME: POTENTIAL BUG: Here we test for symmetry in oldys. If it is symmetric, then we assume that newys should also be symmetric. Now, technically, this might not always be true. However, in my situations, I think it will always be true. And loss of roundoff is a problem. So here we go
    IF(.true.)THEN
        tester=oldys(1:a)-oldys(a:1:-1)
        IF(maxval(abs(tester))<1.0E-9) THEN
            !Enforce symmetry
            DO i=1,floor(a*0.5_dp)
            newys(i)=0.5_dp*(newys(i)+newys(a-i+1))
            newys(a-i+1)=newys(i)
            END DO
        END IF
    END IF

    oldys(2:a-1)=newys(2:a-1) !Update the y values

end subroutine interp3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Ddx_3E(n2, Q,delX,dQdx)
    !Subroutine to calculate the delQ/delX, using an explicit third order method
    !with a flux limiter. It is assumed that Q (the raw flux at point i) contains
    !two boundary values to the left, and two to the right.
    INTEGER, INTENT(IN) :: n2
    REAL(dp), INTENT(IN):: delX,Q((-1):n2+2) !Notice how Q has 2 boundary points to the left and right
    REAL(dp), INTENT(INOUT):: dQdx(n2) !Notice how dQdx doesn't have the boundary points

    INTEGER:: i
    REAL(dp):: limi(0:n2), theta(0:n2+1), Q1(0:n2), mu
    !limi is the flux limiter
    !theta is a parameter for the flux limiter
    !mu us another parameter for the flux limiter
    !Q1 is the 2nd order flux at i+1/2

    mu=1._dp

    !Define limiter, see Hundsdorfer and Verwer
    DO i=0,n2+1
        theta(i)=(Q(i)-Q(i-1))/(Q(i+1)-Q(i))  

        !Catch special cases - theta = nan or 0
        IF(isnan(theta(i))) THEN
            IF(Q(i).ne.Q(i-1)) THEN
                theta(i)=1.0E+12_dp*sign(1._dp, Q(i)-Q(i-1))
            ELSE
                theta(i)=1.0_dp
            END IF
        END IF

        IF(theta(i).eq.0._dp) THEN
            IF(Q(i+1).ne.(Q(i))) THEN
                theta(i)=1.0E-12_dp*sign(1._dp, Q(i+1)-Q(i))
            ELSE
                theta(i)=1.0_dp 
            END IF

        END IF

    END DO

    !Calculate fluxes, with limiting -- 3rd order without limiting
    DO i=0, n2
        IF(0.5_dp*(Q(i)+Q(i+1))>=0._dp) THEN
            limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/6._dp*theta(i) , mu*theta(i)))
            !if(min(abs(Q(i)),abs(Q(i+1)) ).eq.0._dp) limi(i)=0._dp
            !limi(i)= 0.5_dp*(theta(i)+abs(theta(i)))/(1._dp+abs(theta(i)))
            Q1(i)= Q(i)+limi(i)*(Q(i+1)-Q(i)) 
        ELSE
            limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/(6._dp*theta(i+1)) , mu/theta(i+1)))
            !if(min(abs(Q(i)),abs(Q(i+1)) ).eq.0._dp) limi(i)=0._dp
            !limi(i)= 0.5_dp*(1._dp/theta(i+1)+abs(1._dp/theta(i+1)))/(1._dp+abs(1._dp/theta(i+1)))
            Q1(i)= ( Q(i+1)+limi(i)*(Q(i)-Q(i+1)) ) 
        END IF
    END DO

    dQdx=(Q1(1:n2)-Q1(0:n2-1))/delX

    DO i=1,n2
        IF(isnan(dQdx(i))) THEN
            PRINT*, 'dQdx', i, 'is nan'
            PRINT*, Q((i-2):(i+2))
            PRINT*, Q1(max((i-2),0):min((i+2),n2))
            PRINT*, limi( (i-1):(i+1) )
            PRINT*, theta((i-1):(i+1))
        END IF
    END DO

END SUBROUTINE Ddx_3E
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE conservation_tests1(l1,Anew,Aold,delX,Q2,DT,Q1in,Vol1)
    INTEGER, INTENT(IN):: l1 !Length of Anew
    REAL(dp), INTENT(IN):: DT !Time step
    REAL(dp), INTENT(IN):: Anew(l1),Aold(l1), Q2(l1+1), delX !Hydrodynamic variables
    REAL(dp), INTENT(IN OUT):: Q1in, Vol1 !Variables in which we store the conservation metrics

    Q1in=Q1in+ (Q2(1) - Q2(l1+1))*DT !The total discharge added to the reach
    Vol1=sum(Anew)*delX !The volume of water in the reach


END SUBROUTINE conservation_tests1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE conservation_tests2(l1,Anew,wetwidth,delX,Q2,DT,C,C_old,FL1_mouth, FL1_river, E, wset, rhos, QS2in, VolS2, Source2)
    INTEGER, INTENT(IN):: l1 !Length of Anew
    REAL(dp), INTENT(IN):: DT !Time step
    REAL(dp), INTENT(IN):: Anew(l1),wetwidth(l1), Q2(l1), C(l1),C_old(l1), E(l1), delX,wset,rhos, FL1_mouth, FL1_river !Hydrodynamic/sediment variables
    REAL(dp), INTENT(IN OUT):: QS2in, VolS2,Source2 !Variables in which we store the conservation metrics

    !IF(Q2(1)<0._dp) THEN
    QS2in=QS2in+ (FL1_mouth -  FL1_river)*DT !The total sediment mass discharge added to the reach
    !ELSE
    !QS2in=QS2in+ (Q2(1)*Cmouth -  Q2(l1)*C(l1))*DT !The total sediment mass discharge added to the reach
    !END IF
    VolS2=sum(Anew*C)*delX !The mass of sediment in the reach
    Source2= Source2+ sum(E*rhos-min(wset*wetwidth,Anew/DT)*0.5*(C+C_old) )*DT*delX !The mass eroded/deposited

END SUBROUTINE conservation_tests2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE qbh_approx(n,ys,qb,qbnP1,bed, bedl, bedu, ysl, ysu, order)
    !On input, qb has the values of some function qb at ys(i).  
    !On output, qb has the values of qb AT ys(i+1/2), EVALUATED BY FITTING A
    !PARABOLA THROUGH 3 POINTS NEIGHBOURING i+1/2, WITH THE DIRECTION DETERMINED
    !BY SLOPE BASED UPWINDING.  
    !Note the use of the boundary condition qbnP1 (=qb at N+1) which we need to
    !do this
    INTEGER, INTENT(IN):: n, order !Dimension of the other vectors, order of polynomial used for interpolation (1 or 2)
    REAL(dp),INTENT(IN):: ys(n), bed(n), bedu, bedl, ysu, ysl, qbnP1 ! y coordinate and elevation, and the values of bed and ys beyond the left and right boundaries, and the single point value of qb and n+1
    REAL(dp),INTENT(IN OUT):: qb(0:n) !qb. 

    REAL(dp):: qbt(0:n), bedt(0:n+1),yst(0:n+1) !Temporary qb(i), bed(i) and ys(i)
    REAL(dp):: delyP(0:n), delyM(1:n+1), delqbP(0:n), delqbM(1:n+1) !Forward increment of ys, Backward of ys, forward of qb, backward of qb
    REAL(dp):: a1, a2,a0, useme, qb_old(0:n) !Constants.
    INTEGER:: i, ii
    INTEGER, save:: counter=0

    IF(order.ne.1) THEN
        IF(order.ne.2) THEN
            print*, 'Order (in qbh_approx) has an incorrect value'
            stop
        END IF
    END IF
        

    !Temporary variables
    qbt=qb
    qb_old=qb
    bedt(1:n)=bed
    yst(1:n)=ys

    bedt(0)=bedl
    bedt(n+1)=bedu

    yst(0)=ysl
    yst(n+1)=ysu

    !Forward and backward increments of ys
    delyP(0:n)=yst(1:n+1)-yst(0:n)
    delyM(1:n+1)=delyP(0:n)

    !Forward and backward increments of qb
    delqbP(0:n-1)=qb(1:n)-qb(0:n-1)
    delqbP(n)=qbnP1-qb(n) !Use boundary condition
    delqbM(1:n+1)=delqbP(0:n)

    !Calculate qb(i+1/2)
    !We need to choose either a backward, forward or central extrapolation. Here
    !we create a flag which = [i for backward], = [i+1 for forward], and = [-1 for
    !central]
    DO i=0,n

        !IF(order.eq.2) THEN
        SELECT CASE(order)
            ! Order ==2
            CASE(2)
                IF(i.eq.n) THEN
                    ii=i !-1 !i  
                ELSE
                    IF(i.eq.0) THEN
                        ii=i+1 !-1 !i+1  
                    ELSE
                        IF((bedt(i)>bedt(i+1)).and.(bedt(i-1)>bedt(i))) THEN
                            ii=i !Use backward extrap (points i-1,i and i+1)
                        ELSE
                            IF((bedt(i)<bedt(i+1)).and.(bedt(i-1)<bedt(i))) THEN
                                ii= i+1 !Use forward extrap (using points i, i+1, i+2)
                            ELSE
                                ii=-1 !Use central extrap
                            END IF
                        END IF
                    END IF
                END IF
            ! Order ==1
            CASE(1)
                ii = -1 !Make everything a central extrapolation
            
            CASE DEFAULT
                PRINT*, 'ERROR: the variable order should be either 1 or 2 in qbh_approx'
                stop
        END SELECT


        !Here we do the extrapolation        
        IF(ii==-1) THEN !Central approach
            IF(i<n) THEN
                qbt(i)=0.5_dp*(qb(i)+qb(i+1))
            ELSE
                qbt(i)=0.5_dp*(qb(i)+qbnP1)
            END IF
        ELSE !Use either forward or backward extrap
            !Say the polynomial is F(x)= c+b*(x-x(i)) + a*(x-x(i))*(x-x(i+1))
            !Then c=y(i), b=(y(i+1)-y(i))/(x(i+1)-x(i)), a= ( (y(i-1)-y(i))/(x(i-1)-x(i)) -b)/(x(i-1)-x(i+1))
            a0=qb(ii) !c
            a1=delqbP(ii)/delyP(ii) !b
            a2=(-delqbM(ii)/delyM(ii) + a1)/(delyM(ii)+delyP(ii)) !a
            
            IF(ii==i) THEN
                qbt(i)= a0+(delyP(ii)*0.5_dp)*(a1+a2*(-delyM(ii+1)*0.5_dp)) !The estimate of qb(i+1/2)
                !print*, '-', i
            ELSE
                qbt(i)= a0+(-delyM(ii)*0.5_dp)*(a1+a2*(-delyM(ii)*0.5_dp -delyP(ii)))  !The estimate of qb(i+1/2)
                !print*, '+', i
            END IF
            !qbt(i)= min(qbt(i),0._dp) !Not monotonic
            
        END IF
    END DO
    
    DO i=1,n-1
        !Check for monotonicity, and enforce if needed.
        IF((qbt(i)-qb(i))*(qbt(i)-qb(i+1)) >0._dp ) THEN
            qbt(i) = 0.5_dp*(qb(i)+qb(i+1))
        END IF
    END DO

    qb(0:n)=qbt(0:n) !The estimate of qb(i+1/2)
    !print*, 'qber', maxval(abs(qbt(1:n-1)-qbt(n-1:1:-1))), maxval(abs(qbt(1:n-1)))

    !qb_old(0:n)=0.5_dp*((/qb_old(1:n),qbnP1/) + qb_old(0:n))
    !print*,'qbM', maxval(abs(qb(0:n)-qbt(0:n)))/maxval(abs(qbt))

    !counter=counter+1
    !IF((mod(counter,1).EQ.0).and.(maxval(abs(qbt(1:n-1)-qbt(n-1:1:-1)))>0.0_dp)) THEN
    !DO i=0,n
    !WRITE(13,*) -n, qb(0:n),-n,qb_old(0:n), -n
    !END DO
    !END IF

END SUBROUTINE qbh_approx
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dbeddyH_approx(a,ys,bed, dbeddyH, ysl, ysu, bedl, bedu, order)
! This routine takes ys, bed, and calculates the coefficients b1, b2, b3, b4,
! where dhdy_(i+1/2) = b1*bed(i-1) + b2*bed(i) + b3*bed(i+1) + b4*bed(i+2)
! These are stored in dbeddyH
! The coefficients come from fitting a curve to the points, which is either a
! cubic (order =3), or a linear fit (order =2 -- note that a parabola gives the
! same answer though)
    INTEGER, INTENT(IN)::a, order
    REAL(dp), INTENT(IN):: ys(a), bed(a), ysl, ysu, bedl, bedu
    REAL(dp),INTENT(IN OUT):: dbeddyH(4,0:a)

    INTEGER:: i, ii
    REAL(dp):: xtemp, ys_temp(0:a+1), bed_temp(0:a+1), sl1, sl2, limiter
    REAL(dp):: s1, s2, mono_test, dy, dyinv
    
    !IF(order.ne.2) THEN
    !    IF(order.ne.3) THEN
    !        print*, 'Order (in dbeddyH_approx) has an incorrect value'
    !        stop
    !    END IF
    !END IF

    !Define a convenient ys_temp variable including boundaries
    ys_temp(1:a)=ys
    ys_temp(0)=ysl
    ys_temp(a+1)=ysu

    bed_temp(1:a) = bed
    bed_temp(0)=bedl
    bed_temp(a+1)=bedu

    ! Pre-zero matrix that will store coefficients
    dbeddyH=0._dp

    ! Useful constants
    dy = ys_temp(2)-ys_temp(1)
    dyinv = 1.0_dp/dy

    IF(order==3) THEN
        DO i=1,a+1
            IF( ((ys_temp(i) - ys_temp(i-1))-dy) > 1.0e-8_dp*dy) THEN
                print*, 'ERROR: dy does not appear constant, but the 3rd &
                        order approximation is being used in dbeddyH_approx. &
                        Perhaps try the linear (2nd order) approximation instead' 
                stop
            END IF
        END DO
    END IF


    DO i=0,a
        SELECT CASE (order)
            ! Order ==3
            CASE(3)
                ! Treat boundaries first, then interior points
                IF((i.eq.a).or.(i.eq.a-1)) THEN
                    ii = -1 !Use central extrapolation

                ELSEIF((i.eq.0).or.(i.eq.1)) THEN
                    ! Central extrapolation
                        ii = -1
               ELSE
                    ! Interior points
                    IF( (bed_temp(i)>bed_temp(i+1)).and.(bed_temp(i-1)>bed_temp(i))) THEN
                        ii= 1
                    ELSEIF( (bed_temp(i)<bed_temp(i+1)).and.(bed_temp(i-1)<bed_temp(i))) THEN
                        ii = 1
                    ELSE
                        ii=-1 !Use central extrap
                    END IF
                END IF
            ! Order ==1
            CASE(2)
                ii = -1 ! Make everything central extrap
        
            CASE DEFAULT
                print*, 'ERROR: order should be either 3 or 2 in dbedHdy_approx'
                stop
        END SELECT

        
        
        ! Now compute the coefficients for the extrapolation
        SELECT CASE(ii)
            Case(1)
                ! Note that this case requires a constant dy -- as is tested
                ! above
                dbeddyH(1,i) = 1.0_dp/24._dp*dyinv
                dbeddyH(2,i) = -9.0_dp/8.0_dp*dyinv
                dbeddyH(3,i) = 9.0_dp/8.0_dp*dyinv
                dbeddyH(4,i) = -1.0_dp/24._dp*dyinv
                
                ! Check that the result has the same sign as a central slope --
                ! if not, then there are oscillations in the interpolating
                ! polynomial, and we use the central estimate. 
                ! mono_test = central_slope * cubic_slope
                mono_test = (bed_temp(i+1)-bed_temp(i-1))* & 
                        (-bed_temp(i+2) + 27.0_dp*bed_temp(i+1) - 27.0_dp*bed_temp(i) + bed_temp(i-2))
                IF(mono_test<=0._dp) THEN
                    dbeddyH(1,i) = 0._dp
                    dbeddyH(2,i) = -dyinv
                    dbeddyH(3,i) = dyinv
                    dbeddyH(4,i) = 0._dp
                END IF
            Case(-1) 
            !Central estimate - easy
                IF((i<a).AND.(i>0)) THEN
                    dbeddyH(3,i)=1._dp/(ys_temp(i+1)-ys_temp(i))
                    dbeddyH(2,i)=-dbeddyH(3,i)
                ELSE
                    IF(i==a) THEN
                        dbeddyH(3,i)=1._dp/(ysu-ys_temp(i))
                        dbeddyH(2,i)=-dbeddyH(3,i)
                    ELSE !i==0
                        dbeddyH(3,i)=1._dp/(ys_temp(i+1)-ysl)
                        dbeddyH(2,i)=-dbeddyH(3,i)
                    END IF
                END IF

            CASE DEFAULT
                print*, 'ERROR: ii should be +-1 in dbeddyH_approx'
                stop

        END SELECT
    END DO
    

END SUBROUTINE dbeddyH_approx
!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE active_zone(a,height1,height2,morbl,morbu, offset)
    !This subroutine gives the lower and upper indicies of the zone where height1 is
    !different to height2. Then, an 'offset' is applied to these boundaries - so they are extended outward by 'offset' points.
    INTEGER, INTENT(IN)::a, offset
    INTEGER, INTENT(INOUT):: morbl, morbu
    REAL(dp), INTENT(INOUT):: height1(a),height2(a)

    INTEGER:: i

    morbl=0
    morbu=0

    DO i=1,a
        IF((height1(i)-height2(i)).ne.0._dp) THEN
            !print*, 'morbl', i,height1(i)-height2(i)
            morbl=i !max(i-6,1)
            goto 2345
        END IF
    END DO
    2345 CONTINUE

    DO i=a,1,-1
        IF((height1(i)-height2(i)).ne.0._dp) THEN
            !print*, 'morbu', i,height1(i)-height2(i)
            morbu=i !min(i,a)
            GOTO 2346
        END IF
    END DO
    2346 CONTINUE

    !Add an offset if desired
    if(morbl>0) morbl=max(morbl-offset,1)
    if(morbu>0) morbu=min(morbu+offset,a)


END SUBROUTINE active_zone
!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE reset_ys(a,ys,morbl,morbu, minact,maxnoact)
    !This subroutine re-defines ys so that most sufficient points are located
    !between ys(morbl) and ys(morbu), which is the morhologically active
    !zone
    INTEGER, INTENT(IN)::a
    INTEGER, INTENT(INOUT):: morbl, morbu
    REAL(dp), INTENT(IN):: minact, maxnoact !The minimum spacing between points in the active zone will be (average spacing)*minact, while the max spacing between points out of the active zone will be (average spacing)*maxnoact 
    REAL(dp), INTENT(INOUT):: ys(a)

    INTEGER:: i,pts_act,pts_noact
    REAL(dp):: active_w, full_w, av_space, ltmp(a)

    active_w=ys(morbu)-ys(morbl) !Width of the active zone
    full_w=ys(a)-ys(1) !Width of the cross-section

    av_space=full_w/(1._dp*(a-1)) !The average distance between x points in the cross-section
    pts_act= ceiling(active_w/(av_space*minact))+1 !Maximum desired number of points in the active zone
    pts_noact=ceiling( (full_w-active_w)/(av_space*maxnoact)) !Minimum desired number of points outside the active zone

    !Ensure that the latter is even
    IF(mod(pts_noact,2).ne.0) pts_noact=pts_noact+1


    IF(pts_act+pts_noact < a) THEN
        !So here we have the max number of points in the active zone, and we put the rest outside the active zone
        pts_noact=a-pts_act
        !Ensure that pts_noact is even - because for symmetry, we should have equal numbers of points on the left and the right
        IF(mod(pts_noact,2).NE.0) THEN
            pts_noact=pts_noact+1
            pts_act=pts_act-1
        END IF       
    ELSE
        !So here we possibly do not have enough points to have the max number in the active zone
        pts_act=a- pts_noact
    END IF

    !Redefine the locations of ys in a temporary vector ltmp. non-active zone
    IF(pts_noact>1) THEN
        DO i=1,(pts_noact/2)
            ltmp(i)=ys(1)+(1._dp*(i-1))*(ys(morbl)-ys(1))/(1._dp*(pts_noact/2))
            ltmp(a-i+1)=ys(a)+(1._dp*(i-1))*(ys(morbu)-ys(a))/(1._dp*(pts_noact/2))
        END DO
    END IF
    !As above for the active zone, which includes ys(morbl) and ys(morbu)
    DO i=(pts_noact/2+1), a-(pts_noact/2)
        ltmp(i)= ys(morbl) +(1._dp*(i-(pts_noact/2+1)))*(ys(morbu)-ys(morbl))/(1._dp*(a-pts_noact-1))
    END DO

    !Check if symmetry is intended, by checking for near symmetry in ys. If so, we want to retain symmetry, and so try to supress round-off error
    IF((maxval(ys(1:a)+ys(a:1:-1))).le.(minval(ys(1:a)+ys(a:1:-1))+1.0E-12_dp)) THEN !If this is satisfied, then symmetry is clearly intended
        ltmp=0.5_dp*( (ltmp(1:a)-ltmp(1))+(ltmp(a)-ltmp(a:1:-1))) +ltmp(1) !Try to reduce round-off symm break
    END IF
    !Define the new ys, morbu,and morbl
    ys=ltmp
    morbl=pts_noact/2+1
    morbu=a-pts_noact/2



END SUBROUTINE reset_ys
!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE sus_dist_sect(a, ys, bed, water, tau, vel, sconc, wset,lambdacon,rho,rhos, & 
!                            g, d50, susQbal, Qbed,talmon, sllength, bedl,bedu, ysl, ysu, C, integrated_load_flux)
!    INTEGER, INTENT(IN)::a
!    REAL(dp), INTENT(IN):: ys, bed, water, tau, vel, sconc, wset, lambdacon, rho, rhos,g, d50, Qbed, &
!                           sllength, bedl, bedu,ysl,ysu, integrated_load_flux
!    LOGICAL, INTENT(IN):: susQbal, talmon
!    REAL(dp), INTENT(OUT):: C
!    DIMENSION ys(a), bed(a), tau(a),vel(a), Qbed(a), sllength(a), C(a)
!    !PURPOSE: Calculate the distribution of suspended load over the cross-section C (kg/m^3)
!    !OUTCOME: C is calculated
!
!    ! Solves 
!    !        (negative of lateral flux)= 0 
!    ! or 
!    !        (negative of lateral flux) = qbl
!    !
!    ! d(ZETA)/dy - cbed*d(Depth)/dy = 0  (or = qbl) 
!    ! 
!    ! ZETA is the depth integrated sediment concentration, and
!    ! cbed is the near bed suspended sediment concentration, and
!    ! the suspended sediment is distributed exponentially in the vertical with
!    ! ZETA= edvis_z/wset*(1._dp-exp(-(wset/edvis_z)*(water-bed)) )*cbed, and
!    ! qbl is the lateral bedload transport
!    ! Note that there will be a family of solutions - the homogeneous part of the equation is linear - we
!    ! need to avoid getting the zero solution - do this by forcing a central value
!    ! of cbed to be what we want
!    !
!    ! NOTES BELOW MAY BE OUTDATED NOW
!    ! Note also that I have had problems with points to the edge of the cross-section
!    ! that are 'cut-off' from the main channel (i.e. there are points of zero depth
!    ! between the main channel and them). At such locations, this method was
!    ! producing large spikes in the sediment concentration. How to get around this?
!    ! One way is to note that the equation can be described purely in terms of depth,
!    ! by dividing by dDepth/dy.  Thus, C should be a function of depth. So, at points
!    ! outside the main channel where C has spiked, we could replace their value with
!    ! an extrapolated value using the main channel C's. For example, find the min
!    ! depth in the main channel, and linearly extrapolate to zero - use this to
!    ! reassign C to each bit.
!
!    !weighted average of derivative
!    !dZeta/dy = [ (y(i)-y(i-1))*(zeta(i+1)-zeta(i))/(y(i+1)-y(i)) +  (y(i+1)-y(i))*(zeta(i)-zeta(i-1))/(y(i)-y(i-1))  ]/(y(i+1)-y(i-1))
!
!    REAL(dp):: hss2f(0:a+1), edvisy(a), edvis_z(a), zetamult(a),dyf(a), w1, spec(a), homo(a), dum1(a), dy(a)
!    REAL(dp):: csi_flux1, csi_flux2, bedload_flux, shields(a)
!    INTEGER:: b(1), jjj, i
!    
!
!    hss2f=0._dp ! Predefine hss2f, which stores (lateral bedload flux) / (edvisy * lateral bedslope)
!    IF(susQbal) THEN 
!        !In this case, we are assuming a balance between lateral bedload and the lateral suspended flux
!        
!        ! Lateral eddy viscosity
!        IF(lambdacon>0._dp) THEN
!            edvisy= lambdacon*sqrt(abs(tau)/rho)*max(water-bed,0._dp)
!        ELSE
!            edvisy= 0.2_dp*sqrt(abs(tau)/rho)*max(water-bed,0._dp)
!        END IF
!
!        ! Compute downslope bedload transport
!        ! hss2f*dh/dy = (downslope bedload transport) / (lateral eddy viscosity)
!        IF(maxval(abs(Qbed))>0._dp) THEN
!            shields = (abs(tau)/(rho*g*(rhos/rho-1._dp)*d50))
!            DO i=1, a
!                IF((abs(tau(i))>0._dp).and.(edvisy(i)>0._dp)) THEN 
!
!                    hss2f(i)= -abs(Qbed(i))/edvisy(i)*1._dp/(sqrt(shields(i)))
!                    IF(talmon) hss2f(i)=hss2f(i)/9._dp*(max(water-.5_dp*(bed(i)+bed(i)),0._dp)/d50)**.3_dp 
!
!                ELSE 
!                    hss2f(i)=0._dp
!                END IF
!            END DO
!            ! Calculate coef at i+1/2
!            hss2f(1:a-1)=0.5_dp*(hss2f(2:a)+hss2f(1:a-1))
!        END IF
!
!        ! Compute hss2f at boundaries
!        IF(abs(tau(1))>0._dp) THEN
!
!            hss2f(0)= -.5_dp*(abs(Qbed(1)))/edvisy(1)*1._dp/sqrt(shields(1)) 
!            IF(talmon) hss2f(0)=hss2f(0)/9._dp*(max(water-.5_dp*(bed(1)+bed(1)),0._dp)/d50)**.3_dp
!
!        ELSE
!            hss2f(0)=0._dp
!        END IF
!        
!        ! Compute hss2f at boundaries
!        IF(abs(tau(a))>0._dp) THEN
!            hss2f(a) = -.5_dp*abs(Qbed(a))/edvisy(a)*1._dp/sqrt(shields(a))
!            IF(talmon) hss2f(a) = hss2f(a)/9._dp*(max(water-.5_dp*(bed(a)+bed(a)),0._dp)/d50)**.3_dp
!
!        ELSE
!            hss2f(a)=0._dp
!        END IF
!
!    END IF !susQbal
!
!    !Vertical eddy diffusivity
!    edvis_z= 0.1_dp*sqrt(abs(tau)/rho)*max(water-bed,0._dp) 
!    !print*, maxval(edvis_z)
!    !Zetamult*cbed= ZETA= depth integrated sediment concentration
!    DO i=1, a
!        IF((edvis_z(i)>0._dp)) THEN 
!            zetamult(i)= edvis_z(i)/wset*(1._dp-exp(-(wset/edvis_z(i))*max(water-bed(i),0._dp)) )
!        ELSE !Avoid division by zero in this case
!            zetamult(i)=0._dp
!        END IF
!    END DO
!    
!    !Forward difference increment
!    dyf(1:a-1)=ys(2:a)-ys(1:a-1) 
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! Calculate homogenous part of solution up to some arbitrary constant multiple
!    ! and the right hand side term
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!    DO i=1, a
!              
!        IF(i==1) THEN
!            w1=(zetamult(i)+(bed(i)-bedu))
!                IF(w1.NE.0._dp) THEN
!                    spec(i)= rhos*hss2f(i-1)*(bed(i)-bedu)/w1
!                ELSE
!                    spec(i)=0._dp
!                END IF 
!        
!            homo(i)=  1.0E-07_dp !An arbitrary constant to make things non-zero 
!        ELSE 
!            !   if((abs(2._dp*zetamult(i) + bed(i)-bed(i-1))>1.0E-01_dp).or.(bed(i)==bed(i-1))) THEN
!            IF(bed(i)==bed(i-1)) THEN 
!                w1=(zetamult(i)+0.5_dp*(bed(i)-bed(i-1)))
!                    IF(w1.NE.0._dp) THEN 
!                        spec(i)=(rhos*hss2f(i-1)*(bed(i)-bed(i-1)) +spec(i-1)*(zetamult(i-1) -0.5_dp*(bed(i)-bed(i-1))))/w1
!                        homo(i)=homo(i-1)*( zetamult(i-1) -0.5_dp*(bed(i)-bed(i-1)))/w1
!                    ELSE
!                        spec(i)=spec(i-1)
!                        homo(i)=homo(i-1)
!                    END IF
!
!            ELSE !FIXME: Try to reduce round-off errors -- seems to work -- but apparently could be improved
!                w1=(zetamult(i)/(bed(i)-bed(i-1))+0.5_dp)
!                    IF(w1.NE.0._dp) THEN 
!                        spec(i)=(rhos*hss2f(i-1) +spec(i-1)*(zetamult(i-1)/(bed(i)-bed(i-1)) -0.5_dp))/w1
!                        homo(i)=homo(i-1)*( zetamult(i-1)/(bed(i)-bed(i-1)) -0.5_dp)/w1
!                    ELSE
!                        spec(i)=spec(i-1)
!                        homo(i)=homo(i-1)
!                    END IF
!            END IF
!        END IF
!        !Set things to zero if there is no possibility to suspend sediment
!        !if(zetamult(i)==0._dp) THEN
!        !homo(i)=0._dp
!        !spec(i)=0._dp
!        !end if
!    END DO
!
!    ! Error check          
!    DO i=1,a
!        IF((isnan(homo(i))).OR.(isnan(spec(i)))) THEN 
!            print*, "C(",i,") is nan", edvis_z(i), zetamult(i), homo(i), spec(i)!, rhs(i), hsu 
!            open(632,file='DUMPER')
!            write(632,*) C, homo,spec,tau,bed,zetamult,edvis_z
!            close(632)
!            !print*, '##############CCCCCCCC##################'
!            !print*, C, homo, spec
!            !print*, '##############tausssss##################'
!            !print*, tau
!            !print*, '##############depths####################'
!            !print*, water-hs
!            STOP
!        END IF
!        IF(homo(i)<0._dp) homo(i)=abs(homo(i)) !Can get small negative values, or perhaps can we get the arbitrary multiplicative constant being negative?
!    END DO
!    
!    !Define the point with maximum depth. 
!    b=maxloc(water-bed) !ceiling(a*.5_dp)
!    jjj=b(1)
!    !Renormalise the homogeneous solution so that homo(jjj)=sconc
!    IF(homo(jjj).ne.0._dp) THEN
!        homo=homo/homo(jjj)*sconc
!    end if
!
!    ! Compute cross-sectionally integrated flux for homo 
!    dy(2:(a-1))= 0.5_dp*(ys(3:a)-ys(1:(a-2))) 
!    dy(1)=0.5_dp*(ys(2)-ys(1)) +(ys(1)-ysl)*(water-bed(1))/(bedl-bed(1)) 
!    dy(a)=0.5_dp*(ys(a)-ys(a-1)) +(ysu-ys(a))*(water-bed(a))/(bedu-bed(a))
!    csi_flux1 = sum(homo*zetamult*vel*dy) !Cross-sectionally integrated flux
!
!    !print*, 'susinfo', maxval(homo), minval(homo), maxval(spec), minval(spec), jjj,spec(1)
!    
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! Construct the final solution
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    IF(susQbal) THEN !The case when we have bedload is more complex than when we do not.
!        !Here we do some things to ensure that C is non-negative
!        dum1=9.9E+20_dp !Just a really big number
!        !Find the minimum of spec/homo 
!        DO i=1, a
!            IF(abs(homo(i))>0._dp) THEN 
!                dum1(1)= min(spec(i)/abs(homo(i)), dum1(1))
!            END IF
!        END DO
!
!        IF(dum1(1)==9.9E+20_dp) THEN
!            print*, 'Overflow in the suspended sed dist'
!            open(632,file='DUMPER')
!            write(632,*) C, homo,spec,tau,bed,zetamult,edvis_z
!            close(632)
!            stop
!        END IF
!        !print*, "ddd", dum1(1),maxval(homo), minval(homo)
!        C = spec - dum1(1)*abs(homo) !Force the min val of C to be zero
!        !print*, dum1
!        
!        csi_flux2 = sum(C*zetamult*vel*dy) !Cross-sectionally integrated flux for special solution
!        !Now, the boundary condition is sort of ambiguous. We could add more homo
!        !to the above and it would still be a solution. How about we assume that we add
!        !x*homo as well.
!        !C=1._dp*C+0.1_dp*homo
!
!        !Add multiples of the homogeneous solution so that the total flux (suspended load + bedload) is = integrated_load_flux
!        !C = special_solution + x*homogeneous_solution
!        !integrated_load_flux = bedload_flux+ csi_flux2 + x*csi_flux1
!        !Therefore, x = (integrated_load_flux - bedload_flux - csi_flux2)/csi_flux1
!        bedload_flux = sum(Qbed*dy)*rhos ! (m^3/(s*m)*m*)kg/m^3 = kg/s
!        IF(csi_flux1 > 0._dp) THEN
!            C = C+ ( (integrated_load_flux - bedload_flux -csi_flux2)/csi_flux1)*homo
!        END IF
!        !C = max(C,0._dp)
!        IF(integrated_load_flux<=csi_flux2+bedload_flux) THEN
!            !This situation is problematic, because the above approach would
!            !cause negative values of C in some locations.
!            print*, 'Special solution flux + bedload_flux = ', csi_flux2, '+', bedload_flux,'=', csi_flux2+bedload_flux, & 
!                    ' > integrated_load_flux = ', integrated_load_flux, maxval(spec-dum1(1)*abs(homo))
!            ! Perhaps temporarily let the solution pass through this state.
!            C = .50_dp*((spec - dum1(1)*abs(homo))/csi_flux2 + homo/csi_flux1)*(integrated_load_flux-bedload_flux)
!        !    stop
!        END IF
!    
!        !print*, C(jjj)
!    ELSE 
!        !Here, we calculate the homogenous solution, if there is no need to balance bedload and suspended load
!        !C=1._dp*homo
!        bedload_flux = sum(Qbed*dy)*rhos ! (m^3/(s*m)*m*)kg/m^3 = kg/s
!        if(csi_flux1>0._dp) C=homo/csi_flux1*(integrated_load_flux-bedload_flux)
!    END IF !susQbal
!
!    !!!!!!!!!!!!!!!!!!!!!!!
!    !Sanity checks
!    !!!!!!!!!!!!!!!!!!!!!!!
!    DO i = 1, a
!        IF(C(i)<0._dp) THEN
!            !if(C(i)<-(10._dp)**(-7)) THEN
!            print*, "C<0", C(i), i, water-bed(i)
!            !if(water-bed(i)>0.01_dp) print*, "C<0", C(i), i, water-bed(i)
!            !end if
!            C(i)=0._dp
!        END IF
!
!        IF(isnan(C(i))) THEN 
!            print*, "C(",i,") is nan", edvis_z(i), zetamult(i), dum1(1), homo(i), spec(i), csi_flux1, csi_flux2, bedload_flux!, rhs(i) 
!            
!            stop
!        END IF
!    END DO
!
!
!    !IF(maxval(C)>sconc + 10._dp**(-5)) THEN
!    !    !b=maxloc(C)
!    !    !print*, "maxval C > sconc + 10._dp**(-5), = ", maxval(C), b, a, water-bed((b(1)-1):(b(1)+1)), zetamult((b(1)-1):(b(1)+1)), &
!    !    !minval(water-bed(1:b(1))), minval(water-bed(b(1):a)), minval(water-bed(jjj:b(1)))
!
!    !    !Brute force solution to the occasional spikes in C
!    !    DO i=1, a
!    !        IF((water-bed(i)<0.01_dp).and.(C(i)>sconc)) THEN
!    !            print*, 'shallow C spike', i, C(i), water-bed(i)
!    !            C(i)=sconc !min(C(i), sconc)
!    !        ELSE
!    !            !if(C(i)>sconc*1.1_dp) print*, 'Deep C spike', C(i), water-bed(i), i
!    !        END IF
!    !    END DO
!    !END IF
!END SUBROUTINE sus_dist_sect
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dynamic_sus_dist(a, delT, ys, bed, water, waterlast, Q, tau, vel, wset, Qe,lambdacon, &
                                rho,rhos, g, d50, bedl,bedu, ysl, ysu, cb, Cbar, Qbed, &
                                sconc, counter, high_order_Cflux, a_ref, sus_vert_prof)
    ! Calculate the cross-sectional distribution of suspended sediment using
    ! some ideas from /home/gareth/Documents/H_drive_Gareth/Gareth_and_colab
    ! s/Thesis/Hydraulic_morpho_model/channel_cross_section/paper/idea_for_
    ! simplified_cross_sectional_sus
    
    ! We solve the equation
    ! depth d (Cbar) / dt + U*depth dCbar/dx +
    ! V*depth* d(Cbar)/dy - d/dy( eddif_y d(depth Cbar)/dy + eddif_y*cb*dh/dy) -
    ! (Es - ws*cb) = 0.
    ! With suitable approximations for dCbar/dx and V
    
    INTEGER, INTENT(IN)::a, counter
    REAL(dp), INTENT(IN):: delT, ys, bed, water, waterlast, tau, vel, wset, Qe, lambdacon, rho, rhos,g, & 
                                d50, bedl, bedu,ysl,ysu, sconc, Q, Qbed, a_ref
    REAL(dp), INTENT(IN OUT):: cb, Cbar ! Near bed suspended sediment concentration, Depth averaged suspended sediment concentration
    CHARACTER(20), INTENT(IN):: sus_vert_prof
    LOGICAL, INTENT(IN):: high_order_Cflux
    DIMENSION ys(a), bed(a), tau(a),vel(a), Qe(a), cb(a), Cbar(a),Qbed(a), a_ref(a) 

    ! LOCAL VARIABLES
    INTEGER:: i, info
    LOGICAL:: halt
    REAL(dp):: depth(0:a+1), eddif_y(0:a+1), eddif_z(a), zetamult(0:a+1), vd(0:a+1), ys_temp(0:a+1)
    REAL(dp):: M1_lower(a), M1_diag(a), M1_upper(a), M1_upper2(a), M1_lower2(a)
    REAL(dp):: RHS(a), dy_all(a)
    REAL(dp):: tmp1, dy, dy_outer, xlen, tmp2, Cref, z
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
        IF(counter.eq.1) print*, 'WARNING: Using default eddy diffusivity &
                                    of 0.2 in in dynamic_sus_dist, because &
                                    lambdacon = 0.0'
        eddif_y(1:a)= 0.2_dp*sqrt(abs(tau)/rho)*depth(1:a) 
    END IF
    ! Include zero depth boundaries 0 and a+1
    eddif_y(0) = 0._dp
    eddif_y(a+1) = 0._dp
   
    !IF(.FALSE.) THEN
    !    eddif_y(1:a)=maxval(eddif_y(1:a)) !eddif_y(1:a)+0.01_dp
    !    IF(counter.eq.1.) PRINT*, 'WARNING: constant eddy diff'
    !END IF 
    

    !IF(.FALSE.) THEN
    !    eddif_y=0._dp
    !    IF(counter.eq.1) print*, 'WARNING: Zero eddy diffusivity in dynamic_sus_dist'
    !END IF

    ! Calculate the value of 'zetamult', where:
    ! zetamult*cbed = depth integrated sediment concentration = Cbar*d  
    SELECT CASE(sus_vert_prof) 

        CASE('exp')
            ! Exponential suspended sediment distribution

            ! Vertical eddy diffusivity
            !eddif_z= 0.067_dp*sqrt(abs(tau)/rho)*depth(1:a) 
            eddif_z= 0.1_dp*sqrt(abs(tau)/rho)*depth(1:a) 
            !eddif_z= (lambdacon/2.0_dp)*sqrt(abs(tau)/rho)*depth(1:a) 
            
            !IF(.FALSE.) THEN
            !    eddif_z(1:a)=maxval(eddif_z(1:a)) !eddif_y(1:a)+0.01_dp
            !    IF(counter.eq.1.) PRINT*, 'WARNING: constant eddy diff'
            !END IF 
            
            DO i=1, a
                IF((eddif_z(i)>0._dp).and.(depth(i)>0.0e-00_dp)) THEN 
                    zetamult(i)= eddif_z(i)/wset*(1._dp-exp(-(wset/eddif_z(i))*max(depth(i),0._dp)) )
                ELSE 
                    zetamult(i)=1.0e-012_dp !1.0e-04_dp
                END IF
            END DO

        CASE('Rouse')
            ! Rouse vertical suspended sediment distribution
            DO i=1,a

                IF(abs(tau(i))>0.0_dp) THEN
                    ! Rouse number
                    z = wset/(0.4_dp*sqrt(tau(i)/rho)) 
                ELSE
                    zetamult(i) = 0.0_dp
                    CYCLE
                END IF

                ! Compute rouse integral factor using a function
                zetamult(i) = rouse_int(z,a_ref(i)/depth(i))

            END DO

        ! Catch errors
        CASE DEFAULT
            PRINT*,  'ERROR: sus_vert_prof does not have the correct value in dynamic_sus_dist'
            stop
        
    END SELECT

    ! Limit zetamult to a non-zero value to avoid division problems
    zetamult(0)   = 1.0e-012_dp
    zetamult(a+1) = 1.0e-012_dp
    zetamult = max(zetamult, 1.0e-012)


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
        
    !Check if the mesh is constant -- because if it is, we can use high order
    !derivative methods
    IF(maxval(dy_all) - minval(dy_all) < 1.0e-010_dp) THEN
        const_mesh=.TRUE.
    ELSE
        const_mesh=.FALSE.
    END IF 

    ! Advection term -- dC/dx is = (C - C/mean_C*desired_C)/x_length_scale
    dy_all = (ys_temp(2:a+1)-ys_temp(0:a-1))*0.5_dp
    tmp1 = sum(Cbar*vel*depth(1:a)*dy_all) !Cbar flux
    tmp2 =sconc*abs(Q) - sum(Qbed(1:a)*dy_all) !Desired Cbar flux = 'Measure of total load less bedload'
    tmp2 = max(tmp2, 0._dp)

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
            ! 2 point derivative approx
            M1_upper(i) = M1_upper(i) - tmp1*depth(i+1)
            M1_diag(i)  = M1_diag(i)  + tmp1*depth(i)
        END IF
 
        tmp1 = 0.5_dp*(eddif_y(i) + eddif_y(i-1))/((ys_temp(i) - ys_temp(i-1))*dy_outer)
        IF(i>1) THEN
            ! 2 point derivative approx
            M1_diag(i)  = M1_diag(i)  + tmp1*depth(i)
            M1_lower(i) = M1_lower(i) - tmp1*depth(i-1)
        END IF

        ! d/dy ( eddify* dbed/dy * cb)

        ! First compute eddify*dbed/dy at i+1/2
        tmp1 = 0.5_dp*(eddif_y(i+1)+eddif_y(i))
        tmp1 = tmp1*( -(depth(i+1)-depth(i) )/(ys_temp(i+1)-ys_temp(i)))
        tmp1 = tmp1/dy_outer

        IF(i<a) THEN
            !Estimate of 1/dy_outer*(eddify* dbed/dy *cb) at i+1/2
            M1_upper(i) = M1_upper(i) - 0.5_dp*tmp1*(depth(i+1)/zetamult(i+1))  ! Note that depth(i)/zetamult(i)*Cbar = cb
            M1_diag(i)  = M1_diag(i)  - 0.5_dp*tmp1*(depth(i)/zetamult(i))
        END IF

        ! Compute eddify*dbed/dy at i-1/2
        tmp1 = 0.5_dp*(eddif_y(i)+eddif_y(i-1))
        tmp1 = tmp1*( -(depth(i)-depth(i-1) )/(ys_temp(i)-ys_temp(i-1)))
        tmp1 = tmp1/dy_outer

        IF(i>1) THEN
            !Estimate of 1/dy_outer*(eddify*dbed/dy*cb) at i-1/2
            M1_diag(i)   = M1_diag(i)   + 0.5_dp*tmp1*(depth(i)/zetamult(i))  ! Note that depth(i)/zetamult(i)*Cbar = cb
            M1_lower(i)  = M1_lower(i)  + 0.5_dp*tmp1*(depth(i-1)/zetamult(i-1))
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
    IF(.TRUE.) THEN
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
    
    ! Solve matrix equations
    XXX(1:a,1) = RHS(1:a)
    call DGTSV(a, 1, bandmat(4,1:a-1), bandmat(3,1:a), bandmat(2,2:a), XXX(1:a,1), a, info)

    ! New Cbar, converted to kg/m^3
    Cbar = XXX(1:a,1)*rhos

    IF(info.ne.0) THEN
            print*, 'ERROR: info = ', info, ' in DGTSV, dynamic_sus_dist'
        stop
    END IF
    
    ! Sometimes, the near bed diffusive flux can lead to small negative values
    ! of Cbed being predicted right near the channel edge. This is because when
    ! the depth is changing by order of magnitude (e.g. from 1.0e-05 to 1.0e-02)
    ! outgoing flux from the small depth point [say (cbed*eddvis*slope)(i+1/2)
    ! ] can be much greater than the amount of suspended sediment over the small
    ! depth point. In my experience, this often happens at points with a tiny
    ! depth right next to the dry bank. (e.g. depth(i-1)=0.0)
    ! A simple fix is to make Cbed = 0.0 if negative, and remove the required
    ! amount of suspended sediment from the neighbouring point. 
    DO i = 1,a
        IF((Cbar(i)<0._dp).and.(.TRUE.)) THEN
            ! Clip negligably small Cbar values
            IF(abs(Cbar(i))<1.0e-10) THEN
                Cbar(i) = 0._dp
                !goto 21121            
                CYCLE
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
END SUBROUTINE dynamic_sus_dist
!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_friction(friction_type, grain_friction_type, rough_coef, water,&
                            a, bed, vel, man_nveg,d50,veg_ht, rhos, rho, g, f,&
                            vegdrag,f_g, dsand, counter, a_ref) 
    INTEGER, INTENT(IN):: a, counter
    CHARACTER(LEN=20), INTENT(IN):: friction_type, grain_friction_type
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


    
    !f(1) = -9999.0_dp !Trick to catch errors in setting the variable 'friction_type' 
   

    ! Friction factors
    SELECT CASE(friction_type)
    !IF(friction_type == 'manning') THEN
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
    !END IF
        CASE('darcy')
    !IF(friction_type == 'darcy') THEN
            !Darcy-weisbach style friction factor
            f= rough_coef 
    !END IF !m
        CASE('ks') 
    !IF(friction_type == 'ks') THEN
            ! Roughness height style friction factor
            ks = rough_coef

            DO i=1,a
            !f(i) = 8._dp*g/(18.0_dp*log10(12.0_dp*max((water-bed(i)),ks)/ks))**2
            !f(i) = 8._dp*g/(18.0_dp*log10(12.0_dp*max((water-bed(i)),ks)/ks))**2
            f(i) = 8.0_dp*(0.4_dp/log(max(water-bed(i), 3.0_dp*ks)/ks - 1.0_dp))**2
            END DO
    !END IF
        CASE('vanrijn') 
    !IF(friction_type == 'vanrijn') THEN
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
             ! have a type
             f_g =(8._dp*g/(18._dp*log10(12._dp*max(water-bed, 20.0_dp*3.0_dp*d50)/(3._dp*d50)+0.0_dp)+0.0e+00_dp)**2)

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
                        PRINT*,'f_g did not converge in ', j,' iterations' 
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
        !a_ref = max(0.01_dp, 0.5_dp*k_scr) !, 0.01_dp*(water-bed(i))) !Reference level (m) 
        a_ref(i) = min(max(0.5_dp*k_scr, 0.01_dp*(water-bed(i))),0.5_dp*(water-bed(i))) !Reference level (m) 
    END DO
    

END SUBROUTINE calc_friction
!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(dp) FUNCTION rouse_int(z,d_aref)
    ! Suppose that Cbar = cbed*K
    ! Where Cbar is the depth-averaged sediment concentration
    ! cbed is the near bed concentration
    ! and K is an integrating factor, which is determined from the Rouse
    ! distribution for suspended sediment.
    ! Then this function calculates K.
    
    ! INPUTS
    ! z = rouse number = wset/(0.4*ustar)
    ! d_aref = dimensionless reference level 
    !        = (van rijn reference level)/ depth

    ! NOTE -- This uses a method from Guo and Julien (2004) 'Efficient
    ! Algorithms for Computing Einstein Integrals', Journal of Hydraulic
    ! Engineering 130:1198-1201.
    ! They show how to calculate J1=int_{d_aref}^{1} [(1-eps)/eps]^z d(eps)
    ! It can be shown that K=J1*db_const (where db_const is defined below) 
    REAL(dp), INTENT(IN):: z, d_aref

    INTEGER:: i
    REAL(dp):: db_const, F1, J1, j, E2, z2, rouse_int
   
    ! If z>10.0, there is no suspended load, make a quick exit, otherwise proceed with algorithm 
    IF((z>10.0_dp).or.(d_aref>0.3_dp)) THEN
        !FIXME: Check that the 0.3 above is okay
        rouse_int = 0.0_dp 
        
    ELSE

        ! Prevent integer values of z 
        IF(abs(z-anint(z))<0.0005_dp) THEN
            z2 = anint(z)+0.0005_dp
        ELSE
            z2=z
        END IF

        ! Compute 1/((1-deltab)/deltab)^z2
        db_const = ((1.0_dp-d_aref)/d_aref)**(-z2)

        ! Compute F1, eqn 8
        F1 =((1.0_dp-d_aref)**z2)/d_aref**(z2-1.0_dp) 
        E2 = d_aref/(1.0_dp-d_aref)

        DO i=1,10
            j=i*1.0_dp
            F1 = F1 - z2*((-1)**j)/(j-z2)*(E2)**(j-z2)
        END DO

        ! Compute J1, eqn 9
        J1 = z2*pi/sin(z2*pi) -F1

        ! Compute the desired integration factor
        rouse_int=J1*db_const

        IF((rouse_int<0.0_dp).or.(rouse_int>=1.0_dp)) THEN
            PRINT*, ' ERROR in rouse_int: unphysical rouse_int value ', rouse_int, d_aref, z
            stop
        END IF
    END IF

END FUNCTION rouse_int
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE crosssection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

