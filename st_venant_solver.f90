!!!!!!!!!!!!!!!!!! Time-stepping routine for variant of the St Venant Equations
Module st_venant_solver

use global_defs
IMPLICIT NONE

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE hyupdate(delT, delX, bottom, B, A, Q,Q2H_dT,Y, t,l, counter,dWidth_dwaters, rmu,inuc, LAKE, hlim, Qb, tr,&
 mouth_data, mouthread,mouth_data_len, mns,rho, g,cfl, v1coef,v4coef,seabuf)
    ! delT = time-step (computed herein based on the CFL condition)
    ! delX = distance between cross-sections (a constant, though this should be
    !        reasonably easy to hack if variable delX is desired)
    ! bottom    = the mean wetted-bed elevation
    ! B    = wetted Width
    ! A    = Cross-sectional area
    ! Q    = Discharge (m^3/s)
    ! Q2H_dT  = The time-integrated discharge (m^3), computed at spatial points
    !        i+1/2. This satisfies A_new(i) -A_old(i) = (Q2H_dT(i)-Q2H_dT(i-1) )/ delX
    !        It allows us to keep track of the discharge over multiple delT
    !        steps (which will sum to dT).
    ! Y    = Free surface elevation
    ! t    = time
    ! l    = number of cross-sections
    ! counter = a counter (keeps track of how many times this has been called).
    ! dWidth_dwaters = rate of change of width w.r.t. the water elevation. This can hold 2
    !        values for each cross-section (so allowing a different treatment
    !        for rising and falling water levels), but this is not used at
    !        present.
    ! rmu  = (roughness_coefficient/depth). The roughness coefficient and depth
    !        are computed from the cross-sectional hydrodynamics / morphology.
    !        For darcy friction, rmu is basically f/(8*g*depth). 
    !        'depth' appears in rmu because doing so adds stability to the
    !        algorithm
    ! inuc = correction constant to the convective intertial terms, computed
    !        from the cross-sectional models. Accounts for the non-uniformity of
    !        velocity over each cross-section. The convective terms in this code
    !        are equal to (1+inuc)*d(Q^2/A)/dx.
    ! LAKE = True/False -- if True, assume a constant water surface elevation
    !        throughout the channel, and use the continuity equation to solve for
    !        velocities. Should be False to solve the St Venant equations
    ! hlim = small depth, at which we start doing funny things to avoid negative
    !        depths -- a wetting and drying trick
    ! Qb   = River discharge imposed at boundary 'l'
    ! tr   = tidal range, which can be passed to 'mouth_height' if
    !        mouthread=FALSE
    ! mouth_data = a 2 column matrix containing (time, waterlevel) data that can be
    !          used to set the downstream boundary condition
    ! mthread = True/false -- determine if we use mouth_data to define the downstream boundary
    !           condition.
    ! mouth_data_len = number of rows in mouth_data
    ! mns = DEFUNCT PARAMETER
    ! rho = water density
    ! g = gravity
    ! cfl = cfl number (<= 1) which influences the size of the time step
    ! v1coef = coefficient for artificial viscosity
    ! v4coef = DEFUNCT coefficient for artificial viscosity
    ! seabuf = number of 'morphologically frozen' cross-sections at the
    !          downstream boundary. In the driver routine, only cross-sections
    !          'seabuf+1:l' are alowed to evolve. At some stage, I found it
    !          useful to  only apply artificial viscosity to the evolving
    !          cross-sections, I think because in one application, the
    !          cross-sectional area changed rapidly between section seabuf and
    !          seabuf+1.

    INTEGER, INTENT (IN)::l, mouth_data_len,seabuf
    REAL(dp), INTENT(IN):: delX, Qb, tr, mouth_data, hlim, rho, mns, g, cfl, bottom, B
    REAL(dp),INTENT(IN OUT)::t, delT
    REAL(dp), INTENT(IN OUT):: A,Q, Q2H_dT, Y
    DIMENSION:: bottom(l), B(l),A(l), Q(l), Q2H_dT(0:l), Y(l),mouth_data(mouth_data_len,2), mns(l) 
    INTEGER, INTENT (IN):: counter
    REAL(dp), INTENT(IN)::dWidth_dwaters(l,2), rmu(l),inuc(l)
    REAL(dp), INTENT(IN):: v1coef, v4coef
    LOGICAL, INTENT(IN)::LAKE, mouthread

    !Local variables-- this is a mess, you really need to look at the code to see what they mean, I just wacked them in as I needed them
    !logical:: isnan
    REAL(dp)::d0
    REAL(dp):: th, dlim
    INTEGER:: i,  j, swit,lo,ind  
    REAL(dp)::Ashift(l),Qshift(l),Bbshift(l)
    REAL(dp)::Alast(l),Qlast(l), df(l), vf(l), dc(l), vc(l),   zer(l), zer2(l+1)
    REAL(dp)::Apred(l),Qpred(l),Qpred1(l),Qpred2(l),  efd(l),efd2(l)
    REAL(dp)::Apredbshift(l), Qpredbshift(l), Qpredc(l), bottombshift(l)
    REAL(dp)::Acor(l),Qcor(l),Qcor1(l),Qcor2(l),Qcor3(l), Qcor4(l), Af(l),  Acorb(l)
    REAL(dp):: Areafil(l), Qfil(l),bottomshift(l) 
    REAL(dp)::Ylast(l),Ypred(l),Ycor(l), Yshift(l),  YY(l), QQ(l), slop(l), Q2_tmp(0:l) 
    REAL(dp)::Bnew(l), inucshift(l), s(l), up(l), low(l), diag(l),  visc(l),   Qext, w1, w3, useme, useme2(l)
    REAL(dp):: FLUXp(l,1:2), FLUXm(l,1:2), FLUXpP(l,1:2), FLUXmP(l,1:2), viscF(l), viscf4(l)
    LOGICAL::chnge   

    dlim=1._dp  !Whenever I use dlim in the code, it appears as a proportion of hlim -- so 1 means the discharge limit will be hlim

    Alast=A
    Qlast=Q
    Ylast=Y

    IF (LAKE) THEN !!The momentum equation reduces to a flat free surface.
        delT=10._dp 
        GOTO 1209 !Straight to the water elevation and the lake hydrodynamics
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate the time step 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    delT=1000._dp

    call time_step(delT, cfl, delX, Y-bottom, abs(Q/A), rmu, l )
    1209 t=t+delT !time in seconds

    ! Sanity Checks
    IF(abs(delT-1000._dp)<1.1_dp) print*, "delT= 1000."

    IF(delT<1.0E-04_dp) THEN 
        print*, "short time step", delT
        print*, "counter=", counter, "maxloc=", maxloc(abs(Q/A)), "l=", l, "maxvel=", maxval(abs(Q/A))
        print*,  Q/A
        !print*, maxval(bottom), minval(bottom)!maxval(efd), maxval(abs(Q)/A)
        STOP
    END IF


    !Sanity check
    IF(isnan(t)) THEN
        call Dat
        PRINT*, " Error: time is NAN", Y, bottom!, efd
        STOP
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Compute Mouth boundary condition 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call mouth_height(th, t,tr, mouth_data,mouthread, mouth_data_len)

    d0=th 

    IF(LAKE) THEN  !This is supposed to calculate the flow for the theoretical case
        !of a flat free surface
        Y=max(d0+0._dp*bottom, bottom+dlim*hlim)
        swit= l
        DO i= 2, l
            IF(d0<bottom(i)) THEN
                swit = i
                GOTO 22
            END IF
        END DO                        
        22 Apred= (B*(Y-bottom))
        
        !!Solve delQ/delX = -delA/delT with matrix inversion, upwind method. 
        DO i = 2, swit-1
            IF(Q(i)>0._dp) THEN 
                diag(i)= 1._dp
                up(i)= 0._dp
                low(i)= -1._dp
                s(i)=  -(delX/delT)*(Apred(i)-A(i))
            ELSE
                diag(i)= -1._dp
                up(i)=1._dp
                low(i)=0._dp
                s(i)=  -(delX/delT)*(Apred(i)-A(i))
            END IF
        END DO

        low(1)=0._dp
        up(1)=0._dp
        diag(1)=1._dp
        s(1)= sum(Apred-A)*delX/delT

        call dgtsv(swit-1,1, low(1:(swit-1)), diag(1:(swit-1)), up(1:(swit-1)), s(1:(swit-1)),swit-1, j)

        Q(1:swit-1)=s(1:swit-1)
        Q(swit:l)=0._dp
        A=Apred
        
        GOTO 3541 !!go to the End of the program
    END IF

    !various shifted parameters that are helpful. 
    Ashift(1:l-1)= Alast(2:l)
    Qshift(1:l-1)= Qlast(2:l)
    Bbshift(2:l)= B(1:l-1)
    Yshift(1:l-1)=Ylast(2:l)

    Qshift(l)=Qb!Q(l)!Qb!2*Qshift(l-1)  -Qshift(l-2)
    Ashift(l)=A(l)!2*Ashift(l-1)  -Ashift(l-2)
    Yshift(l)= 2.0_dp*Y(l)-Y(l-1) !2*Yshift(l-1)-Yshift(l-2)
    Bbshift(1)=B(1)!2*Bbshift(2)-Bbshift(3)
    bottomshift(1:l-1)=bottom(2:l)
    bottomshift(l)= 2._dp*bottom(l)-bottom(l-1)
    !bottombshift(2:l)=bottom(2:l)
    !bottombshift(1)=2._dp*bottom(1)-bottom(2) 
    inucshift(1:l-1)=inuc(2:l)
    inucshift(l)=inuc(l)

    !Now some halfway averages of the above channel variables
    Af= (Ashift+A)*0.5_dp  !/((Yshift+Y-hshift-h)*.5)  !

    !!!!!


    !!!!!!!!!Predictor water level
    Apred(1:l)= Alast(1:l)-(delT/delX)*(Qshift(1:l)-Q(1:l))

    !For this step, note that
    !Ypred=Ylast + (Apred-Alast)/(B+ 0.5_dp*db/dh*(Ypred-Ylast))
    !Rearranging things,
    !(Ypred-Ylast)*B + 0.5*db/dh*(Ypred-Ylast)^2= (Apred-Alast)
    !Ypred-Ylast= -B+-sqrt( B^2 - 2*db/dh*(Alast-Apred))/(2*0.5*db/dh)
    DO i=1,l
        IF(Apred(i)-Alast(i)>0._dp) THEN
            ind=1
        ELSE
            ind=2
        END IF

        w1= B(i)**2+2._dp*dWidth_dwaters(i,ind)*(Apred(i)-Alast(i))
        
        IF((dWidth_dwaters(i,ind)>0._dp).AND.(w1>0._dp)) THEN
            Ypred(i)=Ylast(i)+ (-B(i)+ sqrt(w1))/(dWidth_dwaters(i,ind))
        ELSE
            Ypred(i)= Ylast(i)+ (Apred(i)-Alast(i))/B(i)
        END IF

    END DO

    DO i = 1, l
        IF(isnan(Ypred(i))) THEN 
            PRINT*, "elev", i," is Nan", Ylast !, Qshift(i), B(i), i, Q(i+1), Q(i),l !, Bshift(i)
            PRINT*, '....Qlast....'
            PRINT*, Q
            PRINT*, '....dWidth_dwaters....'
            PRINT*, dWidth_dwaters
            PRINT*, '....Alast....'
            PRINT*, Alast
            PRINT*, '....B....'
            PRINT*, B
            STOP
        END IF
    END DO

    zer=1._dp

    !A more effective zeroing?
    zer2=1._dp
    DO i=1,l
        IF(Y(i)-bottom(i)<1._dp*hlim) THEN
            zer2(i)=0._dp
        END IF
    END DO
    zer2(l+1)=zer2(l)

    !Define the slope
    DO i=1,l
      slop(i)= (Yshift(i)-Y(i))  
    END DO

    DO i=2,l
        IF(Y(i)-bottom(i)<1._dp*hlim) THEN
            !slop(i-1:i)=0._dp
            slop(i)=0._dp
        END IF
    END DO

    
    !!!MOMENTUM EQUTION WITH IMPLICIT FRICTION
    Qpred1=Qlast-(delT/delX)*((1._dp+inucshift)*Qshift**2._dp/Ashift*zer2(2:l+1)-& 
        (1._dp+inuc)*Qlast**2._dp/Alast*zer2(1:l))*zer-&
        (delT/delX)*g*((Af*(slop))) ! !+A*(bottomsbottomift-bottom)*(1-zer))
    Qpred2= (g*Af*(-sign(1._dp, Qpred1)/(Apred**2._dp))*rmu )!*efd2 !*zer 
    !Explicit friction
    !Qpred2=(g*Af*(-sign(1._dp, Qlast)*Qlast**2/(Alast**2._dp))*rmu )
    !Normal method
    !Qpred= Qpred1+ delT*Qpred2

    visc=0._dp 

    !!Implicit
    !Qpred= Qpred1+visc+ delT*Qpred2*Qpred^2
    DO i=1,l
        IF(Qpred2(i).ne.0._dp) THEN
            Qpred(i)= (1._dp - sqrt(1._dp- 4._dp*delT*Qpred2(i)*(Qpred1(i)+visc(i)) ))/(2._dp*delT*Qpred2(i))
        ELSE
            Qpred(i)=Qpred1(i)
        END IF
    END DO

    !Boundary conditions --Actually if we extrapolate there is no need
    Qpred(l)= Qb 

    !Prevent negative areas later
    DO i=2,l
        useme=max((Alast(i)-dlim*hlim*B(i)),0._dp)*delX/delT
        IF((Qpred(i)-Qpred(i-1))> useme ) THEN
            Qpred(i)=0._dp !Qpred(i)/(Alast(i))*delT/delX
            Qpred(i-1)=0._dp !Qpred(i-1)/(Alast(i))*delT/delX
        END IF
    END DO
    

    !!!!!!!!!!!!!!!!!!
    ! Corrector step
    !!!!!!!!!!!!!!!!!!

    DO i=1,(l-1)
        Apredbshift(i+1)= Apred(i)
        Qpredbshift(i+1)=Qpred(i)
    END DO

    Apredbshift(1)= 2._dp*Apredbshift(2)-Apredbshift(3)  ! Apredbshift(2) !2*Apredbshift(2)-Apredbshift(3)
    Qpredbshift(1)= Apredbshift(1)*(2._dp*Qpredbshift(2)/Apredbshift(2) - Qpredbshift(3)/Apredbshift(3))!2*Qpredbshift(2) -Qpredbshift(3) !Qpredbshift(2) !2*Qpredbshift(2) -Qpredbshift(3)

    Qpredc=(Qpred+Qpredbshift)*0.5_dp
    Acorb=(Apred+Apredbshift)*0.5_dp !/((Ypredbshift+Y-hbshift-h)*.5)

    Bnew(1:l)= B(1:l) 
    !do i= 1, l
    !   if(Ypred(i)>Ylast(i)) THEN
    !   Bnew(i)=Bnew(i)+dWidth_dwaters(i,1)*(Ypred(i)-Ylast(i))
    !   ELSE
    !   Bnew(i)=Bnew(i)+dWidth_dwaters(i,2)*(Ypred(i)-Ylast(i))
    !   END IF
    !   IF(Bnew(i)<0._dp) THEN
    !           Bnew(i)=B(i)
    !           !dWidth_dwaters(i)=0._dp
    !   END IF
    !end do

    Acor(1:l)= Alast(1:l) - delT/delX*(Qpred(1:l)-Qpredbshift(1:l)) !!Note that we account for changes in B 
    DO i=1,l
        IF(Acor(i)-Alast(i)>0._dp) THEN
            ind=1
        ELSE
            ind=2
        END IF
        !To understand this step, see a similar treatment above
        w1= Bnew(i)**2+2._dp*dWidth_dwaters(i,ind)*(Acor(i)-Alast(i))
        
        IF((dWidth_dwaters(i,ind)>0._dp).and.(.true.).and.(w1>0._dp)) THEN
            Ycor(i)=Ylast(i)+ (-Bnew(i)+ sqrt(w1))/(dWidth_dwaters(i,ind))
        ELSE
            Ycor(i)= Ylast(i)+ (Acor(i)-Alast(i))/B(i)
        END IF
    END DO

    !!Slope
    DO i=2, l
        slop(i)= Ypred(i)-Ypred(i-1) 
    END DO
    slop(1)=slop(2)
    
    DO i=1,l-1
        IF(Ypred(i)-bottom(i)<1._dp*hlim) THEN
            slop(i)=0._dp
        END IF
    END DO

    !A more effective zeroing?
    zer2=1._dp
    DO i=1,l
        IF(Y(i)-bottom(i)<1._dp*hlim) THEN
            zer2(i+1)=0._dp
        END IF
    END DO
    zer2(1)=zer2(2)

    !Calculate discharge corrector

    DO i=2,l
        !!!!!!! IMPLICIT FRICTION FORM
        Qcor1(i)= Qlast(i)- (delT/delX)*& 
        ((1._dp+inuc(i))*Qpred(i)**2/Apred(i)*zer2(i+1)-(1._dp+inuc(i-1))*Qpred(i-1)**2/Apred(i-1)*zer2(i) )*zer(i)
        Qcor2(i)=(delT/delX)*g*(Acorb(i)*(slop(i)))! !+Apred(i)*(bottom(i)-bottom(i-1))*(1-zer(i)))
        Qcor4(i)=g*Acorb(i)*(-sign(1._dp,Qcor1(i)-Qcor2(i))/(Acor(i)**2._dp)*rmu(i))
        !Qcor4(i)=g*Acorb(i)*(-sign(1._dp,Qpred(i))*Qpred(i)**2/(Apred(i)**2._dp)*rmu(i))
    END DO
    !Boundary conditions
    Qcor1(1)= Qlast(1) 
    Qcor2(1)= delT/delX*g*Acorb(1)*slop(1)
    !Qcor4(1)= g*Acorb(1)*(-sign(1._dp, Qcor1(1)-Qcor2(1)))/(Acor(1)**2*efd(1))*rmu(1)
    Qcor4(1)=g*Acorb(1)*(-sign(1._dp,Qcor1(1)-Qcor2(1))/(Acor(1)**2._dp)*rmu(1))
    
    !!Normal form
    !Qcor=Qcor1-Qcor2+delT*Qcor4 !-(delT/delX)*g*0.5*Qcor3(2:l)

    visc=0._dp 

    !Implicit friction form -- use quadratic formula
    !Qcor= Qcor1 - Qcor2 +visc + Qcor^2*Qcor4*delT
    DO i=1, l
        IF(Qcor4(i).NE.0._dp) THEN
            Qcor(i)= (1._dp-sqrt(1._dp-4._dp*Qcor4(i)*delT*(Qcor1(i)-Qcor2(i)+visc(i))))/(2._dp*Qcor4(i)*delT)
        ELSE
            Qcor(i)= Qcor1(i)-Qcor2(i)
        END IF
    END DO


    !Limit Qcor to prevent negative depths
    DO i=1,l-1
        useme=max(0.5_dp*(Apred(i)+Acor(i))-dlim*hlim*B(i),0._dp)*delX/delT !A useful variable
        IF(0.5_dp*(Qcor(i+1)-Qcor(i))> useme-0.5_dp*(Qpred(i+1)-Qpred(i)) ) THEN  !We need to do limiting
            Qcor(i+1)=0._dp !-Qpred(i+1)
            Qcor(i)=0._dp !-Qpred(i)
        END IF
    END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Final estimates, prior to addition of viscosity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Q(1:l)=0.5_dp*(Qpred(1:l )+Qcor(1:l))
    A=0.5_dp*(Apred+Acor)

    !!Again, back calculate Y_{i}
    DO i=1,l
        IF(A(i)-Alast(i)>0._dp) THEN
            ind=1
        ELSE
            ind=2
        END IF
            w1= B(i)**2+2._dp*dWidth_dwaters(i,ind)*(A(i)-Alast(i))
        IF((dWidth_dwaters(i,ind)>0._dp).AND.(w1>0._dp)) THEN
            Y(i)=Ylast(i)+ (-B(i)+ sqrt(w1))/(dWidth_dwaters(i,ind))
        ELSE
            Y(i)= Ylast(i)+ (A(i)-Alast(i))/B(i)
        END IF
    END DO
    
    ! Artificial viscosity
    visc(2:l-1)=min(abs(Y(3:l)-2._dp*Y(2:l-1)+Y(1:l-2)),10000.1_dp)*&
    (abs(Qlast(2:l-1)/Alast(2:l-1))+sqrt(g*(Ylast(2:l-1)-bottom(2:l-1))))*delT/delX ! & 
    
    visc(2:l-1)=visc(2:l-1)/(abs(Y(3:l)-bottom(3:l)) &
    +2._dp*abs(Y(2:l-1)-bottom(2:l-1) ) + abs(Y(1:l-2)-bottom(1:l-2)) ) ! &
    
    visc(1)=(abs(Y(2)-Y(1)))/(abs(Y(2)-bottom(2))+abs(Y(1)-bottom(1)))*&
            (abs(Qlast(1)/Alast(1))+sqrt(g*(Ylast(1)-bottom(1))))*delT/delX
    visc(l)=(abs(Y(l)-Y(l-1)))/(abs(Y(l)-bottom(l))+abs(Y(l-1)-bottom(l-1)))*&
            (abs(Qlast(l)/Alast(l))+sqrt(g*(Ylast(l)-bottom(l))))*delT/delX

    !Here we modify the artificial viscosity to reduce the chance of drying.
    viscf=0._dp
    DO i= seabuf+3,l-1
        viscf(i)= min(max(visc(i),visc(i+1))*v1coef, 1000.5_dp)
        
        IF(i.EQ.seabuf+1) viscf(i)=min(visc(i+1)*v1coef,1000.5_dp)

        !viscf4(i)=max(0._dp,v4coef*1.0_dp-viscf(i))
        IF((minval(Y(max(i-1,1):i+1)-bottom(max(i-1,1):i+1))<5.00001_dp*hlim+0.00_dp)) THEN
            viscf(i)=0._dp
            viscf(max(i-1,1))=0._dp
        END IF
    END DO
    
    ! 'Conservative discharge'
    Q2_tmp(1:l)=0.5_dp*(Qshift +Qpred)
    Q2_tmp(0)=((delX/delT)*(A(1)-Alast(1)) + 0.5_dp*(Qshift(1)+Qpred(1)))

    !Apply viscosity to A and Q 
    Qpred2=A(1:l)
    Qpred1=Q
    DO i=2,l-1
        A(i)=A(i)+viscf(i)*(Qpred2(i+1)-Qpred2(i)) -viscf(i-1)*(Qpred2(i)-Qpred2(i-1))
        Q(i)=Q(i)+viscf(i)*(Qpred1(i+1)-Qpred1(i)) -viscf(i-1)*(Qpred1(i)-Qpred1(i-1)) 
    END DO

    A(1)=A(1)+viscf(1)*(Qpred2(1+1)-Qpred2(1)) 
    Q(1)=Q(1)+viscf(1)*(Qpred1(1+1)-Qpred1(1))
    
    A(l)=A(l)-viscf(l-1)*(Qpred2(l)-Qpred2(l-1)) 
    Q(l)=Q(l)-viscf(l-1)*(Qpred1(l)-Qpred1(l-1))

    !Check if we have actually enhanced local variability
    chnge=.false.
    IF(.false.) THEN
        DO j=1,300
            DO i=2,l-1
                IF((viscf(i)>0).and.(viscf(i-1)>0)) THEN
                    IF( (abs(Q(i+1)- 2._dp*Q(i) +Q(i-1) ) > &
                        1.0_dp*abs(Qpred1(i+1)- 2._dp*Qpred1(i)+Qpred1(i-1)) ).OR.& 
                        !(sign(1._dp, Q(i+1)- 2._dp*Q(i) +Q(i-1) ).ne. &
                        !sign(1._dp, Qpred1(i+1)- 2._dp*Qpred1(i)+Qpred1(i-1)) ).or.&
                        ( (abs(A(i+1)- 2._dp*A(i) +A(i-1) ) > &
                        1.0_dp*abs(Qpred2(i+1)- 2._dp*Qpred2(i)+Qpred2(i-1)) )) ) THEN !.or.&
                        !( sign(1._dp, A(i+1)- 2._dp*A(i) +A(i-1) ).ne. &
                        !sign(1._dp, Qpred2(i+1)- 2._dp*Qpred2(i)+Qpred2(i-1)) ))  THEN
                        viscf(i)=0._dp !viscf(i)*0.5_dp
                        viscf(i-1)=0._dp !viscf(i-1)*0.5_dp
                        chnge=.true. !If this is not true, we will break out of the loop
                    END IF
                END IF
            END DO

            DO i=2,l-1
                A(i)=Qpred2(i)+viscf(i)*(Qpred2(i+1)-Qpred2(i)) -viscf(i-1)*(Qpred2(i)-Qpred2(i-1))
                !
                Q(i)=Qpred1(i)+viscf(i)*(Qpred1(i+1)-Qpred1(i)) -viscf(i-1)*(Qpred1(i)-Qpred1(i-1)) 
            END DO
            
            A(1)=A(1)+viscf(1)*(Qpred2(1+1)-Qpred2(1)) 
            Q(1)=Q(1)+viscf(1)*(Qpred1(1+1)-Qpred1(1))
            !
            A(l)=A(l)-viscf(l-1)*(Qpred2(l)-Qpred2(l-1)) 
            Q(l)=Q(l)-viscf(l-1)*(Qpred1(l)-Qpred1(l-1))
            if(chnge.eqv..false.) goto 1449 !Break out of loop
            chnge=.false. !Reset chnge
        END DO
    END IF
    1449 CONTINUE
    
    !Correct for viscosity - this is supposed to ensure that the mass conservation still holds well. Note that Q2H_dT is actually the discharge integrated over the time step, not the raw discharge - so it has units m^3
    do i=1,l-1
    Q2_tmp(i)=Q2_tmp(i) -delX/delT*viscf(i)*(Qpred2(i+1)-Qpred2(i)) 
    end do
    !The above formula is designed to ensure that even with viscosity, A_{i}^{j+1}= A_{i}^{j} - delt/delx*[ Q_{i+1/2}^{j+1/2}-Q_{i-1/2}^{j+1/2} ]
    !As you see will see below, viscosity is implemented as A_{i}^{j+1} <-- A_{i}^{j+1} + visc(i+1/2)*(A_{i+1}^{j+1}-A_{i}^{j+1}) + visc(i-1/2)*(A_{i}^{j+1}-A_{i-1}^{j+1})
    !The correction in the above loop ensures that even for this updated A, the mass conservation with Q2H_dT is satisfied.

    !Make the 'mass conservative Q2H_dT', a discharge integrated over time (and divided
    !by the longer time step in driver2)
    !Numerically, this can be thought of as Q_{i+1/2}^{j+1/2}*delT where i is the spatial
    !increment and j is the temporal one. 
    Q2H_dT(1:l) = Q2H_dT(1:l)+ Q2_tmp(1:l)*delT
    Q2H_dT(0)= Q2H_dT(0)+ Q2_tmp(0)*delT

    !!Again, back calculate Y_{i}
    DO i=1,l
        IF(A(i)-Alast(i)>0._dp) THEN
            ind=1
        ELSE
            ind=2
        END IF
        w1= B(i)**2+2._dp*dWidth_dwaters(i,ind)*(A(i)-Alast(i))
        IF((dWidth_dwaters(i,ind)>0._dp).and.(w1>0._dp)) THEN
            Y(i)=Ylast(i)+ (-B(i)+ sqrt(w1))/(dWidth_dwaters(i,ind))
        ELSE
            Y(i)= Ylast(i)+ (A(i)-Alast(i))/B(i)
        END IF
    END DO
    
    !Condition suggested by Blayo and Debreu (2005). Qext=0 is my decision, since it
    !seems to do good things. - Could be interpreted as uext=0, indicating that the
    !cross sectional area at the boundary is so large that u is effectively 0.
    Qext=Qlast(1)*A(1)/Alast(1)! 0._dp!Q(1)!Q(1)/2._dp!0._dp!A(1)*(2.*Q(1)/A(1)-Q(2)/A(2))!0._dp!Qlast(1)!0._dp!.1*sum(A(2:l)+A(1:l-1) -Alast(2:l)-Alast(1:l-1))*delX/delT
    w1= Qext/A(1) +sqrt(g/(Y(1)-bottom(1)))*(d0-bottom(1))
    w3= Q(1)/A(1)-sqrt(g/(Y(1)-bottom(1)))*(Y(1)-bottom(1)) !2.*(Q(2)/A(2) -sqrt(g/(Y(1)-bottom(1)))*(Y(2)) ) - (Q(3)/A(3) -sqrt(g/(Y(1)-bottom(1)))*(Y(3))) !Note that this could be extrap in different ways.Q(1)/A(1)-sqrt(g/(Y(1)-bottom(1)))*Y(1)
    Q(1)= .5_dp*A(1)*(w1+w3) !.5*A(1)*(Qext/A(1) + 2.*Q(2)/A(2) -1.*Q(3)/A(3)+sqrt(g/(Y(1)-bottom(1)))*(d0 - 2.*Y(2)+1.*Y(3))) 
    Y(1)= sqrt((Y(1)-bottom(1))/g)*(-w3+w1)/2._dp +bottom(1)!.5*(d0+ 2.*Y(2) -1.*Y(3) - sqrt((Y(1)-bottom(1))/g)*(-Qext/A(1) + 2.*Q(2)/A(2) - 1.*Q(3)/A(3))) 

    !Y(1)=d0
  
    3541 CONTINUE

    DO i=1, l
        IF(isnan(Ypred(i))) print*, "Ypred ",i," is nan"
        IF(isnan(Ycor(i))) print*, "Ycor ",i," is nan", dWidth_dwaters(i,1:2), B(i), Bnew(i), Acor(i)-Alast(i)
        IF(isnan(Qpred(i))) print*, "Qpred ",i," is nan"
        IF(isnan(Qcor(i))) print*, "Qcor ",i," is nan", Qcor4(i), Qcor1(i), Qcor2(i), Acorb(i) & 
        , (1._dp-sqrt(1._dp-4._dp*Qcor4(i)*delT*(Qcor1(i)-Qcor2(i)+visc(i))))/(2._dp*Qcor4(i)*delT), &
        sqrt(1._dp-4._dp*Qcor4(i)*delT*(Qcor1(i)-Qcor2(i)+visc(i))), 4._dp*Qcor4(i)*delT*(Qcor1(i)-Qcor2(i)+visc(i))
    END DO


END SUBROUTINE hyupdate
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE time_step(delT, cfl, delX, depth, abvel, rmu, a) 
    INTEGER, INTENT(IN):: a
    REAL(dp), INTENT(IN):: cfl, delX
    REAL(dp), INTENT(IN):: depth(a), abvel(a), rmu(a)
    REAL(dp), INTENT(INOUT):: delT

    REAL(dp) mxvl, mxd, mnd

    !Useful things
    mxd= maxval(depth)
    !mnd= minval(depth)
    mxvl= maxval(abvel)

    !Time step
    !delT= min(delT,cfl*((delX/( ( 9.8_dp*( mxd))**0.5_dp+mxvl)))) !, & 
    !minval(.5_dp*mnd/((.001_dp+ abvel)*(9.8_dp*rmu))) )

    delT= min(delT,cfl*(delX/( maxval(( 9.8_dp*( depth))**0.5_dp+abvel))))

END SUBROUTINE time_step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Dat

INTEGER:: Date_Time(8)
Character::a(3)
Call DATE_AND_TIME(a(1),a(2),a(3), Date_Time)

Print * ,Date_Time

End Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mouth_height(th, t, tr, mouth_data, mouthread, mouth_data_len)
    INTEGER, INTENT(IN):: mouth_data_len
    REAL(dp), INTENT(IN):: t, tr
    REAL(dp), INTENT(IN)::mouth_data(mouth_data_len,2)
    REAL(dp), INTENT(INOUT):: th
    LOGICAL, INTENT(IN):: mouthread

    INTEGER:: lwind
    REAL(dp):: tmax, tscale, dtscale
    !!You could pass another argument here - mouthdta - a matrix with 2 columns -
    !time and the height - for interpolation

    IF(mouthread.eqv..false.) THEN
        !th=  (tr/2._dp)*( sin(2.00_dp*pi*t/(3600.00_dp*12.40_dp)) ) !*atan(counter/4000.)*2./pi
        th=  (tr/2._dp)*( sin(2.00_dp*pi*t/(180.0_dp)) ) + 0.082 !Tambrioni experiment 1 !*atan(counter/4000.)*2./pi
        !if(t< 3260._dp) THEN
        !th= .05_dp*(2._dp/2.00_dp)*( sin(2.00_dp*pi*t/(1.250_dp*3260._dp)))!(3600.00_dp*12.40_dp)) )+0._dp !*atan(counter/4000.)*2./pi
    ELSE 

        !        th=0._dp
        !END IF
        tmax= mouth_data(mouth_data_len, 1) !The maximum time in the input data. 
        !tscale = mod(t, tmax) !t modulo the max time - so we loop after this
        IF(tmax<t) THEN
            print*, 'PROBLEM: t is greater than the maximum time in the &
                    downstream boundary condition input', t, tmax
            stop
        END IF
        !
        dtscale= mouth_data(2,1)-mouth_data(1,1) !The time increment in the data
        !
        lwind= floor(t/dtscale)+1 !The index corresponding to tscale
        !
        th= (( t-mouth_data(lwind,1))*mouth_data(lwind+1,2) + (mouth_data(lwind+1,1)-t)*mouth_data(lwind,2))/& 
        (mouth_data(lwind+1,1)-mouth_data(lwind,1)) !The height at the mouth

    END IF

END SUBROUTINE mouth_height
!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine hy22(delT, delX, h, B, A, Q,Q2H,Y, t,l, counter,dbdh, rmu,inuc, LAKE, hlim, Qb,tr, & 
!mouth_data,mouthread, mouth_data_len, mns, rho, g,cfl)!slL,srL,sLlh,sLrh,C0, qbbl,export)
!INTEGER, Intent (In)::l, mouth_data_len
!REAL(dp), Intent(In):: delX, Qb,tr, mouth_data, rho, mns, g, cfl,Q2H
!REAL(dp),Intent(In out)::t
!REAL(dp), Intent(In out):: h, B,A,Q, Y!, qbbl!,slL,srL,sLlh,sLrh, qbbl,C0
!DIMENSION:: h(l), B(l),A(l),Q(l),Q2H(0:l), Y(l), mouth_data(mouth_data_len,2), mns(l)!, qbbl(l)!,slL(l),srL(l),sLlh(l),sLrh(l), qbbl(l) 
!REAL(dp), Intent(In out)::delT!,export
!Integer, intent (in):: counter
!REAL(dp), Intent(in)::dbdh(l,2), rmu(l),inuc(l)
!Logical, intent(in)::LAKE, mouthread
!
!!Local variables-- this is a mess, you really need to look at the code to see what they mean, I just wacked them in as I needed them
!!logical:: isnan
!REAL(dp) :: Qblast98,Ablast98,one=1.00_dp, morpfact, hlim
!REAL(dp)::d0,bb,cc,tmod,r1,r2,vel,KonstB1p1, KonstB1p2,KonstA1, KonstA100, KonstB1,KonstB100,blastminus
!REAL(dp)::KonstB1p3, KonstB100p1, KonstB100p2, KonstB100p3, KonstB1p4, KonstB100p4, sfirstplus, slastminus, hfirstplusb
!REAL(dp):: Qfirstplusb, hlastminusb, Qlastminusb,tcrit, Qm, t1, t2, D50,th,dlim,fillim, filladd
!INTEGER::s1,s2
!INTEGER:: i,  j, track, swit,lo
!REAL:: upw(l)  
!REAL(dp):: Ab,v, Qml, bconst,savvy,savvy2,bcrit, qs, echar,mxdep,aa,arr(l),derf,eps=epsilon(1.00_dp)
!REAL(dp)::Ashift(l),Abshift(l),Qshift(l),Qbshift(l),Bshift(l),Bbshift(l), Bs(l),d2hdx2(l), entil(l), rat(l), ratlast(l)
!REAL(dp)::Alast(l),Qlast(l), df(l), vf(l), dc(l), vc(l), randomunif(l),Bpred(l), Bcor(l), dscale(l), epscale(l),zer(l)
!REAL(dp)::Apred(l),Qpred(l),Qpred1(l),Qpred2(l), Qpred3(l), Apredc1(l), Apredc2(l), er(l), fact(l),bshear(l), efd(l),efd2(l)
!REAL(dp)::Apredbshift(l), Apredfshift(l), Apredp1(l), Apredp2(l), Qpredbshift(l), Qpredfshift(l), Qpredp(l), Qpredc(l), hbshift(l)
!REAL(dp)::Acor(l),Qcor(l),Qcor1(l),Qcor2(l),Qcor3(l), Qcor4(l), alpha(l), beta(l), Af(l),Ac(l), Ap(l), Acorb(l), Acorb2(l)
!REAL(dp)::delHdelXp(l), delHdelXc(l), delBdelXp(l), delBdelXc(l), Ap1(l), Ac1(l),Ac2(l), Qp(l), Qc(l), Bp(l), Bc(l),efdshift(l) 
!REAL(dp):: taub(l), qbbp(l),qbbc(l), hpred(l),hcor(l),Ablast(l),Qblast(l),Areafil(l), Qfil(l),hshift(l), sih(l), bdist(l), slop(l)
!REAL(dp)::Ylast(l),Ypred(l),Ycor(l), Yshift(l),Ybshift(l),sihmin(l),efdmin(l),Abb(l), tok, fric(l), ccc(l), YY(l), QQ(l), zer2(l)
!REAL(dp)::Bnew(l), inucshift(l), s(l), up(l), low(l), diag(l), Dfil(l),Ypredbshift(l), visc(l), vcoef(l)!, QQ(l)
!
!LOGICAL::Limiter=.false., bound=.false.   !these variables control the use of the limiter and boundery treatment, and the simplified flat water treatment 
!!REAL(dp):: cfl=.95_dp  !Note that we are running off a grid that is half the size of delX, and taking 2 time steps, so the cfl should be the same as normal
!!!!NOTE-- it would be really cool to see if the morphological results were
!!effected by the flood or ebb dominance. Everyone assumes that they are,
!!although in general there might not be such a major difference. 
!
!!!Now, it seems that in some situations they can create morphodynamic problems.
!!Still need to fix the method of characteristics at the landward boundary
!!And at the ocean--I'm still not sure that the retrying for xfrstplus is
!!working-- sometimes it takes a good value but still fails.  
!    
!!!IMPORTANT VARIABLES
!
!!Qb= -160. !landward boundary discharge
!!qs= .00/100000.00*(Qb) !Upstream sediment input
!!tr=4.0 !Tide range
!!tcrit= 1000.047 !Dimensionless critical shear for moving bed sediment
!!bcrit=100000.047 !critical shear for bank erosion, see Julian and Torres (2006)  
!!D50=0.0001
!!morpfact= 400.00 !Morphological factor to accellerate convergence. 
!!bconst=0.000001
!!hlim= 0.005
!!IF(hlim>0._dp) THEN
!!        dlim=max(maxval(.1001_dp*delX*.5_dp*rmu*g), hlim+0.001_dp)/hlim !0.0601/hlim  !discharge limit (=how many times hlim) where discharge is 0 if height is < dlim 
!!ELSE
!        dlim=1._dp
!!END IF
!!print*, dlim*hlim
!!!!!!!!Note the paper Burguete et al (2007) where they provide a 'stability
!!criterion' for wetting and drying that basically says that delX< 2*depth/(f/8)
!!always. Now, in my case, with f/8 ~0.002 or something like that, this means
!!that the depth should be > ~ 5cm. I found indeed that the stability of the
!!drying front improved lot when doing this. Yep, sure does
!
!!print*, dlim
!
!
!!Make sure that the initial state conforms to the limits on Y
!swit=l !Predefine the wet boundary
!
!!         DO j= 2, l
!!         IF(Y(j)<mns(j)+dlim*hlim) Q(j)=0._dp
!!        
!!         IF(Y(j)<mns(j)+hlim) THEN
!!       Y(j)= h(j) + hlim/1.001_dp !(((dlim-1)/(j-(i-1)) +1)*hlim) !-arr(i) !Forcing the elevation seems to
!!        Q(j)=0._dp
!!         END IF
!!         END DO
!      
! efd= Y-h!max(Y-h,hlim+0.*h)
!
!A=B*efd! cross sectional area. Note that this is reset here, and so the input area actually doesn't effect the outcome. Is this really what you want?? 
!
!!IF (t<0.000001_dp) then !set some variables
!!	        counter=0
!!        END IF
!!       
!!       
!!        counter=counter+1
!
!
!
!
!Alast=A
!Qlast=Q
!
!Ylast=Y
!
!!!!!!!!!!!!!!!!!!!!TIME STEP
!delT=1000._dp
!!mxdep=maxval(efd)
!!aa=maxval(abs(Q/A))
!
!!Time step
!!delT= min(delT,cfl*((delX/( ( 9.8_dp*( mxdep ))**0.5+aa+0.000_dp*delX))), & 
!!minval(.5_dp*minval(efd)/((0.001_dp+ abs(Q/A))*(9.8_dp*rmu))) ) !I think the latter limit comes from Burguete?
!
!call time_step(delT,cfl, delX, Y-h, abs(Q/A), rmu, l )
!
!!print*, delT
!!print*, minval(rmu)
!
!!Checks
!IF(abs(delT-1000._dp)<1.1_dp) print*, "delT= 1000."
!
!IF(delT<1.0E-04_dp) THEN 
!        print*, "short time step", delT
!        print*, "counter=", counter, "maxloc=", maxloc(abs(Q/A)), "l=", l, "max=", maxval(abs(Q/A)), 1
!        print*,  Q/A
!        !print*, maxval(h), minval(h)!maxval(efd), maxval(abs(Q)/A)
!stop
!END IF
!
!
!t=t+delT/2._dp !time in seconds
!
!!!!!!!!!Flag for nan time
!        IF (isnan(t)) THEN
!        call Dat
!        print*, " Error: time is NAN", Y, h, efd
!        Stop
!        END IF
!!!!!!!!!!!!!!!
!!TIDE --Now we find the change in water depth at the mouth. 
!
!!th=  (tr/2.00_dp)*( sin(2.00_dp*pi*t/(3600.00_dp*12.40_dp)) )!*atan(counter/4000.)*2./pi  
!call mouth_height(th, t,tr, mouth_data,mouthread, mouth_data_len)
!d0=th !+(tr/2.)*(counter)**(-0.01)! A simple tide !max(th, hlim+h(1)) ! need to fix this to be depth at the mouth. 
!
!!print*, d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!! water level, lax-friedrichs approach
!!!SO here the "pred " 's are the first half step, and the "cor's" are the second
!!half step. See Liska and Wendroff. 
!Ypred=minval(h)-1._dp
!!!First half step. Y(2)= grid position 3/2
!DO i= 2, l
!IF(Y(i)>mns(i)+hlim) THEN  !Here, we don't want dry points with a bedslope to be effected by higher upstream or lower downstream points
!!Ypred(i)= Ylast(i) + (delT/delX)*( .5*(1._dp+sign(1._dp, Q(i)))*(Q(i)-Q(i-1)) &
!!+ .5*(1._dp-sign(1._dp, Q(i)))*(Q(i+1)-Q(i)))/B(i) 
!
!Ypred(i)= .5_dp*(Ylast(i-1)+Ylast(i)) -(delT/(2._dp*delX))*(Q(i)-Q(i-1))/(.5_dp*(B(i)+B(i-1)))
!
!END IF
!        IF((Y(i-1)>(mns(i-1))+hlim).and.(Y(i)<mns(i)+hlim)) THEN
!  
!     Ypred(i)=max(1._dp*Y(i-1)-0._dp*(delT/(2._dp*delX))*(0._dp-Q(i-1))/(.5_dp*(B(i)+B(i-1))),.5_dp*(h(i)+h(i-1))+hlim/1.001_dp) !Set it to the water level at the end point
!        END IF
!
!        IF(Ypred(i)<.5_dp*(mns(i)+mns(i-1))+hlim)THEN
!                
!        
!                Ypred(i)=.5_dp*(h(i)+h(i-1)) +hlim/1.001_dp ! Set it to the default dry water level
!        END IF
!
!
!END DO
!
!
!
!!print*, "yes"
!!!!!!!!!!!Boundary condition
!       ! Ypred(l)=Y(l) 
!        Ypred(1)= d0!+ (g*(Y(1)-h(1)))**.5 -(g*(Y(2)-h(1)))**.5 ! -(d0-Ypred(2))!Y(1)- .5*(Y(2)-Y(1)) ! This is Y(0.5)
!
!
!
!Apred(2:l)=.5_dp*(B(2:l)+B(1:l-1))*(Ypred(2:l) -.5_dp*(h(1:l-1)+h(2:l) ) ) 
!Apred(1)= .5_dp*(B(1)+B(1))*(Ypred(1)-.5_dp*(h(1)+h(1)))!2.*Apred(2)-Apred(3)!Ypred(1)-h(1) !Extrapolate
!
!
!!Define the slope
!  do i=2, l
!  slop(i)= Y(i)-Y(i-1)!sign(1._dp,(Y(i)-Y(i-1)))*min(abs(Y(i)-Y(i-1)),0.03) !Same indexing as Ypred 
!  !slop(l)=slop(l-1)
!  end do
!
!slop(1)=Y(1)-(2._dp*Y(1)-Y(2)) !slop(2)!Y(2)-Y(1)! Slop(1/2)
!!slop(l)=0.
!!!Now, we need some way to regulate the dry points. 
!  DO i=1, l
! IF(Y(i)<mns(i)+dlim*hlim) THEN
!      !  slop(i)=slop(i-1)
!         slop((i):l)=0._dp
! !slop(i:i+1)=slop(i-1:i)        
! !slop(i)=0.!slop(i-1)
! swit=i
! goto 20 
!
! END IF
! 
! END DO
!
!
!!DISCHARGE first half. Indexing as for Ypred
!20 DO i= 2, l
!Qpred1(i)= .5_dp*(Qlast(i-1)+Qlast(i))-(delT/(2._dp*delX))*(Q(i)**2._dp/A(i)-Q(i-1)**2._dp/A(i-1) +inuc(i)-inuc(i-1))-& 
!(delT/(2._dp*delX))*g*(((A(i)+A(i-1))*.5_dp*(slop(i)))) !
!end do
!
!
!!Boundaries-- Note the extrapolation here. 
!Qpred1(1)= .5_dp*(Qlast(1)+Qlast(1))-(delT/(2._dp*delX))*( Q(1)**2._dp/A(1)-(2._dp*Q(1)-Q(2))**2._dp/(2._dp*A(1)-A(2)) )-& 
!(delT/(2._dp*delX))*g*((.5_dp*(A(1)+A(1))*(slop(1))))
!!Qpred1(1)=0.!2.*Qpred1(2)-Qpred1(3)
!!Setting Qpred1 to zero seems to stop the growing oscillations. But it's not
!!very physical
!!Qpred1(1)= .5*Qlast(1) + ((Qpred1(2)-.5*(Qlast(1)+Qlast(2)))/1.) !*(1.5*A(1)-.5*A(2))/(.5*A(2)+.5*A(1)) +& 
!!Apred(1)*.5*(Qlast(1)/Alast(1)+(2.*Qlast(1)/Alast(1)-Qlast(2)/Alast(2)))!Apred(1)*(Qpred1(2)/Apred(2) )!-2.*(g*(Ypred(2)-h(1)))**.5 + 2.*(g*(Ypred(1)-h(1)))**.5)
!
!!!Explicit friction
!
!!DO i= 2, l
!!Qpred2(i)= -(g*.5*(A(i)+A(i-1))*.5*(Q(i)+Q(i-1))*abs(.5*(Q(i-1)+Q(i)))/((.5*(A(i)+A(i-1)))**2. & 
!!*.5*(Y(i)+Y(i-1)-h(i)-h(i-1)) ))*rmu(i) 
!!end do
!!Qpred2(1)= -(g*.5*(A(1)+A(1))*Q(1)*abs(Q(1))/(A(1)**2.*(Y(1)-h(1)) ))*rmu(1) 
!!Qpred2(l)=0.01
!
!!Implicit friction
!DO i= 2, l
!Qpred2(i)= -(g*.5_dp*(A(i)+A(i-1))*sign(1._dp,Qpred1(i))/(Apred(i)**2._dp*(Ypred(i)-.5_dp*(h(i)+h(i-1)) ) )) & 
!*.5_dp*(rmu(i)+rmu(i-1)) 
!end do
!Qpred2(1)= -(g*.5_dp*(A(1)+A(1))*sign(1._dp, Qpred1(1))/ & 
!(Apred(1)**2._dp*(Ypred(1)-.5_dp* (h(1)+h(1))) &
!))*.5_dp*(rmu(1)+ rmu(1)) 
!!Qpred2(l)=0.01
!
!
!
!visc=0._dp !
!
!!!Explicit friction
!!Qpred= Qpred1+(delT/2.)*Qpred2
!
!!Implcit friction
!Qpred= (1._dp - sqrt(1._dp- 4._dp*(delT/2._dp)*Qpred2*(Qpred1+visc) ))/(2._dp*(delT/2._dp)*Qpred2)
!
!!Qpred(1)=2.*(1.*Q(1)-(Qpred(2)))
!!Qpred(1)= Qlast(1)!Q(1)- .5*(Q(2)-Q(1))
!
!!Qpred(1)=Apred(1)*(Qlast(1)/Alast(1) - (g*(Y(1)-h(1)))**.5 +(g*(Ypred(1)-h(1)))**.5)
!!Boundary conditions
!!Qpred(l)= Q(l)!Qb !- 2.*Apred(l)*((g*Ypred(l)-h(l))**0.5- (g*(Ypred(l-1)-h(l-1)))**0.5)!+ Qpred(l-1) !Qlast(l)!2*Qpred(l-1)-Q(l-2)
!
!!Ypred(l)= Y(l)!((Q(l-1)/A(l-1)- Qb/A(l) + (g*(Ypred(l-1)-h(l-1)))**0.5)**2.)/g +h(l)
!
!!print*, Qpred-Q
!!stop
!
!!Special treatment of boundary
!DO i= 2, l
!IF(isnan(Qpred(i))) THEN
!        print*, "Qpred is nan", i, Qpred1(i), Qpred2(i), 1._dp- 4._dp*(delT/2._dp)*Qpred2(i)*(Qpred1(i)+visc(i)) 
!        stop
!end IF
!IF((Ypred(i)<.5_dp*(mns(i)+mns(i-1))+dlim*hlim)) THEN 
!Qpred(i:l)=0._dp
!goto 4122
!END IF
!END DO
!
!
!!4122 !mxdep=maxval(efd)
!!aa=maxval(abs(Qpred/Apred))
!
!!Time step
!!delT= min(delT,cfl*((delX/( ( 9.8*( mxdep ))**0.5+aa+0.000*delX))), minval(.5*hlim/(abs(Qpred/Apred)*(9.8*rmu/C0**2.))) )
!
!4122 t=t+delT/2._dp !Step through the time
!!TIDE --Now we find the change in water depth at the mouth. 
!
!!th=  (tr/2.00_dp)*( sin(2.00_dp*pi*t/(3600.00_dp*12.40_dp)) )!*atan(counter/4000.)*2./pi  
!call mouth_height(th, t,tr, mouth_data,mouthread, mouth_data_len)
!d0=th !+(tr/2.)*(counter)**(-0.01)! A simple tide !max(th, hlim+h(1)) ! need to fix this to be depth at the mouth. 
!
!!Now we update the width, eh. 
!Bnew(2:l-1)= B(2:l-1) + 0._dp*(.25_dp*(B(1:l-2)+2._dp*B(2:l-1)+B(3:l) ) ) &
!+dbdh(2:l-1,1)*(.5_dp*(Ypred(2:l-1)+Ypred(3:l))-Y(2:l-1)) 
!Bnew(l)=  B(l)+dbdh(l,1)*(Ypred(l)-Y(l))
!Bnew(1)= B(1)+dbdh(1,1)*(Ypred(1)-Y(1))
!!!!!!!!!!!!!!!!
!!!!!!!!Elevation corrector
!
!Ycor=minval(h)-1.
! DO i= 2, l-1
! IF(Bnew(i)<0._dp) Bnew(i)=B(i) !Handy fix so we don't push B to zero with the derivative
!IF(Ypred(i+1)>.5_dp*(mns(i)+mns(i+1))+hlim) THEN  !Here, we don't want dry points with a slope to be effected by higher upstream or lower downstream points
!
!Ycor(i)= .5_dp*(Ypred(i)+Ypred(i+1)) - (delT/(2._dp*delX))*(Qpred(i+1)-Qpred(i))/(.25_dp*(0._dp*B(i-1)+4._dp*B(i)+0._dp*B(i+1)))!( Bnew(i) ) !!Note treatment of dbdh. 
!End if
!        IF ((Ypred(i)>.5_dp*(mns(i)+mns(i-1))+hlim).and.(Ypred(i+1)<.5_dp*(mns(i)+mns(i+1))+hlim)) THEN !The water level at i-1/2 is greater than the dry level at the ith point
!	Ycor(i)=max(1._dp*Ypred(i) -0._dp*(delT/(2._dp*delX))*(0._dp-Qpred(i))/(.25_dp*(0._dp*B(i-1)+4._dp*B(i)+0._dp*B(i+1) )),&
!       	h(i)+hlim/1.001_dp)
!        END IF
!       
!! IF((Y(i-1)>(mns(i-1))+hlim).and.(Y(i)<mns(i)+hlim)) THEN
! ! 	Ycor(i)= max(Y(i-1), h(i)+hlim/1.001_dp)
!!	END IF
!
!        IF(Ycor(i)<h(i)+hlim) THEN
!        Ycor(i)=h(i)+hlim/1.001_dp !.5*(h(i)+Ypred(i+1))
!        END IF
!    
!END DO
!
!
!!!Boundary conditions
!Ycor(1)=d0
!
!Ycor(l)= 1.5_dp*Ypred(l)-.5_dp*Ypred(l-1) - (delT/(2._dp*delX))*(Qb-Qpred(l))/Bnew(l) !max(Ypred(l), h(l)+hlim/1.001)
!efd= Ycor-h!
!Acor(1:l-1)=(Bnew(1:l-1))*efd(1:l-1)  !Note treatment of dbdh
!Acor(l)= (Bnew(l))*efd(l)
!
!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!Discharge cor
!
!!Define the slope
!  do i=1, l-1
!  slop(i)= Ypred(i+1)-Ypred(i) !sign(1._dp,(Ypred(i+1)-Ypred(i)))*min(abs(Ypred(i+1)-Ypred(i)),.03)
!  !slop(l)=slop(l-1)
!  end do
!
!slop(l)=slop(l-1)!2.*slop(l-1)-slop(l-2)
!
!!!Now, we need some way to regulate the dry points. 
!  DO i=1, l-1
!IF(Ypred(i+1)<.5_dp*(mns(i+1)+mns(i))+dlim*hlim) THEN
!       ! slop(i)=slop(i-1)
!        slop(i:l)=0._dp
! !slop(i:i+1)=slop(i-1:i)        
! !slop(i)=0.!slop(i-1)
! swit=i
! goto 30
!
! END IF
! 
! END DO
!
!!DISCHARGE second half. Indexing as for Ypred
!30 DO i= 2, l-1
!
!Qcor1(i)= .5_dp*(Qpred(i)+Qpred(i+1))-(delT/(2._dp*delX))*(Qpred(i+1)**2._dp/Apred(i+1)-Qpred(i)**2._dp/Apred(i) &
!+ (inuc(i)-inuc(i-1)))-(delT/(2._dp*delX))*g*((.5_dp*(Apred(i)+Apred(i+1))*(slop(i)))) !
!
!end do
!
!!Needed due to the inuc term
!Qcor1(1)= .5_dp*(Qpred(1)+Qpred(1+1))-(delT/(2._dp*delX))*(Qpred(1+1)**2._dp/Apred(1+1)-Qpred(1)**2._dp/Apred(1)+ & 
!1._dp*(inuc(1+1)-inuc(1)))- (delT/(2._dp*delX))*g*((.5_dp*(Apred(1)+Apred(1+1))*(slop(1)))) !
!
!
!!Qcor(1)= Qcor1(1)-.5*(Qpred(1)-Qpred(2))
!!Qcor1(l)=0.!Qcor1(l-1)
!
!!Qcor1(l)= .5*(Qpred(l)+Qb)-(delT/(2.*delX))*(Qb**2./Apred(l)-Qpred(l)**2./Apred(l))-& 
!!(delT/(2.*delX))*g*((.5*(Apred(l)+Apred(l))*(slop(l)))) !
!
!!!Explicit friction
!
!!DO i= 1, l-1
!!Qcor2(i)= -(g*.5*(Apred(i)+Apred(i+1))*.5*(Qpred(i)+Qpred(i+1))*abs(.5*(Qpred(i+1)+Qpred(i)))/&
!!((.5*(Apred(i)+Apred(i+1)))**2.*.5*(Ypred(i)+Ypred(i+1)-.5*h(i-1) -h(i)-.5*h(i+1)) ))*rmu(i) 
!!end do
!!Qcor2(l)=Qcor2(l-1)!0.01
!
!!!Implciti froction
!DO i= 1, l-1
!Qcor2(i)= -(g*.5_dp*(Apred(i)+Apred(i+1))*sign(1._dp,Qcor1(i))/(Acor(i)**2.*(efd(i) ) ))*rmu(i) 
!end do
!!Qcor2(l)=0.01
!!Qcor2(l)= -(g*.5*(Apred(l)+Apred(l))*sign(1._dp,Qb)/(Acor(l)**2.*(efd(l) ) ))*rmu(l) 
!
!
!!Explicit
!!Qcor=Qcor1+(delT/2.)*Qcor2
!!visc(2:l-1)= 0. !-.01*(Q(1:l-2)-2.*Q(2:l-1)+Q(3:l-1))
!!implicit
!
!
!
!Qcor= (1. - sqrt(1._dp- 4._dp*(delT/2._dp)*Qcor2*(Qcor1+visc) ))/(2._dp*(delT/2._dp)*Qcor2)
!
!
!
!Qcor(l)=Qb 
!
!!QQ=Qcor
!!DO i= 2, l-1
!!Qcor(i)= Qcor(i)+ 0.01*(QQ(i-1)-2.*QQ(i)+QQ(i+1))
!!END DO
!!print*,Ycor-Y
!!print*, "Next"
!!print*, Qpred
!!print*, "Next"
!!print*, Qcor
!!stop
!
!DO i=1, l
!IF(isnan(Qcor(i))) THEN
!        print*, "Qcor is nan", i, Qcor1(i), Qcor2(i), Acor(i), rmu(i)
!        stop
!end IF
!IF(Ycor(i)<mns(i)+dlim*hlim) THEN
!Qcor(i:l)=0.
!goto 333 
!END IF
!END DO
!
!333 DO i= 2, l
!IF( Ypred(i)<.5*(h(i)+h(i-1)) )THEN 
!        print*, "Ypred is neg", i
!        stop
!END IF
!
!IF( (Ycor(i)<h(i)).or.(isnan(Ycor(i))))THEN 
!        print*, "Ycor is a neg depth or nan", i, Ycor(i), h(i)
!        stop
!END IF
!
!
!END DO
!
!Y=Ycor
!Ycor(1)=d0
!A=B*(Ycor-h)
!Q=Qcor
!
!
!
!
!
!
!
!end subroutine hy22
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine hy23(delT, delX, h, B, A, Q,Q2,Y, t,l, counter,dbdh, rmu,inuc, LAKE, hlim, Qb,tr, & 
!mouth_data,mouthread, mouth_data_len, mns, rho, g,cfl)!slL,srL,sLlh,sLrh,C0, qbbl,export)
!
!INTEGER, Intent (In)::l, mouth_data_len
!REAL(dp), Intent(In):: delX, Qb,tr, mouth_data, rho, mns, g, cfl,Q2
!REAL(dp),Intent(In out)::t
!REAL(dp), Intent(In out):: h, B,A,Q, Y!, qbbl!,slL,srL,sLlh,sLrh, qbbl,C0
!DIMENSION:: h(l), B(l),A(l),Q(l),Q2(l), Y(l), mouth_data(mouth_data_len,2), mns(l)!, qbbl(l)!,slL(l),srL(l),sLlh(l),sLrh(l), qbbl(l) 
!REAL(dp), Intent(In out)::delT!,export
!Integer, intent (in):: counter
!REAL(dp), Intent(in)::dbdh(l,2), rmu(l),inuc(l)
!Logical, intent(in)::LAKE, mouthread
!
!!Local variables-- this is a mess, you really need to look at the code to see what they mean, I just wacked them in as I needed them
!!logical:: isnan
!REAL(dp) :: Qblast98,Ablast98,one=1.00_dp, morpfact, hlim
!REAL(dp)::d0,bb,cc,tmod,r1,r2,vel,KonstB1p1, KonstB1p2,KonstA1, KonstA100, KonstB1,KonstB100,blastminus
!REAL(dp)::KonstB1p3, KonstB100p1, KonstB100p2, KonstB100p3, KonstB1p4, KonstB100p4, sfirstplus, slastminus, hfirstplusb
!REAL(dp):: Qfirstplusb, hlastminusb, Qlastminusb,tcrit, Qm, t1, t2, D50,th,dlim,fillim, filladd
!INTEGER::s1,s2
!INTEGER:: i,  j, track, swit,lo, ind
!REAL:: upw(l), w1  
!REAL(dp):: Ab,v, Qml, bconst,savvy,savvy2,bcrit, qs, echar,mxdep,aa,arr(l),derf,eps=epsilon(1.00_dp)
!REAL(dp)::Ashift(l),Abshift(l),Qshift(l),Qbshift(l),Bshift(l),Bbshift(l), Bs(l),d2hdx2(l), entil(l), rat(l), ratlast(l)
!REAL(dp)::Alast(l),Qlast(l), df(l), vf(l), dc(l), vc(l), randomunif(l),Bpred(l), Bcor(l), dscale(l), epscale(l),zer(l)
!REAL(dp)::Apred(l),Qpred(l),Qpred1(l),Qpred2(l), Qpred3(l), Apredc1(l), Apredc2(l), er(l), fact(l),bshear(l), efd(l),efd2(l)
!REAL(dp)::Apredbshift(l), Apredfshift(l), Apredp1(l), Apredp2(l), Qpredbshift(l), Qpredfshift(l), Qpredp(l), Qpredc(l), hbshift(l)
!REAL(dp)::Acor(l),Qcor(l),Qcor1(l),Qcor2(l),Qcor3(l), Qcor4(l), alpha(l), beta(l), Af(l),Ac(l), Ap(l), Acorb(l), Acorb2(l)
!REAL(dp)::delHdelXp(l), delHdelXc(l), delBdelXp(l), delBdelXc(l), Ap1(l), Ac1(l),Ac2(l), Qp(l), Qc(l), Bp(l), Bc(l),efdshift(l) 
!REAL(dp):: taub(l), qbbp(l),qbbc(l), hpred(l),hcor(l),Ablast(l),Qblast(l),Areafil(l), Qfil(l),hshift(l), sih(l), bdist(l), slop(l)
!REAL(dp)::Ylast(l),Ypred(l),Ycor(l), Yshift(l),Ybshift(l),sihmin(l),efdmin(l),Abb(l), tok, fric(l), ccc(l), YY(l), QQ(l), zer2(l)
!REAL(dp)::Bnew(l), inucshift(l), s(l), up(l), low(l), diag(l), Dfil(l),Ypredbshift(l), visc(l), vcoef(l)!, QQ(l)
!
!LOGICAL::Limiter=.false., bound=.false.   !these variables control the use of the limiter and boundery treatment, and the simplified flat water treatment 
!!REAL(dp):: cfl=.95_dp  !Note that we are running off a grid that is half the size of delX, and taking 2 time steps, so the cfl should be the same as normal
!!!!NOTE-- it would be really cool to see if the morphological results were
!!effected by the flood or ebb dominance. Everyone assumes that they are,
!!although in general there might not be such a major difference. 
!
!!!Now, it seems that in some situations they can create morphodynamic problems.
!!Still need to fix the method of characteristics at the landward boundary
!!And at the ocean--I'm still not sure that the retrying for xfrstplus is
!!working-- sometimes it takes a good value but still fails.  
!    
!!!IMPORTANT VARIABLES
!
!!Qb= -160. !landward boundary discharge
!!qs= .00/100000.00*(Qb) !Upstream sediment input
!!tr=4.0 !Tide range
!!tcrit= 1000.047 !Dimensionless critical shear for moving bed sediment
!!bcrit=100000.047 !critical shear for bank erosion, see Julian and Torres (2006)  
!!D50=0.0001
!!morpfact= 400.00 !Morphological factor to accellerate convergence. 
!!bconst=0.000001
!!hlim= 0.005
!!IF(hlim>0._dp) THEN
!!        dlim=max(maxval(.1001_dp*delX*.5_dp*rmu*g), hlim+0.001_dp)/hlim !0.0601/hlim  !discharge limit (=how many times hlim) where discharge is 0 if height is < dlim 
!!ELSE
!        dlim=1._dp
!!END IF
!!print*, dlim*hlim
!!!!!!!!Note the paper Burguete et al (2007) where they provide a 'stability
!!criterion' for wetting and drying that basically says that delX< 2*depth/(f/8)
!!always. Now, in my case, with f/8 ~0.002 or something like that, this means
!!that the depth should be > ~ 5cm. I found indeed that the stability of the
!!drying front improved lot when doing this. Yep, sure does
!
!!print*, dlim
!
!
!!Make sure that the initial state conforms to the limits on Y
!swit=l !Predefine the wet boundary
!
!!         DO j= 2, l
!!         IF(Y(j)<mns(j)+dlim*hlim) Q(j)=0._dp
!!        
!!         IF(Y(j)<mns(j)+hlim) THEN
!!       Y(j)= h(j) + hlim/1.001_dp !(((dlim-1)/(j-(i-1)) +1)*hlim) !-arr(i) !Forcing the elevation seems to
!!        Q(j)=0._dp
!!         END IF
!!         END DO
!      
!! efd= Y-h!max(Y-h,hlim+0.*h)
!
!!A=B*efd! cross sectional area. Note that this is reset here, and so the input area actually doesn't effect the outcome. Is this really what you want?? 
!
!!IF (t<0.000001_dp) then !set some variables
!!	        counter=0
!!        END IF
!!       
!!       
!!        counter=counter+1
!
!
!
!
!Alast=A
!Qlast=Q
!
!Ylast=Y
!
!!!!!!!!!!!!!!!!!!!!TIME STEP
!delT=1000._dp
!!mxdep=maxval(efd)
!!aa=maxval(abs(Q/A))
!
!!Time step
!!delT= min(delT,cfl*((delX/( ( 9.8_dp*( mxdep ))**0.5+aa+0.000_dp*delX))), & 
!!minval(.5_dp*minval(efd)/((0.001_dp+ abs(Q/A))*(9.8_dp*rmu))) ) !I think the latter limit comes from Burguete?
!
!call time_step(delT,cfl, delX, Y-h, abs(Q/A), rmu, l )
!
!!print*, delT
!!print*, minval(rmu)
!
!!Checks
!IF(abs(delT-1000._dp)<1.1_dp) print*, "delT= 1000."
!
!IF(delT<1.0E-04_dp) THEN 
!        print*, "short time step", delT
!        print*, "counter=", counter, "maxloc=", maxloc(abs(Q/A)), "l=", l, "max=", maxval(abs(Q/A)), 1
!        print*,  Q/A
!        !print*, maxval(h), minval(h)!maxval(efd), maxval(abs(Q)/A)
!stop
!END IF
!
!
!t=t+delT!/2._dp !time in seconds
!
!!!!!!!!!Flag for nan time
!        IF (isnan(t)) THEN
!        call Dat
!        print*, " Error: time is NAN", Y, h, efd
!        Stop
!        END IF
!!!!!!!!!!!!!!!
!!TIDE --Now we find the change in water depth at the mouth. 
!
!!th=  (tr/2.00_dp)*( sin(2.00_dp*pi*t/(3600.00_dp*12.40_dp)) )!*atan(counter/4000.)*2./pi  
!call mouth_height(th, t,tr, mouth_data,mouthread, mouth_data_len)
!d0=th !+(tr/2.)*(counter)**(-0.01)! A simple tide !max(th, hlim+h(1)) ! need to fix this to be depth at the mouth. 
!
!!print*, d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!do i=2,l-1
!Apred(i)=0.5_dp*(A(i-1)+A(i+1)) - (delT/delX)*(Q(i+1)-Q(i-1))*0.5_dp
!end do
!Apred(l)=0.5_dp*(A(l-1)+A(l)) -(delT/delX)*(Qb-Q(l-1))*0.5_dp !Boundary cond, assuming zero gradient in A upstream
!Apred(1)=0.5_dp*(A(2)+A(1)) - (delT/delX)*(Q(2)-Q(1))*0.5_dp !Boundary cond, assuming zero gradient in A downstream
!
!!Use Apred to calculate Ypred
!do i=1,l
!        if(Apred(i)-Alast(i)>0._dp) THEN
!        ind=1
!        else
!        ind=2
!        end if
!
!w1= B(i)**2+2._dp*dbdh(i,ind)*(Apred(i)-Alast(i))
!if((dbdh(i,ind)>0._dp).and.(w1>0._dp)) THEN
!Ypred(i)=Ylast(i)+ (-B(i)+ sqrt(w1))/(dbdh(i,ind))
!else
!	Ypred(i)= Ylast(i)+ (Apred(i)-Alast(i))/B(i)
!END IF
!end do
!!Momentum equation, internal points
!do i=2,l-1
!Qpred1(i)=0.5_dp*(Q(i-1)+Q(i+1)) - (delT/delX)*( &
! ((1._dp+inuc(i+1))*Q(i+1)**2/A(i+1) -(1._dp+inuc(i-1))*Q(i-1)**2/A(i-1) )*0.5_dp &
! +g*A(i)*(Y(i+1)-Y(i-1))*0.5_dp &
!)
!end do
!!Upstream boundary, discharge is Qb, area is A(l), water elevation is Y(l)
!Qpred1(l)=0.5_dp*(Q(l-1)+Qb) - (delT/delX)*( &
! ((1._dp+inuc(l))*Qb**2/A(l) -(1._dp+inuc(l-1))*Q(l-1)**2/A(l-1) )*0.5_dp &
! +g*A(l)*(Y(l)-Y(l-1))*0.5_dp &
!)
!Qpred1(1)=0.5_dp*(Q(1)+Q(2)) - (delT/delX)*( &
! ((1._dp+inuc(1+1))*Q(1+1)**2/A(i+1) -(1._dp+inuc(1))*Q(1)**2/A(1) )*0.5_dp &
! +g*A(i)*(Y(1+1)-Y(1))*0.5_dp &
!)
!
!!Implicit friction term
!Qpred2= (g*A*(-sign(1._dp, Qpred1)/(Apred**2))*rmu )!*efd2 !*zer 
!
!!Calculate Qpred
!do i=1,l
!if(Qpred2(i).ne.0._dp) THEN
!Qpred(i)= (1._dp - sqrt(1._dp- 4._dp*delT*Qpred2(i)*(Qpred1(i)+visc(i)) ))/(2._dp*delT*Qpred2(i))
!else
!Qpred(i)=0._dp
!end if
!end do
!
!!print*, "yes"
!!!!!!!!!!!Boundary condition
!       ! Ypred(l)=Y(l) 
!        Ypred(1)= d0!+ (g*(Y(1)-h(1)))**.5 -(g*(Y(2)-h(1)))**.5 ! -(d0-Ypred(2))!Y(1)- .5*(Y(2)-Y(1)) ! This is Y(0.5)
!
!A=Apred
!Y=Ypred
!Q=Qpred
!
!end subroutine hy23
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine fctlim(Y,Q,Ylast,Qlast, delT,delX,l)
!!This hasn't worked well and should be debugged
!
!Integer, intent(in)::l
!REAL(dp), INTENT(in):: Ylast(l), Qlast(l), delT,delX
!REAL(dp),intent(in out):: Y(l), Q(l)
!
!REAL(dp):: D(l,2), AA(l,2 ), LL(l, 2), eta(l,2), B(l,2), one=1.00
!Integer:: i
!
!D(2:l-1, 1)= (Ylast(1:l-2) -2*Ylast(2:l-1)+Ylast(3:l))*.5 
!AA(2:l-1,1)= (-Y(2:l-1)+Y(3:l))*.5
!D(2:l-1, 2)= (Qlast(1:l-2)-2*Qlast(2:l-1)+ Qlast(3:l))*.5 
!AA(2:l-1,2)= (-Q(2:l-1)+Q(3:l))*.5
!
!!print*, D
!
!Y(2:l-1)= Y(2:l-1)+D(2:l-1,1)
!Q(2:l-1)= Q(2:l-1)+D(2:l-1,2)
!
!DO i =2, l-1
!B(i,1)=sign(one,AA(i,1))
!B(i,2)=sign(one,AA(i,2))
!
!IF(Q(i)>=0) THEN
!        eta(i,1)= abs(AA(i,1))
!        eta(i,2)=abs(AA(i,2))
!ELSE
!        eta(i,1)=abs(AA(i-1,1))
!        eta(i,2)=abs(AA(i-1,2))
!END IF
!END DO
!!print*, B(2:l-1,1)
!!print*, B(2:l-1,2)
!
!LL(2:l-2,1)=B(2:l-2,1)*max(0., min( B(2:l-2,1)*(Y(2:l-2)-Y(1:l-3)),B(2:l-2,1)*(Y(4:l)-Y(3:l-1)),eta(2:l-2,1)))!abs(AA(2:l-2,1))))
!LL(2:l-2,2)=B(2:l-2,2)*max(0., min( B(2:l-2,2)*(Q(2:l-2)-Q(1:l-3)),B(2:l-2,2)*(Q(4:l)-Q(3:l-1)),eta(2:l-2,2)))!abs(AA(2:l-2,2))))
!
!!print*,LL(2:l-2,1)
!
!Y(2:l-2)= Y(2:l-2)- LL(2:l-2,1)+LL(1:l-2,1)
!Q(2:l-2)= Q(2:l-2)- LL(2:l-2,2)+LL(1:l-2,2)
!
!end subroutine fctlim
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine limi(l,delT, delX, Areafil, Qfil, Q, A, B, h,Y)
!
!!!DEFINE VARIABLES AS NEEDED
!INTEGER, INTENT(In):: l
!
!REAL(dp), INTENT(In out):: Areafil, Qfil
!REAL(dp), INTENT(In):: Q, A, B, h,delX,delT,Y
!DIMENSION::Areafil(l), Qfil(l), Q(l), A(l), B(l), h(l),Y(l)
!
!REAL(dp):: g=9.8,one=1.00, v, r1, r2,nnn
!
!
!INTEGER::i,s1,s2
!REAL(dp)::lambda1(l),lambda2(l), lam1half(l), lam2half(l), ent1half(l), ent2half(l),alpha1half(l)
!REAL(dp):: alpha2half(l), lim1(l), lim2(l), fun1half(l), fun2half(l)
!REAL(dp):: uh(l-1), ch(l-1), I1(l)
!
!
!
!	!!!!!!!!!!!!!!!!!!!!!Now we need to set up the limiter. Note that it can be based on either the corrector values, or the values from the last time step. ('comparison TVD Schemes paper') 
!	!!!!Hey, might be useful to change the code so the the  limiter is only invoked if it is needed, based on some criterion of the smoothness of the water surface. 
!	!!First calculate eigenvalues
!
!	lambda1= Q/A+(g*A/B)**0.5
!	lambda2= Q/A-(g*A/B)**0.5
!
!	!!!Now the forward half velocity and celerity, averaged as in Tseng. 	
!	uh=(A(1:(l-1))**0.5*Q(1:(l-1))/A(1:(l-1)) + A(2:l)**0.5*Q(2:l)/A(2:l) )/( A(1:(l-1))**0.5 +A(2:l)**0.5 ) 
!
!
!	!Forward half celerity 'ch' needs a loop to get it. This is using the approach of Garcia Navarro et al., 1992. It causes this to fall over commonly.  This did not happen with the simple (g*h)^0.5 method, however I still had conservation problems galore with that. 
!
!	I1= A**2/(2*B)
!	DO i=1,(l-1)
!		IF (  abs(A(i+1)-A(i))<0.00005 .OR. sign(one, (I1(i+1)-I1(i)) )/=sign(one,(A(i+1)-A(i)))  ) THEN
!		ch(i)= 0.5*( (g*A(i+1)/B(i+1))**0.5+(g*A(i)/B(i))**0.5)
!
!        ELSE
!		ch(i)=(g*(I1(i)+ I1(i+1)  )/( A(i+1)+A(i) ) )**0.5 !!Note the trial alteration. This seems sensible to me (sort of a weighted average approximation of the height) and it seems to work, whereas the other approach causes errors and doesn't seem correct to me. 
!		END IF
!	END DO
! 
!	!write(7,*) ch
! 
!	!!!Now define them at the midpoints. Might want a better choice of end point.
!	lam1half(1:(l-1))= uh+ch
!	lam1half(l)= 1.5*lambda1(l)-0.5*lambda1(l-1)
!	lam2half(1:(l-1))= uh-ch
!	lam2half(l)= 1.5*lambda2(l)-0.5*lambda2(l-1)
!
!	!!Now define the entropy fix. Note the final variable in the max-- why
!        !do this?? When blanked out I don't get some of the problems. 
!	DO i=1,(l-1)
!		ent1half(i)=max(0.,lam1half(i)-lambda1(i), lambda1(i+1)-lam1half(i)) !abs(lam1half(i)) )
!		ent2half(i)=max(0.,lam2half(i)-lambda2(i), lambda2(i+1)-lam2half(i)) !abs(lam2half(i)) )
!	END DO
!	ent1half(l)= max(0.,lam1half(l)-lambda1(l))!abs(lam1half(l)) )
!	ent2half(l)= max(0.,lam2half(l)-lambda2(l))!abs(lam2half(l)) )
!
!	!!now define the alphas. Note that the definition of alpha is the key aspect of Tseng's treatment. Here we use a similar approach.  
!
!	DO i=1,(l-1)
!		v=0.5* (0.5*(g*A(i)/B(i))**0.5 +0.5*(g*A(i+1)/B(i+1))**0.5)**(-1)
!                !nnn=((Q(i)/A(i))*(B(i)*h(i))-(Q(i+1)/A(i+1))*(B(i+1)*h(i+1)))
!                
!alpha1half(i)=v*(-lam2half(i)*((B(i+1)*Y(i+1))-(B(i)*Y(i))) +(Q(i+1)-Q(i))) 
!alpha2half(i)=v*(lam1half(i)*((B(i+1)*Y(i+1))-(B(i)*Y(i)))  -(Q(i+1)-Q(i)))
!
!                !alpha1half(i)= v*(-lam2half(i)*((A(i+1)+h(i+1)*B(i+1))-(A(i)+h(i)*B(i)) ) +(Q(i+1)-Q(i) ) )  !  /(0.5*(B(i+1)+B(i)))  ) 
!		!alpha2half(i)=  v*(lam1half(i)*((A(i+1)+h(i+1)*B(i+1))-(A(i)+h(i)*B(i)) ) - (Q(i+1)-Q(i)) )                 !/(0.5*(B(i+1)+B(i)))  ) 
!
!
!	END DO
!		alpha1half(l)= 2*alpha1half(l-1)-alpha2half(l-2)
!		alpha2half(l)=2*alpha2half(l-1)-alpha2half(l-2)
!	
!	
!
!	!!!!!!!!!!Now we define the limiter
!
!	DO i=2,(l-1)
!	s1= sign(one, lam1half(i))
!	s2= sign(one, lam2half(i))
!	
!
!	r1= alpha1half(i-s1)/alpha1half(i)
!	r2= alpha2half(i-s2)/alpha2half(i)
!
!
!	!!!Minmod limiter	
!	IF (r1>0.) then
!		lim1(i)= min(abs(r1),one) 
!		ELSE 
!		lim1(i)=0.
!		END IF
!	IF (r2>0.) then
!		lim2(i)= min(abs(r2),one) 
!		ELSE 
!		lim2(i)=0.
!		END IF
!
!
!	!Superbee limiter.
!
!	lim1(i)=max(0.,min(1.,2.*r1),min(2.,r1))
!	lim2(i)=max(0.,min(1.,2.*r2),min(2.,r2)) 
!
!	END DO
!	
!	!!!For now we just make these things zero at the end points. Hopefully the boundary conditions will make everything work
!	lim1(1)=2*lim1(2)-lim1(3)
!	lim2(1)=2*lim2(2)-lim2(3)
!	lim1(l)=2*lim1(l-1)-lim1(l-1)
!	lim2(l)=2*lim2(l-1)-lim2(l-1)
!
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NOW lets get closer.
!
!	fun1half= ent1half*(1-(delT/delX)*abs(lam1half))*(1-lim1)*alpha1half
!	fun2half= ent2half*(1-(delT/delX)*abs(lam2half))*(1-lim2)*alpha2half
!	
!	!!!Filter on the areas
!	Areafil(2:l)= 0.5*(delT/delX)*(fun1half(2:l)+fun2half(2:l) -fun1half(1:(l-1))-fun2half(1:(l-1)) )
!	Areafil(1)= 2*Areafil(2)-Areafil(3)
!
!	!!!Filter on the discharges. Need to break it up to deal with that line length issue. 
!	Qfil(2:l)=lam1half(2:l)*fun1half(2:l)+lam2half(2:l)*fun2half(2:l) 
!	Qfil(2:l)= Qfil(2:l)-lam1half(1:(l-1))*fun1half(1:(l-1))-lam2half(1:(l-1))*fun2half(1:(l-1))
!	Qfil(2:l)= 0.5*(delT/delX)*Qfil(2:l)
!	Qfil(1)=2*Qfil(2)-Qfil(3)
!
!
!End Subroutine limi
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine methchar(Q,A,B,delT,delX,delHdelX,C0,l,Alast,Qlast)
!        !!!!!!!!!!!!!!!Now let's try Garcia Navarro's (1992, first paper)method
!        !of characteristics (though including the unsteady terms). 
!       
!
!Integer, Intent (In):: l        
!REAL(dp), Intent(In):: delX, C0,delT,delHdelX,B, Alast,Qlast
!REAL(dp), Intent(In out):: A,Q
!DIMENSION:: B(l),A(l),Q(l), C0(l),delHdelX(l), Alast(l), Qlast(l) 
!
!
!LOGICAL:: conv
!REAL(dp):: g=9.8
!Integer:: aaa
!REAL(dp)::xfirstplus,xfirstpluslast,Qml,Qm,aa,hfirstplus,hfirstplusb,Qfirstplus,Qfirstplusb
!REAL(dp)::bfirstplus,sfirstplus,bb,Kb1p1,Kb1p2,KB1p3,KB1p4,KA1, Afirstplus
!REAL(dp):: xlastminus,hlastminus,hlastminusb,Qlastminus,Qlastminusb
!REAL(dp)::blastminus,slastminus, cc, KA100, Kb100p1,Kb100p2,KB100p3,KB100p4, Aml
!
!
!	
!	xfirstplus= 1- (Qlast(1)/Alast(1)- (g*Alast(1)/B(1))**0.5 )*delT/delX
!        
!        conv=.false.  ! a flag to test for convergence. 
!        Qml=Q(1)
!        
!aaa=0
!        DO While (.NOT. conv)
!        aaa=aaa+1
!        IF (aaa.ge.2000) STOP
!        !!Add a check in here
!
!
!	!!!define the constants -- Hey, this shouldn't necessarily be based
!        !around the point 2! -- I wonder if there are errors more generally in
!        !this?? Have I properly accounted for time?? Is my modification crazy!! 
!	
!	hfirstplus= Alast(2)/B(2)- (Alast(2)/B(2)- Alast(1)/B(1))*(2-xfirstplus)
!	!hfirstplusb=Alast(2)/B(2) -(Alast(2)/B(2)- Alast(1)/B(1))*(2-xfirstplus)  !for time derivative
!	Qfirstplus= Qlast(2)-(Qlast(2)-Qlast(1))*(2-xfirstplus)
!        Afirstplus= Alast(2)-(Alast(2)-Alast(1))*(2-xfirstplus)
!	!Qfirstplusb= Qlast(2)-(Qlast(2)-Qlast(1))*(2-xfirstplus) !for time derivative
!	bfirstplus=B(2)-(B(2)-B(1))*(2-xfirstplus)
!	sfirstplus=delHdelX(2)-(delHdelX(2)-delHdelX(1))*(2-xfirstplus) !(xfirstplus-1)*(delHdelX(2))+(2-xfirstplus)*delHdelX(1) 
!	
!        !print*, hfirstplus, Qfirstplus, bfirstplus, sfirstplus
!        !print*, xfirstplus
!	!!One change that might be useful is treating the bank/bed derivatives more carefully with some weighted averaging. 
!	!!Bunch of important constants. The constants with 1 in the label are for the mouth, and those with 100 in the label are for the landward end. 
!
!
!	!!!!So here we are using the approach of Garcia Navarro and Sav.. (1992), but I've used the full characteristic equation in getting the difference eq, whereas they ignored temporal derivative terms.  
!	!the mouth first. Note at present the temporal derivative ignores any morphodynamic changes, probably fine in most cases. 	
!	
!	
!	KB1p2= Qfirstplus/(Afirstplus) !useful shortcuts		
!	bb=delX*(xfirstplus-1) !spatial step
!	!!Now we define useful constants
!	KA1= (KB1p2 + (g*hfirstplus)**0.5)*bfirstplus
!        !KA1=0
!        KB1p1= Qfirstplus-KA1*hfirstplus
!        
!	KB1p3=delT*g*Afirstplus*( sfirstplus-(KB1p2)*abs(KB1p2)/(C0(1)**2*9.8*(hfirstplus)) )
!	
!        !print*,KB1p2, KA1, KB1p1, KB1p3
!
!	aa= KA1*A(1)/B(1)+KB1p1+KB1p3 !2*(g*hfirstplus)**0.5*bfirstplus  !!This term is the characteristic form of discharge at the mouth
!	
!	Qm= 0.5*(aa +Qml) !Average of this estimate and the last one. 
!	!print*, Qm/A(1), Qml/A(1)
!	!!Iterate here to improve convergence
!	
!        IF (abs(Qm-Qml)>0.00001) THEN   
!
!	xfirstplus= 1- 0.5*( (Qm/A(1)-(g*A(1)/B(1))**0.5) + (KB1p2-(g*hfirstplus)**0.5) )*delT/(delX)
!        ! So conv is still false, so we keep iterating.
!        
!        IF (xfirstplus>2) STOP !!Need to fix this up!
!
!        Qml= Qm
!
!        ELSE
!        conv=.true.
!        Q(1)=Qm
!        aaa=0
!        END IF
!        END DO !end of the do while statement. 
!	Q(1)=Qm
!	
!
!	!Now we treat the other boundary
!        
!
!	xlastminus=l-( Qlast(l)/Alast(l)+(g*Alast(l)/B(l))**0.5 )*delT/delX
!        conv=.false.
!        
!
!        IF((xlastminus<l-1).OR.(xlastminus>l)) print*, "problem with landward characteristics", l, xlastminus
!       
!
!        Aml=A(l)
!        aaa=0
!        DO WHILE (.not. conv)
!        IF (aaa.ge.2000) STOP
!        aaa=aaa+1
!	hlastminus= Alast(l-1)/B(l-1)+ (Alast(l-1)/B(l-1)- Alast(l)/B(l))*((l-1)-xlastminus)
!	!hlastminusb= Alast(l-1)/B(l-1)+ (Alast(l-1)/B(l-1)- Alast(l)/B(l))*((l-1)-xlastminus) !for time derivative
!	Qlastminus= Qlast(l-1)+(Qlast(l-1)-Qlast(l))*((l-1)-xlastminus)
!	!Qlastminusb= Qlast(l-1)+(Qlast(l-1)-Qlast(l))*((l-1)-xlastminus)
!	blastminus=B(l-1)+(B(l-1)-B(l))*((l-1)-xlastminus)	
!	slastminus=delHdelX(l-1)+(delHdelX(l-1)-delHdelX(l))*(l-1-xlastminus)
!
!	cc= delX*(l-xlastminus)
!	KB100p2=Qlastminus/(hlastminus*blastminus)
!
!	
!	KB100p1= (KB100p2 -(g*hlastminus)**0.5)*blastminus
!
!	KB100p3=delT*g*hlastminus*blastminus*(slastminus- KB100p2*abs(KB100p2)/(C0(l)**2*9.8*(hlastminus))) 
!	
!	
!        bb=((KB100p3-(Q(l)-Qlastminus))/(-KB100p1)+hlastminus) 
!!bb=bb*(1-(KB100p2+(g*hlastminus)**0.5))+(KB100p2+(g*hlastminus)**0.5)*hlastminus
!
!
!        bb=0.5*(bb*B(l)+Aml)  !Final area
!        IF (abs(bb-Aml)>0.00001) THEN
!        xlastminus= l -0.5*( Q(l)/A(l) +(g*A(l)/B(l))**0.5 + KB100p2 +(g*hlastminus)**0.5)*delT/delX
!        IF(xlastminus<l-1) THEN 
!                print*, "methcar failure"; STOP
!        END IF
!        Aml=bb
!
!ELSE
!        conv=.true.
!        A(l)=bb
!        aaa=0
!END IF
!END DO
!	
!	
!        !END DO
!	
!
!
!        
!End Subroutine methchar
!
!Subroutine reader4(a,b,l)
!
!!Character, Intent(In):: b
!Integer, Intent (In) :: l
!Real(dp), Intent(In out):: a(l), b(l)
!open(66,file= 'depths', status="old")
!read(66,*) a
!close(66)
!
!open(66,file= 'widths', status="old")
!read(66,*) b
!close(66)
!
!
!END subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine hyupdate2(delT, delX, h, B, A, Q,Q2H,Y, t,l, counter,dbdh, rmu,inuc, LAKE, hlim, Qb, tr,&
! mouth_data, mouthread, mouth_data_len, mns,rho, g,cfl, v1coef,v4coef)!slL,srL,sLlh,sLrh,C0, qbbl,export) 
!
!INTEGER, Intent (In)::l, mouth_data_len
!REAL(dp), Intent(In):: delX, Qb, tr, mouth_data, hlim, rho, mns, g, cfl, h, B
!REAL(dp),Intent(In out)::t
!REAL(dp), Intent(In out):: A,Q, Q2H, Y!, qbbl!,slL,srL,sLlh,sLrh, qbbl,C0
!DIMENSION:: h(l), B(l),A(l), Q(l), Q2H(0:l), Y(l),mouth_data(mouth_data_len,2), mns(l) !, qbbl(l)!,slL(l),srL(l),sLlh(l),sLrh(l), qbbl(l) 
!REAL(dp), Intent(In out)::delT!,export
!Integer, intent (in):: counter
!REAL(dp), Intent(in)::dbdh(l,2), rmu(l),inuc(l)
!REAL(dp), Intent(in):: v1coef, v4coef
!Logical, intent(in)::LAKE, mouthread
!
!!Local variables-- this is a mess, you really need to look at the code to see what they mean, I just wacked them in as I needed them
!!logical:: isnan
!REAL(dp)::d0
!REAL(dp):: th, dlim
!INTEGER, save:: local_counter=0
!INTEGER:: i,  j, swit,lo,ind, sde  
!REAL(dp)::Ashift(l),Qshift(l),Bbshift(l)
!REAL(dp)::Alast(l),Qlast(l), df(l), vf(l), dc(l), vc(l),   zer(l)
!REAL(dp)::Apred(l),Qpred(l),Qpred1(l),Qpred2(l),  efd(l),efd2(l)
!REAL(dp)::Apredbshift(l), Qpredbshift(l), Qpredc(l), hbshift(l)
!REAL(dp)::Acor(l),Qcor(l),Qcor1(l),Qcor2(l),Qcor3(l), Qcor4(l), Af(l),  Acorb(l)
!REAL(dp):: Areafil(l), Qfil(l),hshift(l) 
!REAL(dp)::Ylast(l),Ypred(l),Ycor(l), Yshift(l),  YY(l), QQ(l), slop(l) 
!REAL(dp)::Bnew(l), inucshift(l), s(l), up(l), low(l), diag(l),  visc(l),   Qext, w1, w3, useme, useme2(l)
!REAL(dp):: FLUXp(l,1:2), FLUXm(l,1:2), FLUXpP(l,1:2), FLUXmP(l,1:2), viscF(l), viscf4(l)
!
!LOGICAL::Limiter=.false.   !these variables control the use of the limiter and boundery treatment, and the simplified flat water treatment 
!
!!!!NOTE-- it would be really cool to see if the morphological results were
!!effected by the flood or ebb dominance. Everyone assumes that they are,
!!although in general there might not be such a major difference. 
!local_counter=local_counter+1
!IF(mod(local_counter,2).eq.0) THEN
!sde=0
!ELSE
!sde=1
!END IF   
!dlim=1._dp  !Whenever I use dlim in the code, it appears as a proportion of hlim -- so 1 means the discharge limit will be hlim
!
!!!!!!!!Note the paper Burguete et al (2007) where they provide a 'stability
!!criterion' for wetting and drying that basically says that delX< 2*depth/(f/8)
!!always. Now, in my case, with f/8 ~0.002 or something like that, this means
!!that the depth should be > ~ 5cm. I found indeed that the stability of the
!!drying front improved alot when doing this. 
!
!!Make sure that the initial state conforms to the limits on Y
!!swit=l !Predefine the wet boundary
!      
!!efd= Y-h!max(Y-h,dlim*hlim+0.*h)
!
!Alast=A
!Qlast=Q
!
!Ylast=Y
!
!IF (LAKE) THEN !!The momentum equation reduces to a flat free surface.
!        delT=10._dp 
!        goto 1209 !Straight to the water elevation and the lake hydrodynamics
!END IF
!!!here we calculate the time step 
!delT=1000._dp
!
!call time_step(delT, cfl, delX, Y-h, abs(Q/A), rmu, l )
!
!!Checks
!IF(abs(delT-1000._dp)<1.1_dp) print*, "delT= 1000."
!
!IF(delT<1.0E-04_dp) THEN 
!        print*, "short time step", delT
!        print*, "counter=", counter, "maxloc=", maxloc(abs(Q/A)), "l=", l, "maxvel=", maxval(abs(Q/A))
!        print*,  Q/A
!        !print*, maxval(h), minval(h)!maxval(efd), maxval(abs(Q)/A)
!STOP
!END IF
!
!
!1209 t=t+delT !time in seconds
!
!!!!!!!!!Flag for nan time
!        IF (isnan(t)) THEN
!        call Dat
!        print*, " Error: time is NAN", Y, h!, efd
!        Stop
!        END IF
!!!!!!!!!!!!!!!
!!TIDE --Now we find the change in water depth at the mouth. 
!
!call mouth_height(th, t,tr, mouth_data,mouthread, mouth_data_len)
!
!d0=th 
!
!IF(LAKE) THEN  !THis is supposed to calculate the flow for the theoretical case
!        !of a flat free surface
!       
!        !delT=100.
!
!        Y=max(d0+0._dp*h, h+dlim*hlim)
!        swit= l
!        DO i= 2, l
!                IF(d0<h(i)) THEN
!                        swit = i
!                        goto 22
!                END IF
!        END DO                        
!        22 Apred= (B*(Y-h))
!        
!        
!        !!Solve delQ/delX = -delA/delT with matrix inversion, upwind method. 
!        DO i = 2, swit-1
!        IF(Q(i)>0._dp) THEN 
!                diag(i)= 1._dp
!                up(i)= 0._dp
!                low(i)= -1._dp
!                s(i)=  -(delX/delT)*(Apred(i)-A(i))
!        ELSE
!                diag(i)= -1._dp
!                up(i)=1._dp
!                low(i)=0._dp
!                s(i)=  -(delX/delT)*(Apred(i)-A(i))
!        END IF
!        END DO
!        
!
!        low(1)=0._dp
!        up(1)=0._dp
!        diag(1)=1._dp
!        s(1)= sum(Apred-A)*delX/delT
!
!        call dgtsv(swit-1,1, low(1:(swit-1)), diag(1:(swit-1)), up(1:(swit-1)), s(1:(swit-1)),swit-1, j)
!
!        Q(1:swit-1)=s(1:swit-1)
!        Q(swit:l)=0._dp
!        A=Apred
!        
!        
!        goto 3541 !!go to the End of the program
!END IF
!
!!various shifted parameters that are helpful. 
!!Ashift(1:l-1)= Alast(2:l)
!!Qshift(1:l-1)= Qlast(2:l)
!!Bbshift(2:l)= B(1:l-1)
!!Yshift(1:l-1)=Ylast(2:l)
!
!
!!Qshift(l)=Qb!Q(l)!Qb!2*Qshift(l-1)  -Qshift(l-2)
!!Ashift(l)=A(l)!2*Ashift(l-1)  -Ashift(l-2)
!!Yshift(l)= 2*Y(l)-Y(l-1) !2*Yshift(l-1)-Yshift(l-2)
!!Bbshift(1)=B(1)!2*Bbshift(2)-Bbshift(3)
!!hshift(1:l-1)=h(2:l)
!!hshift(l)= 2._dp*h(l)-h(l-1)
!!hbshift(2:l)=h(2:l)
!!hbshift(1)=2._dp*h(1)-h(2) 
!!inucshift(1:l-1)=inuc(2:l)
!!inucshift(l)=inuc(l)
!
!!Now some halfway averages of the above channel variables
!!Af= (Ashift+A)*0.5_dp  !/((Yshift+Y-hshift-h)*.5)  !
!
!!!!!!
!
!
!!!!!!!!!!Predictor water level
!DO i=2,l-1
!Apred(i)= Alast(i)-(delT/delX)*(Qlast(i+1-sde)-Qlast(i-sde))
!END DO
!
!!Boundary conditions
!IF(sde==0) THEN
!Apred(1)=Alast(1)-delT/delX*(Qlast(2)-Qlast(1)) 
!Apred(l)=Alast(l)-delT/delX*(Qb-Qlast(l)) !Upstream imposed discharge
!ELSE
!Apred(1)=Alast(1)-delT/delX*(Qlast(1)-Alast(1)*(2._dp*Qlast(1)/Alast(1)-Qlast(2)/Alast(2))) !Linear extrapolation of velocity
!Apred(l)=Alast(l)-delT/delX*(Qlast(l)-Qlast(l-1))
!END IF
!
!!For this step, note that
!!Ypred=Ylast + (Apred-Alast)/(B+ 0.5_dp*db/dh*(Ypred-Ylast))
!!Rearranging things,
!!(Ypred-Ylast)*B + 0.5*db/dh*(Ypred-Ylast)^2= (Apred-Alast)
!!Ypred-Ylast= -B+-sqrt( B^2 - 2*db/dh*(Alast-Apred))/(2*0.5*db/dh)
!DO i=1,l
!        IF(Apred(i)-Alast(i)>0._dp) THEN
!        ind=1
!        ELSE
!        ind=2
!        END IF
!
!w1= B(i)**2+2._dp*dbdh(i,ind)*(Apred(i)-Alast(i))
!IF((dbdh(i,ind)>0._dp).and.(w1>0._dp)) THEN
!Ypred(i)=Ylast(i)+ (-B(i)+ sqrt(w1))/(dbdh(i,ind))
!ELSE
!	Ypred(i)= Ylast(i)+ (Apred(i)-Alast(i))/B(i)
!END IF
!END DO
!
!DO i = 1, l
!IF(isnan(Ypred(i))) THEN 
!        
!        print*, "elev", i," is Nan", Ylast !, Qshift(i), B(i), i, Q(i+1), Q(i),l !, Bshift(i)
!        print*, '....Qlast....'
!        print*, Q
!        print*, '....dbdh....'
!        print*, dbdh
!        print*, '....Alast....'
!        print*, Alast
!        print*, '....B....'
!        print*, B
!        stop
!END IF
!END DO
!
!
!!Define the slope
!DO i=2,l-1
!slop(i)= (Ylast(i+1-sde)-Ylast(i-sde))  
!END DO
!IF(sde==0) THEN
!slop(1)=Ylast(2)-Ylast(1)
!slop(l)=slop(l-1) !Extrapolate water surface
!ELSE
!slop(l)=Ylast(l)-Ylast(l-1)
!slop(1)=slop(2) !Extrapolate water surface
!END IF
!
!DO i=2,l-1
!IF((Ylast(i-sde)-h(i-sde)<hlim).or.(Ylast(i+1-sde)-h(i+1-sde)<hlim)) THEN
!slop(i)=0._dp
!END IF
!END DO
!IF(Ylast(l)-h(l)<hlim) THEN
!slop(l)=0._dp
!if(sde==0) slop(l-1)=0
!END IF
!
!!!!implicit friction
!DO i=2,l-1
!Qpred1(i)= Qlast(i)-(delT/delX)*( & 
!(1._dp+inuc(i+1-sde))*Qlast(i+1-sde)**2._dp/Alast(i+1-sde)-(1._dp+inuc(i-sde))*Qlast(i-sde)**2._dp/Alast(i-sde))-&
!(delT/delX)*g*(0.5_dp*(Alast(i+1-sde)+Alast(i-sde))*(slop(i)))  
!!Key change
!Qpred2(i)= (g*0.5_dp*(Alast(i+1-sde)+Alast(i-sde))*(-sign(1._dp, Qpred1(i))/(Apred(i)**2._dp))*rmu(i) ) 
!END DO
!
!!Boundary conditions
!IF(sde==0) THEN
!!Upstream
!Qpred1(1)= Qlast(1)-(delT/delX)*( & 
!(1._dp+inuc(2))*Qlast(1+1-sde)**2._dp/Alast(1+1-sde)-(1._dp+inuc(1-sde))*Qlast(1-sde)**2._dp/Alast(1-sde))-&
!(delT/delX)*g*(0.5_dp*(Alast(1+1-sde)+Alast(1-sde))*(slop(1)))  
!!Key change
!Qpred2(1)= (g*0.5_dp*(Alast(1+1-sde)+Alast(1-sde))*(-sign(1._dp, Qpred1(1))/(Apred(1)**2._dp))*rmu(1)) 
!
!!Downstream
!Qpred1(l)=0._dp
!Qpred2(l)=0._dp !Notice how this boundary is taken care of later
!
!ELSE
!
!!Upstream
!Qpred1(l)= Qlast(l)-(delT/delX)*( & 
!(1._dp+inuc(l+1-sde))*Qlast(l+1-sde)**2._dp/Alast(l+1-sde)-(1._dp+inuc(l-sde))*Qlast(l-sde)**2._dp/Alast(l-sde))-&
!(delT/delX)*g*(0.5_dp*(Alast(l+1-sde)+Alast(l-sde))*(slop(l)))  
!!Key change
!Qpred2(l)= (g*0.5_dp*(Alast(l+1-sde)+Alast(l-sde))*(-sign(1._dp, Qpred1(l))/(Apred(l)**2._dp))*rmu(l) )
!
!!Downstream
!Qpred1(1)= Qlast(1)-& 
!(delT/delX)*g*(0.5_dp*(Alast(2)+Alast(1))*(slop(1)))
!Qpred2(1)=(g*0.5_dp*(Alast(2)+Alast(1))*(-sign(1._dp, Qpred1(1))/(Apred(1)**2._dp))*rmu(1) ) 
!END IF
!
!
!
!!!Implicit
!!Qpred= Qpred1+visc+ delT*Qpred2*Qpred^2
!DO i=1,l
!IF(Qpred2(i).ne.0._dp) THEN
!Qpred(i)= (1._dp - sqrt(1._dp- 4._dp*delT*Qpred2(i)*(Qpred1(i)) ))/(2._dp*delT*Qpred2(i))
!ELSE
!Qpred(i)=0._dp
!END IF
!END DO
!
!!Boundary conditions 
!IF(sde==0) THEN
!Qpred(l)= Qb 
!END IF
!
!
!
!!Prevent negative areas later
!DO i=2,l-1
!        useme=max((Alast(i)-dlim*hlim*B(i)),0._dp)*delX/delT
!        IF((Qpred(i+1-(1-sde))-Qpred(i-(1-sde)))> useme ) THEN
!                Qpred(i+1-(1-sde))=0._dp 
!                Qpred(i-(1-sde))=0._dp 
!        END IF
!END DO
!IF(sde==0) THEN
!        useme=max((Alast(l)-dlim*hlim*B(l)),0._dp)*delX/delT
!        IF((Qpred(l+1-(1-sde))-Qpred(l-(1-sde)))> useme ) THEN
!                Qpred(l+1-(1-sde))=0._dp 
!                Qpred(l-(1-sde))=0._dp 
!        END IF
!ELSE
!
!END IF
!
!
!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!Now we move onto the corrector step
!!!!!!!!!!!!!!!!!!!!!
!sde=1-sde !Change the side on which we calculate derivatives
!
!!!!!!!!!!Predictor water level
!DO i=2,l-1
!Acor(i)= Alast(i)-(delT/delX)*(Qpred(i+1-sde)-Qpred(i-sde))
!END DO
!
!!Boundary conditions
!IF(sde==0) THEN
!Acor(1)=Alast(1)-delT/delX*(Qpred(2)-Qpred(1)) 
!Acor(l)=Alast(l)-delT/delX*(Qb-Qpred(l)) !Upstream imposed discharge
!ELSE
!Acor(1)=Alast(1)-delT/delX*(Qpred(1)-Apred(1)*(2._dp*Qpred(1)/Apred(1)-Qpred(2)/Apred(2))) !Linear extrapolation of velocity
!Acor(l)=Alast(l)-delT/delX*(Qpred(l)-Qpred(l-1))
!END IF
!
!!For this step, note that
!!Ycor=Ylast + (Acor-Alast)/(B+ 0.5_dp*db/dh*(Ycor-Ylast))
!!Rearranging things,
!!(Ycor-Ylast)*B + 0.5*db/dh*(Ycor-Ylast)^2= (Acor-Alast)
!!Ycor-Ylast= -B+-sqrt( B^2 - 2*db/dh*(Alast-Acor))/(2*0.5*db/dh)
!DO i=1,l
!        IF(Acor(i)-Acor(i)>0._dp) THEN
!        ind=1
!        ELSE
!        ind=2
!        END IF
!
!w1= B(i)**2+2._dp*dbdh(i,ind)*(Acor(i)-Alast(i))
!if((dbdh(i,ind)>0._dp).and.(w1>0._dp)) THEN
!        Ycor(i)=Ylast(i)+ (-B(i)+ sqrt(w1))/(dbdh(i,ind))
!else
!        Ycor(i)= Ylast(i)+ (Acor(i)-Alast(i))/B(i)
!END IF
!end do
!
!DO i = 1, l
!IF(isnan(Ycor(i))) THEN 
!        PRINT*, "elev", i," is Nan", Ylast 
!        PRINT*, '....Qlast....'
!        PRINT*, Q
!        PRINT*, '....dbdh....'
!        PRINT*, dbdh
!        PRINT*, '....Alast....'
!        PRINT*, Alast
!        PRINT*, '....B....'
!        PRINT*, B
!        STOP
!END IF
!END DO
!
!
!!Define the slope
!DO i=2,l-1
!slop(i)= (Ypred(i+1-sde)-Ypred(i-sde))  
!END DO
!IF(sde==0) THEN
!slop(1)=Ypred(2)-Ypred(1)
!slop(l)=slop(l-1) !Extrapolate water surface
!ELSE
!slop(l)=Ypred(l)-Ypred(l-1)
!slop(1)=slop(2) !Extrapolate water surface
!END IF
!DO i=2,l-1
!IF((Ypred(i-sde)-h(i-sde)<hlim).or.(Ypred(i+1-sde)-h(i+1-sde)<hlim)) THEN
!slop(i)=0._dp
!END IF
!END DO
!IF(Ypred(l)-h(l)<hlim) THEN
!slop(l)=0._dp
!if(sde==0) slop(l-1)=0._dp
!END IF
!
!
!!!!implicit friction
!do i=2,l-1
!Qcor1(i)= Qlast(i)-(delT/delX)*( & 
!(1._dp+inuc(i+1-sde))*Qpred(i+1-sde)**2._dp/Apred(i+1-sde)-(1._dp+inuc(i-sde))*Qpred(i-sde)**2._dp/Apred(i-sde))-&
!(delT/delX)*g*(0.5_dp*(Apred(i+1-sde)+Apred(i-sde))*(slop(i)))  
!!Key change
!Qcor2(i)= (g*0.5_dp*(Acor(i+1-sde)+Acor(i-sde))*(-sign(1._dp, Qcor1(i))/(Acor(i)**2._dp))*rmu(i) ) 
!end do
!
!!Boundary conditions
!if(sde==0) THEN
!!Downstream
!Qcor1(1)= Qlast(1)-(delT/delX)*( & 
!(1._dp+inuc(2))*Qpred(1+1-sde)**2._dp/Apred(1+1-sde)-(1._dp+inuc(1-sde))*Qpred(1-sde)**2._dp/Apred(1-sde))-&
!(delT/delX)*g*(0.5_dp*(Apred(1+1-sde)+Apred(1-sde))*(slop(1)))  
!!Key change
!Qcor2(1)= (g*0.5_dp*(Acor(1+1-sde)+Acor(1-sde))*(-sign(1._dp, Qcor1(1))/(Acor(1)**2._dp))*rmu(1) ) 
!
!!Upstream
!Qcor(l)=0._dp
!Qcor(l)=0._dp !Notice how this boundary is taken care of later
!
!ELSE
!
!!Upstream
!Qcor1(l)= Qlast(l)-(delT/delX)*( & 
!(1._dp+inuc(l+1-sde))*Qpred(l+1-sde)**2._dp/Apred(l+1-sde)-(1._dp+inuc(l-sde))*Qpred(l-sde)**2._dp/Apred(l-sde))-&
!(delT/delX)*g*(0.5_dp*(Apred(l+1-sde)+Apred(l-sde))*(slop(l)))  
!!Key change
!Qcor2(l)= (g*0.5_dp*(Apred(l+1-sde)+Apred(l-sde))*(-sign(1._dp, Qcor1(l))/(Acor(l)**2._dp))*rmu(l) )
!
!!Downstream
!Qcor1(1)= Qlast(1)-& 
!(delT/delX)*g*(0.5_dp*(Apred(2)+Apred(1))*(slop(1)))
!Qcor2(1)=(g*0.5_dp*(Apred(2)+Apred(1))*(-sign(1._dp, Qcor1(1))/(Acor(1)**2._dp))*rmu(1) ) 
!END IF
!
!
!
!!!Implicit
!!Qcor= Qcor1+visc+ delT*Qcor2*Qcor^2
!do i=1,l
!if(Qcor2(i).ne.0._dp) THEN
!Qcor(i)= (1._dp - sqrt(1._dp- 4._dp*delT*Qcor2(i)*(Qcor1(i)) ))/(2._dp*delT*Qcor2(i))
!else
!Qcor(i)=0._dp
!end if
!end do
!
!!Boundary conditions 
!IF(sde==0) THEN
!Qcor(l)= Qb 
!END IF
!
!
!!Limit Qcor to prevent negative depths
!!do i=1,l-1
!!useme=max(0.5_dp*(Apred(i)+Acor(i))-dlim*hlim*B(i),0._dp)*delX/delT !A useful variable
!!If(0.5_dp*(Qcor(i+1)-Qcor(i))> useme-0.5_dp*(Qpred(i+1)-Qpred(i)) ) THEN  !We need to do limiting
!!!        if((Qcor(i)>0._dp).and.(Qcor(i+1)>0._dp)) THEN
!!!        !Reduce the outflow
!!!        Qcor(i+1)=useme+Qcor(i) !Qcor(i+1)/(Alast(i))*delT/delX
!!!        else
!!!                if((Qcor(i)<0._dp).and.(Qcor(i+1)<0._dp)) THEN
!!!                !Reduce the outflow
!!!                Qcor(i)= Qcor(i+1)-useme
!!!                else 
!!!                !Limit both Q's
!!!                useme=useme/abs(Qcor(i+1)-Qcor(i))
!!!                Qcor(i+1)=Qcor(i+1)*useme
!!!                Qcor(i)=Qcor(i)*useme
!!!                end if        
!!!        
!!!        end if
!!!useme=(useme-0.5_dp*(Qpred(i+1)-Qpred(i)))/(0.5_dp*(Qcor(i+1)-Qcor(i)))
!!Qcor(i+1)=-Qpred(i+1)
!!Qcor(i)=-Qpred(i)
!!END IF
!!end do
!
!11234 t=t
!
!
!!!!!!!!!!!!!!!!!!!NOW WE MOVE ON TO THE LIMTER. NOTE THAT YOU CAN TURN IT OF WITH THE LOGICAL VARIABLE 'LIMITER'
!
!Areafil= 0._dp!A*0.
!Qfil= 0._dp!A*0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END OF THE LIMITER
!
!Q=0.5_dp*(Qpred+Qcor)
!A=.5_dp*(Apred+Acor)
!
!!!Again, back calculate Y_{i}
!do i=1,l
!        if(A(i)-Alast(i)>0._dp) THEN
!        ind=1
!        else
!        ind=2
!        end if
!w1= B(i)**2+2._dp*dbdh(i,ind)*(A(i)-Alast(i))
!if((dbdh(i,ind)>0._dp).and.(w1>0._dp)) THEN
!Y(i)=Ylast(i)+ (-B(i)+ sqrt(w1))/(dbdh(i,ind))
!else
!Y(i)= Ylast(i)+ (A(i)-Alast(i))/B(i)
!END IF
!end do
!
!!Artificial viscosity - it's no good! But try this, modified from Chaudry's
!!suggestion. 
!!do i=2,l-1
!visc(2:l-1)=abs(Y(3:l)-2._dp*Y(2:l-1)+Y(1:l-2))*(abs(Qlast(2:l-1)/Alast(2:l-1))+sqrt(g*(Ylast(2:l-1)-h(2:l-1))))*delT/delX ! & 
!!visc(2:l-1)= abs(Ylast(3:l)-2._dp*Ylast(2:l-1)+Ylast(1:l-2))&
!!*min(((abs(Qlast(2:l-1)/Alast(2:l-1))+sqrt(g*(Ylast(2:l-1)-h(2:l-1))))*delT/delX),1.0_dp) ! & 
!!+ abs( Q(3:l)/A(3:l) -2._dp*Q(2:l-1)/A(2:l-1) +Q(1:l-2)/A(1:l-2))
!visc(2:l-1)=visc(2:l-1)/(abs(Y(3:l)-h(3:l)) &
!+2._dp*abs(Y(2:l-1)-h(2:l-1) ) + abs(Y(1:l-2)-h(1:l-2)) ) ! &
!!visc(2:l-1)=visc(2:l-1)/(abs(Ylast(3:l)-h(3:l)) &
!!+2._dp*abs(Ylast(2:l-1)-h(2:l-1) ) + abs(Ylast(1:l-2)-h(1:l-2)) )! &
!!+ abs( abs(Q(3:l)/A(3:l)) +2._dp*abs(Q(2:l-1)/A(2:l-1)) +abs(Q(1:l-2)/A(1:l-2))) )
!!visc(2:l-1)= abs( (Q(3:l)/A(3:l)) - 2._dp*(Q(2:l-1)/A(2:l-1)) + (Q(1:l-2)/A(1:l-2)) )& 
!!/max(abs(Q(3:l)/A(3:l)) + 2._dp*abs(Q(2:l-1)/A(2:l-1)) + abs(Q(1:l-2)/A(1:l-2)), 1.0E-04_dp)
!!end do
!visc(1)=(abs(Y(2)-Y(1)))/(abs(Y(2)-h(2))+abs(Y(1)-h(1)))*(abs(Qlast(1)/Alast(1))+sqrt(g*(Ylast(1)-h(1))))*delT/delX
!visc(l)=(abs(Y(l)-Y(l-1)))/(abs(Y(l)-h(l))+abs(Y(l-1)-h(l-1)))*(abs(Qlast(l)/Alast(l))+sqrt(g*(Ylast(l)-h(l))))*delT/delX
!
!!visc(1)= (abs(Ylast(2)-Ylast(1)))/(abs(Ylast(2)-h(2))+abs(Ylast(1)-h(1))  )&
!!*min( ((abs(Qlast(1)/Alast(1))+sqrt(g*(Ylast(1)-h(1))))*delT/delX), 1.0_dp)
!!visc(l)= (abs(Ylast(l)-Ylast(l-1)))/(abs(Ylast(l)-h(l))+abs(Ylast(l-1)-h(l-1))  )&
!!*min(((abs(Qlast(l)/Alast(l))+sqrt(g*(Ylast(l)-h(l))))*delT/delX),1._dp)
!
!!visc(1)= abs(Q(2)/A(2)-Q(1)/A(1))/max(abs(Q(2)/A(2))+abs(Q(1)/A(1)),1.0E-04_dp)
!!visc(l)= abs( Q(l)/A(l)-Q(l-1)/A(l-1) )/max(abs(Q(l)/A(l))+abs(Q(l-1)/A(l-1)),1.0E-04_dp)
!
!
!!visc=0._dp !TURN OFF ARTIFICIAL VISC
!
!
!!Here we modify the artificial viscosity to reduce the chance of drying.
!do i=1,l-1
!viscf(i)= min(max(visc(i),visc(i+1))*v1coef, 0.5_dp)
!!viscf(i)=viscf(i)
!!viscf(i)= 0.5_dp*(visc(i)+visc(i+1))*1.0_dp
!!
!!        !Trial a fix to prevent sign changes in velocity.
!!        useme=(Q(i+1)+Q(i))*(Alast(i+1)-Alast(i))-(A(i+1)+A(i))*(Qlast(i+1)-Qlast(i))
!!        if(abs(useme)>1.0E-09_dp) THEN
!!        viscf(i)=max( min(viscf(i),1.0_dp*(-Q(i+1)*A(i)+Q(i)*A(i+1) )/useme), 0._dp)
!!        END IF
!
!viscf4(i)=max(0._dp,v4coef*1.0_dp-viscf(i))
!if(min(Y(i), Ylast(i))-h(i)<1.5_dp*hlim+0.00_dp) THEN
!viscf(i)=0._dp
!viscf(i-1)=0._dp
!end if
!end do
!
!!print*, 'v1', maxval(viscf(1:l-1))
!!Make the 'mass conservative Q2H', a discharge integrated over time (and divided
!!by the longer time step in driver2)
!!Numerically, this can be thought of as Q_{i+1/2}^{j+1/2}*delT where i is the spatial
!!increment and j is the temporal one.
!sde=1-sde !Take sde back to the value at the start of this step
!do i=1,l-1
!Q2H(i) = Q2H(i)+ 0.5_dp*( Qlast(i+1-sde) +Qpred(i+1-(1-sde)) )*delT
!end do
!IF(sde==0) THEN
!Q2H(l)=Q2H(l)+0.5_dp*(Qb+Qpred(l))*delT
!Q2H(0)= Q2H(0)+ ((delX/delT)*(A(1)-Alast(1)) + 0.5_dp*(Qlast(2)+Qpred(1)))*delT
!ELSE
!Q2H(l)=Q2H(l)+0.5_dp*(Qlast(l)+Qb)*delT
!Q2H(0)= Q2H(0)+ ((delX/delT)*(A(1)-Alast(1)) + 0.5_dp*(Qlast(1)+Qpred(2)))*delT
!END IF
!
!!print*, 'lc= ', local_counter
!
!!Correct for viscosity - this is supposed to ensure that the mass conservation still holds well. Note that this is actually the discharge integrated over the time step, not the raw discharge - so it has units m^3
!do i=1,l-1
!Q2H(i)=Q2H(i) -delX/delT*viscf(i)*(A(i+1)-A(i))*delT !The last delT is so we are integrating over the time step
!end do
!!The above formula is designed to ensure that even with viscosity, A_{i}^{j+1}= A_{i}^{j} - delt/delx*[ Q_{i+1/2}^{j+1/2}-Q_{i-1/2}^{j+1/2} ]
!!As you see will see below, viscosity is implemented as A_{i}^{j+1} <-- A_{i}^{j+1} + visc(i+1/2)*(A_{i+1}^{j+1}-A_{i}^{j+1}) + visc(i-1/2)*(A_{i}^{j+1}-A_{i-1}^{j+1})
!!The correction in the above loop ensures that even for this updated A, the mass conservation with Q2H is satisfied.
!
!!Apply viscosity to Y and V- these are the variables for which large gradients are not reasonable
!!Here we use Qpred2 for Y_old and Y for Ynew. We use Qpred1 for (Q/A)_old
!Qpred2=A
!Qpred1=Q
!do i=2,l-1
!A(i)=A(i)+viscf(i)*(Qpred2(i+1)-Qpred2(i)) -viscf(i-1)*(Qpred2(i)-Qpred2(i-1))
!!
!Q(i)=(Q(i)+viscf(i)*(Qpred1(i+1)-Qpred1(i)) -viscf(i-1)*(Qpred1(i)-Qpred1(i-1)) )
!end do
!
!A(1)=A(1)+viscf(1)*(Qpred2(1+1)-Qpred2(1)) 
!Q(1)=Q(1)+viscf(1)*(Qpred1(1+1)-Qpred1(1))
!!
!A(l)=A(l)-viscf(l-1)*(Qpred2(l)-Qpred2(l-1)) 
!Q(l)=Q(l)-viscf(l-1)*(Qpred1(l)-Qpred1(l-1))
!
!!A(1)=Qpred2(1)+viscf(1)*(Qpred2(1+1)-Qpred2(1)) 
!!Q(1)=(Qpred1(1)+viscf(1)*(Qpred1(1+1)-Qpred1(1)))
!!!
!!A(l)=Qpred2(l)-viscf(l-1)*(Qpred2(l)-Qpred2(l-1)) 
!!Q(l)=(Qpred1(l)-viscf(l-1)*(Qpred1(l)-Qpred1(l-1)))
!
!!!Again, back calculate Y_{i}
!do i=1,l
!        if(A(i)-Alast(i)>0._dp) THEN
!        ind=1
!        else
!        ind=2
!        end if
!w1= B(i)**2+2._dp*dbdh(i,ind)*(A(i)-Alast(i))
!if((dbdh(i,ind)>0._dp).and.(w1>0._dp)) THEN
!Y(i)=Ylast(i)+ (-B(i)+ sqrt(w1))/(dbdh(i,ind))
!else
!Y(i)= Ylast(i)+ (A(i)-Alast(i))/B(i)
!END IF
!end do
!
!!!!Limit discharge
!do i=2,l-1
!if( (Q(i+1-(1-sde))-Q(i-(1-sde)))>max(A(i)-dlim*hlim*B(i), 0._dp)*delX/delT) THEN
!Q(i+1-(1-sde))=0._dp
!Q(i-(1-sde))=0._dp
!END IF
!end do
!IF((1-sde).eq.0) THEN
!        if( (Qb-Q(l))>max(A(l)-dlim*hlim*B(l), 0._dp)*delX/delT) THEN
!        Q(l)=0._dp
!        !Q(i-(1-sde))=0._dp
!        END IF
!ELSE
!        IF(Q(l)-Q(l-1)>max(A(l)-dlim*hlim*B(l), 0._dp)*delX/delT) THEN
!        Q(l)=0._dp
!        Q(l-1)=0._dp
!        END IF
!END IF
!!Y(swit)= .5*(3.*Ylast(swit)-Ylast(swit-1))
!!Q(swit)= .5*(A(swit))*( 2.*Q(swit)/A(swit) + sqrt(g/(Y(swit)-h(swit)))*(Y(swit)- (2.*Y(swit)-Y(swit-1))) )
!!8665 Q(swit)= (A(swit)-Alast(swit))/delT*delX
!!IF(Q(swit)<0.) THEN
!      !  Q(swit) = Q(swit+1) + (A(swit)-Alast(swit))/delT*delX
!       ! ELSE
!        !        Q(swit)= Q(swit-1) - (A(swit)-Alast(swit))/delT*delX
!        !END IF
!!Y(2:swit-1)= Y(2:swit-1) + delT*.075*(Y(3:swit)-2.*Y(2:swit-1)+Y(1:swit-2))
!!     Q(2:swit-1)= Q(2:swit-1) + delT*0.045*A(2:swit-1)*& 
!!     ( Q(3:swit)/A(3:swit)- 2.*Q(2:swit-1)/A(2:swit-1) +Q(1:swit-2)/A(1:swit-2))  
!
!!!THIS CONDITION IS THE FLATHER CONDITION
!!Q(1)=  A(1)*(sqrt(g/(Y(1)-h(1)))*d0 - sqrt(g/(Y(1)-h(1)))*Y(1))
!!Q(1)= !A(1)*2.*sqrt(g*A(1)/B(1))
!!Try flather boundary condition
!!Y(1)= Ylast(1) - delT*sqrt(g*(Ylast(1)-h(1)))*(Y(1)-d0)/(delX)
!!A=Bnew*(Y-h)!0.5_dp*(Apred+Acor) +Areafil
!!The velocities were quite unstable
!!Q(1)= Q(2) + delX*(A(1)-Alast(1))/delT
!!Q(1)= A(1)*(Q(1)/A(1) + sqrt(g/(Y(1)-h(1)))*(Y(1)-d0))
!!Q(1)= .5*sum(A(1:l-1)+A(2:l)-Alast(1:l-1)-Alast(2:l))*delX/delT - A(1)*(sqrt(g/(Y(1)-h(1)))*(Y(1)-d0))
!
!
!!Y(1)=d0
!!Condition suggested by Blayo and Debreu (2005). Qext=0 is my decision, since it
!!seems to do good things. - Could be interpreted as uext=0, indicating that the
!!cross sectional area at the boundary is so large that u is effectively 0.
!Qext=Qlast(1)*A(1)/Alast(1)! 0._dp!Q(1)!Q(1)/2._dp!0._dp!A(1)*(2.*Q(1)/A(1)-Q(2)/A(2))!0._dp!Qlast(1)!0._dp!.1*sum(A(2:l)+A(1:l-1) -Alast(2:l)-Alast(1:l-1))*delX/delT
!w1= Qext/A(1) +sqrt(g/(Y(1)-h(1)))*(d0-h(1))
!w3= Q(1)/A(1)-sqrt(g/(Y(1)-h(1)))*(Y(1)-h(1)) !2.*(Q(2)/A(2) -sqrt(g/(Y(1)-h(1)))*(Y(2)) ) - (Q(3)/A(3) -sqrt(g/(Y(1)-h(1)))*(Y(3))) !Note that this could be extrap in different ways.Q(1)/A(1)-sqrt(g/(Y(1)-h(1)))*Y(1)
!Q(1)= .5_dp*A(1)*(w1+w3) !.5*A(1)*(Qext/A(1) + 2.*Q(2)/A(2) -1.*Q(3)/A(3)+sqrt(g/(Y(1)-h(1)))*(d0 - 2.*Y(2)+1.*Y(3))) 
!Y(1)= sqrt((Y(1)-h(1))/g)*(-w3+w1)/2._dp +h(1)!.5*(d0+ 2.*Y(2) -1.*Y(3) - sqrt((Y(1)-h(1))/g)*(-Qext/A(1) + 2.*Q(2)/A(2) - 1.*Q(3)/A(3))) 
!
!!Y(1)=d0
!!!This is supposed to solve dc/dt - sqrt(gh)dc/dx = 0 where c is the negatively
!! or positively moving characteristic - change some of the signs to switch. 
!!!This seems to not be reflective for e crazy settings, at the price of totally fucking everything! - For this, make the characteristic in d/dx positive
!!! It is sensible when c= u-2sqrt(gh). But it is reflective and def not correct
!!Q(1)= A(1)*(Qlast(1)/Alast(1) + 2.*( sqrt(g*A(1)/Bnew(1))-sqrt(g*Alast(1)/B(1)) ) + & 
!!delT*sqrt(g*A(1)/Bnew(1))*( -2.*(sqrt(g*Alast(2)/B(2))- sqrt(g*Alast(1)/B(1)) )/delX & 
!!+(Qlast(2)/Alast(2)-Qlast(1)/Alast(1))/delX ) )
!!Y(1) = Ylast(1)- delT*(Q(2)-Q(1))/(delX*Bnew(1))
!!Q(1)= Q(2) + delX*(A(1)-Alast(1))/delT
!!Y(1)= Ylast(1)+delT*(Q(1)-( sqrt(g*(d0-h(1)))- sqrt(g*Ylast(1)-h(1)))*A(1))/(delX*Bnew(1))
!
!!This condition is non-reflective, but doesn't enforce a good mouth tide
!!Actually presently it is more reflective!
!!Y(1)= d0+ h(1)+ (2.*sqrt(g*(abs(d0-h(1))))- Q(1)/A(1))**2/(4.*g)
!
!!Q(1)=!2.*Q(2)-Q(3)
!!Suppose a travelling wave at infinity - 
!
!!Q(l)=0._dp
!!Y(1)= d0 - sqrt(g*(Y(1)-h(1)))*delT/delX*( .5*(Ylast(1)+Y(1))- .5*(Ylast(2)+Y(2))) !-delT*(.5*(Ylast(1)+Y(1)) - d0)/20._dp
!
!!if(counter>100) Y(1)= Y(1) +sqrt(g*(Y(1)-.5*(h(1)+h(2))))-sqrt(g*(Y(2)-.5*(h(2)+h(1))) )
!
!!Q(1)= Q(2)+ delX*(A(1)-Alast(1))/delT
!!Q(1)= A(1)*(Y(1)/sqrt(g*(Y(1)-h(1)))- sqrt(g*Y(1)-h(1)))
!!Q(1)= .5_dp*sum( (A(1:(l-1))+A(2:l))-(Alast(1:(l-1))+Alast(2:l)))*delX/delT + 
!!Q(1)= -sqrt(abs(Y(2)-Y(1))/delX*A(1)**2/rmu(1))*sign(1._dp, Y(2)-Y(1))
!!Q(1)= !+ (Qlast(2)-Qlast(1))!A(1)*( Qlast(2)/Alast(2) -Qlast(1)/Alast(1))
!!Q(1)= (Qlast(1)) + delT*(0.-Q(1))/delX!+ delT*.99*( (A(1)-Alast(1))/delT +(Qlast(2)-Qlast(1))/delX)/2.
!!Q(1)= A(1)*(2.*Qlast(1)/Alast(1) - 4.*sqrt(g*Alast(1)/B(1)) +4.*sqrt(g*A(1)/B(1)))/.5
!
!!usf1= Qlast(1)/Alast(1)- 2.*sqrt(g*Alast(1)/B(1))
!!usf2= Qlast(2)/Alast(2)- 2.*sqrt(g*Alast(2)/B(2))
!!Q(1)= A(1)*( usf1 -(usf2-usf1)/delX*delT+ 2.*sqrt(g*A(1)/B(1)) )
!
!
!
!!swit=l
!!usf1= Qlast(swit)/Alast(swit)+ 2.*sqrt(g*Alast(swit)/B(swit))
!!usf2= Qlast(swit-1)/Alast(swit-1)+ 2.*sqrt(g*Alast(swit-1)/B(swit-1))
!!Q(swit)= A(swit)*( usf1 -(Qlast(swit-1)/Alast(swit-1)+sqrt(g*Alast(swit-1)/B(swit-1)))*(usf1-usf2)/delX*delT- & 
!!2.*sqrt(g*A(swit)/B(swit)) )
!!Q(l)= Q(l-1)
!
!!Q(swit)= A(swit)*(Qlast(swit-1)/Alast(swit-1) + 2.*sqrt(g*Alast(swit-1)/B(swit-1))- 2.*sqrt(g*A(swit)/B(swit))+ & 
!!g*delT*Qlast(swit)**2/(Alast(swit)**2*Ylast(swit)-h(swit))*rmu(swit)*g/sqrt(g*A(swit)/B(swit)) )
!!!
!!Y(l)= Ylast(l-1)*(.5)+ Ylast(l)*(1.-.5)
!!Q(l)= A(l)*( Qlast(l-1)/Alast(l-1)*(.5)+ Qlast(l)/Alast(l)*(1.-.5))
!!Q(l)= -A(l)*(sqrt(g*(A(l)/B(l))) +Q(l-1)/A(l-1)-sqrt(g*(A(l-1)/B(l-1))))
!
!
!!YY=Y
!!QQ=Q
!
!!DO i=1, l
!!IF((Y(i)<mns(i)).or.(isnan(Y(i)) )) THEN
!!       ! print*, "Y<h in McCormack, or Y is nan", i, Y(i), Ypred(i), Ycor(i), Qpred(i), Qcor(i), Ylast(i), & 
!!       ! Qlast(i-1), Qlast(i), Qlast(i+1), B(i), Qpred1(i), Qpred2(i), h(i), Q(i)
!!!Y(i)=h(i)+hlim/1.001_dp
!!        !        print*, h
!!!        stop
!!END IF
!!END DO
!!
!!!Find the 'wet-dry' boundary
!!swit=l !Predefine
!!DO i=1, l
!!IF(Y(i)-mns(i)<dlim*hlim) THEN
!!
!!        !Q(i)=0._dp
!!
!!!       print*,i
!!swit=i
!!       goto 2109
!!END IF
!!END DO
!
!
!
!
!!2109 lo=2!swit-20!swit-20
!
!!DO i=lo, swit-2 !swit-20,swit-2
!
!
!!Q(i)= QQ(i) + .01*(delT)*( (QQ(i+1)/A(i+1) - 1.*QQ(i)/A(i))   &
!!+ (-1.*QQ(i)/A(i) + QQ(i-1)/A(i-1)) )*A(i)
!
!!END DO
!
!!Q(lo-1)= QQ(lo-1) + .01*(delT)*( & 
!!(-1.*QQ(lo-1)/A(lo-1) + QQ(lo)/A(lo)) )*A(lo-1)
!
!!Q(swit-1)= QQ(swit-1) +.01*(delT)*( & 
!!(-1.*QQ(swit-1)/A(swit-1) + QQ(swit-2)/A(swit-2)) )*A(swit-1)
!
!
!
!!update width -- actually I don't need this at all
!3541 continue !b= b!bnew!+(y-ylast)*(dbdh)
!
!
!do i=1, l
!if(isnan(Ypred(i))) print*, "Ypred ",i," is nan"
!if(isnan(Ycor(i))) print*, "Ycor ",i," is nan", dbdh(i,1:2), B(i), Bnew(i), Acor(i)-Alast(i)
!
!if(isnan(Qpred(i))) print*, "Qpred ",i," is nan"
!if(isnan(Qcor(i))) print*, "Qcor ",i," is nan", Qcor4(i), Qcor1(i), Qcor2(i), Acorb(i) & 
!, (1._dp-sqrt(1._dp-4._dp*Qcor4(i)*delT*(Qcor1(i)-Qcor2(i)+visc(i))))/(2._dp*Qcor4(i)*delT), &
!sqrt(1._dp-4._dp*Qcor4(i)*delT*(Qcor1(i)-Qcor2(i)+visc(i))), 4._dp*Qcor4(i)*delT*(Qcor1(i)-Qcor2(i)+visc(i)), &
!efd(i)!,Acorb(i)
!
!end do
!
!
!!Check that the width is not too small
!!DO i=2,l
!
!!IF (B(i)<1.) Bnew(i)= B(i)- (Y(i)-Ylast(i))*dbdh(i)  !Don't let the width fall too low-- it can lead to stability problems
!
!!END DO
!
!
!!update the average depth estimation-- do this with a weighted average
!!(weighting by width, assume the newly inundated areas have average elevation of
!!Ylast+0.5*(Y-Ylast)-- so their depth is 0.5*(Y-Ylast)
!!h= (B*h+(Bnew-B)*(Ylast+ 0.5*(Y-Ylast) ) )/Bnew
!
!!B=Bnew
!
!
!End subroutine hyupdate2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
End Module st_venant_solver
