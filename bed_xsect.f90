MODULE bed_xsect
! Module to compute rates of resuspension, deposition, and bedload transport,
! and to compute the evolution of the bed on a cross-section

! Module with global parameters
USE global_defs
USE util_various, only: minmod
IMPLICIT NONE
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_resus_bedload(a, dT, water, Q, bed,ys,Area, ff,recrd, E, C, wset, a2, tau,tau_g,vel & 
    ,counter,slopes, hlim,mor,taucrit_dep,layers, taucrit_dep_ys, dst, taucrit, rho, Qe, Qbed,qb_G, rhos,& 
    voidf, dsand, d50, g, kvis, norm, alpha, Qbedon,talmon, ysl,ysu,bedl,bedu, resus_type, bedload_type, a_ref) 
    ! Purpose: Calculate the rate of resuspension and bedload transport over a
    ! cross-section

    INTEGER, INTENT(IN)::a,a2,counter,layers
    REAL(dp), INTENT(IN)::water,Q, Area, ff, hlim,mor, dt, rho, rhos, voidf, dsand, d50, g, & 
        kvis, alpha, wset, a_ref, ys,vel,slopes, tau,tau_g, bed,taucrit_dep, taucrit_dep_ys,&
        C, taucrit, dst
    REAL(dp), INTENT(IN):: ysl,ysu,bedl,bedu 
    LOGICAL, INTENT(IN):: norm, Qbedon,talmon
    CHARACTER(char_len), INTENT(IN):: resus_type, bedload_type

    ! Output variables - 'recrd' is an optional output, 'E' is the integrated
    ! resuspension, 'Qe' is the resuspension, 'Qbed' is the bedload, and 'qb_G'
    ! is a factor in the downslope bedload formual: Qb_y = -qb_G*dh/dy
    REAL(dp), INTENT(IN OUT):: recrd, E, Qe, Qbed,qb_G

    DIMENSION bed(a),ys(a), ff(a),recrd(a),tau(a),tau_g(a), vel(a), slopes(a),taucrit_dep(a,layers),C(a),& 
              taucrit_dep_ys(a), dst(a,0:(layers+1)), taucrit(a,0:layers),  Qe(a), Qbed(a), & 
              qb_G(0:a+1), a_ref(a) ! 

    INTEGER::i, j, bgwet, up,  jj,  info,ii, n(a)
    REAL(dp)::wslope,  Qelocal, tt, corfact, useme, dgravel
    REAL(dp):: kkkk(a), tbst(a), f(a)  
    REAL(dp)::sllength(a), slope_f(a),   Qb(a), & 
        bedlast(a), slope_minmod(a), mu_d, Qtemp, useful(a), Ceq(a) 
    REAL(dp)::writout(a2), d_star, c_a, k_scr, f_cs, si, tmp,tmp1
    !logical::writ_tau=.TRUE.  !This will write out the cross sectional taus-- will lead to massive files if you're not careful. 
    LOGICAL::  dry(a)


    ! To direct erosion normal to the bed, adjust the erosion factor accordingly.
    IF(norm) THEN
        slope_f(1:a-1) = (bed(2:a)-bed(1:a-1))/( ys(2:a)-ys(1:a-1))
        ! Here we compute the 'erosion factor' using the minmod of the left and
        ! right slopes 
        DO i=2,a-1
            slope_minmod(i)=minmod( slope_f(i), slope_f(i-1))
        END DO
            slope_minmod(1) = slopes(1)
            slope_minmod(a) = slopes(a)
            !sllength = sqrt(1._dp+ slopes(1:a)**2)
            sllength = sqrt(1._dp+ slope_minmod(1:a)**2)
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
        !dst(i, 0)=0._dp !The distance to the zeroth layer -- always 0.
        !The distance to the shear transition layers, in the vertical
        DO j= 1, layers
        !    dst(i, j)=max((bed(i)-taucrit_dep(i,j) ), 0._dp) !The distance from the bed surface to the next shear layer 
            !Check for consistency
            IF(bed(i)<taucrit_dep(i,j)) THEN 
                PRINT*, 'taucrit_dep prob again', i, j, bed(i), taucrit_dep(i,j)        
            END IF
        END DO

        !dst(i, layers+1)=9.9E+10_dp !The distance to the layer past the last layer -- a convenience 


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
        !IF(taucrit(i,jj)<0._dp) THEN
        !    PRINT*, 'taucrit <0', counter, i,j, slopes(i) 
        !    STOP
        !END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !EROSION -- There are a number of cases to consider
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(mor>0._dp) THEN

          1456  IF ((abs(tau_g(i))>taucrit(i, jj)).and.(.true.)) THEN !At least some erosion through this layer occurs
                    
                    SELECT CASE (resus_type)
                        CASE('cohesive') 
                            !! This is the rate of erosion in this layer, in
                            !meters per second of SOLID material, i.e. not
                            !accounting for porosity, which is accounted for
                            !later. 
                            Qelocal=  (alpha/rhos)*( abs(tau_g(i))-taucrit(i,jj) )**1._dp*& 
                                      (taucrit(i,jj)**(-.5_dp))*sllength(i) 
                                !min(taucrit(i,jj)**(-.5_dp), 50._dp)*sllength(i) 
                        CASE('vanrijn')
                            ! vanrijn 2007 reference concentration method
                            d_star = (d50*((rhos/rho-1._dp)*g/kvis**2)**(1._dp/3._dp))  !Van Rijn d_star parameter
                            c_a = 0.015_dp*max(dsand/d50,1.0_dp)*d50/(a_ref(i)*d_star**0.3_dp)* & 
                                    (max(0._dp,abs(tau_g(i))-taucrit(i,jj))/taucrit(i,jj))**1.5_dp ! Van Rijn reference concentration, in m^3/m^3     
                            !IF(c_a>150.0_dp) THEN
                                 
                            !END IF
                            Qelocal = wset*c_a*sllength(i) !/rhos !Rate of erosion in m/s of SOLID material
                            
                        CASE('smithmac')
                            !Based on description in Garcia and Parker
                            tmp = 2.4e-03_dp*(max(0._dp,abs(tau_g(i))-taucrit(i,jj))/taucrit(i,jj))
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
                        Qe(i)= dst(i,jj)/(mor*dT)*(1._dp-voidf) + Qelocal*(dT-tt)/dT  
                    END IF !jj<layers and qelocal        
                ELSE 
                    !So here, we (may) have eroded, and now we hit a layer that we
                    !couldn't cut through
                    Qe(i)= (dst(i,jj))/(mor*dT)*(1._dp-voidf)  
                    !But, if we were already in a higher shear layer, then we
                    !better move the points down to reflect that fact.

                END IF !tau_g>taucrit 
                    

                IF((Qe(i)<0._dp).OR.(isnan(Qe(i)))) THEN 
                    print*, "erosion < 0 or nan", taucrit(i, jj), tau_g(i), Qe(i)
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
        IF( (abs(tau_g(i))>taucrit(i,0)).AND.(Qbedon)) THEN
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
                            rho**(-.5_dp)*(abs(tau_g(i)))**(.5_dp)*&
                            ( abs(tau_g(i))-taucrit(i,0))/taucrit(i,0)*sign(1._dp, tau_g(i)) &
                            *sllength(i)                 
                CASE('mpm')
                ! Meyer-Peter and Muller
                    Qbed(i) = 8.0_dp*((abs(tau_g(i))-taucrit(i,0))/useme)**(1.5_dp)*sign(1._dp,tau_g(i))&
                                *sqrt(useme/rho*d50**2._dp)*sllength(i)
                !Qbed(i) = 8.0_dp*(abs(tau_g(i))/useme-0.047_dp)**(1.5_dp)*sign(1._dp,tau_g(i))&
                !            *sqrt(useme/rho*d50**2._dp)*sllength(i)
                CASE DEFAULT
                    PRINT*, 'ERROR: bedload_type was not correctly specified in calc_resus_bedload'
            END SELECT
        END IF

    END DO
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! COMPUTE qb_G, a downslope transport coefficient 
    ! Qby = -qb_G*d(bed)/dy
    ! qb_G is evaluated at y+1/2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    qb_G=0._dp
    
    IF(Qbedon) THEN
        ! Calculate downslope bedload transport coefficients
        tmp1 = (rho*g*(rhos/rho-1._dp)*d50) ! A constant in the lateral bedload formula
        DO i=1, a
            IF(abs(tau_g(i))>taucrit(i,0)) THEN
                !Choose downslope bedload relation
                !Note that if tau_g>taucrit, then tau>0, so we will not divide by
                !zero
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
    
    IF(.FALSE.) THEN
        IF(counter.eq.1) print*, 'WARNING: qb_G(i+1/2) set to zero if taug(i) or taug(i+1) < taucrit'
        DO i=1,a-1
           IF((tau_g(i)<=taucrit(i,0)).AND.(tau_g(i+1)<=taucrit(i+1,0))) qb_G(i)=0._dp
        END DO 
    END IF
   
 
    !!Total erosion - useful for the 1d sus sed 
    E=.5_dp*(sum( (Qe(2:a-1)*(ys(3:a)-ys(1:a-2)))) & 
        + Qe(1)*(ys(2)-ys(1) +1._dp*((water-bed(1))/(bedl-bed(1))*(ys(1)-ysl) )) & 
        + Qe(a)*(1._dp*((water-bed(a))/(bedu-bed(a))*(ysu-ys(a)) )+ys(a)-ys(a-1)) )


END SUBROUTINE calc_resus_bedload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_bed(a, dT, water, Q, bed,ys,Area, recrd, E, D,C,a2, tau,taug,&
    counter,slopes, hlim,mor,taucrit_dep,layers, taucrit_dep_ys, nos, taucrit, rho, Qe, Qbed, & 
    qb_G, wset,dqbeddx, rhos, voidf, d50, g, Qbedon, normmov,sus2d,ysl, & 
    ysu,bedl,bedu,iii, bedlast, talmon, high_order_bedload, too_steep)
    ! Purpose: Solve the sediment continuity equation to update the channel bed.

    INTEGER, INTENT(IN)::a,a2,counter,layers, nos, iii, too_steep
    REAL(dp), INTENT(IN)::water,Q, Area, hlim,mor, dt, rho, Qbed, Qe, dqbeddx, &
        rhos, voidf, d50, g, wset, ysl,ysu, bedlast, taug,C, qb_G !QbedI, dQbedI
    REAL(dp), INTENT(IN OUT):: bed, recrd, E, D,tau, ys,taucrit_dep, taucrit_dep_ys, slopes, & 
        taucrit, bedu, bedl
    LOGICAL, INTENT(IN):: Qbedon, normmov, sus2d, talmon, high_order_bedload
    DIMENSION bed(a),ys(a), recrd(0:a),tau(a),taug(a), slopes(a),taucrit_dep(nos,layers),C(a),&
        taucrit_dep_ys(nos),taucrit(nos,0:layers), Qbed(a), qb_G(0:a+1), Qe(a), dqbeddx(a), bedlast(a), too_steep(a) ! 

    INTEGER::i, j, bgwet, up, bfall, jj,dstspl, jjj, minX,maxX, storindex(a), info,ii, indd(a,layers), n(a), b(1)
    REAL(dp):: val, tmp1,dt_on_lambda 
    REAL(dp)::wslope, Qeold, tt, Qelocal, bed_tmp(0:a+1), ys_tmp(0:a+1)!, dqbeddx(a)
    REAL(dp):: kkkk(a), tbst(a), f(a), Qd(a), h_rhs(0:a+1), newxs(a), newys(a), lcslp(a)
    REAL(dp)::vinext,vi, vilast, epsl, epsu , mxslopes(a), erode, newslope(a), slim,p1, w1,w2,eqspc,dsts(a), dsts2(a), taud(a) &
        , sllength(a), hlast1,hlast2, storeys(a), C_diag(a),C_upper(a),C_lower(a), C_out(a),dDdy(a), dst(a,0:(layers+1)), taucinc,& 
        Bdist,p11,p12,taucrit_depnew(nos,layers),edvis(a),dedvisdy(a),edvisd(a),dedvisddy(a),d2hdy2(a),rnd(a),hss2_deriv(a), & 
        h_diag(0:a+1), h_upper(0:a+1), h_lower(0:a+1), h_rhs2(0:a+1), vel(a), sinf(a), bednew(a),&
        hss22(a),bedlast_tmp(a), sinsl(a), mu_d, Qtemp, useful(a), Ceq(a), h_lower2(0:a+1),h_upper2(0:a+1), &
        ecu, zetamult(a), upper(a), lower(a), diag(a), rhs(a), dyf(a), upper2(a), lower2(a), diag2(a),rhs2(a), &
        edvisy(a), dum1(a), dum2(a), dum3(a), dum4(a),spec(a),homo(a), sllength2(a)  

    REAL(dp)::writout(a2), impcon, useme
    REAL(dp)::bandmat(7,0:a+1), dbeddyH(4,0:a), AFB(7,0:a+1), RRR(0:a+1), CCC(0:a+1), XXX(0:a+1,1), rcond, ferr, berr,work(3*(a+2))
    CHARACTER(1)::EQUED
    INTEGER:: IPV(a+2), iwork(a+2)
    LOGICAL:: crazy_check=.false., dry(a), cnan


    !!!!!!!!!!!!Check if we have enough points 
    IF(a<2) THEN 
            RETURN
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate deposition rate
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SELECT CASE(sus2d)
        CASE(.FALSE.)
            ! When we do not use the fully 2D suspended load routine
            WHERE (bed<water)
                !Qd=(min(wset, water-bed)/rhos)*C
                Qd=(wset*too_steep/rhos)*C
            ELSEWHERE
                Qd = 0.0_dp
            END WHERE

        CASE(.TRUE.)
            ! Where we do use the 2D (multiple cross-sections) suspended load routine
            WHERE((abs(tau)>0._dp).and.(wset>0._dp))
                    Qd=(wset*too_steep/rhos)*C*& 
                        min( &
                        (water-bed)/(.1_dp*sqrt(abs(tau)/rho)*(water-bed)/wset*& 
                        (1._dp-exp(-wset/(.1_dp*sqrt(abs(tau)/rho)*(water-bed))*(water-bed)))) &
                        , 20._dp) 
            ELSEWHERE
                    Qd= (wset*too_steep/rhos)*C*20._dp 
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

    !Next we check if qb_G =0 everywhere. If so, then there is no point solving the
    !matrix equation.
    IF(maxval(abs(qb_G)) == 0.0_dp) THEN 
        IF(normmov.eqv..false.) THEN

            !print*, 'b ', mor*dT/(1.0_dp-voidf), maxval(Qe)
            IF(iii==1) bed = bed + (-dqbeddx + Qd - Qe)*mor*dT/(1._dp-voidf)
            IF(iii==2) bed = bedlast + (-dqbeddx + Qd - Qe)*mor*dT/(1._dp-voidf)

        ELSE !!USE POINTS WHICH MOVE IN THE DIRECTION OF THE VECTOR D-E
            print*, 'FIXME: Use of normmov is not consistent with the new definition of Qe &
                             (30/12/2010), which already includes sllength'
            stop
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
        impcon = 0.500_dp !Constant determining 'implicitness' of the bed solver. 0.5= 'Crank Nicholson', except that we lag the non-linear terms
        dt_on_lambda = mor*dT*1._dp/(1._dp-voidf) ! Time evolution constants
        
        !Calculate coefficients for bed slope terms
        IF(high_order_bedload) THEN
            call dbeddyH_approx(a,ys,bed,dbeddyH, ysl, ysu, bedl, bedu, 3)
        ELSE
            call dbeddyH_approx(a,ys,bed,dbeddyH, ysl, ysu, bedl, bedu, 2)
        END IF
        !!Set diagonals for Matrix -- termed 'h'
        !! h stores the implicit terms for [1 / (1-voidf)] dbed/dt = -d/dy (qb_G*dbed/dy) + Ds - Es
    
        DO i=1,a
            tmp1 = dt_on_lambda*impcon*(2._dp/(ys_tmp(i+1)-ys_tmp(i-1)))
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

            !print*,'   ', recrd(i), qb_G(i)
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
        !h_lower(0)=dt_on_lambda*impcon*(2._dp/(ys_tmp(2)-ys_tmp(0)))*(qb_G(0)*dbeddyH(1,0))
        tmp1=dt_on_lambda*impcon*(2._dp/(ys_tmp(2)-ys_tmp(0)))
        h_diag(0)=tmp1*(qb_G(0)*dbeddyH(2,0))
        h_upper(0)=tmp1*(qb_G(0)*dbeddyH(3,0))
        h_upper2(0)=tmp1*(qb_G(0)*dbeddyH(4,0))
        h_rhs2(0)= - (1._dp-impcon)/impcon*( h_diag(0)*bed_tmp(0)+h_upper(0)*bed_tmp(1) +h_upper2(0)*bed_tmp(2))
        !Correct h_diag
        h_diag(0)=1._dp+h_diag(0)

        tmp1=dt_on_lambda*impcon*(2._dp/(ys_tmp(a+1)-ys_tmp(a-1)))
        h_lower2(a+1)=tmp1*( - qb_G(a)*dbeddyH(1,a))
        h_lower(a+1)=tmp1*( - qb_G(a)*dbeddyH(2,a))
        h_diag(a+1)=tmp1*( - qb_G(a)*dbeddyH(3,a))
        h_upper(a+1)=0._dp !This will be zero, because dhdy(4,a) will be zero.
        h_upper2(a+1)=0._dp
        h_rhs2(a+1)=- (1._dp-impcon)/impcon*( h_diag(a+1)*bed_tmp(a+1)+h_lower(a+1)*bed_tmp(a) +h_lower2(a+1)*bed_tmp(a-1))
        !Correct h_diag
        h_diag(a+1)=1._dp+h_diag(a+1)


        IF(iii==1) THEN
            !Now define the right_hand side
            h_rhs(1:a)= (bed_tmp(1:a) + (-dqbeddx + Qd - Qe)*dt_on_lambda) + h_rhs2(1:a)
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
            h_rhs(1)= bedlast(1) + 2._dp/(ys(2)-ysl)*qb_G(0)*(-bedl/(ys(1)-ysl))*dt_on_lambda
            h_rhs(a) = bedlast(a) +2._dp/(ysu-ys(a-1))*qb_G(a)*(-bedu/(ysu-ys(a)))*dt_on_lambda
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
            print*, "bedjump in update_bed", i, bedlast_tmp(i)-bed(i), bedlast_tmp(i), bed(i), Qe(i), C(i),&
            Qd(i), dqbeddx(i)
        END IF
    END DO


END SUBROUTINE update_bed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
! Purpose: This routine takes ys, bed, and calculates the coefficients b1, b2,
! b3, b4, where dhdy_(i+1/2) = b1*bed(i-1) + b2*bed(i) + b3*bed(i+1) +
! b4*bed(i+2) These are stored in dbeddyH
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE basic_slope_limit(nos,ys,bed,failure_slope, remesh,limit_fail)
    ! Purpose: Basic routine to limit the absolute value of the lateral slope to be <= failure_slope
    ! Input: The channel geometry, and the slope at which the bank 'fails', and
    ! limit_fail = a number from (0-1], which can be used to make the failure
    ! happen over more time.
    ! Output: The updated channel geometry
    INTEGER, INTENT(IN):: nos
    REAL(dp), INTENT(IN):: ys(nos), failure_slope, limit_fail
    REAL(dp), INTENT(IN OUT):: bed(nos)
    LOGICAL, INTENT(IN):: remesh

    !Local variables
    REAL(dp):: hss(nos), tmp
    INTEGER:: i

    IF((remesh.eqv..FALSE.)) THEN
        !FIXME: At present, this is only valid(mass conserving) with a uniform mesh
        hss = bed
        
        ! Move from centre of channel to left bank
        !DO i=nos/2,2,-1 !nos/2,2,-1
        ! Move from left bank to the channel centre
        DO i=2,floor(nos*0.5_dp), 1 !nos/2,2,-1
            IF(abs(bed(i)-bed(i-1))>failure_slope*(ys(i)-ys(i-1))) THEN
                !print*, '.'
                IF(bed(i)>bed(i-1)) THEN
                    ! Note the 'explicit' nature of the computation:
                    ! bed(post_failure) = bed(pre_failure) +
                    !                     failure_amount(pre_failure)
                    tmp = (bed(i)-(bed(i-1) + failure_slope*(ys(i)-ys(i-1)) ))*limit_fail
                    hss(i) = hss(i)-tmp*0.5_dp
                    hss(i-1) = hss(i-1) + tmp*0.5_dp    
                ELSE
                    tmp = (bed(i-1)-(bed(i) + failure_slope*(ys(i)-ys(i-1))))*limit_fail
                    hss(i) = hss(i)+tmp*0.5_dp
                    hss(i-1) = hss(i-1) - tmp*0.5_dp    

                END IF
                
            END IF
        END DO
        ! Move from centre of channel to right bank
        !DO i=nos/2,nos-1,1 !nos/2,nos-1,1
        ! Move from right bank to the channel centre
        DO i=nos-1,ceiling(nos*0.5_dp), -1 !nos/2,nos-1,1
            IF(abs(bed(i)-bed(i+1))>failure_slope*(ys(i+1)-ys(i))) THEN
                IF(bed(i)>bed(i+1)) THEN
                    tmp = (bed(i)-(bed(i+1) + failure_slope*(ys(i+1)-ys(i)) ))*limit_fail
                    hss(i) = hss(i)-tmp*0.5_dp
                    hss(i+1) = hss(i+1) + tmp*0.5_dp    
                ELSE
                    tmp = (bed(i+1)-(bed(i) + failure_slope*(ys(i+1)-ys(i))))*limit_fail
                    hss(i) = hss(i)+tmp*0.5_dp
                    hss(i+1) = hss(i+1) - tmp*0.5_dp    
                END IF
            END IF
        END DO

    !Reset the bed
    bed = hss
    END IF
END SUBROUTINE basic_slope_limit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE basic_jump_limit(nos,ys,bed,bedjump, remesh,limit_fail)
    ! Purpose: Basic routine to limit the absolute value of the difference between bed points to be <=bedjump 
    ! Input: The channel geometry, and the bedjump at which the bank 'fails', and
    ! limit_fail = a number from (0-1], which can be used to make the failure
    ! happen over more time.
    ! Output: The updated channel geometry
    INTEGER, INTENT(IN):: nos
    REAL(dp), INTENT(IN):: ys(nos), bedjump, limit_fail
    REAL(dp), INTENT(IN OUT):: bed(nos)
    LOGICAL, INTENT(IN):: remesh

    !Local variables
    REAL(dp):: hss(nos), tmp
    INTEGER:: i

    IF((remesh.eqv..FALSE.)) THEN
        !FIXME: At present, this is only valid(mass conserving) with a uniform mesh
        hss = bed
        
        ! Move from centre of channel to left bank
        !DO i=nos/2,2,-1 !nos/2,2,-1
        ! Move from left bank to the channel centre
        DO i=2,floor(nos*0.5_dp), 1 !nos/2,2,-1
            IF(abs(bed(i)-bed(i-1))>bedjump) THEN
                !print*, '.'
                IF(bed(i)>bed(i-1)) THEN
                    ! Note the 'explicit' nature of the computation:
                    ! bed(post_failure) = bed(pre_failure) +
                    !                     failure_amount(pre_failure)
                    tmp = (bed(i)-(bed(i-1) + bedjump ))*limit_fail
                    hss(i) = hss(i)-tmp*0.5_dp
                    hss(i-1) = hss(i-1) + tmp*0.5_dp    
                ELSE
                    tmp = (bed(i-1)-(bed(i) + bedjump))*limit_fail
                    hss(i) = hss(i)+tmp*0.5_dp
                    hss(i-1) = hss(i-1) - tmp*0.5_dp    

                END IF
                
            END IF
        END DO
        ! Move from centre of channel to right bank
        !DO i=nos/2,nos-1,1 !nos/2,nos-1,1
        ! Move from right bank to the channel centre
        DO i=nos-1,ceiling(nos*0.5_dp), -1 !nos/2,nos-1,1
            IF(abs(bed(i)-bed(i+1))>bedjump) THEN
                IF(bed(i)>bed(i+1)) THEN
                    tmp = (bed(i)-(bed(i+1) + bedjump ))*limit_fail
                    hss(i) = hss(i)-tmp*0.5_dp
                    hss(i+1) = hss(i+1) + tmp*0.5_dp    
                ELSE
                    tmp = (bed(i+1)-(bed(i) + bedjump ))*limit_fail
                    hss(i) = hss(i)+tmp*0.5_dp
                    hss(i+1) = hss(i+1) - tmp*0.5_dp    
                END IF
            END IF
        END DO

    !Reset the bed
    bed = hss
    END IF
END SUBROUTINE basic_jump_limit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE critical_slope_wasting(dT, nos,ys,bed,failure_slope, rate)
    ! Purpose: Basic routine to limit the absolute value of the lateral slope.
    !          This happens by having downslope motion whenever the slope is <= failure_slope
    ! Input:   The channel geometry, the slope at which the bank starts to
    !          have downslope motion, and rate = constant determining rate of
    !          failure over time
    ! Output:  The updated channel geometry
    INTEGER, INTENT(IN):: nos
    REAL(dp), INTENT(IN):: ys(nos), failure_slope, rate, dT
    REAL(dp), INTENT(IN OUT):: bed(nos)

    ! Local variables
    INTEGER:: i
    REAL(dp):: flux(nos), slope, bed_pred(nos), flux_cor(nos)

    !bed =0.5_dp*( bed(1:nos) + bed(nos:1:-1))
    ! Determine rate of mass failure at i+1/2
    DO i=1,nos-1
        slope = (bed(i+1)-bed(i))/(ys(i+1)-ys(i))
        flux(i) = -rate*max( abs(slope) - failure_slope,0.0_dp)*sign(1.0_dp, slope)
    END DO

    ! This routine uses predictor-corrector type timestepping. 

    ! d(bed)/dT = -d(flux)/dy
    bed_pred(2:nos-1) = bed(2:nos-1) - dT*(flux(2:nos-1) - flux(1:nos-2))/(0.5_dp*(ys(3:nos)-ys(1:nos-2)))
    bed_pred(1) = bed(1)             - dT*(flux(1)                      )/(ys(2)-ys(1))
    bed_pred(nos) = bed(nos)         - dT*(              - flux(nos-1)  )/(ys(nos)-ys(nos-1))

    ! Re-calculate the flux using the predictor value of the bed
    DO i=1,nos-1
        slope = (bed_pred(i+1)-bed_pred(i))/(ys(i+1)-ys(i))
        flux_cor(i) = -rate*max( abs(slope) - failure_slope,0.0_dp)*sign(1.0_dp, slope)
    END DO

    ! New value of the flux = 0.5*(predictor + corrector)
    flux = 0.5_dp*(flux + flux_cor) 
   
    !print*, count(abs(flux(1:nos-1))>0.0_dp) 
    ! d(bed)/dT = -d(flux)/dy
    bed(2:nos-1) = bed(2:nos-1) - dT*(flux(2:nos-1) - flux(1:nos-2))/(0.5_dp*(ys(3:nos)-ys(1:nos-2)))
    bed(1) = bed(1)             - dT*(flux(1)                      )/(ys(2)-ys(1))
    bed(nos) = bed(nos)         - dT*(              - flux(nos-1)  )/(ys(nos)-ys(nos-1))
    
END SUBROUTINE critical_slope_wasting
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE critical_bedjump_wasting(dT, nos,ys,bed,failure_jump, rate)
    ! Purpose: Basic routine to limit the absolute value of the difference in
    !          conscutive bed elevations.  This happens by having downslope motion
    !          whenever the difference in two consecutive bed values is >= a constant
    !          value. 
    ! Input:   The channel geometry, the jump at which there starts to
    !          be downslope motion, and rate = constant determining rate of
    !          failure over time
    ! Output:  The updated channel geometry
    INTEGER, INTENT(IN):: nos
    REAL(dp), INTENT(IN):: ys(nos), failure_jump, rate, dT
    REAL(dp), INTENT(IN OUT):: bed(nos)

    ! Local variables
    INTEGER:: i
    REAL(dp):: flux(nos), jump, bed_pred(nos), flux_cor(nos)


    ! Determine rate of mass failure at i+1/2
    DO i=1,nos-1
        jump = (bed(i+1)-bed(i))
        flux(i) = -rate*max( abs(jump) - failure_jump,0.0_dp)*sign(1.0_dp, jump)
    END DO

    ! This routine uses predictor-corrector type timestepping. 

    ! d(bed)/dT = -d(flux) !/dy
    bed_pred(2:nos-1) = bed(2:nos-1) - dT*(flux(2:nos-1) - flux(1:nos-2))/(0.5_dp*(ys(3:nos)-ys(1:nos-2)))
    bed_pred(1) = bed(1)             - dT*(flux(1)                      )/(ys(2)-ys(1))
    bed_pred(nos) = bed(nos)         - dT*(              - flux(nos-1)  )/(ys(nos)-ys(nos-1))

    ! Re-calculate the flux using the predictor value of the bed
    DO i=1,nos-1
        jump = (bed_pred(i+1)-bed_pred(i))
        flux_cor(i) = -rate*max( abs(jump) - failure_jump,0.0_dp)*sign(1.0_dp, jump)
    END DO

    ! New value of the flux = 0.5*(predictor + corrector)
    flux = 0.5_dp*(flux + flux_cor) 
    
    ! d(bed)/dT = -d(flux) !/dy
    bed(2:nos-1) = bed(2:nos-1) - dT*(flux(2:nos-1) - flux(1:nos-2))/(0.5_dp*(ys(3:nos)-ys(1:nos-2)))
    bed(1) = bed(1)             - dT*(flux(1)                      )/(ys(2)-ys(1))
    bed(nos) = bed(nos)         - dT*(              - flux(nos-1)  )/(ys(nos)-ys(nos-1))
    
END SUBROUTINE critical_bedjump_wasting
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE bed_xsect
