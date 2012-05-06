MODULE util_various
! Module with random routines that don't fit so neatly in the single cross-section code

! Global parameters
USE global_defs
IMPLICIT NONE
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!
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

    INTEGER:: wetdry(nos), lenwet(nos), lastdry, width, i


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

    INTEGER:: i,j, bankl, bankr, tmp(1), n1, n2, num_pts(5), mid_chan
    REAL(dp):: slope_c(a),newxs(a), slope_f(a),&
               tmpR, res_pts(6), high_res_width, bank_old, slope_steep
    SAVE res_pts ! This will record the boundaries between zones of different resolutions
    DATA res_pts /6*0.0_dp/
    CHARACTER(char_len):: left='left', right='right'
    
    !Check the input data
    DO i= 1, a
    IF ( (isnan(ys(i))).or.(isnan(bed(i))) ) THEN
            print*, "isnan before refit "
            stop
    END IF
    END DO

    ! Bed slope, forward estimate
    slope_f(1:a-1) = (bed(2:a)-bed(1:a-1))/(ys(2:a)-ys(1:a-1))
    ! Bed slope, centred estimate
    slope_c(2:a-1) = (slope_f(2:a-1)*(ys(2:a-1) - ys(1:a-2)) + &
                      slope_f(1:a-2)*(ys(3:a) - ys(2:a-1)) )/ &
                     (ys(3:a) - ys(1:a-2))
    slope_c(1) = 0._dp
    slope_c(a) = 0._dp

    ! Find location of the maximum slope on the left half of the channel The
    ! method accepts the left most (or right most) point within a tolerence of
    ! the maximum, to avoid floating point round off chaos
    bankl=maxloc2(-slope_c, a, left, 1.0e-05)
    bankr=maxloc2(slope_c, a, right, 1.0e-05)
    !mid_chan=floor(0.5_dp*a)
    !slope_steep = minval(slope_c(2:mid_chan))
    !DO i=2,mid_chan
    !    IF(abs(slope_c(i)-slope_steep)<1.0e-05*abs(slope_steep)) THEN
    !        bankl = i
    !        continue
    !    END IF
    !END DO
    !! Find the steepest slope on the right half of the channel.
    !! The method accepts the right most point within a tolerence of the maximum,
    !! to avoid floating point round off chaos
    !slope_steep = maxval(abs(slope_c(mid_chan:a-1)))
    !DO i=a-1,mid_chan,-1
    !    IF(abs(slope_c(i)-slope_steep)<1.0e-05*abs(slope_steep)) THEN
    !        bankr = i
    !        continue
    !    END IF
    !END DO

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
   
 
    IF((res_pts(2)<ys(1)).or.(res_pts(5)>ys(a))& 
        .or.(res_pts(4)<res_pts(3)) ) THEN
        PRINT*, 'ERROR - the bank region points are not ordered correctly. &
                 Need to recode the refit routine to make this more general'
        PRINT*, res_pts
        stop
    END IF 

    n1 = min(2*floor( a*(high_res_width/(ys(a)-ys(1)))*3.0_dp), floor(0.4_dp*a)) ! Total number of points in the two high res regions
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
    print*, 'res_pts', res_pts
    print*, 'bank indices', bankl, bankr, a
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
    IF(.FALSE.)THEN
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
SUBROUTINE reader(a,l)
    ! Read a vector from a file
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

SUBROUTINE set_geo(bed, ys, waters,fs,fs_g, a_ref, a,b,hlim, read_initial_geo, read_initial_waters, water_m, water_mthick) !initial geometry
    ! A routine to set initial geometry conditions for the quasi-2D model

    INTEGER, INTENT(IN):: a, b 
    REAL(dp), INTENT(IN OUT):: bed, waters,fs,fs_g, a_ref,ys
    REAL(dp), INTENT(IN):: hlim, water_m, water_mthick
    LOGICAL, INTENT(IN):: read_initial_geo, read_initial_waters
    DIMENSION bed(a,b),waters(b),ys(a,b),fs(a,b), fs_g(a,b),a_ref(a,b)
    
    INTEGER:: i, j,m,n
    REAL(dp), ALLOCATABLE:: waters_tmp(:,:)
    CHARACTER(char_len):: waters_input_filename='watersold2'

    IF(read_initial_geo) THEN
        call readcs(bed,ys,a,b)  
    ELSE
        DO j = 1, b
            DO i= 1, a
                ys(i,j)= 1._dp*(i-1)*1000._dp/(1._dp*a-1._dp)
                bed(i,j) = 5.5_dp- 14.5*(((1._dp*(b-j))/(1._dp*b)))**1.2_dp+ &
                            (6._dp)*abs((i*1._dp-a*0.5_dp-0.5_dp)/(a*0.5_dp) )**1.2_dp  
                !4.5_dp- 8.5*(((1._dp*(b-j))/(0.6_dp*b)))**1.2_dp+ &
                !(6._dp)*abs((i/1._dp-a/2._dp-0.5_dp)/(a/2._dp) )**1.2_dp  
                !- 4.*(((1.*(b-j))/(1.*b)))**1.2  !+ (4.)*abs((i-a/2-0.5)/(a/2) )**1.2  
                ! -7.*exp(-j*1.0/(4.*b)) !initial valley waterations  !0.000002*abs((i-a/2))**2 - 5.*(100.-j)/100.
            END DO
        END DO
    END IF
   
    ! Set arbitrary fs, fs_g and a_ref values, which will be reset later 
    fs=0.032_dp
    fs_g=fs*0.1_dp
    a_ref=0.01_dp

    IF(maxval(bed)>5.0E+04_dp) THEN
        print*, 'ERROR: The maximum bed value is > 50000m. This is pretty big! &
                At present the code puts side walls of height 100000m on the &
                edge of each cross-section to prevent water flowing out. These should never &
                become inundated. As your bed is getting close to this number, I thought it &
                appropriate to throw an error in case the flow is going to get near 100000m. &
                If it is not, then change the code to stop this error -- &
                otherwise, increase the height of the side walls in the code &
                below'
                stop
    END IF
    bed(1,:)= 1.0E+05_dp  !Side wall, to prevent water flowing outside of the cross-section
    bed(a,:)= 1.0E+05_dp  !Side Wall

    ! Set the initial water elevation
    IF(read_initial_waters) THEN
        ! Read in the initial water elevation. Note that because waters is a
        ! vector, but this routine reads the data into an array, we have to use
        ! waters_tmp to read the data first
        CALL read_real_table(waters_input_filename, waters_tmp, n, 1, .false.)
        ! Check that the file was of the correct size
        IF(n.ne.b) THEN
            print*, 'ERROR: the water level initial condition file watersold2 &
                    did not have the correct number of values', b, n 
            stop
        END IF
        waters=waters_tmp(:,1)
    ELSE
        DO j= 1, b
            waters(j)= max(bed(a/2,j)+hlim/1.001_dp+water_mthick, water_m) !0.d0 !hs(a/2,j)+2.!0d0 !initial water surface wateration
        END DO 
    END IF

END SUBROUTINE set_geo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE meanvars(bed,ys,waters,fs,a,b,u,l,Width,Area,bottom, dWidth_dwaters, even,hlim) !use this to get suitably averaged variables for input into the 1D code 
    INTEGER, INTENT(IN)::a,b
    INTEGER, INTENT(IN OUT):: l,u
    LOGICAL, INTENT(IN):: even !!To enforce even cross sections. Sometimes the wetting and drying routine can fall over if we require evenness, in situations with say several pools of water 
    REAL(dp), INTENT(IN):: bed,ys,waters,hlim
    REAL(dp), INTENT(INOUT):: Width, Area, bottom, dWidth_dwaters !Averaged variables
    REAL(dp), INTENT(IN):: fs
    DIMENSION bed(a,b), ys(a,b), waters(b), fs(a,b), Width(b), Area(b), bottom(b),  l(b), u(b),dWidth_dwaters(b,2)
    INTEGER:: i,ll,uu,j, wetpts(a)

    REAL(dp):: increm(a) , edgel, edgeu, nom, db, Area_old(b)
    LOGICAL:: alldry, allwet
    ! FIXME: Notice how this assumes that there exists a 'dry' point ll-1 and
    ! uu+1. So we cannot have all points in the cross-section being wet. This
    ! makes for a slight inefficiency in prismatic 1D cases, which you may wish
    ! to remove.

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
                dWidth_dwaters(i,1:2)= -(ys(ll,i) -ys(ll-1,i))/(bed(ll,i)-bed(ll-1,i)) &
                +  (ys(uu+1,i) -ys(uu,i))/(bed(uu+1,i)-bed(uu,i))
            ELSE
                dWidth_dwaters(i,1:2)=0._dp
            END IF
           
            
        ELSE    !!Now we consider the case with mid channel bars
            IF ( (uu>ll).AND.(.not.alldry)) THEN 
                Width(i)=0._dp 
                Area(i)=0._dp
                dWidth_dwaters(i,1:2)=0._dp 

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
                                    dWidth_dwaters(i,1)= dWidth_dwaters(i,1)+ (abs((ys(j,i)-ys(j-1,i))*&
                                                (waters(i)+1.0E-06_dp-bed(j,i))/(bed(j-1,i)-bed(j,i)))-db)/1.0E-06_dp
                                        !Note that the above expression is = ( db(h+delh) - db(h))/delh
                                    dWidth_dwaters(i,2)=dWidth_dwaters(i,1)
                                END IF
                            END IF	
                        ELSE ! so wetpts(j)=0, but point j-1 might be wet and we need to account for that
                            IF (j>1) THEN
                                IF(wetpts(j-1)>0) THEN
                                    db=abs((ys(j,i)-ys(j-1,i))*(waters(i)-bed(j-1,i))/(bed(j,i)-bed(j-1,i)))
                                    Width(i)= Width(i)+ db
                                    Area(i)=Area(i)+db*(waters(i)-bed(j-1,i))*.5_dp
                                    dWidth_dwaters(i,1)= dWidth_dwaters(i,1)+ &
                    ( abs((ys(j,i)-ys(j-1,i))*(waters(i)+1.0E-06_dp-bed(j-1,i))/(bed(j,i)-bed(j-1,i)))-db)/1.0E-06_dp
                                                        dWidth_dwaters(i,2)=dWidth_dwaters(i,1)
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
                    
                    dWidth_dwaters(i,1)= -(ys(ll,i) -ys(ll-1,i))/(bed(ll,i)-bed(ll-1,i)) &
                    +  (ys(uu+1,i) -ys(uu,i))/(bed(uu+1,i)-bed(uu,i)) 
                    dWidth_dwaters(i,2)=0._dp
            END IF !END of the if statement where we deal with cross sections that either have internal dry points or are entirely dry



        END IF !End of the main calculation, i.e. end of IF(minval(wetpts(ll:uu))>0) THEN 

        !dWidth_dwaters=0._dp
         
        IF(dWidth_dwaters(i,1)<0._dp) THEN 
          PRINT*, "dWidth_dwaters<0", i, l(i),u(i), bed(l(i)-1,i),bed(l(i),i), bed(u(i)+1,i), bed(u(i),i) , ys(l(i)-1:l(i),i), &
          ys(u(i):u(i)+1,i), dWidth_dwaters(i,1) 
          STOP
        END IF

        !if(bottom(i)>waters(i)-0.03) print*, 'somethings less is going on'
        IF(Width(i)>1000._dp) THEN
            PRINT*, 'Width(', i, ') > 1000 =', Width(i), ll, uu, edgel, edgeu, alldry, allwet, l(i), u(i)
            STOP
        END IF

    END DO

    !dWidth_dwaters=0._dp
    !if(b>80) write(23,*) Area_old(80)-Area(80)

END SUBROUTINE meanvars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(dp) FUNCTION timeseries_interpolate(timeseries,timeseries_len,timeout)
    ! Purpose: Given a timeseries (array with 2 columns = (time, value)),
    !          interpolate the value at time=timeout
    INTEGER, INTENT(IN):: timeseries_len
    REAL(dp), INTENT(IN):: timeseries(timeseries_len,2)
    REAL(dp), INTENT(IN):: timeout

    REAL(dp):: tmax, dtscale
    INTEGER:: lwind
    

    tmax= timeseries(timeseries_len, 1) !The maximum time in the input data. 
    IF(timeout>tmax) THEN
        print*, 'ERROR in timeseries_interpolate: Requested to interpolate &
                beyond the range of the timeseries, ', tmax, timeout
        stop
    END IF

    dtscale= timeseries(2,1)-timeseries(1,1) !The time increment in the data
    !
    lwind= floor(timeout/dtscale)+1 !The index corresponding to tscale
    !
    timeseries_interpolate= (( timeout-timeseries(lwind,1))*timeseries(lwind+1,2) +&
                             (timeseries(lwind+1,1)-timeout)*timeseries(lwind,2))/& 
                            (timeseries(lwind+1,1)-timeseries(lwind,1)) ! Output value

END FUNCTION timeseries_interpolate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE check_for_uneven_time_increments(timeseries,timeseries_len)
    ! Check that the time spacing is constant in a timeseries array
    INTEGER, INTENT(IN):: timeseries_len
    REAL(dp), INTENT(IN):: timeseries(timeseries_len,2)

    ! Local variables
    REAL(dp):: dt(timeseries_len-1)
    INTEGER:: i

    dt=timeseries(2:timeseries_len,1) - timeseries(1:timeseries_len-1,1)

    DO i=2,timeseries_len-1
        IF(dt(i).ne.dt(i-1)) THEN
            print*, "ERROR: timeseries does not have even time increments"
            stop
        END IF
    END DO
    
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(dp) FUNCTION simpson_integrate(a,dz,f)
    ! Purpose: Numerically integrate f
    ! Assume f(1), f(2), ...f(a) are evenly spaced by dz
    INTEGER, INTENT(IN)::a
    REAL(dp), INTENT(IN):: dz, f(a)
    REAL(dp):: onethird=1.0_dp/3.0_dp

    IF(mod(a,2)==0) THEN
        print*, 'ERROR: Even number of points for simpsons rule in util_various'
        stop
    END IF

    simpson_integrate = f(1) + f(a) + 4.0_dp*sum(f(2:a-1:2)) + 2.0*sum(f(3:a-2:2))

    simpson_integrate = simpson_integrate*dz*onethird

END FUNCTION simpson_integrate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(dp) FUNCTION trapz_integrate(a,dz,f)
    ! Purpose: Numerically integrate f
    ! Assume f(1), f(2), ...f(a) are evenly spaced by dz
    INTEGER, INTENT(IN)::a
    REAL(dp), INTENT(IN):: dz, f(a)


    trapz_integrate = ((f(1) + f(a))*0.5_dp + sum(f(2:a-1)))*dz

END FUNCTION trapz_integrate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


REAL(dp) FUNCTION wedint ( ntab, h, ftab)

!*****************************************************************************80
!
!! WEDINT uses Weddle's rule to integrate data at equally spaced points.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NTAB, is the number of data points.  (NTAB-1) must be
!    divisible by 6.
!
!    Input, real ( kind = 8 ) H, is the spacing between the points at which
!    the data was evaluated.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the tabulated data values.
!
  implicit none

  integer ntab

  real(dp) ftab(ntab)
  real(dp) h
  integer i

  wedint = 0.0_dp
 
  !if ( ntab <= 1 ) then
  !  write ( *, '(a)' ) ' '
  !  write ( *, '(a)' ) 'WEDINT - Fatal error!'
  !  write ( *, '(a)' ) '  NTAB < 2'
  !  write ( *, '(a,i8)' ) '  NTAB = ', ntab
  !  stop
  !end if
 
  if ( mod ( ntab, 6 ) /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEDINT - Fatal error!'
    write ( *, '(a)' ) '  NTAB must equal 6*N+1 for some N!'
    stop
  end if
 
  do i = 1, ntab-6, 6
    wedint = wedint & 
      +           ftab(i)   &
      + 5.0_dp * ftab(i+1) &
      +           ftab(i+2) &
      + 6.0_dp * ftab(i+3) &
      +           ftab(i+4) &
      + 5.0_dp * ftab(i+5) &
      +           ftab(i+6)
  end do
 
  wedint = 3.0_dp * h * wedint / 10.0_dp
 
  return
end function wedint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(dp) FUNCTION newtcotes7( ntab, h, ftab)

!*****************************************************************************80
!
! newtcotes7 uses a seven point  newton-cotes rule to integrate data at equally
! spaced points.
!
!  Parameters:
!
!    Input, integer NTAB, is the number of data points.  (NTAB-1) must be
!    divisible by 6.
!
!    Input, real ( kind = 8 ) H, is the spacing between the points at which
!    the data was evaluated.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the tabulated data values.
!
  implicit none

  integer ntab

  real(dp) ftab(ntab)
  real(dp) h
  integer i

  newtcotes7 = 0.0_dp
 
  !if ( ntab <= 1 ) then
  !  write ( *, '(a)' ) ' '
  !  write ( *, '(a)' ) 'WEDINT - Fatal error!'
  !  write ( *, '(a)' ) '  NTAB < 2'
  !  write ( *, '(a,i8)' ) '  NTAB = ', ntab
  !  stop
  !end if
 
  if ( mod ( ntab, 6 ) /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEDINT - Fatal error!'
    write ( *, '(a)' ) '  NTAB must equal 6*N+1 for some N!'
    stop
  end if
 
  do i = 1, ntab-6, 6
    newtcotes7 = newtcotes7 & 
      + 41.0_dp  * ftab(i)   &
      + 216.0_dp * ftab(i+1) &
      + 27.0_dp  * ftab(i+2) &
      + 272.0_dp * ftab(i+3) &
      + 27.0_dp  * ftab(i+4) &
      + 216.0_dp * ftab(i+5) &
      + 41.0_dp  * ftab(i+6)
  end do
 
  newtcotes7 = h * newtcotes7 / 140.0_dp
 
  return
end function newtcotes7
!!!!!!!!!!!!!!!!!!!!
REAL(dp) FUNCTION minmod(a,b)
    REAL(dp) a, b    

    IF(sign(1.0_dp,a).ne.sign(1.0_dp,b)) THEN
        minmod=0.0_dp
    ELSE
        IF(abs(a)>abs(b)) THEN
            minmod = b
        ELSE
            minmod = a
        END IF
    END IF
        
    return
END FUNCTION minmod
!!!!!!!!!!!!!!!!!!!
REAL(dp) FUNCTION maxmod(a,b)
    REAL(dp) a, b    

    IF(sign(1.0_dp,a).ne.sign(1.0_dp,b)) THEN
        maxmod=0.0_dp
    ELSE
        IF(abs(a)>abs(b)) THEN
            maxmod = a
        ELSE
            maxmod = b
        END IF
    END IF
        
    return
END FUNCTION maxmod
!!!!!!!!!!!!!!!!!!!
REAL(dp) FUNCTION compute_area(nos, water,bed,ys,l,u)
    ! Used to compute the cross-sectional area.
    ! nos = number of points in cross-section (including dry regions)
    ! water = water elevation
    ! ys,bed = y and z coordinates of cross-section
    ! l,u = indices of lower wet point and the upper wet point on the
    ! cross-section
    INTEGER:: l,u,nos
    REAL(dp):: water, bed(nos), ys(nos)

    ! Local variables
    INTEGER:: i
    REAL(dp):: Area

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
        PRINT*, "Area is nan", l, u, water, maxval(bed(l:u)), minval(bed(l:u))
        STOP
    END IF

    !Return the value Area
    compute_area = Area    

END FUNCTION compute_area
!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(dp) FUNCTION compute_wet_width(nos, water,bed,ys,l,u)
    ! Used to compute the cross-sectional area.
    ! nos = number of points in cross-section (including dry regions)
    ! water = water elevation
    ! ys,bed = y and z coordinates of cross-section
    ! l,u = indices of lower wet point and the upper wet point on the
    ! cross-section
    INTEGER:: l,u,nos
    REAL(dp):: water, bed(nos), ys(nos)

    ! Local variables
    REAL(dp):: wet_width

    wet_width=ys(u)-ys(l) 
    IF(l>0) THEN
        IF (u<nos) wet_width = wet_width+ (water-bed(u))/(bed(u+1)-bed(u))*(ys(u+1)-ys(u))
        IF (l>1)   wet_width = wet_width+ (water-bed(l))/(bed(l-1)-bed(l))*(ys(l)-ys(l-1))
    END IF

    !Return the value Area
    compute_wet_width = wet_width    

END FUNCTION compute_wet_width
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE compute_xsect_wetted_boundary_coords(nos,ys,bed,water,l,u,&
                                                  ysl, ysu, bedl, bedu)
    ! Compute the 'exact' y,bed values at the wet/dry edge of the cross-section
    INTEGER, INTENT(IN):: nos, l, u
    REAL(dp), INTENT(IN):: ys(nos), bed(nos), water
    REAL(dp), INTENT(OUT):: ysl, ysu, bedl, bedu

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

END SUBROUTINE compute_xsect_wetted_boundary_coords
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE compute_slope(nos, slopes, bed, ys)
    ! Compute the lateral bed slope, as a distance weighted average of the right
    ! facing and left facing slopes.
    INTEGER, INTENT(IN):: nos
    REAL(dp), INTENT(IN):: bed(nos), ys(nos)
    REAL(dp), INTENT(OUT):: slopes(nos)

    REAL(dp):: slope_f(nos), dy_f(nos)

    ! Forward dy value
    dy_f(1:nos-1) = (ys(2:nos)-ys(1:nos-1))
    ! Forward estimate of the slope
    slope_f(1:nos-1) = (bed(2:nos)-bed(1:nos-1))/dy_f(1:nos-1)

    ! Central estimate of the slope
    slopes(2:nos-1)= (slope_f(2:nos-1)*dy_f(1:nos-2) + & 
                      slope_f(1:nos-2)*dy_f(2:nos-1)) &
                      /(ys(3:nos)-ys(1:nos-2))

    slopes(1)= (bed(2)-bed(1))/(ys(2)-ys(1))
    slopes(nos)= (bed(nos)-bed(nos-1))/(ys(nos)-ys(nos-1))

END SUBROUTINE compute_slope
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE compute_critical_shear(nos, layers, bed, slopes, taucrit_dep, taucrit,&
                                    dst, lincrem, mu, lifttodrag, &
                                    taucrit_slope_reduction, erconst)
    ! Compute the critical shear. This also accounts for multiple layers in the
    ! bed, and can treat the reduction of the critical shear on a slope

    INTEGER, INTENT(IN):: nos, layers
    LOGICAL, INTENT(IN):: taucrit_slope_reduction
    REAL(dp), INTENT(IN):: bed(nos), slopes(nos), lincrem, mu, lifttodrag, erconst
    REAL(dp), INTENT(IN OUT):: taucrit_dep(nos,layers), taucrit(nos, 0:layers), dst(nos,0:layers+1)

    INTEGER:: i, n, jj
    REAL(dp):: multa, aa, bb, cc

    ! taucrit_dep = the value of taucrit in sediment buried at a
    ! certain depth
    DO  i= 1, nos
        DO n= 1, layers
            !Sediments buried at a certain depth will have their critical shear increase.
            taucrit_dep(i,n)= max(min(taucrit_dep(i,n), bed(i)),bed(i)-lincrem*n)  
        END DO !N 
    END DO


    !dst = the distance of any critical shear layer to the surface.
    dst(:,0)= 0._dp
    DO i= 1, nos
        DO jj= 1, layers
            dst(i, jj)=max((bed(i)-taucrit_dep(i,jj) ), 0._dp)
        END DO
    END DO
    dst(:,layers+1)=9.9E+10_dp !The distance to the layer past the last layer -- a convenience 




    ! Update the critical shear layers to account for any slope effects
    DO i = 1, nos
        DO jj= 0, layers

            
            IF(taucrit_slope_reduction.eqv..TRUE.) THEN 
                ! Compute slope-related reduction in critical shear stress
                !IF((j.eq.1).and.(i.eq.1).and.(jj==0)) print*, 'WARNING: Critical shear on a slope is reduced'
                
                aa= mu**2*(mu*lifttodrag-1._dp)/(mu*lifttodrag+1._dp)
                bb= -2._dp*mu**2*lifttodrag*cos(atan(slopes(i)))*mu/(1._dp+mu*lifttodrag) 
                cc= mu**2*cos(atan(slopes(i)))**2 - sin(atan(slopes(i)))**2
                !multa*(critical shear on a flat bed) = critical shear on a slope.
                multa= (-bb - sqrt(bb**2-4._dp*aa*cc))/ (2._dp*aa)  
            ELSE

                multa=1._dp

            END IF

            taucrit(i,jj) = erconst*(1._dp+ jj*1._dp)*max(multa,1.0e-01_dp)

            !if(bed(i)>-1.0_dp) taucrit(i,jj) = taucrit(i,jj)*2.0_dp

            IF( isnan(taucrit(i,jj))) THEN
                PRINT*, "taucrit(", i,",", jj, ") is nan"
                STOP
            END IF

        END DO 
    END DO

END SUBROUTINE compute_critical_shear
!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE create_initial_geometry( ys, bed, taucrit_dep_ys, taucrit_dep,  &
                                  read_geo, nos, layers, Width, max_init_depth)
    ! Create the initial cross-sectional geometry,
    ! either by read_geog from some input files, or using a conveniently defined
    ! function
    INTEGER, INTENT(IN):: nos, layers
    LOGICAL, INTENT(IN):: read_geo
    REAL(dp), INTENT(IN):: Width, max_init_depth
    REAL(dp), INTENT(IN OUT):: ys(nos), bed(nos), & 
                              taucrit_dep_ys(nos),&
                              taucrit_dep(nos, layers)    

    INTEGER:: i

    !Create the channel geometry
    If(read_geo) THEN
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
            bed(i) =min( max(abs(ys(i)-Width/2._dp)-1.5_dp*max_init_depth,-max_init_depth),&
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
        
        
        
    END IF 

END SUBROUTINE create_initial_geometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_real_table(input_file, storage_array, nrows, ncols, check_regular_timeseries)
    ! Routine to read a table (with a given number of columns) into a
    ! real array 
    CHARACTER(char_len), INTENT(in):: input_file
    INTEGER, INTENT(in):: ncols
    LOGICAL, INTENT(in):: check_regular_timeseries
    INTEGER, INTENT(out)::nrows
    REAL(dp), ALLOCATABLE, INTENT(inout):: storage_array(:,:)

    INTEGER:: i, iost, offset
    REAL(dp):: tmp1
    REAL(dp), ALLOCATABLE:: tmp2(:)

    OPEN(77,file=input_file, status="old")

    ! Count the number of lines in the file, and assign to 'nrows'
    i=0 ! The counter
    DO
        READ(77,*,iostat=iost) tmp1
        !print*, 'tmp1 ', tmp1
        IF(iost<0) THEN
            ! We have reached the end of the file, record the number of rows
            nrows=i
            ! Go back to the start of the file
            REWIND(77)
            EXIT
        END IF
        i=i+1
    END DO
    
    ! Read file into temporary vector 
    print*, "READING FILE: ", input_file, '; Record length is ', i
    ALLOCATE(tmp2(nrows*ncols))
    READ(77,*) tmp2

    ! Allocate array holding boundary information
    ALLOCATE(storage_array(nrows,ncols))
    offset=0
    DO i=1,nrows
       storage_array(i,:) = tmp2(1+offset:ncols + offset) 
       offset=offset+ncols
    END DO
    !print*, storage_array(:,1), storage_array(:,2)    

    CLOSE(77)

    IF(check_regular_timeseries) THEN
        call check_for_uneven_time_increments(storage_array, nrows)
    END IF
END SUBROUTINE read_real_table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION maxloc2(quantity,a, side, tolerance)
    ! Function like 'maxloc', except that it finds the first value within
    ! 'tolerance'*maximum of the maximum, approaching from side 'side' 
    ! It is supposed to be less prone to round-off error than the standard max
    ! function
    INTEGER:: a
    REAL(dp):: quantity(a), tolerance
    CHARACTER(char_len):: side
   
    ! Local variables
    INTEGER:: i 
    REAL(dp):: quant_max

    quant_max=maxval(quantity)

    SELECT CASE (side)
        CASE('left')
            DO i=1,a
                IF(abs(quantity(i) - quant_max)<tolerance*quant_max) THEN
                    maxloc2=i
                    return
                END IF
            END DO

        CASE('right')
            DO i=a,1,-1
                IF(abs(quantity(i) - quant_max)<tolerance*quant_max) THEN
                    maxloc2=i
                    return
                END IF
            END DO

        CASE DEFAULT
            print*, "ERROR: quantity 'side' in maxloc2 must be either 'left' or&
                    'right'. Received ", side
            stop

    END SELECT
     

END FUNCTION maxloc2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE util_various
