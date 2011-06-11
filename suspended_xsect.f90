MODULE suspended_xsect

!File of global parameters
use global_defs

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dynamic_sus_dist(a, delT, ys, bed, water, waterlast, Q, tau, vel, wset, Qe,lambdacon, &
                                rho,rhos, g, d50, bedl,bedu, ysl, ysu, cb, Cbar, Qbed, &
                                sed_lag_scale, counter, high_order_Cflux, a_ref, sus_vert_prof, edify_model, &
                                x_len_scale, sconc, lat_sus_flux)
    ! Calculate the cross-sectional distribution of suspended sediment using
    ! some ideas from /home/gareth/Documents/H_drive_Gareth/Gareth_and_colab
    ! s/Thesis/Hydraulic_morpho_model/channel_cross_section/paper/idea_for_
    ! simplified_cross_sectional_sus
    
    ! We solve the equation
    !
    ! depth d (Cbar) / dt + U*depth dCbar/dx +
    ! V*depth* d(Cbar)/dy = - d/dy( Fl ) +
    ! (Es - ws*cb)
    !
    ! With suitable approximations for dCbar/dx and V, and if needed, numerical
    ! integration to calculate Fl (the lateral flux)
    ! Numerically, we use (linear) implicit discretizations, except for the
    ! dCbar/dx term (explicit)
    
    INTEGER, INTENT(IN)::a, counter
    REAL(dp), INTENT(IN):: delT, ys, bed, water, waterlast, tau, vel, wset, Qe, lambdacon, rho, rhos,g, & 
                                d50, bedl, bedu,ysl,ysu, Q, Qbed, a_ref, x_len_scale, sconc
    REAL(dp), INTENT(OUT):: sed_lag_scale, lat_sus_flux  ! sed_lag_scale = used to estimate derivative terms, e.g. dCbar/dx
    ! cb = Near bed suspended sediment concentration, 
    ! Cbar = Depth averaged suspended sediment concentration
    REAL(dp), INTENT(IN OUT):: cb, Cbar
    CHARACTER(20), INTENT(IN):: sus_vert_prof, edify_model
    LOGICAL, INTENT(IN):: high_order_Cflux
    DIMENSION ys(a), bed(a), tau(a),vel(a), Qe(a), cb(a), Cbar(a),Qbed(a), a_ref(a), lat_sus_flux(a+1)

    ! LOCAL VARIABLES
    INTEGER:: i, info, j
    LOGICAL:: halt
    REAL(dp):: depth(0:a+1), eddif_y(0:a+1), eddif_z(a), zetamult(0:a+1), vd(0:a+1), ys_temp(0:a+1)
    REAL(dp):: M1_lower(a), M1_diag(a), M1_upper(a), M1_upper2(a), M1_lower2(a)
    REAL(dp):: RHS(a), dy_all(a)
    REAL(dp):: tmp1, dy, dy_outer, tmp2, Cref, z, cbed_tmp1, cbed_tmp2
    REAL(dp):: dQdx, dhdt, Cbar_old(a), dUd_dx(0:a+1), sus_flux
    REAL(dp):: DLF(a), DF(a), DUF(a), DU2(a),rcond, ferr, berr, work(3*a), XXX(a, 1)
    REAL(dp):: bandmat(5,a), AFB(7,a), RRR(a), CCC(a), int_edif_f(a+1), int_edif_dfdy(a+1)
    INTEGER::  IPV(a), iwork(a)   
    CHARACTER(1):: EQUED

    ! This routine calculates cb, Cbar in units m^3/m^3 --- however, elsewhere
    ! they are in kg/m^3 --- so we convert here, and convert back at the end of
    ! the routine
    cb = cb/rhos
    Cbar = Cbar/rhos 

    Cbar_old = Cbar

    ! Depth, including at zero-depth boundaries ( 0 and a+1 )
    depth(1:a) = max(water -bed, 0._dp)
    depth(0) = 0._dp
    depth(a+1) = 0._dp

    ! y coordinate, including at zero depth boundaries (0 and a+1)
    ys_temp(1:a) = ys
    ys_temp(0) = ys(1) - (ys(1) - ysl)*abs(water-bed(1))/(bedl-bed(1))!ysl
    ys_temp(a+1) = ys(a) + (ysu-ys(a))*abs(water-bed(a))/(bedu-bed(a)) !ysu


    ! Calculate the value of 'zetamult', where:
    ! zetamult*cbed = Cbar = depth integrated sediment concentration
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
                    ! Ikeda and Izumi (1991) relation, divided by depth
                    zetamult(i)= eddif_z(i)/(depth(i)*wset)*(1._dp-exp(-(wset/eddif_z(i))*depth(i)) )
                ELSE 
                    zetamult(i)=1.0e-12_dp
                END IF
            END DO

        CASE('Rouse')
            ! Rouse vertical suspended sediment distribution
            DO i=1,a

                IF(abs(tau(i))>0.0_dp) THEN
                    ! Rouse number
                    z = wset/(0.4_dp*sqrt(tau(i)/rho)) 
                    ! Compute rouse integral factor using a function
                    zetamult(i) = rouse_int(z,a_ref(i)/depth(i))
                ELSE
                    zetamult(i) = 0.0_dp
                END IF

            END DO

        ! Catch errors
        CASE DEFAULT
            PRINT*,  'ERROR: sus_vert_prof does not have the correct value in dynamic_sus_dist'
            stop
        
    END SELECT

    ! Limit zetamult to a non-zero value to avoid division problems
    zetamult(0)   = 1.0_dp
    zetamult(a+1) = 1.0_dp
    zetamult = max(zetamult, 1.0e-12)


    !!! PUSH THE SUSPENDED FLUX TOWARDS THE DESIRED VALUE   
    ! Calculate total sediment flux at time = t.
    ! We will use this to 'correct' the sediment flux towards the desired value 
    sus_flux = sum( & 
            (Qbed+Cbar*abs(vel)*max(water-bed,0._dp) )*& ! Total load
            ( ( (/ ys(2:a), ysu /) - (/ ysl, ys(1:a-1) /) )*0.5_dp) &  ! dy
                  )
    ! Compute the discharge using the same approach
    tmp1 = sum( abs(vel)*max(water-bed,0.0_dp)*& 
              ( ( (/ ys(2:a), ysu /) - (/ ysl, ys(1:a-1) /) )*0.5_dp) &  ! dy
                  )

    IF(sus_flux > 1.0e-12_dp) THEN
        sed_lag_scale = 1.0_dp*((sconc*tmp1)/sus_flux) !Desired flux / actual flux

        ! Prevent very high or low values
        sed_lag_scale = min(max(sed_lag_scale,0.666_dp),1.5_dp) 
        IF(mod(counter,1000).eq.1) PRINT*, 'sed_lag_scale = ', sed_lag_scale

    ELSE
        sed_lag_scale = 1.0_dp

    END IF

    IF(mod(counter,1000).eq.1) print*, 'sus_flux is:', sus_flux, ' desired flux is:', sconc*Q

    ! Here, we try to add a constant to cb across the channel, such that the 
    ! sus_flux is as desired
    ! If we do it by adjusting cb, then the idea is we will not affect Fl
    !tmp1 = sum(1.0_dp*zetamult(1:a)*abs(vel)**max(water-bed,0._dp)*& 
    !            ( ( (/ ys(2:a), ysu /) - (/ ysl, ys(1:a-1) /) )*0.5_dp) &  ! dy
    !              )
    ! Note that if k*tmp1  = desired_extra_flux, then k*zetamult represents the
    ! extra Cbar that we need to add everywhere get the flux correct (and k is
    ! the extra cbed, which is constant, as desired). Note that the
    ! desired_extra_flux = sed_lag_scale*sus_flux - sus_flux  
    !Cbar = Cbar + ( (sed_lag_scale*sus_flux-sus_flux)/tmp1)*zetamult(1:a) !Adding a constant


    ! Now we add the term U*d*dCbar/dx using an operator splitting technique
    ! depth*dCbar/dT + depth*vel*dCbar/dx = 0.0
    ! Implicit 
    Cbar = Cbar*(1.0_dp - 0.0_dp*delT*vel*(1._dp-sed_lag_scale)/x_len_scale)/ &
           (1.0_dp + 1.0_dp*delT*vel*(1.0_dp-sed_lag_scale)/x_len_scale)


    ! Solve initially for the depth - averaged suspended sediment concentration
    ! depth d (Cbar) / dt + U*depth*dCbar/dx +
    ! V*depth* d(Cbar)/dy - d/dy( Fl ) -
    ! (Es - ws*cb) = 0.
    !
    
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
        
    ! Advection term -- dC/dx is = (C - C*sed_lag_scale)/x_length_scale

    !!dy_all = (ys_temp(2:a+1)-ys_temp(0:a-1))*0.5_dp
    !!tmp1 = sum(Cbar*vel*depth(1:a)*dy_all) !Cbar flux
    !!tmp2 =sconc*abs(Q) - sum(Qbed(1:a)*dy_all) !Desired Cbar flux = 'Measure of total load less bedload'
    !!tmp2 = max(tmp2, 0._dp)

    ! Write output to monitor convegence
    !IF(mod(counter,1000).eq.0) PRINT*, 'sus flux is =', tmp1, '; desired flux is', tmp2
    
    !DO i=1,a
    !    !Explicit
    !    Cref = Cbar(i)*sed_lag_scale 
    !    ! C_lag_scale is chosen so that at steady state, the total sediment flux equals the desired sediment flux
    !    RHS(i) = RHS(i) - vel(i)*depth(i)*(Cbar(i) - Cref)/x_len_scale
    !    !RHS(i) = RHS(i) - vel(i)*depth(i)*Cbar(i)*(1.0_dp - sed_lag_scale)*0.5_dp/x_len_scale
    !    !M1_diag(i) = M1_diag(i) + vel(i)*depth(i)*Cbar(i)*(1.0_dp-sed_lag_scale)*0.5_dp/x_len_scale
    !END DO 
             
    DO i = 1, a
        !PRINT*, 'SANITY CHECK a'
        IF(isnan(M1_diag(i))) print*, 'M1_diag(', i,') is NaN a'    
        IF(isnan(M1_lower(i))) print*, 'M1_upper(', i,') is NaN a'    
        IF(isnan(M1_upper(i))) print*, 'M1_upper(', i,') is NaN a'    
    END DO

    ! Vd Advection term -- requires calculation of vd, using an
    ! approximate model. Here we use the 2D continuity equation:
    !
    ! dwater/dt + d(U*d)/dx + d(V*d)/dy = 0
    ! 
    ! in the form:
    !
    ! V*d = int( -dwater/dt - d(U*d)/dx) dy 
    ! 
    ! with appropriate boundary conditions 
    !vd(0) = 0._dp ! =  (lateral velocity)*(depth at zero depth boundary)=0
    !vd(a+1) = 0._dp 

    vd = 0.0_dp

    ! Estimate d(ud)/dx, using
    ! d(U*d)/dx ~ Ud/Q*dQ/dx --- approximate model for longitudinal momentum derivative
    dUd_dx = 0._dp !Predefine

    ! dQ/dx = -B* dh/dt, from the cross-sectionally integrated continuity equation
    dhdt = (water - waterlast)/delT
    dQdx = - (ys_temp(a+1) - ys_temp(0))*dhdt  

    DO i = 1, a
        IF(Q.ne.0._dp) dUd_dx(i) = (vel(i)*depth(i)/Q)*dQdx  

        ! Vd(y) = int_{left_bank}^{y} [- dh/dt -dUd/dx ] dy
        vd(i) = vd(i-1) - (ys_temp(i)-ys_temp(i-1))*(dhdt +0.5_dp*(dUd_dx(i)+dUd_dx(i-1)) )  

        IF(depth(i)==0._dp) vd(i) = 0._dp

        ! Upwind discretization of the advective term -- fill out the matrices
        ! (FIXME -- consider making use of a higher order discretization)
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

        ! Treat the cases with vertically constant / variable lateral eddy
        ! diffusivity differently
        IF(edify_model=='Constant') THEN
            ! Method with vertically constant lateral eddy diffusivity
            ! This does not require numerical integration to calculate the
            ! lateral flux:
            ! Fl = -int ( eddif_y * d(cb*f)/dy ) dz
            ! -- so is much faster than the case with vertically
            ! variable lateral eddy diffusivity. 
            ! Note though that the code for the case with vertically variable
            ! lateral eddy diffusivity can also treat the case with vertically
            ! constant lateral eddy diffusivity -- but it is slower.
   
            IF(i==1) THEN 
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
            END IF 
 
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
                M1_upper(i) = M1_upper(i) - 0.5_dp*tmp1/zetamult(i+1)  ! Note that 1.0/zetamult(i)*Cbar = cb
                M1_diag(i)  = M1_diag(i)  - 0.5_dp*tmp1/zetamult(i)
            END IF

            ! Compute eddify*dbed/dy at i-1/2
            tmp1 = 0.5_dp*(eddif_y(i)+eddif_y(i-1))
            tmp1 = tmp1*( -(depth(i)-depth(i-1) )/(ys_temp(i)-ys_temp(i-1)))
            tmp1 = tmp1/dy_outer

            IF(i>1) THEN
                !Estimate of 1/dy_outer*(eddify*dbed/dy*cb) at i-1/2
                M1_diag(i)   = M1_diag(i)   + 0.5_dp*tmp1/zetamult(i)  ! Note that 1.0/zetamult(i)*Cbar = cb
                M1_lower(i)  = M1_lower(i)  + 0.5_dp*tmp1/zetamult(i-1)
            END IF
 
        ELSE
            ! Case with vertically variable eddy diffusivity. This requires
            ! numerical integration to calculate the flux terms. We calculate
            !
            ! int( Fl) dz =  - int( eddif_y* d(cb*f)/dy ) dz
            !
            ! where eddif_y is the lateral eddy diffusivity, which varies with
            ! z, and f defines the vertical profile of suspended sediment
            !
            ! Practically, the integration is based on the following method:
            !
            ! int( eddif_y*d(cb*f)/dy ) dz = 
            ! int( eddif_y*cb*df/dy + eddif_y*f*dcb/dy) dz = 
            ! cb*[ int( eddif_y * df/dy) dz] + dcb/dy*[ int(eddif_y*f) dz] = 
            ! cb*(int_edif_dfdy) + dcb/dy*(int_edif_f)

            IF(i==1) THEN
                ! Calculate important integration factors int_edif_dfdy and
                ! int_edif_f
                call int_edify_f(edify_model, sus_vert_prof, a, ys, bed, ysl, ysu,&
                                 bedl, bedu, water, sqrt(abs(tau)/rho), wset,a_ref,&
                                 int_edif_f, int_edif_dfdy,100)
                !int_edif_f=0._dp
                !int_edif_dfdy=0._dp
            END IF
            
            ! d/dy [ dcb/dy*(int(eddif_y*f) dz) ]
            tmp1 = 1.0_dp/dy_outer
            tmp2 = 1.0_dp/(ys_temp(i+1)-ys_temp(i))
            
            !IF((depth(i)<10.0_dp*depth(i+1)).and.(depth(i)<10.0_dp*depth(i-1)).and. &
            !  (depth(i+1)<10.0_dp*depth(i)).and.(depth(i-1)<10.0_dp*depth(i))) THEN
                ! When the depth changes rapidly, this term can induce negative
                ! values of Cbar. 
            !IF((depth(i+1)/depth(i)<10.0_dp).and.(depth(i)/depth(i+1)<10.0_dp)) THEN
                IF(i<a) THEN
                    ! 2 point derivative approx -- note cb = Cbar/zetamult
                    M1_upper(i) = M1_upper(i) - tmp1*tmp2*int_edif_f(i+1)/zetamult(i+1)
                    M1_diag(i)  = M1_diag(i)  + tmp1*tmp2*int_edif_f(i+1)/zetamult(i)
                END IF
            !END IF

                tmp2 = 1.0_dp/(ys_temp(i)-ys_temp(i-1))
            !IF((depth(i)/depth(i-1)<10.0_dp).and.(depth(i-1)/depth(i)<10.0_dp)) THEN
                IF(i>1) THEN
                    ! 2 point derivative approx
                    M1_diag(i)  = M1_diag(i)  + tmp1*tmp2*int_edif_f(i)/zetamult(i)
                    M1_lower(i) = M1_lower(i) - tmp1*tmp2*int_edif_f(i)/zetamult(i-1)
                END IF
            !END IF 

            
            !! d/dy ( cb*INT(edify*df/dy )dz  )
            !IF((depth(i+1)/depth(i)<10.0_dp).and.(depth(i)/depth(i+1)<10.0_dp)) THEN
                IF(i<a) THEN
                    !M1_upper(i) = M1_upper(i) - 0.5_dp*tmp1*int_edif_dfdy(i+1)/zetamult(i+1)  ! Note that 1/zetamult(i)*Cbar = cb
                    !M1_diag(i)  = M1_diag(i)  - 0.5_dp*tmp1*int_edif_dfdy(i+1)/zetamult(i)
                    
                    ! Experimenting with this approach to try to avoid negative
                    ! Cbar values
                    IF(Cbar(i)/zetamult(i) <Cbar(i+1)/zetamult(i+1)) & 
                            M1_diag(i)  = M1_diag(i)  - 1.0_dp*tmp1*int_edif_dfdy(i+1)/zetamult(i)
                    IF(Cbar(i)/zetamult(i)>=Cbar(i+1)/zetamult(i+1)) &
                            M1_upper(i)  = M1_upper(i)  - 1.0_dp*tmp1*int_edif_dfdy(i+1)/zetamult(i+1)
                END IF
            !END IF
    
            !IF((depth(i)/depth(i-1)<10.0_dp).and.(depth(i-1)/depth(i)<10.0_dp)) THEN
                IF(i>1) THEN
                    !M1_diag(i)   = M1_diag(i)   + 0.5_dp*tmp1*int_edif_dfdy(i)/zetamult(i)  ! Note that 1/zetamult(i)*Cbar = cb
                    !M1_lower(i)  = M1_lower(i)  + 0.5_dp*tmp1*int_edif_dfdy(i)/zetamult(i-1)
        
                    ! Experimenting with this approach to try to avoid negative
                    ! Cbar values.
                    IF(Cbar(i)/zetamult(i)<Cbar(i-1)/zetamult(i-1)) &
                             M1_diag(i)   = M1_diag(i)   + 1.0_dp*tmp1*int_edif_dfdy(i)/zetamult(i)  ! Note that 1/zetamult(i)*Cbar = cb
                    IF(Cbar(i)/zetamult(i)>=Cbar(i-1)/zetamult(i-1)) & 
                             M1_lower(i)   = M1_lower(i)   + 1.0_dp*tmp1*int_edif_dfdy(i)/zetamult(i-1)  ! Note that 1/zetamult(i)*Cbar = cb

                END IF
            !END IF 
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
            !M1_diag(i) = M1_diag(i) + wset/zetamult(i)  ! Note that 1/zetamult(i)*Cbar = cb


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

    ! New Cbar
    Cbar = XXX(1:a,1)

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
    ! A simple 'fix' is to make Cbed = 0.0 if negative, and remove the required
    ! amount of suspended sediment from the neighbouring point -- equivalent to
    ! assuming that the flux 'stopped' when the sediment ran out during the last
    ! time step.
    DO i = 1,a
        IF((Cbar(i)*depth(i)<0.0).and.(.FALSE.)) THEN
            ! Clip negligably small Cbar values
            !IF(abs(Cbar(i))<1.0e-10) THEN
            !    Cbar(i) = 0._dp
            !    !goto 21121            
            !    CYCLE
            !END IF
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
                    Cbar(i) = 1.e-012_dp
                END IF

            ELSE
                ! Treat boundary values
                IF(i==1) THEN
                    IF(depth(i+1)>0._dp) THEN
                        Cbar(i+1) =(Cbar(i+1)*depth(i+1)+Cbar(i)*depth(i))/depth(i+1)
                        Cbar(i) = 1.0e-012_dp
                    ELSE
                        ! Isolated point
                        Cbar(i) = 1.e-012_dp
                    END IF 
                ELSE IF(i==a) THEN
                    IF(depth(i-1)>0._dp) THEN
                        Cbar(i-1) =(Cbar(i-1)*depth(i-1)+Cbar(i)*depth(i))/depth(i-1)
                        Cbar(i) = 1.0e-012_dp
                    ELSE
                        Cbar(i)=1.e-12_dp
                    END IF
                END IF

            END IF
    21121 CONTINUE
        END IF
    END DO
    
    DO i=1,a
        IF(Cbar(i)<0.0_dp) THEN
            IF(Cbar(i)< -1.0e-012_dp) print*, 'Cbar clip', i, Cbar(i), Cbar_old(i), depth(i)
            Cbar(i) = 0.0e-12_dp
        END IF
    END DO
    !IF((.TRUE.).and.(minval(Cbar)<0._dp)) THEN
    !    print*, 'Cbar clip', minval(Cbar), minval(Cbar)*depth(minloc(Cbar))
    !    Cbar = max(Cbar, 0._dp)
    !    IF(counter.eq.1) print*, 'WARNING: Negative Cbar values are being clipped &
    !                            to zero'
    !    ! IDEA -- try approaching this with an approach along the lines of the
    !    ! one in that paper that Steve Roberts pointed you to.
    !END IF 

   
 
    ! Store the lateral flux of suspended load = -cb*int_edif_dfdy -
    ! dcb/dy*int_edif_f
    lat_sus_flux = 0.0_dp
    DO i=0,a
        ! Near bed suspended sediment concentration at i+1, i
        IF((i<a).and.(i>0)) THEN
            cbed_tmp1 = Cbar(i+1)/zetamult(i+1) - Cbar(i)/zetamult(i) ! diff(cbed)
            !cbed_tmp2 = cbed(i)
            IF(Cbar_old(i)/zetamult(i)<Cbar_old(i+1)/zetamult(i+1)) THEN
                cbed_tmp2 = Cbar(i)/zetamult(i)
            ELSE
                cbed_tmp2 = Cbar(i+1)/zetamult(i+1)
            END IF
        ELSE
            cbed_tmp1 = 0.0_dp
            cbed_tmp2 = 0.0_dp
        END IF
        ! Compute lat sus flux at i+1/2
        ! e.g. lat_sus_flux(1) = flux at 1/2 
        !      lat_sus_flux(a+1) = flux at a+1/2
        lat_sus_flux(i+1) = -int_edif_f(i+1)*(cbed_tmp1)/(ys_temp(i+1)-ys_temp(i))
        lat_sus_flux(i+1) = lat_sus_flux(i+1) -(cbed_tmp2)*int_edif_dfdy(i+1) 
    END DO 

    ! Add deposition here using operator splitting
    ! depth*dCbar/dt +wset*Cbed =  0 
    ! depth/dt*(Cbar_new -Cbar) + wset*(Cbar_new/zetamult) = 0
    ! Cbar_new( depth/dt + wset/zetamult) = 0 + depth/dt*Cbar
    Cbar = (Qe*0.0_dp + depth(1:a)/delT*Cbar - 0.0_dp*wset/zetamult(1:a)*Cbar)/ &
           (depth(1:a)/delT + 1.0_dp*wset/zetamult(1:a))


    

    ! Convert back to kg/m^3
    Cbar = Cbar*rhos

    ! Compute new near bed cb
    DO i = 1, a
        IF(zetamult(i) > 0._dp) THEN
            cb(i) = Cbar(i)/zetamult(i)
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

    ! num_trapz_pts = number of points used in trapezoidal integration
    INTEGER:: i, num_trapz_pts=100 
    REAL(dp):: db_const, F1, J1, j, E2, z2, perturb = 1.0e-05_dp, eps, d_eps
   
    IF((z>20.0_dp).or.(d_aref>1.0_dp)) THEN
        ! If z is large, or aref is larger than the depth, there is no suspended load, make a quick exit
        rouse_int = 0.0_dp 
        
    ELSE
        !Proceed with the Guo and Julien algorithm        

        ! Prevent integer values of z, by 'perturbing' the z value away from an
        ! integer if needed. 
        IF((abs(z-anint(z))<perturb).and.(anint(z)>0)) THEN
            IF(z<anint(z)) THEN 
                z2 = anint(z)-perturb
            ELSE
                z2 = anint(z) + perturb
            END IF
        ELSE
            z2=z
        END IF

        ! Compute 1/((1-deltab)/deltab)^z2
        db_const = ((1.0_dp-d_aref)/d_aref)**(-z2)

        ! Compute F1, eqn 8
        IF(d_aref<0.1_dp) THEN
            ! Here is the Guo and Julien algorithm -- it reportedly only
            ! converges in the infinite sum for d_aref<0.5. 
            ! FIXME: I chose the 0.1 after some bad experiences, although I am
            ! not sure about the convergence behaviour
            F1 =((1.0_dp-d_aref)**z2)/d_aref**(z2-1.0_dp) 
            E2 = d_aref/(1.0_dp-d_aref)

            DO i=1,10
                j=i*1.0_dp
                F1 = F1 - z2*((-1)**j)/(j-z2)*(E2)**(j-z2)
            END DO

            ! Compute J1, eqn 9
            J1 = z2*pi/sin(z2*pi) -F1

        ELSE
            ! Here we compute J1 using a brute-force trapezoidal rule
            J1 = 0.0_dp

            ! Trapezoidal Rule integration
            d_eps = ((1.0_dp-d_aref)/(num_trapz_pts*1._dp))
            DO i=1,num_trapz_pts
                ! Evaluate the function at at (0.5, 1.5, 2.5, ...99.5) / 100 of the
                ! integration domain, 
                eps = d_aref + ((i-0.5_dp))*d_eps
                ! J1 = J1 + f(eps,z)*deps
                J1 = J1 + ( (1.0_dp-eps)/eps)**(z)*d_eps
            END DO       

        END IF

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
SUBROUTINE int_edify_f(edify_model,sus_vert_prof,& 
                      a,ys,bed,ysl, ysu, bedl, bedu, &
                      water, ustar,wset, a_ref, &
                      int_edif_f, int_edif_dfdy, no_subints)
    ! PURPOSE: 
    !   To calculate
    !
    ! Integral ( edify*f) dz
    ! and,
    ! Integral (edify*df/dy) dz
    !
    ! Where z is a vertical coordinate,
    ! edify is the horizontal eddy diffusivity (varying in the vertical), and 
    ! f defines the vertical profile of suspended sediment ( so c(z) = cbed*f(z) )
    !
    ! These integral is useful in computing the vertically integrated lateral flux of
    ! suspended load = 
    ! INT (edify*dc/dy) dz = 
    ! INT(edify*cb*df/dy + edify*f*dcb/dy) dz =
    ! cb*INT(edify*df/dy) dz + dcb/dy*INT(edify*f) dz
    !
    ! OUTPUTS: 
    !   int_edif_f = Integral ( edify*f) dz, evaluated between grid points (i.e. 0.5, 1.5, ...a+0.5)
    !   int_edif_dfdy = Integral ( edify*df/dy) dz, evaluated between grid points (i.e. 0.5, 1.5, ...a+0.5)
    !
    INTEGER, INTENT(IN):: a, no_subints
    CHARACTER(len=20), INTENT(IN):: edify_model, sus_vert_prof
    REAL(dp), INTENT(IN):: ys, bed, ysl, ysu, bedl, bedu, ustar, water, wset, a_ref
    REAL(dp), INTENT(OUT):: int_edif_f, int_edif_dfdy
    DIMENSION ys(a), bed(a), ustar(a), int_edif_f(a+1), int_edif_dfdy(a+1), a_ref(a)

    ! Local variables
    INTEGER:: i, j

    REAL(dp):: d, us,bedh, f(no_subints), edify(no_subints), df_dy(no_subints), z_tmp(no_subints), &
                bed_tmp(0:a+1), ustar_tmp(0:a+1), &
                eps_z, ys_tmp(0:a+1), dbed_dy, depsz_dy, aref_tmp(0:a+1),&
                arefh, daref_dy, dus_dy, df_dbedh(no_subints), df_darefh(no_subints), df_dus(no_subints), &
                tmp2(no_subints), z2surf(no_subints), z2bed(no_subints), z2ratio(no_subints), dz

    ! Predefine bed_tmp, ys, and ustar, including boundary values
    bed_tmp(1:a) = bed
    bed_tmp(0)   = bedl
    bed_tmp(a+1) = bedu

    ys_tmp(1:a) = ys
    ys_tmp(0)   = ysl
    ys_tmp(a+1) = ysu

    ustar_tmp(1:a) = ustar
    ustar_tmp(0)   = ustar(1)
    ustar_tmp(a+1) = ustar(a)
    
    aref_tmp(1:a) = a_ref(1:a)
    aref_tmp(0)   = a_ref(1)
    aref_tmp(a+1) = a_ref(a)

    ! Loop through (0.5, 1.5, ... a+0.5) to compute integrals
    DO i=1,a+1

        ! Define depth, ustar, bed, aref, epsz, at i-1/2
        d = max(water - 0.5_dp*(bed_tmp(i)+bed_tmp(i-1)), 0.0_dp) 
        us= 0.5_dp*(ustar_tmp(i)+ustar_tmp(i-1))
        bedh = 0.5_dp*(bed_tmp(i)+bed_tmp(i-1))
        arefh = 0.5_dp*(aref_tmp(i)+aref_tmp(i-1))

        ! Define daref/dy and dbed/dy, dustar/dy at i-1/2
        daref_dy = (aref_tmp(i) - aref_tmp(i-1))/(ys_tmp(i)-ys_tmp(i-1)) 
        dbed_dy = (bed_tmp(i) - bed_tmp(i-1))/(ys_tmp(i)-ys_tmp(i-1)) 
        dus_dy = (ustar_tmp(i) - ustar_tmp(i-1))/(ys_tmp(i)-ys_tmp(i-1)) 
        

        ! Create vertical suspended sediment profile
        SELECT CASE(sus_vert_prof) 
            CASE('exp') 
                ! Vertical eddy diffusivity
                eps_z =0.5_dp*( 0.1_dp*ustar_tmp(i)*max(water-bed_tmp(i),0.0_dp) + &
                                0.1_dp*ustar_tmp(i-1)*max(water-bed_tmp(i-1),0.0_dp)) 
                eps_z = max(eps_z,1.0e-08_dp) 

                depsz_dy = ( 0.1_dp*ustar_tmp(i)*max(water-bed_tmp(i),0.0_dp) - &
                             0.1_dp*ustar_tmp(i-1)*max(water-bed_tmp(i-1),0.0_dp) &
                             )/(ys_tmp(i)-ys_tmp(i-1))

                !z_tmp = elevation above bed = at 0.5, 1.5, ... 99.5 * depth/no_subints.0 
                z_tmp = bedh + (d/no_subints*1.0_dp)*( (/ (j,j=1,no_subints) /)*1.0_dp-0.5_dp)

                ! Exponential vertical suspended sediment profile
                f = exp(-wset/eps_z*(z_tmp-bedh))
                
                ! df/dy = df/deps_z * deps_z/dy + df/dbed*dbed/dy
                df_dy = wset*(z_tmp-bedh)/eps_z**2*f*depsz_dy + &
                                wset/eps_z*f*dbed_dy

            CASE('Rouse')
                 
                IF(d>arefh)THEN  !((d*0.3_dp>arefh)) THEN !.and.(us>wset)) THEN
                    !FIXME: d*0.3_dp > arefh --- this is to be consistent with
                    !the 'rouse_int' function, which only converges for
                    !d>2arefh (I apparently built in some safety after a bad
                    !experience)

                    !z_tmp = elevation above bed = at 0.5, 1.5, ... 99.5 * depth/no_subints.0,
                    ! adjusted so z>arefh 
                    z_tmp = bedh+arefh+ ((d-arefh)/no_subints*1.0_dp)*( (/ (j,j=1,no_subints) /)*1.0_dp-0.5_dp)

                    ! Useful shorthand variables, which save computation
                    !z2surf= water-z_tmp ! Distance from z_tmp to the surface
                    z2bed = z_tmp-bedh  ! Distance from z_tmp to the bed
                    z2ratio = d/z2bed ! A ratio that comes up a lot

                    ! Calculate vertical profile of suspended sediment
                    f = ( (z2ratio-1.0_dp)/(d/arefh -1.0_dp))**(wset/(0.4_dp*us))
                    !IF(minval(f)<0._dp) THEN
                    !    print*, ' MINVAL f < 0._dp'
                    !    stop
                    !END IF

                    ! Calculate derivative of f. Use this approach:
                    ! df_dy = df/dbedh*dbedh/dy + df/aref*daref/dy + df/dus*dus/dy
                    ! because all the d/dy terms have only a single value at
                    ! i-1/2. 

                    ! Step1: df/dbedh, evaluated using maxima (symbolic algebra) 
                    ! -- see code in the file lat_flux.max
                    ! [2.5E+0*wset*((Y-h)/aref-1)*(((Y-h)/(z-h)-1)/((Y-h)/aref-1))**(2.5
                    !     1   E+0*wset/ustar)*(((Y-h)/(z-h)-1)/(aref*((Y-h)/aref-1)**2)+((Y-h
                    !     2   )/(z-h)**2-1/(z-h))/((Y-h)/aref-1))/(ustar*((Y-h)/(z-h)-1))]
                    df_dbedh = &
                            (wset/0.4_dp)*(d/arefh-1.0_dp)*f* &
                            ((z2ratio-1.0_dp)/(arefh*(d/arefh-1.0_dp)**2)+&
                            (d/z2bed**2-1.0_dp/z2bed)/(d/arefh-1.0_dp))/(us*(z2ratio-1.0_dp))

                    ! Step2: df/darefh, evaluated using maxima (symbolic algebra) 
                    ! -- see code in the file lat_flux.max
                    ! [2.5E+0*wset*(Y-h)*(((Y-h)/(z-h)-1)/((Y-h)/aref-1))**(2.5E+0*wset/
                    !     1   ustar)/(aref**2*ustar*((Y-h)/aref-1))]
                    df_darefh = (wset/0.4_dp)*d*f/(arefh**2*us*(d/arefh-1.0_dp))

                    ! Step3: df/dus, evaluated using maxima (symbolic algebra) 
                    ! -- see code in the file lat_flux.max
                    ! [-2.5E+0*wset*(((Y-h)/(z-h)-1)/((Y-h)/aref-1))**(2.5E+0*wset/ustar
                    !     1   )*log(((Y-h)/(z-h)-1)/((Y-h)/aref-1))/ustar**2]
                    df_dus = -(wset/0.4_dp)*f*log((z2ratio-1.0_dp)/(d/arefh-1.0_dp))/us**2

                    ! Step4: df_dy = df/dbedh*dbedh/dy + df/aref*daref/dy + df/dus*dus/dy
                    df_dy = df_dbedh*dbed_dy + df_darefh*daref_dy + df_dus*dus_dy
                ELSE
                    !PRINT*, 'ERROR - d< aref in suspended_xsect (int_edify_f)', &
                    !        d, arefh, aref_tmp(i), aref_tmp(i-1), bed_tmp(i), bed_tmp(i-1)
                    !stop
                    !z_tmp = elevation above bed = at 0.5, 1.5, ... 99.5 * depth/no_subints.0,
                    
                    ! adjusted so z>arefh 
                    !z_tmp = bedh+arefh+ ((d-arefh)/no_subints*1.0_dp)*( (/ (j,j=1,no_subints) /)*1.0_dp-0.5_dp)

                    z_tmp = bedh+ ((d-bedh)/no_subints*1.0_dp)*( (/ (j,j=1,no_subints) /)*1.0_dp-0.5_dp)
                    ! In these shallow waters, define things so there is no
                    ! lateral flux of suspended load -- hmm, actually, not such
                    ! a good idea?
                    f = 1.0_dp ! Uniform suspended load in very shallow water?
                    df_dy= 0.0_dp ! FIXME: Is this appropriate?
                END IF

            CASE DEFAULT
                print*, 'ERROR: Invalid value of sus_vert_prof in int_edify_f'
                stop

        END SELECT

        ! Create horizontal eddy diffusivity profile     
        SELECT CASE(edify_model)

            CASE('constant')
                IF(d>0.0_dp) THEN
                    edify = 0.24*d*us ! Constant
                ELSE
                    edify = 0.0_dp
                END IF

            CASE('Parabolic')
                IF(d>0.0_dp) THEN
                    edify = 0.4_dp*us*d*(z_tmp-bedh)*(water-z_tmp)/(0.25_dp*d**2) ! Parabolic
                ELSE
                    edify=0.0_dp
                END IF

            CASE('Parabola_const')
                ! Parabolic lower half, constant upper half
                IF(d>0.0_dp) THEN
                    edify = 0.4_dp*us*d*min(z_tmp-bedh, d*0.5_dp)*max(water-z_tmp,d*0.5_dp)/(0.25_dp*d**2)
                ELSE
                    edify=0.0_dp
                END IF

            CASE DEFAULT
                print*, 'ERROR: Invalid value of edify_model in int_edify_f'
                stop
                
        END SELECT    


        !Integral (edify*f) dz
        dz = (z_tmp(no_subints)-z_tmp(1))/(1.0_dp*(no_subints-1))
        int_edif_f(i) = sum(edify*f)*dz
        
        !Integral (edify*df/dy) dz
        int_edif_dfdy(i) = sum(edify*df_dy)*dz

        !print*,'#########', i, int_edif_f(i), int_edif_dfdy(i) 
        !print*, edify, df_dy
    END DO
   

END SUBROUTINE int_edify_f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

END MODULE suspended_xsect
