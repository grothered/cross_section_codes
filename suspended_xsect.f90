MODULE suspended_xsect

!File of global parameters
use global_defs
use util_various
IMPLICIT NONE
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dynamic_sus_dist(a, delT, ys, bed, water, waterlast, Q, tau, vel, wset, Qe,lambdacon, &
                                rho,rhos, g, d50, bedl,bedu, ysl, ysu, cb, Cbar, Qbed,Qbedon, &
                                sed_lag_scale, counter, a_ref, sus_vert_prof, edify_model, &
                                x_len_scale, sconc, lat_sus_flux, bedlast, int_edif_f, int_edif_dfdy, zetamult, &
                                too_steep)
    ! Calculate the cross-sectional distribution of suspended sediment using
    ! some ideas from /home/gareth/Documents/H_drive_Gareth/Gareth_and_colab
    ! s/Thesis/Hydraulic_morpho_model/channel_cross_section/paper/idea_for_
    ! simplified_cross_sectional_sus
    
    ! We solve the equation
    !
    ! d (depth*Cbar) / dt + d(U*depth*Cbar)/dx +
    ! d(V*depth*Cbar)/dy = - d/dy( Fl ) +
    ! (Es - ws*cb)
    !
    ! With suitable approximations for d(U*depth*Cbar)/dx and V*depth, and if needed, numerical
    ! integration to calculate Fl (the lateral flux)
    
    INTEGER, INTENT(IN)::a, counter, too_steep
    LOGICAL, INTENT(IN):: Qbedon
    REAL(dp), INTENT(IN):: delT, ys, bed, water, waterlast, tau, vel, wset, Qe, lambdacon, rho, rhos,g, & 
                                d50, bedl, bedu,ysl,ysu, Q, Qbed, a_ref, x_len_scale, sconc, bedlast
    REAL(dp), INTENT(OUT):: lat_sus_flux 
    REAL(dp), INTENT(INOUT):: int_edif_f, int_edif_dfdy, sed_lag_scale , zetamult
    ! int_edif_f, int_edif_dfdy = 2 integrals that appear in the lateral diffusive flux
    ! sed_lag_scale = used to estimate derivative terms, e.g.:
    ! d(U*d*Cbar)/dx = (U*d*Cbar - (U*d*Cbar)_ideal)/x_len_scale
    ! where (U*d*Cbar)_ideal is computed so that the 'ideal' sediment flux is
    ! equal to a desired value.
    ! cb = Near bed suspended sediment concentration, 
    ! Cbar = Depth averaged suspended sediment concentration
    REAL(dp), INTENT(IN OUT):: cb, Cbar
    CHARACTER(char_len), INTENT(IN):: sus_vert_prof, edify_model
    DIMENSION ys(a), bed(a), tau(a),vel(a), Qe(a), cb(a), Cbar(a),Qbed(a), a_ref(a), lat_sus_flux(a+1), bedlast(a), &
              int_edif_f(a+1), int_edif_dfdy(a+1), zetamult(0:a+1), too_steep(a)

    ! LOCAL VARIABLES
    INTEGER:: i, info, j
    LOGICAL:: halt, xderivative_operator_splitting=.FALSE., erode_deposit_splitting=.FALSE.
    ! NOTE: BE CAREFUL WITH USING OPERATOR SPLITTING IT CAN NEGATIVELY THE
    ! AFFECT ABILITY OF THE CODE TO REACH STEADY STATE / CONVERGE IN TIME. 

    ! 14/2/2012
    ! The following tests suggest that setting both splitting flags to FALSE is
    ! a good idea, and that any other choice will probably require great great
    ! care. 

    ! I ran a simulation (no_22 in the xsect_tests folder) with each combination
    ! of xderivative = True/False, and erode_deposit_splitting=True/False. I
    ! denote tf to be xderivative(True), erode_deposit_splitting (false), and
    ! likewise define ff, tt, and ft. 

    ! I first ran a test with 'Zero' lateral diffusion, rouse vertical profile,
    ! running the code to an equilibrium with every combination. 'tf' got a
    ! different equilibrium answer to each of the other 3 combinations of these
    ! variables (which all agreed with each other in terms of equilibrium
    ! solutions). However, the time-evolution of tt was different to ft and ff,
    ! which seemed the same as each other.

    ! I then re-ran the above with twice as many points in space. All seemed to
    ! converge to the same solution as with half as many points, although tf has
    ! larger errors. These 'fine' runs also seemed to behave the same in time as
    ! with half as many points ===> Little effect of spatial discretization on
    ! equilibrium or evolution. 

    ! WHAT ABOUT IF WE ADD LATERAL DIFFUSION?
    ! I ran the same check as above with lateral diffusion of suspended load,
    ! using the 'parabolic' model. 
    ! With these settings, ff and ft do not agree in time anymore on the
    ! coarsest grid, although based on an initial run, it seemed plausible that
    ! they may be approaching the same equilibrium solution (I didn't run either
    ! for enough time to find equilibrium, although for ff, the channel change
    ! became quite slow) 
    !
    ! If we double the number of points, then ff required a smaller time-step
    ! (selected as *0.5) to be stable enough to generally preserve symmetry. The
    ! evolution in time seemed very similar in this case to the case with less
    ! grid points. I further tried a *0.25 timestep, to check for strong
    ! dependende on the ratio of dt/dx -- however, this solution seemed very
    ! similar to the other fine-grid solution. ==> ff appears to be approx
    ! converged
    !
    ! For ft, I also tried doubling the number of points. The agreement in time
    ! was not so good with the other ft, it was somewhat between the ff
    ! solutions and the coarse ft one. I tried shrinking the time-step to *0.25,
    ! and the tendency was again to be moving closer to ff. Dropping the
    ! time-step it to *0.125, it was pretty close to ff. Using *0.0625, it was
    ! closer still to the ff case. Using *0.03125, it got closer to the ff case
    ! with twice as many points, which is still quite close to the ff case.
    
    ! CONCLUSIONS: This analysis doesn't show any version to be correct,
    ! however, ff seems to have the best convergence properties, requiring much
    ! less effort than ft (which however seems to give the same solution as ff
    ! with a small enough time-step / sufficient number of points). This must be
    ! because of the operator splitting, as nothing else differs between the 2
    ! runs.
    
    ! So operator splitting seems to require great care in this example --
    ! probably because the time-step I take is pretty large compared with the
    ! relevant timescales for diffusion?

    ! In addition, beware large timesteps in general. The above coarse runs were
    ! done with a 30s time-step (and 2000 pts over 150m), halving as described.
    ! No remeshing was used.  Although the ff case behaved pretty well even with
    ! these settings, I have often used much larger time-steps, at my peril.

    ! I will now re-investigate the equilibrium of the no-22 case with 4000pts
    ! and a 10s time-step, as these seem to be moderately conservative values.
    ! .... -- This run does seem to reach a convergent solution. 
    ! - Does it have sed_lag_scale =1.0? {Does seem to be converging to this, e.g.
    ! sed_lag_scale=0.9999989 when the bed_change is around 3.0e-10 ) 
    ! - Is the equilibrium sufficiently grid-independent? {A run with half as many
    ! points and 2x the time step did not converge well. A finer run seems to
    ! be similar, but still significantly different near the banks.}
    ! - Can we use re-meshing to get the fine grid near the banks and still get as
    ! good a solution, and does this require that we keep the time-step small in
    ! keeping with the finer mesh? { Still running, looks as though the
    ! remeshing run with half as many points may get closer to the finer run
    ! with 2x as many points than does the standard run. CHECK}.
    ! THEN ONCE WE UNDERSTAND THIS
    ! -- how does it go on more challenging problems like chi_1?

    REAL(dp):: depth(0:a+1), eddif_y(0:a+1), eddif_z(a), vd(0:a+1), ys_temp(0:a+1)
    REAL(dp):: M1_lower(a), M1_diag(a), M1_upper(a), M1_upper2(a), M1_lower2(a)
    REAL(dp):: RHS(a), dy_all(a), depthlast(0:a+1), zetamult_old(0:a+1),&
                diffuse1(0:a), diffuse2(0:a), cb_weight(1:a-1)
    REAL(dp):: tmp1, dy, dy_outer,dy_outer_inv, tmp2,dy_inner_inv, Cref, z,&
                cbed_tmp1, cbed_tmp2, tmp3, tmp4, discharge
    REAL(dp):: dQdx, dhdt, Cbar_old(a), dUd_dx(0:a+1), sus_flux, impcon, d1, d2
    REAL(dp):: DLF(a), DF(a), DUF(a), DU2(a),rcond, ferr, berr, work(3*a), XXX(a, 1)
    REAL(dp):: bandmat(5,a), AFB(7,a), RRR(a), CCC(a), int_edif_f_old(a+1), int_edif_dfdy_old(a+1)
    INTEGER::  IPV(a), iwork(a), imax   
    CHARACTER(1):: EQUED

    ! This routine calculates cb, Cbar in units m^3/m^3 --- however, elsewhere
    ! they are in kg/m^3 --- so we convert here, and convert back at the end of
    ! the routine
    cb = cb/rhos
    Cbar = Cbar/rhos 


    ! Degree of implicitness of the lateral diffusive flux
    impcon=1.00_dp
    IF(impcon .ne. 1.0_dp) THEN
        print*, 'ERROR: impcon in dynamic_sus_dist is not presently supported for values other than 1.0_dp' 
    END IF

    Cbar_old = Cbar
    int_edif_f_old = int_edif_f
    int_edif_dfdy_old = int_edif_dfdy
    zetamult_old = zetamult

    ! Depth, including at zero-depth boundaries ( 0 and a+1 )
    depth(1:a) = max(water -bed, 0._dp)
    depth(0) = 0._dp
    depth(a+1) = 0._dp

    depthlast(1:a) = max(waterlast-bedlast,0.0_dp)
    depthlast(0)=0._dp
    depthlast(a+1) = 0._dp

    ! y coordinate, including at zero depth boundaries (0 and a+1)
    ys_temp(1:a) = ys
    ys_temp(0) = ys(1) - (ys(1) - ysl)*abs(water-bed(1))/(bedl-bed(1))!ysl
    ys_temp(a+1) = ys(a) + (ysu-ys(a))*abs(water-bed(a))/(bedu-bed(a)) !ysu


    ! Calculate the value of 'zetamult', where:
    ! zetamult*cbed = Cbar = depth averaged sediment concentration
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
                    zetamult(i) = rouse_int(z,a_ref(i)/depth(i), Qbedon)
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
    
    ! Should ensure cb <= ~ 150g/L, and make sure we dont divide by zero
    !zetamult(1:a) = max(zetamult(1:a), Cbar/(150.0_dp/rhos), 1.0e-010_dp) 
    !zetamult_old(1:a) = max(zetamult_old(1:a), Cbar_old/(150.0_dp/rhos), 1.0e-010_dp) 
    zetamult(1:a) = max(zetamult(1:a),1.0e-10_dp)
    IF(any(zetamult(1:a)*150.0_dp/rhos < Cbar(1:a))) THEN
        PRINT*, 'ERROR: cb will apparently be > 150g/L'
        DO i=1,a
            PRINT*, i, zetamult(i)*150.0_dp/rhos, Cbar(i),zetamult(i), tau(i), &
                    a_ref(i), wset/(0.4*sqrt(tau(i)/rho)), depth(i)
        END DO
        STOP
    END IF


    ! Compute the discharge using the same approach as is used to compute
    ! sus_flux
    ! FIXME: Note that this method of computing the discharge is not 100%
    ! consistent with the method in hydro_xsect, because it uses the ysu, ysl
    ! boundaries, instead of the linearly interpolated approximations of the
    ! boundaries. The difference is very minor though. 
    discharge = sum( abs(vel)*max(water-bed,0.0_dp)*& 
              ( ( (/ ys(2:a), ysu /) - (/ ysl, ys(1:a-1) /) )*0.5_dp) &  ! dy
                  )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! COMPUTE A HALF STEP OF d(depth*Cbar)/dt + d/dx(U*depth*Cbar) = 0
    !! This allows the d/dx derivative to be included in the calculations via
    !! operator splitting
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !IF(xderivative_operator_splitting) THEN
    !    imax=1
    !    DO i=1,imax
    !        ! Take 'imax' small time-steps, which in total sum to delT/2

    !        ! PUSH THE SUSPENDED FLUX TOWARDS THE DESIRED VALUE   
    !        ! Calculate total sediment flux at time = t.
    !        ! We will use this to 'correct' the sediment flux towards the desired value 
    !        sus_flux = sum( & 
    !                (Qbed+Cbar*abs(vel)*max(water-bed,0._dp) )*& ! Total load
    !                ( ( (/ ys(2:a), ysu /) - (/ ysl, ys(1:a-1) /) )*0.5_dp) &  ! dy
    !                      )

    !        ! Store old value of sed_lag_scale, to use below
    !        tmp2 = sed_lag_scale

    !        IF(sus_flux > 1.0e-12_dp) THEN
    !            sed_lag_scale = 1.0_dp*((sconc*discharge)/sus_flux) !Desired flux / actual flux

    !            ! Prevent very high or low values
    !            sed_lag_scale = min(max(sed_lag_scale,0.666_dp),1.5_dp) 
    !            !IF(mod(counter,1000).eq.1) PRINT*, 'sed_lag_scale = ', sed_lag_scale

    !        ELSE
    !            IF(sconc*discharge<sus_flux) THEN
    !                sed_lag_scale = 0.666_dp
    !            ELSEIF(sconc*discharge==sus_flux) THEN
    !                sed_lag_scale = 1.0_dp
    !            ELSE
    !                sed_lag_scale = 1.5_dp
    !            END IF

    !       END IF

    !        ! NON-CONSERVATIVE VERSION
    !        !! Now we add the term U*d*dCbar/dx using an operator splitting technique
    !        !! depth*dCbar/dT + depth*vel*dCbar/dx = 0.0
    !        !! Implicit 
    !        !!IF(maxval(Cbar)>0.0_dp) THEN
    !        !    Cbar = Cbar*(1.0_dp - (1.0_dp-impcon)*(delT/(2.0_dp*imax))*vel*(1._dp-tmp2)/x_len_scale)/ &
    !        !          (1.0_dp + impcon*(delT/(2.0_dp*imax))*vel*(1.0_dp-sed_lag_scale)/x_len_scale)
    !        
    !        ! CONSERVATIVE VERSION
    !        ! Now we add the term d(U*d*Cbar)/dx using an operator splitting technique
    !        ! d(depth*Cbar)/dT + d(depth*vel*Cbar)/dx = 0.0
    !        ! Here we are assuming that the spatially lagged value of depth*vel*Cbar 
    !        !  = depth*vel*Cbar/(actual_flux)*desired_flux -- i.e. same shape, but
    !        !  scaled so that the flux is exactly the desired flux.
    !        ! Implicit
    !        WHERE(depth(1:a)>0.0_dp)        
    !            Cbar = (depthlast(1:a)*Cbar_old)/ &
    !                   (depth(1:a) + (vel*depth(1:a)*delT/(2.0_dp*imax))*(1.0_dp-sed_lag_scale)/x_len_scale)
    !        ELSEWHERE
    !            Cbar = 0.0_dp
    !        END WHERE

    !    END DO
    !END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solve initially for the depth - averaged suspended sediment concentration
    !  d (depth*Cbar) / dt + 
    ! d(V*depth*Cbar)/dy - d/dy( Fl ) -
    ! (Es - ws*cb) = 0.
    !
    ! Note that the d/dx term is included via operator spliiting (see above and
    ! at the end of the code).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        ! Non-conservative version = depth*dCbar/dt
        !M1_diag(i) = M1_diag(i) + depth(i)/delT
        !RHS(i)     = RHS(i)     + Cbar_old(i)*depth(i)/delT
        ! Conservative version = d/dt (depth*Cbar)
        M1_diag(i) = M1_diag(i) + depth(i)/delT
        RHS(i)     = RHS(i)     + Cbar_old(i)*depthlast(i)/delT
        
    END DO
 
    IF(xderivative_operator_splitting.eqv..FALSE.) THEN 
        ! PUSH THE SUSPENDED FLUX TOWARDS THE DESIRED VALUE   
        ! Calculate total sediment flux at time = t.
        ! We will use this to 'correct' the sediment flux towards the desired value 
        sus_flux = sum( & 
                (Qbed+Cbar*abs(vel)*max(water-bed,0._dp) )*& ! Total load
                ( ( (/ ys(2:a), ysu /) - (/ ysl, ys(1:a-1) /) )*0.5_dp) &  ! dy
                      )

        IF(sus_flux > 1.0e-12_dp) THEN
            !sed_lag_scale = 1.0_dp*((sconc*discharge)/sus_flux)**5.0 !Desired flux / actual flux
            sed_lag_scale = ((sconc*discharge)/sus_flux) !Desired flux / actual flux

            ! Prevent very high or low values
            sed_lag_scale = min(max(sed_lag_scale,0.666_dp),1.5_dp) 
            !IF(mod(counter,1000).eq.1) PRINT*, 'sed_lag_scale = ', sed_lag_scale

        ELSE
            IF(sconc*discharge<sus_flux) THEN
                sed_lag_scale = 0.666_dp
            ELSEIF(sconc*discharge==sus_flux) THEN
                sed_lag_scale = 1.0_dp
            ELSE
                sed_lag_scale = 1.5_dp
            END IF

        END IF
       
        DO i=1,a
            M1_diag(i) = M1_diag(i) + vel(i)*depth(i)*(1.0_dp-sed_lag_scale)/x_len_scale
        END DO

    END IF
 
    DO i = 1, a
        !PRINT*, 'SANITY CHECK a'
        IF(isnan(M1_diag(i))) print*, 'M1_diag(', i,') is NaN a'    
        IF(isnan(M1_lower(i))) print*, 'M1_upper(', i,') is NaN a'    
        IF(isnan(M1_upper(i))) print*, 'M1_upper(', i,') is NaN a'    
        IF(isnan(RHS(i))) print*, 'RHS(', i,') is NaN a'    
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
    END DO

    DO i = 1,a
        ! Upwind discretization of the advective term -- fill out the matrices
        ! (FIXME -- consider making use of a higher order discretization)
        IF((i.ne.1).and.(i.ne.a)) THEN
            ! NON-CONSERVATIVE VERSION
            ! v*depth*d/dy(Cbar)
            !IF(vd(i) < 0.0) THEN 
            !    dy = ys_temp(i+1) - ys_temp(i)
            !    M1_upper(i) = M1_upper(i) + vd(i)/dy
            !    M1_diag(i)  = M1_diag(i)  - vd(i)/dy
            !ELSE
            !    dy = ys_temp(i) - ys_temp(i-1)
            !    M1_diag(i)  = M1_diag(i)  + vd(i)/dy
            !    M1_lower(i) = M1_lower(i) - vd(i)/dy
            !END IF
            
            ! CONSERVATIVE VERSION
            ! d/dy (v*depth*Cbar)
            IF(vd(i) < 0.0_dp) THEN 
                dy = ys_temp(i+1) - ys_temp(i)
                M1_upper(i) = M1_upper(i) + vd(i+1)/dy
                M1_diag(i)  = M1_diag(i)  - vd(i)/dy
            ELSE
                dy = ys_temp(i) - ys_temp(i-1)
                M1_diag(i)  = M1_diag(i)  + vd(i)/dy
                M1_lower(i) = M1_lower(i) - vd(i-1)/dy
            END IF
        END IF

    END DO

    DO i = 1, a
        !PRINT*, 'SANITY CHECK b'
        IF(isnan(vd(i))) print*, 'vd(', i,') is NAN'
        IF(isnan(M1_diag(i))) print*, 'M1_diag(', i,') is NaN b'    
        IF(isnan(M1_lower(i))) print*, 'M1_upper(', i,') is NaN b'    
        IF(isnan(M1_upper(i))) print*, 'M1_upper(', i,') is NaN b'    
        IF(isnan(RHS(i))) print*, 'RHS(', i,') is NaN b'    

    END DO

    ! Diffusion term
    DO i = 1, a
        ! Calculate a centred dy
        dy_outer = 0.5_dp*(ys_temp(i+1) - ys_temp(i-1))

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
                call int_edify_f(edify_model, sus_vert_prof, a, ys, bed, ys_temp(0), &
                                 ys_temp(a+1), water, water, water, sqrt(abs(tau)/rho),&
                                 wset,a_ref, Qbedon, int_edif_f, int_edif_dfdy) !, 205)

                ! Compute pointwise value of int(eddif_y*f)dz at (i+1/2)
                ! (denoted 'diffuse1') -- this is reused later
                diffuse1(0) = 0.0_dp
                diffuse1(a) = 0.0_dp
                !diffuse1(1:a-1) = 0.5_dp*max(int_edif_f(2:a) + int_edif_f(3:a+1), &
                !                            int_edif_f(2:a) + int_edif_f(1:a-1))
                diffuse1(1:a-1) = int_edif_f(2:a) 
                ! Here, we try to weight 'diffuse1' by the depth, biasing towards
                ! higher depths
                !diffuse1(1:a-1) = ((int_edif_f(2:a) + int_edif_f(3:a+1))*depth(2:a)**1 + &
                !                   (int_edif_f(2:a) +int_edif_f(1:a-1))*depth(1:a-1)**1  &
                !                  )/(depth(2:a)**1+depth(1:a-1)**1)
                
                ! Compute pointwise value of int(eddif_y*dfdy)dz at (i+1/2)
                ! (denoted 'diffuse2') -- this is reused later
                diffuse2(0) = 0.0_dp
                diffuse2(a) = 0.0_dp
                DO j=1,a-1
                    !diffuse2(j) = 0.5_dp*minmod(int_edif_dfdy(j+1)+int_edif_dfdy(j+2), &
                    !                              int_edif_dfdy(j+1)+int_edif_dfdy(j))
                    diffuse2(j)=int_edif_dfdy(j+1)
                    ! Here we try to weight diffuse2 by the depth, biasing
                    ! towards lower depths
                    !diffuse2(j) =( (int_edif_dfdy(j+1)+int_edif_dfdy(j+2))*depth(j)**1 + &
                    !               (int_edif_dfdy(j+1)+int_edif_dfdy(j))*depth(j+1)**1   &
                    !             )/(depth(j+1)**1+depth(j)**1)
 
                END DO

                ! Compute the weighting of cb(i+1/2) = (1-cb_weight)*cb(i+1) + cb_weight*cb(i)
                DO j=1,a-1
                    !IF(bed(j)>bed(j+1)) THEN
                    !    cb_weight(j) = 1.0_dp
                    !ELSEIF(bed(j)<bed(j+1)) THEN
                    !    cb_weight(j) = 0.0_dp
                    !ELSE
                    !    cb_weight(j) = 0.5_dp
                    !END IF
                    !cb_weight(j)=0.5_dp

                    ! Try to weight toward the cbed associated with the shallower depth
                    !cb_weight(j) = 1.0_dp - depth(j+1)**1/(depth(j)**1+depth(j+1)**1)

                    ! Try to weight toward the lower cbed value 

                    ! Special treatment of boundaries
                    !IF((j==1).OR.(j==a-1)) THEN
                    !    IF(cb(j+1)<cb(j)) THEN
                    !        cb_weight(j)=0.0_dp
                    !    ELSE
                    !        cb_weight(j)=1.0_dp 
                    !    END IF

                    !ELSE
                    cb_weight(j) = 0.50_dp + 0.5*(cb(j+1)-cb(j))/(cb(j)**1+cb(j+1)**1+1.0e-10_dp)
                        ! 15/2/2012
                        ! The above term would be the standard second order
                        ! approx if it were just 0.5_dp. The additional part
                        ! means we end up adding the following term to the
                        ! discretized pde (ignoring the 1.0e-10 which prevents
                        ! zero division):
                        ! 1/dy_outer*( 0.5*(cb(j+1)-cb(j))/(2*cb(j+1/2))*[-cb(j+1)+cb(j)]*int_edif_dfdy(i+1/2) - 
                        !              0.5*(cb(j)-cb(j-1))/(2*cb(j-1/2))*[-cb(j)+cb(j-1)]*int_edif_dfdy(i-1/2)
                        !            ),
                        ! where the part in square brackets follows from
                        ! factorising. This extra term helps prevent negative cb
                        ! values in the solution.
                        ! Noting that the (cb(j+1)-cb(j)) terms are proportional
                        ! to dy_inner*(dcb/dy)@(i+1/2), the extra term scales at
                        ! most with 1/dy_outer*(dy_inner)**2 ~= dy_outer.
                        ! Actually it probably scales with dy_outer**2 because
                        ! of the subtraction of similar terms at i+1/2 and
                        ! i-1/2), so goes to zero as dy_outer goes to zero.
                        ! Hence this should not ruin the convergence, at least
                        ! in smooth areas. 
                    !END IF
                END DO

            END IF
           
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            ! Add 
            ! d/dy [ dcb/dy*(int(eddif_y*f) dz) ]
            ! to the equations
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            dy_outer_inv = 1.0_dp/dy_outer
            dy_inner_inv = 1.0_dp/(ys_temp(i+1)-ys_temp(i))
           
            IF(i<a) THEN
                
                !diffuse1(i) = 0.5_dp*max(int_edif_f(i+1)+int_edif_f(i+2), int_edif_f(i+1)+int_edif_f(i))
                tmp4 = impcon*dy_outer_inv*dy_inner_inv*diffuse1(i)                

                ! Compute derivatives involving i+1 and i
                M1_upper(i) = M1_upper(i) - tmp4/zetamult(i+1)
                M1_diag(i)  = M1_diag(i)  + tmp4/zetamult(i)
        
                !tmp4 = (1.0_dp-impcon)*dy_outer_inv*dy_inner_inv*int_edif_f_old(i+1)
                !tmp4 = (1.0_dp-impcon)*dy_outer_inv*dy_inner_inv*&
                !        0.5_dp*minmod(int_edif_f_old(i+1)+int_edif_f_old(i+2), int_edif_f_old(i+1)+int_edif_f_old(i))
                !RHS(i)      = RHS(i) + tmp4/zetamult_old(i+1)*Cbar(i+1)
                !RHS(i)      = RHS(i) - tmp4/zetamult_old(i)*Cbar(i)

            END IF
            
            ! Redefine dy_inner_inv for computations of the derivatives
            ! involving i and i-1
            dy_inner_inv = 1.0_dp/(ys_temp(i)-ys_temp(i-1))
            IF(i>1) THEN

                tmp4 = impcon*dy_outer_inv*dy_inner_inv*diffuse1(i-1) 
                M1_diag(i)  = M1_diag(i)  + tmp4/zetamult(i)
                M1_lower(i) = M1_lower(i) - tmp4/zetamult(i-1)

                !tmp4 = (1.0_dp-impcon)*dy_outer_inv*dy_inner_inv*int_edif_f_old(i)
                !tmp4 = (1.0_dp-impcon)*dy_outer_inv*dy_inner_inv*&
                !        minmod(int_edif_f_old(i) + int_edif_f_old(i+1), int_edif_f_old(i)+int_edif_f_old(i-1))
                !RHS(i) = RHS(i) - tmp4/zetamult_old(i)*Cbar(i)
                !RHS(i) = RHS(i) + tmp4/zetamult_old(i-1)*Cbar(i-1)
                    
            END IF

            ! Add:
            ! d/dy ( cb*INT(edify*df/dy )dz )
            ! to the equations

            IF(i<a) THEN
                ! Try weighting cb (i+1/2) based on the depths at i and i+1 ---
                ! a potenial compromise between the centred approach, and just
                ! picking the one with the smallest depth.
                !tmp3=1.0_dp
                !tmp3 = (water-bed(i+1))**2/((water-bed(i+1))**2+(water-bed(i))**2)
                !if(tmp3>0.5_dp) tmp3 = min(tmp3, & 
                !                    (water-bedlast(i+1))**2/((water-bedlast(i+1))**2 + (water-bedlast(i))**2   ) &
                !                         )
                !if(tmp3<0.1_dp) tmp3 = 0.0_dp
                !if(tmp3>0.9_dp) tmp3 = 1.0_dp
                !!d1 = water-bed(i)
                !!d2 = water-bed(i+1)
                !! cb(i+1/2) = tmp3*(cb(i+1) ) + (1-tmp3)*cb(i)
                !! To ensure the positivity of cb, we may need to set tmp3 to 0
                !! or 1 where the depth changes rapidly, although where this is
                !! not an issue we expect 0.5 to be more accurate. 
                !IF((d1>1.5_dp*d2).or.(1.5_dp*d1<d2)) THEN
                    !tmp3=0.5_dp
                !ELSE
                !    tmp3 = 0.5_dp
                !END IF                

                tmp4 = impcon*dy_outer_inv*diffuse2(i)

                M1_upper(i) = M1_upper(i) - (1.0_dp-cb_weight(i))*tmp4/zetamult(i+1)  ! Note that 1/zetamult(i)*Cbar = cb
                M1_diag(i)  = M1_diag(i)  - cb_weight(i)*tmp4/zetamult(i)
               
                !tmp4 = (1.0_dp-impcon)*dy_outer_inv*int_edif_dfdy_old(i+1) 
                !RHS(i) = RHS(i) +  (1.0_dp-tmp3)*tmp4/zetamult_old(i+1)*Cbar(i+1)
                !RHS(i) = RHS(i) +  tmp3*tmp4/zetamult_old(i)*Cbar(i)
            END IF
    
            IF(i>1) THEN
                ! Try weighting cb (i-1/2) based on the depths at i and i-1 ---
                ! a potenial compromise between the centred approach, and just
                ! picking the one with the smallest depth.
                !tmp3 = (water-bed(i))**2/( (water-bed(i))**2+(water-bed(i-1))**2)
                !tmp3 = min(tmp3, & 
                !          (water-bedlast(i))**2/((water-bedlast(i))**2 + (water-bedlast(i-1))**2   ) &
                !        )
                !if(tmp3<0.1_dp) tmp3 = 0.0_dp
                !if(tmp3>0.9_dp) tmp3 = 1.0_dp
                !d1 = water-bed(i)
                !d2 = water-bed(i-1)
               
                !IF((d1>1.5_dp*d2).or.(1.5_dp*d1<d2)) THEN
                    !tmp3=0.5_dp
                !ELSE
                !    tmp3 = 0.5_dp
                !END IF

                tmp4 = impcon*dy_outer_inv*diffuse2(i-1)

                M1_diag(i)   = M1_diag(i)   + (1.0_dp-cb_weight(i-1))*tmp4/zetamult(i)  ! Note that 1/zetamult(i)*Cbar = cb
                M1_lower(i)  = M1_lower(i)  + cb_weight(i-1)*tmp4/zetamult(i-1)
    
                !tmp4 = (1.0_dp-impcon)*dy_outer_inv*int_edif_dfdy_old(i)
                !RHS(i) = RHS(i) -  (1.0_dp-tmp3)*tmp4/zetamult_old(i)*Cbar(i)
                !RHS(i) = RHS(i) -  tmp3*tmp4/zetamult_old(i-1)*Cbar(i-1)
            END IF
        END IF
    END DO

    !Sanity-check
    DO i = 1, a
    !    PRINT*, 'SANITY CHECK c'
        IF(isnan(vd(i))) print*, 'vd(', i,') is NAN, c'
        IF(isnan(M1_diag(i))) print*, 'M1_diag(', i,') is NaN, c', M1_diag(i)   
        IF(isnan(M1_lower(i))) print*, 'M1_lower(', i,') is NaN, c' , M1_lower(i), int_edif_dfdy(i), zetamult(i)  
        IF(isnan(M1_upper(i))) print*, 'M1_upper(', i,') is NaN, c' , M1_upper(i), int_edif_dfdy(i+1), zetamult(i+1)   
        IF(isnan(RHS(i))) print*, 'RHS(', i,') is NaN, c', RHS(i), Cbar(i),zetamult(i), Cbar_old(i),zetamult_old(i)   
    END DO

    ! Erosion and deposition
    IF(.NOT. erode_deposit_splitting) THEN
        DO i = 1, a

            !RHS(i) = RHS(i) +Qe(i) - (1.0_dp-impcon)*min(wset/zetamult_old(i), depthlast(i)/delT)*Cbar(i)
            !M1_diag(i) = M1_diag(i) + impcon*min(wset/zetamult(i), depth(i)/delT)  ! Note that 1/zetamult(i)*Cbar = cb
            RHS(i) = RHS(i) +Qe(i) - (1.0_dp-impcon)*wset*too_steep(i)/zetamult_old(i)*Cbar(i)
            M1_diag(i) = M1_diag(i) + impcon*wset*too_steep(i)/zetamult(i)  ! Note that 1/zetamult(i)*Cbar = cb

        END DO
    END IF

    !Sanity-check
    DO i = 1, a
        !PRINT*, 'SANITY CHECK d'
        IF(isnan(vd(i))) print*, 'vd(', i,') is NAN'
        IF(isnan(M1_diag(i))) print*, 'M1_diag(', i,') is NaN d' , zetamult(i)   
        IF(isnan(M1_lower(i))) print*, 'M1_upper(', i,') is NaN d'    
        IF(isnan(M1_upper(i))) print*, 'M1_upper(', i,') is NaN d'    
        IF(isnan(RHS(i))) print*, 'RHS(', i,') is NaN d'    

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

    ! Check for nan
    DO i=1,a
        IF((isnan(Cbar(i)))) THEN
            print*, 'Cbar', i, ' is nan right after matrix solve'
            DO j=1,a
                print*, j, M1_lower(j), M1_diag(j), M1_upper(j), RHS(j)
            END DO        
            stop
        END IF
    END DO
    
    !Check for lapack errors
    IF(info.ne.0) THEN
            print*, 'ERROR: info = ', info, ' in DGTSV, dynamic_sus_dist'
        stop
    END IF
   

    ! Check for negative Cbar values -- clip, and warn if they are not very
    ! small 
    DO i=1,a
        IF(Cbar(i)<0.0_dp) THEN
            IF(Cbar(i)< -1.0e-012_dp) THEN
                print*, 'Cbar clip', i, Cbar(i) !, int_edif_f(i-1:i+2), int_edif_dfdy(i-1:i+2)
                !IF((i>1).and.(i<a)) print*, depth(i-1),depth(i+1)
            END IF
            Cbar(i) = 0.0e-12_dp
        END IF
    END DO
   
 
    ! STORE THE LATERAL FLUX OF SUSPENDED LOAD = 
    ! -cb*int_edif_dfdy - dcb/dy*int_edif_f
    lat_sus_flux = 0.0_dp
    DO i=0,a
        ! Near bed suspended sediment concentration at i+1, i
        IF((i<a).and.(i>0)) THEN
            cbed_tmp1 = Cbar(i+1)/zetamult(i+1) - Cbar(i)/zetamult(i) ! diff(cbed)
            ! cbed(i+1/2)
            cbed_tmp2 = cb_weight(i)*Cbar(i)/zetamult(i) + (1.0_dp-cb_weight(i))*Cbar(i+1)/zetamult(i+1)
        ELSE
            cbed_tmp1 = 0.0_dp
            cbed_tmp2 = 0.0_dp
        END IF
        ! Compute lat sus flux at i+1/2
        ! e.g. lat_sus_flux(1) = flux at 1/2 
        !      lat_sus_flux(a+1) = flux at a+1/2
        lat_sus_flux(i+1) = -diffuse1(i)*(cbed_tmp1)/(ys_temp(i+1)-ys_temp(i))
        lat_sus_flux(i+1) = lat_sus_flux(i+1) -(cbed_tmp2)*diffuse2(i)
    END DO 

    
    !IF(xderivative_operator_splitting) THEN
   
    !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    ! Take a half time step of remaining terms 
    !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !    ! Compute the discharge using the same approach as is used to compute
    !    ! sus_flux
    !    ! FIXME: Note that this method of computing the discharge is not 100%
    !    ! consistent with the method in hydro_xsect, because it uses the ysu, ysl
    !    ! boundaries, instead of the linearly interpolated approximations of the
    !    ! boundaries. The difference is very minor though. 
    !    discharge = sum( abs(vel)*max(water-bed,0.0_dp)*& 
    !              ( ( (/ ys(2:a), ysu /) - (/ ysl, ys(1:a-1) /) )*0.5_dp) &  ! dy
    !                  )
    !    imax=1
    !    DO i=1,imax
    !        ! Take 'imax' small time-steps, which in total sum to delT/2
    !        !Cbar = (Qe + depth(1:a)/(delT/(2.0*imax))*Cbar - 0.5_dp*wset/zetamult_old(1:a)*Cbar)/ &
    !        !       (depth(1:a)/(delT/(2.0*imax)) + 0.5_dp*wset/zetamult(1:a))
    !    
    !        ! PUSH THE SUSPENDED FLUX TOWARDS THE DESIRED VALUE   
    !        ! Calculate total sediment flux at time = t.
    !        ! We will use this to 'correct' the sediment flux towards the desired value 
    !        sus_flux = sum( & 
    !                (Qbed+Cbar*abs(vel)*max(water-bed,0._dp) )*& ! Total load
    !                ( ( (/ ys(2:a), ysu /) - (/ ysl, ys(1:a-1) /) )*0.5_dp) &  ! dy
    !                      )

    !        ! Store old value of sed_lag_scale, to use below
    !        tmp2 = sed_lag_scale

    !        IF(sus_flux > 1.0e-12_dp) THEN
    !            sed_lag_scale = ((sconc*discharge)/sus_flux) !Desired flux / actual flux

    !            ! Prevent very high or low values
    !            sed_lag_scale = min(max(sed_lag_scale,0.666_dp),1.5_dp) 

    !        ELSE
    !            IF(sconc*discharge<sus_flux) THEN
    !                sed_lag_scale = 0.666_dp
    !            ELSEIF(sconc*discharge==sus_flux) THEN
    !                sed_lag_scale = 1.0_dp
    !            ELSE
    !                sed_lag_scale = 1.5_dp
    !            END IF

    !        END IF

    !        ! NON-CONSERVATIVE VERSION
    !        ! Now we add the term U*d*dCbar/dx using an operator splitting technique
    !        ! depth*dCbar/dT + depth*vel*dCbar/dx = 0.0
    !        ! Implicit 
    !        !    Cbar = Cbar*(1.0_dp - (1.0_dp-impcon)*(delT/(2.0_dp*imax))*vel*(1._dp-tmp2)/x_len_scale)/ &
    !        !           (1.0_dp + impcon*(delT/(2.0_dp*imax))*vel*(1.0_dp-sed_lag_scale)/x_len_scale)

    !        ! CONSERVATIVE VERSION
    !        ! Now we add the term d(U*d*Cbar)/dx using an operator splitting technique
    !        ! d(depth*Cbar)/dT + d(depth*vel*Cbar)/dx = 0.0
    !        ! Here we are assuming that the spatially lagged value of depth*vel*Cbar 
    !        !  = depth*vel*Cbar/(actual_flux)*desired_flux -- i.e. same shape, but
    !        !  scaled so that the flux is exactly the desired flux.
    !        ! Implicit 
    !        WHERE(depth(1:a)>0.0_dp)
    !            Cbar = (depth(1:a)*Cbar)/ &
    !                   (depth(1:a) + (vel*depth(1:a)*delT/(2.0_dp*imax))*(1.0_dp-sed_lag_scale)/x_len_scale)
    !        ELSEWHERE
    !            Cbar = 0.0_dp
    !        END WHERE
    !    END DO
    !  
    !END IF
    !
    !IF(erode_deposit_splitting) THEN
    !    ! Add deposition and erosion here using operator splitting
    !    ! depth*dCbar/dt +wset*Cbed =  Qe 
    !    ! depth/dt*(Cbar_new -Cbar) + wset*(Cbar_new/zetamult) = Qe
    !    ! Cbar_new( depth/dt + wset/zetamult) = Qe + depth/dt*Cbar
    !    !FIXME: Doing this separate to lateral diffusion may lead to negative Cbar
    !    !values, except with a fully implicit approach (impcon =1.0_dp)
    !    Cbar = 1.0_dp*(Qe*1.0_dp + depth(1:a)/delT*Cbar - 0.0_dp*wset/zetamult(1:a)*Cbar)/ &
    !           (depth(1:a)/delT + 1.0_dp*wset/zetamult(1:a))
    !END IF
    !    
    !! Check for negative Cbar, clip, and warn if it is not small 
    !DO i=1,a
    !    IF(Cbar(i)<0.0_dp) THEN
    !        IF(Cbar(i)< -1.0e-012_dp) print*, 'Cbar clip', i, int_edif_f(i:(i+2)), int_edif_dfdy(i:(i+2))
    !        Cbar(i) = 0.0e-12_dp
    !    END IF
    !END DO

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
REAL(dp) FUNCTION rouse_int(z,d_aref, Qbedon)
    ! Suppose that Cbar = cbed*(K +a_ref/d) [without bedload]
    !         or
    !              Cbar = cbed*K [with bedload]    
    ! Where Cbar is the depth-averaged sediment concentration
    ! cbed is the near bed concentration
    ! and K is an integrating factor, which is determined from the Rouse
    ! distribution for suspended sediment.
    ! Then this function calculates K + a_ref/d [without bedload]
    !                               or just K [with bedload]
    !
    !
    ! INPUTS
    ! z = rouse number = wset/(0.4*ustar)
    ! d_aref = dimensionless reference level 
    !        = (van rijn reference level)/ depth
    ! Qbedon = TRUE/FALSE [bedload 'on' or 'off']

    ! NOTE -- This uses a method from Guo and Julien (2004) 'Efficient
    ! Algorithms for Computing Einstein Integrals', Journal of Hydraulic
    ! Engineering 130:1198-1201.
    ! They show how to calculate J1=int_{d_aref}^{1} [(1-eps)/eps]^z d(eps)
    ! It can be shown that K=J1*db_const (where db_const is defined below) 
    REAL(dp), INTENT(IN):: z, d_aref
    LOGICAL, INTENT(IN):: Qbedon

    ! num_trapz_pts = number of points used in trapezoidal integration
    INTEGER:: i, num_trapz_pts=400 
    REAL(dp):: db_const, F1, J1, j, E2, z2, perturb = 1.0e-05_dp, eps, d_eps
   
    IF((d_aref>1.0_dp).or.(z>50.0_dp)) THEN
        ! If aref is larger than the depth, there is no suspended load, make a quick exit
        ! FIXME: Assuming that if z is large enough, we can neglect rouse_int
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

            ! Trapezoidal Rule integration -- FIXME: should really do Gauss quadrature.
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
        ! Include the contribution between [0, a_ref]
        ! FIXME: This should only be included when bedload is OFF
        IF(Qbedon.eqv..FALSE.) THEN
            rouse_int=rouse_int+d_aref 
        END IF

        IF((rouse_int<0.0_dp).or.(rouse_int>=1.0_dp).or.(rouse_int.ne.rouse_int)) THEN
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
                      water, ustar,wset, a_ref,Qbedon, &
                      int_edif_f, int_edif_dfdy) !, no_subints)
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
    !   int_edif_f = Integral_{0}^{water_surface} ( edify*f) dz
    !   int_edif_dfdy = Integral_{0}^{water_surface} ( edify*df/dy) dz
    !
    ! EVALUATED AT i+1/2
    !
    ! Note: In practice we integrate over [a_ref ,  water_surface] first, then
    ! add on the contribution from [0, a_ref] if Qbedon=FALSE . Otherwise, the
    ! latter zone is associated with bedload -- I think this is standard in other
    ! models (e.g. Delft3D)

    INTEGER, INTENT(IN):: a 
    LOGICAL, INTENT(IN)::Qbedon
    CHARACTER(char_len), INTENT(IN):: edify_model, sus_vert_prof
    REAL(dp), INTENT(IN):: ys, bed, ysl, ysu, bedl, bedu, ustar, water, wset, a_ref
    REAL(dp), INTENT(OUT):: int_edif_f, int_edif_dfdy
    DIMENSION ys(a), bed(a), ustar(a), int_edif_f(a+1), int_edif_dfdy(a+1), a_ref(a)

    ! Local variables
    INTEGER:: i, j

    REAL(dp):: d, us,bedh, f(64), edify(64), df_dy(64), z_tmp(64), &
                bed_tmp(0:a+1), ustar_tmp(0:a+1), &
                eps_z, ys_tmp(0:a+1), dbed_dy, depsz_dy, aref_tmp(0:a+1),&
                arefh, daref_dy, dus_dy, df_dbedh(64), df_darefh(64), df_dus(64), &
                z2ratio(64), dz, dyinv, d_on_aref_les1_inv, z2bed_inv(64), arefh_inv, &
                parabola(64), rouseno, tmp1(a+1), tmp2(a+1)
                !z2surf(64), z2bed(64)

    ! Note -- 64 point gaussian quadrature is used for vertical integration.
    ! Experiments suggest that this works very well, compared with e.g. 800
    ! points using simpsons rule or another higher order method -- get the same
    ! answer -- and it is way cheaper.  I got the coefficients from the web --
    ! see my 'math notes' folder for the website, with a folder on
    ! gaussian_quadrature there.
    REAL(dp):: gauss_weights(64)= (/ 0.0486909570091397_dp,0.0486909570091397_dp,0.0485754674415034_dp,0.0485754674415034_dp,&
                                     0.0483447622348030_dp,0.0483447622348030_dp,0.0479993885964583_dp,0.0479993885964583_dp,&
                                     0.0475401657148303_dp,0.0475401657148303_dp,0.0469681828162100_dp,0.0469681828162100_dp,&
                                     0.0462847965813144_dp,0.0462847965813144_dp,0.0454916279274181_dp,0.0454916279274181_dp,&
                                     0.0445905581637566_dp,0.0445905581637566_dp,0.0435837245293235_dp,0.0435837245293235_dp,&
                                     0.0424735151236536_dp,0.0424735151236536_dp,0.0412625632426235_dp,0.0412625632426235_dp,&
                                     0.0399537411327203_dp,0.0399537411327203_dp,0.0385501531786156_dp,0.0385501531786156_dp,&
                                     0.0370551285402400_dp,0.0370551285402400_dp,0.0354722132568824_dp,0.0354722132568824_dp,&
                                     0.0338051618371416_dp,0.0338051618371416_dp,0.0320579283548516_dp,0.0320579283548516_dp,&
                                     0.0302346570724025_dp,0.0302346570724025_dp,0.0283396726142595_dp,0.0283396726142595_dp,&
                                     0.0263774697150547_dp,0.0263774697150547_dp,0.0243527025687109_dp,0.0243527025687109_dp,&
                                     0.0222701738083833_dp,0.0222701738083833_dp,0.0201348231535302_dp,0.0201348231535302_dp,&
                                     0.0179517157756973_dp,0.0179517157756973_dp,0.0157260304760247_dp,0.0157260304760247_dp,&
                                     0.0134630478967186_dp,0.0134630478967186_dp,0.0111681394601311_dp,0.0111681394601311_dp,&
                                     0.0088467598263639_dp,0.0088467598263639_dp,0.0065044579689784_dp,0.0065044579689784_dp,&
                                     0.0041470332605625_dp,0.0041470332605625_dp,0.0017832807216964_dp,0.0017832807216964_dp /)

    REAL(dp)::gauss_abscissae(64)= (/-0.0243502926634244_dp,0.0243502926634244_dp,-0.0729931217877990_dp,0.0729931217877990_dp,&
                                     -0.1214628192961206_dp,0.1214628192961206_dp,-0.1696444204239928_dp,0.1696444204239928_dp,&
                                     -0.2174236437400071_dp,0.2174236437400071_dp,-0.2646871622087674_dp,0.2646871622087674_dp,&
                                     -0.3113228719902110_dp,0.3113228719902110_dp,-0.3572201583376681_dp,0.3572201583376681_dp,&
                                     -0.4022701579639916_dp,0.4022701579639916_dp,-0.4463660172534641_dp,0.4463660172534641_dp,&
                                     -0.4894031457070530_dp,0.4894031457070530_dp,-0.5312794640198946_dp,0.5312794640198946_dp,&
                                     -0.5718956462026340_dp,0.5718956462026340_dp,-0.6111553551723933_dp,0.6111553551723933_dp,&
                                     -0.6489654712546573_dp,0.6489654712546573_dp,-0.6852363130542333_dp,0.6852363130542333_dp,&
                                     -0.7198818501716109_dp,0.7198818501716109_dp,-0.7528199072605319_dp,0.7528199072605319_dp,&
                                     -0.7839723589433414_dp,0.7839723589433414_dp,-0.8132653151227975_dp,0.8132653151227975_dp,&
                                     -0.8406292962525803_dp,0.8406292962525803_dp,-0.8659993981540928_dp,0.8659993981540928_dp,&
                                     -0.8893154459951141_dp,0.8893154459951141_dp,-0.9105221370785028_dp,0.9105221370785028_dp,&
                                     -0.9295691721319396_dp,0.9295691721319396_dp,-0.9464113748584028_dp,0.9464113748584028_dp,&
                                     -0.9610087996520538_dp,0.9610087996520538_dp,-0.9733268277899110_dp,0.9733268277899110_dp,&
                                     -0.9833362538846260_dp,0.9833362538846260_dp,-0.9910133714767443_dp,0.9910133714767443_dp,&
                                     -0.9963401167719553_dp,0.9963401167719553_dp,-0.9993050417357722_dp,0.9993050417357722_dp /)

    IF(edify_model=='Zero') THEN
        int_edif_f=0.0_dp
        int_edif_dfdy=0.0_dp
        return
    END IF

    ! Predefine bed_tmp, ys, and ustar, including boundary values
    bed_tmp(1:a) = bed
    bed_tmp(0)   = bedl
    bed_tmp(a+1) = bedu

    ys_tmp(1:a) = ys
    ys_tmp(0)   = ysl
    ys_tmp(a+1) = ysu

    ustar_tmp(1:a) = ustar
    ustar_tmp(0)   = max(0.0_dp, 2.0_dp*ustar(1)-ustar(2)) !ustar(1)
    ustar_tmp(a+1) = max(0.0_dp, 2.0_dp*ustar(a)-ustar(a-1)) !ustar(a)
    
    aref_tmp(1:a) = a_ref(1:a)
    aref_tmp(0)   = max(0.0_dp, 2.0_dp*a_ref(1)-a_ref(2)) !a_ref(1)
    aref_tmp(a+1) = max(0.0_dp, 2.0_dp*a_ref(a)-a_ref(a-1)) !a_ref(a)

    ! Loop through (0.5, 1.5, ... a+0.5) to compute integrals
    DO i=1,a+1

        ! Define depth, ustar, bed, aref, epsz, at i-1/2
        !d = 0.5_dp*(max( (water - bed_tmp(i)), 0.0_dp) + max( (water - bed_tmp(i-1)), 0.0_dp))
        us= 0.5_dp*(ustar_tmp(i)+ustar_tmp(i-1))
        bedh = 0.5_dp*(bed_tmp(i)+bed_tmp(i-1))
        arefh = 0.5_dp*(aref_tmp(i)+aref_tmp(i-1))
        d = max(water - bedh, 0.0_dp)
        ! Define daref/dy and dbed/dy, dustar/dy at i-1/2
        dyinv = 1.0_dp/(ys_tmp(i)-ys_tmp(i-1)) 
        daref_dy = (aref_tmp(i) - aref_tmp(i-1))*dyinv
        dbed_dy = (bed_tmp(i) - bed_tmp(i-1))*dyinv
        dus_dy = (ustar_tmp(i) - ustar_tmp(i-1))*dyinv 
       
        ! I experienced numerical problems if the derivative of aref changes
        ! sign -- try setting the derivative to zero in this case (which is
        ! reasonable)
        ! FIXME: Check if this is still needed.
        IF((i>2).and.(i<a)) THEN
            IF( ((aref_tmp(i)-aref_tmp(i-1))*(aref_tmp(i+1)-aref_tmp(i))<0.0_dp)&
              .or.((aref_tmp(i)-aref_tmp(i-1))*(aref_tmp(i-1)-aref_tmp(i-2))<0.0_dp) ) THEN
                daref_dy = 0.0_dp                
            END IF
        END IF

        ! Create vertical suspended sediment profile
        SELECT CASE(sus_vert_prof) 
            CASE('exp') 
                PRINT*, "ERROR: sus_vert_prof = 'exp' is not really supported now."
                PRINT*, "The routine has been coded with 'Rouse' in mind"
                ! Vertical eddy diffusivity
                eps_z =0.5_dp*( 0.1_dp*ustar_tmp(i)*max(water-bed_tmp(i),0.0_dp) + &
                                0.1_dp*ustar_tmp(i-1)*max(water-bed_tmp(i-1),0.0_dp)) 
                eps_z = max(eps_z,1.0e-08_dp) 

                depsz_dy = ( 0.1_dp*ustar_tmp(i)*max(water-bed_tmp(i),0.0_dp) - &
                             0.1_dp*ustar_tmp(i-1)*max(water-bed_tmp(i-1),0.0_dp) &
                             )/(ys_tmp(i)-ys_tmp(i-1))

                !z_tmp = elevation above bed = at 0.5, 1.5, ... 99.5 * depth/no_subints.0 
                !z_tmp = bedh + (d/((no_subints-1)*1.0_dp))*( (/ (j,j=0,no_subints-1) /))
                z_tmp = bedh+arefh+ (d-arefh)/(2.0_dp)*(gauss_abscissae +1.0_dp)

                ! Exponential vertical suspended sediment profile
                f = exp(-wset/eps_z*(z_tmp-bedh))
                
                ! df/dy = df/deps_z * deps_z/dy + df/dbed*dbed/dy
                df_dy = wset*(z_tmp-bedh)/eps_z**2*f*depsz_dy + &
                                wset/eps_z*f*dbed_dy

            CASE('Rouse')
                 
                IF((d>arefh).and.(wset/us<1.0e+3_dp)) THEN  !((d*0.3_dp>arefh)) THEN !.and.(us>wset)) THEN

                    !z_tmp = elevation above bed at gauss integration points 
                    z_tmp = bedh+arefh+ (d-arefh)/(2.0_dp)*(gauss_abscissae +1.0_dp)
                    ! Useful shorthand variables, which save computation
                    ! This routine is computationally demanding, so it is worth
                    ! making some effort.
                    z2bed_inv = 1.0_dp/(z_tmp-bedh) ! Inverse of above, reuse below
                    z2ratio = d*z2bed_inv ! A ratio that comes up a lot
                    rouseno = wset/(0.4_dp*us) ! Rouse number
                    arefh_inv = 1.0_dp/arefh
                    d_on_aref_les1_inv = 1.0_dp/(d/arefh -1.0_dp) ! Useful term
                    
                    ! Calculate vertical profile of suspended sediment
                    f = ((z2ratio-1.0_dp)*(d_on_aref_les1_inv))**rouseno
                    !print*, maxval(f)
                    IF(maxval(f)>1.0_dp) THEN
                        print*, 'maxval(f)>1', maxval(f), d, maxloc(f), z_tmp(maxloc(f)), arefh+bedh
                        stop
                    END IF
                    ! Calculate derivative of f. Use this approach:
                    ! df_dy = df/dbedh*dbedh/dy + df/aref*daref/dy + df/dus*dus/dy
                    ! because all the d/dy terms have only a single value at
                    ! i-1/2. 

                    ! Step1: df/dbedh, evaluated using maxima (symbolic algebra) 
                    ! -- see code in the file lat_flux.max
                    ! [2.5E+0*wset*((Y-h)/aref-1)*(((Y-h)/(z-h)-1)/((Y-h)/aref-1))**(2.5
                    !     1   E+0*wset/ustar)*(((Y-h)/(z-h)-1)/(aref*((Y-h)/aref-1)**2)+((Y-h
                    !     2   )/(z-h)**2-1/(z-h))/((Y-h)/aref-1))/(ustar*((Y-h)/(z-h)-1))]
                    df_dbedh =rouseno*f*(d_on_aref_les1_inv*arefh_inv+z2bed_inv)

                    ! Step2: df/darefh, evaluated using maxima (symbolic algebra) 
                    ! -- see code in the file lat_flux.max
                    ! [2.5E+0*wset*(Y-h)*(((Y-h)/(z-h)-1)/((Y-h)/aref-1))**(2.5E+0*wset/
                    !     1   ustar)/(aref**2*ustar*((Y-h)/aref-1))]
                    df_darefh = rouseno*d*f*d_on_aref_les1_inv*arefh_inv**2

                    ! Step3: df/dus, evaluated using maxima (symbolic algebra) 
                    ! -- see code in the file lat_flux.max
                    ! [-2.5E+0*wset*(((Y-h)/(z-h)-1)/((Y-h)/aref-1))**(2.5E+0*wset/ustar
                    !     1   )*log(((Y-h)/(z-h)-1)/((Y-h)/aref-1))/ustar**2]
                    !df_dus = -(wset/0.4_dp)*f*log((z2ratio-1.0_dp)/(d_on_aref_les1))/us**2
                    !df_dus = -(f/us)*log(f)
                    ! Note f*log(f) --> 0 as f--> 0
                    !df_dus(1:(no_subints-1)) = -(f(1:(no_subints-1))/us)*log(f(1:(no_subints-1)))
                    !df_dus(no_subints)=0.0_dp
                    !df_dus = -(f/us)*log(f)

                    ! Treat special case where f=0, which can go singular. But,
                    ! f*log(f) --> 0 from below as f--> 0 
                    DO j=1,64                    
                        IF(f(j)>1.0e-32) THEN
                            df_dus(j) = -(f(j)/us)*log(f(j))
                        ELSE
                            df_dus(j)=0.0_dp
                        END IF
                    END DO 

                    ! Step4: df_dy = df/dbedh*dbedh/dy + df/aref*daref/dy + df/dus*dus/dy
                    df_dy = df_dbedh*dbed_dy + df_darefh*daref_dy + df_dus*dus_dy

                    ! Note that (analytically), df_dy =0 for z<a_ref+bedh,
                    ! because f=1 there in this zone. If we perturb y by a sufficiently
                    ! small amount below z=a_ref+bedh, then f will also be 1 in
                    ! the neighbouring point. Hence df_dy =0 there
                ELSE
                    !PRINT*, 'ERROR - d< aref in suspended_xsect (int_edify_f)', &
                    !        d, arefh, aref_tmp(i), aref_tmp(i-1), bed_tmp(i), bed_tmp(i-1)
                    !stop
                    !z_tmp = elevation above bed = at 0.5, 1.5, ... 99.5 * depth/no_subints.0,
                    z_tmp = bedh+arefh+ (d-arefh)/(2.0_dp)*(gauss_abscissae +1.0_dp)
                    
                    ! In these shallow waters, define things so there is no
                    ! lateral flux of suspended load -- hmm, actually, not such
                    ! a good idea?
                    f = 0.0_dp !1.0_dp ! Uniform suspended load in very shallow water?
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
                    edify = 0.24_dp*d*us ! Constant
                ELSE
                    edify = 0.0_dp
                END IF

            CASE('Parabolic')
                IF(d>0.0_dp) THEN
                    ! Note this is parabolic, with peak lateral eddy viscosity = 0.4_dp*us*d
                    !edify = (0.4_dp*us*d)*[ (z_tmp-bedh)*(water-z_tmp)/(0.25_dp*d**2)] ! Parabolic
                    edify = (1.6_dp*us/d)*(z_tmp-bedh)*(water-z_tmp) ! Parabolic
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
        !dz = (z_tmp(no_subints)-z_tmp(1))/(1.0_dp*(no_subints-1))
        !int_edif_f(i) = sum(edify*f)*dz
        !int_edif_f(i) = trapz_integrate(no_subints, dz, edify*f)
        !int_edif_f(i) = simpson_integrate(no_subints, dz, edify*f)
        !int_edif_f(i) = wedint(no_subints, dz, edify*f)
        !int_edif_f(i) = newtcotes7(no_subints, dz, edify*f)
        int_edif_f(i) = sum(gauss_weights*edify*f)*(d-arefh)/2.0_dp ! gaussian quadrature

        IF( (sus_vert_prof=='Rouse').and.(Qbedon.eqv..FALSE.)) THEN
            ! Try adding in near-bed portion [a_ref >= z >= bed]. Note that here f =
            ! 1, while the eddy viscosity profile is still parabolic. Integrating
            ! this eddy viscosity profile from [zero , a_ref] gives us the extra constant to add
            int_edif_f(i) = int_edif_f(i) + ( 1.6_dp*us/d*( (arefh)**3)/3)
        END IF
        
        IF(isnan(int_edif_f(i))) THEN
            print*, 'Error: int_edif_f(',i, ') is nan, '  
            stop
        END IF
     
        !Integral (edify*df/dy) dz
        !int_edif_dfdy(i) = sum(edify*df_dy)*dz
        !int_edif_dfdy(i) = trapz_integrate(no_subints, dz, edify*df_dy)
        !int_edif_dfdy(i) = simpson_integrate(no_subints, dz, edify*df_dy)
        !int_edif_dfdy(i) = wedint(no_subints, dz, edify*df_dy)
        !int_edif_dfdy(i) = newtcotes7(no_subints, dz, edify*df_dy)
        int_edif_dfdy(i) = sum(gauss_weights*edify*df_dy)*(d-arefh)/2.0_dp !Gaussian quadrature   
        ! Note -- since df_dy is zero below a_refh, we do not need to add
        ! anything here to cover the region [0, arefh]
        IF(isnan(int_edif_dfdy(i))) THEN
            print*, 'Error: int_edif_dfdy(',i, ') is nan, ', us,'$$$$', f,'$$$$', df_dus
            stop
        END IF
        !print*,'#########', i, int_edif_f(i), int_edif_dfdy(i) 
        !print*, edify, df_dy
    END DO

    ! EXPERIMENTAL LIMITING 
    !DO i=2,a
    !    tmp1(i) = 0.5_dp*maxmod( (int_edif_f(i)+int_edif_f(i-1)), (int_edif_f(i)+int_edif_f(i+1)) )
    !    tmp2(i) = 0.5_dp*maxmod( (int_edif_dfdy(i)+int_edif_dfdy(i-1)), (int_edif_dfdy(i)+int_edif_dfdy(i+1)) )
    !END DO 
    !    tmp1(1) = 0.5_dp*(int_edif_f(1) + int_edif_f(2))
    !    tmp1(a+1) = 0.5_dp*(int_edif_f(a)+int_edif_f(a+1))
    !    tmp2(1) = 0.5_dp*(int_edif_dfdy(1) + int_edif_dfdy(2))
    !    tmp2(a+1) = 0.5_dp*(int_edif_dfdy(a) + int_edif_dfdy(a+1))

    !int_edif_f=tmp1
    !int_edif_dfdy=tmp2
    ! END EXPERIMENTAL LIMITING

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
