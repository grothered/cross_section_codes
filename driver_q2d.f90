Program driver
! Main program to run the quasi-2D model.

use global_defs ! Various important parameters
!use util, only: read_real_table ! Various utility routines
use util_various, only: set_geo, reset_ys, meanvars, compute_slope, compute_critical_shear, Ddx_3E, &
                        interp3, active_zone, conservation_tests1, conservation_tests2, read_real_table
use hydro_xsect, only: calc_friction, calc_shear ! Hydrodynamic and grain friction and shear
use bed_xsect, only: calc_resus_bedload, update_bed ! Compute rates of sediment transport, and evolve bed
use st_venant_solver, only: hyupdate !The longitudinal hydrodynamic solver
use sus, only: susconc_up35, susconc2dup3rdV2  !The longitudinal suspended sediment solver

IMPLICIT  NONE

!!Define variables - This is quite disorganised, but cleaning this is not a priority.
INTEGER:: a, b, i,j, m, p, o,d2,d1, ii, n, i1, jj1, dlength, jmax, writfreq,writfreq1, jj, writcount
REAL(dp):: longt,longt1, cfl1, tmp1
REAL(dp):: bed,bed_old,bed_Vold,bed_oldrefit, ys,ys_oldrefit, mxdeps, fs,fs_g,a_ref, ws, waters, rough_coef 
REAL(dp):: Width, Area, bottom, bottom_old,wid_old, Q,Q2, Q2_old,& 
           dBdh,inuc,NN,taus,taus_g, Qsav, A2, &
           waters_avg,waters_avg_old, waters_old, tsav,C_old, &
           Area_old, U2_old,acUdlast, NN_old, NN_old2, taus_old2, & 
           NN_last, bedl, bedu, ysl, ysu, diff1D, diff1D_old, Q_old, Q2H, Q2H_dT, useme1,useme2, & 
           wset_tmp, Q2_geo, recrd
REAL(dp)::delT,delX, wset, wetwidth, wetwidth_old, DT_old, t, DT, tlast, U2, vels,vels_old, & 
            velslast,taucrit_dep, Qe,Qe_old, Qbed, QbedI, dQbedI, dQbedIp, filler, dqbeddx
REAL(dp)::  R, E,E_old, D, C, q1, rmu,bt,el,x,w,slopes, wt, wt2, vegdrag, taucrit, Cdist, Cdist_old,& 
            Cdist_in, taucrit_dep_ys, dst, qb_G
INTEGER:: l, u,k,kk, incount, count2, seabuf, LF, layers, bedwrite,&
          remeshfreq,morbl,morbu,morbl_old,morbu_old, iost, too_steep
REAL(dp):: Q1in, Vol1, QS2in, VolS2, Source2, pars_out, xxx, visc_bedp, visc_bedm, visc_bed 
REAL(dp):: hlim , Qb, tr, mor,mor1,  mu, erconst, multa, aa,bb, cc, lifttodrag, & 
    rho, mthdta, z0, rhos, burnin, lfkick ,&
    voidf, dsand, d50, g, kvis,  lambdacon, alpha, cfl,man_nveg, Cmouth,& 
    Criver, water_m, water_mthick, veg_ht, &
    v1coef,v4coef, eddis1D,lincrem
LOGICAL:: susdist=.false., sus2d, LAKE, mouthread, norm, vertical, tbston, normmov, readin, & 
    remesh, Qbedon, talmon, susQbal=.false., manning, printall, taucrit_slope_reduction=.false.
CHARACTER(char_len):: boundary_downstream_file, friction_type, grain_friction_type, resus_type, &
                      bedload_type

!Variables that are read in from the inputdata file
NAMELIST /inputdata2/ a, b, jmax, writfreq, t,longt, delX, wset, seabuf, LF,  hlim, &
     Qb, tr, mor, mu, erconst,lifttodrag, rho, rhos, burnin, lfkick,sus2d, LAKE, mouthread, & 
    voidf, dsand, d50, g, kvis, norm, vertical, lambdacon, tbston, alpha, readin, cfl, &
     rough_coef, man_nveg, Cmouth, Criver, layers, bedwrite, remesh, remeshfreq, normmov,& 
     water_m, water_mthick, veg_ht, Qbedon, talmon, manning, v1coef,v4coef,eddis1D,lincrem, &
     boundary_downstream_file, friction_type, grain_friction_type, resus_type, bedload_type

ALLOCATABLE bed(:,:),bed_old(:,:),bed_Vold(:,:),bed_oldrefit(:,:), ys(:,:),ys_oldrefit(:,:),& 
            mxdeps(:,:), fs(:,:), fs_g(:,:),a_ref(:,:), waters(:), waters_old(:),&
            l(:), u(:), morbl(:),morbu(:), recrd(:), & 
            morbl_old(:),morbu_old(:), dbdh(:,:),rmu(:), slopes(:,:), Qsav(:), A2(:),& 
            waters_avg(:),waters_avg_old(:), taucrit_dep(:,:,:),dst(:,:,:), wt(:),wt2(:),&
            vegdrag(:,:), taucrit(:,:,:), wset_tmp(:),& 
            taucrit_dep_ys(:,:), mthdta(:,:), C_old(:), Area_old(:),&
            U2_old(:), wetwidth(:),wetwidth_old(:),& 
            QbedI(:), dQbedI(:), dqbeddx(:,:),dQbedIp(:), qb_G(:,:), &
            filler(:), diff1D(:), diff1D_old(:), Q_old(:), Q1in(:),&
            Vol1(:),pars_out(:),  Width(:),Area(:),bottom(:), bottom_old(:),wid_old(:),&
            Q(:),Q2(:),Q2_old(:), R(:),E(:),E_old(:), D(:),C(:),U2(:),inuc(:),NN(:,:), &
            taus(:,:), taus_g(:,:), vels(:,:), vels_old(:,:), acUdlast(:,:), velslast(:,:), &
            NN_old(:,:), NN_old2(:,:),taus_old2(:,:), NN_last(:,:), Q2H(:), Q2H_dT(:), Q2_geo(:), &
             x(:), Cdist(:,:), Cdist_old(:,:), Cdist_in(:,:), Qe(:,:), Qe_old(:,:), Qbed(:,:), ws(:), & 
            bedl(:), bedu(:), ysl(:), ysu(:), visc_bedp(:), visc_bedm(:), visc_bed(:,:), too_steep(:)

!Read the input parameters
!open(1001, file='inputdata2.modin')
!read(1001,nml=inputdata2)
!close(1001)
READ(*, nml=inputdata2)

PRINT inputdata2

ALLOCATE( bed(a,b),bed_old(a,b),bed_Vold(a,b),bed_oldrefit(a,b), ys(a,b),ys_oldrefit(a,b), & 
          mxdeps(a,b),fs(a,b),fs_g(a,b),a_ref(a,b), ws(b), waters(b), waters_old(b),l(b), u(b), &
          morbl(b),morbu(b) ,morbl_old(b),& 
          morbu_old(b) ,taus(a,b), taus_g(a,b), vels(a,b),vels_old(a,b), Cdist(a,b),Cdist_old(a,b),Cdist_in(a,b), & 
          C_old(b), Area_old(b), U2_old(b), Qe(a,b),Qe_old(a,b), Qbed(a,b), QbedI(b), dQbedI(b),dQbedIp(0:b),& 
          filler(-1:b+2), dqbeddx(a,b) , acUdlast(a,b),velslast(a,b), diff1D(b), diff1D_old(b),  &
          Width(b), Area(b), bottom(b),bottom_old(b),wid_old(b), Q(b), Q2(b), Q2_old(b), R(a), D(b),&
          E(b),E_old(b), C(b), U2(b),dbdh(b,2) ,rmu(b), inuc(b), NN(a,b),NN_old(a,b), NN_old2(a,b),& 
          taus_old2(a,b), NN_last(a,b), x(a),slopes(a,b),wetwidth(b), wetwidth_old(b), Q2H(0:b),Q2H_dT(0:b), Q2_geo(b),& 
          wset_tmp(b), A2(b), Qsav(b), waters_avg(b),waters_avg_old(b), taucrit_dep(a,b,layers), dst(a,b,0:layers+1), &
          wt(b), wt2(b), vegdrag(a,b), qb_G(0:a+1,b), recrd(0:a), & 
          taucrit(a,b,0:layers), taucrit_dep_ys(a,b) , bedl(b), bedu(b), ysl(b), ysu(b), Q_old(b), & 
          visc_bedp(a), visc_bedm(a), visc_bed(a,b), Q1in(1), Vol1(1),pars_out(10), too_steep(a) ) 

!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!BEGIN THE ROUTINE

!!Read in the mouth boundary condition -- or if mouthread=.false., it should be in an equation in the
!tidalmod_fg31.f95 subroutine 'mouth_height'
IF(mouthread) THEN
    call read_real_table(boundary_downstream_file, mthdta ,dlength,2)   
ELSE
    print*, "Analytical water level mouth boundary needed"
END IF
!!!!!!!!!!!!!!!!!!!

!!Initialise the geometry
CALL set_geo(bed,ys,waters,mxdeps,fs,a,b, hlim, readin, water_m,water_mthick) 
! FIXME: Temporarily initialise fs_g and a_ref here -- later, lump with fs
fs_g=fs*0.1_dp ! Initialise Grain friction factor
a_ref=0.01_dp ! Initialise Reference height for suspended load

!!Calculate the averages needed for the 1d hydro
CALL meanvars(bed,ys,waters,fs,a,b,u,l,Width, Area, bottom,dbdh, .false., hlim ) 

!Predefine starting values for a whole host of variables. ANYTHING IMPORTANT IS REDEFINED DEEPER IN THE CODE
delT=0._dp
DT=0._dp
cfl1=cfl
Q=Qb+0._dp*Area
E=0._dp+0._dp*Area
D=0._dp+0._dp*Area
inuc=0._dp
C=0._dp+ 0._dp*Q
Cdist=0._dp
Cdist_old=Criver
too_steep=1

taucrit=1._dp
taucrit_dep=-9.99E+08_dp
taucrit_dep_ys=ys
dst(:,:,0)=0._dp
DO i=1,layers
    dst(:,:,i)=bed(:,:)-i*lincrem
END DO
dst(:,:,layers+1)=9.9E+10_dp

taus=0._dp
taus_g=0._dp
Area_old=Area
C_old=Criver
U2_old=0._dp+0._dp*Q
E_old=0._dp
waters_old=-99999._dp !waters
slopes=0._dp
acUdlast=0._dp
NN_old=0._dp
vels=0._dp
vels_old=0._dp
velslast=0._dp
diff1D=0._dp
diff1D_old=0._dp
longt1=longt
Q2_old=0._dp
Q2H=0._dp
Q2H_dT=0._dp
Q2_geo=0._dp
wetwidth=0._dp
wetwidth_old=0._dp
writfreq1=9E+08
Q1in=0._dp
Vol1=0._dp
QS2in=0._dp 
VolS2=0._dp
Source2=0._dp
bed_old=bed
bed_Vold=bed
bed_oldrefit=bed
ys_oldrefit=ys
morbl=0
morbu=0
morbl_old=0
morbu_old=0

writcount=9E+8

!!!If we have a 'seabuf' (a region of the model where no sediment or morphological changes happen - i.e. a hydrodynamic buffer zone), then fix the suspended sediment concentration in this zone
IF(seabuf>0) THEN
	C(1:seabuf)=Cmouth
	Cdist(:,1:seabuf)=Cmouth
END IF

!!The roughness multiplier
rmu= fs(floor(a/2._dp),:)/(8._dp*g) !Roughness multiplier, a preliminary value. This used subsequently as the 'real' (properly integrated) roughness. So it is updated as velocity profiles become available

!Some more random definitions
NN=bed*0._dp !At the moment this is unused, but it flows through all the routines. Very inefficient



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Define files to write output to. 
!   The file extensions have the following info:
! .1DO = Longitudinal 1D variable (e.g. Discharge)
! .2DO = 2D variable (e.g. the Bed array)
! .0DO = 0D variable (e.g. the Time)
OPEN(1,file="Area.1DO")
OPEN(2,file="Discharge.1DO")
OPEN(42,file='Discharge_halftime_lim.1DO')
OPEN(23,file="Area_halftime.1DO")
OPEN(4, file="Width.1DO")
OPEN(5,file="Bottom.1DO")
OPEN(7,file="Water.1DO")
OPEN(8,file="Susconc.1DO")
OPEN(9,file="Times.ODO")
OPEN(10,file="Resuspension.1DO")
OPEN(11, file="BedloadI.1DO")
OPEN(15,file="Vels.1DO")
OPEN(24,file="dWidth_dWater.1DO")
OPEN(12,file="Taus.2DO")
OPEN(13,file="Ys.2DO")
OPEN(3,file="Bed.2DO")
OPEN(33,file="Susconc.2DO")
OPEN(25,file="dropout.0DO")
OPEN(26,file="dropout4")
OPEN(35,file='Randomstuff')
!OPEN(14,file="taucrit_dep")
!OPEN(15,file="tausss")
!OPEN(22,file="bndry")
!OPEN(34,file="Nfile.2DO")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!MAIN LOOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO j= 1, jmax

    !Here we define 3334, so that we can come back to the start of an iteration if
    !needed -- this is done to reduce the time step after the rate of morphological
    !evolution becomes too rapid.
    3334 CONTINUE        
    !!Define variables that will serve as the previous time step variables
    tlast=t
    DT_old=DT
    C_old=C
    Area_old=Area
    U2_old=U2
    E_old=E
    diff1D_old=diff1D
    Q2_old=Q2
    Cdist_old=Cdist
    vels_old=vels
    waters_old=waters
    waters_avg=waters
    Q_old=Q
    Qe_old=Qe
    bottom_old=bottom
    wid_old=Width
    wetwidth_old=wetwidth
    bed_old=bed
    waters_avg_old=waters_avg

    !Print out info so we can track the progress
    IF(mod(j-1,writfreq*5).eq.0) THEN
                ! PRINT*,'status', t/3600., j, longt1, delT, cfl1, cfl!, "hours"!, C(1:5)
                PRINT*, '#################'
                PRINT*, '# Output Step:', j
                PRINT*, '#   Time (hrs):', t/3600._dp 
                PRINT*, '#   Long timestep: ', longt1, ' Hydrodynamic timestep: ', delT
                PRINT*, '#   cfl_local: ', cfl1, 'cfl: ', cfl 
                PRINT*, '################'
    END IF

    !!!Pre-set some variables
    incount= 0 !A convenient flag to determine how often we recalculate the geometry within the small hydrodynamic time step
    Q2=0._dp !Q2 is the half time step discharge.
    Q2H=0._dp !Q2H is the half time step discharge, centred at i+1/2
    Q2H_dT=0._dp !Q2H_dT is the time-integrated half time step discharge, centred at i+1/2


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!LONGITUDINAL HYDRODYNAMIC TIME STEP
    ! FROM tlast TO tlast+DT
    !(actually DT will be defined afterwoulds, but is approximately=
    !longt1)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO WHILE ((t-tlast<longt1)) 

        incount=incount+1 !update the flag for counting how many hydrodynamic time steps I do in this loop
        m=b !m can be used to include less cross-sections in the x direction (m<b), but use m=b for now

    !!Recalculate the mean variables, if needed
        IF(incount>1) THEN
            CALL meanvars(bed(:,1:m),ys(:,1:m),waters(1:m),fs(:,1:m),a,m,u(1:m),l(1:m),Width(1:m), & 
            Area(1:m), bottom(1:m), dbdh(1:m,1:2), .false., hlim ) 
        END IF

            !Normal MacCormack
            CALL hyupdate(delT,delX,bottom(1:m),Width(1:m),Area(1:m),Q(1:m), Q2H_dT(0:m),&
                          waters(1:m),t,m,j,dbdh(1:m,1:2), &
                          rmu(1:m),inuc(1:m),LAKE, hlim, Qb,tr, mthdta,mouthread, &
                          dlength, bottom(1:m), rho, g,cfl1,v1coef,v4coef,seabuf)

    END DO 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!END OF LONGITUDINAL HYDRODYNAMIC LOOP 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!The time step
    DT= t-tlast  
    !print*, fs(5,300)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!BEGIN THE CROSS SECTIONAL HYDRODYNAMICS/MORPHODYNAMICS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Check whether we should be updating the morphology
    IF(t<3600._dp*burnin) THEN  !This allows for a hydrodynamic 'burn-in' where we don't update the morphology - useful to stabilise the hydrodynamics prior to morphodynamic updating 
        mor1=1.0E-16_dp
    ELSE
        mor1=mor
    END IF

    !!Evaluate half time step variables
    waters_avg=0.5_dp*(waters+waters_old) !The half time step waters
    m=b
    !More half time-step variables
    CALL meanvars(bed(:,1:m),ys(:,1:m),waters_avg(1:m),fs(:,1:m),a,m,u(1:m),l(1:m),Width(1:m), & 
            A2(1:m), bottom(1:m),dbdh(1:m,1:2),.false.,hlim) !Figure out the average geometry
    Q2H=Q2H_dT/DT !The conservative 'mean' discharge throughout the last time step, evaluated spatially at 1/2, 3/2, 5/2, ... i.e. halfway between the cross-sections
    Q2H(0)= delX/DT*(A2(1)-Area_old(1))*2._dp + Q2H(1) !Note the use of (A2(1)-Area_old(1))*2 here--This should be Area(1)-Area_old(1) - except that Area(1) is not correct, because it has not yet adjusted for the change in Y(1) due to the boundary condition. However, A2(1) has adjusted, and should be halfway between Area_old and Area. So this works quite well.
    Q2=0.5_dp*(Q2H(1:b)+Q2H(0:(b-1))) !Q2H, evaluated spatially at 1, 2, 3, 4, ... i.e. at the cross-sections
    
    !!Define 'limited' version of Q2
    Q2_geo=min(abs(Q2H(1:b)), abs(Q2H(0:b-1)), 0.5_dp*(abs(Q(1:b))+abs(Q_old(1:b)) )) & 
            *0.5_dp*(sign(1._dp,Q2H(1:b))+sign(1._dp,Q2H(0:b-1)) ) !sqrt(abs(Q2H(1:b))*abs(Q2H(0:(b-1))))*0.5_dp*(sign(1._dp,Q2H(1:b))+sign(1._dp,Q2H(0:b-1)) )
    ! The 'limited' version is only used if it is sufficiently different to the
    ! ordinary Q2
    DO i=1,b
        IF(abs(Q2_geo(i)-Q2(i))/(min(abs(Q2(i)), abs(Q2_geo(i)))+1.0E-07_dp)<= .3_dp) THEN
            !!!Note - not sure if this is useful - but the difference between [
            !Q2^{k+1/2}_{j+1/2} + Q2^{k+1/2}_{j-1/2} ] and [ Q^{k+1}_{j} +
            !Q^{k}_{j}] (which are both reasonable estimates of the half time
            !step Q at j) is equal to 1/4( del^2 Q^{k}_{j} + (Qpred^{k+1}_{j}  -
            !Qcor^{k+1}_{j})). Perhaps it could be useful to design a method
            !that detects when there are large differences in the velocities and
            !adjusts accordingly. I've no idea how though!
            Q2_geo(i)=Q2(i)
        END IF
    END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Loop over every cross section, calculate shear and rates of erosion. We
    ! can use this to do sus sed and bedload, and after that we can update the
    ! morphology.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO i =(seabuf+1), b
  
        IF(l(i)>0) THEN !This just checks that the cross-section is wet.
            call calc_friction(friction_type, grain_friction_type,rough_coef, waters_avg(i),&
                               u(i)-l(i)+1, bed(l(i):u(i),i), vels(l(i):u(i),i), man_nveg, &
                               d50,veg_ht, rhos, rho, g, fs(l(i):u(i), i),&
                               vegdrag(l(i):u(i), i) ,fs_g(l(i):u(i),i), dsand, j, a_ref(l(i):u(i),i))           
 
          
            IF((l(i)>1).and.(u(i)<a)) THEN
            ! Use central slope estimate at l(i) and u(i)
                call compute_slope(u(i)+1-(l(i)-1)+1,slopes(l(i)-1:u(i)+1,i), bed(l(i)-1:u(i)+1,i), ys(l(i)-1:u(i)+1,i))
            ELSE
            ! Use inner slope estimate at l(i) and u(i)
                call compute_slope(u(i)-l(i)+1,slopes(l(i):u(i),i), bed(l(i):u(i),i), ys(l(i):u(i),i))
            END IF    
            
            ! Figure out the value for heights and ys at the wetted edge of the
            ! cross-section. This is useful for the bed solver
            IF(l(i)>1) THEN
                ysl(i)=ys(l(i)-1, i)
                bedl(i)=bed(l(i)-1, i)
            ELSE
                ysl(i)=ys(l(i), i)- 0.001_dp !2._dp*ys(l)-ys(l+1)
                bedl(i)=waters_avg(i) !2._dp*heights(l)-heights(l+1)
            END IF
            IF(u(i)<a) THEN
                ysu(i)=ys(u(i)+1, i)
                bedu(i)=bed(u(i)+1, i)
            ELSE
                ysu(i)=ys(u(i), i)+0.001_dp!2._dp*ys(u)-ys(u-1)
                bedu(i)=waters_avg(i) !2._dp*heights(u)-heights(u-1)
            END IF
            !!!!!

            !!!Calculate critical shear stress

            call compute_critical_shear(u(i)-l(i)+1, layers, bed(l(i):u(i),i),&
                                         slopes(l(i):u(i),i), taucrit_dep(l(i):u(i),i,1:layers),&
                                         taucrit(l(i):u(i),i, 0:layers),&
                                         dst(l(i):u(i),i, 0:layers+1), lincrem, mu, lifttodrag, &
                                         taucrit_slope_reduction, erconst)
            !! ADD WARNING ABOUT BED LAYERS
            DO i1=l(i), u(i)
                IF(dst(i1,i,1).eq.0._dp) THEN
                    print*, 'FIXME: Cutting down into harder bed layers. This is not really &
                    supported at present, you need to think about it / test first.'
                    stop
                END IF
            END DO

            ! Define the lateral variation of the suspended sediment
            ! concentration - this was needed to distinguish the cases in which
            ! we were using a lateral sediment distribution equation 
            IF((.NOT.susdist).AND.(.NOT.sus2d)) Cdist(l(i):u(i),i)=C(i)

            !Predefine some more things
            Qe(:,i)=0._dp
            Qbed(:,i)=0._dp
            qb_G(:,i)=0._dp
            taus(:,i)=0._dp !Predefine hydrodynamic shear stress
            taus_g(:,i)=0._dp ! Predefine the grain shear stress
            E(i)=0._dp  !Total rate of erosion over the cross section. 
            D(i)=0._dp !Total rate of deposition over the cross section. 
        END IF !l(i)>0

    END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!Calculate the shear over every cross-section
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO i=(seabuf+1),b
        
        CALL calc_shear(u(i)-l(i)+1, DT, waters_avg(i), Q2_geo(i), bed(l(i):u(i),i), &
                        ys(l(i):u(i),i),A2(i), bottom(i), fs(l(i):u(i),i),rmu(i),inuc(i), &
                        taus(l(i):u(i),i),NN(l(i):u(i),i), j &
                        ,slopes(l(i):u(i),i), hlim,u(i)-l(i)+1, vegdrag(l(i):u(i),i), rho, & 
                        rhos, voidf, d50, g, kvis, vertical, lambdacon, tbston, ysl(i),ysu(i),bedl(i),bedu(i), & 
                        .FALSE.) 

        ! Compute velocity, and the grain shear stress
        vels(:,i)= sqrt(8._dp*abs(taus(:,i))/(fs(:,i)*rho))*sign(1._dp+0._dp*taus(:,i), taus(:,i))
        taus_g(l(i):u(i),i) = rho*vels(l(i):u(i),i)**2*(fs_g(l(i):u(i),i)/8._dp)*&
                              sign(1._dp+0._dp*taus(l(i):u(i),i), taus(l(i):u(i),i))
    END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
    !Calculate rates of resuspension and bedload transport over each cross-section
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO i=(seabuf+1),b
        call calc_resus_bedload(u(i)-l(i)+1, DT, waters_avg(i), Q2_geo(i), bed(l(i):u(i),i),&
                                      ys(l(i):u(i),i),A2(i), fs(l(i):u(i),i),recrd(l(i):u(i)),&
                                      E(i), Cdist(l(i):u(i),i), wset, 1, taus(l(i):u(i),i), &
                                      taus_g(l(i):u(i),i),vels(l(i):u(i),i), j, slopes(l(i):u(i),i), &
                                      hlim,mor1,taucrit_dep(l(i):u(i),i,1:layers),layers, &
                                      taucrit_dep_ys(l(i):u(i),i), dst(l(i):u(i),i, 0:layers+1), &
                                      taucrit(l(i):u(i),i,0:layers), rho, Qe(l(i):u(i),i), &
                                      Qbed(l(i):u(i),i),qb_G(l(i)-1:u(i)+1,i), rhos, & 
                                      voidf, dsand, d50, g, kvis, norm, alpha, Qbedon,talmon,&
                                      ysl(i),ysu(i),bedl(i),bedu(i), resus_type, bedload_type, &
                                      a_ref(l(i):u(i),i)) 
    END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! SUSPENDED SEDIMENT / BEDLOAD ROUTINE -Choice of 1 or 2d
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(.NOT.sus2d) THEN
    !!Calculate velocity
        DO i = 1,b
            IF(A2(i)>0._dp) THEN
                U2(i)= Q2(i)/A2(i)
                wetwidth(i)= ys(u(i),i)-ys(l(i),i)
            ELSE
                U2(i)=0._dp
                wetwidth(i)=0._dp
            END IF
        END DO
        !Calculate 1D dispersion coeff. 
        diff1D=0._dp+ eddis1D*sqrt(abs(taus(a/2,:))*max(waters_avg-bottom, 0._dp ) )
        !In the seabuf region, set 1D dispersion coef to the value at the channel mouth
        diff1d(1:seabuf)=diff1d(seabuf+1)
        
        !if(j==400) THEN
        !Criver=0.05_dp !INJECT C_RIVER 
        !print*, 'INJECT C_RIVER'
        !END IF
        !Cmouth=E(seabuf+1)*rhos/wset

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Update suspended sediment concentration in 1D.  Do the update from tlast to t
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        wset_tmp(seabuf+1:b)=wset !Settling velocity
        wset_tmp(1:seabuf)=0._dp !No settling in the seabuf region

        ! NOTE: C, Cmouth, Criver etc are in g/L (or kg/m^3), NOT (m^3/m^3)
        CALL susconc_up35(b,DT, Area, Q2, Q2_old, delX,C, U2,E*rhos, E_old*rhos, &
                          D*rhos, Cmouth , C_old, Area_old, &
                          Criver, wset_tmp, wetwidth,wetwidth_old,&
                          diff1D, diff1D_old,pars_out)

        !!BEDLOAD IN THE 1D CASE
        QbedI=0._dp !Cross-sectionally integrated bedload flux
        DO i=seabuf+1,b
            IF(norm) THEN !Include the slope factor in this case
                QbedI(i)=0.5_dp*(sum( Qbed((l(i)+1):(u(i)-1),i)*sqrt(1._dp+slopes((l(i)+1):(u(i)+1),i)**2._dp )*&
                        (ys((l(i)+2):u(i), i)-ys(l(i):(u(i)-2), i))) &
                        +Qbed(l(i),i)*sqrt(1._dp+slopes(l(i),i)**2._dp)*(ys(l(i)+1,i)-ys(l(i),i)) &
                        +Qbed(u(i),i)*sqrt(1._dp+slopes(u(i),i)**2._dp)*(ys(u(i),i)-ys(u(i)-1,i))  )
            ELSE !Don't include the slope factor
                QbedI(i)= 0.5_dp*(sum( Qbed((l(i)+1):(u(i)-1),i)*&
                        (ys((l(i)+2):u(i), i)-ys(l(i):(u(i)-2), i))) &
                        +Qbed(l(i),i)*(ys(l(i)+1,i)-ys(l(i),i)) &
                        +Qbed(u(i),i)*(ys(u(i),i)-ys(u(i)-1,i))  )
            END IF
        END DO
        !Now estimate the bedload derivatives. Filler is just a useful thing to put the boundary conditions in
        filler(1:b)=QbedI(1:b)
        !filler(-1:0)=QbedI(1) !Boundary condition - QbedI(mouth)=QbedI(1)
        filler(0)=0._dp !QbedI(1)-(QbedI(2)-QbedI(1))
        filler(-1)=0._dp !QbedI(1)-2._dp*(QbedI(2)-QbedI(1))
        !filler(b+1:b+2)=QbedI(b)!Boundary condition - QbedI(upstream)=QbedI(b)
        filler(b+1)=0._dp !QbedI(b)+(QbedI(b)-QbedI(b-1))
        filler(b+2)=0._dp !QbedI(b)+2.0_dp*(QbedI(b)-QbedI(b-1))
        !Calculate the third order flux limited bedload derivative
        CALL Ddx_3E(b,filler,delX, dQbedI) 
        !dQbedI=0._dp

        !Now we distribute that bedload derivative over each cross-section
        DO i=1, b
            !Here we include a trick to stop division by zero.
            IF(abs(QbedI(i))>1.0E-12_dp) THEN
                dqbeddx(:,i)= abs(Qbed(:,i))/abs(QbedI(i))*dQbedI(i)
            ELSE
                dqbeddx(:,i)=dQbedI(i)/(ys(u(i),i)-ys(l(i),i))+0._dp*Qbed(:,i)
            END IF
        END DO
        !!!!END SUSPENDED SEDIMENT / BEDLOAD IN THE 1D CASE

    ELSE !Here we have the sediment calculations where the suspended sediment is calculated using a fully 2D model

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Compute suspended sediment in 2d
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! So note that now this is time-stepping from tlast-DT_old/2 to
        ! tlast+DT/2 -- could be a problem if DT changes quickly (or at all?)?
        CALL susconc2dup3rdV2(a,b-seabuf, 0.5_dp*(DT+DT_old),delX,ys(:,seabuf+1:b),waters_avg(seabuf+1:b), &
        waters_avg_old(seabuf+1:b),vels(:,seabuf+1:b),vels_old(:,seabuf+1:b),bed(:,seabuf+1:b),&
         Cdist(:,seabuf+1:b),Cdist_old(:,seabuf+1:b), Qe(:,seabuf+1:b)*rhos,&
        Qe_old(:,seabuf+1:b)*rhos, Cmouth+0._dp*ys(:,1), Criver+0._dp*ys(:,1), wset,fs(:,seabuf+1:b), &
        slopes(:,seabuf+1:b), l(seabuf+1:b), u(seabuf+1:b),lambdacon)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!BEDLOAD -- 2D CASE
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !print*, 'NEED TO CHANGE BEDLOAD TO THIRD ORDER'
        !stop
        dqbeddx(:,2:b-1)= (Qbed(:,3:b)-Qbed(:,1:b-2))/(2._dp*delX)
        dqbeddx(:,1)= (Qbed(:,2)-Qbed(:,1))/delX
        dqbeddx(:,b)= (Qbed(:,b)-Qbed(:,b-1))/delX
    END IF 

    !!!!!!!!!!!!!!!!!!
    !!!!!END SUSPENDED-SEDIMENT / BEDLOAD CALCULATIONS
    !!!!!!!!!!!!!!!!!!
     
    !!!!!!!!!!!!!!!!!!!!
    !!!MORPHOLOGICAL UPDATING
    !!!!!!!!!!!!!!!!!!!
    DO i=seabuf+1,b 
           
        p=l(i) !save space 
        o=u(i) !save space
        x(p:o)=NN(p:o,i)

        !Here we define C, depending on whether we have sus2d, susdist, or constant C
        IF((.NOT.susdist).AND.(.NOT.sus2d)) THEN
            Cdist(:,i)=C(i)
            Cdist_in(:,i)=0.5_dp*(C(i)+Cdist_old(:,i)) !Half time step useage
        ELSE
            Cdist_in(:,i)=Cdist(:,i) !Note that this is automatically half time step usage, since the 2Dsusload routine steps from tlast-DT_old/2 to tlast+DT/2
        END IF
        
        !Morphological evolution routine
        CALL update_bed(o-p+1, DT, waters_avg(i), Q2(i), bed(p:o,i),ys(p:o,i),A2(i), &
                        recrd(p-1:o), E(i),&
                        D(i), Cdist_in(p:o,i) ,1 , taus(p:o,i),taus_g(p:o,i),&
                        j,slopes(p:o,i), hlim,mor1,taucrit_dep(p:o,i,1:layers),layers,&
                        taucrit_dep_ys(p:o,i), o-p+1, taucrit(p:o,i,0:layers), rho, & 
                        Qe(p:o,i), Qbed(p:o,i), qb_G(p-1:o+1,i), wset,dqbeddx(p:o,i),&
                        rhos, voidf, d50, g, Qbedon, normmov,sus2d,ysl(i), & 
                        ysu(i),bedl(i),bedu(i),1, bed(p:o,i), talmon, .false., too_steep)

        ! Test for changes in the bank values, which shouldn't happen with the
        ! new routine
        IF(l(i)>1) THEN
            IF(bed(l(i)-1,i)>bedl(i)) THEN
               print*, 'ERROR: bedl has changed'
               stop
            END IF
        END IF
        IF(u(i)<a) THEN
            IF(bed(u(i)+1,i)>bedu(i)) THEN
               print*, 'ERROR: bedu has changed'
               stop
            END IF    
        END IF

        IF(.TRUE.) THEN
        !   A version of the Delft bank erosion model. 
        !   If erosion is occuring at the channel margins,
        !   then assign it to the neighbouring dry bed point
            IF((bed(l(i),i)<bed_old(l(i),i)).AND.(l(i)>1)) THEN
                    bed(l(i)-1,i) = bed(l(i)-1,i) - (bed_old(l(i),i) - bed(l(i),i))
                    bed(l(i),i) = bed_old(l(i),i)
            END IF
            IF((bed(u(i),i)<bed_old(u(i),i)).AND.(u(i)<a)) THEN
                    bed(u(i)+1,i) = bed(u(i)+1,i) - (bed_old(u(i),i) - bed(u(i),i))
                    bed(u(i),i) = bed_old(u(i),i)
            END IF
        END IF


    END DO 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!END OF MORPHODYNAMIC UPDATE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !NEW DELT-AIM TIME STEP - here we calculate an idea of what we would like the time
    !step to be - basically we don't want the channel movement to be too large, or
    !the suspended sediment CFL number to be exceeded.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF((maxval(abs(bed-bed_old))> 0.0E-03_dp).and.(.FALSE.)) THEN
        !Here we have 2 cases. Firstly, if the rate of evolution is really
        !really fast, then we go back to the last time step with a short longt
        IF(maxval(abs(bed-bed_old))>0.003_dp) THEN
                !Things are evolving too fast! We try to 'undo' this time step, and go
                !back, with a lower longt
                PRINT*, 'undoing a time step', t, delT, maxval(abs(bed-bed_old)), j, longt1, cfl1
                longt1=min(longt1*1.0E-03_dp/maxval(abs(bed-bed_old)), longt1/2._dp) !This is short enough that we will just do one hydrodynamic time step
                cfl1=min(cfl1,0.5_dp*longt1/(delX/( maxval(( g*(waters-bottom))**0.5_dp+abs(Q/Area)))) ) !And this will make sure that even one hydrodynamic time step will not be too much, eventually
                PRINT*, '--and--', longt1, cfl1 
                t=tlast
                DT=DT_old
                C=C_old
                Q=Q_old
                Area=Area_old
                U2=U2_old
                E=E_old
                diff1D=diff1D_old
                Q2=Q2_old
                Cdist=Cdist_old
                vels=vels_old
                waters=waters_old
                Qe=Qe_old
                bottom=bottom_old
                Width=wid_old
                bed=bed_old
                waters_avg=waters_avg_old 
                GOTO 3334  !Go back to the start of the time step
        
        ELSE
                ! In this case, delt-aim is chosen so the erosion is not too
                ! fast, and also the suspended sediment cfl is respected
                ! (too-much, for safety purposes!)
                longt1=min(longt1*1.0E-03_dp/maxval(abs(bed-bed_old)), longt, 0.6_dp*minval(delX/(abs(Q/Area)+0.01_dp)))
                cfl1=min(cfl,0.5_dp*longt1/(delX/(maxval(( g*(waters-bottom))**0.5_dp+abs(Q/Area)))) )
        END IF
    ELSE
        ! In this case, the time step is limited by longt and the susconc CFL
        ! condition
        longt1=min(longt, minval(0.6_dp*delX/(abs(Q/Area)+0.01_dp))) 
        cfl1=min(cfl,1.5_dp*cfl1)
    END IF
    !!!!!!!!!!!!END DELT-AIM

    !!!Conservation_measurements
    !Discharge
    CALL conservation_tests1( b-1,Area(2:b),Area_old(2:b),delX,Q2H(1:b),DT,Q1in(1),Vol1(1) )
    !1D susconc
    IF((.NOT.susdist).AND.(.NOT.sus2d)) THEN
        CALL conservation_tests2(b,Area(1:b),wetwidth(1:b),delX,Q2(1:b),DT,C(1:b),C_old(1:b), & 
            pars_out(1), pars_out(2), E(1:b), wset, rhos, QS2in, VolS2, Source2)
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!Save variables
    IF((mod(j-1, writfreq1).EQ.0).AND.(mor1.EQ.mor)) THEN !.or.((mod(j,2).eq.0).AND.(t>72.*3600.))) THEN
        WRITE(1,*) Area
        WRITE(23,*) A2 !QbedI
        WRITE(24,*) dbdh(:,1)
        WRITE(25,*) Q1in(1), Vol1(1), QS2in, VolS2, Source2
        WRITE(2,*) Q
        WRITE(42,*) Q2_geo
        WRITE(4,*) Width
        WRITE(5,*) bottom
        WRITE(7,*) waters
                IF(sus2d) THEN
                    WRITE(8,*), Cdist(a/2,:)
                ELSE
                    WRITE(8,*) C 
                END IF
        WRITE(9,*) t
        WRITE(10,*) E
        WRITE(11,*) QbedI !dqbeddx(a/2,:)
        WRITE(15,*) sqrt(abs(taus(a/2,:))/(rho*fs(a/2,:)/8._dp))*sign(1._dp,taus(a/2,:))
        !WRITE(22,*) bed(a/2,75:76) !rmu!QbedI
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Save larger '2D' variables less frequently
    IF((mod(j-1,bedwrite).EQ.0).AND.(mor1.EQ.mor)) THEN
        WRITE(3,*) bed
        
        IF(remesh) write(13,*) ys
        
        WRITE(12,*) taus
        !WRITE(34,*) NN
        !write(12,*) vels
        WRITE(33,*) Cdist

        !BURST OF WRITING 
        !if(writfreq1.ne.writfreq) THEN
        writfreq1=writfreq !So we will write every writfreq steps following this, until a time determined by the next 'BURST OF WRITING' code
        writcount=-1
        !end if
    END IF

    !BURST OF WRITING
    !These pieces of code control writing of hydrodynamic output
    writcount=writcount+1
    IF(writcount==1000) THEN
        writfreq1=9E+8 !So we will never write with this value of writfreq.
    END IF

    !Smoothing routine - here we stop the banks from being too steep to
    !erode, by smoothing if the height difference between points is too
    !large. For symmetry preservation purposes, we approach from the banks
    !toward the centre, starting from both the left and right banks. If we
    !smooth at a point on the left side, then we smooth at a corresponding
    !point on the right side, without doing a test to see if we should.
    !Although by symmetry we obviously should, sometimes round-off error can
    !mess with this.
    !It's interesting to note that when using the linear re-interpolation
    !(which is very diffusive), I concluded that I didn't need to use this -
    !whereas now I am using cubic interpolation, and I can see that it's
    !probably important to prevent steep banks which lead to irregular
    !changes in the channel.
    !IF(mod(j,20).eq.-1) THEN
    !    IF(.true.) THEN
    !        DO i=1,b
    !            DO ii=max(5-1,2), floor(a*0.5_dp) !min(u(i),a-1)
    !                !!Approach from left 
    !                w=(bed(ii,i)-bed(ii+1,i))
    !                IF(abs(w)>abs(ys(ii,i)-ys(ii+1,i))) THEN
    !                    PRINT*, 'smoothing ', i
    !                    bed(ii,i)=bed(ii,i)-0.33333_dp*w !-(w-sign(1._dp,w))  !(w-1._dp)
    !                    bed(ii+1,i)=bed(ii+1,i)+0.33333_dp*w !(w-sign(1._dp,w)) !+(w-1._dp)             
    !                    !END IF
    !                    !!Approach from right
    !                    !w=(bed(a-ii+1,i)-bed(a-ii,i))
    !                    !if(abs(w)>1._dp) THEN
    !                    !print*, 'smoothing ', i
    !                    bed(a-ii+1,i)=bed(a-ii+1,i)-0.33333_dp*w !-(w-sign(1._dp,w))  !(w-1._dp)
    !                    bed(a-ii,i)=bed(a-ii,i)+0.33333_dp*w !(w-sign(1._dp,w)) !+(w-1._dp)             
    !                END IF
    !            END DO
    !        END DO
    !    END IF
    !END IF

    !!!!!!!Occasionally we remesh the cross-section    
    IF(remesh) THEN
        IF(sus2d) THEN
                PRINT*, "Remeshing not supported - and note susconc2d needs all the ys(:,i) to be identical"
                STOP
        END IF

        IF(mod(j,remeshfreq).eq.0) THEN
            DO i=seabuf+1, b
                morbl(i)=0
                morbu(i)=0
                CALL active_zone(a-2,bed(2:a-1,i),bed_oldrefit(2:a-1,i),morbl(i),morbu(i), 6) !Fine the zone where heights( is different to sectionsold 
                IF(((morbl(i).ne.0).AND.(morbu(i).ne.0)).AND. & 
                        ((abs(morbl(i)-morbl_old(i))>1).OR.(abs(morbu(i)-morbu_old(i))>1))) THEN
            
                    IF(morbl(i)+morbu(i).ne.a-1) THEN
                        !Sometimes symmetry breaking can happen in the active_zone function. We work against that trend here
                        morbl(i)=min(morbl(i), a-1-morbu(i))
                        morbu(i)=a-1-morbl(i)
                    END IF
            
                    PRINT*, 'Refitting section', i, morbl(i), morbu(i), bed(morbl(i)+1,i), bed(morbu(i)+1,i)

                    CALL reset_ys(a-2,ys(2:a-1,i),morbl(i),morbu(i), 0.1_dp, 5.0_dp)

                    !Find the heights associated with the new ys
                    CALL interp3(ys_oldrefit(:,i),bed(:,i),ys(:,i),a)
                    
                    DO jj=1, layers

                        CALL interp3(ys_oldrefit(:,i),  taucrit_dep(:,i,jj),ys(:,i), a) !Beware the risk that this could cause 'leakage'
                        taucrit_dep(:,i, jj)= min(bed(:,i), taucrit_dep(:,i,jj)) !So if we interpolate the critical shear layer, it should of course not be above the bed level.

                    END DO

                    morbl_old(i)=morbl(i)
                    morbu_old(i)=morbu(i)
                END IF

                ys_oldrefit(:,i)=ys(:,i) 
                bed_oldrefit(:,i)=bed(:,i) 
            END DO
        END IF
    END IF

    !!Recalculate mean variables
    CALL meanvars(bed,ys, waters,fs,a,b,u,l,Width, Area, bottom, dbdh(:,1:2), .false. ,hlim) 

    !Sometimes too rapid accretion can push points out of the water - try this to reduce that
    printall=.false.
    DO i=1,b
        IF(waters(i)<bottom(i)) THEN
            !printall=.true.
            PRINT*, 'Point', i, 'accreted or drained out', C(i), bottom(i), bottom_old(i), & 
                                    waters(i), waters_old(i), minval(bed(l(i):u(i),i))
            PRINT*, '........taus.......'
            PRINT*, taus(l(i):u(i),i)
            PRINT*, '........bed.......'
            PRINT*, bed(l(i):u(i),i)
            PRINT*, '.........Qbed......'
            PRINT*, Qbed(l(i):u(i),i)
            PRINT*, '.........QbedI.....'
            PRINT*, QbedI(i-1:i+1)
            waters(i)=minval(bed(:,i))+hlim/1.0001_dp
            Q(i)=0._dp
        END IF
    END DO

    IF(printall) THEN
        WRITE(1,*) Area !, Area_old
        WRITE(23,*) QbedI
        WRITE(24,*) dbdh(:,1)
        WRITE(25,*) dbdh(:,2)
        WRITE(2,*) Q  !, Q_old
        WRITE(42,*) Q2H
        WRITE(4,*) Width
        WRITE(5,*) bottom !, bottom_old
        WRITE(7,*) waters !, waters_old
                IF(sus2d) THEN
                WRITE(8,*), Cdist(a/2,:)
                ELSE
                WRITE(8,*) C 
                END IF
        WRITE(9,*) t
        WRITE(10,*) E
        WRITE(11,*) D
        WRITE(15,*) taus(a/2,:) !sqrt(abs(taus(a/2,:))/(rho*fs(a/2,:)/8._dp))*sign(1._dp,taus(a/2,:))
        WRITE(22,*) bed(a/2,75:76) !rmu!QbedI
    END IF

    !Trial to improve mass conservation. Notice that the any change between
    !bottom and bottom_old must be due to sedimentation/erosion
    !waters=waters+min(0.02_dp, abs(bottom-bottom_old))*sign(1._dp, bottom-bottom_old) ! However, sometimes when the flow goes overbank, this is badly behaved. Changes fast, need to limit it.
    !waters=waters+min(0.02_dp, abs(bottom-bottom_old))*sign(1._dp, bottom-bottom_old) ! However, sometimes when the flow goes overbank, this is badly behaved. Changes fast, need to limit it.

END DO !j=1, jmax 

PRINT*, "Finished", " time =", t, "j=", j 

END PROGRAM driver



