&inputdata

!!!!!!!!!!!!!!!!!
!Spatial and temporal variables
!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!
!Initial conditions
!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!
!Roughness
!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!
!Bed variables
!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!
!Algorithm options
!!!!!!!!!!!!!!!!!!!!!

nos = 2000 ! Number of spatial grid points
writfreq = 2000 !24*5*5 ! The output is written every writfreq 'th timestep
jmax= 1200000 !1382400 !24*5*2*45*4*4*4*4*5 ! The number of time steps
t = 0.00 ! Starting time (s)
dT = 1.0 !Time step (s)
variable_timestep=.TRUE. !Do we change the timestepping for high sediment concentrations? -- this is presently inconsistent with bed layers

waterM = 0.0 !Initial water elevation (m) and mean water elevation
TR=0.0 !Tidal range. If it is set to zero then the discharge is constant, otherwise the continuity based method is used.
Width = 1.0 !Width of computational domain (m)

no_discharges=1
Discharges = 0.00738
susconcs = 1.0e-04
friction_type = 'manning'!, 'darcy', 'vanrijn', 'manning', 'ks'
rough_coef=0.03
grain_friction_type='vanrijn' !'colebrook'!'vanrijn', 'one', 'onethird'
man_nveg = 0.3 !Mannings n for vegetated regions
veg_ht = 9.0e+20 !Height at which vegetation is assumed to occur.
lambdacon=0.24 !Dimensionless eddy viscosity constant
rho = 1026.0 ! Density of water (kg/m^3)
tbston=.true. !When true (false) this term switches on (off) the sqrt(1+slopes^2) in the bed shear equation: tau*sqrt(1+slopes^2) = rho g Sf h + d/dy ... If .false., then it is replaced with 1. This effects both the shear calculation, and the 'roughmult' friction factor estimation

layers=1 !The number of bed layers
lincrem = 1000.031 ! The distance between bed layers (m). Set it to a very high number to avoid the multi bed layers having any influence.
mu = 10.60 !Angle of repose - this can be used to influence the critical shear stress if the code is adjusted
failure_slope = 0.6 ! Slope at which mass failure occurs
tauinc = 0.00 ! A DEFUNCT constant (Pa)
erconst = 0.13  ! The constant determining the min critical shear and the critical shear increment
taucrit_slope_reduction=.TRUE. ! Does taucrit reduce on a slope
wset = 0.014 ! Settling velocity of sediment in m/s
voidf = 0.4 ! Void fraction (non sediment fraction) of bed = porosity
lifttodrag = 0.0 ! Lift to drag coefficient ratio
hlim = 0.01 ! If the mean depth of the cross section is < hlim (m), then it is treated as dry - so we don't calculate the shear stress
mor = 1.0 !Morphological factor
rhos = 2600.0 !Density of solid sediment (kg/m^3)
dsand = 0.000062 ! dsand from van Rijn (2004) bedload formula (m)
d50 = 0.000149 !0.000062 ! median grain size for bedload formula (m)
g = 9.8 !Gravity (m/s^2)
kvis = 1.0E-06 !Kinematic viscosity of water. 
alpha=0.000228 !The constant for the erosion formula E= alpha*(tau-taue)/sqrt(taue)
Qbedon=.TRUE. !Is bedload active
bedload_type='mpm' ! 'vanrijn', 'mpm'
talmon=.FALSE. !Do we use a talmon lateral bedload closure?
resus_type = 'vanrijn'! ! 'cohesive', 'vanrijn', 'smithmac'

susdist = .TRUE. !Do we have a laterally variable suspended load? This can ONLY treat the case of steady state. The total load flux (bedload + spsuended load) is assumed to be integrated_load_flux
sus_vert_prof='exp' !'exp', 'Rouse'
edify_model='Parabolic' ! 'Constant', 'Parabolic', 'Parabola_const'
x_len_scale=1.0 ! x length scale. In dynamic_sus_dist dC/dx ~= (C -k*C)/x_len_scale
susQbal= .FALSE. ! DEPRECIATED: Is there a balance between the lateral flux of suspended load and bedload? Only relevant if susdist=.true.
integrated_load_flux=-1.0 !DEPRECIATED: The total flux (suspended load + bedload) through the cross-section, in kg/s =  (kg/m^3)*m^2*m/s. Only used if susdist=.true.
sus2d = .false. !Do we use a fully 2d suspended sediment - this is only applicable to the case with many cross sections strung together - the quasi 2d model.
norm=.false. !Is erosion to be directed normal to the bed?
vertical=.true. ! POORLY TESTED .false. CASE: Is the vertical shear method (SKM) to be used (support for Pizzuto method may not be complete, and in this case it should be .true.

readin = .FALSE. !Do we read the initial conditions from a file?
geo = .false. ! RARELY USE .true. ANYMORE: Do we use the 'geotech' algorithm to prevent steep slopes? This is also a bit outdated
smax= 200.0 ! DEPRECIATED: The max slope when geo=.true. A bit outdated
remesh=.FALSE. !Do we remesh
remesh_freq= 5000 !How many time steps before we consider remesh (only active if remesh=.true.)
normmov=.false. !Do the bed points actually shift with the D-E vector? This is only supported for pure suspended load without bed layers.

high_order_shear=.FALSE. !LITTLE EXPERIENCE WITH .true. CASE: Do we try to use a higher order estimate of derivatives in the shear approximation
high_order_bedload=.FALSE. !LITTLE EXPERIENCE WITH .true CASE: Do we try to use a higher order estimate of derivatives in the downslope bedload approximation
high_order_Cflux=.FALSE. !DEPRECIATED: Do we try to use a higher order estimate of derivatives in the dynamic_sus_dist routine
/

