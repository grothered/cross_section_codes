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
writfreq = 12000 !24*5*5 ! The output is written every writfreq 'th timestep
jmax= 2400000 !1382400 !24*5*2*45*4*4*4*4*5 ! The number of time steps
t = 0.00 ! Starting time (s)
dT = 100.00 !Time step (s)
variable_timestep=.FALSE. !Do we change the timestepping for high sediment concentrations? -- this is presently inconsistent with bed layers

waterM = 0.0 !Initial water elevation (m) and mean water elevation
TR=0.0 !Tidal range. If it is set to zero then the discharge is constant, otherwise the continuity based method is used.
Width = 200.0 !Width of computational domain (m)

no_discharges=3
Discharges = 30.0 24.0 18.0
susconcs =   0.1  0.1  0.1
friction_type = 'darcy'
rough_coef = 0.07056 !Mannings n for the basic bed
grain_friction_type='one' ! 'colebrook, vanrijn, onethird
man_nveg = 0.3 !Mannings n for vegetated regions
veg_ht = 0.3 !Height at which vegetation is assumed to occur.
lambdacon=0.0 !Dimensionless eddy viscosity constant
rho = 1026.0 ! Density of water (kg/m^3)
tbston=.false. !When true (false) this term switches on (off) the sqrt(1+slopes^2) in: tau*sqrt(1+slopes^2) = rho g Sf h + d/dy ... If .false., then it is replaced with 1. This effects both the shear calculation, and the 'roughmult' friction factor estimation

layers=1 !The number of bed layers
lincrem = 100.031 ! The distance between bed layers (m)
mu = 100.0 !Angle of repose - this influences the critical shear stress
erconst = 0.4  ! The constant determining the min critical shear and the critical shear increment
wset = 0.0004 ! Settling velocity of sediment in m/s
voidf = 0.4 ! Void fraction (non sediment fraction) of bed = porosity
lifttodrag = 0.0 ! Lift to drag coefficient ratio
hlim = 0.01 ! If the mean depth of the cross section is < hlim (m), then it is treated as dry - so we don't calculate the shear stress
mor = 45.0 !Morphological factor
rhos = 2600.0 !Density of solid sediment (kg/m^3)
dsand = 0.000062 ! dsand from van Rijn (2004) bedload formula (m)
d50 = 0.000062 !0.000062 ! median grain size for bedload formula (m)
g = 9.8 !Gravity (m/s^2)
kvis = 1.0E-06 !Kinematic viscosity of water. 
alpha=0.000228 !The constant for the erosion formula E= alpha*(tau-taue)/sqrt(taue)
Qbedon=.true. !Is bedload active
bedload_type='vanrijn'
talmon=.false.
resus_type='cohesive'

susdist = .false. !Do we have a laterally variable suspended load? This can ONLY treat the case of steady state, with suspended load only. The near bed concentration at the centre is assumed to be = sconc
susQbal= .false. !Is there a balance between the lateral flux of suspended load and bedload? Only relevant if susdist=.true.
integrated_load_flux= -1.0 !The total flux (suspended load + bedload) through the cross-section, in kg/s =  (kg/m^3)*m^2*m/s. Only used if susdist=.true., but NOT in dynamic_sus_dist, which is presently in favour. Hence I set it negative.
sus2d = .false. !Do we use a fully 2d suspended sediment - this is only applicable to the case with many cross sections strung together - the quasi 2d model.
norm=.false. !Is erosion to be directed normal to the bed?
vertical=.true. !Is the vertical shear method (SKM) to be used (support for Pizzuto method may not be complete, and in this case it should be .true.
evolve_bed=.TRUE. ! Do we evolve the bed? If FALSE, hydrodynamics, erosion and deposition are computed, but no change to the bed occurs.

readin = .false. !Do we read the initial conditions from a file?
geo = .false. !Do we use the 'geotech' algorithm to prevent steep slopes? This is also a bit outdated
smax= 200.0 !The max slope when geo=.true. A bit outdated
remesh=.TRUE. !Do we remesh
remesh_freq= 5000 !How many time steps before we consider remesh (only active if remesh=.true.)
normmov=.false. !Do the bed points actually shift with the D-E vector? This is only supported for pure suspended load without bed layers.

high_order_shear=.FALSE. ! Do we try to use a higher order estimate of derivatives in the shear approximation
high_order_bedload=.FALSE. !Do we try to use a higher order estimate of derivatives in the downslope bedload approximation
high_order_Cflux=.FALSE. ! Do we try to use a higher order estimate of derivatives in the dynamic_sus_dist routine
/

