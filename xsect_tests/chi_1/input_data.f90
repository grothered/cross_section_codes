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

nos = 800 ! Number of spatial grid points
writfreq = 8000 !24*5*5 ! The output is written every writfreq 'th timestep
jmax= 2000001 !2000001 !1382400 !24*5*2*45*4*4*4*4*5 ! The number of time steps
t = 0.00 ! Starting time (s)
dT = 300.00 !Time step (s)
variable_timestep=.FALSE. !Do we change the timestepping for high sediment concentrations? -- this is presently inconsistent with bed layers

waterM = 0.00 !Initial water elevation (m) and mean water elevation
TR=0.0 !Tidal range. If it is set to zero then the discharge is constant, otherwise the continuity based method is used.
Width = 800.0 !Width of computational domain (m)

num_simulations= 1
Discharges = 2100.0  !2100.0  2100.0  2100.0  2100.0  2100.0  2100.0 
susconcs =   1.0e-03 !1.0e-04 8.0e-04 4.0e-04 1.0e-05 6.0e-05  
friction_type = 'manning' !'manning' !, 'darcy', 'vanrijn'
rough_coef = 0.039 !0.039 !Mannings n for the basic bed
grain_friction_type  = 'vanrijn' ! 'vanrijn', 'colebrook', 'onethird' -- model for bed shear acting on grains 
man_nveg = 0.3 !Mannings n for vegetated regions
veg_ht = 1000.1 !Height at which vegetation is assumed to occur.
lambdacon=0.24 !Dimensionless eddy viscosity constant
rho = 1026.0 ! Density of water (kg/m^3)
tbston=.TRUE. !When true (false) this term switches on (off) the sqrt(1+slopes^2) in the bed shear equation: tau*sqrt(1+slopes^2) = rho g Sf h + d/dy ... If .false., then it is replaced with 1. This effects both the shear calculation, and the 'roughmult' friction factor estimation

layers=1 !The number of bed layers
lincrem = 1000.031 ! The distance between bed layers (m). Set it to a very high number to avoid the multi bed layers having any influence.
mu = 10000.60 !Angle of repose - this can be used to influence the critical shear stress if the code is adjusted

bank_erosion_type='Delft' !'Delft_if_too_steep' !'Delft' !'Delft_if_big_jump' !'Delft_if_too_steep'
failure_slope = 0.50 ! Slope at which mass failure occurs.
erconst = 0.14  ! The constant determining the min critical shear and the critical shear increment
taucrit_slope_reduction=.FALSE. ! Is taucrit reduced on lateral slopes?
wset = 0.016 ! Settling velocity of sediment in m/s
voidf = 0.4 ! Void fraction (non sediment fraction) of bed = porosity
lifttodrag = 0.0 ! Lift to drag coefficient ratio
hlim = 0.01 ! If the mean depth of the cross section is < hlim (m), then it is treated as dry - so we don't calculate the shear stress
mor = 1.0 !Morphological factor
rhos = 2600.0 !Density of solid sediment (kg/m^3)
dsand = 0.000062 ! dsand from van Rijn (2004) bedload formula (m)
d50 = 0.00015 !0.000062 ! median grain size for bedload formula (m)
g = 9.8 !Gravity (m/s^2)
kvis = 1.0E-06 !Kinematic viscosity of water. 
alpha=0.000228 !The constant for the erosion formula E= alpha*(tau-taue)/sqrt(taue)
Qbedon=.FALSE. !Is bedload active
bedload_type='vanrijn' ! 'mpm', 'vanrijn'
talmon=.TRUE. !Do we use a talmon lateral bedload closure?
resus_type = 'vanrijn' ! 'cohesive', 'vanrijn'

susdist = .TRUE. !Do we have a laterally variable suspended load? 
sus_vert_prof='Rouse' !'exp', 'Rouse'
edify_model='Parabolic' ! 'Constant', 'Parabolic', 'Parabola_const'
x_len_scale=1000.0 ! x length scale. In dynamic_sus_dist dC/dx ~= (C -k*C)/x_len_scale
sus2d = .false. !Do we use a fully 2d suspended sediment - this is only applicable to the case with many cross sections strung together - the quasi 2d model.
norm=.FALSE. !Is erosion to be directed normal to the bed?
vertical=.TRUE. !Is the vertical shear method (SKM) to be used (support for Pizzuto method may not be complete, and in this case it should be .true.
evolve_bed=.TRUE. ! Do we evolve the bed? If FALSE, hydrodynamics, erosion and deposition are computed, but no change to the bed occurs.

readin = .FALSE. !Do we read the initial conditions from a file?
geo = .false. !Do we use the 'geotech' algorithm to prevent steep slopes? This is also a bit outdated
smax= 200.0 !The max slope when geo=.true. A bit outdated
remesh=.TRUE. !Do we remesh
remesh_freq= 3000 !How many time steps before we consider remesh (only active if remesh=.true.)
normmov=.false. !Do the bed points actually shift with the D-E vector? This is only supported for pure suspended load without bed layers.

high_order_shear=.FALSE. ! Do we try to use a higher order estimate of derivatives in the shear approximation
high_order_bedload=.FALSE. !Do we try to use a higher order estimate of derivatives in the downslope bedload approximation
/

