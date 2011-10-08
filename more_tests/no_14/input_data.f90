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
dT = 300.00 !Time step (s)
variable_timestep=.FALSE. !Do we use a variable time-step? Can cause problems

waterM = 0.0 !Initial water elevation (m) and mean water elevation
TR=0.0 !Tidal range. If it is set to zero then the discharge is constant, oth
Width = 100.0 !Width of computational domain (m)

num_simulations=1 ! Number of different discharge simulations
Discharges = 29.4 ! Discharge values
susconcs = 1.0e-03 !suspended sediment concentration(m^3/m^3) for each discha
friction_type = 'manning'!, 'darcy', 'vanrijn'
rough_coef=0.027 ! Friction coefficient corresponding to the friction_type mo
grain_friction_type='vanrijn' !'colebrook'!'vanrijn'
man_nveg = 0.3 !Mannings n for vegetated regions
veg_ht = 9.0e+20 !Height at which vegetation is assumed to occur.
lambdacon=0.24 !Dimensionless eddy viscosity constant
rho = 1026.0 ! Density of water (kg/m^3)
tbston=.true. !When true (false) this term switches on (off) the sqrt(1+slope

layers=1 !The number of bed layers -- layers>1 is not presently compatible wi
lincrem = 1000.031 ! The distance between bed layers (m). Set it to a very hi
mu = .60 !Angle of repose - this can be used to influence the critical shear 
failure_slope = 1.0 ! Slope at which mass failure occurs
erconst = 0.13  ! The constant determining the min critical shear and the cri
taucrit_slope_reduction=.FALSE. ! Does taucrit reduce on a lateral slope?
wset = 0.014 ! Settling velocity of sediment in m/s
voidf = 0.4 ! Void fraction (non sediment fraction) of bed = porosity
lifttodrag = 0.0 ! Lift to drag coefficient ratio
hlim = 0.01 ! If the mean depth of the cross section is < hlim (m), then it i
mor = 1.0 !Morphological factor
rhos = 2600.0 !Density of solid sediment (kg/m^3)
dsand = 0.000062 ! dsand from van Rijn (2004) bedload formula (m)
d50 = 0.000149 !0.000062 ! median grain size for bedload formula (m)
g = 9.8 !Gravity (m/s^2)
kvis = 1.0E-06 !Kinematic viscosity of water. 
alpha=0.000228 !The constant for the erosion formula E= alpha*(tau-taue)/sqrt
Qbedon=.TRUE. !Is bedload active
bedload_type='vanrijn'
talmon=.TRUE. !Do we use a talmon lateral bedload closure?
resus_type = 'vanrijn'! ! 'cohesive', 'vanrijn', 'smithmac'

susdist = .TRUE. !Do we have a laterally variable suspended load? This can ON
sus_vert_prof='Rouse' !'exp', 'Rouse'
edify_model='Parabolic' ! 'Constant', 'Parabolic', 'Parabola_const'
x_len_scale=10000.0 ! x length scale. In dynamic_sus_dist dC/dx ~=(C-k*C)/x_l
susQbal= .FALSE. ! DEPRECIATED: Is there a balance between the lateral flux o
integrated_load_flux=-1.0 !DEPRECIATED: The total flux (suspended load + bedl
sus2d = .false. !Do we use a fully 2d suspended sediment - this is only appli
norm=.false. !Is erosion to be directed normal to the bed?
vertical=.true. ! POORLY TESTED .false. CASE: Is the vertical shear method (S
evolve_bed=.TRUE. ! Do we evolve the bed? If FALSE, hydrodynamics, erosion an

readin = .FALSE. !Do we read the initial conditions from a file?
geo = .false. ! RARELY USE .true. ANYMORE: Do we use the 'geotech' algorithm 
smax= 200.0 ! DEPRECIATED: The max slope when geo=.true. A bit outdated
remesh=.FALSE. !Do we remesh
remesh_freq= 5000 !How many time steps before we consider remesh (only active
normmov=.false. !Do the bed points actually shift with the D-E vector? This i

high_order_shear=.FALSE. !LITTLE EXPERIENCE WITH .true. CASE: Do we try to us
high_order_bedload=.FALSE. !LITTLE EXPERIENCE WITH .true CASE: Do we try to u
high_order_Cflux=.FALSE. !DEPRECIATED: Do we try to use a higher order estima
/

