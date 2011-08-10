Module Pizzutotry
! Module to compute the shear using Pizzuto's approach


!Module of important constants
use global_defs

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



Subroutine bigB(Dn, ks, nn, B)
      !evaluate the term 'B' in the equations by analytical means. 

Integer, intent(in):: nn      !number of elements along the cross section
      
REAL(dp), intent(in):: Dn, ks  !'depth' along normals and roughness height. Dn(2)=Dn(1.5)

REAL(dp),intent(out):: B      !the big B where the output will go B(2) = B(1.5), 

DIMENSION:: Dn(nn), ks(nn), B(nn)


B= (0.5/(Dn-ks))*( Dn**2*log(Dn/ks)*(Dn/6._dp-0.5*ks)-(5._dp/36._dp)*(Dn**3-ks**3)+(0.75_dp)*Dn*(Dn*ks-ks**2)) 

!B=0._dp
!This is an analytical evaluation of that integral.
!DO i=1, nn
!IF(B(i)<0.) THEN 
!        B(i)=0.
!print*, "B(",i,")<0."
!END IF

!END DO
!print*, maxval(B), minval(B)
!stop
      
end subroutine bigB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dAdPP(nn, ys,hs, elev, Dn, dAP,dPer, slps, Ps, ks, useful)

!Calculate that annoying dA/dP term that appears in the equation, and a few
!other geometric things

Integer, intent(in)::nn
REAL(dp), intent(in):: ys, hs,elev, slps, ks !So Ps is a vector of y coordinates for the bed points, hs is their heights, depth is the water depth
REAL(dp),intent(out)::Dn, dAP, dPer, Ps, useful ! the distance along the normals, and the final derivative term (area centered around each point), and the wetted perimeter increment (b!etween each point

DIMENSION:: Dn(nn),dAP(nn),ys(nn),hs(nn), dPer(nn), slps(nn), Ps(nn), useful(nn), ks(nn)

REAL(dp):: nbp(nn,2) !normal "base" points
REAL(dp):: ndir(nn,2) !normal direction vectors
REAL(dp):: tptx(nn) ! water surface point of each normal,x value only
REAL(dp):: btx(nn) ! bottom of trapezium (shifted point on normal), x values only
REAL(dp)::DA(nn), W0(nn), Wt(nn), useful1(nn) !Area, and useful temp variables
REAL(dp)::rnd
Logical:: bbb(nn)
!nbp(2:nn,1)= 0.5*(Ps(1:nn-1)+Ps(2:nn)) !!So this array holds the base points (x,y) of the normal vectors
!nbp(2:nn,2)=0.5*(hs(1:nn-1)+hs(2:nn))
bbb=.false.

nbp(1:nn,1)= ys(1:nn)!0.5*(Ps(1:nn-2)+Ps(2:nn)) !!So this array holds the base points (x,y) of the normal vectors
nbp(1:nn,2)=hs(1:nn) !0.5*(hs(1:nn-1)+hs(2:nn))


ndir(1:nn,1)= -slps(1:nn)!-(hs(2:nn)-hs(1:nn-1)) !And this array holds the direction vectors (r1,r2) of the normals. Note they are prependicular to the difference of y1,h1 and y2,h2
ndir(1:nn,2)= 1._dp !Ps(2:nn)-Ps(1:nn-1)

!Now we compute the trapezium areas, as the first part of getting the areas. 
!Note that the depth at each normal is just the (elev - the y coord of the base point)

tptx(1:nn)= nbp(1:nn,1)+ ((elev-nbp(1:nn,2))/ndir(1:nn,2))*ndir(1:nn,1) !So this is the x value of the point where the normal intersects that water surface (obviously the y value there is elev)

!call random_number(rand)
!if(rand>.998) write(8,*) tptx


!Check that the normals don't intersect -- this is equivalent to checking that
!tptx is increasing.
IF (minval(tptx(2:nn)-tptx(1:nn-1))<0.) THEN
        print*, "normals intersect"
     !  bbb=.false.
      !   bbb(2:nn) = tptx(2:nn)-tptx(1:nn-1)<0._dp 
       ! bbb(1)=.false.
     !   do i= 1, 10
      !  print*, tptx(i), ys(i), slps(i)
       ! end do
        !print*, tptx
      ! STOP
       !TRY NORMAL DEPTH METHOD -- This is often a problem
       dAP(1:nn)= (elev-hs)/cos(atan(slps))! .5*( ( (tptx(2:nn-1)-nbp(2:nn-1,1))**2.+ (elev-nbp(2:nn-1,2))**2.)**.5 +& 
      
       ! ( (tptx(3:nn)-nbp(3:nn,1))**2.+ (elev-nbp(3:nn,2))**2.)**.5 ) !Normal depth, average of the depth at each normal  !(elev-hs(2:nn-1))*(1.+ ((hs(3:nn)-hs(1:nn-2))/(Ps(3:nn)-Ps(1:nn-2)))**2.)**.5!cos(atan(abs( (hs(3:nn)-hs(2:nn-1))/(Ps(3:nn)-Ps(2:nn-1)) )))
      ! dAP(1)=elev-hs(1)
       !dAP(nn)=elev-hs(nn)
       B=0.
       dPer= 1. ! This won't matter since B =0.
       Dn= -1. !(This will ensure that the B term is kept to 0.
      Ps=-1. !This also shouldn't be used
      useful= 0._dp !This will signal to the other code that we are in a special case

      IF ((minval(dAP)<0.).or.(minval(dAP)+1.== minval(dAP))) THEN
        print*, "minval dAP<0.; minval=", minval(dAP), maxval(dAP)
        print*, (1.+((hs(3:nn)-hs(1:nn-2))/(Ps(3:nn)-Ps(1:nn-2)))**2.)**.5
        stop
      end if


      goto 1111  !End of the subroutine

END IF


btx(1:nn-1)= nbp(1:nn-1,1)+((nbp(2:nn,2)-nbp(1:nn-1,2))/ndir(1:nn-1,2))*ndir(1:nn-1,1)
!btx is an important for calculating the area between normals. For any given
!point in the numerical scheme, the area between the normals either side of it
!which contains water is a value that we need to find. 
!We do this by first finding the area of a trapezium with both bases parallel to
!the water surface, and then subtracting or adding as approp the triangle which makes the area equal to the value we want.
!The parallelogram calculation requires taking one of the normal points (the left most in this case), and projecting
!it along its normal until it is at the same height as the other normal point.
!btx stores this x value

!Now we find the area of the trapezium and later subtract the area of the
!triangle

DA(1:nn-1)= 0.5*( (tptx(2:nn)-tptx(1:nn-1))+ (nbp(2:nn,1)-btx(1:nn-1)))*(elev-nbp(2:nn,2)) !Area of the parallelogram with both bases parallel to the water surface

!Now add the triangle to it

DA(1:nn-1)= DA(1:nn-1) +0.5*( nbp(2:nn,1)- btx(1:nn-2))*(nbp(2:nn,2)-nbp(1:nn-1,2)) !This expresssion is "proper area= area of full parallellogram +  .5*base*vertical height

!Note that these area increments are really the increments in the direction of the index increasing. We centre them in dAP 

!First the perimeter increments around the data points.
dPer(2:nn-1)= .5*( ((nbp(3:nn,1)-nbp(2:nn-1,1))**2+(nbp(3:nn,2)-nbp(2:nn-1,2))**2)**0.5 + & 
((nbp(2:nn-1,1)-nbp(1:nn-2,1))**2+(nbp(2:nn-1,2)-nbp(1:nn-2,2))**2)**0.5)  ! dper is actually centred at the channel point (rather than the normal point).

dPer(1)=dPer(2)
dPer(nn)=dPer(nn-1)


dAP(2:nn-1)= .5*(DA(1:nn-2)+ DA(2:nn-1))/dPer(2:nn-1) !So this is the area increment divided by the wetted perimeter increment
dAP(1)= .5*DA(1)/dPer(1)
dAP(nn)= .5*DA(nn-1)/dPer(nn)


!The following variable will serve the role of 'lengths' in the depth averaged
!model
Ps(1)=0._dp
do i=2, nn
Ps(i)= Ps(i-1)+ sqrt( (hs(i)-hs(i-1))**2+(ys(i)-ys(i-1))**2)
end do

!dAP(1)= elev-hs(1) !A crap guess!
!dAP(nn)= elev-hs(nn) !But at least the same degree of crapness
!now compute the lengths of the normals eh

Dn(1:nn)= (elev-hs)/cos(atan(slps)) !!((nbp(1:nn,1)-tptx(1:nn))**2 +(nbp(1:nn,2)-elev)**2)**0.5  !The length of the normal

W0(2:nn)= Ps(2:nn)-Ps(1:nn-1) !Distance along the bed between i and i-1 normal
Wt(2:nn)= max(0._dp, (tptx(2:nn)-tptx(1:nn-1)))*cos(atan((hs(2:nn)-hs(1:nn-1))/(ys(2:nn)-ys(1:nn-1)) )) !Distance at the top of the ith normal between the i and i-1 normal. The max(0,..) is there for the case when the normals intersect, and we use the normal depth method- Crude I know, but our assumptions fail massively here anyway. Hence why I didn't want to use this method. 
!print*, maxval(W0(2:nn)), minval(W0(2:nn)), maxval(Wt(2:nn)), minval(Wt(2:nn))
!Calculate the associated area/velocity shape function for i between i and i-1 -
!this is needed to correct the discharge later - see the appendix of cros4.tex
!for an explanation 
useful1(2:nn)= -.5*((Ks(2:nn)**2-4*Dn(2:nn)*Ks-2*Dn(2:nn)**2*log(Dn(2:nn)/Ks(2:nn))+3*Dn(2:nn)**2)*W0(2:nn)-Wt(2:nn)*Ks(2:nn)**2 & 
-2*Dn(2:nn)**2*Wt(2:nn)*log(Dn(2:nn)/Ks(2:nn))+Dn(2:nn)**2*Wt(2:nn))/(4._dp*Dn(2:nn))
useful(1)= useful1(2)
useful(nn)=useful1(nn)
useful(2:nn-1)= useful1(3:nn)+useful1(2:nn-1)

!print*, maxval(useful), minval(useful), sum(useful)/nn
!So Dn(2) = Dn(1.5) really
!Dn(1)=Dn(2)
!dAP(1)= Dn(1) !A guess. 
!dAP(nn)=Dn(nn) 

!do i=1, nn
!print*, elev-hs(i), Dn(i), slps(i)
!end do

!do i= 2, nn-1
!if ((bbb(i)).or.(bbb(i+1)).or.(bbb(i-1))) THEN
!        Dn(i)=ks !1This will ensure B there is zero
!        dAP(i)= (elev-hs(i))/cos(atan(slps(i)))
!end if
!if((bbb(1)).or.(bbb(2))) THEN
!        Dn(1)= ks
!        dAP(1)= (elev-hs(1))/cos(atan(slps(1)))
!END IF
!if((bbb(nn)).or.(bbb(nn-1))) THEN
!        Dn(nn)= ks
!        dAP(nn)= (elev-hs(nn))/cos(atan(slps(nn)))
!END IF
!
!end do

1111 End Subroutine dAdPP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine shearP(nn,ys,hs,elev,wslope,tau,ks, f, NNN,slps,counter, Q, vegdrag, rho, useful)
!!At the moment some of the arguments are dummies that are useful to makethe
!function compatible with other code I have 

Integer, intent(in)::nn, counter
REAL(dp), intent(in):: ys, hs,elev,wslope, ks,f,NNN, slps, Q, vegdrag, rho !So Ps is a vector of y coordinates for the bed points, hs is their heights, depth is the water depth
REAL(dp), intent(out):: tau !shear
REAL(dp), intent(in out):: useful

Dimension:: ys(nn),hs(nn),tau(nn), ks(nn), f(nn), NNN(nn), slps(nn), vegdrag(nn), useful(nn) 

Integer:: info
REAL(dp)::g=9.8_dp
REAL(dp):: dAP(nn), B(nn),Dn(nn), dPer(nn) !derivative term, B integral term, distance along normals term, wetted perimiter increment term
REAL(dp)::alpht(nn),alphb(nn), s(nn)  !first 2 are used in the matrix to represent the diagonals (top and bottom), and last is the right hand side
REAL(dp):: diag(nn), Ps(nn), p0, pmax, Bderiv(nn), l(nn), u(nn), ds(nn), dyf(nn)  !the main diagonal of the matrix
!Calculate distance along normals, and the derivative term

Dn=ks

call dAdPP(nn, ys,hs, elev, Dn, dAP,dPer, slps, Ps, ks, useful)

!Calculate the B term
IF((nn>3).and.(maxval(Dn)>0.)) THEN
        call bigB(Dn, ks, nn, B)
ELSE
        B=0._dp !Ignore lateral momentum exchange
END IF



IF((maxval(B)>0._dp).and.(maxval(useful)>0._dp)) THEN !Solve with lateral momentum exchange
       ! B=0.

p0= 2._dp*Ps(1)- Ps(2)
pmax=  2._dp*Ps(nn)-Ps(nn-1) !Assume that the dry point is here. 


dyf(1:nn-1)= Ps(2:nn)-Ps(1:nn-1) !Forward difference increment along the wetted perimeter

!!Here we use the compact difference strategy
alpht(2:nn-1) = - 1._dp/(.5_dp*(dyf(2:nn-1)+dyf(1:nn-2)))*( .5_dp*(B(2:nn-1)+B(3:nn))*1._dp/dyf(2:nn-1)) !Upper diagonal
alphb(2:nn-1) = - 1._dp/(.5_dp*(dyf(2:nn-1)+dyf(1:nn-2)))*(.5_dp*(B(1:nn-2)+B(2:nn-1))*1._dp/dyf(1:nn-2)) !Lower diagonal
diag(2:nn-1)=  1._dp - alphb(2:nn-1) -alpht(2:nn-1) !Main diagonal

!!The friction slope term, on the right hand side

s(1:nn)= rho*g*dAP*max(wslope,1._dp) !!Note that with the way we do it now, the water slope is irrelevent -- the shear is corrected below to ensure that the discharge is equal to the imposed discharge. We could have any constant, but for historical reasons I keep using wslope (!). 

!print*, dAP(1), dAP(2),dAP(3)
!stop

!!Boundary conditions
alpht(1)=  0._dp ! - 1._dp/(.5_dp*(dyf(1)+dyf(1)))*( .5_dp*(B(2)+B(1))*1._dp/dyf(1)) 
alphb(1)= 0._dp
diag(1)= 1._dp !-alpht(1) + .5_dp*B(1)*1._dp/dyf(1)

alpht(nn)=0._dp
alphb(nn)= 0._dp ! - 1._dp/(.5_dp*(dyf(nn-1)+dyf(nn-1)))*( .5_dp*(B(nn)+B(nn-1))*1._dp/dyf(nn-1))
diag(nn) =  1._dp ! -alphb(nn) + .5_dp*B(nn)*1._dp/dyf(nn-1)




tau=s !Predefine this for matrix solver
l=alphb!/diag 
ds=diag!/diag
u=alpht!/diag

!!Matrix solver
call DGTSV(nn,1,l(2:nn),ds,u(1:nn-1), tau,1, info)

       
ELSE !This is a fallback calculation for when we do not have enough cross sectional points to solve the shear distribution
        useful=0._dp
        tau=rho*g*dAP*wslope!max(wslope,1._dp)
       
END IF

DO i=1, nn
If ((tau(i)+1.==tau(i)).or.(tau(i)<0.)) THEN

        print*, "tau is nan or <0 in pizzutotry", tau(i)
        print*, dAP(i), maxval(B), elev-hs
        print*, abs(diag)-abs(alphb)-abs(alpht)

        print*, " ....... "
        print*, elev -hs
        print*, "........."
        print*, (1.+ ((hs(3:nn)-hs(1:nn-2))/(Ps(3:nn)-Ps(1:nn-2)))**2.)**.5
        print*, "...."
        print*, B

        stop
END IF
END DO



END Subroutine shearP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Pizzutotry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!

!!!!!!!!!!!!
!Program crosser
!use pizzutotry

!implicit none
!Integer:: nos=200, i,j
!REAL(dp)::wslope= 0.5**2/(80**2*2), ar
!REAL(dp):: lengths,heights, water, dists, tau,ks
!Allocatable lengths(:), heights(:), dists(:), tau(:), ks(:)

!open(1,file="data") 

!Allocate (lengths(nos),heights(nos),dists(nos),tau(nos),ks(nos))

!water=1.1  !Water elevation
!
!Assign cross section
!
!DO i =1,nos
!lengths(i)= (40./nos)*i
!
!heights(i)= 0. !+0.02*(lengths(i)-lengths(1))*(lengths(i)-nos*lengths(1))
!
!
!ks(i)=0.003
!END DO

!write(1,*) tau
!write(1,*), water-heights

!DO j=1,1000000

!DO i=1,nos
!IF(heights(i)<(water-0.001)) heights(i)=heights(i)+0.000002
!END DO

!ar= sum(water-heights)*lengths(1) !area of section
!wslope= (20/ar)**2/(80**2*(ar/(lengths(1)*nos))) !water surface slope using Chezy

!call shearP(nos, lengths,heights,water, wslope ,tau,ks, 0.*ks, 0.*ks)

!DO i = 1, nos
  
!IF ((tau(i)>0.4).AND.((water-heights(i))< 500)) heights(i)= heights(i)- 0.0005*(tau(i)-0.4)
!END DO

!IF(mod(j,10000).eq.0) write(1,*), water-heights

!END DO

!write(1,*) tau
!write(1,*), water-heights
!close(1)

!End program      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

