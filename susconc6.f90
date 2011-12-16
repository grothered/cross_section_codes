module sus

use global_defs
!use crosssection

IMPLICIT NONE
!INTEGER, PARAMETER, private  :: dp = SELECTED_REAL_KIND(12, 60)

CONTAINS
SUBROUTINE susconc(n,DT, A,delX,C, U, Qe, Qd, Cmouth, C_old, Area_old, UU_old, Criver, wset, wetwidth, D)
    INTEGER, intent(in):: n
    REAL(dp), intent(in):: DT, delX, Qe, Qd, Cmouth, A, U, C_old, Area_old, UU_old,Criver, wset, wetwidth, D
    REAL(dp), intent(in out):: C


    !!Note - here Qe is different to elsewhere. 
    Dimension C(n), A(n), U(n), Qe(n), Qd(n), C_old(n), Area_old(n), UU_old(n), wetwidth(n), D(n)

    INTEGER:: info, i
    REAL(dp):: diag(n), up(n), lo(n)  !Matrix diagonals
    REAL(dp):: r(n) !Right hand side of equation
    logical:: flag
    !REAL(dp):: D(n) !dispersion constant-- actually not constant, need to fix.
    REAL(dp):: impcon=.5_dp
    REAL(dp):: Qf1(n), Qf0(n), Af1(n), Af0(n), Df1(n), Df0(n)

    !D= 2._dp!*abs(U)
    !C(n)=0.
    !!!
    !!##Solves (del AC/ del T) +  (del CQ / delX) =   d/dx (A D del C / delX) + (E-D)  

    !!Note that E and D are the total rate of erosion/deposition over each section. 
    !!Solved with a Crank-Nicholson Method. Although this is a common reference, I
    !found Vreugdenhill (1989:59) a useful reference, and the way the code is
    !written reflects that.  
    !!At the mouth boundary, a given value is imposed
    !! At the landward boundary, a zero gradient condition is enforced if there is
    !no river inflow. Otherwise, the river concentration is enforced, or used to
    !provide other boundary values, depending on the version of this code you have.  

    !!!START
    !Calculate the diagonals for the matrix inversion

    !!Forward average of some important variables
    Qf1(1:n-1)=.5_dp*(A(2:n)*U(2:n)+A(1:n-1)*U(1:n-1))
    Af1(1:n-1)= .5_dp* (A(2:n)+A(1:n-1))
    Df1(1:n-1)= .5_dp*(D(1:n-1)+D(2:n))

    !!Forward average of some important old variables
    Qf0(1:n-1)=.5_dp*(Area_old(2:n)*UU_old(2:n)+Area_old(1:n-1)*UU_old(1:n-1))
    Af0(1:n-1)= .5_dp* (Area_old(2:n)+Area_old(1:n-1))
    Df0(1:n-1)= .5_dp*(D(1:n-1)+D(2:n))


    !!Upper diagonal
    up(1:n-1)=(impcon/delX)*( Qf1(1:n-1)*.5_dp-Af1(1:n-1)*Df1(1:n-1)/delX)
    up(n)=0._dp
    !!Lower diagonal
    lo(2:n)=(impcon/delX)*(-Qf1(1:n-1)*.5_dp-Af1(1:n-1)*Df1(1:n-1)/delX)!-up(1:n-1) !-(impcon/delX)*& ((.5_dp*(A(2:n)*U(2:n)+A(1:n-1)*U(1:n-1)))*.5_dp+.5_dp*(A(2:n)+A(1:n-1))*.5_dp*(D(1:n-1)+D(2:n))/delX)
    lo(1)=0._dp

    !Main diagonal
    diag(2:n-1)= A(2:n-1)/DT + wset*wetwidth(2:n-1)+ &
    (impcon/delX)*(Qf1(2:n-1)*.5_dp - Qf1(1:n-2)*.5_dp  + Af1(2:n-1)*Df1(2:n-1)/delX + Af1(1:n-2)*Df1(1:n-2)/delX) !up(2:n-1) +(impcon/delX)*.5_dp*(A(2:nn-1)+A(3:nn))*(D(2:n-1)+D(3:nn))/delX + lo(2:n-1) &
    !+ .5_dp*(A(2:nn-1)+A(1:nn-2))*(D(1:n-2)+D(2:n-1))/delX

    !!Right hand side
    r(2:n-1) = Qe(2:n-1) +Area_old(2:n-1)*C_old(2:n-1)/DT - (1._dp-impcon)/delX*( & 
    ( Qf0(2:n-1)*.5_dp*(C_old(3:n)+C_old(2:n-1)) & !Advective
    - Af0(2:n-1)*Df0(2:n-1)*(C_old(3:n)-C_old(2:n-1))/delX ) & !diffusive
    - ( Qf0(1:n-2)*.5_dp*(C_old(2:n-1)+C_old(1:n-2)) & !Advective
    - Af0(1:n-2)*Df0(1:n-2)*(C_old(2:n-1)-C_old(1:n-2))/delX ) ) !diffusive

    !!!Boundary conditions

    !!CASE 1 - At the mouth, ebbing
    ! Assume the discharge and area at the mouth are the same as at 1, and the sediment
    !concentration beyond the boundary is the same at at 1. 
    If(U(1)<=0._dp) THEN
    diag(1)= A(1)/DT+ wset*wetwidth(1)+ (impcon/delX)*(Qf1(1)*.5_dp - Qf1(1)*.5_dp*2._dp  + Af1(1)*Df1(1)/delX + 0._dp ) !So the  Qf(1)*.5_dp*2._dp reflects the addition of the value at the mouth, and the final 0._dp reflects the Af(0)*Df(1)/delX*(C(1)-C(mouth))/delX which is zero in this situation

    r(1) = Qe(1) + Area_old(1)*C_old(1)/DT - (1._dp-impcon)/delX*( &
    (Qf0(1)*.5_dp*(C_old(1)+C_old(2)) - Af0(1)*Df0(1)*(C_old(2)-C_old(1))/delX) &
    -( Qf0(1)*.5_dp*(C_old(1)+C_old(1)) - Af0(1)*Df0(1)*(0._dp) )) !!Note that here, C_old(1) takes the place of C_mouth
    END IF

    !!! CASE 2- At the mouth, flooding
    If(U(1)>=0._dp) THEN
    diag(1)= A(1)/DT+ wset*wetwidth(1)+ (impcon/delX)*(Qf1(1)*.5_dp - Qf1(1)*.5_dp + Af1(1)*Df1(1)/delX + Af1(1)*Df1(1)/delX ) !Here the c_mouth terms go on the RHS

    !Note that the beginning of the RHS includes the Cmouth terms, which are at time
    !1 rather than 0
    r(1) = impcon/delX*(Qf1(1)*.5_dp*Cmouth +Af1(1)*Df1(1)/delX*Cmouth) + Qe(1) + Area_old(1)*C_old(1)/DT - & 
    (1._dp-impcon)/delX*( (Qf0(1)*.5_dp*(C_old(1)+C_old(2)) - Af0(1)*Df0(1)*(C_old(2)-C_old(1))/delX) &
    -( Qf0(1)*.5_dp*(C_old(1)+ Cmouth) - Af0(1)*Df0(1)*(C_old(1)-Cmouth)/delX))
    END IF



    !!!!Next case - the behaviour at the landward end - I suggest a zero gradient
    !condition is reasonable if there is no river input

    if(Criver==0._dp) THEN
    diag(n)= A(n)/DT+ wset*wetwidth(n)+ impcon/delX*(A(n)*U(n)*.5_dp*2._dp - Qf1(n-1)*.5_dp +0._dp + Af1(n-1)*Df1(n-1)/delX) !Here the A(n)*U(n)*.5_dp*2._dp accounts for the zero gradient in the upper advective part, and the  +0._dp also accounts for the zero gradient in the upper diffusive part
    r(n)= Qe(n)+ Area_old(n)*C_old(n)/DT - (1._dp-impcon)/delX*( &
    (Area_old(n)*UU_old(n)*C_old(n) - 0._dp)- (Qf0(n-1)*.5_dp*(C_old(n) +C_old(n-1)) &
    -Af0(n-1)*Df0(n-1)*(C_old(n)-C_old(n-1))/delX))  
    ELSE
    diag(n)= A(n)/DT+wset*wetwidth(n)+ impcon/delX*( A(n)*U(n)*.5_dp - Qf1(n-1)*.5_dp +A(n)*D(n)/delX + Af1(n-1)*Df1(n-1)/delX) 

    r(n)= Qe(n)+ impcon/delX*(-A(n)*U(n)*.5_dp*Criver + A(n)*D(n)/delX*Criver) + Area_old(n)*C_old(n)/DT & 
    - (1._dp-impcon)/delX*((Area_old(n)*UU_old(n)*.5_dp*(C_old(n)+Criver) -Area_old(n)*D(n)*(Criver-C_old(n))/delX )- &
    (Qf0(n-1)*.5_dp*(C_old(n) +C_old(n-1)) -Af0(n-1)*Df0(n-1)*(C_old(n)-C_old(n-1))/delX))  

    !diag(n)=1._dp
    !r(n)=Criver
    !lo(n)=0._dp
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!Solve

    !print*, maxval(C_old), minval(C_old)
    !print*, lo, diag, up, r
    !stop
    call DGTSV(n,1, lo(2:n), diag, up(1:n-1), r,n, info)

    !The solver writes the new C to r- so we change it back.
    C=r

    !print*, info
    !print*, maxval(diag), minval(diag), maxval(lo), minval(lo), maxval(up), minval(up)
    !print*, "vel", maxval(U), minval(U)

    IF(minval(C)<0._dp) THEN 
             print*, "min Sediment conc <0", minval(C),  maxval(C_old), minval(C_old), maxval(U), minval(U)
             print*, "..........minloc C is", minloc(C), "............"
    !         print*, r
    !         print*, "..........."
             print*, maxval(Qe), minval(Qe), maxval(Qd), minval(Qd)
         !print*, ">>>>>>>>>>>>>"
         !print*, C
    !         !print*, A !diag-(abs(up)+abs(lo))
    do i=1, n
    IF(C(i)<0._dp) THEN

        IF(abs(C(i))<10._dp**(-10)) THEN !This could just be due to round off error in the matrix solver - fix it.
        C(i)=0._dp
        ELSE
        print*, "C is < 0; violation is", abs(C(i))
    !	stop
        END IF
    END IF	
    END DO
    !         
    !        stop
    END IF

    !DO i=1, n
    !IF(C(i)<0._dp) C(i)=0._dp
    !END DO


    DO i = 1, n
    Flag=isnan(C(i))
    if(Flag) THEN
           print*, "sedconc is NAN"
           print*, U
           print*, "................"
           print*, A
           print*, Qd
           print*, Qe
           
            

            stop
    END IF
    END DO

end subroutine susconc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine susconc_up3(n2,DT, A2,delX,C2, U2, Qe2, Qd2, Cmouth, C_old2, Area_old2, UU_old2, Criver, wset, wetwidth2,&
 D2,D2_old, third)
    INTEGER, intent(in):: n2
    REAL(dp), intent(in):: DT, delX, Qe2, Qd2, Cmouth, A2, U2, C_old2, Area_old2, UU_old2,Criver, wset, wetwidth2, D2, D2_old
    REAL(dp), intent(in out):: C2
    LOGICAL, intent(in):: third


    !!Note - here Qe is different to elsewhere. 
    Dimension C2(n2), A2(n2), U2(n2), Qe2(n2), Qd2(n2), C_old2(n2), Area_old2(n2), UU_old2(n2), wetwidth2(n2), D2(n2), D2_old(n2)

    INTEGER:: info, i, KL=2, KU=2, IPV(n2+4), ii !Note KL,KU are the number of upper and lower diagonals of the banded matrix
    REAL(dp):: r(n2) !Right hand side of equation
    logical:: flag
    !REAL(dp):: D(n) !dispersion constant-- actually not constant, need to fix.
    REAL(dp):: impcon=.50_dp, impconU=.00_dp
    !REAL(dp):: Qf1(n), Qf0(n), Af1(n), Af0(n), Df1(n), Df0(n)
    REAL(dp):: band(7,n2+4), rhs(n2+4)
    REAL(dp):: A(n2+4), C(n2+4), U(n2+4), Qe(n2+4), Qd(n2+4), C_old(n2+4), Area_old(n2+4), UU_old(n2+4), wetwidth(n2+4), &
    D(n2+4), D_old(n2+4)

    !!Define new variables with 'ghost'points on the edges, useful for implementing
    !boundary conditions
    A(3:(n2+2))=A2
    C(3:(n2+2))=C2
    U(3:(n2+2))=U2
    Qe(3:(n2+2))=Qe2
    Qd(3:(n2+2))=Qd2
    C_old(3:(n2+2))=C_old2
    Area_old(3:(n2+2))=Area_old2
    UU_old(3:(n2+2))=UU_old2
    wetwidth(3:(n2+2))=wetwidth2
    D(3:(n2+2))=D2
    D_old(3:(n2+2))=D2_old

    !Boundary conditions
    A(1:2)=A2(1)
    A((n2+3):(n2+4)) = A2(n2)
    Area_old(1:2)=Area_old2(1)
    Area_old((n2+3):(n2+4)) = Area_old2(n2)
    U(1:2)=U2(1)
    U((n2+3):(n2+4)) = U2(n2)
    UU_old(1:2)=U2(1)
    UU_old((n2+3):(n2+4)) = UU_old2(n2)
    Qe(1:2)=0._dp
    Qe((n2+3):(n2+4)) = 0._dp 
    Qd(1:2)=0._dp
    Qd((n2+3):(n2+4)) = 0._dp 
    wetwidth(1:2)=wetwidth2(1)
    wetwidth((n2+3):(n2+4)) = wetwidth2(n2)
    D(1:2)=D2(1)
    D((n2+3):(n2+4))=D2(n2)
    D_old(1:2)=D2_old(1)
    D_old((n2+3):(n2+4))=D2_old(n2)

    if(U2(1)>0._dp) THEN
    C_old(1:2)=Cmouth
    ELSE
    C_old(1:2)=C_old2(1)
    END IF

    if(U2(n2)>0._dp) THEN
    C_old((n2+3):(n2+4)) = C_old2(n2) !0._dp !Criver 
    ELSE
    C_old((n2+3):(n2+4)) = Criver 
    END IF

    !D= 2._dp!*abs(U)
    !C(n)=0.
    !!!
    !!##Solves (del AC/ del T) +  (del CQ / delX) =   d/dx (A D del C / delX) + (E-D)  

    !!Note that E and D are the total rate of erosion/deposition over each section. 
    !!Solved with a Crank-Nicholson Method, using a third order upwind
    !discretization of the advection term. 

    !Although this is a common reference, I
    !found Vreugdenhill (1989:59) a useful reference, and the way the code is
    !written reflects that. 
    !For the third order advection, I have used 'ghost' cells. There are 2
    !downstream, and 2 upstream
    !!At the mouth boundary, a given value is imposed if the flow is inward, and
    !otherwise a zero gradient condition is imposed.
    !! At the landward boundary, a zero gradient condition is enforced if there is
    !outflow. Otherwise, the river concentration is enforced, or used to
    !provide other boundary values, depending on the version of this code you have.  
    !An elementary reference, which treats the boundaries a bit differently to me,
    !is Karahan (2007). 


    !Predefine
    band=0._dp
    rhs=0._dp
    DO i=3,n2+2
    IF ((U(i)>0._dp)) THEN 
    !Note the unusual banded storage that lapack does
    !!Upper-upper diagonal
            !IF(sign((C_old(i+1)-C_old(i))*(C_old(i)-C_old(i-1)), 1._dp).gt.0._dp) THEN
            !IF((C_old(i).gt.0.0E-7_dp).and.(C_old(i)>0.5_dp*max(C_old(i+1),C_old(i-1)))) THEN
            if(third.eqv..TRUE.) THEN
            band(KL+KU+1+i-(i+2),i+2)= 0._dp !Advection

            !!Upper-diagonal
            band(KL+KU+1+i-(i+1),i+1) = impcon/delX*1._dp/6._dp*2._dp*(U(i+1)*A(i+1)) & !Advection
            -impcon/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX   ) !Diffusion

            !!Main diagonal
            band(KL+KU+1+i-i,i) = A(i)/DT +wset*wetwidth(i) + & !Unsteady
            impcon/delX*1._dp/6._dp*3._dp*(U(i)*A(i)) & !Advection
            +impcon/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX  ) & !Diffusion
            +impcon/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-diagonal
            band(KL+KU+1+i-(i-1),i-1)= -impcon/delX*1._dp/6._dp*6._dp*(U(i-1)*A(i-1)) & !Advection
            -impcon/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-lower-diagonal
            band(KL+KU+1+i-(i-2),i-2)= impcon/delX*1._dp/6._dp*1._dp*(U(i-2)*A(i-2)) !Advection

            !Right hand side
            rhs(i)= Qe(i) + Area_old(i)*C_old(i)/DT -(1._dp-impcon)/delX*( &
            +1._dp/6._dp*2._dp*(UU_old(i+1)*Area_old(i+1))*C_old(i+1) & !Advective
            + 1._dp/6._dp*3._dp*(UU_old(i)*Area_old(i))*C_old(i) & !Advective
            -0.5_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i))*(C_old(i+1)-C_old(i))/delX & !Diffusive
            -1._dp/6._dp*6._dp*(UU_old(i-1)*Area_old(i-1))*C_old(i-1) & !Advective
            -0.5_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1))*(C_old(i)-C_old(i-1))/delX  &!Diffusive
            + 1._dp/6._dp*1._dp*(UU_old(i-2)*Area_old(i-2))*C_old(i-2) ) !Advective
            ELSE
            !Here we use first order upwind

            band(KL+KU+1+i-(i+2),i+2)= 0._dp !Advection

            !!Upper-diagonal
            band(KL+KU+1+i-(i+1),i+1) = 0._dp & !impconU/delX*1._dp/6._dp*2._dp*(U(i+1)*A(i+1)) & !Advection
            -impconU/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX   ) !Diffusion

            !!Main diagonal
            band(KL+KU+1+i-i,i) = A(i)/DT +wset*wetwidth(i) + & !Unsteady
            impconU/delX*(U(i)*A(i)) & !Advection
            +impconU/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX  ) & !Diffusion
            +impconU/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-diagonal
            band(KL+KU+1+i-(i-1),i-1)= -impconU/delX*(U(i-1)*A(i-1)) & !Advection
            -impconU/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-lower-diagonal
            band(KL+KU+1+i-(i-2),i-2)= 0._dp!1._dp/delX*1._dp/6._dp*1._dp*(U(i-2)*A(i-2)) !Advection

            !Right hand side
            rhs(i)= Qe(i) + Area_old(i)*C_old(i)/DT -(1._dp-impconU)/delX*( &
            + (UU_old(i)*Area_old(i))*C_old(i) & !Advective
            -0.5_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i))*(C_old(i+1)-C_old(i))/delX & !Diffusive
            -(UU_old(i-1)*Area_old(i-1))*C_old(i-1) & !Advective
            -0.5_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1))*(C_old(i)-C_old(i-1))/delX ) !Diffusive
            END IF

    ELSE

    !Note the unusual banded storage that lapack does
            !IF((C_old(i).gt.0.0E-7_dp).and.(C_old(i)>0.5_dp*max(C_old(i+1),C_old(i-1)))) THEN
            !IF(sign((C_old(i+1)-C_old(i))*(C_old(i)-C_old(i-1)), 1._dp).gt.0._dp) THEN
            if(third.eqv..true.) THEN 
            !!Upper-upper diagonal
            band(KL+KU+1+i-(i+2),i+2)=  -impcon/delX*1._dp/6._dp*1._dp*(U(i+2)*A(i+2))!Advection

            !!Upper-diagonal
            band(KL+KU+1+i-(i+1),i+1) = impcon/delX*1._dp/6._dp*6._dp*(U(i+1)*A(i+1)) & !Advection
            -impcon/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX   ) !Diffusion

            !!Main diagonal
            band(KL+KU+1+i-i,i) = A(i)/DT + wset*wetwidth(i) & !Unsteady
            -impcon/delX*1._dp/6._dp*3._dp*(U(i)*A(i)) & !Advection
            +impcon/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX  ) & !Diffusion
            +impcon/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-diagonal
            band(KL+KU+1+i-(i-1),i-1)= -impcon/delX*1._dp/6._dp*2._dp*(U(i-1)*A(i-1)) & !Advection
            -impcon/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-lower-diagonal
            band(KL+KU+1+i-(i-2),i-2)= 0._dp !Advection

            rhs(i)= Qe(i) + Area_old(i)*C_old(i)/DT -(1._dp-impcon)/delX*( &
            -1._dp/6._dp*2._dp*(UU_old(i-1)*Area_old(i-1))*C_old(i-1) & !Advective
            - 1._dp/6._dp*3._dp*(UU_old(i)*Area_old(i))*C_old(i) & !Advective
            -0.5_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i))*(C_old(i+1)-C_old(i))/delX & !Diffusive
            +1._dp/6._dp*6._dp*(UU_old(i+1)*Area_old(i+1))*C_old(i+1) & !Advective
            -0.5_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1))*(C_old(i)-C_old(i-1))/delX  &!Diffusive
            - 1._dp/6._dp*1._dp*(UU_old(i+2)*Area_old(i+2))*C_old(i+2) ) !Advective
            ELSE
            !Use first order upwind
            
            !!Upper-upper diagonal
            band(KL+KU+1+i-(i+2),i+2)= 0._dp !-1._dp/delX*1._dp/6._dp*1._dp*(U(i+2)*A(i+2))!Advection

            !!Upper-diagonal
            band(KL+KU+1+i-(i+1),i+1) = impconU/delX*(U(i+1)*A(i+1)) & !Advection
            -impconU/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX   ) !Diffusion

            !!Main diagonal
            band(KL+KU+1+i-i,i) = A(i)/DT + wset*wetwidth(i) & !Unsteady
            -impconU/delX*(U(i)*A(i)) & !Advection
            +impconU/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX  ) & !Diffusion
            +impconU/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-diagonal
            band(KL+KU+1+i-(i-1),i-1)= 0._dp & !-1._dp/delX*1._dp/6._dp*2._dp*(U(i-1)*A(i-1)) & !Advection
            -impconU/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-lower-diagonal
            band(KL+KU+1+i-(i-2),i-2)= 0._dp !Advection

            rhs(i)= Qe(i) + Area_old(i)*C_old(i)/DT -(1._dp-impconU)/delX*( &
            - (UU_old(i)*Area_old(i))*C_old(i) & !Advective
            -0.5_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i))*(C_old(i+1)-C_old(i))/delX & !Diffusive
            +(UU_old(i+1)*Area_old(i+1))*C_old(i+1) & !Advective
            -0.5_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1))*(C_old(i)-C_old(i-1))/delX ) !Diffusive

            END IF
    END IF 

    END DO 

    !Boundary conditions - we say that C(1)=C(2)=C2(1) (=C(3)) if the flow is
    !going outward, and otherwise we set the values to the same as the boundary
    !conditions
    !DOWNSTREAM BOUNDARY
    IF(U2(1)<0._dp) THEN
    !This says that C(1)-C(3) =0
    i=1
    band(KL+KU+1+i-(i+2),i+2) = -1._dp !Coefficient of 1 for C(3)
    band(KL+KU+1+i- (i),i)=1._dp !Coefficient of 1 for C(1)

    !This says that C(2)-C(3) =0
    i=2
    band(KL+KU+1+i-(i+1),i+1) = -1._dp !Coefficient of 1 for C(3)
    band(KL+KU+1+i- (i),i)=1._dp !Coefficient of 1 for C(2)

    ELSE
    !This says that C(1)=Cmouth
    i=1
    band(KL+KU+1+i- (i),i)=1._dp !Coefficient of 1 for C(1)
    rhs(i)=Cmouth

    !This says that C(2)=Cmouth
    i=2
    band(KL+KU+1+i- (i),i)=1._dp !Coefficient of 1 for C(1)
    rhs(i)=Cmouth
    END IF

    !UPSTREAM BOUNDARY 
    IF(U2(n2)>0._dp) THEN
    i=n2+3
    band(KL+KU+1+i-(i-1),i-1)=-1._dp !Coefficient for C(n2+2)
    band(KL+KU+1+i-(i),i)=1._dp !Coefficient for C(n2+3)


    i=n2+4
    band(KL+KU+1+i-(i-2),i-2)=-1._dp !Coefficient for C(n2+2)
    band(KL+KU+1+i-(i),i)=1._dp !Coefficient for C(n2+4)

    ELSE
    i=n2+3
    band(KL+KU+1+i-(i),i)=1._dp !Coefficient for C(n2+3)
    rhs(n2+3)=Criver

    i=n2+4
    band(KL+KU+1+i-(i),i)=1._dp !Coefficient for C(n2+4)
    rhs(n2+4)=Criver

    END IF
    !!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!Solve

    !call DGTSV(n,1, lo(2:n), diag, up(1:n-1), r,n, info)
    call DGBSV(n2+4, 2, 2, 1, band, 2*KL+KU+1, IPV, rhs, n2+4, info)
    !print*, 'info=', info
    !The solver writes the new C to r- so we change it back.
    C2(1:n2)=rhs(3:n2+2)

    !print*, info
    !print*, maxval(diag), minval(diag), maxval(lo), minval(lo), maxval(up), minval(up)
    !print*, "vel", maxval(U), minval(U)
    do i=1, n2
            if((i>1).and.(i<n2)) THEN
                    if(C2(i)<0._dp) THEN
                    ii=int(sign(1._dp, U2(i)))
                    C2(i-ii)=C2(i-ii)-A(i)*C2(i)/A(i-ii)
                    C2(i)=0._dp
                    END IF
            else
                    if(i==1) THEN
                            if(C2(i)<0._dp) THEN
                            ii=-1
                            C2(i-ii)=C2(i-ii)-A(i)*C2(i)/A(i-ii)
                            C2(i)=0._dp
                            end if 
                    end if
                    if(i==n2) THEN
                            if(C2(i)<0._dp) THEN
                            ii=1
                            C2(i-ii)=C2(i-ii)-A(i)*C2(i)/A(i-ii)
                            C2(i)=0._dp
                            end if 
                    end if
            end if
    END DO


    IF(minval(C2)<0._dp) THEN 
             print*, "min Sediment conc <0", minval(C2),  maxval(C_old2), minval(C_old2), maxval(U2), minval(U2)
             print*, "..........minloc C is", minloc(C2), "............ depth is =",A2(minloc(C2))/wetwidth2(minloc(C2))
    !         print*, r
    !         print*, "..........."
             print*, C_old(minloc(C2)+2), C_old(minloc(C2)-1+2), C_old(minloc(C2)+1+2), maxval(Qe2), minval(Qe2),&
     maxval(Qd2), minval(Qd2)
         !print*, ">>>>>>>>>>>>>"
         !print*, C
    !         !print*, A !diag-(abs(up)+abs(lo))
    do i=1, n2
    IF(C2(i)<0._dp) THEN

        IF(abs(C2(i))<10._dp**(-10)) THEN !This could just be due to round off error in the matrix solver - fix it.
        C2(i)=0._dp
        ELSE
        print*, "C2 is < 0; violation is", abs(C2(i))
    !	stop
        END IF
    END IF	
    END DO
    !         
    !        stop
    END IF

    !DO i=1, n
    !IF(C(i)<0._dp) C(i)=0._dp
    !END DO


    DO i = 1, n2
    Flag=isnan(C2(i))
    if(Flag) THEN
           print*, "sedconc is NAN"
           print*, U2
           print*, "................"
           print*, A2
           print*, Qd2
           print*, Qe2
           
            

            stop
    END IF
    END DO

end subroutine susconc_up3
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine susconc_up32(n2,DT, A2, QH2, QH2_old, delX,C2, U2, Qe2,Qe2_old, Qd2, Cmouth, C_old2, Area_old2, UU_old2, & 
Criver, wset, wetwidth2,wetwidth_old2, D2,D2_old, third)
    INTEGER, intent(in):: n2
    REAL(dp), intent(in):: DT, delX, Qe2, Qe2_old, Qd2, Cmouth, A2, QH2, QH2_old, U2, C_old2, Area_old2, UU_old2,Criver,&
     wset, wetwidth2,wetwidth_old2, D2, D2_old
    REAL(dp), intent(in out):: C2
    LOGICAL, intent(in):: third


    !!Note - here Qe is different to elsewhere. 
    Dimension C2(n2), A2(n2), QH2(n2), QH2_old(n2), U2(n2), Qe2(n2), Qe2_old(n2), Qd2(n2),C_old2(n2),Area_old2(n2),UU_old2(n2),&
    wetwidth2(n2),wetwidth_old2(n2), D2(n2), D2_old(n2)

    INTEGER:: info, i, KL=2, KU=2, IPV(n2+4), ii !Note KL,KU are the number of upper and lower diagonals of the banded matrix
    REAL(dp):: r(n2) !Right hand side of equation
    logical:: flag
    !REAL(dp):: D(n) !dispersion constant-- actually not constant, need to fix.
    REAL(dp):: impcon=.50_dp, impconU=.00_dp
    !REAL(dp):: Qf1(n), Qf0(n), Af1(n), Af0(n), Df1(n), Df0(n)
    REAL(dp):: band(7,n2+4), rhs(n2+4)
    REAL(dp):: A(n2+4), QH(n2+4), C(n2+4), U(n2+4), Qe(n2+4), Qd(n2+4), C_old(n2+4), Area_old(n2+4), UU_old(n2+4), wetwidth(n2+4), &
    D(n2+4), D_old(n2+4)

    !!Define new variables with 'ghost'points on the edges, useful for implementing
    !boundary conditions
    A(3:(n2+2))=A2
    QH(3:(n2+2))=QH2
    C(3:(n2+2))=C2
    U(3:(n2+2))=U2
    Qe(3:(n2+2))=Qe2
    Qd(3:(n2+2))=Qd2
    C_old(3:(n2+2))=C_old2
    Area_old(3:(n2+2))=Area_old2
    UU_old(3:(n2+2))=UU_old2
    wetwidth(3:(n2+2))=wetwidth2
    D(3:(n2+2))=D2
    D_old(3:(n2+2))=D2_old

    !Boundary conditions
    A(1:2)=A2(1)
    A((n2+3):(n2+4)) = A2(n2)
    QH(1:2)=QH2(1)
    QH((n2+3):(n2+4)) = QH2(n2)
    Area_old(1:2)=Area_old2(1)
    Area_old((n2+3):(n2+4)) = Area_old2(n2)
    U(1:2)=U2(1)
    U((n2+3):(n2+4)) = U2(n2)
    UU_old(1:2)=U2(1)
    UU_old((n2+3):(n2+4)) = UU_old2(n2)
    Qe(1:2)=0._dp
    Qe((n2+3):(n2+4)) = 0._dp 
    Qd(1:2)=0._dp
    Qd((n2+3):(n2+4)) = 0._dp 
    wetwidth(1:2)=wetwidth2(1)
    wetwidth((n2+3):(n2+4)) = wetwidth2(n2)
    D(1:2)=D2(1)
    D((n2+3):(n2+4))=D2(n2)
    D_old(1:2)=D2_old(1)
    D_old((n2+3):(n2+4))=D2_old(n2)

    if(U2(1)>0._dp) THEN
    C_old(1:2)=Cmouth
    ELSE
    C_old(1:2)=C_old2(1)
    END IF

    if(U2(n2)>0._dp) THEN
    C_old((n2+3):(n2+4)) = C_old2(n2) !0._dp !Criver 
    ELSE
    C_old((n2+3):(n2+4)) = Criver 
    END IF

    !D= 2._dp!*abs(U)
    !C(n)=0.
    !!!
    !!##Solves (del AC/ del T) +  (del CQ / delX) =   d/dx (A D del C / delX) + (E-D)  

    !!Note that E and D are the total rate of erosion/deposition over each section. 
    !!Solved with a Crank-Nicholson Method, using a third order upwind
    !discretization of the advection term. 

    !Although this is a common reference, I
    !found Vreugdenhill (1989:59) a useful reference, and the way the code is
    !written reflects that. 
    !For the third order advection, I have used 'ghost' cells. There are 2
    !downstream, and 2 upstream
    !!At the mouth boundary, a given value is imposed if the flow is inward, and
    !otherwise a zero gradient condition is imposed.
    !! At the landward boundary, a zero gradient condition is enforced if there is
    !outflow. Otherwise, the river concentration is enforced, or used to
    !provide other boundary values, depending on the version of this code you have.  
    !An elementary reference, which treats the boundaries a bit differently to me,
    !is Karahan (2007). 

    !Flux form of third order advection. dF/dx = [ F(i+1/2)-F(i-1/2) ] / dx
    !with F[i+1/2] = (1/6)[ - F(i-1) +5F(i) +2F(i+1) ]     if (Q(i+1/2) > 0)
    !F[i+1/2] = (1/6)[ 2F(i) +5F(i+1) -F(i+2) ] if ( Q(i+1/2) <0 )

    !Predefine
    band=0._dp
    rhs=0._dp
    DO i=3,n2+2
    IF (0.5_dp*(QH(i)+QH(i+1))>0._dp) THEN 
            IF(0.5_dp*(QH(i)+QH(i-1))>0._dp) THEN
            !Q(i+1/2)> 0; Q(i-1/2) > 0
            
    !Note the unusual banded storage that lapack does
    !!Upper-upper diagonal
            !IF(sign((C_old(i+1)-C_old(i))*(C_old(i)-C_old(i-1)), 1._dp).gt.0._dp) THEN
            !IF((C_old(i).gt.0.0E-7_dp).and.(C_old(i)>0.5_dp*max(C_old(i+1),C_old(i-1)))) THEN
            band(KL+KU+1+i-(i+2),i+2)= 0._dp !Advection

            !!Upper-diagonal
            band(KL+KU+1+i-(i+1),i+1) = impcon/delX*1._dp/6._dp*2._dp*(QH(i+1)) & !Advection
            -impcon/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX   ) !Diffusion

            !!Main diagonal
            band(KL+KU+1+i-i,i) = A(i)/DT +wset*wetwidth(i) + & !Unsteady
            impcon/delX*1._dp/6._dp*3._dp*(QH(i)) & !Advection
            +impcon/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX  ) & !Diffusion
            +impcon/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-diagonal
            band(KL+KU+1+i-(i-1),i-1)= -impcon/delX*1._dp/6._dp*6._dp*(QH(i-1)) & !Advection
            -impcon/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-lower-diagonal
            band(KL+KU+1+i-(i-2),i-2)= impcon/delX*1._dp/6._dp*1._dp*(QH(i-2)) !Advection

            !Right hand side
            rhs(i)= Qe(i) + Area_old(i)*C_old(i)/DT -(1._dp-impcon)/delX*( &
            +1._dp/6._dp*2._dp*(QH(i+1))*C_old(i+1) & !Advective
            + 1._dp/6._dp*3._dp*(QH(i))*C_old(i) & !Advective
            -0.5_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i))*(C_old(i+1)-C_old(i))/delX & !Diffusive
            -1._dp/6._dp*6._dp*(QH(i-1))*C_old(i-1) & !Advective
            +0.5_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1))*(C_old(i)-C_old(i-1))/delX  &!Diffusive
            + 1._dp/6._dp*1._dp*(QH(i-2))*C_old(i-2) ) !Advective
            ELSE
            !Q(i+1/2)>0, Q(i-1/2) <=0        
     
            band(KL+KU+1+i-(i+2),i+2)= 0._dp !Advection

            !!Upper-diagonal
            band(KL+KU+1+i-(i+1),i+1) = impcon/delX*1._dp/6._dp*3._dp*(QH(i+1)) & !Advection
            -impcon/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX   ) !Diffusion

            !!Main diagonal
            band(KL+KU+1+i-i,i) = A(i)/DT +wset*wetwidth(i) - & !Unsteady
            impcon/delX*1._dp/6._dp*0._dp*(QH(i)) & !Advection
            +impcon/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX  ) & !Diffusion
            +impcon/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-diagonal
            band(KL+KU+1+i-(i-1),i-1)= -impcon/delX*1._dp/6._dp*3._dp*(QH(i-1)) & !Advection
            -impcon/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-lower-diagonal
            band(KL+KU+1+i-(i-2),i-2)= 0._dp !impcon/delX*1._dp/6._dp*1._dp*(U(i-2)*A(i-2)) !Advection

            !Right hand side
            rhs(i)= Qe(i) + Area_old(i)*C_old(i)/DT -(1._dp-impcon)/delX*( &
            +1._dp/6._dp*3._dp*(QH(i+1))*C_old(i+1) & !Advective
            - 1._dp/6._dp*0._dp*(QH(i))*C_old(i) & !Advective
            -0.5_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i))*(C_old(i+1)-C_old(i))/delX & !Diffusive
            -1._dp/6._dp*3._dp*(QH(i-1))*C_old(i-1) & !Advective
            +0.5_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1))*(C_old(i)-C_old(i-1))/delX  &!Diffusive
            + 0._dp )!1._dp/6._dp*1._dp*(UU_old(i-2)*Area_old(i-2))*C_old(i-2) ) !Advective


            END IF

    ELSE

    !Note the unusual banded storage that lapack does
            !IF((C_old(i).gt.0.0E-7_dp).and.(C_old(i)>0.5_dp*max(C_old(i+1),C_old(i-1)))) THEN
            !IF(sign((C_old(i+1)-C_old(i))*(C_old(i)-C_old(i-1)), 1._dp).gt.0._dp) THEN
            !if(third.eqv..true.) THEN 
            if(0.5_dp*(QH(i)+QH(i-1))<0._dp) THEN
            !Q(i+1/2)<0 ; Q(i-1/2) <0  
            !!Upper-upper diagonal
            band(KL+KU+1+i-(i+2),i+2)=  -impcon/delX*1._dp/6._dp*1._dp*(QH(i+2))!Advection

            !!Upper-diagonal
            band(KL+KU+1+i-(i+1),i+1) = impcon/delX*1._dp/6._dp*6._dp*(QH(i+1)) & !Advection
            -impcon/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX   ) !Diffusion

            !!Main diagonal
            band(KL+KU+1+i-i,i) = A(i)/DT + wset*wetwidth(i) & !Unsteady
            -impcon/delX*1._dp/6._dp*3._dp*(QH(i)) & !Advection
            +impcon/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX  ) & !Diffusion
            +impcon/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-diagonal
            band(KL+KU+1+i-(i-1),i-1)= -impcon/delX*1._dp/6._dp*2._dp*(QH(i-1)) & !Advection
            -impcon/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-lower-diagonal
            band(KL+KU+1+i-(i-2),i-2)= 0._dp !Advection

            rhs(i)= Qe(i) + Area_old(i)*C_old(i)/DT -(1._dp-impcon)/delX*( &
            -1._dp/6._dp*2._dp*(QH(i-1))*C_old(i-1) & !Advective
            - 1._dp/6._dp*3._dp*(QH(i))*C_old(i) & !Advective
            -0.5_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i))*(C_old(i+1)-C_old(i))/delX & !Diffusive
            +1._dp/6._dp*6._dp*(QH(i+1))*C_old(i+1) & !Advective
            +0.5_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1))*(C_old(i)-C_old(i-1))/delX  &!Diffusive
            - 1._dp/6._dp*1._dp*(QH(i+2))*C_old(i+2) ) !Advective
            ELSE
            !Q(i+1/2)<0 ; Q(i-1/2) >=0  
            !!Upper-upper diagonal
            band(KL+KU+1+i-(i+2),i+2)=  -impcon/delX*1._dp/6._dp*1._dp*(QH(i+2))!Advection

            !!Upper-diagonal
            band(KL+KU+1+i-(i+1),i+1) = impcon/delX*1._dp/6._dp*5._dp*(QH(i+1)) & !Advection
            -impcon/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX   ) !Diffusion

            !!Main diagonal
            band(KL+KU+1+i-i,i) = A(i)/DT + wset*wetwidth(i) & !Unsteady
            -0._dp & !Advection
            +impcon/delX*( 0.5_dp*(A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX  ) & !Diffusion
            +impcon/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-diagonal
            band(KL+KU+1+i-(i-1),i-1)= -impcon/delX*1._dp/6._dp*5._dp*(QH(i-1)) & !Advection
            -impcon/delX*( 0.5_dp*(A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX  ) !Diffusion

            !!Lower-lower-diagonal
            band(KL+KU+1+i-(i-2),i-2)= +impcon/delX*1._dp/6._dp*1._dp*(QH(i-2))  !Advection

            rhs(i)= Qe(i) + Area_old(i)*C_old(i)/DT -(1._dp-impcon)/delX*( &
            +1._dp/6._dp*1._dp*(QH(i-2))*C_old(i-2) & !Advection
            -1._dp/6._dp*5._dp*(QH(i-1))*C_old(i-1) & !Advective
            -0.5_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i))*(C_old(i+1)-C_old(i))/delX & !Diffusive
            +1._dp/6._dp*5._dp*(QH(i+1))*C_old(i+1) & !Advective
            +0.5_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1))*(C_old(i)-C_old(i-1))/delX  &!Diffusive
            - 1._dp/6._dp*1._dp*(QH(i+2))*C_old(i+2) ) !Advective
            END IF
    END IF 

    END DO 

    !Boundary conditions - we say that C(1)=C(2)=C2(1) (=C(3)) if the flow is
    !going outward, and otherwise we set the values to the same as the boundary
    !conditions
    !DOWNSTREAM BOUNDARY
    IF(U2(1)<0._dp) THEN
    !This says that C(1)-C(3) =0
    i=1
    band(KL+KU+1+i-(i+2),i+2) = -1._dp !Coefficient of 1 for C(3)
    band(KL+KU+1+i- (i),i)=1._dp !Coefficient of 1 for C(1)

    !This says that C(2)-C(3) =0
    i=2
    band(KL+KU+1+i-(i+1),i+1) = -1._dp !Coefficient of 1 for C(3)
    band(KL+KU+1+i- (i),i)=1._dp !Coefficient of 1 for C(2)

    ELSE
    !This says that C(1)=Cmouth
    i=1
    band(KL+KU+1+i- (i),i)=1._dp !Coefficient of 1 for C(1)
    rhs(i)=Cmouth

    !This says that C(2)=Cmouth
    i=2
    band(KL+KU+1+i- (i),i)=1._dp !Coefficient of 1 for C(1)
    rhs(i)=Cmouth
    END IF

    !UPSTREAM BOUNDARY 
    IF(U2(n2)>0._dp) THEN
    i=n2+3
    band(KL+KU+1+i-(i-1),i-1)=-1._dp !Coefficient for C(n2+2)
    band(KL+KU+1+i-(i),i)=1._dp !Coefficient for C(n2+3)


    i=n2+4
    band(KL+KU+1+i-(i-2),i-2)=-1._dp !Coefficient for C(n2+2)
    band(KL+KU+1+i-(i),i)=1._dp !Coefficient for C(n2+4)

    ELSE
    i=n2+3
    band(KL+KU+1+i-(i),i)=1._dp !Coefficient for C(n2+3)
    rhs(n2+3)=Criver

    i=n2+4
    band(KL+KU+1+i-(i),i)=1._dp !Coefficient for C(n2+4)
    rhs(n2+4)=Criver

    END IF
    !!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!Solve

    !call DGTSV(n,1, lo(2:n), diag, up(1:n-1), r, n,info)
    call DGBSV(n2+4, 2, 2, 1, band, 2*KL+KU+1, IPV, rhs, n2+4, info)
    !print*, 'info=', info
    !The solver writes the new C to r- so we change it back.
    C2(1:n2)=rhs(3:n2+2)

    !print*, info
    !print*, maxval(diag), minval(diag), maxval(lo), minval(lo), maxval(up), minval(up)
    !print*, "vel", maxval(U), minval(U)
    do i=1, n2
            if((i>1).and.(i<n2)) THEN
                    if(C2(i)<0._dp) THEN
                    ii=int(sign(1._dp, U2(i)))
                    C2(i-ii)=C2(i-ii)-A(i)*C2(i)/A(i-ii)
                    C2(i)=0._dp
                    END IF
            else
                    if(i==1) THEN
                            if(C2(i)<0._dp) THEN
                            ii=-1
                            C2(i-ii)=C2(i-ii)-A(i)*C2(i)/A(i-ii)
                            C2(i)=0._dp
                            end if 
                    end if
                    if(i==n2) THEN
                            if(C2(i)<0._dp) THEN
                            ii=1
                            C2(i-ii)=C2(i-ii)-A(i)*C2(i)/A(i-ii)
                            C2(i)=0._dp
                            end if 
                    end if
            end if
    END DO


    IF(minval(C2)<0._dp) THEN 
             print*, "min Sediment conc <0", minval(C2),  maxval(C_old2), minval(C_old2), maxval(U2), minval(U2)
             print*, "..........minloc C is", minloc(C2), "............ depth is =",A2(minloc(C2))/wetwidth2(minloc(C2))
    !         print*, r
    !         print*, "..........."
             print*, C_old(minloc(C2)+2), C_old(minloc(C2)-1+2), C_old(minloc(C2)+1+2), maxval(Qe2), minval(Qe2),&
     maxval(Qd2), minval(Qd2)
         !print*, ">>>>>>>>>>>>>"
         !print*, C
    !         !print*, A !diag-(abs(up)+abs(lo))
    do i=1, n2
    IF(C2(i)<0._dp) THEN

        IF(abs(C2(i))<10._dp**(-10)) THEN !This could just be due to round off error in the matrix solver - fix it.
        C2(i)=0._dp
        ELSE
        print*, "C2 is < 0; violation is", abs(C2(i))
    !	stop
        END IF
    END IF	
    END DO
    !         
    !        stop
    END IF

    !DO i=1, n
    !IF(C(i)<0._dp) C(i)=0._dp
    !END DO


    DO i = 1, n2
    Flag=isnan(C2(i))
    if(Flag) THEN
           print*, "sedconc is NAN"
           print*, U2
           print*, "................"
           print*, A2
           print*, Qd2
           print*, Qe2
           
            

            stop
    END IF
    END DO

end subroutine susconc_up32
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine susconc_up33(n2,DT, A2, QH2, QH2_old, delX,C2, U2, Qe2,Qe2_old, Qd2, Cmouth, C_old2, Area_old2, UU_old2, & 
Criver, wset, wetwidth2,wetwidth_old2, D2,D2_old, third)
    INTEGER, intent(in):: n2
    REAL(dp), intent(in):: DT, delX, Qe2, Qe2_old, Qd2, Cmouth, A2, QH2, QH2_old, U2, C_old2, Area_old2, UU_old2,Criver,&
     wset, wetwidth2,wetwidth_old2, D2, D2_old
    REAL(dp), intent(in out):: C2
    LOGICAL, intent(in):: third


    !!Note - here Qe is different to elsewhere. 
    Dimension C2(n2), A2(n2), QH2(n2), QH2_old(n2), U2(n2), Qe2(n2), Qe2_old(n2), Qd2(n2),C_old2(n2),Area_old2(n2),UU_old2(n2),&
    wetwidth2(n2),wetwidth_old2(n2), D2(n2), D2_old(n2)

    INTEGER:: info, i, KL=2, KU=2, IPV(n2+4), ii !Note KL,KU are the number of upper and lower diagonals of the banded matrix
    REAL(dp):: r(n2) !Right hand side of equation
    logical:: flag
    !REAL(dp):: D(n) !dispersion constant-- actually not constant, need to fix.
    REAL(dp):: impcon=.50_dp, impconU=.00_dp
    !REAL(dp):: Qf1(n), Qf0(n), Af1(n), Af0(n), Df1(n), Df0(n)
    REAL(dp):: band(7,n2+4), rhs(n2+4)
    REAL(dp):: A(n2+4), QH(n2+4), QH_old(n2+4), C(n2+4), U(n2+4), Qe(n2+4),Qe_old(n2+4), Qd(n2+4), C_old(n2+4), Area_old(n2+4), &
     UU_old(n2+4), wetwidth(n2+4),wetwidth_old(n2+4), D(n2+4), D_old(n2+4), Fl1(n2+4), Q_old(n2+4), limi(n2+4), FL_old(n2+4), & 
    theta(n2+4), Cpred(n2+4)
    REAL(dp):: mu=1._dp, eeps=1.0E-10_dp
    !!##Solves (del AC/ del T) +  (del CQ / delX) =   d/dx (A Diffuse del C / delX) + (E-D)  

    !!Note that E and D are the total rate of erosion/deposition over each section. 
    !Solved Using an approach developed from Hundsdorfer's notes, probably also in their book. 
    !Note that it is easy to have this explicit (stability requirement is much less
    !stringent than for 1D St Venant -- this has a stability determined by the velocity and delX only, no gravity waves)

    !!The numerical time stepping scheme is described by Hundsdorfer's lecture notes
    !(page 49) as 'the implicit midpoint rule with Euler predictor', although it is
    !not implicit. In the book Hundsdorfer and Verwer (2003) it is called the 'one
    !step explicit midpoint rule' (page 142)
    !Basically we take a predictor half-step
    !Cpred = Clast + (1/2 delt)/delx *F(tlast,Clast)
    !And then a corrector full step
    !C=Clast+ delt/delX * F(t+1/2 delT, Cpred)
    !On the first step, the discharge is evaluated as QH2, which is
    !the conservative discharge estimates from the McCormack method. Note
    !that this is a reasonable estimate of the discharge at the old time level if we
    !assume that the discharge is constant between the last step and the next one -
    !While this is a crude assumption, more complex things I tried (e.g Q =
    !0.5*(QH2+QH2_old) lead to greater overshoot near the mouth than the present
    !approach. (NOT UP TO DATE - PRESENTLY I AM USING THE LATTER - HOPEFULLY FINE).
    !Note that while we could just use the straight discharge output from McCormack, it would not be
    !conservative. I have checked this by running a constant discharge case with no
    !erosion or deposition - the present algorithm predicts a constant sediment
    !concentration (and the QH2 discharge is constant), while if we use the straight
    !output from McCormack, then both the Discharge and sediment conc show slight
    !variation from constant.
    ! The diffusion coefs are evaluated at the old time level, as is the rate of erosion. 
    !On the second step, the discharge is evaluated as QH (a good half-time step
    !estimate), as are the rate of erosion and the diffusion coef, by suitable
    !averaging of the input variables 

    !The spatial discretization is a standard central scheme for diffusion, and a
    !flux limited third order upwind scheme for advection. The latter is described
    !in Hundsdorfer's notes (page 38)

    !Note -- the discharge issues. With the McCormack Scheme, we can show that
    !A(t+1)-A(t) = (dT/dX)*(0.5*(Q(t)_[i+1]+Q(Pred)_[i])-0.5*(Q(t)_[i]+Q(Pred)_[i-1]) )
    !This means that at steady state (constant discharge), it is actually (Q(t)_[i+1]+Q(Pred)_[i]) which
    !will not be changing either in time or in space.

    !Although this is a common reference, I
    !found Vreugdenhill (1989:59) a useful reference, and the way the code USED TO be
    !written reflects that. 


    !Flux form of third order advection. dF/dx = [ F(i+1/2)-F(i-1/2) ] / dx
    !with F[i+1/2] = (1/6)[ - F(i-1) +5F(i) +2F(i+1) ]     if (Q(i+1/2) > 0)
    !F[i+1/2] = (1/6)[ 2F(i) +5F(i+1) -F(i+2) ] if ( Q(i+1/2) <0 )



    !For the third order advection, I have used 'ghost' cells. There are 2
    !downstream, and 2 upstream
    !!At the mouth boundary, a given value is imposed if the flow is inward, and
    !otherwise a zero gradient condition is imposed.
    !! At the landward boundary, a zero gradient condition is enforced if there is
    !outflow. Otherwise, the river concentration is enforced, or used to
    !provide other boundary values, depending on the version of this code you have.  

    !Note that if we force advection to first order, presently I can find weird
    !long-term behaviour - accumulation of sediment in upstream zones -odd.

    !!Define new variables with 'ghost'points on the edges, useful for implementing
    !boundary conditions
    A(3:(n2+2))=A2
    QH(3:(n2+2))=QH2
    QH_old(3:(n2+2))=QH2_old
    C(3:(n2+2))=C2
    U(3:(n2+2))=U2
    Qe(3:(n2+2))=Qe2
    Qe_old(3:(n2+2))=Qe2_old
    Qd(3:(n2+2))=Qd2
    C_old(3:(n2+2))=C_old2
    Area_old(3:(n2+2))=Area_old2
    UU_old(3:(n2+2))=UU_old2
    wetwidth(3:(n2+2))=wetwidth2
    wetwidth_old(3:(n2+2))=wetwidth_old2
    D(3:(n2+2))=D2
    D_old(3:(n2+2))=D2_old

    !Boundary conditions
    A(1:2)=A2(1)
    A((n2+3):(n2+4)) = A2(n2)
    QH(1:2)=QH2(1)
    QH((n2+3):(n2+4)) = QH2(n2)
    QH_old(1:2)=QH2_old(1)
    QH_old((n2+3):(n2+4)) = QH2_old(n2)
    Area_old(1:2)=Area_old2(1)
    Area_old((n2+3):(n2+4)) = Area_old2(n2)
    U(1:2)=U2(1)
    U((n2+3):(n2+4)) = U2(n2)
    UU_old(1:2)=U2(1)
    UU_old((n2+3):(n2+4)) = UU_old2(n2)
    Qe(1:2)=0._dp
    Qe((n2+3):(n2+4)) = 0._dp 
    Qe_old(1:2)=0._dp
    Qe_old((n2+3):(n2+4)) = 0._dp 
    Qd(1:2)=0._dp
    Qd((n2+3):(n2+4)) = 0._dp 
    wetwidth(1:2)=wetwidth2(1)
    wetwidth((n2+3):(n2+4)) = wetwidth2(n2)
    wetwidth_old(1:2)=wetwidth_old2(1)
    wetwidth_old((n2+3):(n2+4)) = wetwidth_old2(n2)
    D(1:2)=D2(1)
    D((n2+3):(n2+4))=D2(n2)
    D_old(1:2)=D2_old(1)
    D_old((n2+3):(n2+4))=D2_old(n2)

    if(U2(1)>0._dp) THEN
    C_old(1)=Cmouth
    C_old(2)=Cmouth !0.5_dp*(Cmouth+C_old(3))
    ELSE
    C_old(2)=max(2._dp*C_old(3)-C_old(4),0._dp)
    C_old(1)=max(2._dp*C_old(2) - C_old(3),0._dp )
    END IF

    if(U2(n2)>0._dp) THEN
    C_old((n2+3):(n2+4)) = C_old2(n2) !0._dp !Criver 
    ELSE
    C_old((n2+3):(n2+4)) = Criver 
    END IF


    !The new flux
    FL1=0._dp
    !The old Q
    Q_old= 0.5_dp*(QH+QH_old )!This is like the discharge at the beginning of the time step 
    !The old Flux
    FL_old=Q_old*C_old
    !Useful for limiter, Hundsdorfer
    do i=2,n2+3
        if((FL_old(i+1).ne.FL_old(i)).and.((FL_old(i).ne.FL_old(i-1)))) THEN
            theta(i)=(FL_old(i)-FL_old(i-1))/(FL_old(i+1)-FL_old(i))  
        ELSE
            theta(i)=(FL_old(i)-FL_old(i-1)+eeps)/(FL_old(i+1)-FL_old(i) +eeps)  
        END IF
    end do

    !Calculate fluxes, with limiting -- 3rd order without limiting
    do i=2, n2+2
            if(0.5_dp*(Q_old(i)+Q_old(i+1))>=0._dp) THEN
            
            limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/6._dp*theta(i) , mu*theta(i)))
            !if((i>n2).or.(i<4).or.(.false.)) limi(i)=0._dp !Force the boundaries to use an upwind flux
            FL1(i)= FL_old(i)+limi(i)*(FL_old(i+1)-FL_old(i)) !(1._dp/6._dp)*( 2._dp*Q_old(i+1)*C_old(i+1) + 5._dp*Q_old(i)*C_old(i) -Q_old(i-1)*C_old(i-1)) 
            ELSE
            limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/(6._dp*theta(i+1)) , mu/theta(i+1)))
            !if((i>n2).or.(i<4).or.(.false.)) limi(i)=0._dp !Force the boundaries to use an upwind flux
            FL1(i)= ( FL_old(i+1)+limi(i)*(FL_old(i)-FL_old(i+1)) ) ! (1._dp/6._dp)*( 2._dp*Q_old(i)*C_old(i) + 5._dp*Q_old(i+1)*C_old(i+1) -Q_old(i+2)*C_old(i+2)) 
            END IF
    end do

    !Calculate C at the next half time step, predictor style
    do i=3, n2+2
    !C(i) = (0.5_dp*DT)/(0.5_dp*(A(i)+A(i)))*(Area_old(i)*C_old(i)/(DT*0.5_dp) + Qe_old(i) -wset*C_old(i)*wetwidth_old(i) & !Source Terms
    C(i) = ((0.5_dp*(A(i)+Area_old(i)))/(0.5_dp*DT)+wset*0.5_dp*(wetwidth_old(i)+wetwidth(i)))**(-1._dp)*&  !Constant of C(i), including implicit settling
    (Area_old(i)*C_old(i)/(DT*0.5_dp) &  !Unsteady
    + Qe(i) & !Erosion at half time (Source Terms)
    - 1._dp/delX*( (FL1(i)-FL1(i-1)) - &  !Advection
    (0.5_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i))*(C_old(i+1)-C_old(i))/delX - & !Diffusion
    0.5_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1))*(C_old(i)-C_old(i-1))/delX )) ) !Diffusion
    end do

    !Mouth boundary conditions for the half step -- we impose zero gradient if the flow is
    !outward, and a mouth boundary condition otherwise
    IF(Q_old(3)>0._dp) THEN
    C(1)=Cmouth
    C(2)=Cmouth !0.5_dp*(Cmouth+C(3))
    ELSE
    C(2)=max(2._dp*C(3) -C(4), 0._dp)
    C(1)=max(2._dp*C(2)-C(3), 0._dp)
    END IF

    !River boundary conditions for the half step -- we impose zero gradient if the flow is
    !outward, and a river boundary condition otherwise
    IF(Q_old(n2+2)>0._dp) THEN
    C(n2+3)=C(n2+2)
    C(n2+4)=C(n2+2)
    ELSE
    C(n2+3)=Criver
    C(n2+4)=Criver
    END IF

    Cpred=C

    !write(23,*) Cpred
    !!Second half step --reuse lots of code



    !The second half Flux - evaluated at centred position - This is based on the
    !lecture notes from Hundsdorfer, page 
    Fl_old=QH*Cpred
    !Useful for limiter
    do i=2,n2+3
        if((FL_old(i+1).ne.FL_old(i)).and.((FL_old(i).ne.FL_old(i-1)))) THEN
            theta(i)=(FL_old(i)-FL_old(i-1))/(FL_old(i+1)-FL_old(i))  
        ELSE
            theta(i)=(FL_old(i)-FL_old(i-1)+eeps)/(FL_old(i+1)-FL_old(i) +eeps)  
        END IF

        if(theta(i).eq.0._dp) print*, 'zero theta'
        if(isnan(theta(i))) print*, 'nan theta'
    end do

    !Calculate fluxes, with limiting -- 3rd order without limiting
    do i=2, n2+2
            IF(0.5_dp*(QH(i)+QH(i+1))>=0._dp) THEN
            limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/6._dp*theta(i) , mu*theta(i)))
            !if((i>n2).or.(i<4).or.(.false.)) limi(i)=0._dp !Force the boundaries to use an upwind flux
            FL1(i)= FL_old(i)+limi(i)*(FL_old(i+1)-FL_old(i)) 
            ELSE
            limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/(6._dp*theta(i+1)) , mu/theta(i+1)))
            !if((i>n2).or.(i<4).or.(.false.)) limi(i)=0._dp !Force the boundaries to use an upwind flux
            FL1(i)= ( FL_old(i+1)+limi(i)*(FL_old(i)-FL_old(i+1)) )  
            END IF
    end do

    !Calculate C at the next time step, corrector style
    do i=3, n2+2
    C(i) = DT/A(i)*(Area_old(i)*C_old(i)/DT + Qe(i) -wset*Cpred(i)*0.5_dp*(wetwidth(i) +wetwidth_old(i)) & !Source Terms
    - 1._dp/delX*( (FL1(i)-FL1(i-1)) - &  !Advection
    (0.25_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i)+A(i+1)*D(i+1)+A(i)*D(i))*(Cpred(i+1)-Cpred(i))/delX - & !Diffusion
    0.25_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1)+ A(i)*D(i)+A(i-1)*D(i-1))*(Cpred(i)-Cpred(i-1))/delX )) ) !Diffusion
    end do

    !Update C2 - C2 holds the sediment concentration.
    C2(1:n2)=C(3:n2+2) 

    !Hundsdorfer states that the use of eeps in theta (for the limiter) can allow
    !negative values of order eps. Let's get rid of them 
    do i=1, n2
            if((i>1).and.(i<n2)) THEN
                    if(C2(i)<0._dp) THEN
                    ii=int(sign(1._dp, U2(i)))
                    C2(i-ii)=C2(i-ii)-A(i)*C2(i)/A(i-ii)
                    C2(i)=0._dp
                    END IF
            else
                    if(i==1) THEN
                            if(C2(i)<0._dp) THEN
                            ii=-1
                            C2(i-ii)=C2(i-ii)-A(i)*C2(i)/A(i-ii)
                            C2(i)=0._dp
                            end if 
                    end if
                    if(i==n2) THEN
                            if(C2(i)<0._dp) THEN
                            ii=1
                            C2(i-ii)=C2(i-ii)-A(i)*C2(i)/A(i-ii)
                            C2(i)=0._dp
                            end if 
                    end if
            end if
    END DO


    IF(minval(C2)<0._dp) THEN 
             print*, "min Sediment conc <0", minval(C2),  maxval(C_old2), minval(C_old2), maxval(U2), minval(U2)
             print*, "..........minloc C is", minloc(C2), "............ depth is =",A2(minloc(C2))/wetwidth2(minloc(C2))
    !         print*, r
    !         print*, "..........."
             print*, C_old(minloc(C2)+2), C_old(minloc(C2)-1+2), C_old(minloc(C2)+1+2), maxval(Qe2), minval(Qe2),&
     maxval(Qd2), minval(Qd2)
         !print*, ">>>>>>>>>>>>>"
         !print*, C
    !         !print*, A !diag-(abs(up)+abs(lo))
    do i=1, n2
    IF(C2(i)<0._dp) THEN

        IF(abs(C2(i))<10._dp**(-10)) THEN !This could just be due to round off error in the matrix solver - fix it.
        C2(i)=0._dp
        ELSE
        print*, "C2 is < 0; violation is", abs(C2(i))
    !	stop
        END IF
    END IF	
    END DO
    !         
    !        stop
    END IF

    !DO i=1, n
    !IF(C(i)<0._dp) C(i)=0._dp
    !END DO


    DO i = 1, n2
    Flag=isnan(C2(i))
    if(Flag) THEN
           print*, "sedconc is NAN"
           print*, U2
           print*, "................"
           print*, A2
           print*, Qd2
           print*, Qe2
           
            

            stop
    END IF
    END DO

end subroutine susconc_up33
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine susconc_up34(n2,DT, A2, QH2, QH2_old, delX,C2, U2, Qe2,Qe2_old, Qd2, Cmouth, C_old2, Area_old2, UU_old2, & 
Criver, wset, wetwidth2,wetwidth_old2, D2,D2_old, third)
    INTEGER, intent(in):: n2
    REAL(dp), intent(in):: DT, delX, Qe2, Qe2_old, Qd2, Cmouth, A2, QH2, QH2_old, U2, C_old2, Area_old2, UU_old2,Criver,&
     wset, wetwidth2,wetwidth_old2, D2, D2_old
    REAL(dp), intent(in out):: C2
    LOGICAL, intent(in):: third


    !!Note - here Qe is different to elsewhere. 
    Dimension C2(n2), A2(n2), QH2(n2), QH2_old(n2), U2(n2), Qe2(n2), Qe2_old(n2), Qd2(n2),C_old2(n2),Area_old2(n2),UU_old2(n2),&
    wetwidth2(n2),wetwidth_old2(n2), D2(n2), D2_old(n2)

    INTEGER:: info, i, KL=2, KU=2, IPV(n2+4), ii !Note KL,KU are the number of upper and lower diagonals of the banded matrix
    REAL(dp):: r(n2) !Right hand side of equation
    logical:: flag
    !REAL(dp):: D(n) !dispersion constant-- actually not constant, need to fix.
    REAL(dp):: impcon=.50_dp, impconU=.00_dp
    !REAL(dp):: Qf1(n), Qf0(n), Af1(n), Af0(n), Df1(n), Df0(n)
    REAL(dp):: band(7,n2+4), rhs(n2+4)
    REAL(dp):: A(n2+4), QH(n2+4), QH_old(n2+4), C(n2+4), U(n2+4), Qe(n2+4),Qe_old(n2+4), Qd(n2+4), C_old(n2+4), Area_old(n2+4), &
     UU_old(n2+4), wetwidth(n2+4), D(n2+4), D_old(n2+4), Fl1(n2+4), Q_old(n2+4), limi(n2+4), FL_old(n2+4), & 
    theta(n2+4), Cpred(n2+4)
    REAL(dp):: mu=1._dp, eeps=1.0E-10_dp
    !!##Solves (del AC/ del T) +  (del CQ / delX) =   d/dx (A Diffuse del C / delX) + (E-D)  

    !!Note that E and D are the total rate of erosion/deposition over each section. 
    !Solved Using an approach developed from Hundsdorfer's notes, probably also in their book. 
    !Note that it is easy to have this explicit (stability requirement is much less
    !stringent than for 1D St Venant -- this has a stability determined by the velocity and delX only, no gravity waves)

    !!The numerical time stepping scheme is described by Hundsdorfer's lecture notes
    !(page 49) as 'the implicit midpoint rule with Euler predictor', although it is
    !not implicit. Basically we take a predictor half-step
    !Cpred = Clast + (1/2 delt)/delx *F(tlast,Clast)
    !And then a corrector full step
    !C=Clast+ delt/delX * F(t(1/2), Cpred)
    !On the first step, the discharge is evaluated as QH2, which is
    !the conservative discharge estimates from the McCormack method. Note
    !that this is a reasonable estimate of the discharge at the old time level if we
    !assume that the discharge is constant between the last step and the next one -
    !While this is a crude assumption, more complex things I tried (e.g Q =
    !0.5*(QH2+QH2_old) lead to greater overshoot near the mouth than the present
    !approach.
    !Note that while we could just use the straight discharge output from McCormack, it would not be
    !conservative. I have checked this by running a constant discharge case with no
    !erosion or deposition - the present algorithm predicts a constant sediment
    !concentration (and the QH2 discharge is constant), while if we use the straight
    !output from McCormack, then both the Discharge and sediment conc show slight
    !variation from constant.
    ! The diffusion coefs are evaluated at the old time level, as is the rate of erosion. 
    !On the second step, the discharge is evaluated as QH (a good half-time step
    !estimate), as are the rate of erosion and the diffusion coef, by suitable
    !averaging of the input variables 

    !The spatial discretization is a standard central scheme for diffusion, and a
    !flux limited third order upwind scheme for advection. The latter is described
    !in Hundsdorfer's notes (page 38)

    !Note -- the discharge issues. With the McCormack Scheme, we can show that
    !A(t+1)-A(t) = (dT/dX)*(0.5*(Q(t)_[i+1]+Q(Pred)_[i])-0.5*(Q(t)_[i]+Q(Pred)_[i-1]) )
    !This means that at steady state (constant discharge), it is actually (Q(t)_[i+1]+Q(Pred)_[i]) which
    !will not be changing either in time or in space.

    !Although this is a common reference, I
    !found Vreugdenhill (1989:59) a useful reference, and the way the code USED TO be
    !written reflects that. 


    !Flux form of third order advection. dF/dx = [ F(i+1/2)-F(i-1/2) ] / dx
    !with F[i+1/2] = (1/6)[ - F(i-1) +5F(i) +2F(i+1) ]     if (Q(i+1/2) > 0)
    !F[i+1/2] = (1/6)[ 2F(i) +5F(i+1) -F(i+2) ] if ( Q(i+1/2) <0 )



    !For the third order advection, I have used 'ghost' cells. There are 2
    !downstream, and 2 upstream
    !!At the mouth boundary, a given value is imposed if the flow is inward, and
    !otherwise a zero gradient condition is imposed.
    !! At the landward boundary, a zero gradient condition is enforced if there is
    !outflow. Otherwise, the river concentration is enforced, or used to
    !provide other boundary values, depending on the version of this code you have.  

    !Note that if we force advection to first order, presently I can find weird
    !long-term behaviour - accumulation of sediment in upstream zones -odd.

    !!Define new variables with 'ghost'points on the edges, useful for implementing
    !boundary conditions
    A(3:(n2+2))=A2
    QH(3:(n2+2))=QH2
    QH_old(3:(n2+2))=QH2_old
    C(3:(n2+2))=C2
    U(3:(n2+2))=U2
    Qe(3:(n2+2))=Qe2
    Qe_old(3:(n2+2))=Qe2_old
    Qd(3:(n2+2))=Qd2
    C_old(3:(n2+2))=C_old2
    Area_old(3:(n2+2))=Area_old2
    UU_old(3:(n2+2))=UU_old2
    wetwidth(3:(n2+2))=wetwidth2
    D(3:(n2+2))=D2
    D_old(3:(n2+2))=D2_old

    !Boundary conditions
    A(1:2)=A2(1)
    A((n2+3):(n2+4)) = A2(n2)
    QH(1:2)=QH2(1)
    QH((n2+3):(n2+4)) = QH2(n2)
    QH_old(1:2)=QH2_old(1)
    QH_old((n2+3):(n2+4)) = QH2_old(n2)
    Area_old(1:2)=Area_old2(1)
    Area_old((n2+3):(n2+4)) = Area_old2(n2)
    U(1:2)=U2(1)
    U((n2+3):(n2+4)) = U2(n2)
    UU_old(1:2)=U2(1)
    UU_old((n2+3):(n2+4)) = UU_old2(n2)
    Qe(1:2)=0._dp
    Qe((n2+3):(n2+4)) = 0._dp 
    Qe_old(1:2)=0._dp
    Qe_old((n2+3):(n2+4)) = 0._dp 
    Qd(1:2)=0._dp
    Qd((n2+3):(n2+4)) = 0._dp 
    wetwidth(1:2)=wetwidth2(1)
    wetwidth((n2+3):(n2+4)) = wetwidth2(n2)
    D(1:2)=D2(1)
    D((n2+3):(n2+4))=D2(n2)
    D_old(1:2)=D2_old(1)
    D_old((n2+3):(n2+4))=D2_old(n2)

    if(U2(1)>0._dp) THEN
    C_old(1:2)=Cmouth
    ELSE
    C_old(1:2)=Cmouth !C_old2(1)
    END IF

    if(U2(n2)>0._dp) THEN
    C_old((n2+3):(n2+4)) = Criver !C_old2(n2) !0._dp !Criver 
    ELSE
    C_old((n2+3):(n2+4)) = Criver 
    END IF


    !The new flux
    FL1=0._dp
    !The old Q
    Q_old= 0.5_dp*(QH+QH_old )  !This is like the discharge at the beginning of the time step 
    !The old Flux
    FL_old=Q_old*C_old
    !Useful for limiter, Hundsdorfer
    do i=2,n2+3
        if((FL_old(i+1).ne.FL_old(i)).and.((FL_old(i).ne.FL_old(i-1)))) THEN
            theta(i)=(FL_old(i)-FL_old(i-1))/(FL_old(i+1)-FL_old(i))  
        ELSE
            theta(i)=(FL_old(i)-FL_old(i-1)+eeps)/(FL_old(i+1)-FL_old(i) +eeps)  
        END IF
    end do

    !Calculate fluxes, with limiting -- 3rd order without limiting
    do i=2, n2+2
            !IF(sign(1._dp, Q_old(i-1)*Q_old(i+1))>0._dp) THEN
                    if(0.5_dp*(Q_old(i)+Q_old(i+1))>=0._dp) THEN
                    
                    limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/6._dp*theta(i) , mu*theta(i)))
                    if((i>n2).or.(i<4).or.(.false.)) limi(i)=0._dp !Force the boundaries to use an upwind flux
                    FL1(i)= FL_old(i)+limi(i)*(FL_old(i+1)-FL_old(i)) !(1._dp/6._dp)*( 2._dp*Q_old(i+1)*C_old(i+1) + 5._dp*Q_old(i)*C_old(i) -Q_old(i-1)*C_old(i-1)) 
                    ELSE
                    limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/(6._dp*theta(i+1)) , mu/theta(i+1)))
                    if((i>n2).or.(i<4).or.(.false.)) limi(i)=0._dp !Force the boundaries to use an upwind flux
                    FL1(i)=  FL_old(i+1) +limi(i)*(FL_old(i)-FL_old(i+1))  ! (1._dp/6._dp)*( 2._dp*Q_old(i)*C_old(i) + 5._dp*Q_old(i+1)*C_old(i+1) -Q_old(i+2)*C_old(i+2)) 
                    END IF
            !ELSE
            ! 
            !        FL1(i)=0._dp
            !END IF
    end do

    !Calculate C at the next half time step, predictor style
    do i=3, n2+2
    C(i) = (DT)/(0.5_dp*(A(i)+A(i)))*(Area_old(i)*C_old(i)/(DT) + Qe_old(i) -wset*C_old(i)*wetwidth(i) & !Source Terms
    - 1._dp/delX*( (FL1(i)-FL1(i-1)) - &  !Advection
    (0.5_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i))*(C_old(i+1)-C_old(i))/delX - & !Diffusion
    0.5_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1))*(C_old(i)-C_old(i-1))/delX )) ) !Diffusion
    end do

    C2(1:n2)=C(3:n2+2)

    !Mouth boundary conditions for the half step -- we impose zero gradient if the flow is
    !outward, and a mouth boundary condition otherwise
    !IF(Q_old(3)>0._dp) THEN
    !C(1)=Cmouth
    !C(2)=Cmouth
    !ELSE
    !C(1)=Cmouth !C(3)
    !C(2)=Cmouth !C(3)
    !END IF
    !
    !!River boundary conditions for the half step -- we impose zero gradient if the flow is
    !!outward, and a river boundary condition otherwise
    !IF(Q_old(n2+2)>0._dp) THEN
    !C(n2+3)=Criver !C(n2+2)
    !C(n2+4)=Criver !C(n2+2)
    !ELSE
    !C(n2+3)=Criver
    !C(n2+4)=Criver
    !END IF


    !write(23,*) Cpred
    !!Second half step --reuse lots of code



    !The second half Flux - evaluated at centred position - This is based on the
    !lecture notes from Hundsdorfer, page 
    !Fl_old=QH*Cpred
    !!Useful for limiter
    !do i=2,n2+3
    !if((FL_old(i+1).ne.FL_old(i)).and.((FL_old(i).ne.FL_old(i-1)))) THEN
    !theta(i)=(FL_old(i)-FL_old(i-1))/(FL_old(i+1)-FL_old(i))  
    !ELSE
    !theta(i)=(FL_old(i)-FL_old(i-1)+eeps)/(FL_old(i+1)-FL_old(i) +eeps)  
    !END IF
    !
    !if(theta(i).eq.0._dp) print*, 'zero theta'
    !if(isnan(theta(i))) print*, 'nan theta'
    !end do
    !
    !!Calculate fluxes, with limiting -- 3rd order without limiting
    !do i=3, n2+2
    !        if(0.5_dp*(QH(i)+QH(i+1))>=0._dp) THEN
    !        
    !        limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/6._dp*theta(i) , mu*theta(i)))
    !        if((i>n2).or.(i<4).or.(.false.)) limi(i)=0._dp !Force the boundaries to use an upwind flux
    !        FL1(i)= FL_old(i)+limi(i)*(FL_old(i+1)-FL_old(i)) 
    !        ELSE
    !        limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/(6._dp*theta(i+1)) , mu/theta(i+1)))
    !        if((i>n2).or.(i<4).or.(.false.)) limi(i)=0._dp !Force the boundaries to use an upwind flux
    !        FL1(i)= ( FL_old(i+1)+limi(i)*(FL_old(i)-FL_old(i+1)) )  
    !        END IF
    !end do
    !
    !!Calculate C at the next time step, corrector style
    !do i=3, n2+2
    !C(i) = DT/A(i)*(Area_old(i)*C_old(i)/DT + Qe(i) -wset*Cpred(i) & !Source Terms
    !- 1._dp/delX*( (FL1(i)-FL1(i-1)) - &  !Advection
    !(0.25_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i)+A(i+1)*D(i+1)+A(i)*D(i))*(Cpred(i+1)-Cpred(i))/delX - & !Diffusion
    !0.25_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1)+ A(i)*D(i)+A(i-1)*D(i-1))*(Cpred(i)-Cpred(i-1))/delX )) ) !Diffusion
    !end do
    !
    !!Update C2 - C2 holds the sediment concentration.
    !C2(1:n2)=C(3:n2+2) 

    !Hundsdorfer states that the use of eeps in theta (for the limiter) can allow
    !negative values of order eps. Let's get rid of them 
    !do i=1, n2
    !        if((i>1).and.(i<n2)) THEN
    !                if(C2(i)<0._dp) THEN
    !                ii=int(sign(1._dp, U2(i)))
    !                C2(i-ii)=C2(i-ii)-A(i)*C2(i)/A(i-ii)
    !                C2(i)=0._dp
    !                END IF
    !        else
    !                if(i==1) THEN
    !                        if(C2(i)<0._dp) THEN
    !                        ii=-1
    !                        C2(i-ii)=C2(i-ii)-A(i)*C2(i)/A(i-ii)
    !                        C2(i)=0._dp
    !                        end if 
    !                end if
    !                if(i==n2) THEN
    !                        if(C2(i)<0._dp) THEN
    !                        ii=1
    !                        C2(i-ii)=C2(i-ii)-A(i)*C2(i)/A(i-ii)
    !                        C2(i)=0._dp
    !                        end if 
    !                end if
    !        end if
    !END DO
    !
    !
    IF(minval(C2)<0._dp) THEN 
    !         print*, "min Sediment conc <0", minval(C2),  maxval(C_old2), minval(C_old2), maxval(U2), minval(U2)
    !         print*, "..........minloc C is", minloc(C2), "............ depth is =",A2(minloc(C2))/wetwidth2(minloc(C2))
    !!         print*, r
    !!         print*, "..........."
    !         print*, C_old(minloc(C2)+2), C_old(minloc(C2)-1+2), C_old(minloc(C2)+1+2), maxval(Qe2), minval(Qe2),&
    ! maxval(Qd2), minval(Qd2)
    !	 !print*, ">>>>>>>>>>>>>"
    !	 !print*, C
    !!         !print*, A !diag-(abs(up)+abs(lo))
    do i=1, n2
    IF(C2(i)<0._dp) THEN
            IF(abs(C2(i))<10._dp**(-10)) THEN !This could just be due to round off error in the matrix solver - fix it.
            C2(i)=0._dp
            ELSE
            print*, "C2 is < 0; violation is", abs(C2(i)), i, A2(i)
            stop
            END IF
    END IF
    END DO
    !!         
    !!        stop
    END IF

    !DO i=1, n
    !IF(C(i)<0._dp) C(i)=0._dp
    !END DO


    DO i = 1, n2
    Flag=isnan(C2(i))
    if(Flag) THEN
           print*, "sedconc is NAN"
           print*, U2
           print*, "................"
           print*, A2
           print*, Qd2
           print*, Qe2
           
            

            stop
    END IF
    END DO

end subroutine susconc_up34
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine susconc_up35(n2,DT, A2, QH2, QH2_old, delX,C2, U2, Qe2,Qe2_old, Qd2, Cmouth, C_old2, Area_old2, UU_old2, & 
                        Criver, wset2, wetwidth2,wetwidth_old2, D2,D2_old, third, pars_out)
    INTEGER, INTENT(IN):: n2
    REAL(dp), INTENT(IN):: DT, delX, Qe2, Qe2_old, Qd2, Cmouth, A2, QH2, QH2_old, U2, C_old2, Area_old2, UU_old2,Criver,&
                     wset2, wetwidth2,wetwidth_old2, D2, D2_old
    REAL(dp), INTENT(IN OUT):: C2, pars_out(10)
    LOGICAL, INTENT(IN):: third

    !!Note - here Qe is different to elsewhere. 
    DIMENSION C2(n2), A2(n2), QH2(n2), QH2_old(n2), U2(n2), Qe2(n2), Qe2_old(n2), Qd2(n2),C_old2(n2),Area_old2(n2),UU_old2(n2),&
              wetwidth2(n2),wetwidth_old2(n2), D2(n2), D2_old(n2),wset2(n2)

    INTEGER:: info, i, KL=2, KU=2, IPV(n2+4), ii !Note KL,KU are the number of upper and lower diagonals of the banded matrix
    REAL(dp):: r(n2) !Right hand side of equation
    LOGICAL:: flag
    !REAL(dp):: D(n) !dispersion constant-- actually not constant, need to fix.
    REAL(dp):: impcon=.50_dp, impconU=.00_dp
    !REAL(dp):: Qf1(n), Qf0(n), Af1(n), Af0(n), Df1(n), Df0(n)
    !REAL(dp):: band(7,n2+4), rhs(n2+4)
    REAL(dp):: A(n2+4), QH(n2+4), QH_old(n2+4), C(n2+4), U(n2+4), Qe(n2+4),Qe_old(n2+4), Qd(n2+4), C_old(n2+4), Area_old(n2+4), &
               UU_old(n2+4), wetwidth(n2+4),wetwidth_old(n2+4), D(n2+4), D_old(n2+4), Fl1(n2+4), Q_old(n2+4), limi(n2+4), & 
               FL_old(n2+4), theta(n2+4), Cpred(n2+4), diag(n2+4), lower(n2+4), upper(n2+4), rhs(n2+4),wset(n2+4)
    REAL(dp):: usef1(n2+4), usef2(n2+4), usef3(n2+4), usef4(n2+4)
    REAL(dp):: mu=1._dp, eeps=1.0E-10_dp
    !!##Solves (del AC/ del T) +  (del CQ / delX) =   d/dx (A Diffuse del C / delX) + (E-D)  

    !! See a test with an analytical advection-diffusion solution in:
    !!/home/gareth/Doc_Win/My_Documents/H_drive_Gareth/Maths and bits of Code/fortran code/Hydrodynamic model/full model/2009jan-2009date/good_version_without_N_support/port_with_namespace/most_updated_nov2309/third_order_sussed/with_bedload/lower_Q/bound/even_smallerQ/iforttry/vels1/simple_geo/1d_sedconcheck/analytical

    !!And with pure diffusion in:

    !!Note that E and D are the total rate of erosion/deposition over each section. 
    !Solved Using an approach developed from Hundsdorfer's notes, probably also in their book. 
    !Note that it is easy to have this explicit (stability requirement is much less
    !stringent than for 1D St Venant -- this has a stability determined by the velocity and delX only, no gravity waves)

    !!The numerical time stepping scheme is modified based n something described by Hundsdorfer's lecture notes
    !(page 49) as 'the implicit midpoint rule with Euler predictor', although it is
    !not implicit (but I use implicit diffusion on the first half step, which seems more stable).
    ! In the book Hundsdorfer and Verwer (2003) report the 'one
    !step explicit midpoint rule' (page 142)
    !Which is the basis of this method, the only difference being that I treat
    !diffusion implicitly.
    !Basically we take a predictor half-step
    !Cpred = Clast + (1/2 delt)/delx*( ADVECTION(tlast,Clast) + DIFFUSION(CPRED,tlast+1/2delT) )
    !And then a corrector full step
    !C=Clast+ delt/delX*(ADVECTION(t+1/2 delT, Cpred) + DIFFUSION(C,tlast+delT) )
    !On the first step, the discharge is evaluated as QH2, which is
    !the conservative discharge estimates from the McCormack method. Note
    !that this is a reasonable estimate of the discharge at the old time level if we
    !assume that the discharge is constant between the last step and the next one -
    !While this is a crude assumption, more complex things I tried (e.g Q =
    !0.5*(QH2+QH2_old) lead to greater overshoot near the mouth than the present
    !approach. (NOT UP TO DATE - PRESENTLY I AM USING THE LATTER - HOPEFULLY FINE).
    !Note that while we could just use the straight discharge output from McCormack, it would not be
    !conservative. I have checked this by running a constant discharge case with no
    !erosion or deposition - the present algorithm predicts a constant sediment
    !concentration (and the QH2 discharge is constant), while if we use the straight
    !output from McCormack, then both the Discharge and sediment conc show slight
    !variation from constant.
    ! The diffusion coefs are evaluated implicitly (more stable and accurate than
    ! explicit), as is the rate of erosion. 
    !On the second step, the discharge is evaluated as QH (a good half-time step
    !estimate), as are the rate of erosion and the diffusion coef, by suitable
    !averaging of the input variables 

    !The spatial discretization is a standard central scheme for diffusion, and a
    !flux limited third order upwind scheme for advection. The latter is described
    !in Hundsdorfer's notes (page 38)

    !Note -- the discharge issues. With the McCormack Scheme, we can show that
    !A(t+1)-A(t) = (dT/dX)*(0.5*(Q(t)_[i+1]+Q(Pred)_[i])-0.5*(Q(t)_[i]+Q(Pred)_[i-1]) )
    !This means that at steady state (constant discharge), it is actually (Q(t)_[i+1]+Q(Pred)_[i]) which
    !will not be changing either in time or in space.

    !Although this is a common reference, I
    !found Vreugdenhill (1989:59) a useful reference, and the way the code USED TO be
    !written reflects that. 

    !Flux form of third order advection. dF/dx = [ F(i+1/2)-F(i-1/2) ] / dx
    !with F[i+1/2] = (1/6)[ - F(i-1) +5F(i) +2F(i+1) ]     if (Q(i+1/2) > 0)
    !F[i+1/2] = (1/6)[ 2F(i) +5F(i+1) -F(i+2) ] if ( Q(i+1/2) <0 )

    !For the third order advection, I have used 'ghost' cells. There are 2
    !downstream, and 2 upstream
    !!At the mouth boundary, a given value is imposed if the flow is inward, and
    !otherwise a zero gradient condition is imposed.
    !! At the landward boundary, a zero gradient condition is enforced if there is
    !outflow. Otherwise, the river concentration is enforced, or used to
    !provide other boundary values, depending on the version of this code you have.  

    !Note that if we force advection to first order, presently I can find weird
    !long-term behaviour - accumulation of sediment in upstream zones -odd.

    !!Define new variables with 'ghost'points on the edges, useful for implementing
    !boundary conditions
    A(3:(n2+2))=A2
    QH(3:(n2+2))=QH2
    QH_old(3:(n2+2))=QH2_old
    C(3:(n2+2))=C2
    U(3:(n2+2))=U2
    Qe(3:(n2+2))=Qe2
    Qe_old(3:(n2+2))=Qe2_old
    Qd(3:(n2+2))=Qd2
    C_old(3:(n2+2))=C_old2
    Area_old(3:(n2+2))=Area_old2
    UU_old(3:(n2+2))=UU_old2
    wetwidth(3:(n2+2))=wetwidth2
    wetwidth_old(3:(n2+2))=wetwidth_old2
    D(3:(n2+2))=D2
    D_old(3:(n2+2))=D2_old
    wset(3:n2+2)=wset2

    !Boundary conditions
    A(1:2)=A2(1)
    A((n2+3):(n2+4)) = A2(n2)
    QH(1:2)=QH2(1)
    QH((n2+3):(n2+4)) = QH2(n2)
    QH_old(1:2)=QH2_old(1)
    QH_old((n2+3):(n2+4)) = QH2_old(n2)
    Area_old(1:2)=Area_old2(1)
    Area_old((n2+3):(n2+4)) = Area_old2(n2)
    U(1:2)=U2(1)
    U((n2+3):(n2+4)) = U2(n2)
    UU_old(1:2)=U2(1)
    UU_old((n2+3):(n2+4)) = UU_old2(n2)
    Qe(1:2)=Qe(3)!0._dp
    Qe((n2+3):(n2+4)) = Qe(n2+2) !0._dp 
    Qe_old(1:2)= Qe_old(3) !0._dp
    Qe_old((n2+3):(n2+4))= Qe_old(n2+2) ! 0._dp 
    Qd(1:2)=0._dp
    Qd((n2+3):(n2+4)) = 0._dp 
    wetwidth(1:2)=wetwidth2(1)
    wetwidth((n2+3):(n2+4)) = wetwidth2(n2)
    wetwidth_old(1:2)=wetwidth_old2(1)
    wetwidth_old((n2+3):(n2+4)) = wetwidth_old2(n2)
    D(1:2)=D2(1)
    D((n2+3):(n2+4))=D2(n2)
    D_old(1:2)=D2_old(1)
    D_old((n2+3):(n2+4))=D2_old(n2)
    wset(1:2)=wset2(1)
    wset(n2+3:n2+4)=wset2(n2)
    !!If I cut diffusion at the mouth, mass conservation seems much better - weird, might be a bug deeper in my code?
    !D(1:3)=0._dp
    !D_old(1:3)=0._dp


    IF(QH2(1)>0._dp) THEN
        C_old(1)=Cmouth
        C_old(2)=Cmouth !0.5_dp*(Cmouth+C_old(3))
    ELSE
        C_old(2)=C_old(3) !max(2._dp*C_old(3)-C_old(4),0._dp)
        C_old(1)=C_old(3) !max(2._dp*C_old(2) - C_old(3),0._dp )
    END IF

    IF(QH2(n2)>0._dp) THEN
        C_old((n2+3):(n2+4)) = C_old2(n2) !0._dp !Criver 
    ELSE
        C_old((n2+3):(n2+4)) = Criver 
    END IF


    !The new flux
    FL1=0._dp
    !The old Q
    Q_old= 0.5_dp*(QH+QH_old )!This is like the discharge at the beginning of the time step 
    !The old Flux
    FL_old=Q_old*C_old
    !Useful for limiter, Hundsdorfer
    DO i=2,n2+3
        !if(abs(FL_old(i+1)-FL_old(i))>1.0E-12_dp) THEN
        theta(i)=(FL_old(i)-FL_old(i-1))/(FL_old(i+1)-FL_old(i))  
        !ELSE
        !        if(abs(FL_old(i)-FL_old(i-1))>1.0E-12_dp) THEN
        !!        theta(i)=1.0E+9_dp !(FL_old(i)-FL_old(i-1))/(FL_old(i+1)-FL_old(i) +eeps)  
        !        else
        !        theta(i)= 1._dp !1._dp !(FL_old(i)-FL_old(i-1)+eeps)/(FL_old(i+1)-FL_old(i) +eeps)  
        !        end if
        !END IF
        !Catch special cases - theta = nan or 0
        IF(isnan(theta(i))) THEN
            IF(FL_old(i).NE.FL_old(i-1)) THEN
                theta(i)=1.0E+12_dp*sign(1._dp, FL_old(i)-FL_old(i-1))
            ELSE
                theta(i)=1.0_dp
            END IF
        END IF

        IF(theta(i).eq.0._dp) THEN
            IF(FL_old(i+1).ne.(FL_old(i))) THEN
                theta(i)=1.0E-12_dp*sign(1._dp, FL_old(i+1)-FL_old(i))
            ELSE
                theta(i)=1.0_dp 
            END IF
        END IF
    END DO

    !Calculate fluxes, with limiting -- 3rd order without limiting
    DO i=2, n2+2
        IF(0.5_dp*(Q_old(i)*C_old(i)+Q_old(i+1)*C_old(i+1))>=0._dp) THEN
            limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/6._dp*theta(i) , mu*theta(i)))
            !limi(i)= 0.5_dp*(theta(i)+abs(theta(i)))/(1._dp+abs(theta(i)))
            !if((i>n2).or.(i<4).or.(.false.)) limi(i)=0._dp !Force the boundaries to use an upwind flux
            FL1(i)= FL_old(i)+limi(i)*(FL_old(i+1)-FL_old(i)) !(1._dp/6._dp)*( 2._dp*Q_old(i+1)*C_old(i+1) + 5._dp*Q_old(i)*C_old(i) -Q_old(i-1)*C_old(i-1)) 
        ELSE
            limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/(6._dp*theta(i+1)) , mu/theta(i+1)))
            !limi(i)= 0.5_dp*(1._dp/theta(i+1)+abs(1._dp/theta(i+1)))/(1._dp+abs(1._dp/theta(i+1)))
            !if((i>n2).or.(i<4).or.(.false.)) limi(i)=0._dp !Force the boundaries to use an upwind flux
            FL1(i)= ( FL_old(i+1)+limi(i)*(FL_old(i)-FL_old(i+1)) ) ! (1._dp/6._dp)*( 2._dp*Q_old(i)*C_old(i) + 5._dp*Q_old(i+1)*C_old(i+1) -Q_old(i+2)*C_old(i+2)) 
        END IF
    END DO

    upper=0._dp
    lower=0._dp
    diag=0._dp
    rhs=0._dp
    !Calculate C at the next half time step, predictor style
    DO i=3, n2+2
        !C(i) = (0.5_dp*DT)/(0.5_dp*(A(i)+A(i)))*(Area_old(i)*C_old(i)/(DT*0.5_dp) + Qe_old(i) -wset*C_old(i)*wetwidth_old(i) & !Source Terms
        !C(i) = ((0.5_dp*(A(i)+Area_old(i)))/(0.5_dp*DT)+wset*0.5_dp*(wetwidth_old(i)+wetwidth(i)))**(-1._dp)*&  !Constant of C(i), including implicit settling
        !(Area_old(i)*C_old(i)/(DT*0.5_dp) &  !Unsteady
        !+ Qe(i) & !Erosion at half time (Source Terms)
        !- 1._dp/delX*( (FL1(i)-FL1(i-1)) - &  !Advection
        !(0.5_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i))*(C_old(i+1)-C_old(i))/delX - & !Diffusion
        !0.5_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1))*(C_old(i)-C_old(i-1))/delX )) ) !Diffusion

        upper(i)=-0.25_dp*(Area_old(i+1)*D(i+1)+Area_old(i)*D(i)+A(i+1)*D(i+1)+A(i)*D(i))*1._dp/delX**2 !Diffusion
        lower(i)=-0.25_dp*(Area_old(i)*D(i)+Area_old(i-1)*D(i-1)+A(i)*D(i)+A(i-1)*D(i-1))*1._dp/delX**2 !Diffusion
        diag(i)=(0.5_dp*(A(i)+Area_old(i)))/(0.5_dp*DT)+& !Time derivative
                wset(i)*0.5_dp*(wetwidth(i)+wetwidth(i)) & !Settling
                -upper(i) - lower(i) !Diffusion
        rhs(i)= Area_old(i)*C_old(i)/(DT*0.5_dp) & !Unsteady 
                + Qe(i)*8._dp/8._dp +(Qe(i+1)+Qe(i-1))*0._dp/8._dp & !Erosion
                -1._dp/delX*(FL1(i)-FL1(i-1)) !Advection

    END DO




    !Mouth boundary conditions for the half step -- we impose zero gradient if the flow is
    !outward, and a mouth boundary condition otherwise
    IF(Q_old(3)>0._dp) THEN
        !C(1)=Cmouth
        diag(1)=1._dp
        rhs(1)=Cmouth

        !C(2)=Cmouth !0.5_dp*(Cmouth+C(3))
        diag(2)=1._dp
        rhs(2)=Cmouth
    ELSE
        !C(2)=C(3) !max(2._dp*C(3) -C(4), 0._dp)
        diag(2)=1._dp
        upper(2)=-1._dp

        !C(1)=C(2) !max(2._dp*C(2)-C(3), 0._dp)
        diag(1)=1._dp
        upper(1)=-1._dp
    END IF

    !River boundary conditions for the half step -- we impose zero gradient if the flow is
    !outward, and a river boundary condition otherwise
    IF(Q_old(n2+2)>0._dp) THEN
        !C(n2+3)=C(n2+2)
        diag(n2+3)=1._dp
        lower(n2+3)=-1._dp

        !C(n2+4)=C(n2+3)
        diag(n2+4)=1._dp
        lower(n2+4)=-1._dp
    ELSE
        !C(n2+3)=Criver
        diag(n2+3)=1._dp
        rhs(n2+3)=Criver

        !C(n2+4)=Criver
        diag(n2+4)=1._dp
        rhs(n2+4)=Criver
    END IF

    usef1=lower
    usef2=diag
    usef3=upper
    usef4=rhs
    !Solve it
    call dgtsv(n2+4,1, lower(2:n2+4), diag(1:n2+4), upper(1:n2+3), rhs(1:n2+4),n2+4, info)
    Cpred=rhs
    IF(info.ne.0) THEN
        PRINT*, 'info .ne. 0 in susconc_up35 , pred step', info
    END IF

    DO i = 1, n2+4
        Flag=isnan(Cpred(i))
        IF(Flag) THEN
            PRINT*, "Cpred is NAN"
            PRINT*, Q_old
            PRINT*, "................"
            PRINT*, A2
            PRINT*, Qd2
            PRINT*, Qe2
            PRINT*,"........MATRIX DIAGONALS.."
            PRINT*, '.......Lower.......'
            PRINT*, usef1(2:n2+4)
            PRINT*, '.......Main.......'
            PRINT*, usef2(1:n2+4)
            PRINT*, '.......Upper.......'
            PRINT*, usef3(1:n2+3)
            PRINT*, '.......RHS.......'
            PRINT*, usef4(1:n2+4)
            PRINT*, '.......FL1.......'
            PRINT*, FL1
            PRINT*, '.......theta.......'
            PRINT*, theta
            EXIT
        END IF
    END DO
    !write(23,*) Cpred
    !!Second half step --reuse lots of code

    !The second half Flux - evaluated at centred position - This is based on the
    !lecture notes from Hundsdorfer, page 
    Fl_old=QH*Cpred
    !Useful for limiter
    DO i=2,n2+3
        theta(i)=(FL_old(i)-FL_old(i-1))/(FL_old(i+1)-FL_old(i))  
        
        !Catch special cases - theta = nan or 0
        IF(isnan(theta(i))) THEN
            IF(FL_old(i).NE.FL_old(i-1)) THEN
                theta(i)=1.0E+12_dp*sign(1._dp, FL_old(i)-FL_old(i-1))
            ELSE
                theta(i)=1.0_dp
            END IF
        END IF

        IF(theta(i).eq.0._dp) THEN
            IF(FL_old(i+1).ne.(FL_old(i))) THEN
                theta(i)=1.0E-12_dp*sign(1._dp, FL_old(i+1)-FL_old(i))
            ELSE
                theta(i)=1.0_dp 
            END IF
        END IF
    END DO

    !Calculate fluxes, with limiting -- 3rd order without limiting
    DO i=2, n2+2
        IF(0.5_dp*(QH(i)*Cpred(i)+QH(i+1)*Cpred(i+1))>=0._dp) THEN
            limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/6._dp*theta(i) , mu*theta(i)))
            !limi(i)= 0.5_dp*(theta(i)+abs(theta(i)))/(1._dp+abs(theta(i)))
            !if((i>n2).or.(i<4).or.(.false.)) limi(i)=0._dp !Force the boundaries to use an upwind flux
            FL1(i)= FL_old(i)+limi(i)*(FL_old(i+1)-FL_old(i)) 
        ELSE
            limi(i)=max(0._dp, min(1._dp, 1._dp/3._dp+1._dp/(6._dp*theta(i+1)) , mu/theta(i+1)))
            !limi(i)= 0.5_dp*(1._dp/theta(i+1)+abs(1._dp/theta(i+1)))/(1._dp+abs(1._dp/theta(i+1)))
            !if((i>n2).or.(i<4).or.(.false.)) limi(i)=0._dp !Force the boundaries to use an upwind flux
            FL1(i)= ( FL_old(i+1)+limi(i)*(FL_old(i)-FL_old(i+1)) )  
        END IF
    END DO

    !Calculate C at the next time step, corrector style

    upper=0._dp
    diag=0._dp
    lower=0._dp
    rhs=0._dp

    DO i=3, n2+2
        !C(i) = DT/A(i)*(Area_old(i)*C_old(i)/DT + Qe(i) -wset*Cpred(i)*0.5_dp*(wetwidth(i) +wetwidth_old(i)) & !Source Terms
        !- 1._dp/delX*( (FL1(i)-FL1(i-1)) - &  !Advection
        !(0.25_dp*(Area_old(i+1)*D_old(i+1)+Area_old(i)*D_old(i)+A(i+1)*D(i+1)+A(i)*D(i))*(Cpred(i+1)-Cpred(i))/delX - & !Diffusion
        !0.25_dp*(Area_old(i)*D_old(i)+Area_old(i-1)*D_old(i-1)+ A(i)*D(i)+A(i-1)*D(i-1))*(Cpred(i)-Cpred(i-1))/delX )) ) !Diffusion
        upper(i)= -0.5_dp*0.5_dp*(0.5_dp*(A(i+1)+Area_old(i+1))*D(i+1) +0.5_dp*(A(i)+Area_old(i))*D(i))/delX**2 !Diffusion
        lower(i)= -0.5_dp*0.5_dp*(0.5_dp*(A(i)+Area_old(i))*D(i) +0.5_dp*(A(i-1)+Area_old(i-1))*D(i-1))/delX**2 !Diffusion
        diag(i)= A(i)/DT & !Unsteady
                -upper(i) - lower(i) & !Diffusion
                +0.5_dp*min(wset(i)*0.5_dp*(wetwidth(i) +wetwidth(i)), 0.5_dp*(Area_old(i)+A(i))/DT) !Deposition
        rhs(i) = Area_old(i)*C_old(i)/DT & !Unsteady
                !+ Qe(i) -min(wset(i)*0.5_dp*(wetwidth(i) +wetwidth(i)), 0.5_dp*(Area_old(i)+A(i))/DT)*Cpred(i) & !Erosion and deposition
                + Qe(i)*6._dp/8._dp +1._dp/8._dp*(Qe(i+1)+Qe(i-1)) & !Erosion 
                -0.5_dp*min(wset(i)*0.5_dp*(wetwidth(i) +wetwidth(i)), 0.5_dp*(Area_old(i)+A(i))/DT)*C_old(i) & ! deposition
                - 1._dp/delX*(FL1(i)-FL1(i-1)) & !Advection
                -upper(i)*C_old(i+1) +upper(i)*C_old(i) +lower(i)*C_old(i) - lower(i)*C_old(i-1) ! Diffusion
    END DO

    IF(QH(3)>0._dp) THEN
        !C(1)=Cmouth
        diag(1)=1._dp
        rhs(1)=Cmouth
        !C(2)=Cmouth !0.5_dp*(Cmouth+C(3))
        diag(2)=1._dp
        rhs(2)=Cmouth
    ELSE
        !C(2)=C(3) !max(2._dp*C(3) -C(4), 0._dp)
        diag(2)=1._dp
        upper(2)=-1._dp
        !C(1)=C(2) !max(2._dp*C(2)-C(3), 0._dp)
        diag(1)=1._dp
        upper(1)=-1._dp
    END IF

    !River boundary conditions for the half step -- we impose zero gradient if the flow is
    !outward, and a river boundary condition otherwise
    IF(QH(n2+2)>0._dp) THEN
        !C(n2+3)=C(n2+2)
        diag(n2+3)=1._dp
        lower(n2+3)=-1._dp

        !C(n2+4)=C(n2+3)
        diag(n2+4)=1._dp
        lower(n2+4)=-1._dp
    ELSE
        !C(n2+3)=Criver
        diag(n2+3)=1._dp
        rhs(n2+3)=Criver

        !C(n2+4)=Criver
        diag(n2+4)=1._dp
        rhs(n2+4)=Criver
    END IF

    usef1=lower
    usef2=diag
    usef3=upper
    usef4=rhs

    call dgtsv(n2+4,1, lower(2:n2+4), diag(1:n2+4), upper(1:n2+3), rhs(1:n2+4),n2+4, info)
    !Update C2 - C2 holds the sediment concentration.
    C2(1:n2)=rhs(3:n2+2) 

    !Record the boundary fluxes - useful for checking conservation.
    pars_out(1)=FL1(2) - &  !Advection
            0.25_dp*((A(3)+Area_old(3))*D(3) +(A(3-1)+Area_old(3-1))*D(3-1))/delX*0.5_dp* & !Diffusion coeff
            ((C_old(3)-C_old(2))+(rhs(3)-rhs(2))) !Sus gradient
    pars_out(2)=FL1(n2+2) - &  !Advection
            0.25_dp*((A(n2+3)+Area_old(n2+3))*D(n2+3) +(A(n2+3-1)+Area_old(n2+3-1))*D(n2+3-1))/delX*0.5_dp* & !Diffusion coef
            ((C_old(n2+3)-C_old(n2+2))+(rhs(n2+3)-rhs(n2+2))) !Sus gradient


    IF(info.ne.0) THEN
        PRINT*, 'info .ne. 0 in susconc_up35 , cor step', info
    END IF

    !Hundsdorfer states that the use of eeps in theta (for the limiter) can allow
    !negative values of order eps. Let's get rid of them 
    DO i=1, n2
        IF((i>1).AND.(i<n2)) THEN
            IF(C2(i)<0._dp) THEN
                ii=int(sign(1._dp, U2(i)))
                C2(i-ii)=C2(i-ii)-A(i+2)*C2(i)/A(i-ii+2)
                C2(i)=0._dp
            END IF
        ELSE
            IF(i==1) THEN
                IF(C2(i)<0._dp) THEN
                    ii=-1
                    C2(i-ii)=C2(i-ii)-A(i+2)*C2(i)/A(i-ii+2)
                    C2(i)=0._dp
                END IF 
            END IF
            IF(i==n2) THEN
                IF(C2(i)<0._dp) THEN
                    ii=1
                    C2(i-ii)=C2(i-ii)-A(i+2)*C2(i)/A(i-ii+2)
                    C2(i)=0._dp
                END IF 
            END IF
        END IF
    END DO

    IF(minval(C2)<0._dp) THEN 
        PRINT*, "min Sediment conc <0", minval(C2),  maxval(C_old2), minval(C_old2), maxval(U2), minval(U2)
        PRINT*, "..........minloc C is", minloc(C2), "............ depth is =",A2(minloc(C2))/wetwidth2(minloc(C2))
        PRINT*, C_old(minloc(C2)+2), C_old(minloc(C2)-1+2), C_old(minloc(C2)+1+2), maxval(Qe2), minval(Qe2),&
                maxval(Qd2), minval(Qd2)
        !print*, ">>>>>>>>>>>>>"
        !print*, C
        !print*, A !diag-(abs(up)+abs(lo))
        DO i=1, n2
            IF(C2(i)<0._dp) THEN
                IF(abs(C2(i))<10._dp**(-10)) THEN !This could just be due to round off error in the matrix solver - fix it.
                    C2(i)=0._dp
                ELSE
                    print*, "C2 is < 0; violation is", abs(C2(i))
                    !stop
                END IF
            END IF
        END DO
        !         
        !        stop
    END IF

    !DO i=1, n
    !IF(C(i)<0._dp) C(i)=0._dp
    !END DO

    DO i = 1, n2
    Flag=isnan(C2(i))
        IF(Flag) THEN
            PRINT*, "sedconc is NAN"
            PRINT*, U2
            PRINT*, "................"
            PRINT*, A2
            PRINT*, Qd2
            PRINT*, Qe2
           
            PRINT*,"........MATRIX DIAGONALS.."
            PRINT*, '.......Lower.......'
            PRINT*, usef1(2:n2+4)
            PRINT*, '.......Main.......'
            PRINT*, usef2(1:n2+4)
            PRINT*, '.......Upper.......'
            PRINT*, usef3(1:n2+3)
            PRINT*, '.......RHS.......'
            PRINT*, usef4(1:n2+4)
            !EXIT
            
            STOP
        END IF
    END DO

end subroutine susconc_up35
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine susconc_constD(n,DT, A,delX,C, U, Qe, Qd, Cmouth, C_old, Area_old, UU_old, Criver, wset, wetwidth)
    INTEGER, intent(in):: n
    REAL(dp), intent(in):: DT, delX, Qe, Qd, Cmouth, A, U, C_old, Area_old, UU_old,Criver, wset, wetwidth
    REAL(dp), intent(in out):: C


    !!Note - here Qe is different to elsewhere. 
    Dimension C(n), A(n), U(n), Qe(n), Qd(n), C_old(n), Area_old(n), UU_old(n), wetwidth(n)

    INTEGER:: info, i
    REAL(dp):: diag(n), up(n), lo(n)  !Matrix diagonals
    REAL(dp):: r(n) !Right hand side of equation
    logical:: flag
    REAL(dp):: D(n) !dispersion constant-- actually not constant, need to fix.
    REAL(dp):: impcon=.5_dp
    REAL(dp):: Qf1(n), Qf0(n), Af1(n), Af0(n), Df1(n), Df0(n)

    D= 2._dp!*abs(U)
    !C(n)=0.
    !!!
    !!##Solves (del AC/ del T) +  (del CQ / delX) =   d/dx (A D del C / delX) + (E-D)  

    !!Note that E and D are the total rate of erosion/deposition over each section. 
    !!Solved with a Crank-Nicholson Method. Although this is a common reference, I
    !found Vreugdenhill (1989:59) a useful reference, and the way the code is
    !written reflects that.  
    !!At the mouth boundary, a given value is imposed
    !! At the landward boundary, a zero gradient condition is enforced if there is
    !no river inflow. Otherwise, the river concentration is enforced, or used to
    !provide other boundary values, depending on the version of this code you have.  

    !!!START
    !Calculate the diagonals for the matrix inversion

    !!Forward average of some important variables
    Qf1(1:n-1)=.5_dp*(A(2:n)*U(2:n)+A(1:n-1)*U(1:n-1))
    Af1(1:n-1)= .5_dp* (A(2:n)+A(1:n-1))
    Df1(1:n-1)= .5_dp*(D(1:n-1)+D(2:n))

    !!Forward average of some important old variables
    Qf0(1:n-1)=.5_dp*(Area_old(2:n)*UU_old(2:n)+Area_old(1:n-1)*UU_old(1:n-1))
    Af0(1:n-1)= .5_dp* (Area_old(2:n)+Area_old(1:n-1))
    Df0(1:n-1)= .5_dp*(D(1:n-1)+D(2:n))


    !!Upper diagonal
    up(1:n-1)=(impcon/delX)*( Qf1(1:n-1)*.5_dp-Af1(1:n-1)*Df1(1:n-1)/delX)
    up(n)=0._dp
    !!Lower diagonal
    lo(2:n)=(impcon/delX)*(-Qf1(1:n-1)*.5_dp-Af1(1:n-1)*Df1(1:n-1)/delX)!-up(1:n-1) !-(impcon/delX)*& ((.5_dp*(A(2:n)*U(2:n)+A(1:n-1)*U(1:n-1)))*.5_dp+.5_dp*(A(2:n)+A(1:n-1))*.5_dp*(D(1:n-1)+D(2:n))/delX)
    lo(1)=0._dp

    !Main diagonal
    diag(2:n-1)= A(2:n-1)/DT + wset*wetwidth(2:n-1)+ &
    (impcon/delX)*(Qf1(2:n-1)*.5_dp - Qf1(1:n-2)*.5_dp  + Af1(2:n-1)*Df1(2:n-1)/delX + Af1(1:n-2)*Df1(1:n-2)/delX) !up(2:n-1) +(impcon/delX)*.5_dp*(A(2:nn-1)+A(3:nn))*(D(2:n-1)+D(3:nn))/delX + lo(2:n-1) &
    !+ .5_dp*(A(2:nn-1)+A(1:nn-2))*(D(1:n-2)+D(2:n-1))/delX

    !!Right hand side
    r(2:n-1) = Qe(2:n-1) +Area_old(2:n-1)*C_old(2:n-1)/DT - (1._dp-impcon)/delX*( & 
    ( Qf0(2:n-1)*.5_dp*(C_old(3:n)+C_old(2:n-1)) & !Advective
    - Af0(2:n-1)*Df0(2:n-1)*(C_old(3:n)-C_old(2:n-1))/delX ) & !diffusive
    - ( Qf0(1:n-2)*.5_dp*(C_old(2:n-1)+C_old(1:n-2)) & !Advective
    - Af0(1:n-2)*Df0(1:n-2)*(C_old(2:n-1)-C_old(1:n-2))/delX ) ) !diffusive

    !!!Boundary conditions

    !!CASE 1 - At the mouth, ebbing
    ! Assume the discharge and area at the mouth are the same as at 1, and the sediment
    !concentration beyond the boundary is the same at at 1. 
    If(U(1)<=0._dp) THEN
    diag(1)= A(1)/DT+ wset*wetwidth(1)+ (impcon/delX)*(Qf1(1)*.5_dp - Qf1(1)*.5_dp*2._dp  + Af1(1)*Df1(1)/delX + 0._dp ) !So the  Qf(1)*.5_dp*2._dp reflects the addition of the value at the mouth, and the final 0._dp reflects the Af(0)*Df(1)/delX*(C(1)-C(mouth))/delX which is zero in this situation

    r(1) = Qe(1) + Area_old(1)*C_old(1)/DT - (1._dp-impcon)/delX*( &
    (Qf0(1)*.5_dp*(C_old(1)+C_old(2)) - Af0(1)*Df0(1)*(C_old(2)-C_old(1))/delX) &
    -( Qf0(1)*.5_dp*(C_old(1)+C_old(1)) - Af0(1)*Df0(1)*(0._dp) )) !!Note that here, C_old(1) takes the place of C_mouth
    END IF

    !!! CASE 2- At the mouth, flooding
    If(U(1)>=0._dp) THEN
    diag(1)= A(1)/DT+ wset*wetwidth(1)+ (impcon/delX)*(Qf1(1)*.5_dp - Qf1(1)*.5_dp + Af1(1)*Df1(1)/delX + Af1(1)*Df1(1)/delX ) !Here the c_mouth terms go on the RHS

    !Note that the beginning of the RHS includes the Cmouth terms, which are at time
    !1 rather than 0
    r(1) = impcon/delX*(Qf1(1)*.5_dp*Cmouth +Af1(1)*Df1(1)/delX*Cmouth) + Qe(1) + Area_old(1)*C_old(1)/DT - & 
    (1._dp-impcon)/delX*( (Qf0(1)*.5_dp*(C_old(1)+C_old(2)) - Af0(1)*Df0(1)*(C_old(2)-C_old(1))/delX) &
    -( Qf0(1)*.5_dp*(C_old(1)+ Cmouth) - Af0(1)*Df0(1)*(C_old(1)-Cmouth)/delX))
    END IF



    !!!!Next case - the behaviour at the landward end - I suggest a zero gradient
    !condition is reasonable if there is no river input

    if(Criver==0._dp) THEN
    diag(n)= A(n)/DT+ wset*wetwidth(n)+ impcon/delX*(A(n)*U(n)*.5_dp*2._dp - Qf1(n-1)*.5_dp +0._dp + Af1(n-1)*Df1(n-1)/delX) !Here the A(n)*U(n)*.5_dp*2._dp accounts for the zero gradient in the upper advective part, and the  +0._dp also accounts for the zero gradient in the upper diffusive part
    r(n)= Qe(n)+ Area_old(n)*C_old(n)/DT - (1._dp-impcon)/delX*( &
    (Area_old(n)*UU_old(n)*C_old(n) - 0._dp)- (Qf0(n-1)*.5_dp*(C_old(n) +C_old(n-1)) &
    -Af0(n-1)*Df0(n-1)*(C_old(n)-C_old(n-1))/delX))  
    ELSE
    diag(n)= A(n)/DT+wset*wetwidth(n)+ impcon/delX*( A(n)*U(n)*.5_dp - Qf1(n-1)*.5_dp +A(n)*D(n)/delX + Af1(n-1)*Df1(n-1)/delX) 

    r(n)= Qe(n)+ impcon/delX*(-A(n)*U(n)*.5_dp*Criver + A(n)*D(n)/delX*Criver) + Area_old(n)*C_old(n)/DT & 
    - (1._dp-impcon)/delX*((Area_old(n)*UU_old(n)*.5_dp*(C_old(n)+Criver) -Area_old(n)*D(n)*(Criver-C_old(n))/delX )- &
    (Qf0(n-1)*.5_dp*(C_old(n) +C_old(n-1)) -Af0(n-1)*Df0(n-1)*(C_old(n)-C_old(n-1))/delX))  

    !diag(n)=1._dp
    !r(n)=Criver
    !lo(n)=0._dp
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!Solve

    !print*, maxval(C_old), minval(C_old)
    !print*, lo, diag, up, r
    !stop
    call DGTSV(n,1, lo(2:n), diag, up(1:n-1), r,n, info)

    !call DGBSV(n,2,2)

    !The solver writes the new C to r- so we change it back.
    C=r

    !print*, info
    !print*, maxval(diag), minval(diag), maxval(lo), minval(lo), maxval(up), minval(up)
    !print*, "vel", maxval(U), minval(U)

    IF(minval(C)<0._dp) THEN 
             print*, "min Sediment conc <0", minval(C),  maxval(C_old), minval(C_old), maxval(U), minval(U)
             print*, "..........minloc C is", minloc(C), "............"
    !         print*, r
    !         print*, "..........."
             print*, maxval(Qe), minval(Qe), maxval(Qd), minval(Qd)
         !print*, ">>>>>>>>>>>>>"
         !print*, C
    !         !print*, A !diag-(abs(up)+abs(lo))
    do i=1, n
    IF(C(i)<0._dp) THEN

        IF(abs(C(i))<10._dp**(-10)) THEN !This could just be due to round off error in the matrix solver - fix it.
        C(i)=0._dp
        ELSE
        print*, "C is < 0; violation is", abs(C(i))
    !	stop
        END IF
    END IF	
    END DO
    !         
    !        stop
    END IF

    !DO i=1, n
    !IF(C(i)<0._dp) C(i)=0._dp
    !END DO


    DO i = 1, n
    Flag=isnan(C(i))
    if(Flag) THEN
           print*, "sedconc is NAN"
           print*, U
           print*, "................"
           print*, A
           print*, Qd
           print*, Qe
           
            

            stop
    END IF
    END DO

end subroutine susconc_constD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!Attempt at a 2d susconc solver, using a simplified version of the
!equations - an unsteady term, longitudinal advection, lateral diffusion, and
!erosion and deposition.
!!Numerical method is Peaceman-Rachford ADI, see Colub and Ortega (1981:281)
subroutine susconc2d(a,b, DT,delX, lengths,elevs,elevs_old,vels,vels_old,hs, Cdist,Cdist_old, Qe, Qe_old, Cmouth, Criver,&
wset,fs,slopes)
    !use matrix_solvers
    !use crosssection

    !implicit none
    INTEGER, intent(in):: a, b
    REAL(dp), intent(in):: DT,delX, wset
    REAL(dp), intent(in):: elevs, elevs_old, vels, vels_old,hs, Cdist_old, Qe,Qe_old, Cmouth, Criver,fs, lengths, slopes
    REAL(dp), intent(inout):: Cdist
    DIMENSION elevs(b), elevs_old(b), vels(a,b), vels_old(a,b), hs(a,b), Cdist(a,b), Cdist_old(a,b), Qe(a,b),Qe_old(a,b), & 
    fs(a,b), lengths(a,b),slopes(a,b), Criver(a), Cmouth(a)

    REAL(dp):: depthN(a,b), depthH(a,b), depthF(a,b), elevsH, velh(a,b),sus_h(a,0:b+1), eddN(a,b), veldH(a,0:b+1)
    REAL(dp):: sus_N(0:a+1,b), lnths(0:a+1,b), eddN_yh(0:a,b) , eddH(a,b), eddH_xh(a,0:b), eddF(a,b), eddF_yh(0:a,b)
    REAL(dp):: upperb(b), diagb(b),lowerb(b), rhsb(b)
    REAL(dp):: uppera(a), diaga(a),lowera(a), rhsa(a), dyc(a), dy(0:a), sus_F(a,b)
    REAL(dp):: ded=.2_dp !Dimensionless eddy diffusivity
    INTEGER:: i, j, info


    !!!First define water depths
    do i= 1,b
    !Depth at previous time step
    depthN(:,i)= max(elevs_old(b)-hs(:,i),0._dp)
    !Depth at half time step
    elevsH=(elevs(b)+elevs_old(b))*.5_dp
    depthH(:,i)= max(elevsH-hs(:,i), 0._dp)
    !Depth at next time step
    depthF(:,i)= max(elevs(b)-hs(:,i),0._dp)
    end do

    !Now calculate half time step velocity
    velh=.5_dp*(vels+vels_old)
    !Halfway depth times velocity - a useful efficiency device
    veldH(:,1:b)=depthH*velH
    veldH(:,0)= veldH(:,1) !B
    veldH(:,b+1)=veldH(:,b) !B


    !!Eddy diffusivity*depth at last time step
    eddN= ded*abs(vels_old)*sqrt(fs)*sqrt(0.125_dp)*depthN*depthN !+10._dp*depthN !Note that the sqrt(0.125_dp) is sqrt(1/8)


    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial y average of eddy diffusivity times the depth, with boundaries
    eddN_yh(1:a-1,:)= .5_dp*(eddN(2:a,:)+eddN(1:a-1,:))
    eddN_yh(0,:)=eddN(1,:) !B
    eddN_yh(a,:)=eddN(a,:) !B


    !!Eddy diffusivity*depth at half time step
    eddH= ded*abs(velh)*sqrt(fs)*sqrt(0.125_dp)*depthH*depthH !Note that the sqrt(0.125_dp) is sqrt(1/8)
    !Halfway spatial x average of eddydiff * depth
    eddH_xh(: , 1:b-1 )= .5_dp*(eddH(:,2:b)+eddH(:,1:b-1)) +.0*DT!+0._dp*.5_dp*(depthH(:,1:b-1)+depthH(:,2:b)) !Note the addittion of an extra diffusivity here - this is needed for the method to be stable - and also, perhaps is justified because of the vertical velocity profiles?
    eddH_xh(:,0)= eddH_xh(:,1) !B
    eddH_xh(:,b)= eddH_xh(:,b-1) !B


    !!Susconc at previous time step, with suitable boundary conditions attached
    sus_N(1:a,:)= Cdist_old
    !Zero susconc along the lateral boundaries
    sus_N(0,:)=sus_N(1,:) !0._dp !B
    sus_N(a+1,:)=sus_N(a,:)!0._dp !B

    !Lengths, with suitable boundaries attached
    lnths(1:a,:)= lengths
    lnths(0,:)= 2._dp*lengths(2,:)-lengths(1,:) !B
    lnths(a+1,:)= 2._dp*lengths(a,:)-lengths(a-1,:) !B

    !!!!!!First half time step - implicit in the longitudinal coord - note how the lateral boundary condition is
    !automatically included - need to include the longitudinal boundary conditions
    !on the RHS
    do j=1,a

    !!Calculate matrix diagonals and RHS - the diagonals have the components of the 2nd deriv term on
    !the end of the expression
    upperb= veldH(j,2:b+1)*.5_dp/delX -1._dp/delX*(eddH_xh(j,1:b)/delX)
    diagb= depthH(j,:)*2._dp/DT +wset +1._dp/delX*(eddH_xh(j,1:b)/delX) + 1._dp/delX*(eddH_xh(j,0:b-1)/delX)
    lowerb= -veldH(j,0:b-1)*.5_dp/delX -1._dp/delX*(eddH_xh(j,0:b-1)/delX)
    rhsb= Qe_old(j,:)*sqrt(1._dp+slopes(j,:)**2)+ depthN(j,:)*sus_N(j,:)*2._dp/DT + 2._dp/(lnths(j+1,:)-lnths(j-1,:))*&  !!
    ( eddN_yh(j,:)*(sus_N(j+1,:)-sus_N(j,:) )/(lnths(j+1,:)-lnths(j,:)) - & 
    eddN_yh(j-1,:)*((sus_N(j,:)-sus_N(j-1,:) )/(lnths(j,:)-lnths(j-1,:)) ) ) 

    !Boundary conditions at the mouth
    if(veldH(j,1)>0._dp) THEN
        rhsb(1)= rhsb(1)+ veldH(j,0)*.5_dp/delX*Cmouth(j) +1._dp/delX*eddH_xh(j,0)*Cmouth(j)/delX !B - this is when the sediment inflows from the 'ocean' 
    else
    !	print*, "outflow condition"
    !	diagb(1)= diagb(1)-(2._dp*veldH(j,1)-veldH(j,2))*.5_dp/delX !B - this is when the flow is outward directed, and we assume a zero gradient in veldH*susconc at the mouth - perhaps this is foolish, but lets try it
        rhsb(1)= rhsb(1) +veldH(j,0)*.5_dp/delX*Cdist_old(j,1) +1._dp/delX*eddH_xh(j,0)*Cdist_old(j,1)/delX
    END IF
        
    !if(veldH(j,b)<0._dp) THEN
    rhsb(b)= rhsb(b) - veldH(j,b+1)*.5_dp/delX*Criver(j) + 1._dp/delX*eddH_xh(j,b)*Criver(j)/delX !B - this is where there is river inflow

    call DGTSV(b,1, lowerb(2:b), diagb, upperb(1:b-1), rhsb, b,info)

    sus_h(j,1:b)= rhsb

    if(info.ne.0) print*, "matrix problem in susconc2d"

    end do

    !print*, "half", sus_h(175,1), sus_h(175,40)
     
    !!!!!!!!!!!!!!End of first half time step

    !!!These boundary conditions are needed for the next half step
    sus_h(:,0)= Cmouth !B
    sus_h(:,b+1)=Criver !B

    !!!Begin second half time step

    !!Eddy diffusivity*depth at last time step
    eddF= ded*abs(vels)*sqrt(fs)*sqrt(0.125_dp)*depthF*depthF !+10._dp*depthF !Note that the sqrt(0.125_dp) is sqrt(1/8) - note also how we ignore the change in the friction factor here


    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times depth, with boundaries
    eddF_yh(1:a-1,:)=.5_dp*(eddF(2:a,:)+eddF(1:a-1,:))
    eddF_yh(0,:)=eddF(1,:) !B
    eddF_yh(a,:)=eddF(a,:) !B

    do i=1,b

    !Define some dy type variables
    dyc(2:a-1)= .5*(lengths(3:a,i)-lengths(1:a-2,i))
    dyc(1)=dyc(2) !B
    dyc(a)=dyc(a-1) !B

    dy(1:a-1)=lengths(2:a,i)-lengths(1:a-1,i)
    dy(a)=dy(a-1)
    dy(0)=dy(1)

    !!Matrix diagonals
    uppera= -1._dp/dyc*(eddF_yh(1:a,i)/dy(1:a))
    lowera= -1._dp/dyc*(eddF_yh(0:a-1,i)/dy(0:a-1))
    diaga= depthF(:,i)/DT*2._dp -uppera - lowera +wset !!

    rhsa= depthH(:,i)*sus_h(:,i)/DT*2._dp + Qe(:,i)*sqrt(1._dp+slopes(:,i)**2) - & 
    (veldH(:,i+1)*sus_h(:,i+1) - veldH(:,i-1)*sus_h(:,i-1))&
    *.5_dp/delX + 1._dp/delX*( eddH_xh(:,i)*(sus_h(:,i+1)-sus_h(:,i))/delX - eddH_xh(:,i-1)*(sus_h(:,i)-sus_h(:,i-1))/delX )
    !For boundary conditions here, it is appropriate that the sus is zero on the
    !lateral boundaries-  which equates to us not having to enforce anything
    !further. 

    call DGTSV(a,1, lowera(2:a), diaga, uppera(1:a-1), rhsa,a, info)

    Cdist(:,i)= rhsa

    if(info.ne.0) print*, "matrix problem in susconc2d"

    if(minval(Cdist(:,i))<0._dp) THEN
    !	if(abs(minval(Cdist(:,i)))< 1.00E-12) THEN
            Cdist=max(Cdist,0._dp)
    !	ElSE
    !	print*, "Cdist negative with non-roundoff type magnitude"
    !	print*, "min Cdist", i, "<0", minval(Cdist(:,i)),minloc(Cdist(:,i)),depthF(minloc(Cdist(:,i)),max((i-1),1):min((i+1),b))& 
    !	, vels(minloc(Cdist(:,i)), max((i-1),1):min((i+1),b)) 
    !	END IF

    end if

    end do



end subroutine susconc2d


!!Numerical method is upwind in the x derivative advective term, but still ADI
!otherwise
subroutine susconc2dup(a,b, DT,delX, lengths,elevs,elevs_old,vels,vels_old,hs, Cdist,Cdist_old, Qe,Qe_old, Cmouth, Criver, & 
wset,fs, slopes)
    !use matrix_solvers
    !use crosssection

    !implicit none
    INTEGER, intent(in):: a, b
    REAL(dp), intent(in):: DT,delX, wset
    REAL(dp), intent(in):: elevs, elevs_old, vels, vels_old,hs, Cdist_old, Qe,Qe_old, Cmouth, Criver,fs, lengths, slopes
    REAL(dp), intent(inout):: Cdist
    DIMENSION elevs(b), elevs_old(b), vels(a,b), vels_old(a,b), hs(a,b), Cdist(a,b), Cdist_old(a,b), Qe(a,b),fs(a,b), lengths(a,b) &
    ,Criver(a), Cmouth(a), Qe_old(a,b), slopes(a,b)

    REAL(dp):: depthN(a,b), depthH(a,b), depthF(a,b), elevsH, velh(a,b),sus_h(a,0:b+1), eddN(a,b), veldh(a,0:b+1)
    REAL(dp):: sus_N(0:a+1,b), lnths(0:a+1,b), eddN_yh(0:a,b) , eddF(a,b), eddF_yh(0:a,b)
    REAL(dp):: upperb(b), diagb(b),lowerb(b), rhsb(b)
    REAL(dp):: uppera(a), diaga(a),lowera(a), rhsa(a), dyc(a), dy(0:a), sus_F(a,b)
    REAL(dp):: ded=.20_dp !Dimensionless eddy diffusivity
    INTEGER:: i, j, info


    !!!!So this solves the equation
    ! d/dt(depth*C) + d/dx(vel*depth*C)= d/dy(eddn*depth*dC/dy)
    !By the ADI method, except that the convective x derivative is treated as upwind
    !(seems that it has to be to get stable results). 


    !!!First define water depths
    do i= 1,b
    !Depth at previous time step
    depthN(:,i)= max(elevs_old(b)-hs(:,i),0._dp)
    !Depth at half time step
    elevsH=(elevs(b)+elevs_old(b))*.5_dp
    depthH(:,i)= max(elevsH-hs(:,i), 0._dp)
    !Depth at next time step
    depthF(:,i)= max(elevs(b)-hs(:,i),0._dp)
    end do

    !Now calculate half time step velocity
    velH=.5_dp*(vels+vels_old)
    !Halfway depth times velocity - a useful efficiency device
    veldH(:,1:b)=depthH*velH
    veldH(:,0)= veldH(:,1) !B
    veldH(:,b+1)=veldH(:,b) !B


    !!Eddy diffusivity*depth at last time step
    eddN= ded*abs(vels_old)*sqrt(fs)*sqrt(0.125_dp)*depthN*depthN !Note that the sqrt(0.125_dp) is sqrt(1/8)


    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times the depth, with boundaries
    eddN_yh(1:a-1,:)= .5_dp*(eddN(2:a,:)+eddN(1:a-1,:))
    eddN_yh(0,:)=.5_dp*eddN(1,:) !B
    eddN_yh(a,:)=.5_dp*eddN(a,:) !B

    !!Susconc at previous time step, with suitable boundary conditions attached
    sus_N(1:a,:)= Cdist_old
    !Zero susconc along the lateral boundaries
    sus_N(0,:)=0._dp !B
    sus_N(a+1,:)=0._dp !B

    !Lengths, with suitable boundaries attached
    lnths(1:a,:)= lengths
    lnths(0,:)= 2._dp*lengths(2,:)-lengths(1,:) !B
    lnths(a+1,:)= 2._dp*lengths(a,:)-lengths(a-1,:) !B

    !!!!!!First half time step - implicit in the longitudinal coord - note how the lateral boundary condition is
    !automatically included - need to include the longitudinal boundary conditions
    !on the RHS
    do j=1,a

    !!Calculate matrix diagonals and RHS - upwind method
    do i=1,b
    if(veldH(j,i)>0._dp) THEN
        upperb(i)= 0._dp !veldH(j,2:b+1)*.5_dp/delX
    diagb(i)= depthH(j,i)*2._dp/DT +wset + veldh(j,i)/delX !!
    lowerb(i)= -veldh(j,i-1)/delX !veldh(j,0:b-1)*.5_dp/delX
    else
    upperb(i)= veldh(j,i+1)/delX !veldh(j,2:b+1)*.5_dp/delX
    diagb(i)= depthH(j,i)*2._dp/DT +wset - veldh(j,i)/delX !!
    lowerb(i)= 0._dp !veldh(j,0:b-1)*.5_dp/delX
    end if
    end do

    rhsb= Qe_old(j,:)*sqrt(1._dp+slopes(j,:)**2)+ depthN(j,:)*sus_N(j,:)*2._dp/DT + 2._dp/(lnths(j+1,:)-lnths(j-1,:))*&  !!
    ( eddN_yh(j,:)*(sus_N(j+1,:)-sus_N(j,:) )/(lnths(j+1,:)-lnths(j,:)) - & 
    eddN_yh(j-1,:)*((sus_N(j,:)-sus_N(j-1,:) )/(lnths(j,:)-lnths(j-1,:)) ) ) 

    !Boundary conditions at the mouth
    if(veldh(j,1)>0._dp) THEN
        rhsb(1)= rhsb(1)+ veldh(j,0)/delX*Cmouth(j) !B - this is when the sediment inflows from the 'ocean' 
    !else
    !	diagb(1)= diagb(1)-veldh(j,1)/delX !B - this is when the flow is outward directed, and we assume a zero gradient in veldH*susconc at the mouth - perhaps this is foolish, but lets try it
    END IF
        
    if(veldh(j,b)<0._dp) THEN
    rhsb(b)= rhsb(b) - veldh(j,b+1)/delX*Criver(j) !B - this is where there is river inflow
    end if

    call DGTSV(b,1, lowerb(2:b), diagb, upperb(1:b-1), rhsb,b, info)

    sus_h(j,1:b)= rhsb

    if(info.ne.0) print*, "matrix problem in susconc2d"

    end do

    !print*, "half", sus_h(175,1), sus_h(175,40)
     
    !!!!!!!!!!!!!!End of first half time step

    !!!These boundary conditions are needed for the next half step
    sus_h(:,0)= Cmouth !B
    sus_h(:,b+1)=Criver !B

    !!!Begin second half time step

    !!Eddy diffusivity*depth at last time step
    eddF= ded*abs(vels)*sqrt(fs)*sqrt(0.125_dp)*depthF*depthF !Note that the sqrt(0.125_dp) is sqrt(1/8) - note also how we ignore the change in the friction factor here


    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times depth, with boundaries
    eddF_yh(1:a-1,:)= .5_dp*(eddF(2:a,:)+eddF(1:a-1,:))
    eddF_yh(0,:)=.5_dp*eddF(1,:) !B
    eddF_yh(a,:)=.5_dp*eddF(a,:) !B

    do i=1,b

    !Define some dy type variables
    dyc(2:a-1)= .5*(lengths(3:a,i)-lengths(1:a-2,i))
    dyc(1)=dyc(2) !B
    dyc(a)=dyc(a-1) !B

    dy(1:a-1)=lengths(2:a,i)-lengths(1:a-1,i)
    dy(a)=dy(a-1)
    dy(0)=dy(1)

    !!Matrix diagonals
    uppera= -1._dp/dyc*(eddF_yh(1:a,i)/dy(1:a))
    lowera= -1._dp/dyc*(eddF_yh(0:a-1,i)/dy(0:a-1))
    diaga= depthF(:,i)/DT*2._dp -uppera - lowera +wset !!

    !Right hand side, upwind method
    do j=1,a
    if(veldh(j,i)>0._dp) THEN
    rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
    (veldh(j,i)*sus_h(j,i) - veldh(j,i-1)*sus_h(j,i-1))/delX !!
    ELSE
    rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
    (veldh(j,i+1)*sus_h(j,i+1) - veldh(j,i)*sus_h(j,i))/delX !!
    END IF
    END DO
    !For boundary conditions here, it is appropriate that the sus is zero on the
    !lateral boundaries-  which equates to us not having to enforce anything
    !further. 

    call DGTSV(a, 1,lowera(2:a), diaga, uppera(1:a-1), rhsa,a, info)

    Cdist(:,i)= rhsa

    if(info.ne.0) print*, "matrix problem in susconc2d"

    if(minval(Cdist(:,i))<0._dp) THEN
            Cdist=max(Cdist,0._dp)
        !	if(abs(minval(Cdist(:,i)))< 1.00E-12) THEN

    !	ElSE
        !print*, "Cdist negative with non-roundoff type magnitude"
        !print*, "min Cdist", i, "<0", minval(Cdist(:,i)), minloc(Cdist(:,i)), depthF(minloc(Cdist(:,i)), (i-1):(i+1))& 
        !, vels(minloc(Cdist(:,i)), (i-1):(i+1)) 
        !END IF

    end if

    end do



end subroutine susconc2dup

!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!




!!Numerical method is upwind in the x derivative advective term, but still ADI
!otherwise
subroutine susconc2duprev(a,b, DT,delX, lengths,elevs,elevs_old,vels,vels_old,hs, Cdist,Cdist_old, Qe,Qe_old, Cmouth, Criver, & 
wset,fs, slopes)
    !use matrix_solvers
    !use crosssection

    !implicit none
    INTEGER, intent(in):: a, b
    REAL(dp), intent(in):: DT,delX, wset
    REAL(dp), intent(in):: elevs, elevs_old, vels, vels_old,hs, Cdist_old, Qe,Qe_old, Cmouth, Criver,fs, lengths, slopes
    REAL(dp), intent(inout):: Cdist
    DIMENSION elevs(b), elevs_old(b), vels(a,b), vels_old(a,b), hs(a,b), Cdist(a,b), Cdist_old(a,b), Qe(a,b),fs(a,b), lengths(a,b) &
    ,Criver(a), Cmouth(a), Qe_old(a,b), slopes(a,b)

    REAL(dp):: depthN(a,b), depthH(a,b), depthF(a,b), elevsH, velh(a,b),sus_h(0:a+1,b), eddH(a,b), veldF(a,0:b+1), sus_N(a,0:b+1)
    REAL(dp):: lnths(0:a+1,b), eddH_yh(0:a,b) , eddF(a,b), eddF_yh(0:a,b), veldN(a,0:b+1)
    REAL(dp):: upperb(b), diagb(b),lowerb(b), rhsb(b)
    REAL(dp):: uppera(a), diaga(a),lowera(a), rhsa(a), dyc(a), dy(0:a), sus_F(a,b)
    REAL(dp):: ded=.20_dp !Dimensionless eddy diffusivity
    INTEGER:: i, j, info
    !!!First define water depths
    do i= 1,b
    !Depth at previous time step
    depthN(:,i)= max(elevs_old(b)-hs(:,i),0._dp)
    !Depth at half time step
    elevsH=(elevs(b)+elevs_old(b))*.5_dp
    depthH(:,i)= max(elevsH-hs(:,i), 0._dp)
    !Depth at next time step
    depthF(:,i)= max(elevs(b)-hs(:,i),0._dp)
    end do


    veldN(:,1:b)=vels_old*depthN
    veldN(:,0)= veldN(:,1) !B
    veldN(:,b+1)= veldN(:,b) !B

    !Now calculate half time step velocity
    velh=.5_dp*(vels+vels_old)
    veldF(:,1:b)= vels*depthF
    veldF(:,b+1)=veldF(:,b) !B
    veldF(:,0)=veldF(:,1) !B

    !!Eddy diffusivity*depth at half time step
    eddH= ded*abs(velh)*sqrt(fs)*sqrt(0.125_dp)*depthH*depthH !Note that the sqrt(0.125_dp) is sqrt(1/8)

    !Halfway spatial average of halfway eddy diffusivity times the depth, with boundaries
    eddH_yh(1:a-1,:)= .5_dp*(eddH(2:a,:)+eddH(1:a-1,:))
    eddH_yh(0,:)=.5_dp*eddH(1,:) !B
    eddH_yh(a,:)=.5_dp*eddH(a,:) !B

    !!Susconc at previous time step, with suitable boundary conditions attached
    sus_N(:,1:b)= Cdist_old
    !Imposed susconc at the boundaries
    sus_N(:,1)=Cmouth !B
    sus_N(:,b+1)=Criver !B

    !Lengths, with suitable boundaries attached
    lnths(1:a,:)= lengths
    lnths(0,:)= 2._dp*lengths(2,:)-lengths(1,:) !B
    lnths(a+1,:)= 2._dp*lengths(a,:)-lengths(a-1,:) !B

    !!!!!!First half time step - implicit in the longitudinal coord - note how the lateral boundary condition is
    !automatically included - need to include the longitudinal boundary conditions
    !on the RHS
    do i=1,b

    !Define some dy type variables
    dyc(2:a-1)= .5*(lengths(3:a,i)-lengths(1:a-2,i))
    dyc(1)=dyc(2) !B
    dyc(a)=dyc(a-1) !B

    dy(1:a-1)=lengths(2:a,i)-lengths(1:a-1,i)
    dy(a)=dy(a-1)
    dy(0)=dy(1)

    !!Matrix diagonals
    uppera= -1._dp/dyc*(eddH_yh(1:a,i)/dy(1:a))
    lowera= -1._dp/dyc*(eddH_yh(0:a-1,i)/dy(0:a-1))
    diaga= depthH(:,i)/DT*2._dp -uppera - lowera +wset !!

    !Right hand side, upwind method
    do j=1,a
    if(veldN(j,i)>0._dp) THEN
    rhsa(j)= depthN(j,i)*sus_N(j,i)/DT*2._dp + Qe_old(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
    (veldN(j,i)*sus_N(j,i) - veldN(j,i-1)*sus_N(j,i-1))/delX !!
    ELSE
    rhsa(j)= depthN(j,i)*sus_N(j,i)/DT*2._dp + Qe_old(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
    (veldN(j,i+1)*sus_N(j,i+1) - veldN(j,i)*sus_N(j,i))/delX !!
    END IF
    END DO
    !For boundary conditions here, it is appropriate that the sus is zero on the
    !lateral boundaries-  which equates to us not having to enforce anything
    !further. 

    call DGTSV(a,1, lowera(2:a), diaga, uppera(1:a-1), rhsa,a, info)

    sus_h(1:a,i)= rhsa

    if(info.ne.0) print*, "matrix problem in susconc2d"

    if(minval(sus_h(:,i))<0._dp) THEN
    !	print*, "halfway prob in susconc2duprev", minval(sus_h(:,i))
        sus_h(1:a,i)= max(sus_h(1:a,i), 0._dp)	
    END IF

    end do

    sus_h(0,:)=sus_h(1,:)
    sus_h(a+1,:)=sus_h(a,:)


    !!!!!!!!!!
    do j=1,a

    !!Calculate matrix diagonals and RHS - upwind method
    do i=1,b
    if(veldF(j,i)>0._dp) THEN
        upperb(i)= 0._dp !veldh(j,2:b+1)*.5_dp/delX
    diagb(i)= depthF(j,i)*2._dp/DT +wset + veldF(j,i)/delX !!
    lowerb(i)= -veldF(j,i-1)/delX !veldh(j,0:b-1)*.5_dp/delX
    else
    upperb(i)= veldF(j,i+1)/delX !veldh(j,2:b+1)*.5_dp/delX
    diagb(i)= depthF(j,i)*2._dp/DT +wset - veldF(j,i)/delX !!
    lowerb(i)= 0._dp !veldh(j,0:b-1)*.5_dp/delX
    end if
    end do

    rhsb= Qe(j,:)*sqrt(1._dp+slopes(j,:)**2)+ depthH(j,:)*sus_H(j,:)*2._dp/DT + 2._dp/(lnths(j+1,:)-lnths(j-1,:))*&  !!
    ( eddH_yh(j,:)*(sus_H(j+1,:)-sus_H(j,:) )/(lnths(j+1,:)-lnths(j,:)) - & 
    eddH_yh(j-1,:)*((sus_H(j,:)-sus_H(j-1,:) )/(lnths(j,:)-lnths(j-1,:)) ) ) 

    !Boundary conditions at the mouth
    if(veldF(j,1)>0._dp) THEN
        rhsb(1)= rhsb(1)+ veldF(j,1)/delX*Cmouth(j) !B - this is when the sediment inflows from the 'ocean' 
    !else
    !	diagb(1)= diagb(1)-veldh(j,1)/delX !B - this is when the flow is outward directed, and we assume a zero gradient in veldH*susconc at the mouth - perhaps this is foolish, but lets try it
    END IF
        
    if(veldF(j,b)<0._dp) THEN
    rhsb(b)= rhsb(b) - veldF(j,b)/delX*Criver(j) !B - this is where there is river inflow
    end if

    call DGTSV(b,1, lowerb(2:b), diagb, upperb(1:b-1), rhsb,b, info)

    if(info.ne.0) print*, "matrix problem in susconc2d"


    Cdist(j,1:b)= rhsb

    END DO


    !!!Just check things
    do i=1,b

    if(minval(Cdist(:,i))<0._dp) THEN

        !	if(abs(minval(Cdist(:,i)))< 1.00E-12) THEN

    !	ElSE
     !     print*, "Cdist negative with non-roundoff type magnitude"
          !print*, "min Cdist", i, "<0", minval(Cdist(:,i)), minloc(Cdist(:,i)),depthF(minloc(Cdist(:,i)),max((i-1),1):min((i+1),b))& 
          !, vels(minloc(Cdist(:,i)), (i-1):(i+1)) 
        !END IF
            Cdist=max(Cdist,0._dp)
    end if

    end do



end subroutine susconc2duprev


!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine susconc2dOS(a,b, DT,delX, lengths,elevs,elevs_old,vels,vels_old,hs, Cdist,Cdist_old, Qe,Qe_old, Cmouth, Criver, & 
wset,fs, slopes)
    !use matrix_solvers
    !use crosssection

    !implicit none
    INTEGER, intent(in):: a, b
    REAL(dp), intent(in):: DT,delX, wset
    REAL(dp), intent(in):: elevs, elevs_old, vels, vels_old,hs, Cdist_old, Qe,Qe_old, Cmouth, Criver,fs, lengths, slopes
    REAL(dp), intent(inout):: Cdist
    DIMENSION elevs(b), elevs_old(b), vels(a,b), vels_old(a,b), hs(a,b), Cdist(a,b), Cdist_old(a,b), Qe(a,b),fs(a,b), lengths(a,b) &
    ,Criver(a), Cmouth(a), Qe_old(a,b), slopes(a,b)

    REAL(dp):: depthN(a,b), depthH(a,b), depthF(a,b), elevsH, velh(a,b),sus_h(a,0:b+1), eddN(a,b), veldF(a,0:b+1)
    REAL(dp):: sus_N(0:a+1,b), lnths(0:a+1,b), eddN_yh(0:a,b) , eddF(a,b), eddF_yh(0:a,b)
    REAL(dp):: upperb(b), diagb(b),lowerb(b), rhsb(b)
    REAL(dp):: uppera(a), diaga(a),lowera(a), rhsa(a), dyc(a), dy(0:a), sus_F(a,b)
    REAL(dp):: ded=.0_dp !Dimensionless eddy diffusivity
    INTEGER:: i, j, info, indx(a), pt

    !!!First define water depths
    do i= 1,b
    !Depth at previous time step
    depthN(:,i)= max(elevs_old(b)-hs(:,i),0._dp)
    !Depth at half time step
    !elevsH=(elevs(b)+elevs_old(b))*.5_dp
    !depthH(:,i)= max(elevsH-hs(:,i), 0._dp)
    !Depth at next time step
    depthF(:,i)= max(elevs(b)-hs(:,i),0._dp)
    end do

    !Now calculate half time step velocity
    !velH=.5_dp*(vels+vels_old)
    !Halfway depth times velocity - a useful efficiency device
    !veldH(:,1:b)=depthH*velH
    !veldH(:,0)= veldH(:,1) !B
    !veldH(:,b+1)=veldH(:,b) !B

    veldF(:,1:b)=depthF*vels
    veldF(:,0)= veldF(:,1) !B
    veldF(:,b+1)=veldF(:,b) !B


    !!Eddy diffusivity*depth at last time step
    !eddN= ded*vels_old*sqrt(fs)*sqrt(0.125_dp)*depthN*depthN !Note that the sqrt(0.125_dp) is sqrt(1/8)


    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times the depth, with boundaries
    !eddN_yh(1:a-1,:)= .5_dp*(eddN(2:a,:)+eddN(1:a-1,:))
    !eddN_yh(0,:)=.5_dp*eddN(1,:) !B
    !eddN_yh(a,:)=.5_dp*eddN(a,:) !B

    !!Susconc at previous time step, with suitable boundary conditions attached
    sus_N(1:a,:)= Cdist_old
    !Zero susconc along the lateral boundaries
    sus_N(0,:)=0._dp !B
    sus_N(a+1,:)=0._dp !B

    !Lengths, with suitable boundaries attached
    lnths(1:a,:)= lengths
    lnths(0,:)= 2._dp*lengths(2,:)-lengths(1,:) !B
    lnths(a+1,:)= 2._dp*lengths(a,:)-lengths(a-1,:) !B

    !!!!!!Advection time step - implicit in the longitudinal coord - note how the lateral boundary condition is
    !automatically included - need to include the longitudinal boundary conditions
    !on the RHS
    do j=1,a

    !!Calculate matrix diagonals and RHS - upwind method
    do i=1,b
    if(veldF(j,i)>0._dp) THEN
        upperb(i)= 0._dp !veldH(j,2:b+1)*.5_dp/delX
    diagb(i)= depthF(j,i)/DT +wset + veldF(j,i)/delX !!
    lowerb(i)= -veldF(j,i-1)/delX !veldh(j,0:b-1)*.5_dp/delX
    else
    upperb(i)= veldF(j,i+1)/delX !veldh(j,2:b+1)*.5_dp/delX
    diagb(i)= depthF(j,i)/DT +wset - veldF(j,i)/delX !!
    lowerb(i)= 0._dp !veldh(j,0:b-1)*.5_dp/delX
    end if
    end do

    rhsb= depthN(j,:)*sus_N(j,:)/DT + Qe(j,:)*sqrt(1._dp+slopes(j,:)**2)  ! + 2._dp/(lnths(j+1,:)-lnths(j-1,:))*&  !!
    !( eddN_yh(j,:)*(sus_N(j+1,:)-sus_N(j,:) )/(lnths(j+1,:)-lnths(j,:)) - & 
    !eddN_yh(j-1,:)*((sus_N(j,:)-sus_N(j-1,:) )/(lnths(j,:)-lnths(j-1,:)) ) ) 

    !Boundary conditions at the mouth
    if(veldF(j,1)>0._dp) THEN
        rhsb(1)= rhsb(1)+ veldF(j,0)/delX*Cmouth(j) !B - this is when the sediment inflows from the 'ocean' 
    !else
    !	diagb(1)= diagb(1)-veldh(j,1)/delX !B - this is when the flow is outward directed, and we assume a zero gradient in veldH*susconc at the mouth - perhaps this is foolish, but lets try it
    END IF
        
    if(veldF(j,b)<0._dp) THEN
    rhsb(b)= rhsb(b) - veldF(j,b+1)/delX*Criver(j) !B - this is where there is river inflow
    end if

    do i=1,b
    if(diagb(i)==0._dp) print*, "diag ", i," eq 0"
    END DO

    call DGTSV(b,1, lowerb(2:b), diagb, upperb(1:b-1), rhsb,b, info)

    sus_h(j,1:b)= rhsb

    if(info.ne.0) print*, info, "matrix problem in susconc2d longitudinal"

    end do


     
    !!!!!!!!!!!!!!End of advection time step

    !!!These boundary conditions are needed for the lateral time step
    sus_h(:,0)= Cmouth !B
    sus_h(:,b+1)=Criver !B

    !!!Begin second half time step

    !!Eddy diffusivity*depth at last time step
    eddF= ded*abs(vels)*sqrt(fs)*sqrt(0.125_dp)*depthF*depthF !Note that the sqrt(0.125_dp) is sqrt(1/8) - note also how we ignore the change in the friction factor here


    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times depth, with boundaries
    eddF_yh(1:a-1,:)= .5_dp*(eddF(2:a,:)+eddF(1:a-1,:)) +0.001_dp
    eddF_yh(0,:)=.5_dp*eddF(1,:) !B
    eddF_yh(a,:)=.5_dp*eddF(a,:) !B

    do i=1,b

    !Define some dy type variables
    dyc(2:a-1)= .5*(lengths(3:a,i)-lengths(1:a-2,i))
    dyc(1)=dyc(2) !B
    dyc(a)=dyc(a-1) !B

    dy(1:a-1)=lengths(2:a,i)-lengths(1:a-1,i)
    dy(a)=dy(a-1)
    dy(0)=dy(1)

    !!Matrix diagonals
    uppera= -1._dp/dyc*(eddF_yh(1:a,i)/dy(1:a))
    lowera= -1._dp/dyc*(eddF_yh(0:a-1,i)/dy(0:a-1))
    diaga= depthF(:,i)/DT -uppera - lowera +0._dp*wset !!

    !Right hand side, upwind method
    do j=1,a
    !if(veldF(j,i)>0._dp) THEN
    rhsa(j)= depthF(j,i)*sus_h(j,i)/DT !+ Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) !- & 
    !(veldh(j,i)*sus_h(j,i) - veldh(j,i-1)*sus_h(j,i-1))/delX !!
    !ELSE
    !rhsa(j)= depthF(j,i)*sus_h(j,i)/DT + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) !-&
    !(veldh(j,i+1)*sus_h(j,i+1) - veldh(j,i)*sus_h(j,i))/delX !!
    !END IF
    END DO
    !For boundary conditions here, it is appropriate that the sus is interpolated on the
    !lateral boundaries-  which equates to us not having to enforce anything
    !further. 

    !!!!!!!!!!!!!!!!!
    !NOW, here we need to remove zero rows and columns, and keep track of those we
    !remove, so that everything is okay. This is a little involved
    indx=0 !Predefine the indices of non-zero diagonal points
    pt=0 !Predefine counter which records the number of entries we go through
    do j=1,a
    if(diaga(j).ne.0._dp) THEN !Add an entry to indx
        pt=pt+1	
        indx(pt)=j
    END IF	
    END DO

    Cdist(:,i)=0._dp
    if(pt>1) THEN
    !	print*, pt

    lowera(2:pt)= lowera(indx(2:pt))
    diaga(1:pt)= diaga(indx(1:pt))
    uppera(1:pt-1)= uppera(indx(1:pt-1))
    rhsa(1:pt)= rhsa(indx(1:pt))

    call DGTSV(pt,1, lowera(2:pt), diaga(1:pt), uppera( 1:pt-1), rhsa(1:pt),pt, info)
    Cdist(indx(1:pt),i)= rhsa(1:pt)

    if(info.ne.0) print*, info, "matrix problem in susconc2d lateral"
    ELSE
        if(pt==1) Cdist(indx(1),i)=sus_h(indx(1),i)
    END IF


    if(minval(Cdist(:,i))<0._dp) THEN
            Cdist=max(Cdist,0._dp)
        !	if(abs(minval(Cdist(:,i)))< 1.00E-12) THEN

    !	ElSE
        !print*, "Cdist negative with non-roundoff type magnitude"
        !print*, "min Cdist", i, "<0", minval(Cdist(:,i)), minloc(Cdist(:,i)), depthF(minloc(Cdist(:,i)), (i-1):(i+1))& 
        !, vels(minloc(Cdist(:,i)), (i-1):(i+1)) 
        !END IF

    end if

    end do



end subroutine susconc2dOS

!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!

!!Numerical method is upwind in the x derivative advective term, but still ADI
!otherwise
subroutine susconc2dup2(a,b, DT,delX, lengths,elevs,elevs_old,vels,vels_old,hs, Cdist,Cdist_old, Qe,Qe_old, Cmouth, Criver, & 
wset,fs, slopes)
    !use matrix_solvers
    !use crosssection

    !implicit none
    INTEGER, intent(in):: a, b
    REAL(dp), intent(in):: DT,delX, wset
    REAL(dp), intent(in):: elevs, elevs_old, vels, vels_old,hs, Qe,Qe_old, Cmouth, Criver,fs, lengths, slopes
    REAL(dp), intent(in out):: Cdist, Cdist_old !!Cdist_old will output random interesting information
    DIMENSION elevs(b), elevs_old(b), vels(a,b), vels_old(a,b), hs(a,b), Cdist(a,b), Cdist_old(a,b), Qe(a,b),fs(a,b), lengths(a,b) &
    ,Criver(a), Cmouth(a), Qe_old(a,b), slopes(a,b)

    REAL(dp):: depthN(a,b), depthH(a,b), depthF(a,b), elevsH, velh(a,b),sus_h(a,0:b+1), edN(0:a+1,b), eddN(a,b), veldh(a,0:b+1)
    REAL(dp):: sus_N(0:a+1,b), lnths(0:a+1,b), eddN_yh(0:a,b) ,edF(0:a+1,b), eddF(a,b), eddF_yh(0:a,b)
    REAL(dp):: upperb(b), diagb(b),lowerb(b), rhsb(b)
    REAL(dp):: uppera(a), diaga(a),lowera(a), rhsa(a), dyc(a), dy(0:a), sus_F(a,b), cbedS(0:a+1,b), slopesS(0:a+1,b)
    REAL(dp):: ded=.20_dp !Dimensionless eddy diffusivity
    REAL(dp):: tmpa(a), tmpb(b),tmp
    INTEGER:: i, j, info, ii
    logical:: periodicb=.true.

    !!!!So this solves the equation
    ! d/dt(depth*C) + d/dx(vel*depth*C)= d/dy(eddn*depth*dC/dy + eddn*Cbed*dh/dy)
    ! where I plan to add the last term soon.
    !By the ADI method, except that the convective x derivative is treated as upwind
    !(seems that it has to be to get stable results). 
    !!Note that assuming a constant vertical eddy diffusivity, 
    !Cbed= Caverage*depth/((ez/ws)*(1-exp(-vs/ez*(depth) ) ))
    !!And I have been assuming that ez=0.5_dp*ey


    !!!First define water depths
    do i= 1,b
    !Depth at previous time step
    depthN(:,i)= max(elevs_old(i)-hs(:,i),0._dp)
    !Depth at half time step
    elevsH=(elevs(i)+elevs_old(i))*.5_dp
    depthH(:,i)= max(elevsH-hs(:,i), 0._dp)
    !Depth at next time step
    depthF(:,i)= max(elevs(i)-hs(:,i),0._dp)
    end do

    !Now calculate half time step velocity
    velH=.5_dp*(vels+vels_old)
    !Halfway depth times velocity - a useful efficiency device
    veldH(:,1:b)=depthH*velH
    veldH(:,0)= veldH(:,1) !B
    veldH(:,b+1)=veldH(:,b) !B

    !Eddy diffusivity
    edN(1:a,:)= ded*abs(vels_old)*sqrt(fs)*sqrt(0.125_dp)*depthN

    if(periodicb.eqv..false.) THEN
    edN(0,:)=0._dp
    edN(a+1,:)=0._dp
    ELSE
    edN(0,:)=edN(a,:)
    edN(a+1,:)=edN(1,:)
    end if

    !!Eddy diffusivity*depth at last time step
    eddN= edN(1:a,1:b)*depthN !Note that the sqrt(0.125_dp) is sqrt(1/8)



    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            
            if(edN(j,i)>0._dp) THEN
           tmp= ((.5_dp*edN(j,i))*(1._dp-exp(-wset/(.5_dp*edN(j,i))*depthN(j,i)) ) )
            ELSE
            tmp=0._dp
            END IF

        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthN(j,i)*wset/tmp, 20._dp)
        else
            cbedS(j,i)=20._dp  !So this is the maximum scaling factor we offer
            if(wset==0._dp) cbedS(j,i)=1._dp  !So the idea is that this applies  when wset=0
        
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp


    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times the depth, with boundaries
    eddN_yh(1:a-1,:)= .5_dp*(eddN(2:a,:)+eddN(1:a-1,:))

    if(periodicb.eqv..false.) THEN
    eddN_yh(0,:)=.5_dp*eddN(1,:) !B
    eddN_yh(a,:)=.5_dp*eddN(a,:) !B
    else
    eddN_yh(0,:)=.5_dp*(eddN(1,:)+eddN(a,:)) !B
    eddN_yh(a,:)=.5_dp*(eddN(a,:)+eddN(1,:)) !B
    end if



    !!Susconc at previous time step, with suitable boundary conditions attached
    sus_N(1:a,:)= Cdist_old
    !Zero susconc along the lateral boundaries -- might not be true - zero
    !integrated concentration, but the concentration might not drop off entirely
    sus_N(0,:)=sus_N(1,:) !0._dp !B
    sus_N(a+1,:)=sus_N(a,:) !0._dp !B

    !print*, '0th', sus_N(1,110), sus_N(a,110)


    !Lengths, with suitable boundaries attached
    lnths(1:a,:)= lengths
    lnths(0,:)= 2._dp*lengths(1,:)-lengths(2,:) !B
    lnths(a+1,:)= 2._dp*lengths(a,:)-lengths(a-1,:) !B

    slopesS(1:a,:)=slopes
    slopesS(0,:)=slopes(1,:) !0._dp
    slopesS(a+1,:)=slopes(a,:) !0._dp

    !!!!!!First half time step - implicit in the longitudinal coord - note how the lateral boundary condition is
    !automatically included - need to include the longitudinal boundary conditions
    !on the RHS
    do j=1,a

    !!Calculate matrix diagonals and RHS - upwind method
    do i=1,b
    if(veldH(j,i)>0._dp) THEN
            upperb(i)= 0._dp !veldH(j,2:b+1)*.5_dp/delX
            diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) + veldh(j,i)/delX !!
            lowerb(i)= -veldh(j,i-1)/delX !veldh(j,0:b-1)*.5_dp/delX
    else
            upperb(i)= veldh(j,i+1)/delX !veldh(j,2:b+1)*.5_dp/delX
            diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) - veldh(j,i)/delX !!
            lowerb(i)= 0._dp !veldh(j,0:b-1)*.5_dp/delX
    end if
    end do

    rhsb= Qe_old(j,:)*sqrt(1._dp+slopes(j,:)**2)+ depthN(j,:)*sus_N(j,:)*2._dp/DT + 2._dp/(lnths(j+1,:)-lnths(j-1,:))*&  !!
    ( eddN_yh(j,:)*(sus_N(j+1,:)-sus_N(j,:) )/(lnths(j+1,:)-lnths(j,:)) - & 
    eddN_yh(j-1,:)*((sus_N(j,:)-sus_N(j-1,:) )/(lnths(j,:)-lnths(j-1,:)) ) +&
    ((edN(j+1,:)+edN(j,:))*.5_dp*(cbedS(j+1,:)*sus_N(j+1,:)+cbedS(j,:)*sus_N(j,:))*.5_dp*(slopesS(j+1,:)+slopesS(j,:))*.5_dp - &
     (edN(j,:)+edN(j-1,:))*.5_dp*(cbedS(j,:)*sus_N(j,:)+cbedS(j-1,:)*sus_N(j-1,:))*.5_dp*(slopesS(j,:)+slopesS(j-1,:))*.5_dp) &
    )

    !(eddN_yh(j,:)*cbedS(j,:) ) )

    !Boundary conditions at the mouth
    if(veldh(j,1)>0._dp) THEN
            rhsb(1)= rhsb(1)+ veldh(j,0)/delX*Cmouth(j) !B - this is when the sediment inflows from the 'ocean' 
    !else
    !	diagb(1)= diagb(1)-veldh(j,1)/delX !B - this is when the flow is outward directed, and we assume a zero gradient in veldH*susconc at the mouth - perhaps this is foolish, but lets try it
    END IF
        
    if(veldh(j,b)<0._dp) THEN
    rhsb(b)= rhsb(b) - veldh(j,b+1)/delX*Criver(j) !B - this is where there is river inflow
    end if

    !if((j==1).or.(j==350)) print*, lnths(j+1,110)-lnths(j,110),lnths(j,110)-lnths(j-1,110) ! eddN_yh(j,110)*(sus_N(j+1,110)-sus_N(j,110) )/(lnths(j+1,110)-lnths(j,110)) - & 
    !eddN_yh(j-1,110)*((sus_N(j,110)-sus_N(j-1,110) )/(lnths(j,110)-lnths(j-1,110)) ) 
    ! ((edN(j+1,110)+edN(j,110))*.5_dp*(cbedS(j+1,110)*sus_N(j+1,110)+cbedS(j,110)*sus_N(j,110))& 
    !*.5_dp*(slopesS(j+1,110)+slopesS(j,110))*.5_dp - &
    ! (edN(j,110)+edN(j-1,110))*.5_dp*(cbedS(j,110)*sus_N(j,110)+cbedS(j-1,110)*sus_N(j-1,110))& 
    !*.5_dp*(slopesS(j,110)+slopesS(j-1,110))*.5_dp) !sus_N(j+1,110)-sus_N(j,110), sus_N(j,110)-sus_N(j-1,110)

    !!Prevent main diagonal from being zero
    do ii=1,b
    if(diagb(ii)==0._dp) THEN
            if( (ii==1).or.(ii==b)) THEN
            diagb(ii)=1._dp
            lowerb(ii)=0._dp
            upperb(ii)=0._dp
            rhsb(ii)=0._dp
            else
            diagb(ii)=1._dp
            lowerb(ii)=-0.5_dp
            upperb(ii)=-0.5_dp
            rhsb(ii)=0._dp
            end if
    end if
    end do

    call DGTSV(b,1, lowerb(2:b), diagb, upperb(1:b-1), rhsb,b, info)

    sus_h(j,1:b)= rhsb

    if(info.ne.0) print*, "matrix problem in susconc2d, first bit", info

    end do

    !print*, sus_h(1,110), sus_h(350,110)

    !print*, "half", sus_h(175,1), sus_h(175,40)
     
    !!!!!!!!!!!!!!End of first half time step

    !!!These boundary conditions are needed for the next half step
    sus_h(:,0)= Cmouth !B
    sus_h(:,b+1)=Criver !B

    !!!Begin second half time step

    !!Eddy diffusivity
    edF(1:a,:)=ded*abs(vels)*sqrt(fs)*sqrt(0.125_dp)*depthF

    if(periodicb.eqv..false.) THEN
    edF(0,:)=0._dp
    edF(a+1,:)=0._dp
    ELSE
    edF(0,:)=edF(1,:)
    edF(a+1,:)=edF(a,:)
    END IF
    !!Eddy diffusivity*depth at last time step
    eddF= edF(1:a,:)*depthF !Note that the sqrt(0.125_dp) is sqrt(1/8) - note also how we ignore the change in the friction factor here

    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !cbedS= depthN*wset/((eddN/depthN)*(1._dp-exp(-wset/(eddN/depthN)*depthN) ) )
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            if(wset>0._dp) THEN
            tmp = ((.5_dp*edF(j,i))*(1._dp-exp(-wset/(.5_dp*edF(j,i))*depthF(j,i) ) ))
            else
            tmp=0._dp
            end if
        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthF(j,i)*wset/tmp, 20._dp)
        else
        cbedS(j,i)=20._dp
        if(wset==0) cbedS(j,i)=1._dp
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp


    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times depth, with boundaries
    eddF_yh(1:a-1,:)= .5_dp*(eddF(2:a,:)+eddF(1:a-1,:))
    if(periodicb.eqv..false.) THEN
    eddF_yh(0,:)=.5_dp*eddF(1,:) !B
    eddF_yh(a,:)=.5_dp*eddF(a,:) !B
    else
    eddF_yh(0,:)=.5_dp*(eddF(1,:)+eddF(a,:)) !B
    eddF_yh(a,:)=.5_dp*(eddF(a,:)+eddF(1,:)) !B
    end if


    do i=1,b

    !Define some dy type variables
    dyc(2:a-1)= .5*(lengths(3:a,i)-lengths(1:a-2,i))
    dyc(1)=dyc(2) !B
    dyc(a)=dyc(a-1) !B

    dy(1:a-1)=lengths(2:a,i)-lengths(1:a-1,i)
    dy(a)=dy(a-1)
    dy(0)=dy(1)

    !!Matrix diagonals
    uppera= -1._dp/dyc*(eddF_yh(1:a,i)/dy(1:a)+.5_dp*edF(2:a+1,i)*cbedS(2:a+1,i)*slopesS(2:a+1,i)  )
    lowera= -1._dp/dyc*(eddF_yh(0:a-1,i)/dy(0:a-1)-.5_dp*edF(0:a-1,i)*cbedS(0:a-1,i)*slopesS(0:a-1,i))
    diaga= depthF(:,i)/DT*2._dp +&
    1._dp/dyc*(eddF_yh(1:a,i)/dy(1:a)-.5_dp*edF(1:a,i)*cbedS(1:a,i)*slopesS(1:a,i)  ) +&
    1._dp/dyc*(eddF_yh(0:a-1,i)/dy(0:a-1)+.5_dp*edF(1:a,i)*cbedS(1:a,i)*slopesS(1:a,i)) + &
    wset*cbedS(1:a,i) !!

    !Right hand side, upwind method
    do j=1,a
    if(veldh(j,i)>0._dp) THEN
    rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
    (veldh(j,i)*sus_h(j,i) - veldh(j,i-1)*sus_h(j,i-1))/delX !!
    ELSE
    rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
    (veldh(j,i+1)*sus_h(j,i+1) - veldh(j,i)*sus_h(j,i))/delX !!
    END IF
    END DO
    !For boundary conditions here, it is appropriate that the sus is zero on the
    !lateral boundaries-  which equates to us not having to enforce anything
    !further. Unless we have periodic boundaries, in which case the following might
    !be a practical alternative
    if(periodicb) THEN
    rhsa(1)=rhsa(1) -lowera(a)*sus_h(1,b)
    rhsa(a)= rhsa(a)-uppera(a)*sus_h(a,b)
    END IF
     

    !!Prevent main diagonal from being zero
    do ii=1,a
    if(diaga(ii)==0._dp) THEN
            if( (ii==1).or.(ii==a)) THEN
            diaga(ii)=1._dp
            lowera(ii)=0._dp
            uppera(ii)=0._dp
            rhsa(ii)=0._dp
            else
            diaga(ii)=1._dp
            lowera(ii)=-0.5_dp
            uppera(ii)=-0.5_dp
            rhsa(ii)=0._dp
            end if
    end if
    end do

    call DGTSV(a,1, lowera(2:a), diaga, uppera(1:a-1), rhsa,a, info)

    Cdist(:,i)= rhsa

    if(info.ne.0) print*, "matrix problem in susconc2d, second bit", info

    if(minval(Cdist(:,i))<0._dp) THEN
            Cdist=max(Cdist,0._dp)
        !	if(abs(minval(Cdist(:,i)))< 1.00E-12) THEN

    !	ElSE
        !print*, "Cdist negative with non-roundoff type magnitude"
        !print*, "min Cdist", i, "<0", minval(Cdist(:,i)), minloc(Cdist(:,i)), depthF(minloc(Cdist(:,i)), (i-1):(i+1))& 
        !, vels(minloc(Cdist(:,i)), (i-1):(i+1)) 
        !END IF

    end if


    Cdist_old=cbedS(1:a,:) !At the moment I am just using this to check the output. 

    end do


    !print*, '2nd', Cdist(1,110), Cdist(350,110)

    !stop

end subroutine susconc2dup2




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!Numerical method is upwind in the x derivative advective term, but still ADI
!otherwise
subroutine susconc2dup3(a,b, DT,delX, lengths,elevs,elevs_old,vels,vels_old,hs, Cdist,Cdist_old, Qe,Qe_old, Cmouth, Criver, & 
wset,fs, slopes, l, u)
    !use matrix_solvers
    !use crosssection

    !implicit none
    INTEGER, intent(in):: a, b, l, u
    REAL(dp), intent(in):: DT,delX, wset
    REAL(dp), intent(in):: elevs, elevs_old, vels, vels_old,hs, Qe,Qe_old, Cmouth, Criver,fs, lengths, slopes
    REAL(dp), intent(in out):: Cdist, Cdist_old !!Cdist_old will output random interesting information
    DIMENSION elevs(b), elevs_old(b), vels(a,b), vels_old(a,b), hs(a,b), Cdist(a,b), Cdist_old(a,b), Qe(a,b),fs(a,b), lengths(a,b) &
    ,Criver(a), Cmouth(a), Qe_old(a,b), slopes(a,b), l(b), u(b)

    REAL(dp):: depthN(a,b), depthH(a,b), depthF(a,b), elevsH, velh(a,b),sus_h(a,0:b+1), edN(0:a+1,b), eddN(a,b), veldh(a,0:b+1)
    REAL(dp):: sus_N(0:a+1,b), lnths(0:a+1,b), eddN_yh(0:a,b) ,edF(0:a+1,b), eddF(a,b), eddF_yh(0:a,b)
    REAL(dp):: upperb(b), diagb(b),lowerb(b), rhsb(b), depthN2(0:a+1,b)
    REAL(dp):: uppera(a), diaga(a),lowera(a), rhsa(a), dyc(a), dy(0:a), sus_F(a,b), cbedS(0:a+1,b), slopesS(0:a+1,b)
    REAL(dp):: ded=.20_dp !Dimensionless eddy diffusivity
    REAL(dp):: tmpa(a), tmpb(b),tmp
    INTEGER:: i, j, info, ii
    logical:: periodicb=.false.!false.

    !!!!So this solves the equation
    ! d/dt(depth*C) + d/dx(vel*depth*C)= d/dy(eddn*depth*dC/dy + eddn*Cbed*dh/dy)
    ! where I plan to add the last term soon.
    !By the ADI method, except that the convective x derivative is treated as upwind
    !(seems that it has to be to get stable results). 
    !!Note that assuming a constant vertical eddy diffusivity, 
    !Cbed= Caverage*depth/((ez/ws)*(1-exp(-vs/ez*(depth) ) ))
    !!And I have been assuming that ez=0.5_dp*ey


    !!!First define water depths
    do i= 1,b
    !Depth at previous time step
    depthN(:,i)= max(elevs_old(i)-hs(:,i),0._dp)
    !Depth at half time step
    elevsH=(elevs(i)+elevs_old(i))*.5_dp
    depthH(:,i)= max(elevsH-hs(:,i), 0._dp)
    !Depth at next time step
    depthF(:,i)= max(elevs(i)-hs(:,i),0._dp)
    end do

    !!Useful thing
    depthN2(1:a,:)=depthN
    depthN2(0,:)=0._dp
    depthN2(a+1,:)=0._dp
    !!


    !Now calculate half time step velocity
    velH=.5_dp*(vels+vels_old)
    !Halfway depth times velocity - a useful efficiency device
    veldH(:,1:b)=depthH*velH
    veldH(:,0)= veldH(:,1) !B
    veldH(:,b+1)=veldH(:,b) !B

    !Eddy diffusivity
    edN(1:a,:)= ded*abs(vels_old)*sqrt(fs)*sqrt(0.125_dp)*depthN

    if(periodicb.eqv..false.) THEN
    edN(0,:)=0._dp
    edN(a+1,:)=0._dp

    !Make sure things cancel at the boundaries appropriately
    !do ii=1,b
    !edN(l(ii)-1, ii)=-edN(l(ii),ii)
    !edN(u(ii)+1,ii)=-edN(u(ii),ii)
    !end do
    !edN(0,:)=-edN(1,:)
    !edN(a+1,:)=-edN(a,:)
    ELSE
    edN(0,:)=edN(a,:)
    edN(a+1,:)=edN(1,:)
    end if

    !!Eddy diffusivity*depth at last time step
    eddN= edN(1:a,1:b)*depthN !Note that the sqrt(0.125_dp) is sqrt(1/8)



    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            
            if(edN(j,i)>0._dp) THEN
           tmp= ((.5_dp*edN(j,i))*(1._dp-exp(-wset/(.5_dp*edN(j,i))*depthN(j,i)) ) )
            ELSE
            tmp=0._dp
            END IF

        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthN(j,i)*wset/tmp, 200._dp)
        else
            cbedS(j,i)=200._dp  !So this is the maximum scaling factor we offer
            if((wset==0._dp).and.(edN(j,i)> 0._dp)) cbedS(j,i)=1._dp  !So the idea is that this applies  when wset=0
        
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp


    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times the depth, with boundaries
    eddN_yh(1:a-1,:)= .5_dp*(eddN(2:a,:)+eddN(1:a-1,:))

    if(periodicb.eqv..false.) THEN
    eddN_yh(0,:)=.5_dp*eddN(1,:) !B
    eddN_yh(a,:)=.5_dp*eddN(a,:) !B
    else
    eddN_yh(0,:)=.5_dp*(eddN(1,:)+eddN(a,:)) !B
    eddN_yh(a,:)=.5_dp*(eddN(a,:)+eddN(1,:)) !B
    end if



    !!Susconc at previous time step, with suitable boundary conditions attached
    sus_N(1:a,:)= Cdist_old
    !Zero susconc along the lateral boundaries -- might not be true - zero
    !integrated concentration, but the concentration might not drop off entirely
    sus_N(0,:)=sus_N(1,:) !0._dp !B
    sus_N(a+1,:)=sus_N(a,:) !0._dp !B

    !print*, '0th', sus_N(1,110), sus_N(a,110)


    !Lengths, with suitable boundaries attached
    lnths(1:a,:)= lengths
    lnths(0,:)= 2._dp*lengths(1,:)-lengths(2,:)!lengths(1,:)-depthN(1,:)/max(abs(slopes(1,:)), 0.01_dp)!2._dp*lengths(1,:)-lengths(2,:) !B
    lnths(a+1,:)=  2._dp*lengths(a,:)-lengths(a-1,:)!lengths(a,:) + depthN(a,:)/max(abs(slopes(a,:)), 0.01_dp)!2._dp*lengths(a,:)-lengths(a-1,:) !B

    slopesS(1:a,:)=slopes
    slopesS(0,:)=slopes(1,:) !0._dp
    slopesS(a+1,:)=slopes(a,:) !0._dp

    !!!!!!First half time step - implicit in the longitudinal coord - note how the lateral boundary condition is
    !automatically included - need to include the longitudinal boundary conditions
    !on the RHS
    do j=1,a
    !print*, j, size(Qe_old)
    !!Calculate matrix diagonals and RHS - upwind method
    do i=1,b
    if(veldH(j,i)>0._dp) THEN
            upperb(i)= 0._dp !veldH(j,2:b+1)*.5_dp/delX
            diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) + veldh(j,i)/delX !!
            lowerb(i)= -veldh(j,i-1)/delX !veldh(j,0:b-1)*.5_dp/delX
    else
            upperb(i)= veldh(j,i+1)/delX !veldh(j,2:b+1)*.5_dp/delX
            diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) - veldh(j,i)/delX !!
            lowerb(i)= 0._dp !veldh(j,0:b-1)*.5_dp/delX
    end if
    end do
    !print*, j
    rhsb= Qe_old(j,:)*sqrt(1._dp+slopes(j,:)**2)+ depthN(j,:)*sus_N(j,:)*2._dp/DT + 2._dp/(lnths(j+1,:)-lnths(j-1,:))*&  !!
    ( 0.5_dp*(edN(j,:)+edN(j+1,:))*(depthN2(j+1,:)*sus_N(j+1,:)-depthN2(j,:)*sus_N(j,:) )/(lnths(j+1,:)-lnths(j,:)) - & 
    0.5_dp*(edN(j,:)+edN(j-1,:))*(depthN2(j,:)*sus_N(j,:)-depthN2(j-1,:)*sus_N(j-1,:) )/(lnths(j,:)-lnths(j-1,:))  +&
    ((edN(j+1,:)+edN(j,:))*.5_dp*(cbedS(j+1,:)*sus_N(j+1,:)+cbedS(j,:)*sus_N(j,:))*.5_dp*(slopesS(j+1,:)+slopesS(j,:))*.5_dp - &
     (edN(j,:)+edN(j-1,:))*.5_dp*(cbedS(j,:)*sus_N(j,:)+cbedS(j-1,:)*sus_N(j-1,:))*.5_dp*(slopesS(j,:)+slopesS(j-1,:))*.5_dp) &
    !(.5_dp*(edN(j+1,:)*cbedS(j+1,:)*sus_N(j+1,:)*slopesS(j+1,:)+edN(j,:)*cbedS(j,:)*sus_N(j,:)*slopesS(j,:)) - &
    ! (.5_dp*(edN(j,:)*cbedS(j,:)*sus_N(j,:)*slopesS(j,:)+edN(j-1,:)*cbedS(j-1,:)*sus_N(j-1,:))*slopesS(j-1,:))) &
    )

    !(eddN_yh(j,:)*cbedS(j,:) ) )
    !Boundary conditions at the mouth
    if(veldh(j,1)>0._dp) THEN
            rhsb(1)= rhsb(1)+ veldh(j,0)/delX*Cmouth(j) !B - this is when the sediment inflows from the 'ocean' 
    !else
    !	diagb(1)= diagb(1)-veldh(j,1)/delX !B - this is when the flow is outward directed, and we assume a zero gradient in veldH*susconc at the mouth - perhaps this is foolish, but lets try it
    END IF
        
    if(veldh(j,b)<0._dp) THEN
    rhsb(b)= rhsb(b) - veldh(j,b+1)/delX*Criver(j) !B - this is where there is river inflow
    end if

    tmp=rhsb(b)
    !if((j==1).or.(j==350)) print*, lnths(j+1,110)-lnths(j,110),lnths(j,110)-lnths(j-1,110) ! eddN_yh(j,110)*(sus_N(j+1,110)-sus_N(j,110) )/(lnths(j+1,110)-lnths(j,110)) - & 
    !eddN_yh(j-1,110)*((sus_N(j,110)-sus_N(j-1,110) )/(lnths(j,110)-lnths(j-1,110)) ) 
    ! ((edN(j+1,110)+edN(j,110))*.5_dp*(cbedS(j+1,110)*sus_N(j+1,110)+cbedS(j,110)*sus_N(j,110))& 
    !*.5_dp*(slopesS(j+1,110)+slopesS(j,110))*.5_dp - &
    ! (edN(j,110)+edN(j-1,110))*.5_dp*(cbedS(j,110)*sus_N(j,110)+cbedS(j-1,110)*sus_N(j-1,110))& 
    !*.5_dp*(slopesS(j,110)+slopesS(j-1,110))*.5_dp) !sus_N(j+1,110)-sus_N(j,110), sus_N(j,110)-sus_N(j-1,110)

    !!Prevent main diagonal from being zero
    do ii=1,b
    if((diagb(ii)==0._dp).or.(depthH(j,ii)==0._dp)) THEN
            if( (ii==1).or.(ii==b).or.(depthH(j,ii)==0._dp)) THEN
            diagb(ii)=1._dp
            lowerb(ii)=0._dp
            upperb(ii)=0._dp
            rhsb(ii)=0._dp
            else
            diagb(ii)=1._dp
            lowerb(ii)=-0.5_dp
            upperb(ii)=-0.5_dp
            rhsb(ii)=0._dp
            end if
    end if
    end do

    call DGTSV(b, 1,lowerb(2:b), diagb, upperb(1:b-1), rhsb, b, info)

    sus_h(j,1:b)= rhsb


    if(info.ne.0) print*, "matrix problem in susconc2d, first bit", info

            do ii=1,b
            if(isnan(sus_h(j,ii))) THEN
            print*, 'sus_h(',j,ii,') is NAN' 
            stop
            end if

            if(depthH(j,ii)==0) THEN
                    if(sus_h(j,ii)>0._dp) THEN
                           if(sus_h(j,ii)>0.01_dp) THEN
                            print*, 'zero depth, high sus_h', sus_h(j,ii), j, ii,& 
     depthH(j,ii)*2._dp/DT +wset*cbedS(j,ii) + veldh(j,ii)/delX, veldh(j,i+1)/delX, veldh(j,i-1)/delX 
                            stop
                            end if

                            sus_h(j,ii)=0._dp
                   end if
            end if
            end do

            if((sus_h(j,b).ne.Criver(j)).and.(depthH(j,b)>0._dp)) THEN
    !               print*, 'sus_h(',j,"b) not equal to Criver", sus_h(j,b),Criver(j), depthH(j,b), depthN(j,b)&
    !,depthH(j,b)*2._dp/DT +wset*cbedS(j,b) - veldh(j,b)/delX,veldh(j,b)/delX, tmp, sus_N(j,b)
            end if
    end do
    !stop
    !print*, sus_h(1,110), sus_h(350,110)

    !print*, "half", sus_h(175,1), sus_h(175,40)
     
    !!!!!!!!!!!!!!End of first half time step

    !!!These boundary conditions are needed for the next half step
    sus_h(:,0)= Cmouth !B
    sus_h(:,b+1)=Criver !B

    !!!Begin second half time step

    !!Eddy diffusivity
    edF(1:a,:)=ded*abs(vels)*sqrt(fs)*sqrt(0.125_dp)*depthF

    if(periodicb.eqv..false.) THEN
    edF(0,:)=0._dp
    edF(a+1,:)=0._dp
    !Make things cancel at the boundary
    !do ii=1,b
    !edF(l(ii)-1,ii)=-edF(l(ii),ii)
    !edF(u(ii)+1,ii)=-edF(u(ii),ii)
    !end do

    !edF(0,:)=-edF(1,:) !0._dp
    !edF(a+1,:)=-edF(a,:)!0._dp
    ELSE
    edF(0,:)=edF(1,:)
    edF(a+1,:)=edF(a,:)
    END IF
    !!Eddy diffusivity*depth at last time step
    eddF= edF(1:a,:)*depthF !Note that the sqrt(0.125_dp) is sqrt(1/8) - note also how we ignore the change in the friction factor here

    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !cbedS= depthN*wset/((eddN/depthN)*(1._dp-exp(-wset/(eddN/depthN)*depthN) ) )
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            if(wset>0._dp) THEN
            tmp = ((.5_dp*edF(j,i))*(1._dp-exp(-wset/(.5_dp*edF(j,i))*depthF(j,i) ) ))
            else
            tmp=0._dp
            end if
        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthF(j,i)*wset/tmp, 200._dp)
        else
        cbedS(j,i)=200._dp
        if((wset==0).and.(edF(j,i)> 0._dp)) cbedS(j,i)=1._dp
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp

    !print*, cbedS(175,75), cbedS(175,75)*sus_h(175,75)
    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times depth, with boundaries
    eddF_yh(1:a-1,:)= .5_dp*(eddF(2:a,:)+eddF(1:a-1,:))
    if(periodicb.eqv..false.) THEN
    eddF_yh(0,:)=.5_dp*eddF(1,:) !B
    eddF_yh(a,:)=.5_dp*eddF(a,:) !B
    else
    eddF_yh(0,:)=.5_dp*(eddF(1,:)+eddF(a,:)) !B
    eddF_yh(a,:)=.5_dp*(eddF(a,:)+eddF(1,:)) !B
    end if


    do i=1,b

    !Define some dy type variables
    dyc(2:a-1)= .5_dp*(lengths(3:a,i)-lengths(1:a-2,i))
    dyc(1)=dyc(2) !0.5_dp*(lengths(2,i)-(lengths(1,i)-depthF(1,i)/max(abs(slopes(1,i)), 0.01_dp)))!dyc(2) !B
    dyc(a)=dyc(a-1)!0.5_dp*((lengths(a,i)+depthF(a,i)/max(abs(slopes(a,i)), 0.01_dp))-lengths(a-1,i))!dyc(a-1) !B

    dy(1:a-1)=lengths(2:a,i)-lengths(1:a-1,i)
    dy(a)=dy(a-1) !depthF(a,i)/max(abs(slopes(a,i)), 0.01_dp) !dy(a-1)
    dy(0)=dy(1) !depthF(1,i)/max(abs(slopes(1,i)), 0.01_dp) !dy(1)

    !!Matrix diagonals
    uppera(1:a-1)= -1._dp/dyc(1:a-1)*(0.5_dp*(edF(1:a-1,i)+edF(2:a,i))*depthF(2:a,i)/dy(1:a-1)&
    +.5_dp*edF(2:a,i)*cbedS(2:a,i)*slopesS(2:a,i))
    lowera(2:a)= -1._dp/dyc(2:a)*(0.5_dp*(edF(1:a-1,i)+edF(2:a,i))*depthF(1:a-1,i)/dy(1:a-1)& 
    -.5_dp*edF(1:a-1,i)*cbedS(1:a-1,i)*slopesS(1:a-1,i))
    diaga(2:a-1)= depthF(2:a-1,i)/DT*2._dp +&
    1._dp/dyc(2:a-1)*(0.5_dp*( edF(3:a,i)+edF(2:a-1,i) )*depthF(2:a-1, i)/dy(2:a-1)&
    -.5_dp*edF(2:a-1,i)*cbedS(2:a-1,i)*slopesS(2:a-1,i))+&
    1._dp/dyc(1:a-2)*(0.5_dp*(edF(1:a-2,i)+edF(2:a-1,i))*depthF(2:a-1, i)/dy(1:a-2)&
    +.5_dp*edF(2:a-1,i)*cbedS(2:a-1,i)*slopesS(2:a-1,i)) + &
    wset*cbedS(2:a-1,i) !!

    diaga(1)= depthF(1,i)/DT*2._dp + &
    1._dp/dyc(1)*(0.5_dp*(edF(2,i)+edF(1,i))*depthF(1,i)/dy(1)- 0.5_dp*edF(1,i)*cbedS(1,i)*slopesS(1,i))
    diaga(a)= depthF(a,i)/DT*2._dp +&
    1._dp/dyc(a)*(0.5_dp*(edF(a,i)+edF(a-1,i))*depthF(a,i)/dy(a-1) + 0.5_dp*edF(a,i)*cbedS(a,i)*slopesS(a,i))

    !Right hand side, upwind method
    do j=1,a
    if(veldh(j,i)>0._dp) THEN
    rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
    (veldh(j,i)*sus_h(j,i) - veldh(j,i-1)*sus_h(j,i-1))/delX !!
    ELSE
    rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
    (veldh(j,i+1)*sus_h(j,i+1) - veldh(j,i)*sus_h(j,i))/delX !!
    END IF
    END DO
    !For boundary conditions here, it is appropriate that the sus is zero on the
    !lateral boundaries-  which equates to us not having to enforce anything
    !further. Unless we have periodic boundaries, in which case the following might
    !be a practical alternative
    !if(periodicb) THEN
    !rhsa(1)=rhsa(1) -lowera(1)*sus_h(1,b)
    !rhsa(a)= rhsa(a)-uppera(a)*sus_h(a,b)
    !END IF
     

    !!Prevent main diagonal from being zero
    do ii=1,a
    if(diaga(ii)==0._dp) THEN
            if( (ii==1).or.(ii==a).or.(depthF(ii,i)==0._dp)) THEN
            diaga(ii)=1._dp
            lowera(ii)=0._dp
            uppera(ii)=0._dp
            rhsa(ii)=0._dp
            else
            diaga(ii)=1._dp
            lowera(ii)=-0.5_dp
            uppera(ii)=-0.5_dp
            rhsa(ii)=0._dp
            end if
    end if
    end do

    call DGTSV(a,1, lowera(2:a), diaga, uppera(1:a-1), rhsa,a, info)

    Cdist(:,i)= rhsa

    if(info.ne.0) print*, "matrix problem in susconc2d, second bit", info

    if(minval(Cdist(:,i))<0._dp) THEN
            Cdist=max(Cdist,0._dp)
        !	if(abs(minval(Cdist(:,i)))< 1.00E-12) THEN

    !	ElSE
        !print*, "Cdist negative with non-roundoff type magnitude"
        !print*, "min Cdist", i, "<0", minval(Cdist(:,i)), minloc(Cdist(:,i)), depthF(minloc(Cdist(:,i)), (i-1):(i+1))& 
        !, vels(minloc(Cdist(:,i)), (i-1):(i+1)) 
        !END IF

    end if


    Cdist_old=cbedS(1:a,:) !At the moment I am just using this to check the output. 

    end do


    !print*, '2nd', Cdist(1,110), Cdist(350,110)

    !stop

end subroutine susconc2dup3

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!Numerical method is upwind in the x derivative advective term, but still ADI
!otherwise
subroutine susconc2dup4(a,b, DT,delX, lengths,elevs,elevs_old,vels,vels_old,hs, Cdist,Cdist_old, Qe,Qe_old, Cmouth, Criver, & 
wset,fs, slopes, l, u)
    !use matrix_solvers
    !use crosssection

    !implicit none
    INTEGER, intent(in):: a, b, l, u
    REAL(dp), intent(in):: DT,delX, wset
    REAL(dp), intent(in):: elevs, elevs_old, vels, vels_old,hs, Qe,Qe_old, Cmouth, Criver,fs, lengths, slopes
    REAL(dp), intent(in out):: Cdist, Cdist_old !!Cdist_old will output random interesting information
    DIMENSION elevs(b), elevs_old(b), vels(a,b), vels_old(a,b), hs(a,b), Cdist(a,b), Cdist_old(a,b), Qe(a,b),fs(a,b), lengths(a,b) &
    ,Criver(a), Cmouth(a), Qe_old(a,b), slopes(a,b), l(b), u(b)

    REAL(dp):: depthN(a,b), depthH(a,b), depthF(a,b), elevsH, velh(a,b),sus_h(a,0:b+1), edN(0:a+1,b), eddN(a,b), veldh(a,0:b+1)
    REAL(dp):: sus_N(0:a+1,b), lnths(0:a+1,b), eddN_yh(0:a,b) ,edF(0:a+1,b), eddF(a,b), eddF_yh(0:a,b)
    REAL(dp):: upperb(b), diagb(b),lowerb(b), rhsb(b), depthN2(0:a+1,b)
    REAL(dp):: uppera(a), diaga(a),lowera(a), rhsa(a), dyc(a), dy(0:a), sus_F(a,b), cbedS(0:a+1,b), slopesS(0:a+1,b)
    REAL(dp):: ded=.20_dp !Dimensionless eddy diffusivity
    REAL(dp):: tmpa(a), tmpb(b),tmp
    INTEGER:: i, j, info, ii
    logical:: periodicb=.false.!false.

    !!!!So this solves the equation
    ! d/dt(depth*C) + d/dx(vel*depth*C)= d/dy(eddn*depth*dC/dy + eddn*Cbed*dh/dy)
    ! where I plan to add the last term soon.
    !By the ADI method, except that the convective x derivative is treated as upwind
    !(seems that it has to be to get stable results). 
    !!Note that assuming a constant vertical eddy diffusivity, 
    !Cbed= Caverage*depth/((ez/ws)*(1-exp(-vs/ez*(depth) ) ))
    !!And I have been assuming that ez=0.5_dp*ey


    !!!First define water depths
    do i= 1,b
    !Depth at previous time step
    depthN(:,i)= max(elevs_old(i)-hs(:,i),0._dp)
    !Depth at half time step
    elevsH=(elevs(i)+elevs_old(i))*.5_dp
    depthH(:,i)= max(elevsH-hs(:,i), 0._dp)
    !Depth at next time step
    depthF(:,i)= max(elevs(i)-hs(:,i),0._dp)
    end do

    !!Useful thing
    depthN2(1:a,:)=depthN
    depthN2(0,:)=0._dp
    depthN2(a+1,:)=0._dp
    !!


    !Now calculate half time step velocity
    velH=.5_dp*(vels+vels_old)
    !Halfway depth times velocity - a useful efficiency device
    veldH(:,1:b)=depthH*velH
    veldH(:,0)= veldH(:,1) !B
    veldH(:,b+1)=veldH(:,b) !B

    !Eddy diffusivity
    edN(1:a,:)= ded*abs(vels_old)*sqrt(fs)*sqrt(0.125_dp)*depthN

    if(periodicb.eqv..false.) THEN
    edN(0,:)=0._dp
    edN(a+1,:)=0._dp

    !Make sure things cancel at the boundaries appropriately
    !do ii=1,b
    !edN(l(ii)-1, ii)=-edN(l(ii),ii)
    !edN(u(ii)+1,ii)=-edN(u(ii),ii)
    !end do
    !edN(0,:)=-edN(1,:)
    !edN(a+1,:)=-edN(a,:)
    ELSE
    edN(0,:)=edN(a,:)
    edN(a+1,:)=edN(1,:)
    end if

    !!Eddy diffusivity*depth at last time step
    eddN= edN(1:a,1:b)*depthN !Note that the sqrt(0.125_dp) is sqrt(1/8)



    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            
            if(edN(j,i)>0._dp) THEN
           tmp= ((.5_dp*edN(j,i))*(1._dp-exp(-wset/(.5_dp*edN(j,i))*depthN(j,i)) ) )
            ELSE
            tmp=0._dp
            END IF

        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthN(j,i)*wset/tmp, 200._dp)
        else
            cbedS(j,i)=200._dp  !So this is the maximum scaling factor we offer
            if((wset==0._dp).and.(edN(j,i)> 0._dp)) cbedS(j,i)=1._dp  !So the idea is that this applies  when wset=0
        
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp


    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times the depth, with boundaries
    eddN_yh(1:a-1,:)= .5_dp*(eddN(2:a,:)+eddN(1:a-1,:))

    if(periodicb.eqv..false.) THEN
    eddN_yh(0,:)=.5_dp*eddN(1,:) !B
    eddN_yh(a,:)=.5_dp*eddN(a,:) !B
    else
    eddN_yh(0,:)=.5_dp*(eddN(1,:)+eddN(a,:)) !B
    eddN_yh(a,:)=.5_dp*(eddN(a,:)+eddN(1,:)) !B
    end if



    !!Susconc at previous time step, with suitable boundary conditions attached
    sus_N(1:a,:)= Cdist_old
    !Zero susconc along the lateral boundaries -- might not be true - zero
    !integrated concentration, but the concentration might not drop off entirely
    sus_N(0,:)=sus_N(1,:) !0._dp !B
    sus_N(a+1,:)=sus_N(a,:) !0._dp !B

    !print*, '0th', sus_N(1,110), sus_N(a,110)


    !Lengths, with suitable boundaries attached
    lnths(1:a,:)= lengths
    lnths(0,:)= 2._dp*lengths(1,:)-lengths(2,:)!lengths(1,:)-depthN(1,:)/max(abs(slopes(1,:)), 0.01_dp)!2._dp*lengths(1,:)-lengths(2,:) !B
    lnths(a+1,:)=  2._dp*lengths(a,:)-lengths(a-1,:)!lengths(a,:) + depthN(a,:)/max(abs(slopes(a,:)), 0.01_dp)!2._dp*lengths(a,:)-lengths(a-1,:) !B

    slopesS(1:a,:)=slopes
    slopesS(0,:)=slopes(1,:) !0._dp
    slopesS(a+1,:)=slopes(a,:) !0._dp

    !!!!!!First half time step - implicit in the longitudinal coord - note how the lateral boundary condition is
    !automatically included - need to include the longitudinal boundary conditions
    !on the RHS
    do j=1,a
    !print*, j, size(Qe_old)
    !!Calculate matrix diagonals and RHS - upwind method
    do i=2,b-1
    !if(veldH(j,i)>0._dp) THEN
            upperb(i)= 0.5_dp*veldh(j,i+1)/delX !0._dp !veldH(j,2:b+1)*.5_dp/delX
            diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) !+ veldh(j,i)/delX !!
            lowerb(i)= -0.5_dp*veldh(j,i-1)/delX !veldh(j,0:b-1)*.5_dp/delX
    !else
    !        upperb(i)= veldh(j,i+1)/delX !veldh(j,2:b+1)*.5_dp/delX
    !        diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) - veldh(j,i)/delX !!
    !        lowerb(i)= 0._dp !veldh(j,0:b-1)*.5_dp/delX
    !end if
    end do
            upperb(1)= 0.5_dp*veldh(j,1+1)/delX !0._dp !veldH(j,2:b+1)*.5_dp/delX
            diagb(1)= depthH(j,1)*2._dp/DT +wset*cbedS(j,1) !+ veldh(j,i)/delX !!
            lowerb(1)= 0._dp !-0.5_dp*veldh(j,i-1)/delX !veldh(j,0:b-1)*.5_dp/delX
            upperb(b)= 0._dp!0.5_dp*veldh(j,i+1)/delX !0._dp !veldH(j,2:b+1)*.5_dp/delX
            diagb(b)= depthH(j,b)*2._dp/DT +wset*cbedS(j,b) !+ veldh(j,i)/delX !!
            lowerb(b)= -0.5_dp*veldh(j,b-1)/delX !veldh(j,0:b-1)*.5_dp/delX
    !print*, j
    rhsb= Qe_old(j,:)*sqrt(1._dp+slopes(j,:)**2)+ depthN(j,:)*sus_N(j,:)*2._dp/DT + 2._dp/(lnths(j+1,:)-lnths(j-1,:))*&  !!
    ( 0.5_dp*(edN(j,:)+edN(j+1,:))*(depthN2(j+1,:)*sus_N(j+1,:)-depthN2(j,:)*sus_N(j,:) )/(lnths(j+1,:)-lnths(j,:)) - & 
    0.5_dp*(edN(j,:)+edN(j-1,:))*(depthN2(j,:)*sus_N(j,:)-depthN2(j-1,:)*sus_N(j-1,:) )/(lnths(j,:)-lnths(j-1,:))  +&
    ((edN(j+1,:)+edN(j,:))*.5_dp*(cbedS(j+1,:)*sus_N(j+1,:)+cbedS(j,:)*sus_N(j,:))*.5_dp*(slopesS(j+1,:)+slopesS(j,:))*.5_dp - &
     (edN(j,:)+edN(j-1,:))*.5_dp*(cbedS(j,:)*sus_N(j,:)+cbedS(j-1,:)*sus_N(j-1,:))*.5_dp*(slopesS(j,:)+slopesS(j-1,:))*.5_dp) &
    !(.5_dp*(edN(j+1,:)*cbedS(j+1,:)*sus_N(j+1,:)*slopesS(j+1,:)+edN(j,:)*cbedS(j,:)*sus_N(j,:)*slopesS(j,:)) - &
    ! (.5_dp*(edN(j,:)*cbedS(j,:)*sus_N(j,:)*slopesS(j,:)+edN(j-1,:)*cbedS(j-1,:)*sus_N(j-1,:))*slopesS(j-1,:))) &
    )

    !(eddN_yh(j,:)*cbedS(j,:) ) )
    !Boundary conditions at the mouth
    if(veldh(j,1)>0._dp) THEN
            rhsb(1)= rhsb(1)+ 0.5_dp*veldh(j,0)/delX*Cmouth(j) !B - this is when the sediment inflows from the 'ocean' 
    else
        upperb(1)= upperb(1)-0.5_dp*veldh(j,1+1)/delX !B - this is when the flow is outward directed, and we assume a zero gradient in veldH*susconc at the mouth - perhaps this is foolish, but lets try it
    END IF
        
    !if(veldh(j,b)<0._dp) THEN
    rhsb(b)= rhsb(b) - 0.5_dp*veldh(j,b+1)/delX*Criver(j) !B - this is where there is river inflow
    !end if

    tmp=rhsb(b)
    !if((j==1).or.(j==350)) print*, lnths(j+1,110)-lnths(j,110),lnths(j,110)-lnths(j-1,110) ! eddN_yh(j,110)*(sus_N(j+1,110)-sus_N(j,110) )/(lnths(j+1,110)-lnths(j,110)) - & 
    !eddN_yh(j-1,110)*((sus_N(j,110)-sus_N(j-1,110) )/(lnths(j,110)-lnths(j-1,110)) ) 
    ! ((edN(j+1,110)+edN(j,110))*.5_dp*(cbedS(j+1,110)*sus_N(j+1,110)+cbedS(j,110)*sus_N(j,110))& 
    !*.5_dp*(slopesS(j+1,110)+slopesS(j,110))*.5_dp - &
    ! (edN(j,110)+edN(j-1,110))*.5_dp*(cbedS(j,110)*sus_N(j,110)+cbedS(j-1,110)*sus_N(j-1,110))& 
    !*.5_dp*(slopesS(j,110)+slopesS(j-1,110))*.5_dp) !sus_N(j+1,110)-sus_N(j,110), sus_N(j,110)-sus_N(j-1,110)

    !!Prevent main diagonal from being zero
    do ii=1,b
    if(diagb(ii)==0._dp) THEN
            if( (ii==1).or.(ii==b).or.(depthH(j,ii)==0._dp)) THEN
            diagb(ii)=1._dp
            lowerb(ii)=0._dp
            upperb(ii)=0._dp
            rhsb(ii)=0._dp
            else
            diagb(ii)=1._dp
            lowerb(ii)=-0.5_dp
            upperb(ii)=-0.5_dp
            rhsb(ii)=0._dp
            end if
    end if
    end do

    call DGTSV(b,1, lowerb(2:b), diagb, upperb(1:b-1), rhsb,b, info)

    sus_h(j,1:b)= rhsb


    if(info.ne.0) print*, "matrix problem in susconc2d, first bit", info

            do ii=1,b
            if(isnan(sus_h(j,ii))) THEN
            print*, 'sus_h(',j,ii,') is NAN' 
            stop
            end if

            if(depthH(j,ii)==0) THEN
                    if(sus_h(j,ii)>0._dp) THEN
                           if(sus_h(j,ii)>0.01_dp) THEN
                            print*, 'zero depth, high sus_h', sus_h(j,ii), j, ii,& 
     depthH(j,ii)*2._dp/DT +wset*cbedS(j,ii) + veldh(j,ii)/delX
                            stop
                            end if

                            sus_h(j,ii)=0._dp
                   end if
            end if
            end do

            if((sus_h(j,b).ne.Criver(j)).and.(depthH(j,b)>0._dp)) THEN
    !               print*, 'sus_h(',j,"b) not equal to Criver", sus_h(j,b),Criver(j), depthH(j,b), depthN(j,b)&
    !,depthH(j,b)*2._dp/DT +wset*cbedS(j,b) - veldh(j,b)/delX,veldh(j,b)/delX, tmp, sus_N(j,b)
            end if
    end do
    !stop
    !print*, sus_h(1,110), sus_h(350,110)

    !print*, "half", sus_h(175,1), sus_h(175,40)
     
    !!!!!!!!!!!!!!End of first half time step

    !!!These boundary conditions are needed for the next half step
    do ii=1, j
    if(veldH(ii,1)>0._dp) THEN
    sus_h(:,0)= Cmouth !aB
    ELSE
    sus_h(:,0)=sus_h(:,1)
    end if
    sus_h(:,b+1)=Criver !B
    end do
    !!!Begin second half time step

    !!Eddy diffusivity
    edF(1:a,:)=ded*abs(vels)*sqrt(fs)*sqrt(0.125_dp)*depthF

    if(periodicb.eqv..false.) THEN
    edF(0,:)=0._dp
    edF(a+1,:)=0._dp
    !Make things cancel at the boundary
    !do ii=1,b
    !edF(l(ii)-1,ii)=-edF(l(ii),ii)
    !edF(u(ii)+1,ii)=-edF(u(ii),ii)
    !end do

    !edF(0,:)=-edF(1,:) !0._dp
    !edF(a+1,:)=-edF(a,:)!0._dp
    ELSE
    edF(0,:)=edF(1,:)
    edF(a+1,:)=edF(a,:)
    END IF
    !!Eddy diffusivity*depth at last time step
    eddF= edF(1:a,:)*depthF !Note that the sqrt(0.125_dp) is sqrt(1/8) - note also how we ignore the change in the friction factor here

    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !cbedS= depthN*wset/((eddN/depthN)*(1._dp-exp(-wset/(eddN/depthN)*depthN) ) )
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            if(wset>0._dp) THEN
            tmp = ((.5_dp*edF(j,i))*(1._dp-exp(-wset/(.5_dp*edF(j,i))*depthF(j,i) ) ))
            else
            tmp=0._dp
            end if
        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthF(j,i)*wset/tmp, 200._dp)
        else
        cbedS(j,i)=200._dp
        if((wset==0).and.(edF(j,i)> 0._dp)) cbedS(j,i)=1._dp
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp


    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times depth, with boundaries
    eddF_yh(1:a-1,:)= .5_dp*(eddF(2:a,:)+eddF(1:a-1,:))
    if(periodicb.eqv..false.) THEN
    eddF_yh(0,:)=.5_dp*eddF(1,:) !B
    eddF_yh(a,:)=.5_dp*eddF(a,:) !B
    else
    eddF_yh(0,:)=.5_dp*(eddF(1,:)+eddF(a,:)) !B
    eddF_yh(a,:)=.5_dp*(eddF(a,:)+eddF(1,:)) !B
    end if


    do i=1,b

    !Define some dy type variables
    dyc(2:a-1)= .5_dp*(lengths(3:a,i)-lengths(1:a-2,i))
    dyc(1)=dyc(2) !0.5_dp*(lengths(2,i)-(lengths(1,i)-depthF(1,i)/max(abs(slopes(1,i)), 0.01_dp)))!dyc(2) !B
    dyc(a)=dyc(a-1)!0.5_dp*((lengths(a,i)+depthF(a,i)/max(abs(slopes(a,i)), 0.01_dp))-lengths(a-1,i))!dyc(a-1) !B

    dy(1:a-1)=lengths(2:a,i)-lengths(1:a-1,i)
    dy(a)=dy(a-1) !depthF(a,i)/max(abs(slopes(a,i)), 0.01_dp) !dy(a-1)
    dy(0)=dy(1) !depthF(1,i)/max(abs(slopes(1,i)), 0.01_dp) !dy(1)

    !!Matrix diagonals
    uppera(1:a-1)= -1._dp/dyc(1:a-1)*(0.5_dp*(edF(1:a-1,i)+edF(2:a,i))*depthF(2:a,i)/dy(1:a-1)&
    +.5_dp*edF(2:a,i)*cbedS(2:a,i)*slopesS(2:a,i))
    lowera(2:a)= -1._dp/dyc(2:a)*(0.5_dp*(edF(1:a-1,i)+edF(2:a,i))*depthF(1:a-1,i)/dy(1:a-1)& 
    -.5_dp*edF(1:a-1,i)*cbedS(1:a-1,i)*slopesS(1:a-1,i))
    diaga(2:a-1)= depthF(2:a-1,i)/DT*2._dp +&
    1._dp/dyc(2:a-1)*(0.5_dp*( edF(3:a,i)+edF(2:a-1,i) )*depthF(2:a-1, i)/dy(2:a-1)&
    -.5_dp*edF(2:a-1,i)*cbedS(2:a-1,i)*slopesS(2:a-1,i))+&
    1._dp/dyc(1:a-2)*(0.5_dp*(edF(1:a-2,i)+edF(2:a-1,i))*depthF(2:a-1, i)/dy(1:a-2)&
    +.5_dp*edF(2:a-1,i)*cbedS(2:a-1,i)*slopesS(2:a-1,i)) + &
    wset*cbedS(2:a-1,i) !!

    diaga(1)= depthF(1,i)/DT*2._dp + &
    1._dp/dyc(1)*(0.5_dp*(edF(2,i)+edF(1,i))*depthF(1,i)/dy(1)- 0.5_dp*edF(1,i)*cbedS(1,i)*slopesS(1,i))
    diaga(a)= depthF(a,i)/DT*2._dp +&
    1._dp/dyc(a)*(0.5_dp*(edF(a,i)+edF(a-1,i))*depthF(a,i)/dy(a-1) + 0.5_dp*edF(a,i)*cbedS(a,i)*slopesS(a,i))

    !Right hand side, upwind method
    do j=1,a
    !if(veldh(j,i)>0._dp) THEN
    rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
    0.5_dp*(veldh(j,i+1)*sus_h(j,i+1) - veldh(j,i-1)*sus_h(j,i-1))/delX !!
    !ELSE
    !rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
    !(veldh(j,i+1)*sus_h(j,i+1) - veldh(j,i)*sus_h(j,i))/delX !!
    !END IF
    END DO
    !For boundary conditions here, it is appropriate that the sus is zero on the
    !lateral boundaries-  which equates to us not having to enforce anything
    !further. Unless we have periodic boundaries, in which case the following might
    !be a practical alternative
    !if(periodicb) THEN
    !rhsa(1)=rhsa(1) -lowera(1)*sus_h(1,b)
    !rhsa(a)= rhsa(a)-uppera(a)*sus_h(a,b)
    !END IF
     

    !!Prevent main diagonal from being zero
    do ii=1,a
    if(diaga(ii)==0._dp) THEN
            if( (ii==1).or.(ii==a).or.(depthF(ii,i)==0._dp)) THEN
            diaga(ii)=1._dp
            lowera(ii)=0._dp
            uppera(ii)=0._dp
            rhsa(ii)=0._dp
            else
            diaga(ii)=1._dp
            lowera(ii)=-0.5_dp
            uppera(ii)=-0.5_dp
            rhsa(ii)=0._dp
            end if
    end if
    end do

    call DGTSV(a,1, lowera(2:a), diaga, uppera(1:a-1), rhsa,a, info)

    Cdist(:,i)= rhsa

    if(info.ne.0) print*, "matrix problem in susconc2d, second bit", info

    if(minval(Cdist(:,i))<0._dp) THEN
            Cdist=max(Cdist,0._dp)
        !	if(abs(minval(Cdist(:,i)))< 1.00E-12) THEN

    !	ElSE
        !print*, "Cdist negative with non-roundoff type magnitude"
        !print*, "min Cdist", i, "<0", minval(Cdist(:,i)), minloc(Cdist(:,i)), depthF(minloc(Cdist(:,i)), (i-1):(i+1))& 
        !, vels(minloc(Cdist(:,i)), (i-1):(i+1)) 
        !END IF

    end if


    Cdist_old=cbedS(1:a,:) !At the moment I am just using this to check the output. 

    end do


    !print*, '2nd', Cdist(1,110), Cdist(350,110)

    !stop

end subroutine susconc2dup4
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!Numerical method is upwind in the x derivative advective term, but still ADI
!otherwise
subroutine susconc2dup3V(a,b, DT,delX, lengths,elevs,elevs_old,vels,vels_old,hs, Cdist,Cdist_old, Qe,Qe_old, Cmouth, Criver, & 
wset,fs, slopes, l, u)
    !use matrix_solvers
    !use crosssection

    !implicit none
    INTEGER, intent(in):: a, b, l, u
    REAL(dp), intent(in):: DT,delX, wset
    REAL(dp), intent(in):: elevs, elevs_old, vels, vels_old,hs, Qe,Qe_old, Cmouth, Criver,fs, lengths, slopes
    REAL(dp), intent(in out):: Cdist, Cdist_old !!Cdist_old will output random interesting information
    DIMENSION elevs(b), elevs_old(b), vels(a,b), vels_old(a,b), hs(a,b), Cdist(a,b), Cdist_old(a,b), Qe(a,b),fs(a,b), lengths(a,b) &
    ,Criver(a), Cmouth(a), Qe_old(a,b), slopes(a,b), l(b), u(b)

    REAL(dp):: depthN(a,b), depthH(a,b), depthF(a,b), elevsH, velh(a,b),sus_h(a,0:b+1), edN(0:a+1,b), eddN(a,b), veldh(a,0:b+1)
    REAL(dp):: sus_N(0:a+1,b), lnths(0:a+1,b), eddN_yh(0:a,b) ,edF(0:a+1,b), eddF(a,b), eddF_yh(0:a,b)
    REAL(dp):: upperb(b), diagb(b),lowerb(b), rhsb(b), depthN2(0:a+1,b), VDN(0:a+1,b), VDF(a,b)
    REAL(dp):: uppera(a), diaga(a),lowera(a), rhsa(a), dyc(a), dy(0:a), sus_F(a,b), cbedS(0:a+1,b), slopesS(0:a+1,b)
    REAL(dp):: dedy=0.20_dp, dedz=0.10_dp !Dimensionless eddy diffusivity
    REAL(dp):: tmpa(a), tmpb(b),tmp
    INTEGER:: i, j, info, ii
    logical:: periodicb=.false.!false.

    !!!!So this solves the equation
    ! d/dt(depth*C) + d/dx(vel*depth*C)= d/dy(eddn*depth*dC/dy + eddn*Cbed*dh/dy)
    ! where I plan to add the last term soon.
    !By the ADI method, except that the convective x derivative is treated as upwind
    !(seems that it has to be to get stable results). 
    !!Note that assuming a constant vertical eddy diffusivity, 
    !Cbed= Caverage*depth/((ez/ws)*(1-exp(-vs/ez*(depth) ) ))
    !!And I have been assuming that ez=0.5_dp*ey


    !!!First define water depths
    do i= 1,b
    !Depth at previous time step
    depthN(:,i)= max(elevs_old(i)-hs(:,i),0._dp)
    !Depth at half time step
    elevsH=(elevs(i)+elevs_old(i))*.5_dp
    depthH(:,i)= max(elevsH-hs(:,i), 0._dp)
    !Depth at next time step
    depthF(:,i)= max(elevs(i)-hs(:,i),0._dp)
    end do

    !!Useful thing
    depthN2(1:a,:)=depthN
    depthN2(0,:)=0._dp
    depthN2(a+1,:)=0._dp
    !!


    !Now calculate half time step velocity
    velH=.5_dp*(vels+vels_old)
    !Halfway depth times velocity - a useful efficiency device
    veldH(:,1:b)=depthH*velH
    veldH(:,0)= veldH(:,1) !B
    veldH(:,b+1)=veldH(:,b) !B

    !Eddy diffusivity
    edN(1:a,:)= dedy*abs(vels_old)*sqrt(fs)*sqrt(0.125_dp)*depthN

    if(periodicb.eqv..false.) THEN
    edN(0,:)=0._dp
    edN(a+1,:)=0._dp

    !Make sure things cancel at the boundaries appropriately
    !do ii=1,b
    !edN(l(ii)-1, ii)=-edN(l(ii),ii)
    !edN(u(ii)+1,ii)=-edN(u(ii),ii)
    !end do
    !edN(0,:)=-edN(1,:)
    !edN(a+1,:)=-edN(a,:)
    ELSE
    edN(0,:)=edN(a,:)
    edN(a+1,:)=edN(1,:)
    end if

    !!Eddy diffusivity*depth at last time step
    eddN= edN(1:a,1:b)*depthN !Note that the sqrt(0.125_dp) is sqrt(1/8)



    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            
            if(edN(j,i)>0._dp) THEN
           tmp= ((edN(j,i)/dedy*dedz)*(1._dp-exp(-wset/(edN(j,i)/dedy*dedz)*depthN(j,i)) ) )
            ELSE
            tmp=0._dp
            END IF

        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthN(j,i)*wset/tmp, 200._dp)
        else
            cbedS(j,i)=200._dp  !So this is the maximum scaling factor we offer
            if((wset==0._dp).and.(edN(j,i)> 0._dp)) cbedS(j,i)=1._dp  !So the idea is that this applies  when wset=0
        
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp


    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times the depth, with boundaries
    eddN_yh(1:a-1,:)= .5_dp*(eddN(2:a,:)+eddN(1:a-1,:))

    if(periodicb.eqv..false.) THEN
    eddN_yh(0,:)=.5_dp*eddN(1,:) !B
    eddN_yh(a,:)=.5_dp*eddN(a,:) !B
    else
    eddN_yh(0,:)=.5_dp*(eddN(1,:)+eddN(a,:)) !B
    eddN_yh(a,:)=.5_dp*(eddN(a,:)+eddN(1,:)) !B
    end if



    !!Susconc at previous time step, with suitable boundary conditions attached
    sus_N(1:a,:)= Cdist_old
    !Zero susconc along the lateral boundaries -- might not be true - zero
    !integrated concentration, but the concentration might not drop off entirely
    sus_N(0,:)=sus_N(1,:) !0._dp !B
    sus_N(a+1,:)=sus_N(a,:) !0._dp !B

    !print*, '0th', sus_N(1,110), sus_N(a,110)


    !Lengths, with suitable boundaries attached
    lnths(1:a,:)= lengths
    lnths(0,:)= 2._dp*lengths(1,:)-lengths(2,:)!lengths(1,:)-depthN(1,:)/max(abs(slopes(1,:)), 0.01_dp)!2._dp*lengths(1,:)-lengths(2,:) !B
    lnths(a+1,:)=  2._dp*lengths(a,:)-lengths(a-1,:)!lengths(a,:) + depthN(a,:)/max(abs(slopes(a,:)), 0.01_dp)!2._dp*lengths(a,:)-lengths(a-1,:) !B

    slopesS(1:a,:)=slopes
    slopesS(0,:)=slopes(1,:) !0._dp
    slopesS(a+1,:)=slopes(a,:) !0._dp

    !Calculate VDN, which is the lateral discharge V*D at each point.
    !Based on the idea that
    ! VD = int( -dY/dt - dUD/dx ) dy
    !To ensure conservation, we calculate dY/dt to ensure that over the whole
    !cross-section, int(VD)=0
    VDN=0._dp 
    do i=2, b-1
    !Calculate the integral of dUd/dx. Note that we use u(i)+1 as a trick
    do j=max(l(i),2), min(u(i)+1,a)
    VDN(j,i) = VDN(j-1,i) + 0.5_dp*(-vels_old(j,i+1)*max(elevs_old(i+1)-hs(j,i+1),0._dp) - &
    vels_old(j-1,i+1)*max(elevs_old(i+1)-hs(j-1,i+1),0._dp) + & 
    vels_old(j,i-1)*max(elevs_old(i-1)-hs(j,i-1), 0._dp) + &
    vels_old(j-1,i-1)*max(elevs_old(i-1)-hs(j-1,i-1), 0._dp) )& 
    /(2._dp*delX)*(lengths(j,i)-lengths(j-1,i))
    end do
    !Now calculate dY/dT and correct - note the use of VDN(u(i)+1,i), which is later
    !set to zero
    VDN(l(i):u(i),i) = VDN(l(i):u(i),i) - &
    (lengths(l(i):u(i),i)-(lengths(max(l(i)-1,1),i) ))*&
    (VDN(min(u(i)+1,a),i)/((lengths(min(u(i)+1,a),i))-(lengths(max(l(i)-1,1),i) ) )) !Note that the last (VDN(u(i).. ) is dY/dt 

    VDN(min(u(i)+1,a),i)=0._dp
    !VDN(1:l(i)-1,i) = 999._dp
    !VDN(u(i)+1:a,i)=999._dp

    end do
    !stop
    !print*, 'BUG'
    !VDN=0._dp
    !write(100,*) VDN(:,20) !, (-vels_old(:,21)*max(elevs_old(21)-hs(:,21),0._dp) + & 
    !vels_old(:,19)*max(elevs_old(19)-hs(:,19), 0._dp))

    !!!!!!First half time step - implicit in the longitudinal coord - note how the lateral boundary condition is
    !automatically included - need to include the longitudinal boundary conditions
    !on the RHS

    do j=1,a
    !print*, j, size(Qe_old)
    !!Calculate matrix diagonals and RHS - upwind method
    do i=1,b
    if(veldH(j,i)>0._dp) THEN
            upperb(i)= 0._dp !veldH(j,2:b+1)*.5_dp/delX
            diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) + veldh(j,i)/delX !!
            lowerb(i)= -veldh(j,i-1)/delX !veldh(j,0:b-1)*.5_dp/delX
    else
            upperb(i)= veldh(j,i+1)/delX !veldh(j,2:b+1)*.5_dp/delX
            diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) - veldh(j,i)/delX !!
            lowerb(i)= 0._dp !veldh(j,0:b-1)*.5_dp/delX
    end if
    end do
    !print*, j
    rhsb= Qe_old(j,:)*sqrt(1._dp+slopes(j,:)**2)+ depthN(j,:)*sus_N(j,:)*2._dp/DT + 2._dp/(lnths(j+1,:)-lnths(j-1,:))*&  !!
    ( 0.5_dp*(edN(j,:)+edN(j+1,:))*(depthN2(j+1,:)*sus_N(j+1,:)-depthN2(j,:)*sus_N(j,:) )/(lnths(j+1,:)-lnths(j,:)) - & 
    0.5_dp*(edN(j,:)+edN(j-1,:))*(depthN2(j,:)*sus_N(j,:)-depthN2(j-1,:)*sus_N(j-1,:) )/(lnths(j,:)-lnths(j-1,:))  +&
    ((edN(j+1,:)+edN(j,:))*.5_dp*(cbedS(j+1,:)*sus_N(j+1,:)+cbedS(j,:)*sus_N(j,:))*.5_dp*(slopesS(j+1,:)+slopesS(j,:))*.5_dp - &
     (edN(j,:)+edN(j-1,:))*.5_dp*(cbedS(j,:)*sus_N(j,:)+cbedS(j-1,:)*sus_N(j-1,:))*.5_dp*(slopesS(j,:)+slopesS(j-1,:))*.5_dp) &
    !(.5_dp*(edN(j+1,:)*cbedS(j+1,:)*sus_N(j+1,:)*slopesS(j+1,:)+edN(j,:)*cbedS(j,:)*sus_N(j,:)*slopesS(j,:)) - &
    ! (.5_dp*(edN(j,:)*cbedS(j,:)*sus_N(j,:)*slopesS(j,:)+edN(j-1,:)*cbedS(j-1,:)*sus_N(j-1,:))*slopesS(j-1,:))) &
    ) - (VDN(j+1,:)*sus_N(j+1,:) - VDN(j-1,:)*sus_N(j-1,:))/(lnths(j+1,:)-lnths(j-1,:))

    !(eddN_yh(j,:)*cbedS(j,:) ) )
    !Boundary conditions at the mouth
    if(veldh(j,1)>0._dp) THEN
            rhsb(1)= rhsb(1)+ veldh(j,0)/delX*Cmouth(j) !B - this is when the sediment inflows from the 'ocean' 
    !else
    !	diagb(1)= diagb(1)-veldh(j,1)/delX !B - this is when the flow is outward directed, and we assume a zero gradient in veldH*susconc at the mouth - perhaps this is foolish, but lets try it
    END IF
        
    if(veldh(j,b)<0._dp) THEN
    rhsb(b)= rhsb(b) - veldh(j,b+1)/delX*Criver(j) !B - this is where there is river inflow
    end if

    tmp=rhsb(b)
    !if((j==1).or.(j==350)) print*, lnths(j+1,110)-lnths(j,110),lnths(j,110)-lnths(j-1,110) ! eddN_yh(j,110)*(sus_N(j+1,110)-sus_N(j,110) )/(lnths(j+1,110)-lnths(j,110)) - & 
    !eddN_yh(j-1,110)*((sus_N(j,110)-sus_N(j-1,110) )/(lnths(j,110)-lnths(j-1,110)) ) 
    ! ((edN(j+1,110)+edN(j,110))*.5_dp*(cbedS(j+1,110)*sus_N(j+1,110)+cbedS(j,110)*sus_N(j,110))& 
    !*.5_dp*(slopesS(j+1,110)+slopesS(j,110))*.5_dp - &
    ! (edN(j,110)+edN(j-1,110))*.5_dp*(cbedS(j,110)*sus_N(j,110)+cbedS(j-1,110)*sus_N(j-1,110))& 
    !*.5_dp*(slopesS(j,110)+slopesS(j-1,110))*.5_dp) !sus_N(j+1,110)-sus_N(j,110), sus_N(j,110)-sus_N(j-1,110)

    !!Prevent main diagonal from being zero
    do ii=1,b
    if((diagb(ii)==0._dp).or.(depthH(j,ii)==0._dp)) THEN
            if( (ii==1).or.(ii==b).or.(depthH(j,ii)==0._dp)) THEN
            diagb(ii)=1._dp
            lowerb(ii)=0._dp
            upperb(ii)=0._dp
            rhsb(ii)=0._dp
            else
            diagb(ii)=1._dp
            lowerb(ii)=-0.5_dp
            upperb(ii)=-0.5_dp
            rhsb(ii)=0._dp
            end if
    end if
    end do

    call DGTSV(b,1, lowerb(2:b), diagb, upperb(1:b-1), rhsb,b, info)

    sus_h(j,1:b)= rhsb


    if(info.ne.0) print*, "matrix problem in susconc2d, first bit", info
            !CHECKS
            do ii=1,b
            if(isnan(sus_h(j,ii))) THEN
            print*, 'sus_h(',j,ii,') is NAN' 
            stop
            end if
            
            !if(sus_h(j,ii)<0._dp) THEN
            !        if(sus_h(j,ii)< - 1E-04_dp) print*, 'significant negative sus_h', j, ii, sus_h(j,ii), sus_h(700,ii)
            !        sus_h(j,ii)=0._dp
            !From limited experimentation, it seems that even if this goes negative,
            !we don't end up with significant negative values in the final Cdist-
            !thus, maybe best to leave it alone.
            !end if
            
            if(depthH(j,ii)==0) THEN
                    if(sus_h(j,ii)>0._dp) THEN
                           if(sus_h(j,ii)>0.01_dp) THEN
                            print*, 'zero depth, high sus_h', sus_h(j,ii), j, ii,& 
     depthH(j,ii)*2._dp/DT +wset*cbedS(j,ii) + veldh(j,ii)/delX, veldh(j,i+1)/delX, veldh(j,i-1)/delX 
                            stop
                            end if

                            sus_h(j,ii)=0._dp
                   end if
            end if
            end do

    !        if((sus_h(j,b).ne.Criver(j)).and.(depthH(j,b)>0._dp)) THEN
    !               print*, 'sus_h(',j,"b) not equal to Criver", sus_h(j,b),Criver(j), depthH(j,b), depthN(j,b)&
    !,depthH(j,b)*2._dp/DT +wset*cbedS(j,b) - veldh(j,b)/delX,veldh(j,b)/delX, tmp, sus_N(j,b)
    !        end if


    end do
    !stop
    !print*, sus_h(1,110), sus_h(350,110)
    !print*, maxval(sus_h), minval(sus_h), maxval(Cdist_old), minval(Cdist_old)
    !print*, "half", sus_h(175,1), sus_h(175,40)
    !open(123,file='sus_h')
    !write(123,*) sus_h
    !close(123) 
    !!!!!!!!!!!!!!End of first half time step


    !do i=1, b
    !if((sus_h(l(i),i)>3._dp*sus_h(l(i)+1,i)).and.(abs(sus_h(l(i),i))>1.0E-10_dp) ) print*, 'sus_h spike', i,l(i), & 
    !elevs_old(i)-hs(l(i)-1:l(i)+1,i),"##", sus_h(l(i):l(i)+1,i), sus_h(700,i)
    !end do


    !!!These boundary conditions are needed for the next half step
    sus_h(:,0)= Cmouth !B
    sus_h(:,b+1)=Criver !B

    !!!Begin second half time step

    !!Eddy diffusivity
    edF(1:a,:)=dedy*abs(vels)*sqrt(fs)*sqrt(0.125_dp)*depthF

    if(periodicb.eqv..false.) THEN
    edF(0,:)=0._dp
    edF(a+1,:)=0._dp
    !Make things cancel at the boundary
    !do ii=1,b
    !edF(l(ii)-1,ii)=-edF(l(ii),ii)
    !edF(u(ii)+1,ii)=-edF(u(ii),ii)
    !end do

    !edF(0,:)=-edF(1,:) !0._dp
    !edF(a+1,:)=-edF(a,:)!0._dp
    ELSE
    edF(0,:)=edF(1,:)
    edF(a+1,:)=edF(a,:)
    END IF
    !!Eddy diffusivity*depth at last time step
    eddF= edF(1:a,:)*depthF !Note that the sqrt(0.125_dp) is sqrt(1/8) - note also how we ignore the change in the friction factor here

    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !cbedS= depthN*wset/((eddN/depthN)*(1._dp-exp(-wset/(eddN/depthN)*depthN) ) )
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            if(wset>0._dp) THEN
            tmp = ((edF(j,i)/dedy*dedz)*(1._dp-exp(-wset/(edF(j,i)/dedy*dedz)*depthF(j,i) ) ))
            else
            tmp=0._dp
            end if
        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthF(j,i)*wset/tmp, 200._dp)
        else
        cbedS(j,i)=200._dp
        if((wset==0).and.(edF(j,i)> 0._dp)) cbedS(j,i)=1._dp
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp

    !print*, cbedS(175,75), cbedS(175,75)*sus_h(175,75)
    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times depth, with boundaries
    eddF_yh(1:a-1,:)= .5_dp*(eddF(2:a,:)+eddF(1:a-1,:))
    if(periodicb.eqv..false.) THEN
    eddF_yh(0,:)=.5_dp*eddF(1,:) !B
    eddF_yh(a,:)=.5_dp*eddF(a,:) !B
    else
    eddF_yh(0,:)=.5_dp*(eddF(1,:)+eddF(a,:)) !B
    eddF_yh(a,:)=.5_dp*(eddF(a,:)+eddF(1,:)) !B
    end if



    !Calculate VDF, which is the lateral discharge V*D at each point.
    !Based on the idea that
    ! VD = int( -dY/dt - dUD/dx ) dy
    !To ensure conservation, we calculate dY/dt to ensure that over the whole
    !cross-section, int(VD)=0
    VDF=0._dp 
    do i=2, b-1
    !Calculate the integral of dUd/dx. Note that we use u(i)+1 as a trick
    do j=max(l(i),2), min(u(i)+1,a)
    VDF(j,i) = VDF(j-1,i) + 0.5_dp*(-vels(j,i+1)*max(elevs(i+1)-hs(j,i+1),0._dp) - &
     vels(j-1,i+1)*max(elevs(i+1)-hs(j-1,i+1),0._dp) + & 
    vels(j,i-1)*max(elevs(i-1)-hs(j,i-1), 0._dp) + &
    vels(j-1,i-1)*max(elevs(i-1)-hs(j-1,i-1), 0._dp) )& 
    /(2._dp*delX)*(lengths(j,i)-lengths(j-1,i))
    end do
    !Now calculate dY/dT and correct - note the use of VDF(u(i)+1,i), which is later
    !set to zero
    VDF(l(i):u(i),i) = VDF(l(i):u(i),i) - &
    (lengths(l(i):u(i),i)-(lengths(max(l(i)-1,1),i)) )*&
    (VDF(min(u(i)+1,a),i)/((lengths(min(u(i)+1,a),i))-(lengths(max(l(i)-1,1),i) ) )) !Note that the last (VDN(u(i).. ) is dY/dt 

    VDF(min(u(i)+1,a),i)=0._dp
    !VDN(1:l(i)-1,i) = 999._dp
    !VDN(u(i)+1:a,i)=999._dp
    end do
    !write(100,*) VDF(:,20)
    !print*, 'BUG'
    !VDF=0._dp
    !open(123,file='VDN')
    !write(123,*) VDF
    !print*, size(VDF), size(VDF(:,1))


    do i=1,b

    !Define some dy type variables
    dyc(2:a-1)= .5_dp*(lengths(3:a,i)-lengths(1:a-2,i))
    dyc(1)=dyc(2) !0.5_dp*(lengths(2,i)-(lengths(1,i)-depthF(1,i)/max(abs(slopes(1,i)), 0.01_dp)))!dyc(2) !B
    dyc(a)=dyc(a-1)!0.5_dp*((lengths(a,i)+depthF(a,i)/max(abs(slopes(a,i)), 0.01_dp))-lengths(a-1,i))!dyc(a-1) !B

    dy(1:a-1)=lengths(2:a,i)-lengths(1:a-1,i)
    dy(a)=dy(a-1) !depthF(a,i)/max(abs(slopes(a,i)), 0.01_dp) !dy(a-1)
    dy(0)=dy(1) !depthF(1,i)/max(abs(slopes(1,i)), 0.01_dp) !dy(1)

    !!Matrix diagonals
    uppera(1:a-1)= -1._dp/dyc(1:a-1)*(0.5_dp*(edF(1:a-1,i)+edF(2:a,i))*depthF(2:a,i)/dy(1:a-1)&
    +.5_dp*edF(2:a,i)*cbedS(2:a,i)*slopesS(2:a,i)) &
    + VDF(2:a,i)/(2._dp*dyc(1:a-1))
    lowera(2:a)= -1._dp/dyc(2:a)*(0.5_dp*(edF(1:a-1,i)+edF(2:a,i))*depthF(1:a-1,i)/dy(1:a-1)& 
    -.5_dp*edF(1:a-1,i)*cbedS(1:a-1,i)*slopesS(1:a-1,i)) &
    -VDF(1:a-1,i)/(2._dp*dyc(2:a))
    diaga(2:a-1)= depthF(2:a-1,i)/DT*2._dp +&
    1._dp/dyc(2:a-1)*(0.5_dp*( edF(3:a,i)+edF(2:a-1,i) )*depthF(2:a-1, i)/dy(2:a-1)&
    -.5_dp*edF(2:a-1,i)*cbedS(2:a-1,i)*slopesS(2:a-1,i))+&
    1._dp/dyc(1:a-2)*(0.5_dp*(edF(1:a-2,i)+edF(2:a-1,i))*depthF(2:a-1, i)/dy(1:a-2)&
    +.5_dp*edF(2:a-1,i)*cbedS(2:a-1,i)*slopesS(2:a-1,i)) + &
    wset*cbedS(2:a-1,i) !!

    diaga(1)= depthF(1,i)/DT*2._dp + &
    1._dp/dyc(1)*(0.5_dp*(edF(2,i)+edF(1,i))*depthF(1,i)/dy(1)- 0.5_dp*edF(1,i)*cbedS(1,i)*slopesS(1,i))
    diaga(a)= depthF(a,i)/DT*2._dp +&
    1._dp/dyc(a)*(0.5_dp*(edF(a,i)+edF(a-1,i))*depthF(a,i)/dy(a-1) + 0.5_dp*edF(a,i)*cbedS(a,i)*slopesS(a,i))

    !Right hand side, upwind method
    do j=1,a
    if(veldh(j,i)>0._dp) THEN
    rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
    (veldh(j,i)*sus_h(j,i) - veldh(j,i-1)*sus_h(j,i-1))/delX !!
    ELSE
    rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
    (veldh(j,i+1)*sus_h(j,i+1) - veldh(j,i)*sus_h(j,i))/delX !!
    END IF
    END DO
    !For boundary conditions here, it is appropriate that the sus is zero on the
    !lateral boundaries-  which equates to us not having to enforce anything
    !further. Unless we have periodic boundaries, in which case the following might
    !be a practical alternative
    !if(periodicb) THEN
    !rhsa(1)=rhsa(1) -lowera(1)*sus_h(1,b)
    !rhsa(a)= rhsa(a)-uppera(a)*sus_h(a,b)
    !END IF
     

    !!Prevent main diagonal from being zero
    do ii=1,a
    if((diaga(ii)==0._dp).or.(depthF(ii,i)==0._dp)) THEN
            if( (ii==1).or.(ii==a).or.(depthF(ii,i)==0._dp)) THEN
            diaga(ii)=1._dp
            lowera(ii)=0._dp
            uppera(ii)=0._dp
            rhsa(ii)=0._dp
            else
            diaga(ii)=1._dp
            lowera(ii)=-0.5_dp
            uppera(ii)=-0.5_dp
            rhsa(ii)=0._dp
            end if
    end if
    end do

    call DGTSV(a,1, lowera(2:a), diaga, uppera(1:a-1), rhsa,a, info)

    Cdist(:,i)= rhsa

    if(info.ne.0) print*, "matrix problem in susconc2d, second bit", info


    do ii=1, a
    if(Cdist(ii,i)<0._dp) THEN
            if(Cdist(ii,i)< - 1E-08_dp) print*, 'significant negative Cdist', j, ii, Cdist(ii,i), Cdist(a/2,i)
            !Cdist(j,ii)=0._dp
    end if
    end do

    if(minval(Cdist(:,i))<0._dp) THEN
            Cdist(:,i)=max(Cdist(:,i),0._dp)
        !	if(abs(minval(Cdist(:,i)))< 1.00E-12) THEN

    !	ElSE
        !print*, "Cdist negative with non-roundoff type magnitude"
        !print*, "min Cdist", i, "<0", minval(Cdist(:,i)), minloc(Cdist(:,i)), depthF(minloc(Cdist(:,i)), (i-1):(i+1))& 
        !, vels(minloc(Cdist(:,i)), (i-1):(i+1)) 
        !END IF

    end if


    end do


    Cdist_old=cbedS(1:a,:) !At the moment I am just using this to check the output. 
    !do i=1, b
    !if((Cdist(l(i),i)>3._dp*Cdist(l(i)+1,i)).and.(abs(Cdist(l(i),i))>1.0E-10_dp ) ) print*, 'Cdist spike end', i,l(i), & 
    !elevs(i)-hs(l(i)-1:l(i)+1,i), '##', Cdist(l(i):l(i)+1,i), Cdist(700,i)  
    !end do

    !print*, '2nd', Cdist(1,110), Cdist(350,110)

    !stop

end subroutine susconc2dup3V
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!Numerical method is upwind in the x derivative advective term, but still ADI
!otherwise
subroutine susconc2dup3rdV(a,b, DT,delX, lengths,elevs,elevs_old,vels,vels_old,hs, Cdist,Cdist_old, Qe,Qe_old, Cmouth, Criver, & 
wset,fs, slopes, l, u)
    !use matrix_solvers
    !use crosssection

    !implicit none
    INTEGER, intent(in):: a, b, l, u
    REAL(dp), intent(in):: DT,delX, wset
    REAL(dp), intent(in):: elevs, elevs_old, vels, vels_old,hs, Qe,Qe_old, Cmouth, Criver,fs, lengths, slopes
    REAL(dp), intent(in out):: Cdist, Cdist_old !!Cdist_old will output random interesting information
    DIMENSION elevs(b), elevs_old(b), vels(a,b), vels_old(a,b), hs(a,b), Cdist(a,b), Cdist_old(a,b), Qe(a,b),fs(a,b),& 
            lengths(a,b), Criver(a), Cmouth(a), Qe_old(a,b), slopes(a,b), l(b), u(b)

    INTEGER:: i, j, info, ii, KL=2, KU=2, IPV((-1):(b+2))
    logical:: periodicb=.false.!false.
    REAL(dp):: depthN(a,b), depthH(a,b), depthF(a,b), elevsH, velh(a,b),sus_h(a,(-1):(b+2)),edN(0:a+1,b),eddN(a,b),& 
            veldh(a,(-1):(b+2))
    REAL(dp):: sus_N(0:a+1,b), lnths(0:a+1,b), eddN_yh(0:a,b) ,edF(0:a+1,b), eddF(a,b), eddF_yh(0:a,b)
    REAL(dp):: upperb(b), diagb(b),lowerb(b), rhsb((-1):(b+2)), depthN2(0:a+1,b), VDN(0:a+1,b), VDF(a,b)
    REAL(dp):: uppera(a), diaga(a),lowera(a), rhsa(a), dyc(a), dy(0:a), sus_F(a,b), cbedS(0:a+1,b), slopesS(0:a+1,b)
    REAL(dp):: dedy=0.20_dp, dedz=0.10_dp !Dimensionless eddy diffusivity
    REAL(dp):: tmpa(a), tmpb(b),tmp, band(7,b+4)

    !!!!So this solves the equation
    ! d/dt(depth*C) + d/dx(vel*depth*C)= d/dy(eddn*depth*dC/dy + eddn*Cbed*dh/dy)
    ! where I plan to add the last term soon.
    !By the ADI method, except that the convective x derivative is treated as upwind
    !(seems that it has to be to get stable results). 
    !!Note that assuming a constant vertical eddy diffusivity, 
    !Cbed= Caverage*depth/((ez/ws)*(1-exp(-vs/ez*(depth) ) ))
    !!And I have been assuming that ez=0.5_dp*ey

    !There are ghost points on the left and right edges (i.e. the margins of each
    !cross-section -- 1 each), and on the north
    !and south edges (i.e. the upstream and downstream boundaries -- 2 each)


    !!!First define water depths
    do i= 1,b
    !Depth at previous time step
    depthN(:,i)= max(elevs_old(i)-hs(:,i),0._dp)
    !Depth at half time step
    elevsH=(elevs(i)+elevs_old(i))*.5_dp
    depthH(:,i)= max(elevsH-hs(:,i), 0._dp)
    !Depth at next time step
    depthF(:,i)= max(elevs(i)-hs(:,i),0._dp)
    end do

    !!Useful thing
    depthN2(1:a,:)=depthN
    depthN2(0,:)=0._dp
    depthN2(a+1,:)=0._dp
    !!


    !Now calculate half time step velocity
    velH=.5_dp*(vels+vels_old)
    !Halfway depth times velocity - a useful efficiency device
    veldH(:,1:b)=depthH*velH
    veldH(:,0)= veldH(:,1) !B
    veldH(:,b+1)=veldH(:,b) !B
    veldH(:,-1)=veldH(:,0)
    veldH(:,b+2)=veldH(:,b) !B

    !Eddy diffusivity

    edN(1:a,:)= dedy*abs(vels_old)*sqrt(fs)*sqrt(0.125_dp)*depthN

    if(periodicb.eqv..false.) THEN
    edN(0,:)=0._dp
    edN(a+1,:)=0._dp

    !Make sure things cancel at the boundaries appropriately
    !do ii=1,b
    !edN(l(ii)-1, ii)=-edN(l(ii),ii)
    !edN(u(ii)+1,ii)=-edN(u(ii),ii)
    !end do
    !edN(0,:)=-edN(1,:)
    !edN(a+1,:)=-edN(a,:)
    ELSE
    edN(0,:)=edN(a,:)
    edN(a+1,:)=edN(1,:)
    end if

    !!Eddy diffusivity*depth at last time step
    eddN= edN(1:a,1:b)*depthN !Note that the sqrt(0.125_dp) is sqrt(1/8)



    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            
            if(edN(j,i)>0._dp) THEN
           tmp= ((edN(j,i)/dedy*dedz)*(1._dp-exp(-wset/(edN(j,i)/dedy*dedz)*depthN(j,i)) ) )
            ELSE
            tmp=0._dp
            END IF

        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthN(j,i)*wset/tmp, 200._dp)
        else
            cbedS(j,i)=200._dp  !So this is the maximum scaling factor we offer
            if((wset==0._dp).and.(edN(j,i)> 0._dp)) cbedS(j,i)=1._dp  !So the idea is that this applies  when wset=0
        
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp

    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times the depth, with boundaries
    eddN_yh(1:a-1,:)= .5_dp*(eddN(2:a,:)+eddN(1:a-1,:))

    if(periodicb.eqv..false.) THEN
    eddN_yh(0,:)=.5_dp*eddN(1,:) !B
    eddN_yh(a,:)=.5_dp*eddN(a,:) !B
    else
    eddN_yh(0,:)=.5_dp*(eddN(1,:)+eddN(a,:)) !B
    eddN_yh(a,:)=.5_dp*(eddN(a,:)+eddN(1,:)) !B
    end if

    !!Susconc at previous time step, with suitable boundary conditions attached
    sus_N(1:a,:)= Cdist_old
    !Zero susconc along the lateral boundaries -- might not be true - zero
    !integrated concentration, but the concentration might not drop off entirely
    sus_N(0,:)=sus_N(1,:) !0._dp !B
    sus_N(a+1,:)=sus_N(a,:) !0._dp !B


    !Lengths, with suitable boundaries attached
    lnths(1:a,:)= lengths
    lnths(0,:)= 2._dp*lengths(1,:)-lengths(2,:)!lengths(1,:)-depthN(1,:)/max(abs(slopes(1,:)), 0.01_dp)!2._dp*lengths(1,:)-lengths(2,:) !B
    lnths(a+1,:)=  2._dp*lengths(a,:)-lengths(a-1,:)!lengths(a,:) + depthN(a,:)/max(abs(slopes(a,:)), 0.01_dp)!2._dp*lengths(a,:)-lengths(a-1,:) !B

    slopesS(1:a,:)=slopes
    slopesS(0,:)=slopes(1,:) !0._dp
    slopesS(a+1,:)=slopes(a,:) !0._dp

    !Calculate VDN, which is the lateral discharge V*D at each point.
    !Based on the idea that
    ! VD = int( -dY/dt - dUD/dx ) dy
    !To ensure conservation, we calculate dY/dt to ensure that over the whole
    !cross-section, int(VD)=0
    VDN=0._dp 
    do i=2, b-1
    !Calculate the integral of dUd/dx. Note that we use u(i)+1 as a trick
    do j=max(l(i),2), min(u(i)+1,a)
    !do j=l(i), u(i)+1
    VDN(j,i) = VDN(j-1,i) + 0.5_dp*(-vels_old(j,i+1)*max(elevs_old(i+1)-hs(j,i+1),0._dp) - &
    vels_old(j-1,i+1)*max(elevs_old(i+1)-hs(j-1,i+1),0._dp) + & 
    vels_old(j,i-1)*max(elevs_old(i-1)-hs(j,i-1), 0._dp) + &
    vels_old(j-1,i-1)*max(elevs_old(i-1)-hs(j-1,i-1), 0._dp) )& 
    /(2._dp*delX)*(lengths(j,i)-lengths(j-1,i))
    end do
    !Now calculate dY/dT and correct - note the use of VDN(u(i)+1,i), which is later
    !set to zero
    VDN(l(i):u(i),i) = VDN(l(i):u(i),i) - &
    (lengths(l(i):u(i),i)-(lengths(max(l(i)-1,1),i) ))*&
    (VDN(min(u(i)+1,a),i)/((lengths(min(u(i)+1,a),i))-(lengths(max(l(i)-1,1),i) ) )) !Note that the last (VDN(u(i).. ) is dY/dt 

    VDN(min(u(i)+1,a),i)=0._dp
    !VDN(1:l(i)-1,i) = 999._dp
    !VDN(u(i)+1:a,i)=999._dp

    end do
    !print*, 'BUG'
    !VDN=0._dp
    !write(100,*) VDN(:,20) !, (-vels_old(:,21)*max(elevs_old(21)-hs(:,21),0._dp) + & 
    !vels_old(:,19)*max(elevs_old(19)-hs(:,19), 0._dp))

    !open(133,file='VDN')
    !write(133,*) VDN

    !!!!!!First half time step - implicit in the longitudinal coord - note how the lateral boundary condition is
    !automatically included - need to include the longitudinal boundary conditions
    !on the RHS
    band=0._dp
    do j=1,a
    !print*, j, size(Qe_old)
    !!Calculate matrix diagonals and RHS - upwind method
    do i=1,b
    if(veldH(j,i)>0._dp) THEN
    !Note the storage issues - because of the way I have indexed veldh(-1:b+2), I
    !need to add +2 to each of the index's

    !Upper-upper-diagonal
    band(KL+KU+1+i+2-(i+2 +2),i+2+2)= 0._dp !Advection

    !!Upper-diagonal
    band(KL+KU+1+i+2-(i+1+2),i+1+2) = 1._dp/delX*1._dp/6._dp*2._dp*(veldH(j,i+1))  !Advection

    !!Main diagonal
    band(KL+KU+1+i+2-(i+2),i+2) = depthH(j,i)*2._dp/DT +wset*cbedS(j,i)  & !Unsteady
    +1._dp/delX*1._dp/6._dp*3._dp*(veldH(j,i)) !Advection

    !!Lower-diagonal
    band(KL+KU+1+i+2-(i-1+2),i-1+2)= -1._dp/delX*1._dp/6._dp*6._dp*(veldH(j,i-1))  !Advection

    !!Lower-lower-diagonal
    band(KL+KU+1+i+2-(i-2+2),i-2+2)= 1._dp/delX*1._dp/6._dp*1._dp*(veldH(j,i-2)) !Advection

    !        upperb(i)= 0._dp !veldH(j,2:b+1)*.5_dp/delX
    !        diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) + veldh(j,i)/delX !!
    !        lowerb(i)= -veldh(j,i-1)/delX !veldh(j,0:b-1)*.5_dp/delX
    else
    if(depthH(j,i)>0._dp) THEN
    !!Upper-upper diagonal
    band(KL+KU+1+i+2-(i+2+2),i+2+2)=  -1._dp/delX*1._dp/6._dp*1._dp*(veldH(j,i+2)) !Advection

    !!Upper-diagonal
    band(KL+KU+1+i+2-(i+1+2),i+1+2) = 1._dp/delX*1._dp/6._dp*6._dp*(veldH(j,i+1))  !Advection

    !!Main diagonal
    band(KL+KU+1+i+2-(i+2),i+2) = depthH(j,i)*2._dp/DT +wset*cbedS(j,i) & !Unsteady
    -1._dp/delX*1._dp/6._dp*3._dp*(veldH(j,i))  !Advection

    !!Lower-diagonal
    band(KL+KU+1+i+2-(i-1+2),i-1+2)= -1._dp/delX*1._dp/6._dp*2._dp*(veldH(j,i-1))  !Advection

    !!Lower-lower-diagonal
    band(KL+KU+1+i+2-(i-2+2),i-2+2)= 0._dp !Advection

    !        upperb(i)= veldh(j,i+1)/delX !veldh(j,2:b+1)*.5_dp/delX
    !        diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) - veldh(j,i)/delX !!
    !        lowerb(i)= 0._dp !veldh(j,0:b-1)*.5_dp/delX
    ELSE
    !Set sconc to 0
    band(KL+KU+1+i+2-(i+2),i+2) = 1._dp !depthN(j,i)*2._dp/dT !depthH(j,i)*2._dp/DT +wset*cbedS(j,i) & !Unsteady
    !-1._dp/delX*1._dp/6._dp*3._dp*(veldH(j,i))  !Advection
    END IF

    end if
    end do
    rhsb=0._dp
    rhsb(1:b)= Qe_old(j,:)*sqrt(1._dp+slopes(j,:)**2)+ depthN(j,:)*sus_N(j,:)*2._dp/DT + 2._dp/(lnths(j+1,:)-lnths(j-1,:))*&  !!
    (0.5_dp*(edN(j,:)+edN(j+1,:))*(depthN2(j+1,:)*sus_N(j+1,:)-depthN2(j,:)*sus_N(j,:) )/(lnths(j+1,:)-lnths(j,:)) - & 
    0.5_dp*(edN(j,:)+edN(j-1,:))*(depthN2(j,:)*sus_N(j,:)-depthN2(j-1,:)*sus_N(j-1,:) )/(lnths(j,:)-lnths(j-1,:))  +&
    ((edN(j+1,:)+edN(j,:))*.5_dp*(cbedS(j+1,:)*sus_N(j+1,:)+cbedS(j,:)*sus_N(j,:))*.5_dp*(slopesS(j+1,:)+slopesS(j,:))*.5_dp - &
     (edN(j,:)+edN(j-1,:))*.5_dp*(cbedS(j,:)*sus_N(j,:)+cbedS(j-1,:)*sus_N(j-1,:))*.5_dp*(slopesS(j,:)+slopesS(j-1,:))*.5_dp) &
    !(.5_dp*(edN(j+1,:)*cbedS(j+1,:)*sus_N(j+1,:)*slopesS(j+1,:)+edN(j,:)*cbedS(j,:)*sus_N(j,:)*slopesS(j,:)) - &
    ! (.5_dp*(edN(j,:)*cbedS(j,:)*sus_N(j,:)*slopesS(j,:)+edN(j-1,:)*cbedS(j-1,:)*sus_N(j-1,:))*slopesS(j-1,:))) &
    ) - (VDN(j+1,:)*sus_N(j+1,:) - VDN(j-1,:)*sus_N(j-1,:))/(lnths(j+1,:)-lnths(j-1,:))

    do i=1,b
    if(depthH(j,i)==0._dp) rhsb(i)=0._dp
    end do

    !Boundary conditions at the mouth
    if(veldH(j,1)<0._dp) THEN
    !Here we say that sus_h(j,-1)=sus_h(j,0)=sus_h(j,1)
    i=-1
    band(KL+KU+1+i+2-(i+2+2),i+2+2) = -1._dp !Coefficient of sus_h(j,1) 
    band(KL+KU+1+i+2-(i+2),i+2) = 1._dp !Coefficient of sus_h(j,-1)

    i=0
    band(KL+KU+1+i+2-(i+1+2),i+1+2) = -1._dp !Coefficient of sus_h(j,1) 
    band(KL+KU+1+i+2-(i+2),i+2) = 1._dp !Coefficient of sus_h(j,0)
    else
    !Here we set the boundary points to Cmouth(j)
    i=-1
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp
    rhsb(i)=Cmouth(j)

    i=0
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp
    rhsb(i)=Cmouth(j)

    end if

    if(veldH(j,b)>0._dp) THEN
    !Here we set the boundary points sus_h(j,b+1), sus_h(j,b+2) == sus_h(j,b)
    i=b+1
    band(KL+KU+1+i+2-(i-1+2),i-1+2)=-1._dp !Coefficient for sus_h(j,b)
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp !Coefficient for sus_h(j,b+1)

    i=b+2
    band(KL+KU+1+i+2-(i-2+2),i-2+2)=-1._dp !Coefficient for sus_h(j,b)
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp !Coefficient for sus_h(j,b21)

    else
    !Here we set the boundaries to Criver
    i=b+1
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp !Coefficient for sus_h(j,b+1)
    rhsb(i)=Criver(j)

    i=b+2
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp !Coefficient for sus_h(j,b+2)
    rhsb(i)=Criver(j)

    end if

    !if(j==a/2)THEN 
    !       IF(minval(sus_h(j-1,:))<0._dp) THEN
    !        open(131,file='banded')
    !        open(132,file='rhs')
    !        write(131,*) band
    !        close(131)
    !        write(132,*) rhsb
    !        close(132)
    !        stop
    !        end if
    !end if

    call DGBSV(b+4, 2, 2, 1, band, 2*KL+KU+1, IPV, rhsb((-1):(b+2)), b+4, info)
    !print*, 'info=', info
    !The solver writes the new C to r- so we change it back.
    !sus_h(j,1:b)=rhsb(1:b)
    sus_h(j,:)=rhsb(:)



    if(info.ne.0) print*, "matrix problem in susconc2d, first bit", info
            
            !CHECKS
            do ii=1,b
            if(isnan(sus_h(j,ii))) THEN
            print*, 'sus_h(',j,ii,') is NAN' 
            stop
            end if
            
            !if(sus_h(j,ii)<0._dp) THEN
            !        if(sus_h(j,ii)< - 1E-04_dp) print*, 'significant negative sus_h', j, ii, sus_h(j,ii), sus_h(700,ii)
            !        sus_h(j,ii)=0._dp
            !From limited experimentation, it seems that even if this goes negative,
            !we don't end up with significant negative values in the final Cdist-
            !thus, maybe best to leave it alone.
            !end if
            
           ! if(depthH(j,ii)==0) THEN
           !         if(sus_h(j,ii)>0._dp) THEN
           !                if(sus_h(j,ii)>0.1_dp) THEN
           !                 print*, 'zero depth, high sus_h', sus_h(j,ii), j, ii,& 
     !depthH(j,ii)*2._dp/DT +wset*cbedS(j,ii) ! + veldh(j,ii)/delX, veldh(j,i+1)/delX, veldh(j,i-1)/delX 
           !                 stop
           !                 end if

           !                 sus_h(j,ii)=0._dp
           !        end if
           ! end if
            end do

    !        if((sus_h(j,b).ne.Criver(j)).and.(depthH(j,b)>0._dp)) THEN
    !               print*, 'sus_h(',j,"b) not equal to Criver", sus_h(j,b),Criver(j), depthH(j,b), depthN(j,b)&
    !,depthH(j,b)*2._dp/DT +wset*cbedS(j,b) - veldh(j,b)/delX,veldh(j,b)/delX, tmp, sus_N(j,b)
    !        end if


    end do
    !stop
    !print*, sus_h(1,110), sus_h(350,110)

    !print*, "half", sus_h(175,1), sus_h(175,40)
     
    !!!!!!!!!!!!!!End of first half time step


    !do i=1, b
    !if((sus_h(l(i),i)>3._dp*sus_h(l(i)+1,i)).and.(abs(sus_h(l(i),i))>1.0E-10_dp) ) print*, 'sus_h spike', i,l(i), & 
    !elevs_old(i)-hs(l(i)-1:l(i)+1,i),"##", sus_h(l(i):l(i)+1,i), sus_h(700,i)
    !end do


    !!!These boundary conditions are needed for the next half step
    sus_h(:,0)= Cmouth !B
    sus_h(:,b+1)=Criver !B
    sus_h(:,-1)= Cmouth !B
    sus_h(:,b+2)=Criver !B

    !!!Begin second half time step

    !!Eddy diffusivity
    edF(1:a,:)=dedy*abs(vels)*sqrt(fs)*sqrt(0.125_dp)*depthF

    if(periodicb.eqv..false.) THEN
    edF(0,:)=0._dp
    edF(a+1,:)=0._dp
    !Make things cancel at the boundary
    !do ii=1,b
    !edF(l(ii)-1,ii)=-edF(l(ii),ii)
    !edF(u(ii)+1,ii)=-edF(u(ii),ii)
    !end do

    !edF(0,:)=-edF(1,:) !0._dp
    !edF(a+1,:)=-edF(a,:)!0._dp
    ELSE
    edF(0,:)=edF(1,:)
    edF(a+1,:)=edF(a,:)
    END IF
    !!Eddy diffusivity*depth at last time step
    eddF= edF(1:a,:)*depthF !Note that the sqrt(0.125_dp) is sqrt(1/8) - note also how we ignore the change in the friction factor here

    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !cbedS= depthN*wset/((eddN/depthN)*(1._dp-exp(-wset/(eddN/depthN)*depthN) ) )
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            if(wset>0._dp) THEN
            tmp = ((edF(j,i)/dedy*dedz)*(1._dp-exp(-wset/(edF(j,i)/dedy*dedz)*depthF(j,i) ) ))
            else
            tmp=0._dp
            end if
        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthF(j,i)*wset/tmp, 200._dp)
        else
        cbedS(j,i)=200._dp
        if((wset==0).and.(edF(j,i)> 0._dp)) cbedS(j,i)=1._dp
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp

    !print*, cbedS(175,75), cbedS(175,75)*sus_h(175,75)
    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times depth, with boundaries
    eddF_yh(1:a-1,:)= .5_dp*(eddF(2:a,:)+eddF(1:a-1,:))
    if(periodicb.eqv..false.) THEN
    eddF_yh(0,:)=.5_dp*eddF(1,:) !B
    eddF_yh(a,:)=.5_dp*eddF(a,:) !B
    else
    eddF_yh(0,:)=.5_dp*(eddF(1,:)+eddF(a,:)) !B
    eddF_yh(a,:)=.5_dp*(eddF(a,:)+eddF(1,:)) !B
    end if



    !Calculate VDF, which is the lateral discharge V*D at each point.
    !Based on the idea that
    ! VD = int( -dY/dt - dUD/dx ) dy
    !To ensure conservation, we calculate dY/dt to ensure that over the whole
    !cross-section, int(VD)=0
    VDF=0._dp 
    do i=2, b-1
    !Calculate the integral of dUd/dx. Note that we use u(i)+1 as a trick
    do j=max(l(i),2), min(u(i)+1,a)
    !do j=l(i), u(i)+1
    VDF(j,i) = VDF(j-1,i) + 0.5_dp*(-vels(j,i+1)*max(elevs(i+1)-hs(j,i+1),0._dp) - &
     vels(j-1,i+1)*max(elevs(i+1)-hs(j-1,i+1),0._dp) + & 
    vels(j,i-1)*max(elevs(i-1)-hs(j,i-1), 0._dp) + &
    vels(j-1,i-1)*max(elevs(i-1)-hs(j-1,i-1), 0._dp) )& 
    /(2._dp*delX)*(lengths(j,i)-lengths(j-1,i))
    end do
    !Now calculate dY/dT and correct - note the use of VDF(u(i)+1,i), which is later
    !set to zero
    VDF(l(i):u(i),i) = VDF(l(i):u(i),i) - &
    (lengths(l(i):u(i),i)-(lengths(max(l(i)-1,1),i)) )*&
    (VDF(min(u(i)+1,a),i)/((lengths(min(u(i)+1,a),i))-(lengths(max(l(i)-1,1),i) ) )) !Note that the last (VDN(u(i).. ) is dY/dt 

    VDF(min(u(i)+1,a),i)=0._dp
    !VDN(1:l(i)-1,i) = 999._dp
    !VDN(u(i)+1:a,i)=999._dp

    end do
    !print*, 'BUG'
    !VDF=0._dp



    do i=1,b

    !Define some dy type variables
    dyc(2:a-1)= .5_dp*(lengths(3:a,i)-lengths(1:a-2,i))
    dyc(1)=dyc(2) !0.5_dp*(lengths(2,i)-(lengths(1,i)-depthF(1,i)/max(abs(slopes(1,i)), 0.01_dp)))!dyc(2) !B
    dyc(a)=dyc(a-1)!0.5_dp*((lengths(a,i)+depthF(a,i)/max(abs(slopes(a,i)), 0.01_dp))-lengths(a-1,i))!dyc(a-1) !B

    dy(1:a-1)=lengths(2:a,i)-lengths(1:a-1,i)
    dy(a)=dy(a-1) !depthF(a,i)/max(abs(slopes(a,i)), 0.01_dp) !dy(a-1)
    dy(0)=dy(1) !depthF(1,i)/max(abs(slopes(1,i)), 0.01_dp) !dy(1)

    !!Matrix diagonals
    uppera(1:a-1)= -1._dp/dyc(1:a-1)*(0.5_dp*(edF(1:a-1,i)+edF(2:a,i))*depthF(2:a,i)/dy(1:a-1)&
    +.5_dp*edF(2:a,i)*cbedS(2:a,i)*slopesS(2:a,i)) &
    + VDF(2:a,i)/(2._dp*dyc(1:a-1))
    lowera(2:a)= -1._dp/dyc(2:a)*(0.5_dp*(edF(1:a-1,i)+edF(2:a,i))*depthF(1:a-1,i)/dy(1:a-1)& 
    -.5_dp*edF(1:a-1,i)*cbedS(1:a-1,i)*slopesS(1:a-1,i)) &
    -VDF(1:a-1,i)/(2._dp*dyc(2:a))
    diaga(2:a-1)= depthF(2:a-1,i)/DT*2._dp +&
    1._dp/dyc(2:a-1)*(0.5_dp*( edF(3:a,i)+edF(2:a-1,i) )*depthF(2:a-1, i)/dy(2:a-1)&
    -.5_dp*edF(2:a-1,i)*cbedS(2:a-1,i)*slopesS(2:a-1,i))+&
    1._dp/dyc(1:a-2)*(0.5_dp*(edF(1:a-2,i)+edF(2:a-1,i))*depthF(2:a-1, i)/dy(1:a-2)&
    +.5_dp*edF(2:a-1,i)*cbedS(2:a-1,i)*slopesS(2:a-1,i)) + &
    wset*cbedS(2:a-1,i) !!

    diaga(1)= depthF(1,i)/DT*2._dp + &
    1._dp/dyc(1)*(0.5_dp*(edF(2,i)+edF(1,i))*depthF(1,i)/dy(1)- 0.5_dp*edF(1,i)*cbedS(1,i)*slopesS(1,i))
    diaga(a)= depthF(a,i)/DT*2._dp +&
    1._dp/dyc(a)*(0.5_dp*(edF(a,i)+edF(a-1,i))*depthF(a,i)/dy(a-1) + 0.5_dp*edF(a,i)*cbedS(a,i)*slopesS(a,i))

    !Right hand side, upwind method
    do j=1,a
    if(veldh(j,i)>0._dp) THEN
    rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
    (2._dp*veldh(j,i+1)*sus_h(j,i+1)+3._dp*veldh(j,i)*sus_h(j,i) - 6._dp*veldh(j,i-1)*sus_h(j,i-1) &
    +  veldh(j,i-2)*sus_h(j,i-2))/(6._dp*delX) !!

    !rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
    !(veldh(j,i)*sus_h(j,i) - veldh(j,i-1)*sus_h(j,i-1))/delX !!
    ELSE
    rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
    (-veldh(j,i+2)*sus_h(j,i+2) +6._dp*veldh(j,i+1)*sus_h(j,i+1) -3._dp*veldh(j,i)*sus_h(j,i) &
     -2._dp*veldh(j,i-1)*sus_h(j,i-1) )/(6._dp*delX) !!

    !rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
    !(veldh(j,i+1)*sus_h(j,i+1) - veldh(j,i)*sus_h(j,i))/delX !!
    END IF
    END DO
    !For boundary conditions here, it is appropriate that the sus is zero on the
    !lateral boundaries-  which equates to us not having to enforce anything
    !further. Unless we have periodic boundaries, in which case the following might
    !be a practical alternative
    !if(periodicb) THEN
    !rhsa(1)=rhsa(1) -lowera(1)*sus_h(1,b)
    !rhsa(a)= rhsa(a)-uppera(a)*sus_h(a,b)
    !END IF
     

    !!Prevent main diagonal from being zero
    do ii=1,a
    if((diaga(ii)==0._dp).or.(depthF(ii,i)==0._dp)) THEN
            if( (ii==1).or.(ii==a).or.(depthF(ii,i)==0._dp)) THEN
            diaga(ii)=1._dp
            lowera(ii)=0._dp
            uppera(ii)=0._dp
            rhsa(ii)=0._dp
            else
            diaga(ii)=1._dp
            lowera(ii)=-0.5_dp
            uppera(ii)=-0.5_dp
            rhsa(ii)=0._dp
            end if
    end if
    end do

    call DGTSV(a,1, lowera(2:a), diaga, uppera(1:a-1), rhsa,a, info)

    Cdist(:,i)= rhsa

    if(info.ne.0) print*, "matrix problem in susconc2d, second bit", info


    do ii=1, a
    !if(sus_h(ii,i)<0._dp) THEN
    !        if(sus_h(ii,i)< - 1E-08_dp) print*, 'significant negative sus_h', j, ii, sus_h(ii,i), sus_h(a/2,i)
    !        !Cdist(j,ii)=0._dp
    !end if
    if(Cdist(ii,i)<0._dp) THEN
    !        if(Cdist(ii,i)< - 1E-08_dp) print*, 'significant negative Cdist', j, ii, Cdist(ii,i), Cdist(a/2,i), &
    !depthF(ii,i),depthH(ii,i)
            Cdist(ii,i)=0._dp
    end if
    end do

    !if(minval(Cdist(:,i))<0._dp) THEN
    !		!Cdist(:,i)=max(Cdist(:,i),0._dp)
    !	if(abs(minval(Cdist(:,i)))< 1.00E-12) THEN
    !
    !	ElSE
        !print*, "Cdist negative with non-roundoff type magnitude"
        !print*, "min Cdist", i, "<0", minval(Cdist(:,i)), minloc(Cdist(:,i)), depthF(minloc(Cdist(:,i)), (i-1):(i+1))& 
        !, vels(minloc(Cdist(:,i)), (i-1):(i+1)) 
        !END IF
    !
    !end if


    end do


    !if(minval(Cdist)<0._dp) print*,'Cdist<0 end', minloc(Cdist), minval(Cdist)

    Cdist_old=cbedS(1:a,:) !At the moment I am just using this to check the output. 
    !do i=1, b
    !if((Cdist(l(i),i)>3._dp*Cdist(l(i)+1,i)).and.(abs(Cdist(l(i),i))>1.0E-10_dp ) ) print*, 'Cdist spike end', i,l(i), & 
    !elevs(i)-hs(l(i)-1:l(i)+1,i), '##', Cdist(l(i):l(i)+1,i), Cdist(700,i)  
    !end do

    !print*, '2nd', Cdist(1,110), Cdist(350,110)

    !stop

end subroutine susconc2dup3rdV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Numerical method is upwind in the x derivative advective term, but still ADI
!otherwise
subroutine susconc2dup3rdV2(a,b, DT,delX, lengths,elevs,elevs_old,vels,vels_old,hs, Cdist,Cdist_old, Qe,Qe_old, Cmouth, Criver, & 
wset,fs, slopes, l, u,lambdacon)
    !use matrix_solvers
    !use crosssection

    !implicit none
    INTEGER, intent(in):: a, b, l, u
    REAL(dp), intent(in):: DT,delX, wset
    REAL(dp), intent(in):: elevs, elevs_old, vels, vels_old,hs, Qe,Qe_old, Cmouth, Criver,fs, lengths, slopes,lambdacon
    REAL(dp), intent(in out):: Cdist, Cdist_old !!Cdist_old will output random interesting information
    DIMENSION elevs(b), elevs_old(b), vels(a,b), vels_old(a,b), hs(a,b), Cdist(a,b), Cdist_old(a,b), Qe(a,b),fs(a,b), lengths(a,b) &
    ,Criver(a), Cmouth(a), Qe_old(a,b), slopes(a,b), l(b), u(b)

    INTEGER:: i, j, info, ii, KL=2, KU=2, IPV((-1):(b+2))
    logical:: periodicb=.false.!false.
    REAL(dp):: depthN(a,b), depthH(a,b), depthF(a,b), elevsH, velh(a,b),sus_h(a,(-1):(b+2)),edN(0:a+1,b),eddN(a,b),& 
                veldh(a,(-1):(b+2))
    REAL(dp):: sus_N(0:a+1,b), lnths(0:a+1,b), eddN_yh(0:a,b) ,edF(0:a+1,b), eddF(a,b), eddF_yh(0:a,b)
    REAL(dp):: upperb(b), diagb(b),lowerb(b), rhsb((-1):(b+2)), depthN2(0:a+1,b), VDN(0:a+1,b), VDF(a,b)
    REAL(dp):: uppera(a), diaga(a),lowera(a), rhsa(a), dyc(a), dy(0:a), sus_F(a,b), cbedS(0:a+1,b), slopesS(0:a+1,b)
    REAL(dp):: dedy, dedz=0.10_dp !Dimensionless eddy diffusivity
    REAL(dp):: tmpa(a), tmpb(b),tmp, band(7,b+4)

    !!!!So this solves the equation
    ! d/dt(depth*C) + d/dx(vel*depth*C) +d/dy(vel_lateral*depth*C)= d/dy(eddn*depth*dC/dy + eddn*Cbed*dh/dy)
    ! where I plan to add the last term soon.
    !By the ADI method, except that the convective x derivative is treated as upwind
    !(seems that it has to be to get stable results). 
    !!Note that assuming a constant vertical eddy diffusivity, 
    !Cbed= Caverage*depth/((ez/ws)*(1-exp(-vs/ez*(depth) ) ))
    !!And I have been assuming that ez=0.5_dp*ey

    !There are ghost points on the left and right edges (i.e. the margins of each
    !cross-section -- 1 each), and on the north
    !and south edges (i.e. the upstream and downstream boundaries -- 2 each)

    !Dimensionless eddy diffusivity
    dedy=lambdacon

    !!!First define water depths
    do i= 1,b
    !Depth at previous time step
    depthN(:,i)= max(elevs_old(i)-hs(:,i),0._dp)
    !Depth at half time step
    elevsH=(elevs(i)+elevs_old(i))*.5_dp
    depthH(:,i)= max(elevsH-hs(:,i), 0._dp)
    !Depth at next time step
    depthF(:,i)= max(elevs(i)-hs(:,i),0._dp)
    end do

    !!Useful thing
    depthN2(1:a,:)=depthN
    depthN2(0,:)=0._dp
    depthN2(a+1,:)=0._dp
    !!


    !Now calculate half time step velocity
    velH=.5_dp*(vels+vels_old)
    !Halfway depth times velocity - a useful efficiency device
    veldH(:,1:b)=depthH*velH
    veldH(:,0)= veldH(:,1) !B
    veldH(:,b+1)=veldH(:,b) !B
    veldH(:,-1)=veldH(:,0)
    veldH(:,b+2)=veldH(:,b) !B

    !Eddy diffusivity

    edN(1:a,:)= dedy*abs(vels_old)*sqrt(fs)*sqrt(0.125_dp)*depthN

    if(periodicb.eqv..false.) THEN
    edN(0,:)=0._dp
    edN(a+1,:)=0._dp

    !Make sure things cancel at the boundaries appropriately
    !do ii=1,b
    !edN(l(ii)-1, ii)=-edN(l(ii),ii)
    !edN(u(ii)+1,ii)=-edN(u(ii),ii)
    !end do
    !edN(0,:)=-edN(1,:)
    !edN(a+1,:)=-edN(a,:)
    ELSE
    edN(0,:)=edN(a,:)
    edN(a+1,:)=edN(1,:)
    end if

    !!Eddy diffusivity*depth at last time step
    eddN= edN(1:a,1:b)*depthN !Note that the sqrt(0.125_dp) is sqrt(1/8)



    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            
            if(edN(j,i)>0._dp) THEN
           tmp= ((edN(j,i)/dedy*dedz)*(1._dp-exp(-wset/(edN(j,i)/dedy*dedz)*depthN(j,i)) ) )
            ELSE
            tmp=0._dp
            END IF

        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthN(j,i)*wset/tmp, 200._dp)
        else
            cbedS(j,i)=200._dp  !So this is the maximum scaling factor we offer
            if((wset==0._dp).and.(edN(j,i)> 0._dp)) cbedS(j,i)=1._dp  !So the idea is that this applies  when wset=0
        
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp

    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times the depth, with boundaries
    eddN_yh(1:a-1,:)= .5_dp*(eddN(2:a,:)+eddN(1:a-1,:))

    if(periodicb.eqv..false.) THEN
    eddN_yh(0,:)=.5_dp*eddN(1,:) !B
    eddN_yh(a,:)=.5_dp*eddN(a,:) !B
    else
    eddN_yh(0,:)=.5_dp*(eddN(1,:)+eddN(a,:)) !B
    eddN_yh(a,:)=.5_dp*(eddN(a,:)+eddN(1,:)) !B
    end if

    !!Susconc at previous time step, with suitable boundary conditions attached
    sus_N(1:a,:)= Cdist_old
    !Zero susconc along the lateral boundaries -- might not be true - zero
    !integrated concentration, but the concentration might not drop off entirely
    sus_N(0,:)=sus_N(1,:) !0._dp !B
    sus_N(a+1,:)=sus_N(a,:) !0._dp !B


    !Lengths, with suitable boundaries attached
    lnths(1:a,:)= lengths
    lnths(0,:)= 2._dp*lengths(1,:)-lengths(2,:)!lengths(1,:)-depthN(1,:)/max(abs(slopes(1,:)), 0.01_dp)!2._dp*lengths(1,:)-lengths(2,:) !B
    lnths(a+1,:)=  2._dp*lengths(a,:)-lengths(a-1,:)!lengths(a,:) + depthN(a,:)/max(abs(slopes(a,:)), 0.01_dp)!2._dp*lengths(a,:)-lengths(a-1,:) !B

    slopesS(1:a,:)=slopes
    slopesS(0,:)=slopes(1,:) !0._dp
    slopesS(a+1,:)=slopes(a,:) !0._dp

    !Calculate VDN, which is the lateral discharge V*D at each point.
    !Based on the idea that
    ! VD = int( -dY/dt - dUD/dx ) dy
    !To ensure conservation, we calculate dY/dt to ensure that over the whole
    !cross-section, int(VD)=0
    VDN=0._dp 
    do i=2, b-1
    !Calculate the integral of dUd/dx. Note that we use u(i)+1 as a trick
    do j=max(l(i),2), min(u(i)+1,a)
    !do j=l(i), u(i)+1
    VDN(j,i) = VDN(j-1,i) + 0.5_dp*(-vels_old(j,i+1)*max(elevs_old(i+1)-hs(j,i+1),0._dp) - &
    vels_old(j-1,i+1)*max(elevs_old(i+1)-hs(j-1,i+1),0._dp) + & 
    vels_old(j,i-1)*max(elevs_old(i-1)-hs(j,i-1), 0._dp) + &
    vels_old(j-1,i-1)*max(elevs_old(i-1)-hs(j-1,i-1), 0._dp) )& 
    /(2._dp*delX)*(lengths(j,i)-lengths(j-1,i))
    end do
    !Now calculate dY/dT and correct - note the use of VDN(u(i)+1,i), which is later
    !set to zero
    VDN(l(i):u(i),i) = VDN(l(i):u(i),i) - &
    (lengths(l(i):u(i),i)-(lengths(max(l(i)-1,1),i) ))*&
    (VDN(min(u(i)+1,a),i)/((lengths(min(u(i)+1,a),i))-(lengths(max(l(i)-1,1),i) ) )) !Note that the last (VDN(u(i).. ) is dY/dt 

    VDN(min(u(i)+1,a),i)=0._dp
    !VDN(1:l(i)-1,i) = 999._dp
    !VDN(u(i)+1:a,i)=999._dp

    end do
    !print*, 'BUG'
    !VDN=0._dp
    !write(100,*) VDN(:,20) !, (-vels_old(:,21)*max(elevs_old(21)-hs(:,21),0._dp) + & 
    !vels_old(:,19)*max(elevs_old(19)-hs(:,19), 0._dp))

    !open(133,file='VDN')
    !write(133,*) VDN

    !!!!!!First half time step - implicit in the longitudinal coord - note how the lateral boundary condition is
    !automatically included - need to include the longitudinal boundary conditions
    !on the RHS
    do j=1,a
    band=0._dp
    !print*, j, size(Qe_old)
    !!Calculate matrix diagonals and RHS - upwind method
    DO i=1,b
    if(depthH(j,i)>0._dp) THEN
            if(0.5_dp*(veldH(j,i)+veldH(j,i+1))>0._dp) THEN
                    if(0.5_dp*(veldH(j,i)+veldH(j,i-1))>0._dp) THEN
                    !veldH(j,i+1/2)>0, veldH(j,i-1/2)>0
                    !Note the storage issues - because of the way I have indexed veldh(-1:b+2), I
                    !need to add +2 to each of the index's

                    !Upper-upper-diagonal
                    band(KL+KU+1+i+2-(i+2 +2),i+2+2)= 0._dp !Advection

                    !!Upper-diagonal
                    band(KL+KU+1+i+2-(i+1+2),i+1+2) = 1._dp/delX*1._dp/6._dp*2._dp*(veldH(j,i+1))  !Advection

                    !!Main diagonal
                    band(KL+KU+1+i+2-(i+2),i+2) = depthH(j,i)*2._dp/DT +wset*cbedS(j,i)  & !Unsteady
                    +1._dp/delX*1._dp/6._dp*3._dp*(veldH(j,i)) !Advection

                    !!Lower-diagonal
                    band(KL+KU+1+i+2-(i-1+2),i-1+2)= -1._dp/delX*1._dp/6._dp*6._dp*(veldH(j,i-1))  !Advection

                    !!Lower-lower-diagonal
                    band(KL+KU+1+i+2-(i-2+2),i-2+2)= 1._dp/delX*1._dp/6._dp*1._dp*(veldH(j,i-2)) !Advection

                    !        upperb(i)= 0._dp !veldH(j,2:b+1)*.5_dp/delX
                    !        diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) + veldh(j,i)/delX !!
                    !        lowerb(i)= -veldh(j,i-1)/delX !veldh(j,0:b-1)*.5_dp/delX
                    else
                    !veldH(j,i+1/2)>0, veldH(j,i-1/2)<=0
                    
                    band(KL+KU+1+i+2-(i+2 +2),i+2+2)= 0._dp !Advection

                    !!Upper-diagonal
                    band(KL+KU+1+i+2-(i+1+2),i+1+2) = 1._dp/delX*1._dp/6._dp*3._dp*(veldH(j,i+1))  !Advection

                    !!Main diagonal
                    band(KL+KU+1+i+2-(i+2),i+2) = depthH(j,i)*2._dp/DT +wset*cbedS(j,i)  & !Unsteady
                    +1._dp/delX*1._dp/6._dp*0._dp*(veldH(j,i)) !Advection

                    !!Lower-diagonal
                    band(KL+KU+1+i+2-(i-1+2),i-1+2)= -1._dp/delX*1._dp/6._dp*3._dp*(veldH(j,i-1))  !Advection

                    !!Lower-lower-diagonal
                    band(KL+KU+1+i+2-(i-2+2),i-2+2)= 1._dp/delX*1._dp/6._dp*0._dp*(veldH(j,i-2)) !Advection
                    
                    end if
            else !veldH(j,i+1/2)<=0._dp
                    
                    if(0.5_dp*(veldH(j,i)+veldH(j,i-1))<0._dp) THEN
                    !veldH(j,i+1/2)<=0, veldH(j,i-1/2)<0
                    !!Upper-upper diagonal
                    band(KL+KU+1+i+2-(i+2+2),i+2+2)=  -1._dp/delX*1._dp/6._dp*1._dp*(veldH(j,i+2)) !Advection

                    !!Upper-diagonal
                    band(KL+KU+1+i+2-(i+1+2),i+1+2) = 1._dp/delX*1._dp/6._dp*6._dp*(veldH(j,i+1))  !Advection

                    !!Main diagonal
                    band(KL+KU+1+i+2-(i+2),i+2) = depthH(j,i)*2._dp/DT +wset*cbedS(j,i) & !Unsteady
                    -1._dp/delX*1._dp/6._dp*3._dp*(veldH(j,i))  !Advection

                    !!Lower-diagonal
                    band(KL+KU+1+i+2-(i-1+2),i-1+2)= -1._dp/delX*1._dp/6._dp*2._dp*(veldH(j,i-1))  !Advection

                    !!Lower-lower-diagonal
                    band(KL+KU+1+i+2-(i-2+2),i-2+2)= 0._dp !Advection

                    !        upperb(i)= veldh(j,i+1)/delX !veldh(j,2:b+1)*.5_dp/delX
                    !        diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) - veldh(j,i)/delX !!
                    !        lowerb(i)= 0._dp !veldh(j,0:b-1)*.5_dp/delX
                    else
                    !veldH(j,i+1/2)<=0, veldH(j,i-1/2)>=0
                    band(KL+KU+1+i+2-(i+2+2),i+2+2)=  -1._dp/delX*1._dp/6._dp*1._dp*(veldH(j,i+2)) !Advection

                    !!Upper-diagonal
                    band(KL+KU+1+i+2-(i+1+2),i+1+2) = 1._dp/delX*1._dp/6._dp*5._dp*(veldH(j,i+1))  !Advection

                    !!Main diagonal
                    band(KL+KU+1+i+2-(i+2),i+2) = depthH(j,i)*2._dp/DT +wset*cbedS(j,i) & !Unsteady
                    -1._dp/delX*1._dp/6._dp*0._dp*(veldH(j,i))  !Advection

                    !!Lower-diagonal
                    band(KL+KU+1+i+2-(i-1+2),i-1+2)= -1._dp/delX*1._dp/6._dp*5._dp*(veldH(j,i-1))  !Advection

                    !!Lower-lower-diagonal
                    band(KL+KU+1+i+2-(i-2+2),i-2+2)= 1._dp/delX*1._dp/6._dp*1._dp*(veldH(j,i-2))  

                    !        upperb(i)= veldh(j,i+1)/delX !veldh(j,2:b+1)*.5_dp/delX
                    !        diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) - veldh(j,i)/delX !!
                    !        lowerb(i)= 0._dp !veldh(j,0:b-1)*.5_dp/delX
                    end if
            end if
    ELSE !depth=0
    !Set sconc to 0
    band(KL+KU+1+i+2-(i+2),i+2) = 1._dp !depthN(j,i)*2._dp/dT !depthH(j,i)*2._dp/DT +wset*cbedS(j,i) & !Unsteady
    !-1._dp/delX*1._dp/6._dp*3._dp*(veldH(j,i))  !Advection
    END IF
    END DO

    rhsb=0._dp
    rhsb(1:b)= Qe_old(j,:)*sqrt(1._dp+slopes(j,:)**2)+ depthN(j,:)*sus_N(j,:)*2._dp/DT + 2._dp/(lnths(j+1,:)-lnths(j-1,:))*&  !!
    (0.5_dp*(edN(j,:)+edN(j+1,:))*(depthN2(j+1,:)*sus_N(j+1,:)-depthN2(j,:)*sus_N(j,:) )/(lnths(j+1,:)-lnths(j,:)) - & 
    0.5_dp*(edN(j,:)+edN(j-1,:))*(depthN2(j,:)*sus_N(j,:)-depthN2(j-1,:)*sus_N(j-1,:) )/(lnths(j,:)-lnths(j-1,:))  +&
    ((edN(j+1,:)+edN(j,:))*.5_dp*(cbedS(j+1,:)*sus_N(j+1,:)+cbedS(j,:)*sus_N(j,:))*.5_dp*(slopesS(j+1,:)+slopesS(j,:))*.5_dp - &
     (edN(j,:)+edN(j-1,:))*.5_dp*(cbedS(j,:)*sus_N(j,:)+cbedS(j-1,:)*sus_N(j-1,:))*.5_dp*(slopesS(j,:)+slopesS(j-1,:))*.5_dp) &
    !(.5_dp*(edN(j+1,:)*cbedS(j+1,:)*sus_N(j+1,:)*slopesS(j+1,:)+edN(j,:)*cbedS(j,:)*sus_N(j,:)*slopesS(j,:)) - &
    ! (.5_dp*(edN(j,:)*cbedS(j,:)*sus_N(j,:)*slopesS(j,:)+edN(j-1,:)*cbedS(j-1,:)*sus_N(j-1,:))*slopesS(j-1,:))) &
    ) - (VDN(j+1,:)*sus_N(j+1,:) - VDN(j-1,:)*sus_N(j-1,:))/(lnths(j+1,:)-lnths(j-1,:))

    do i=1,b
    if(depthH(j,i)==0._dp) rhsb(i)=0._dp
    end do

    !Boundary conditions at the mouth
    if(veldH(j,1)<0._dp) THEN
    !Here we say that sus_h(j,-1)=sus_h(j,0)=sus_h(j,1)
    i=-1
    band(KL+KU+1+i+2-(i+2+2),i+2+2) = -1._dp !Coefficient of sus_h(j,1) 
    band(KL+KU+1+i+2-(i+2),i+2) = 1._dp !Coefficient of sus_h(j,-1)

    i=0
    band(KL+KU+1+i+2-(i+1+2),i+1+2) = -1._dp !Coefficient of sus_h(j,1) 
    band(KL+KU+1+i+2-(i+2),i+2) = 1._dp !Coefficient of sus_h(j,0)
    else
    !Here we set the boundary points to Cmouth(j)
    i=-1
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp
    rhsb(i)=Cmouth(j)

    i=0
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp
    rhsb(i)=Cmouth(j)

    end if

    if(veldH(j,b)>0._dp) THEN
    !Here we set the boundary points sus_h(j,b+1), sus_h(j,b+2) == sus_h(j,b)
    i=b+1
    band(KL+KU+1+i+2-(i-1+2),i-1+2)=-1._dp !Coefficient for sus_h(j,b)
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp !Coefficient for sus_h(j,b+1)

    i=b+2
    band(KL+KU+1+i+2-(i-2+2),i-2+2)=-1._dp !Coefficient for sus_h(j,b)
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp !Coefficient for sus_h(j,b21)

    else
    !Here we set the boundaries to Criver
    i=b+1
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp !Coefficient for sus_h(j,b+1)
    rhsb(i)=Criver(j)

    i=b+2
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp !Coefficient for sus_h(j,b+2)
    rhsb(i)=Criver(j)

    end if

    !if(j==a/2)THEN 
    !       IF(minval(sus_h(j-1,:))<0._dp) THEN
    !        open(131,file='banded')
    !        open(132,file='rhs')
    !        write(131,*) band
    !        close(131)
    !        write(132,*) rhsb
    !        close(132)
    !        stop
    !        end if
    !end if

    call DGBSV(b+4, 2, 2, 1, band, 2*KL+KU+1, IPV, rhsb((-1):(b+2)), b+4, info)
    !print*, 'info=', info
    !The solver writes the new C to r- so we change it back.
    !sus_h(j,1:b)=rhsb(1:b)
    sus_h(j,:)=rhsb(:)



    if(info.ne.0) print*, "matrix problem in susconc2d, first bit", info, j, depthH(j,info-2), depthH(j,info)
            
            !CHECKS
            do ii=1,b
            if(isnan(sus_h(j,ii))) THEN
            print*, 'sus_h(',j,ii,') is NAN' 
            stop
            end if
            
            !if(sus_h(j,ii)<0._dp) THEN
            !        if(sus_h(j,ii)< - 1E-04_dp) print*, 'significant negative sus_h', j, ii, sus_h(j,ii), sus_h(700,ii)
            !        sus_h(j,ii)=0._dp
            !From limited experimentation, it seems that even if this goes negative,
            !we don't end up with significant negative values in the final Cdist-
            !thus, maybe best to leave it alone.
            !end if
            
           ! if(depthH(j,ii)==0) THEN
           !         if(sus_h(j,ii)>0._dp) THEN
           !                if(sus_h(j,ii)>0.1_dp) THEN
           !                 print*, 'zero depth, high sus_h', sus_h(j,ii), j, ii,& 
     !depthH(j,ii)*2._dp/DT +wset*cbedS(j,ii) ! + veldh(j,ii)/delX, veldh(j,i+1)/delX, veldh(j,i-1)/delX 
           !                 stop
           !                 end if

           !                 sus_h(j,ii)=0._dp
           !        end if
           ! end if
            end do

    !        if((sus_h(j,b).ne.Criver(j)).and.(depthH(j,b)>0._dp)) THEN
    !               print*, 'sus_h(',j,"b) not equal to Criver", sus_h(j,b),Criver(j), depthH(j,b), depthN(j,b)&
    !,depthH(j,b)*2._dp/DT +wset*cbedS(j,b) - veldh(j,b)/delX,veldh(j,b)/delX, tmp, sus_N(j,b)
    !        end if


    end do
    !stop
    !print*, sus_h(1,110), sus_h(350,110)

    !print*, "half", sus_h(175,1), sus_h(175,40)
     
    !!!!!!!!!!!!!!End of first half time step


    !do i=1, b
    !if((sus_h(l(i),i)>3._dp*sus_h(l(i)+1,i)).and.(abs(sus_h(l(i),i))>1.0E-10_dp) ) print*, 'sus_h spike', i,l(i), & 
    !elevs_old(i)-hs(l(i)-1:l(i)+1,i),"##", sus_h(l(i):l(i)+1,i), sus_h(700,i)
    !end do


    !!!These boundary conditions are needed for the next half step
    sus_h(:,0)= Cmouth !B
    sus_h(:,b+1)=Criver !B
    sus_h(:,-1)= Cmouth !B
    sus_h(:,b+2)=Criver !B

    !!!Begin second half time step

    !!Eddy diffusivity
    edF(1:a,:)=dedy*abs(vels)*sqrt(fs)*sqrt(0.125_dp)*depthF

    if(periodicb.eqv..false.) THEN
    edF(0,:)=0._dp
    edF(a+1,:)=0._dp
    !Make things cancel at the boundary
    !do ii=1,b
    !edF(l(ii)-1,ii)=-edF(l(ii),ii)
    !edF(u(ii)+1,ii)=-edF(u(ii),ii)
    !end do

    !edF(0,:)=-edF(1,:) !0._dp
    !edF(a+1,:)=-edF(a,:)!0._dp
    ELSE
    edF(0,:)=edF(1,:)
    edF(a+1,:)=edF(a,:)
    END IF
    !!Eddy diffusivity*depth at last time step
    eddF= edF(1:a,:)*depthF !Note that the sqrt(0.125_dp) is sqrt(1/8) - note also how we ignore the change in the friction factor here

    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !cbedS= depthN*wset/((eddN/depthN)*(1._dp-exp(-wset/(eddN/depthN)*depthN) ) )
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            if(wset>0._dp) THEN
            tmp = ((edF(j,i)/dedy*dedz)*(1._dp-exp(-wset/(edF(j,i)/dedy*dedz)*depthF(j,i) ) ))
            else
            tmp=0._dp
            end if
        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthF(j,i)*wset/tmp, 200._dp)
        else
        cbedS(j,i)=200._dp
        if((wset==0).and.(edF(j,i)> 0._dp)) cbedS(j,i)=1._dp
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp

    !print*, cbedS(175,75), cbedS(175,75)*sus_h(175,75)
    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times depth, with boundaries
    eddF_yh(1:a-1,:)= .5_dp*(eddF(2:a,:)+eddF(1:a-1,:))
    if(periodicb.eqv..false.) THEN
    eddF_yh(0,:)=.5_dp*eddF(1,:) !B
    eddF_yh(a,:)=.5_dp*eddF(a,:) !B
    else
    eddF_yh(0,:)=.5_dp*(eddF(1,:)+eddF(a,:)) !B
    eddF_yh(a,:)=.5_dp*(eddF(a,:)+eddF(1,:)) !B
    end if



    !Calculate VDF, which is the lateral discharge V*D at each point.
    !Based on the idea that
    ! VD = int( -dY/dt - dUD/dx ) dy
    !To ensure conservation, we calculate dY/dt to ensure that over the whole
    !cross-section, int(VD)=0
    VDF=0._dp 
    do i=2, b-1
    !Calculate the integral of dUd/dx. Note that we use u(i)+1 as a trick
    do j=max(l(i),2), min(u(i)+1,a)
    !do j=l(i), u(i)+1
    VDF(j,i) = VDF(j-1,i) + 0.5_dp*(-vels(j,i+1)*max(elevs(i+1)-hs(j,i+1),0._dp) - &
     vels(j-1,i+1)*max(elevs(i+1)-hs(j-1,i+1),0._dp) + & 
    vels(j,i-1)*max(elevs(i-1)-hs(j,i-1), 0._dp) + &
    vels(j-1,i-1)*max(elevs(i-1)-hs(j-1,i-1), 0._dp) )& 
    /(2._dp*delX)*(lengths(j,i)-lengths(j-1,i))
    end do
    !Now calculate dY/dT and correct - note the use of VDF(u(i)+1,i), which is later
    !set to zero
    VDF(l(i):u(i),i) = VDF(l(i):u(i),i) - &
    (lengths(l(i):u(i),i)-(lengths(max(l(i)-1,1),i)) )*&
    (VDF(min(u(i)+1,a),i)/((lengths(min(u(i)+1,a),i))-(lengths(max(l(i)-1,1),i) ) )) !Note that the last (VDN(u(i).. ) is dY/dt 

    VDF(min(u(i)+1,a),i)=0._dp
    !VDN(1:l(i)-1,i) = 999._dp
    !VDN(u(i)+1:a,i)=999._dp

    end do
    !print*, 'BUG'
    !VDF=0._dp



    do i=1,b

    !Define some dy type variables
    dyc(2:a-1)= .5_dp*(lengths(3:a,i)-lengths(1:a-2,i))
    dyc(1)=dyc(2) !0.5_dp*(lengths(2,i)-(lengths(1,i)-depthF(1,i)/max(abs(slopes(1,i)), 0.01_dp)))!dyc(2) !B
    dyc(a)=dyc(a-1)!0.5_dp*((lengths(a,i)+depthF(a,i)/max(abs(slopes(a,i)), 0.01_dp))-lengths(a-1,i))!dyc(a-1) !B

    dy(1:a-1)=lengths(2:a,i)-lengths(1:a-1,i)
    dy(a)=dy(a-1) !depthF(a,i)/max(abs(slopes(a,i)), 0.01_dp) !dy(a-1)
    dy(0)=dy(1) !depthF(1,i)/max(abs(slopes(1,i)), 0.01_dp) !dy(1)

    !!Matrix diagonals
    uppera(1:a-1)= -1._dp/dyc(1:a-1)*(0.5_dp*(edF(1:a-1,i)+edF(2:a,i))*depthF(2:a,i)/dy(1:a-1)&
    +.5_dp*edF(2:a,i)*cbedS(2:a,i)*slopesS(2:a,i)) &
    + VDF(2:a,i)/(2._dp*dyc(1:a-1))
    lowera(2:a)= -1._dp/dyc(2:a)*(0.5_dp*(edF(1:a-1,i)+edF(2:a,i))*depthF(1:a-1,i)/dy(1:a-1)& 
    -.5_dp*edF(1:a-1,i)*cbedS(1:a-1,i)*slopesS(1:a-1,i)) &
    -VDF(1:a-1,i)/(2._dp*dyc(2:a))
    diaga(2:a-1)= depthF(2:a-1,i)/DT*2._dp +&
    1._dp/dyc(2:a-1)*(0.5_dp*( edF(3:a,i)+edF(2:a-1,i) )*depthF(2:a-1, i)/dy(2:a-1)&
    -.5_dp*edF(2:a-1,i)*cbedS(2:a-1,i)*slopesS(2:a-1,i))+&
    1._dp/dyc(1:a-2)*(0.5_dp*(edF(1:a-2,i)+edF(2:a-1,i))*depthF(2:a-1, i)/dy(1:a-2)&
    +.5_dp*edF(2:a-1,i)*cbedS(2:a-1,i)*slopesS(2:a-1,i)) + &
    wset*cbedS(2:a-1,i) !!

    diaga(1)= depthF(1,i)/DT*2._dp + &
    1._dp/dyc(1)*(0.5_dp*(edF(2,i)+edF(1,i))*depthF(1,i)/dy(1)- 0.5_dp*edF(1,i)*cbedS(1,i)*slopesS(1,i))
    diaga(a)= depthF(a,i)/DT*2._dp +&
    1._dp/dyc(a)*(0.5_dp*(edF(a,i)+edF(a-1,i))*depthF(a,i)/dy(a-1) + 0.5_dp*edF(a,i)*cbedS(a,i)*slopesS(a,i))

    !Right hand side, upwind method
    do j=1,a
    if(0.5_dp*(veldh(j,i)+veldh(j,i+1))>0._dp) THEN
            if(0.5_dp*(veldh(j,i)+veldh(j,i-1))>0._dp) THEN
            rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
            (2._dp*veldh(j,i+1)*sus_h(j,i+1)+3._dp*veldh(j,i)*sus_h(j,i) - 6._dp*veldh(j,i-1)*sus_h(j,i-1) &
            +  veldh(j,i-2)*sus_h(j,i-2))/(6._dp*delX) !!

    !rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
    !(veldh(j,i)*sus_h(j,i) - veldh(j,i-1)*sus_h(j,i-1))/delX !!
            else
            rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
            (3._dp*veldh(j,i+1)*sus_h(j,i+1)+0._dp*veldh(j,i)*sus_h(j,i) - 3._dp*veldh(j,i-1)*sus_h(j,i-1) &
            +  0._dp*veldh(j,i-2)*sus_h(j,i-2))/(6._dp*delX) !!
            end if
    ELSE
            if(0.5_dp*(veldh(j,i)+veldh(j,i-1))<0._dp) THEN
            rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
            (-veldh(j,i+2)*sus_h(j,i+2) +6._dp*veldh(j,i+1)*sus_h(j,i+1) -3._dp*veldh(j,i)*sus_h(j,i) &
             -2._dp*veldh(j,i-1)*sus_h(j,i-1) )/(6._dp*delX) !!

            !rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
            !(veldh(j,i+1)*sus_h(j,i+1) - veldh(j,i)*sus_h(j,i))/delX !!
            else
            rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
            (-veldh(j,i+2)*sus_h(j,i+2) +5._dp*veldh(j,i+1)*sus_h(j,i+1) -0._dp*veldh(j,i)*sus_h(j,i) &
             -5._dp*veldh(j,i-1)*sus_h(j,i-1) +veldh(j,i-2)*sus_h(j,i-2) )/(6._dp*delX) !!
            end if
    END IF
    END DO
    !For boundary conditions here, it is appropriate that the sus is zero on the
    !lateral boundaries-  which equates to us not having to enforce anything
    !further. Unless we have periodic boundaries, in which case the following might
    !be a practical alternative
    !if(periodicb) THEN
    !rhsa(1)=rhsa(1) -lowera(1)*sus_h(1,b)
    !rhsa(a)= rhsa(a)-uppera(a)*sus_h(a,b)
    !END IF
     

    !!Prevent main diagonal from being zero
    do ii=1,a
    if((diaga(ii)==0._dp).or.(depthF(ii,i)==0._dp)) THEN
            if( (ii==1).or.(ii==a).or.(depthF(ii,i)==0._dp)) THEN
            diaga(ii)=1._dp
            lowera(ii)=0._dp
            uppera(ii)=0._dp
            rhsa(ii)=0._dp
            else
            diaga(ii)=1._dp
            lowera(ii)=-0.5_dp
            uppera(ii)=-0.5_dp
            rhsa(ii)=0._dp
            end if
    end if
    end do

    call DGTSV(a,1, lowera(2:a), diaga, uppera(1:a-1), rhsa,a, info)

    Cdist(:,i)= rhsa

    if(info.ne.0) print*, "matrix problem in susconc2d, second bit", info


    do ii=1, a
    !if(sus_h(ii,i)<0._dp) THEN
    !        if(sus_h(ii,i)< - 1E-08_dp) print*, 'significant negative sus_h', j, ii, sus_h(ii,i), sus_h(a/2,i)
    !        !Cdist(j,ii)=0._dp
    !end if
    if(Cdist(ii,i)<0._dp) THEN
    !        if(Cdist(ii,i)< - 1E-08_dp) print*, 'significant negative Cdist', j, ii, Cdist(ii,i), Cdist(a/2,i), &
    !depthF(ii,i),depthH(ii,i)
            Cdist(ii,i)=0._dp
    end if
    end do

    !if(minval(Cdist(:,i))<0._dp) THEN
    !		!Cdist(:,i)=max(Cdist(:,i),0._dp)
    !	if(abs(minval(Cdist(:,i)))< 1.00E-12) THEN
    !
    !	ElSE
        !print*, "Cdist negative with non-roundoff type magnitude"
        !print*, "min Cdist", i, "<0", minval(Cdist(:,i)), minloc(Cdist(:,i)), depthF(minloc(Cdist(:,i)), (i-1):(i+1))& 
        !, vels(minloc(Cdist(:,i)), (i-1):(i+1)) 
        !END IF
    !
    !end if


    end do


    !if(minval(Cdist)<0._dp) print*,'Cdist<0 end', minloc(Cdist), minval(Cdist)

    Cdist_old=cbedS(1:a,:) !At the moment I am just using this to check the output. 
    !do i=1, b
    !if((Cdist(l(i),i)>3._dp*Cdist(l(i)+1,i)).and.(abs(Cdist(l(i),i))>1.0E-10_dp ) ) print*, 'Cdist spike end', i,l(i), & 
    !elevs(i)-hs(l(i)-1:l(i)+1,i), '##', Cdist(l(i):l(i)+1,i), Cdist(700,i)  
    !end do

    !print*, '2nd', Cdist(1,110), Cdist(350,110)

    !stop

end subroutine susconc2dup3rdV2
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Numerical method is upwind in the x derivative advective term, but still ADI
!otherwise
subroutine susconc2dup3rdV3(a,b, DT,delX, lengths,elevs,elevs_old,vels,vels_old,hs, Cdist,Cdist_old, Qe,Qe_old, Cmouth, Criver, & 
wset,fs, slopes, l, u,lambdacon)
    !use matrix_solvers
    !use crosssection

    !implicit none
    INTEGER, intent(in):: a, b, l, u
    REAL(dp), intent(in):: DT,delX, wset
    REAL(dp), intent(in):: elevs, elevs_old, vels, vels_old,hs, Qe,Qe_old, Cmouth, Criver,fs, lengths, slopes,lambdacon
    REAL(dp), intent(in out):: Cdist, Cdist_old !!Cdist_old will output random interesting information
    DIMENSION elevs(b), elevs_old(b), vels(a,b), vels_old(a,b), hs(a,b), Cdist(a,b), Cdist_old(a,b), Qe(a,b),fs(a,b), lengths(a,b) &
    ,Criver(a), Cmouth(a), Qe_old(a,b), slopes(a,b), l(b), u(b)

    INTEGER:: i, j, info, ii, KL=2, KU=2, IPV((-1):(b+2))
    logical:: periodicb=.false.!false.
    REAL(dp):: depthN(a,b), depthH(a,b), depthF(a,b), elevsH, velh(a,b),sus_h(a,(-1):(b+2)),edN(0:a+1,b),eddN(a,b), & 
            veldh(a,(-1):(b+2))
    REAL(dp):: sus_N(0:a+1,b), lnths(0:a+1,b), eddN_yh(0:a,b) ,edF(0:a+1,b), eddF(a,b), eddF_yh(0:a,b)
    REAL(dp):: upperb(b), diagb(b),lowerb(b), rhsb((-1):(b+2)), depthN2(0:a+1,b), VDN(0:a+1,b), VDF(a,b)
    REAL(dp):: uppera(a), diaga(a),lowera(a), rhsa(a), dyc(a), dy(0:a), sus_F(a,b), cbedS(0:a+1,b), slopesS(0:a+1,b)
    REAL(dp):: dedy, dedz=0.10_dp !Dimensionless eddy diffusivity
    REAL(dp):: tmpa(a), tmpb(b),tmp, band(7,b+4)

    !!!!So this solves the equation
    ! d/dt(depth*C) + d/dx(vel*depth*C) +d/dy(vel_lateral*depth*C)= d/dy(eddn*depth*dC/dy + eddn*Cbed*dh/dy)
    ! where I plan to add the last term soon.
    !By the ADI method, except that the convective x derivative is treated as upwind
    !(seems that it has to be to get stable results). 
    !!Note that assuming a constant vertical eddy diffusivity, 
    !Cbed= Caverage*depth/((ez/ws)*(1-exp(-vs/ez*(depth) ) ))
    !!And I have been assuming that ez=0.5_dp*ey

    !There are ghost points on the left and right edges (i.e. the margins of each
    !cross-section -- 1 each), and on the north
    !and south edges (i.e. the upstream and downstream boundaries -- 2 each)

    !Dimensionless eddy diffusivity
    dedy=lambdacon

    !!!First define water depths
    do i= 1,b
    !Depth at previous time step
    depthN(:,i)= max(elevs_old(i)-hs(:,i),0._dp)
    !Depth at half time step
    elevsH=(elevs(i)+elevs_old(i))*.5_dp
    depthH(:,i)= max(elevsH-hs(:,i), 0._dp)
    !Depth at next time step
    depthF(:,i)= max(elevs(i)-hs(:,i),0._dp)
    end do

    !!Useful thing
    depthN2(1:a,:)=depthN
    depthN2(0,:)=0._dp
    depthN2(a+1,:)=0._dp
    !!


    !Now calculate half time step velocity
    velH=.5_dp*(vels+vels_old)
    !Halfway depth times velocity - a useful efficiency device
    veldH(:,1:b)=depthH*velH
    veldH(:,0)= veldH(:,1) !B
    veldH(:,b+1)=veldH(:,b) !B
    veldH(:,-1)=veldH(:,0)
    veldH(:,b+2)=veldH(:,b) !B

    !Eddy diffusivity

    edN(1:a,:)= dedy*abs(vels_old)*sqrt(fs)*sqrt(0.125_dp)*depthN

    if(periodicb.eqv..false.) THEN
    edN(0,:)=0._dp
    edN(a+1,:)=0._dp

    !Make sure things cancel at the boundaries appropriately
    !do ii=1,b
    !edN(l(ii)-1, ii)=-edN(l(ii),ii)
    !edN(u(ii)+1,ii)=-edN(u(ii),ii)
    !end do
    !edN(0,:)=-edN(1,:)
    !edN(a+1,:)=-edN(a,:)
    ELSE
    edN(0,:)=edN(a,:)
    edN(a+1,:)=edN(1,:)
    end if

    !!Eddy diffusivity*depth at last time step
    eddN= edN(1:a,1:b)*depthN !Note that the sqrt(0.125_dp) is sqrt(1/8)



    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            
            if(edN(j,i)>0._dp) THEN
           tmp= ((edN(j,i)/dedy*dedz)*(1._dp-exp(-wset/(edN(j,i)/dedy*dedz)*depthN(j,i)) ) )
            ELSE
            tmp=0._dp
            END IF

        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthN(j,i)*wset/tmp, 200._dp)
        else
            cbedS(j,i)=200._dp  !So this is the maximum scaling factor we offer
            if((wset==0._dp).and.(edN(j,i)> 0._dp)) cbedS(j,i)=1._dp  !So the idea is that this applies  when wset=0
        
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp

    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times the depth, with boundaries
    eddN_yh(1:a-1,:)= .5_dp*(eddN(2:a,:)+eddN(1:a-1,:))

    if(periodicb.eqv..false.) THEN
    eddN_yh(0,:)=.5_dp*eddN(1,:) !B
    eddN_yh(a,:)=.5_dp*eddN(a,:) !B
    else
    eddN_yh(0,:)=.5_dp*(eddN(1,:)+eddN(a,:)) !B
    eddN_yh(a,:)=.5_dp*(eddN(a,:)+eddN(1,:)) !B
    end if

    !!Susconc at previous time step, with suitable boundary conditions attached
    sus_N(1:a,:)= Cdist_old
    !Zero susconc along the lateral boundaries -- might not be true - zero
    !integrated concentration, but the concentration might not drop off entirely
    sus_N(0,:)=sus_N(1,:) !0._dp !B
    sus_N(a+1,:)=sus_N(a,:) !0._dp !B


    !Lengths, with suitable boundaries attached
    lnths(1:a,:)= lengths
    lnths(0,:)= 2._dp*lengths(1,:)-lengths(2,:)!lengths(1,:)-depthN(1,:)/max(abs(slopes(1,:)), 0.01_dp)!2._dp*lengths(1,:)-lengths(2,:) !B
    lnths(a+1,:)=  2._dp*lengths(a,:)-lengths(a-1,:)!lengths(a,:) + depthN(a,:)/max(abs(slopes(a,:)), 0.01_dp)!2._dp*lengths(a,:)-lengths(a-1,:) !B

    slopesS(1:a,:)=slopes
    slopesS(0,:)=slopes(1,:) !0._dp
    slopesS(a+1,:)=slopes(a,:) !0._dp

    !Calculate VDN, which is the lateral discharge V*D at each point.
    !Based on the idea that
    ! VD = int( -dY/dt - dUD/dx ) dy
    !To ensure conservation, we calculate dY/dt to ensure that over the whole
    !cross-section, int(VD)=0
    VDN=0._dp 
    do i=2, b-1
    !Calculate the integral of dUd/dx. Note that we use u(i)+1 as a trick
    do j=max(l(i),2), min(u(i)+1,a)
    !do j=l(i), u(i)+1
    VDN(j,i) = VDN(j-1,i) + 0.5_dp*(-vels_old(j,i+1)*max(elevs_old(i+1)-hs(j,i+1),0._dp) - &
    vels_old(j-1,i+1)*max(elevs_old(i+1)-hs(j-1,i+1),0._dp) + & 
    vels_old(j,i-1)*max(elevs_old(i-1)-hs(j,i-1), 0._dp) + &
    vels_old(j-1,i-1)*max(elevs_old(i-1)-hs(j-1,i-1), 0._dp) )& 
    /(2._dp*delX)*(lengths(j,i)-lengths(j-1,i))
    end do
    !Now calculate dY/dT and correct - note the use of VDN(u(i)+1,i), which is later
    !set to zero
    VDN(l(i):u(i),i) = VDN(l(i):u(i),i) - &
    (lengths(l(i):u(i),i)-(lengths(max(l(i)-1,1),i) ))*&
    (VDN(min(u(i)+1,a),i)/((lengths(min(u(i)+1,a),i))-(lengths(max(l(i)-1,1),i) ) )) !Note that the last (VDN(u(i).. ) is dY/dt 

    VDN(min(u(i)+1,a),i)=0._dp
    !VDN(1:l(i)-1,i) = 999._dp
    !VDN(u(i)+1:a,i)=999._dp

    end do
    !print*, 'BUG'
    !VDN=0._dp
    !write(100,*) VDN(:,20) !, (-vels_old(:,21)*max(elevs_old(21)-hs(:,21),0._dp) + & 
    !vels_old(:,19)*max(elevs_old(19)-hs(:,19), 0._dp))

    !open(133,file='VDN')
    !write(133,*) VDN

    !!!!!!First half time step - implicit in the longitudinal coord - note how the lateral boundary condition is
    !automatically included - need to include the longitudinal boundary conditions
    !on the RHS
    do j=1,a
    band=0._dp
    !print*, j, size(Qe_old)
    !!Calculate matrix diagonals and RHS - upwind method
    DO i=1,b
    if(depthH(j,i)>0._dp) THEN
            if(0.5_dp*(veldH(j,i)+veldH(j,i+1))>0._dp) THEN
                    if(0.5_dp*(veldH(j,i)+veldH(j,i-1))>0._dp) THEN
                    !veldH(j,i+1/2)>0, veldH(j,i-1/2)>0
                    !Note the storage issues - because of the way I have indexed veldh(-1:b+2), I
                    !need to add +2 to each of the index's

                    !Upper-upper-diagonal
                    band(KL+KU+1+i+2-(i+2 +2),i+2+2)= 0._dp !Advection

                    !!Upper-diagonal
                    band(KL+KU+1+i+2-(i+1+2),i+1+2) = 1._dp/delX*1._dp/6._dp*2._dp*(veldH(j,i+1))  !Advection

                    !!Main diagonal
                    band(KL+KU+1+i+2-(i+2),i+2) = depthH(j,i)*2._dp/DT +wset*cbedS(j,i)  & !Unsteady
                    +1._dp/delX*1._dp/6._dp*3._dp*(veldH(j,i)) !Advection

                    !!Lower-diagonal
                    band(KL+KU+1+i+2-(i-1+2),i-1+2)= -1._dp/delX*1._dp/6._dp*6._dp*(veldH(j,i-1))  !Advection

                    !!Lower-lower-diagonal
                    band(KL+KU+1+i+2-(i-2+2),i-2+2)= 1._dp/delX*1._dp/6._dp*1._dp*(veldH(j,i-2)) !Advection

                    !        upperb(i)= 0._dp !veldH(j,2:b+1)*.5_dp/delX
                    !        diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) + veldh(j,i)/delX !!
                    !        lowerb(i)= -veldh(j,i-1)/delX !veldh(j,0:b-1)*.5_dp/delX
                    else
                    !veldH(j,i+1/2)>0, veldH(j,i-1/2)<=0
                    
                    band(KL+KU+1+i+2-(i+2 +2),i+2+2)= 0._dp !Advection

                    !!Upper-diagonal
                    band(KL+KU+1+i+2-(i+1+2),i+1+2) = 1._dp/delX*1._dp/6._dp*3._dp*(veldH(j,i+1))  !Advection

                    !!Main diagonal
                    band(KL+KU+1+i+2-(i+2),i+2) = depthH(j,i)*2._dp/DT +wset*cbedS(j,i)  & !Unsteady
                    +1._dp/delX*1._dp/6._dp*0._dp*(veldH(j,i)) !Advection

                    !!Lower-diagonal
                    band(KL+KU+1+i+2-(i-1+2),i-1+2)= -1._dp/delX*1._dp/6._dp*3._dp*(veldH(j,i-1))  !Advection

                    !!Lower-lower-diagonal
                    band(KL+KU+1+i+2-(i-2+2),i-2+2)= 1._dp/delX*1._dp/6._dp*0._dp*(veldH(j,i-2)) !Advection
                    
                    end if
            else !veldH(j,i+1/2)<=0._dp
                    
                    if(0.5_dp*(veldH(j,i)+veldH(j,i-1))<0._dp) THEN
                    !veldH(j,i+1/2)<=0, veldH(j,i-1/2)<0
                    !!Upper-upper diagonal
                    band(KL+KU+1+i+2-(i+2+2),i+2+2)=  -1._dp/delX*1._dp/6._dp*1._dp*(veldH(j,i+2)) !Advection

                    !!Upper-diagonal
                    band(KL+KU+1+i+2-(i+1+2),i+1+2) = 1._dp/delX*1._dp/6._dp*6._dp*(veldH(j,i+1))  !Advection

                    !!Main diagonal
                    band(KL+KU+1+i+2-(i+2),i+2) = depthH(j,i)*2._dp/DT +wset*cbedS(j,i) & !Unsteady
                    -1._dp/delX*1._dp/6._dp*3._dp*(veldH(j,i))  !Advection

                    !!Lower-diagonal
                    band(KL+KU+1+i+2-(i-1+2),i-1+2)= -1._dp/delX*1._dp/6._dp*2._dp*(veldH(j,i-1))  !Advection

                    !!Lower-lower-diagonal
                    band(KL+KU+1+i+2-(i-2+2),i-2+2)= 0._dp !Advection

                    !        upperb(i)= veldh(j,i+1)/delX !veldh(j,2:b+1)*.5_dp/delX
                    !        diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) - veldh(j,i)/delX !!
                    !        lowerb(i)= 0._dp !veldh(j,0:b-1)*.5_dp/delX
                    else
                    !veldH(j,i+1/2)<=0, veldH(j,i-1/2)>=0
                    band(KL+KU+1+i+2-(i+2+2),i+2+2)=  -1._dp/delX*1._dp/6._dp*1._dp*(veldH(j,i+2)) !Advection

                    !!Upper-diagonal
                    band(KL+KU+1+i+2-(i+1+2),i+1+2) = 1._dp/delX*1._dp/6._dp*5._dp*(veldH(j,i+1))  !Advection

                    !!Main diagonal
                    band(KL+KU+1+i+2-(i+2),i+2) = depthH(j,i)*2._dp/DT +wset*cbedS(j,i) & !Unsteady
                    -1._dp/delX*1._dp/6._dp*0._dp*(veldH(j,i))  !Advection

                    !!Lower-diagonal
                    band(KL+KU+1+i+2-(i-1+2),i-1+2)= -1._dp/delX*1._dp/6._dp*5._dp*(veldH(j,i-1))  !Advection

                    !!Lower-lower-diagonal
                    band(KL+KU+1+i+2-(i-2+2),i-2+2)= 1._dp/delX*1._dp/6._dp*1._dp*(veldH(j,i-2))  

                    !        upperb(i)= veldh(j,i+1)/delX !veldh(j,2:b+1)*.5_dp/delX
                    !        diagb(i)= depthH(j,i)*2._dp/DT +wset*cbedS(j,i) - veldh(j,i)/delX !!
                    !        lowerb(i)= 0._dp !veldh(j,0:b-1)*.5_dp/delX
                    end if
            end if
    ELSE !depth=0
    !Set sconc to 0
    band(KL+KU+1+i+2-(i+2),i+2) = 1._dp !depthN(j,i)*2._dp/dT !depthH(j,i)*2._dp/DT +wset*cbedS(j,i) & !Unsteady
    !-1._dp/delX*1._dp/6._dp*3._dp*(veldH(j,i))  !Advection
    END IF
    END DO

    rhsb=0._dp
    rhsb(1:b)= Qe_old(j,:)*sqrt(1._dp+slopes(j,:)**2)+ depthN(j,:)*sus_N(j,:)*2._dp/DT + 2._dp/(lnths(j+1,:)-lnths(j-1,:))*&  !!
    (0.5_dp*(edN(j,:)+edN(j+1,:))*(depthN2(j+1,:)*sus_N(j+1,:)-depthN2(j,:)*sus_N(j,:) )/(lnths(j+1,:)-lnths(j,:)) - & 
    0.5_dp*(edN(j,:)+edN(j-1,:))*(depthN2(j,:)*sus_N(j,:)-depthN2(j-1,:)*sus_N(j-1,:) )/(lnths(j,:)-lnths(j-1,:))  +&
    ((edN(j+1,:)+edN(j,:))*.5_dp*(cbedS(j+1,:)*sus_N(j+1,:)+cbedS(j,:)*sus_N(j,:))*.5_dp*(slopesS(j+1,:)+slopesS(j,:))*.5_dp - &
     (edN(j,:)+edN(j-1,:))*.5_dp*(cbedS(j,:)*sus_N(j,:)+cbedS(j-1,:)*sus_N(j-1,:))*.5_dp*(slopesS(j,:)+slopesS(j-1,:))*.5_dp) &
    !(.5_dp*(edN(j+1,:)*cbedS(j+1,:)*sus_N(j+1,:)*slopesS(j+1,:)+edN(j,:)*cbedS(j,:)*sus_N(j,:)*slopesS(j,:)) - &
    ! (.5_dp*(edN(j,:)*cbedS(j,:)*sus_N(j,:)*slopesS(j,:)+edN(j-1,:)*cbedS(j-1,:)*sus_N(j-1,:))*slopesS(j-1,:))) &
    ) - (VDN(j+1,:)*sus_N(j+1,:) - VDN(j-1,:)*sus_N(j-1,:))/(lnths(j+1,:)-lnths(j-1,:))

    do i=1,b
    if(depthH(j,i)==0._dp) rhsb(i)=0._dp
    end do

    !Boundary conditions at the mouth
    if(veldH(j,1)<0._dp) THEN
    !Here we say that sus_h(j,-1)=sus_h(j,0)=sus_h(j,1)
    i=-1
    band(KL+KU+1+i+2-(i+2+2),i+2+2) = -1._dp !Coefficient of sus_h(j,1) 
    band(KL+KU+1+i+2-(i+2),i+2) = 1._dp !Coefficient of sus_h(j,-1)

    i=0
    band(KL+KU+1+i+2-(i+1+2),i+1+2) = -1._dp !Coefficient of sus_h(j,1) 
    band(KL+KU+1+i+2-(i+2),i+2) = 1._dp !Coefficient of sus_h(j,0)
    else
    !Here we set the boundary points to Cmouth(j)
    i=-1
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp
    rhsb(i)=Cmouth(j)

    i=0
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp
    rhsb(i)=Cmouth(j)

    end if

    if(veldH(j,b)>0._dp) THEN
    !Here we set the boundary points sus_h(j,b+1), sus_h(j,b+2) == sus_h(j,b)
    i=b+1
    band(KL+KU+1+i+2-(i-1+2),i-1+2)=-1._dp !Coefficient for sus_h(j,b)
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp !Coefficient for sus_h(j,b+1)

    i=b+2
    band(KL+KU+1+i+2-(i-2+2),i-2+2)=-1._dp !Coefficient for sus_h(j,b)
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp !Coefficient for sus_h(j,b21)

    else
    !Here we set the boundaries to Criver
    i=b+1
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp !Coefficient for sus_h(j,b+1)
    rhsb(i)=Criver(j)

    i=b+2
    band(KL+KU+1+i+2-(i+2),i+2)=1._dp !Coefficient for sus_h(j,b+2)
    rhsb(i)=Criver(j)

    end if

    !if(j==a/2)THEN 
    !       IF(minval(sus_h(j-1,:))<0._dp) THEN
    !        open(131,file='banded')
    !        open(132,file='rhs')
    !        write(131,*) band
    !        close(131)
    !        write(132,*) rhsb
    !        close(132)
    !        stop
    !        end if
    !end if

    call DGBSV(b+4, 2, 2, 1, band, 2*KL+KU+1, IPV, rhsb((-1):(b+2)), b+4, info)
    !print*, 'info=', info
    !The solver writes the new C to r- so we change it back.
    !sus_h(j,1:b)=rhsb(1:b)
    sus_h(j,:)=rhsb(:)



    if(info.ne.0) print*, "matrix problem in susconc2d, first bit", info, j, depthH(j,info-2), depthH(j,info)
            
            !CHECKS
            do ii=1,b
            if(isnan(sus_h(j,ii))) THEN
            print*, 'sus_h(',j,ii,') is NAN' 
            stop
            end if
            
            !if(sus_h(j,ii)<0._dp) THEN
            !        if(sus_h(j,ii)< - 1E-04_dp) print*, 'significant negative sus_h', j, ii, sus_h(j,ii), sus_h(700,ii)
            !        sus_h(j,ii)=0._dp
            !From limited experimentation, it seems that even if this goes negative,
            !we don't end up with significant negative values in the final Cdist-
            !thus, maybe best to leave it alone.
            !end if
            
           ! if(depthH(j,ii)==0) THEN
           !         if(sus_h(j,ii)>0._dp) THEN
           !                if(sus_h(j,ii)>0.1_dp) THEN
           !                 print*, 'zero depth, high sus_h', sus_h(j,ii), j, ii,& 
     !depthH(j,ii)*2._dp/DT +wset*cbedS(j,ii) ! + veldh(j,ii)/delX, veldh(j,i+1)/delX, veldh(j,i-1)/delX 
           !                 stop
           !                 end if

           !                 sus_h(j,ii)=0._dp
           !        end if
           ! end if
            end do

    !        if((sus_h(j,b).ne.Criver(j)).and.(depthH(j,b)>0._dp)) THEN
    !               print*, 'sus_h(',j,"b) not equal to Criver", sus_h(j,b),Criver(j), depthH(j,b), depthN(j,b)&
    !,depthH(j,b)*2._dp/DT +wset*cbedS(j,b) - veldh(j,b)/delX,veldh(j,b)/delX, tmp, sus_N(j,b)
    !        end if


    end do
    !stop
    !print*, sus_h(1,110), sus_h(350,110)

    !print*, "half", sus_h(175,1), sus_h(175,40)
     
    !!!!!!!!!!!!!!End of first half time step


    !do i=1, b
    !if((sus_h(l(i),i)>3._dp*sus_h(l(i)+1,i)).and.(abs(sus_h(l(i),i))>1.0E-10_dp) ) print*, 'sus_h spike', i,l(i), & 
    !elevs_old(i)-hs(l(i)-1:l(i)+1,i),"##", sus_h(l(i):l(i)+1,i), sus_h(700,i)
    !end do


    !!!These boundary conditions are needed for the next half step
    sus_h(:,0)= Cmouth !B
    sus_h(:,b+1)=Criver !B
    sus_h(:,-1)= Cmouth !B
    sus_h(:,b+2)=Criver !B

    !!!Begin second half time step

    !!Eddy diffusivity
    edF(1:a,:)=dedy*abs(vels)*sqrt(fs)*sqrt(0.125_dp)*depthF

    if(periodicb.eqv..false.) THEN
    edF(0,:)=0._dp
    edF(a+1,:)=0._dp
    !Make things cancel at the boundary
    !do ii=1,b
    !edF(l(ii)-1,ii)=-edF(l(ii),ii)
    !edF(u(ii)+1,ii)=-edF(u(ii),ii)
    !end do

    !edF(0,:)=-edF(1,:) !0._dp
    !edF(a+1,:)=-edF(a,:)!0._dp
    ELSE
    edF(0,:)=edF(1,:)
    edF(a+1,:)=edF(a,:)
    END IF
    !!Eddy diffusivity*depth at last time step
    eddF= edF(1:a,:)*depthF !Note that the sqrt(0.125_dp) is sqrt(1/8) - note also how we ignore the change in the friction factor here

    !!A factor converting the bottom concentration to the depth averaged concentration
    ! Cbed=cbedS*Cave
    !cbedS= depthN*wset/((eddN/depthN)*(1._dp-exp(-wset/(eddN/depthN)*depthN) ) )
    !!Simplifying some of those depth terms, it simplifies to
    do i=1,b
    do j=1,a
            if(wset>0._dp) THEN
            tmp = ((edF(j,i)/dedy*dedz)*(1._dp-exp(-wset/(edF(j,i)/dedy*dedz)*depthF(j,i) ) ))
            else
            tmp=0._dp
            end if
        if(tmp.ne.0._dp) THEN
        cbedS(j,i)= min(depthF(j,i)*wset/tmp, 200._dp)
        else
        cbedS(j,i)=200._dp
        if((wset==0).and.(edF(j,i)> 0._dp)) cbedS(j,i)=1._dp
            end if
    end do
    end do
    cbedS(0,:)=cbedS(1,:) !0._dp
    cbedS(a+1,:)=cbedS(a,:) !0._dp

    !print*, cbedS(175,75), cbedS(175,75)*sus_h(175,75)
    !!!FROM BELOW, WE NOTE BOUNDARY CONDITIONS WITH A !B at the end

    !Halfway spatial average of eddy diffusivity times depth, with boundaries
    eddF_yh(1:a-1,:)= .5_dp*(eddF(2:a,:)+eddF(1:a-1,:))
    if(periodicb.eqv..false.) THEN
    eddF_yh(0,:)=.5_dp*eddF(1,:) !B
    eddF_yh(a,:)=.5_dp*eddF(a,:) !B
    else
    eddF_yh(0,:)=.5_dp*(eddF(1,:)+eddF(a,:)) !B
    eddF_yh(a,:)=.5_dp*(eddF(a,:)+eddF(1,:)) !B
    end if



    !Calculate VDF, which is the lateral discharge V*D at each point.
    !Based on the idea that
    ! VD = int( -dY/dt - dUD/dx ) dy
    !To ensure conservation, we calculate dY/dt to ensure that over the whole
    !cross-section, int(VD)=0
    VDF=0._dp 
    do i=2, b-1
    !Calculate the integral of dUd/dx. Note that we use u(i)+1 as a trick
    do j=max(l(i),2), min(u(i)+1,a)
    !do j=l(i), u(i)+1
    VDF(j,i) = VDF(j-1,i) + 0.5_dp*(-vels(j,i+1)*max(elevs(i+1)-hs(j,i+1),0._dp) - &
     vels(j-1,i+1)*max(elevs(i+1)-hs(j-1,i+1),0._dp) + & 
    vels(j,i-1)*max(elevs(i-1)-hs(j,i-1), 0._dp) + &
    vels(j-1,i-1)*max(elevs(i-1)-hs(j-1,i-1), 0._dp) )& 
    /(2._dp*delX)*(lengths(j,i)-lengths(j-1,i))
    end do
    !Now calculate dY/dT and correct - note the use of VDF(u(i)+1,i), which is later
    !set to zero
    VDF(l(i):u(i),i) = VDF(l(i):u(i),i) - &
    (lengths(l(i):u(i),i)-(lengths(max(l(i)-1,1),i)) )*&
    (VDF(min(u(i)+1,a),i)/((lengths(min(u(i)+1,a),i))-(lengths(max(l(i)-1,1),i) ) )) !Note that the last (VDN(u(i).. ) is dY/dt 

    VDF(min(u(i)+1,a),i)=0._dp
    !VDN(1:l(i)-1,i) = 999._dp
    !VDN(u(i)+1:a,i)=999._dp

    end do
    !print*, 'BUG'
    !VDF=0._dp



    do i=1,b

    !Define some dy type variables
    dyc(2:a-1)= .5_dp*(lengths(3:a,i)-lengths(1:a-2,i))
    dyc(1)=dyc(2) !0.5_dp*(lengths(2,i)-(lengths(1,i)-depthF(1,i)/max(abs(slopes(1,i)), 0.01_dp)))!dyc(2) !B
    dyc(a)=dyc(a-1)!0.5_dp*((lengths(a,i)+depthF(a,i)/max(abs(slopes(a,i)), 0.01_dp))-lengths(a-1,i))!dyc(a-1) !B

    dy(1:a-1)=lengths(2:a,i)-lengths(1:a-1,i)
    dy(a)=dy(a-1) !depthF(a,i)/max(abs(slopes(a,i)), 0.01_dp) !dy(a-1)
    dy(0)=dy(1) !depthF(1,i)/max(abs(slopes(1,i)), 0.01_dp) !dy(1)

    !!Matrix diagonals
    uppera(1:a-1)= -1._dp/dyc(1:a-1)*(0.5_dp*(edF(1:a-1,i)+edF(2:a,i))*depthF(2:a,i)/dy(1:a-1)&
    +.5_dp*edF(2:a,i)*cbedS(2:a,i)*slopesS(2:a,i)) &
    + VDF(2:a,i)/(2._dp*dyc(1:a-1))
    lowera(2:a)= -1._dp/dyc(2:a)*(0.5_dp*(edF(1:a-1,i)+edF(2:a,i))*depthF(1:a-1,i)/dy(1:a-1)& 
    -.5_dp*edF(1:a-1,i)*cbedS(1:a-1,i)*slopesS(1:a-1,i)) &
    -VDF(1:a-1,i)/(2._dp*dyc(2:a))
    diaga(2:a-1)= depthF(2:a-1,i)/DT*2._dp +&
    1._dp/dyc(2:a-1)*(0.5_dp*( edF(3:a,i)+edF(2:a-1,i) )*depthF(2:a-1, i)/dy(2:a-1)&
    -.5_dp*edF(2:a-1,i)*cbedS(2:a-1,i)*slopesS(2:a-1,i))+&
    1._dp/dyc(1:a-2)*(0.5_dp*(edF(1:a-2,i)+edF(2:a-1,i))*depthF(2:a-1, i)/dy(1:a-2)&
    +.5_dp*edF(2:a-1,i)*cbedS(2:a-1,i)*slopesS(2:a-1,i)) + &
    wset*cbedS(2:a-1,i) !!

    diaga(1)= depthF(1,i)/DT*2._dp + &
    1._dp/dyc(1)*(0.5_dp*(edF(2,i)+edF(1,i))*depthF(1,i)/dy(1)- 0.5_dp*edF(1,i)*cbedS(1,i)*slopesS(1,i))
    diaga(a)= depthF(a,i)/DT*2._dp +&
    1._dp/dyc(a)*(0.5_dp*(edF(a,i)+edF(a-1,i))*depthF(a,i)/dy(a-1) + 0.5_dp*edF(a,i)*cbedS(a,i)*slopesS(a,i))

    !Right hand side, upwind method
    do j=1,a
    if(0.5_dp*(veldh(j,i)+veldh(j,i+1))>0._dp) THEN
            if(0.5_dp*(veldh(j,i)+veldh(j,i-1))>0._dp) THEN
            rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
            (2._dp*veldh(j,i+1)*sus_h(j,i+1)+3._dp*veldh(j,i)*sus_h(j,i) - 6._dp*veldh(j,i-1)*sus_h(j,i-1) &
            +  veldh(j,i-2)*sus_h(j,i-2))/(6._dp*delX) !!

    !rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
    !(veldh(j,i)*sus_h(j,i) - veldh(j,i-1)*sus_h(j,i-1))/delX !!
            else
            rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) - & 
            (3._dp*veldh(j,i+1)*sus_h(j,i+1)+0._dp*veldh(j,i)*sus_h(j,i) - 3._dp*veldh(j,i-1)*sus_h(j,i-1) &
            +  0._dp*veldh(j,i-2)*sus_h(j,i-2))/(6._dp*delX) !!
            end if
    ELSE
            if(0.5_dp*(veldh(j,i)+veldh(j,i-1))<0._dp) THEN
            rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
            (-veldh(j,i+2)*sus_h(j,i+2) +6._dp*veldh(j,i+1)*sus_h(j,i+1) -3._dp*veldh(j,i)*sus_h(j,i) &
             -2._dp*veldh(j,i-1)*sus_h(j,i-1) )/(6._dp*delX) !!

            !rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
            !(veldh(j,i+1)*sus_h(j,i+1) - veldh(j,i)*sus_h(j,i))/delX !!
            else
            rhsa(j)= depthH(j,i)*sus_h(j,i)/DT*2._dp + Qe(j,i)*sqrt(1._dp+slopes(j,i)**2) -&
            (-veldh(j,i+2)*sus_h(j,i+2) +5._dp*veldh(j,i+1)*sus_h(j,i+1) -0._dp*veldh(j,i)*sus_h(j,i) &
             -5._dp*veldh(j,i-1)*sus_h(j,i-1) +veldh(j,i-2)*sus_h(j,i-2) )/(6._dp*delX) !!
            end if
    END IF
    END DO
    !For boundary conditions here, it is appropriate that the sus is zero on the
    !lateral boundaries-  which equates to us not having to enforce anything
    !further. Unless we have periodic boundaries, in which case the following might
    !be a practical alternative
    !if(periodicb) THEN
    !rhsa(1)=rhsa(1) -lowera(1)*sus_h(1,b)
    !rhsa(a)= rhsa(a)-uppera(a)*sus_h(a,b)
    !END IF
     

    !!Prevent main diagonal from being zero
    do ii=1,a
    if((diaga(ii)==0._dp).or.(depthF(ii,i)==0._dp)) THEN
            if( (ii==1).or.(ii==a).or.(depthF(ii,i)==0._dp)) THEN
            diaga(ii)=1._dp
            lowera(ii)=0._dp
            uppera(ii)=0._dp
            rhsa(ii)=0._dp
            else
            diaga(ii)=1._dp
            lowera(ii)=-0.5_dp
            uppera(ii)=-0.5_dp
            rhsa(ii)=0._dp
            end if
    end if
    end do

    call DGTSV(a,1, lowera(2:a), diaga, uppera(1:a-1), rhsa,a, info)

    Cdist(:,i)= rhsa

    if(info.ne.0) print*, "matrix problem in susconc2d, second bit", info


    do ii=1, a
    !if(sus_h(ii,i)<0._dp) THEN
    !        if(sus_h(ii,i)< - 1E-08_dp) print*, 'significant negative sus_h', j, ii, sus_h(ii,i), sus_h(a/2,i)
    !        !Cdist(j,ii)=0._dp
    !end if
    if(Cdist(ii,i)<0._dp) THEN
    !        if(Cdist(ii,i)< - 1E-08_dp) print*, 'significant negative Cdist', j, ii, Cdist(ii,i), Cdist(a/2,i), &
    !depthF(ii,i),depthH(ii,i)
            Cdist(ii,i)=0._dp
    end if
    end do

    !if(minval(Cdist(:,i))<0._dp) THEN
    !		!Cdist(:,i)=max(Cdist(:,i),0._dp)
    !	if(abs(minval(Cdist(:,i)))< 1.00E-12) THEN
    !
    !	ElSE
        !print*, "Cdist negative with non-roundoff type magnitude"
        !print*, "min Cdist", i, "<0", minval(Cdist(:,i)), minloc(Cdist(:,i)), depthF(minloc(Cdist(:,i)), (i-1):(i+1))& 
        !, vels(minloc(Cdist(:,i)), (i-1):(i+1)) 
        !END IF
    !
    !end if


    end do


    !if(minval(Cdist)<0._dp) print*,'Cdist<0 end', minloc(Cdist), minval(Cdist)

    Cdist_old=cbedS(1:a,:) !At the moment I am just using this to check the output. 
    !do i=1, b
    !if((Cdist(l(i),i)>3._dp*Cdist(l(i)+1,i)).and.(abs(Cdist(l(i),i))>1.0E-10_dp ) ) print*, 'Cdist spike end', i,l(i), & 
    !elevs(i)-hs(l(i)-1:l(i)+1,i), '##', Cdist(l(i):l(i)+1,i), Cdist(700,i)  
    !end do

    !print*, '2nd', Cdist(1,110), Cdist(350,110)

    !stop

end subroutine susconc2dup3rdV3
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!
end module sus



