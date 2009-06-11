PROGRAM link_radar_obs
  USE fund_const
  USE orbit_elements
  USE output_control
  USE reference_systems, ONLY: roteceq
  IMPLICIT NONE
! orbital elements
  DOUBLE PRECISION,DIMENSION(6) :: eqtmp1,eqtmp2 ! auxiliary
  TYPE(orbit_elem) :: elem1,elem2 ! EQU elements of the Earth
  TYPE(orbit_elem) :: ecpl1,ecpl2 ! CAR elements of the Earth
  TYPE(orbit_elem) :: elkep   ! Keplerian elements of the asteroid
! ATT orbital elements
  TYPE(orbit_elem) :: attel1t1,attel1t2 ! 1st choice for ATT elems at t1,t2
  TYPE(orbit_elem) :: kepel1t1,kepel1t2 ! 1st choice for KEP elems at t1,t2
  TYPE(orbit_elem) :: attel2t1,attel2t2 ! 2nd choice for ATT elems at t1,t2
  TYPE(orbit_elem) :: kepel2t1,kepel2t2 ! 2nd choice for KEP elems at t1,t2
  DOUBLE PRECISION,DIMENSION(2) :: alpha,delta,rho,rhodot ! radar observations
  DOUBLE PRECISION,DIMENSION(2,2) :: alphadot,deltadot ! unknowns
  DOUBLE PRECISION :: t1,t2 ! time of the observations
!
  DOUBLE PRECISION,DIMENSION(3,2) :: rhohat,rhohat_a,rhohat_d
  DOUBLE PRECISION :: detM,detJ1,detK1,detJ2,detK2,detJ3,detK3
  DOUBLE PRECISION,DIMENSION(3,2) :: pe,pde ! pos and vel of the Earth
  DOUBLE PRECISION,DIMENSION(2) :: c0,c1,c2,c3,c31,c32,c4,c4tilde,c5
  DOUBLE PRECISION :: E1,F1,G1,H1,I1,E2,F2,G2,H2,I2
  DOUBLE PRECISION :: AA,BB,CC,DD,EE,FF
  DOUBLE PRECISION :: calA,calB,calC,discrim ! coeffs of the 2nd degree eq.
! auxiliary --------------------------------------------------------
  DOUBLE PRECISION,DIMENSION(3,2) :: pe_pde,pe_rhat
  DOUBLE PRECISION,DIMENSION(3,2) :: pe_rhata,pe_rhatd,rhat_rhata,rhat_rhatd
  DOUBLE PRECISION,DIMENSION(3,2) :: rhat_pde
  DOUBLE PRECISION,DIMENSION(3) :: coef
  DOUBLE PRECISION :: q,ecc,inc,omnod,omeg,tperi
!
  INTEGER :: i,j ! loop indexes
  INTEGER :: fail_flag
! ==============OPTIONS============================
  CHARACTER(LEN=6) :: progna
  CHARACTER(LEN=80) :: run
! =========================================================

! options
  progna='compmo'
  run='ang_mom'
  CALL compop

  OPEN(1,file='orbital_elements',status='unknown') ! output file
  OPEN(2,file='radar_obs',status='old') ! input file
  READ(2,*) alpha(1),delta(1),rho(1),rhodot(1),t1 ! first attributable
  READ(2,*) alpha(2),delta(2),rho(2),rhodot(2),t2 ! second attributable
  write(*,*)alpha(1),delta(1),rho(1),rhodot(1),t1

  DO j =1,2
     rhohat(1,j) = cos(alpha(j))*cos(delta(j))
     rhohat(2,j) = sin(alpha(j))*cos(delta(j))
     rhohat(3,j) = sin(delta(j))
     rhohat_a(1,j) = -sin(alpha(j))*cos(delta(j))
     rhohat_a(2,j) = cos(alpha(j))*cos(delta(j))
     rhohat_a(3,j) = 0.d0
     rhohat_d(1,j) = -cos(alpha(j))*sin(delta(j))
     rhohat_d(2,j) = -sin(alpha(j))*sin(delta(j))
     rhohat_d(3,j) = cos(delta(j))
  ENDDO

! Earth elements at time t1
  elem1=undefined_orbit_elem
  CALL earth(t1,eqtmp1)
  elem1%coo='EQU'
  elem1%t=t1
  elem1%coord(1:6)=eqtmp1
  call coo_cha(elem1,'CAR',ecpl1,fail_flag)
  pe(1:3,1) = MATMUL(roteceq,ecpl1%coord(1:3))
  pde(1:3,1) = MATMUL(roteceq,ecpl1%coord(4:6))
!  write(*,*)'Earth position at t1:',  pe(1:3,1)

! Earth elements at time t2
  elem2=undefined_orbit_elem
  CALL earth(t2,eqtmp2)
  elem2%coo='EQU'
  elem2%t=t2
  elem2%coord(1:6)=eqtmp2
  call coo_cha(elem2,'CAR',ecpl2,fail_flag)
  pe(1:3,2) = MATMUL(roteceq,ecpl2%coord(1:3))
  pde(1:3,2) = MATMUL(roteceq,ecpl2%coord(4:6))
!  write(*,*)'Earth position at t2:',  pe(1:3,2)

!  DO j=1,2
!     CALL cross_prod(pe(1:3,j),pde(1:3,j),pe_pde(1:3,j))
!     CALL cross_prod(pe(1:3,j),rhohat_a(1:3,j),pe_rhata(1:3,j))
!     CALL cross_prod(pe(1:3,j),rhohat_d(1:3,j),pe_rhatd(1:3,j))
!     CALL cross_prod(rhohat(1:3,j),pde(1:3,j),rhat_pde(1:3,j))
!     CALL cross_prod(rhohat(1:3,j),rhohat_a(1:3,j),rhat_rhata(1:3,j))
!     CALL cross_prod(rhohat(1:3,j),rhohat_d(1:3,j),rhat_rhatd(1:3,j))
!     CALL cross_prod(pe(1:3,j),rhohat(1:3,j),pe_rhat(1:3,j))
!  ENDDO

! ******* ADMISSIBLE REGION CURVE ********
  DO j = 1,2
     c0(j) = DOT_PRODUCT(pe(1:3,j),pe(1:3,j))
     c1(j) = 2.d0*DOT_PRODUCT(pde(1:3,j),rhohat(1:3,j))
!     c2(j) = (ad(j)*cos(delta(j)))**2 + dd(j)**2
     c31(j) = 2.d0*DOT_PRODUCT(pde(1:3,j),rhohat_a(1:3,j))
     c32(j) = 2.d0*DOT_PRODUCT(pde(1:3,j),rhohat_d(1:3,j))
!     c3(j) = ad(j)*cc31(j) + dd(j)*cc32(j);
     c4(j) = DOT_PRODUCT(pde(1:3,j),pde(1:3,j))
     c5(j) = 2.d0*DOT_PRODUCT(pe(1:3,j),rhohat(1:3,j))
! this is used for the admissible region, not for the def. of energy
!     c4tilde(j) = c4(j) + gms/(4*100)
!     c4(j) = c4tilde(j)
  ENDDO

! -----------------------------------------------------------------
! coefficient of the equal ang. momentum curve
  E1 = cos(delta(1))**2*rho(1)**2
  F1 = rho(1)**2
  G1 = c31(1)*rho(1)
  H1 = c32(1)*rho(1)
  I1 = rhodot(1)**2 + c1(1)*rhodot(1) + c4(1) - &
       & 2.d0*gms/sqrt(rho(1)**2+c5(1)*rho(1)+c0(1))

  E2 = cos(delta(2))**2*rho(2)**2
  F2 = rho(2)**2
  G2 = c31(2)*rho(2)
  H2 = c32(2)*rho(2)
  I2 = rhodot(2)**2 + c1(2)*rhodot(2) + c4(2) - &
       & 2.d0*gms/sqrt(rho(2)**2+c5(2)*rho(2)+c0(2))

! eliminate ad1,dd1,ad2 by Cramer's rule
  CALL determinants(pe, pde, rho, rhohat, rhohat_a, rhohat_d, &
       & rhodot, detM, detJ1, detK1, detJ2, detK2, detJ3, detK3)
  write(*,*)'det(M) =',detM

! ad1 = AA*dd2 + BB
  AA = detK1/detM
  BB = detJ1/detM
! dd1 = CC*dd2 + DD
  CC = detK2/detM
  DD = detJ2/detM
! ad2 = EE*dd2 + FF
  EE = detK3/detM
  FF = detJ3/detM

! coefficients of the 2nd degree eq. for deltadot(2)
  calA = E1*AA**2+CC**2*F1-EE**2*E2-F2
  calB = 2*E1*AA*BB+2*DD*CC*F1+G1*AA+CC*H1-2*FF*EE*E2-EE*G2-H2
  calC = E1*BB**2+F1*DD**2+G1*BB+H1*DD+I1-E2*FF**2-G2*FF-I2

! --- compute attributable orbital elements ---
! possible values of deltadot(2)
  discrim=calB**2 - 4.d0*calA*calC
  IF(discrim.ge.0-epsilon(1.d0)*1.d3) THEN
     deltadot(2,1) = (-calB + sqrt(discrim))/(2.d0*calA) 
     deltadot(2,2) = (-calB - sqrt(discrim))/(2.d0*calA) 

     alphadot(1,1:2) = AA*deltadot(2,1:2)+BB
     deltadot(1,1:2) = CC*deltadot(2,1:2)+DD
     alphadot(2,1:2) = EE*deltadot(2,1:2)+FF
  ELSE
     WRITE(*,*)'negative discriminant! stopping computation'
     STOP
  ENDIF

! write possible orbits at time t_1
  write(*,*)'----------------------------------------'
  write(*,*)'t=',t1
  WRITE(*,*) alphadot(1,1),deltadot(1,1)
  WRITE(*,*) alphadot(1,2),deltadot(1,2)
  write(*,*)'t=',t2
  WRITE(*,*) alphadot(2,1),deltadot(2,1)
  WRITE(*,*) alphadot(2,2),deltadot(2,2)

  WRITE(1,100) t1,alpha(1),delta(1),alphadot(1,1), &
       & deltadot(1,1),rho(1),rhodot(1)
  WRITE(1,100) t1,alpha(1),delta(1),alphadot(1,2), &
       & deltadot(1,2),rho(1),rhodot(1)

100 FORMAT(6(e14.7,2x))

! check with omega and elle
! .... STILL TO BE WRITTEN ...

! ATT elements at time t1
  attel1t1=undefined_orbit_elem
  attel1t1%coo='ATT'
  attel1t1%t=t1
  attel1t1%coord(1)= alpha(1)
  attel1t1%coord(2)= delta(1)
  attel1t1%coord(3)= alphadot(1,1)
  attel1t1%coord(4)= deltadot(1,1)
  attel1t1%coord(5)= rho(1)
  attel1t1%coord(6)= rhodot(1)
  call coo_cha(attel1t1,'KEP',kepel1t1,fail_flag)

  attel2t1=undefined_orbit_elem
  attel2t1%coo='ATT'
  attel2t1%t=t1
  attel2t1%coord(1)= alpha(1)
  attel2t1%coord(2)= delta(1)
  attel2t1%coord(3)= alphadot(1,2)
  attel2t1%coord(4)= deltadot(1,2)
  attel2t1%coord(5)= rho(1)
  attel2t1%coord(6)= rhodot(1)
  call coo_cha(attel2t1,'KEP',kepel2t1,fail_flag)

! ATT elements at time t2
  attel1t2=undefined_orbit_elem
  attel1t2%coo='ATT'
  attel1t2%t=t2
  attel1t2%coord(1)= alpha(2)
  attel1t2%coord(2)= delta(2)
  attel1t2%coord(3)= alphadot(2,1)
  attel1t2%coord(4)= deltadot(2,1)
  attel1t2%coord(5)= rho(2)
  attel1t2%coord(6)= rhodot(2)
  call coo_cha(attel1t2,'KEP',kepel1t2,fail_flag)

  attel2t2=undefined_orbit_elem
  attel2t2%coo='ATT'
  attel2t2%t=t2
  attel2t2%coord(1)= alpha(2)
  attel2t2%coord(2)= delta(2)
  attel2t2%coord(3)= alphadot(2,2)
  attel2t2%coord(4)= deltadot(2,2)
  attel2t2%coord(5)= rho(2)
  attel2t2%coord(6)= rhodot(2)
  call coo_cha(attel2t2,'KEP',kepel2t2,fail_flag)

  write(*,*)'------------- KEPLERIAN ELEMENTS -------------'
  WRITE(*,101)kepel1t1%coord
  WRITE(*,101)kepel2t1%coord
  WRITE(*,101)kepel1t2%coord
  WRITE(*,101)kepel2t2%coord
101 FORMAT(6(F12.5,2X))
  CLOSE(2)
  CLOSE(1)

CONTAINS
! =========================================================
  SUBROUTINE compop
    IMPLICIT NONE
    INTEGER            :: le,iunout
    CHARACTER(LEN=100) :: file
    LOGICAL            :: ireq,found
    CHARACTER(LEN=60)  :: comment
! =============================
    CALL initopt(progna,run,'mop')
! read option for physical model and integration method
    CALL rmodel
  ! Output files: for control
    file=run//'.mou'
    CALL rmsp(file,le)
    call filopn(iunout,file(1:le),'UNKNOWN')
! =============WHERE TO FIND ASTEROID ELEMENTS===================
!    ireq=.true.
!    comment='asteroid elements directory'
!    CALL input_cha_opt(progna,'eledir',eledir,ireq,found,comment,iunout)
  END SUBROUTINE compop
END PROGRAM link_radar_obs

SUBROUTINE cross_prod(v1,v2,v3)
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN),DIMENSION(3) :: v1,v2
  DOUBLE PRECISION,INTENT(OUT),DIMENSION(3) :: v3
  v3(1) = v1(2)*v2(3)-v1(3)*v2(2)
  v3(2) = v1(3)*v2(1)-v1(1)*v2(3)
  v3(3) = v1(1)*v2(2)-v1(2)*v2(1)
END SUBROUTINE cross_prod
