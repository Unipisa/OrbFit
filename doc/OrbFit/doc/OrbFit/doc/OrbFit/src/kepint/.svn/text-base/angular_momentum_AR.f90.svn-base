PROGRAM angular_momentum_AR
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

  DOUBLE PRECISION,DIMENSION(2) :: alpha,delta,ad,dd ! attributables
  DOUBLE PRECISION :: t1,t2 ! time of the attributables
  DOUBLE PRECISION,DIMENSION(3,2) :: rhohat,rhohat_a,rhohat_d
  DOUBLE PRECISION,DIMENSION(3,2) :: pe,pde
  DOUBLE PRECISION :: co00,co01,co02,co10,co11,co20,co21,co30,co40
  DOUBLE PRECISION :: C1,C2,C3,C4,C5,C6
  DOUBLE PRECISION :: D1,D2,D3,D4,D5,D6
!  DOUBLE PRECISION :: rd1 ! radial velocity
  DOUBLE PRECISION :: cc0,cc1,cc2,cc3,cc31,cc32,cc4,cc5
! auxiliary --------------------------------------------------------
  DOUBLE PRECISION,DIMENSION(3,2) :: pe_pde,pe_rhat
  DOUBLE PRECISION,DIMENSION(3,2) :: pe_rhata,pe_rhatd,rhat_rhata,rhat_rhatd
  DOUBLE PRECISION,DIMENSION(3,2) :: rhat_pde
  DOUBLE PRECISION,DIMENSION(3) :: coef
  DOUBLE PRECISION :: aux,fact1,fact2
  DOUBLE PRECISION :: q,ecc,inc,omnod,omeg,tperi
  INTEGER :: i1,i2
!
  INTEGER :: isel ! selected index
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

  OPEN(1,file='coeff',status='unknown') ! output file

  OPEN(2,file='attributables',status='old') ! input file
  READ(2,*) alpha(1),delta(1),ad(1),dd(1),t1 ! first attributable
  READ(2,*) alpha(2),delta(2),ad(2),dd(2),t2 ! second attributable
 ! write(*,*)alpha(1),delta(1),ad(1),dd(1),t1

  OPEN(3,file='admis_reg',status='unknown') 

  rhohat(1,1) = cos(alpha(1))*cos(delta(1))
  rhohat(2,1) = sin(alpha(1))*cos(delta(1))
  rhohat(3,1) = sin(delta(1))
  rhohat_a(1,1) = -sin(alpha(1))*cos(delta(1))
  rhohat_a(2,1) = cos(alpha(1))*cos(delta(1))
  rhohat_a(3,1) = 0.d0
  rhohat_d(1,1) = -cos(alpha(1))*sin(delta(1))
  rhohat_d(2,1) = -sin(alpha(1))*sin(delta(1))
  rhohat_d(3,1) = cos(delta(1))

  rhohat(1,2) = cos(alpha(2))*cos(delta(2))
  rhohat(2,2) = sin(alpha(2))*cos(delta(2))
  rhohat(3,2) = sin(delta(2))
  rhohat_a(1,2) = -sin(alpha(2))*cos(delta(2))
  rhohat_a(2,2) = cos(alpha(2))*cos(delta(2))
  rhohat_a(3,2) = 0.d0
  rhohat_d(1,2) = -cos(alpha(2))*sin(delta(2))
  rhohat_d(2,2) = -sin(alpha(2))*sin(delta(2))
  rhohat_d(3,2) = cos(delta(2))

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

  DO j=1,2
     CALL cross_prod(pe(1:3,j),pde(1:3,j),pe_pde(1:3,j))
     CALL cross_prod(pe(1:3,j),rhohat_a(1:3,j),pe_rhata(1:3,j))
     CALL cross_prod(pe(1:3,j),rhohat_d(1:3,j),pe_rhatd(1:3,j))
     CALL cross_prod(rhohat(1:3,j),pde(1:3,j),rhat_pde(1:3,j))
     CALL cross_prod(rhohat(1:3,j),rhohat_a(1:3,j),rhat_rhata(1:3,j))
     CALL cross_prod(rhohat(1:3,j),rhohat_d(1:3,j),rhat_rhatd(1:3,j))
     CALL cross_prod(pe(1:3,j),rhohat(1:3,j),pe_rhat(1:3,j))
  ENDDO

! ******* ADMISSIBLE REGION CURVE ********
  cc0 = DOT_PRODUCT(pe(1:3,1),pe(1:3,1))
  cc1 = 2.d0*DOT_PRODUCT(pde(1:3,1),rhohat(1:3,1))
  cc2 = (ad(1)*cos(delta(1)))**2 + dd(1)**2
  cc31 = 2.d0*DOT_PRODUCT(pde(1:3,1),rhohat_a(1:3,1))
  cc32 = 2.d0*DOT_PRODUCT(pde(1:3,1),rhohat_d(1:3,1))
  cc3 = ad(1)*cc31 + dd(1)*cc32;
  cc4 = DOT_PRODUCT(pde(1:3,1),pde(1:3,1))
  cc5 = 2.d0*DOT_PRODUCT(pde(1:3,1),rhohat(1:3,1))
!  cc4tilde = cc4 + gms/(4*100);
!  cc4 = c4tilde;
  WRITE(3,112) cc0,cc1,cc2,cc3,cc4,cc5,gms
112 FORMAT(7(e14.7,1x))

! ******* ANGULAR MOMENTUM CURVE *******  
  IF(DOT_PRODUCT(pe_rhat(1:3,2),pe_rhat(1:3,2)).gt.epsilon(1.d0)*1000)THEN
     ! can perform elimination of rd2     
     aux=0.d0
     DO i =1,3
        coef(i)=abs(pe_rhat(i,2))
        IF(coef(i).gt.aux) THEN
           isel=i
           aux=coef(i)
        ENDIF
     ENDDO
     IF(isel.eq.1)THEN
        i1=2;i2=3
     ELSEIF(isel.eq.2)THEN
        i1=1;i2=3
     ELSEIF(isel.eq.3)THEN
        i1=1;i2=2
     ENDIF
  ELSE
     WRITE(*,*)'ERROR! cannot eliminate rd2!'
     WRITE(*,*)'STOPPING PROGRAM'
     STOP
  ENDIF

  fact1 = pe_rhat(i1,2)/pe_rhat(isel,2)
  C1 =(pe_pde(i1,1)-pe_pde(i1,2))- (pe_pde(isel,1)-pe_pde(isel,2))*fact1
  C2 =ad(1)*pe_rhata(i1,1)+dd(1)*pe_rhatd(i1,1)+rhat_pde(i1,1) - &
    & (ad(1)*pe_rhata(isel,1)+dd(1)*pe_rhatd(isel,1)+rhat_pde(isel,1))*fact1
  C3 = ad(1)*rhat_rhata(i1,1)+dd(1)*rhat_rhatd(i1,1) - &
       & (ad(1)*rhat_rhata(isel,1)+dd(1)*rhat_rhatd(isel,1))*fact1
  C4 = pe_rhat(i1,1)-pe_rhat(isel,1)*fact1
  C5 =-(ad(2)*pe_rhata(i1,2)+dd(2)*pe_rhatd(i1,2)+rhat_pde(i1,2)) + &
    & (ad(2)*pe_rhata(isel,2)+dd(2)*pe_rhatd(isel,2)+rhat_pde(isel,2))*fact1
  C6 =-(ad(2)*rhat_rhata(i1,2)+dd(2)*rhat_rhatd(i1,2)) + &
       & (ad(2)*rhat_rhata(isel,2)+dd(2)*rhat_rhatd(isel,2))*fact1

  fact2 = pe_rhat(i2,2)/pe_rhat(isel,2)
  D1 =(pe_pde(i2,1)-pe_pde(i2,2))- (pe_pde(isel,1)-pe_pde(isel,2))*fact2
  D2 =ad(1)*pe_rhata(i2,1)+dd(1)*pe_rhatd(i2,1)+rhat_pde(i2,1) - &
    & (ad(1)*pe_rhata(isel,1)+dd(1)*pe_rhatd(isel,1)+rhat_pde(isel,1))*fact2
  D3 = ad(1)*rhat_rhata(i2,1)+dd(1)*rhat_rhatd(i2,1) - &
       & (ad(1)*rhat_rhata(isel,1)+dd(1)*rhat_rhatd(isel,1))*fact2
  D4= pe_rhat(i2,1)-pe_rhat(isel,1)*fact2
  D5 =-(ad(2)*pe_rhata(i2,2)+dd(2)*pe_rhatd(i2,2)+rhat_pde(i2,2)) + &
    & (ad(2)*pe_rhata(isel,2)+dd(2)*pe_rhatd(isel,2)+rhat_pde(isel,2))*fact2
  D6 =-(ad(2)*rhat_rhata(i2,2)+dd(2)*rhat_rhatd(i2,2)) + &
       & (ad(2)*rhat_rhata(isel,2)+dd(2)*rhat_rhatd(isel,2))*fact2

! -----------------------------------------------------------------
! coefficient of the equal ang. momentum curve
  co00 = D6**2*C1**2+C6**2*D1**2+D1*D6*C5**2-D5*C5*C6*D1-2.D0*D6*C1*&
     &C6*D1-D5*D6*C5*C1+C6*D5**2*C1                                     
                                                                        
  co01 = -2.D0*D6*C1*C6*D4-C4*C5*D6*D5-D4*C6*C5*D5+C4*D5**2*C6-2.D0*&
     &D6*C4*C6*D1+2.D0*D1*D4*C6**2+2.D0*C4*C1*D6**2+C5**2*D6*D4         
                                                                        
  co02 = (-C6*D4+D6*C4)**2.D0 
                                                                        
  co10 = -2.D0*D6*C1*C6*D2-2.D0*D6*C2*C6*D1+2.D0*D2*D1*C6**2-D5*D6*C&
     &5*C2-D5*C5*C6*D2+C6*D5**2*C2+D2*D6*C5**2+2.D0*D6**2*C1*C2         
                                                                        
  co11 = 2.D0*(-C6*D4+D6*C4)*(D6*C2-C6*D2) 
                                                                        
  co20 = -2.D0*D6*C1*C6*D3-2.D0*D6*C2*C6*D2-2.D0*D6*C3*C6*D1+D6**2*C&
     &2**2+2.D0*D3*D1*C6**2-D5*D6*C5*C3-D5*C5*C6*D3+C6*D5**2*C3+D3*D6*C5&
     &**2+2.D0*D6**2*C1*C3+D2**2*C6**2                                  
                                                                        
  co21 = 2.D0*(-C6*D4+D6*C4)*(-C6*D3+D6*C3) 
                                                                        
  co30 = 2.D0*(-C6*D3+D6*C3)*(D6*C2-C6*D2) 
                                                                        
  co40 = (-C6*D3+D6*C3)**2.D0 
                                                                        

!  co0 = D6**2*C1**2+2.D0*D6**2*C1*C4*rd1-2.D0*D6*C1*C6*D1+C6**2*D4**&
!       &2*rd1**2-D5*D6*C5*C1+C6*D5**2*C1+D4*rd1*D6*C5**2+D6**2*C4**2*rd1**&
!       &2+C6**2*D1**2+D1*D6*C5**2-2.D0*D6*C1*C6*D4*rd1-2.D0*D6*C4*rd1**2*C&
!       &6*D4-2.D0*D6*C4*rd1*C6*D1+2.D0*C6**2*D4*rd1*D1-D5*D6*C5*C4*rd1-D5*&
!       &C5*C6*D4*rd1-D5*C5*C6*D1+C6*D5**2*C4*rd1                          
!  co1 = 2.D0*D6**2*C1*C2+D2*D6*C5**2-2.D0*D6*C2*C6*D4*rd1-2.D0*D6*C1&
!       &*C6*D2+2.D0*D6**2*C2*C4*rd1-2.D0*D6*C2*C6*D1-2.D0*D6*C4*rd1*C6*D2+&
!       &2.D0*C6**2*D4*rd1*D2+2.D0*C6**2*D1*D2-D5*D6*C5*C2-D5*C5*C6*D2+C6*D&
!       &5**2*C2                                                           
!  co2 = 2.D0*D6**2*C1*C3+C6**2*D2**2+D3*D6*C5**2+D6**2*C2**2-2.D0*D6&
!       &*C1*C6*D3-2.D0*D6*C2*C6*D2+2.D0*D6**2*C3*C4*rd1-2.D0*D6*C3*C6*D4*r&
!       &d1-2.D0*D6*C3*C6*D1-2.D0*D6*C4*rd1*C6*D3+2.D0*C6**2*D4*rd1*D3+2.D0&
!       &*C6**2*D1*D3-D5*D6*C5*C3-D5*C5*C6*D3+C6*D5**2*C3                  
!  co3 = 2.D0*(-D6*C3+C6*D3)*(D2*C6-D6*C2)  
!  co4 = (-D6*C3+C6*D3)**2.D0 

  WRITE(1,100) co00,co01,co02,co10,co11,co20,co21,co30,co40
100 FORMAT(9(e14.7,2x))

  CLOSE(3)
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
END PROGRAM angular_momentum_AR

SUBROUTINE cross_prod(v1,v2,v3)
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN),DIMENSION(3) :: v1,v2
  DOUBLE PRECISION,INTENT(OUT),DIMENSION(3) :: v3
  v3(1) = v1(2)*v2(3)-v1(3)*v2(2)
  v3(2) = v1(3)*v2(1)-v1(1)*v2(3)
  v3(3) = v1(1)*v2(2)-v1(2)*v2(1)
END SUBROUTINE cross_prod

