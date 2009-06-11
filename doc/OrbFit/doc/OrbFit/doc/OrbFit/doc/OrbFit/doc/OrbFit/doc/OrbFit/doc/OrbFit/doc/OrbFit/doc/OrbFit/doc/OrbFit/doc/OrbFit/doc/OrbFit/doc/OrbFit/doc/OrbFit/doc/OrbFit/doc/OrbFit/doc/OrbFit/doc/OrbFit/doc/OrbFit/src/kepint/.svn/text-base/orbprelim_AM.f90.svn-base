
PROGRAM orbprelim_AM
!  USE fund_const
  USE orbit_elements
!  USE output_control
!  USE reference_systems, ONLY: roteceq
  USE attributable
  USE kepint, ONLY: twobodyint,nprelx
  IMPLICIT NONE
  INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision

!  REAL(KIND=8),DIMENSION(4) :: att1,att2
  TYPE(attrib) :: att1,att2
  INTEGER :: nprel ! number of preliminary orbits accepted
  TYPE(orbit_elem) :: orbprel(nprelx)

  CHARACTER*3 :: obscod1,obscod2 ! observatory code
  REAL(KIND=qkind),DIMENSION(2) :: alpha,delta,ad,dd ! attributables
  REAL(KIND=8) :: t1,t2 ! time of the attributables
  INTEGER :: i !loop indexes
! ==============OPTIONS============================
  CHARACTER(LEN=6) :: progna
  CHARACTER(LEN=80) :: run
! =========================================================
! options
  progna='compmo'
  run='ang_mom'
  CALL compop

  att1=undefined_attrib
  att2=undefined_attrib

  OPEN(1,file='p_coeff',status='unknown')    
  OPEN(3,file='q_coeff',status='unknown')    

  OPEN(2,file='attributables',status='old')  ! input file
  READ(2,*) att1%angles(1:4),att1%tdtobs,att1%obscod ! first attributable
  READ(2,*) att2%angles(1:4),att2%tdtobs,att2%obscod ! second attributable
 
! READ(2,*) alpha(1),delta(1),ad(1),dd(1),t1,obscod1 ! first attributable
! READ(2,*) alpha(2),delta(2),ad(2),dd(2),t2,obscod2 ! second attributable
! write(*,*)alpha(1),delta(1),ad(1),dd(1),t1
! write(*,*)alpha(2),delta(2),ad(2),dd(2),t2
! att1(1)=alpha(1)
! att1(2)=delta(1)
! att1(3)=ad(1)
! att1(4)=dd(1)
! att2(1)=alpha(2)
! att2(2)=delta(2)
! att2(3)=ad(2)
! att2(4)=dd(2)
!  CALL twobodyint(att1,t1,obscod1,att2,t2,obscod2)

  CALL twobodyint(att1,att2,nprel,orbprel,CHECKDER=.false.)
  CLOSE(3)
  CLOSE(2)
  CLOSE(1)

  DO i=1,nprel
     write(*,*) '---------------------------------------------'
     write(*,*) orbprel(i)%coord(1:6)
  ENDDO

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
END PROGRAM orbprelim_AM

