MODULE angular_momentum
!  USE fund_const
!  USE orbit_elements
!  USE output_control
  IMPLICIT NONE
  PRIVATE

! public routines
  PUBLIC :: compsylv48
CONTAINS

! *********************************************************
! ***** COMPUTE COEFFICIENTS OF SYLVESTER MATRIX **********
! *********************************************************
  SUBROUTINE compsylv48(AA(0:20),BB(0:2),SYLV) 
    IMPLICIT NONE 
! Sylvester matrix elements                                         
    COMPLEX*16,INTENT(IN) :: AA(0:20)
    COMPLEX*16,INTENT(IN) :: BB(0:2)
    COMPLEX*16,INTENT(OUT) :: SYLV(22,22) 
! ======== end interface ==================================
    INTEGER :: i,j
    SYLV=(0.d0,0.d0)
    DO i = 1,21
       SSYLV(i,1) = AA(21-i)
       SSYLV(i+1,2) = AA(21-i)
    ENDDO
    DO j = 3,22
       SSYLV(j-2,j) = BB(2)
       SSYLV(j-1,j) = BB(1)
       SSYLV(j,j) = BB(0)
    ENDDO
  END SUBROUTINE compsylv48

! ******************************************************************
! *** get the coeff's of the polynomials of the Sylvester matrix ***
! ******************************************************************
  SUBROUTINE matrixdat
    IMPLICIT NONE 

A0()=

A20()

B0(0) = 
B0(1) = 
B0(2) = 

B1 = 

B2 = 

  END SUBROUTINE matrixdat
  
END MODULE angular_momentum
