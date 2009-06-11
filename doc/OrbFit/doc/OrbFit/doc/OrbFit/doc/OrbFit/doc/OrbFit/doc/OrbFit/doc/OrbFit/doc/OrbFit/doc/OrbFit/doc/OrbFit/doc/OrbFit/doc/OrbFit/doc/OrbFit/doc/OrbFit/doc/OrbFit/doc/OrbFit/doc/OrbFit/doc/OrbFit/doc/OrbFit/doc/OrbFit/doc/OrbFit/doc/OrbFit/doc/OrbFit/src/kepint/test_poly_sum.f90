PROGRAM test_poly_sum
  IMPLICIT NONE
  INTEGER,PARAMETER :: amod=2,bmod=2
  DOUBLE PRECISION ::  calpha(0:amod,0:amod),dbeta(0:bmod,0:bmod)

  INTEGER :: gmod ! modulus of the multi-index
!  DOUBLE PRECISION :: egamma(0:amod+bmod,0:amod+bmod)
  DOUBLE PRECISION :: egamma(0:max(amod,bmod),0:max(amod,bmod))
  INTEGER :: i,j ! loop indexes

  calpha=0.d0;dbeta=0.d0

  calpha(0,0)= 0.d0
  calpha(0,1)= 0.d0
  calpha(0,2)= 1.d0
  calpha(1,0)= 0.d0
  calpha(1,1)= 0.d0
  calpha(2,0)= 1.d0

  dbeta(0,0)= 0.d0
  dbeta(0,1)= 0.d0
  dbeta(0,2)= -1.d0
  dbeta(1,0)= 0.d0
  dbeta(1,1)= 0.d0
  dbeta(2,0)= 1.d0

  CALL poly_sum(amod,bmod,calpha,dbeta,gmod,egamma)

  DO i=0,gmod
!     write(*,*)'i=',i
     DO j=0,gmod-i
        WRITE(*,*)'egamma(',i,j,')=',egamma(i,j)
     ENDDO
  ENDDO

END PROGRAM test_poly_sum
