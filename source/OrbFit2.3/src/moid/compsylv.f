c     ******************************************************************
c     ===========================================================
c     COMPUTE COEFFICIENTS OF SYLVESTER MATRIX
c     ===========================================================
      SUBROUTINE compsylv(aleph,beth,ghimel,A,B,C,D,E,SSYLV)
      IMPLICIT NONE
c     =========== Sylvester matrix elements =================
      DOUBLE PRECISION aleph,beth,ghimel
      DOUBLE PRECISION A,B,C,D,E
      DOUBLE PRECISION SSYLV(6,6)
      DOUBLE PRECISION SYLV(6)
c     loop index
      INTEGER h,k
c     ******************************************************************

      SSYLV(1,1) = aleph
      SSYLV(1,2) = 0
      SSYLV(1,3) = 0
      SSYLV(1,4) = 0
      SSYLV(1,5) = A
      SSYLV(1,6) = 0
c
      SSYLV(2,1) = beth
      SSYLV(2,2) = aleph
      SSYLV(2,3) = 0
      SSYLV(2,4) = 0
      SSYLV(2,5) = B
      SSYLV(2,6) = A
c
      SSYLV(3,1) = ghimel
      SSYLV(3,2) = beth
      SSYLV(3,3) = aleph
      SSYLV(3,4) = 0
      SSYLV(3,5) = C
      SSYLV(3,6) = B
c
      SSYLV(4,1) = 0
      SSYLV(4,2) = ghimel
      SSYLV(4,3) = beth
      SSYLV(4,4) = aleph
      SSYLV(4,5) = D
      SSYLV(4,6) = C
c
      SSYLV(5,1) = 0
      SSYLV(5,2) = 0
      SSYLV(5,3) = ghimel
      SSYLV(5,4) = beth
      SSYLV(5,5) = E
      SSYLV(5,6) = D
c
      SSYLV(6,1) = 0
      SSYLV(6,2) = 0
      SSYLV(6,3) = 0
      SSYLV(6,4) = ghimel
      SSYLV(6,5) = 0
      SSYLV(6,6) = E


      RETURN
      END

