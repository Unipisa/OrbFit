*
*  ***************************************************************
*  *                                                             *
*  *                          S O L V 20                         *
*  *                                                             *
*  *              Calcolo delle radici reali                     *
*  *              di un polinomio di grado 20                    *
*  *                                                             *
*  ***************************************************************
*
*
* INPUT:    COEF(i)   -  Coefficiente di z^i (i=0,20)
*
* OUTPUT:   ROOTS(k)  -  Radici reali positive (k=1,nroots)
*           NROOTS    -  Numero delle radici reali positive
*
*
      SUBROUTINE solv20(coef,roots,nroots,radius,zr1)
      IMPLICIT NONE

      DOUBLE PRECISION coef(0:20),roots(20),radius(20)
      INTEGER nroots

      COMPLEX*16 poly(0:20),zr1(20)
      DOUBLE PRECISION epsm,big,small,rad(20),a1(21),a2(21),rim,rre
      INTEGER nit,j
      LOGICAL err(21)

      DO 10 j=0,20
      poly(j)=coef(j)
 10   CONTINUE
      epsm=2.D0**(-53)
      big=2.D0**1023
      small=2.D0**(-1022)
      CALL polzeros(20,poly,epsm,big,small,50,zr1,rad,err,nit,a1,a2)

      nroots=0
c     control
c      WRITE(*,*)'CHECK OF ALL THE ROOTS OF THE RESULTANT'

      DO 11 j=1,20
c     control
c         WRITE(*,*)'ROOT(',j,')=',zr1(j)
c         WRITE(*,*)'======================================='

c      rim=IMAG(zr1(j))
         rim=DIMAG(zr1(j))
c      IF(ABS(rim).LT.10.d-4) THEN
c
      IF(ABS(rim).LT.rad(j)) THEN
c
          rre=DBLE(zr1(j))
c     ===================================================
c     uncomment if you want only positive real roots 
c          IF(rre.LT.0.D0) GOTO 11
c     ===================================================
          nroots=nroots+1
          roots(nroots)=rre
c     error estimate
          radius(nroots)=rad(j)

      END IF
 11   CONTINUE

      END
