! *************************************
!        S  O  L  V  P  O  L  Y        
! *************************************
! find the real roots of a real polynomial with degree poldeg
 SUBROUTINE solvpoly(poldeg,coef,roots,nroots,hzflag,multfl)
   USE output_control
   IMPLICIT NONE 
   INTEGER, INTENT(IN) :: poldeg ! polynomial degree
   DOUBLE PRECISION, INTENT(IN) :: coef(0:poldeg) ! z^i coeff. (i=0,poldeg)   
   DOUBLE PRECISION, INTENT(OUT) :: roots(poldeg) !real roots
   INTEGER, INTENT(OUT) :: nroots ! number of real roots
   LOGICAL,INTENT(INOUT) :: hzflag ! hzflag = .true.  OK!
                                   ! hzflag = .false. abs(root)>10^5
   LOGICAL,INTENT(INOUT) :: multfl ! multfl = .true.  OK!
                                   ! multfl = .false. 0 has multiplicity > 4
! =============== end interface =======================================
   DOUBLE PRECISION, PARAMETER :: toll=1.d-10
   COMPLEX*16 :: zr1(poldeg) 
   DOUBLE PRECISION :: radius(poldeg) 
   COMPLEX*16 :: poly(0:poldeg)
   DOUBLE PRECISION :: epsm,big,small,rad(poldeg),a1(poldeg+1), &
        & a2(poldeg+1),rim,rre 
   INTEGER :: nit,j,i,pdeg 
   INTEGER :: numzerosol,numzerocoe
   LOGICAL err(poldeg+1) 
! =====================================================================

! initialization
   numzerosol = 0
   numzerocoe = 0
   roots(1:poldeg) = 0.d0
   pdeg = poldeg

! for polzeros
   epsm=2.D0**(-53) 
   big=2.D0**1023 
   small=2.D0**(-1022) 

!   write(*,*)'*resultant coefficients*', coef(0:pdeg)

   poly(0:pdeg)=coef(0:pdeg) 
   
! *************************
! look for zero solutions   
! *************************
   IF(abs(poly(0)).lt.toll) THEN
      WRITE(*,*)'CONSTANT TERM is ZERO!',poly(0),poly(1)
      numzerosol = 1
      pdeg = pdeg-1
      IF(abs(poly(1)).lt.toll) THEN
         WRITE(*,*)'FIRST DEGREE TERM is ZERO!',poly(1)
         numzerosol = 2
         pdeg = pdeg-1
         IF(abs(poly(2)).lt.toll) THEN
            WRITE(*,*)'SECOND DEGREE TERM is ZERO!',poly(2),poly(3)
            numzerosol = 3
            pdeg = pdeg-1
            IF(abs(poly(3)).lt.toll) THEN
               WRITE(*,*)'THIRD DEGREE TERM is ZERO!',poly(3)
               numzerosol = 4
               pdeg = pdeg-1
               IF(abs(poly(4)).lt.toll) THEN
                  WRITE(*,*)'solvpoly: ERROR! &
                       & zero solution has multiplicity > 4'
                  multfl = .false.
               ENDIF
            ENDIF
         ENDIF
      ENDIF
   ENDIF
   
! *******************************************
! check the decrease of the degree up to 16  
! *******************************************
   IF(poly(poldeg).eq.0.d0) THEN
      WRITE(*,*)'20th DEGREE TERM is ZERO!',poly(0)
      pdeg = pdeg-1
      numzerocoe = 1
      IF(poly(poldeg-1).eq.0.d0) THEN
         WRITE(*,*)'19th DEGREE TERM is ZERO!',poly(1)
         pdeg = pdeg-1
         numzerocoe = 2
         IF(poly(poldeg-2).eq.0.d0) THEN
            WRITE(*,*)'18th DEGREE TERM is ZERO!',poly(2)
            pdeg = pdeg-1
            numzerocoe = 3
            IF(poly(poldeg-3).eq.0.d0) THEN
               WRITE(*,*)'17th DEGREE TERM is ZERO!',poly(3)
               pdeg = pdeg-1
               numzerocoe = 4
            ENDIF
         ENDIF
      ENDIF
   ENDIF

! *********************************************************************   
   CALL polzeros(pdeg,poly(numzerosol:poldeg-numzerocoe),epsm,big,small, &
        & 50,zr1,rad,err,nit,a1,a2) 
! *********************************************************************   
   
   nroots=0 
   
   DO 11 j = 1,pdeg
! WRITE(*,*)'ROOT(',j,')=',zr1(j)
!      rim=IMAG(zr1(j))                                                 
      rim=DIMAG(zr1(j)) 
!      IF(ABS(rim).LT.10.d-4) THEN                                      
      IF(ABS(rim).LT.rad(j)) THEN 
         rre=DBLE(zr1(j)) 
! uncomment if you want only positive real roots                    
!      IF(rre.LT.0.D0) GOTO 11                                      
         nroots=nroots+1 
         roots(nroots)=rre 
         radius(nroots)=rad(j) ! error estimate            

! ***************************
! check if some root is big
! ***************************
         IF(abs(roots(nroots)).gt.1.d5) THEN
            if(verb_moid.ge.20) then
               write(*,*)'solvpoly: large value for a root',roots(nroots)
            endif
            hzflag=.false.
!            write(*,*)'hzflag',hzflag
         ENDIF

      END IF
11 END DO
   
!   write(*,*)'roots:',roots(1:nroots)

   IF(numzerosol.gt.0) THEN
      DO i = 1,numzerosol
         nroots=nroots+1 
         roots(nroots)=0.d0 
         radius(nroots)=0.d0 ! dummy
!         write(*,*)'roots(',nroots,')=',roots(nroots)
      ENDDO
   ENDIF
   
 END SUBROUTINE solvpoly
