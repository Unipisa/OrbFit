c     ********************************************************************
c     ********** SELECT AMONG MINIMA, MAXIMA and SADLE POINTS ************
c     ********************************************************************
c     *********** written by GIOVANNI F. GRONCHI (1999) ******************
c     ********** Department of Mathematics, UNIVERSITY of PISA ***********
c     ====================================================================
      SUBROUTINE hess(u,upl,ans)
      IMPLICIT NONE
c     u,upl are passed in radians
      DOUBLE PRECISION u,upl
      INTEGER ans
c     ------------------------------------ end interface -----------------
      DOUBLE PRECISION H11,H12,H21,H22,trH,detH
c     eigenvalues
      DOUBLE PRECISION lam1,lam2
c     
      DOUBLE PRECISION A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      DOUBLE PRECISION A11,A12,A13,A14,A15
      COMMON/Aj1to15/ A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15
c     trigonometric constants
      INCLUDE 'trig.h'
c     ===================================================================

c     HESSIAN MATRIX
      H11 = 2.d0*A1*cos(u)**2-2.d0*A1*sin(u)**2+2.d0*A3*sin(u)**2-
     *     2.d0*A3*cos(u)**2-A7*sin(u)*sin(upl)-A8*sin(u)*cos(upl)-
     *     A9*cos(u)*sin(upl)-A10*cos(u)*cos(upl)-A11*sin(u)-A12*cos(u)      
      H12 = A7*cos(u)*cos(upl)-A8*cos(u)*sin(upl)-A9*sin(u)*cos(upl)+
     *     A10*sin(u)*sin(upl)      
      H21 = H12      
      H22 = 2.d0*A4*cos(upl)**2-2.d0*A4*sin(upl)**2+2.d0*A6*sin(upl)**2-
     *     2.d0*A6*cos(upl)**2-A7*sin(u)*sin(upl)-A8*sin(u)*cos(upl)-
     *     A9*cos(u)*sin(upl)-A10*cos(u)*cos(upl)-A13*sin(upl)-
     *     A14*cos(upl)
c     TRACE OF THE HESSIAN
      trH = H11 + H22
c     DETERMINANT OF THE HESSIAN
      detH = H11*H22 - H12*H21
c     check
      IF (abs(detH).le.1.d-10) THEN
         WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(*,*)'!!!!!!!! det(H) = 0 !!!!!!!'
         WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         GOTO 10
      ENDIF
c     compute eigenvalues
      IF ((trH**2-4.d0*detH).ge.0.d0) THEN
         lam1 = (trH + dsqrt(trH**2 - 4.d0*detH))/2.d0
         lam2 = (trH - dsqrt(trH**2 - 4.d0*detH))/2.d0
         IF ((lam1.gt.0.d0).and.(lam2.gt.0.d0)) THEN
            ans = -1
         ELSEIF ((lam1.lt.0.d0).and.(lam2.lt.0.d0)) THEN
            ans = 1
         ELSEIF ((lam1.gt.0.d0).and.(lam2.lt.0.d0)) THEN
            ans = 0            
         ELSEIF ((lam1.lt.0.d0).and.(lam2.gt.0.d0)) THEN
            ans = 0
         ELSE
            ans = 2
         ENDIF    
      ELSE
         WRITE(*,*)'DISCRIMINANT = 0'
      ENDIF
      
 10   CONTINUE
      
      RETURN
      END
