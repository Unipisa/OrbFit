c     ********************************************************************
c     ******** EVALUATE THE SQUARED DISTANCE IN THE POINTS u,upl *********
c     ********************************************************************
c     *********** written by GIOVANNI F. GRONCHI (2001) ******************
c     ******** Department of Mathematics, UNIVERSITY of PISA *************
c     ====================================================================
      SUBROUTINE D2eval(u,upl,D2)     
      IMPLICIT NONE
c     eccentric anomalies (in radians)
      DOUBLE PRECISION u,upl
c     SQUARED DISTANCE function
      DOUBLE PRECISION D2
c     --------------------------------- end interface --------------------
c     functions of the orbital elements
      DOUBLE PRECISION A1,A2,A3,A4,A5,A6,A7
      DOUBLE PRECISION A8,A9,A10,A11,A12,A13,A14,A15
      COMMON/Aj1to15/ A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15
c     trig constants
      INCLUDE 'trig.h'
c     ====================================================================
      
c     ============= SQUARED DISTANCE function ================
      D2=A1*(sin(u)**2) + A3*(cos(u)**2) + A4*(sin(upl)**2) +
     *     A6*(cos(upl)**2) + A7*sin(u)*sin(upl) + A8*sin(u)*cos(upl)
     *     + A9*cos(u)*sin(upl) + A10*cos(u)*cos(upl) + A11*sin(u) +
     *     A12*cos(u) + A13*sin(upl) + A14*cos(upl) + A15

      RETURN
      END
