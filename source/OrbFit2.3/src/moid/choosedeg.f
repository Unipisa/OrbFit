c     ===========================================================
c     subroutine to obtain angles in [0,360]
      SUBROUTINE choosedeg(lambda,chl)
      IMPLICIT NONE
c     INPUT
      DOUBLE PRECISION lambda
c     OUTPUT
      DOUBLE PRECISION chl
c     -------------------------------- end interface ------------
c     sin(lambda),cos(lambda)
      DOUBLE PRECISION x,y
c     trig. constants
      INCLUDE 'trig.h'
c     ============================================================

      x=sin(radeg*lambda)
      y=cos(radeg*lambda)

c     chl in [-180,180]
      chl=degrad*atan2(x,y)
c     chl in [0,360]
      IF((chl.ge.0).and.(chl.le.180)) THEN 
         chl = chl
      ELSEIF((chl.lt.0).and.(chl.ge.-180)) THEN
         chl = 360 + chl
      ELSE
         WRITE(*,*)'ERROR!!!!'
      ENDIF
c     control
      if(chl.gt.360) then
         write(*,*)'ERROR: IT SHOULD NEVER GET IT!'
      endif

      RETURN
      END
