! ****************************************
! subroutine to obtain angles in [0,360]
! input and output in degrees
! written by G.F.Gronchi Nov.2004
! last modified 21/10/2005 GFG
! ****************************************
 SUBROUTINE choosedeg(lambda,chl)
   USE output_control
   USE fund_const
   IMPLICIT NONE
   INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
   REAL(KIND=qkind), INTENT(IN) :: lambda
   REAL(KIND=qkind), INTENT(OUT) :: chl
! ------------- end interface ------------
   REAL(KIND=qkind) :: x,y !sin(lambda),cos(lambda)
! ===================================================   
   x=sin(radeg*lambda)
   y=cos(radeg*lambda)
! chl in [-180,180]                                                 
   chl=degrad*atan2(x,y)
! chl in [0,360]                                                    
   IF((chl.ge.0.d0).and.(chl.le.180.d0)) THEN 
      chl = chl 
   ELSEIF((chl.lt.0.d0).and.(chl.ge.-180.d0)) THEN 
      chl = 360.d0 + chl 
   ELSE 
      WRITE(iun_log,*)'choosedeg: input angle outside its established range!'&
           & ,chl 
!      STOP
   ENDIF
!     control                                                           
   if((chl.gt.360.d0).or.(chl.lt.0.d0)) then 
      WRITE(iun_log,*)'choosedeg: output angle outside its established range!'&
           & ,chl 
!      STOP
   endif
 END SUBROUTINE choosedeg
