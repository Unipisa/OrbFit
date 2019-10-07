! *******************************************
! subroutine to obtain angles in [-180,180]
! written by G.F.Gronchi Nov.2004
! *******************************************
 SUBROUTINE choosedeg(lambda,chl) 
   USE fund_const
   IMPLICIT NONE                                       
   REAL(KIND=8), INTENT(IN) :: lambda 
   REAL(KIND=8), INTENT(OUT) :: chl 
! ------------- end interface ------------                   
   REAL(KIND=8) :: x,y !     sin(lambda),cos(lambda)     
! ===================================================      
   
   x=sin(radeg*lambda) 
   y=cos(radeg*lambda) 
                                                                        
!     chl in [-180,180]                                                 
   chl=degrad*atan2(x,y)
!     chl in [0,360]                                                    
   IF((chl.ge.0).and.(chl.le.180)) THEN 
      chl = chl 
   ELSEIF((chl.lt.0).and.(chl.ge.-180)) THEN 
      chl = 360 + chl 
   ELSE 
      WRITE(*,*)'ERROR!!!!' 
   ENDIF
!     control                                                           
   if(chl.gt.360) then 
      write(*,*)'ERROR: IT SHOULD NEVER GET IT!' 
   endif
   
   RETURN 
 END SUBROUTINE choosedeg
 
