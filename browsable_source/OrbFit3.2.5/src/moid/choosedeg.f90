!     ===========================================================       
!     subroutine to obtain angles in [0,360]                            
!     last modified June 2004 (GFG)
      SUBROUTINE choosedeg(lambda,chl) 
      USE fund_const
      IMPLICIT NONE 
      DOUBLE PRECISION, INTENT(IN) :: lambda                                     
      DOUBLE PRECISION, INTENT(OUT) :: chl 
!     -------------------------------- end interface ------------                   
      DOUBLE PRECISION :: x,y ! sin(lambda),cos(lambda) 
!     ============================================================      
                                               
      x=sin(radeg*lambda)
      y=cos(radeg*lambda)

!     chl in [-180,180]                                    
      chl=degrad*atan2(x,y) 
!     chl in [0,360]                                                    
      IF((chl.ge.0.d0).and.(chl.le.180.d0)) THEN 
         chl = chl 
      ELSEIF((chl.lt.0.d0).and.(chl.ge.-180.d0)) THEN 
         chl = 360.d0 + chl 
      ELSE 
         WRITE(*,*)'choosedeg: ERROR!' 
      ENDIF 
!     control                                                           
      if((chl.gt.360.d0).or.(chl.lt.0.d0)) then 
         write(*,*)'choosedeg: ERROR! IT SHOULD NEVER GET IT!' 
      endif 
                                                                        
      RETURN 
      END SUBROUTINE choosedeg                                          
