!     ******************************************************************
!     OUTPUT CONVERSION FROM THE STANDARD FORM of DTF TO THE USUAL ONE  
!     ******************************************************************
!     *********** written by GIOVANNI F. GRONCHI (2001) ****************
!     ********** Department of Mathematics, UNIVERSITY of PISA *********
!     ==================================================================
      SUBROUTINE decode_out(N,dftout1,dftout2,dfteval) 
      IMPLICIT NONE 
      INTEGER N 
      DOUBLE PRECISION dftout1( * ),dftout2( * ) 
      COMPLEX*16 dfteval( * ) 
!     ----------------------------------- end interface ----------------
!     loop indexes                                                      
      INTEGER j 
!     ==================================================================
      dfteval(1) = DCMPLX(dftout1(1),0.d0) 
      DO j = 2,N/2 
         dfteval(j) = DCMPLX(dftout1(j),dftout1(N+2-j)) 
         dfteval(N/2+j) = DCMPLX(dftout2(N/2+j),dftout2(N/2+2-j)) 
      ENDDO 
      dfteval(N/2+1) = DCMPLX(dftout1(N/2+1),0.d0) 
                                                                        
      RETURN 
      END                                           
