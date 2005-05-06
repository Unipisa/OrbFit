!     ==================================================================
!     **************  S O R T I N G   D2   P R O G R A M ***************
!     ******************************************************************
!     ************** written by GIOVANNI F. GRONCHI (2001) *************
!     ********** Department of Mathematics, UNIVERSITY of PISA *********
!     ==================================================================
      SUBROUTINE sortd2(nstat,u,upl,D2,answer) 
      IMPLICIT NONE 
!     input                                                             
      INTEGER nstat 
      DOUBLE PRECISION u(20),upl(20),D2(20) 
      INTEGER answer(20) 
!     --------------------------------- end interface ------------------
!     loop indexes                                                      
      INTEGER i,j 
!     ==================================================================
                                                                        
      DO j = 2,nstat 
         DO i=1,j-1 
                                                                        
            CALL confr(D2(i),u(i),upl(i),answer(i),                     &
     &           D2(j),u(j),upl(j),answer(j))                           
         ENDDO 
      ENDDO 
                                                                        
      RETURN 
      END                                           
                                                                        
                                                                        
!     ================================================================  
      SUBROUTINE confr(d1,u1,upl1,ans1,d2,u2,upl2,ans2) 
      IMPLICIT NONE 
      DOUBLE PRECISION d1,u1,upl1,d2,u2,upl2 
      INTEGER ans1,ans2 
!     ------------------------------ end interface -------------------  
      DOUBLE PRECISION tmpd2,tmpu2,tmpupl2 
      INTEGER tmpans2 
!     ===============================================================   
                                                                        
      IF (d1.lt.d2) THEN 
!     the order is correct                                              
                                                                        
      ELSEIF (d1.ge.d2) THEN 
         tmpd2 = d2 
         d2 = d1 
         d1 = tmpd2 
!                                                                       
         tmpu2 = u2 
         u2 = u1 
         u1 = tmpu2 
!                                                                       
         tmpupl2 = upl2 
         upl2 = upl1 
         upl1 = tmpupl2 
!                                                                       
         tmpans2 = ans2 
         ans2 = ans1 
         ans1 = tmpans2 
      ENDIF 
                                                                        
      RETURN 
      END                                           
