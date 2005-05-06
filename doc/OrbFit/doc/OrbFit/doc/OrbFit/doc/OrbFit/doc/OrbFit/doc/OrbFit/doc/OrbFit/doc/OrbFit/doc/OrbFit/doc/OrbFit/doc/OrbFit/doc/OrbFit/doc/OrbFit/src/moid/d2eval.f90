
! ******************************************************************
! ******** EVALUATE THE SQUARED DISTANCE IN THE POINTS u,upl *******
! ******************************************************************
! *********** written by GIOVANNI F. GRONCHI (2001) ****************
! ******** Department of Mathematics, UNIVERSITY of PISA ***********
! last modified May 2003
! ==================================================================
  SUBROUTINE D2eval(u,upl,D2) 
    USE fund_const
    IMPLICIT NONE                            
    DOUBLE PRECISION,INTENT(IN) :: u,upl ! eccentric anomalies (in rad.)
    DOUBLE PRECISION D2 ! SQUARED DISTANCE function 
!   ------------- end interface --------------------------------------
!   functions of the orbital elements                                 
    DOUBLE PRECISION A1,A2,A3,A4,A5,A6,A7 
    DOUBLE PRECISION A8,A9,A10,A11,A12,A13,A14,A15 
    COMMON/Aj1to15/ A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15 
!   ==================================================================
                                                                        

    D2=A1*(sin(u)**2) + A3*(cos(u)**2) + A4*(sin(upl)**2) +            &
         & A6*(cos(upl)**2) + A7*sin(u)*sin(upl) + A8*sin(u)*cos(upl)  &
         & + A9*cos(u)*sin(upl) + A10*cos(u)*cos(upl) + A11*sin(u) +   &
         & A12*cos(u) + A13*sin(upl) + A14*cos(upl) + A15               
                                                                        
!   check if D2 i positive; otherwise set it to zero                  
    IF (D2.ge.0) THEN 
!   do nothing                                                        
    ELSEIF (D2.lt.0) THEN 
       D2 = 0.d0 
    ENDIF
    
    RETURN 
  END SUBROUTINE D2eval
