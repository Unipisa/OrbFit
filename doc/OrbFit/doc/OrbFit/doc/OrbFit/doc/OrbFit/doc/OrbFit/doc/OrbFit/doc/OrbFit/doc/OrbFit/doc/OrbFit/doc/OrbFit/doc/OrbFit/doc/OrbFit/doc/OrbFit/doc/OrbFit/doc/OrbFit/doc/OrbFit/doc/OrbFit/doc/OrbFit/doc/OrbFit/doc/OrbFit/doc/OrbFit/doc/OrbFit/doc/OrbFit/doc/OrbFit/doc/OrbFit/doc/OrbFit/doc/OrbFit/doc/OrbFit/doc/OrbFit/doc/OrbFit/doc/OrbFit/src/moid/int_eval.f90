
! ******************************************************************
! ********** SELECT AMONG MINIMA, MAXIMA and SADLE POINTS **********
! **** by EVALUATING D2 in a NEIGHBORHOOD of the CRITICAL POiNT ****
! ******************************************************************
! ************* written by GIOVANNI F. GRONCHI (2003) **************
! ********** Department of Mathematics, UNIVERSITY of PISA *********
! ==================================================================
  SUBROUTINE int_eval(u,upl,ans) 
    USE fund_const
    IMPLICIT NONE 
!   u,upl are passed in radians                                       
    DOUBLE PRECISION,INTENT(IN) :: u,upl 
    INTEGER,INTENT(OUT) :: ans !  ans = 1: maximum 
!                                 ans = 0: saddle 
!                                 ans= -1: minimum 
!                                 ans= -2: cannot decide
! ----------- end interface ---------------------------------------
    INTEGER,PARAMETER :: n = 30 
    DOUBLE PRECISION :: v,vpl 
    DOUBLE PRECISION :: uD2
    DOUBLE PRECISION :: vD2
    INTEGER :: count ! counter
    INTEGER :: j
! =================================================================

    count = 0

    CALL D2eval(u,upl,uD2) 
!    write(*,*)'uD2',uD2

    DO j = 0,n-1
       v = u + radeg*cos(j*dpig/n)
       vpl = upl + radeg*sin(j*dpig/n)
       CALL D2eval(v,vpl,vD2) 
!    write(*,*)'vD2',vD2

       IF(uD2.lt.vD2) THEN
          count = count - 1
       ELSEIF(uD2.gt.vD2) THEN
          count = count + 1
       ELSE
          WRITE(*,*)'int_eval: equal values of D2'
       ENDIF
    ENDDO

    IF(count.eq.-n) THEN
       ans = -1
    ELSEIF(count.eq.n) THEN
       ans = 1
    ELSE
       ans = 0
    ENDIF

    RETURN
  END SUBROUTINE int_eval
