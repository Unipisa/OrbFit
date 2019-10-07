c     ===================================================================
c     **************  S O R T I N G   D2   P R O G R A M ****************
c     *******************************************************************
c     ************** written by GIOVANNI F. GRONCHI (2001) **************
c     ********** Department of Mathematics, UNIVERSITY of PISA **********
c     ===================================================================
      SUBROUTINE sortd2(nstat,u,upl,D2,answer)
      IMPLICIT NONE
c     input
      INTEGER nstat
      DOUBLE PRECISION u(20),upl(20),D2(20)
      INTEGER answer(20)
c     --------------------------------- end interface -------------------
c     loop indexes
      INTEGER i,j
c     ===================================================================
 
      DO j = 2,nstat
         DO i=1,j-1
            
            CALL confr(D2(i),u(i),upl(i),answer(i),
     *           D2(j),u(j),upl(j),answer(j))             
         ENDDO
      ENDDO
      
      RETURN
      END


c     ================================================================
      SUBROUTINE confr(d1,u1,upl1,ans1,d2,u2,upl2,ans2)
      IMPLICIT NONE
      DOUBLE PRECISION d1,u1,upl1,d2,u2,upl2
      INTEGER ans1,ans2
c     ------------------------------ end interface -------------------
      DOUBLE PRECISION tmpd2,tmpu2,tmpupl2
      INTEGER tmpans2
c     ===============================================================

      IF (d1.lt.d2) THEN
c     the order is correct

      ELSEIF (d1.ge.d2) THEN
         tmpd2 = d2
         d2 = d1
         d1 = tmpd2
c
         tmpu2 = u2
         u2 = u1
         u1 = tmpu2
c
         tmpupl2 = upl2
         upl2 = upl1
         upl1 = tmpupl2
c
         tmpans2 = ans2
         ans2 = ans1
         ans1 = tmpans2
      ENDIF
      
      RETURN
      END
