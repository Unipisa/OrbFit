c     ********************************************************************
c     ********** CODING of the n evaluations of a polynomial *************
c     ************* p: p(om(1)), p(om(n)) to apply idft ******************
c     ********************************************************************
c     *********** written by GIOVANNI F. GRONCHI (2001) ******************
c     ********** Department of Mathematics, UNIVERSITY of PISA ***********
c     ====================================================================
      SUBROUTINE code_inp(N,decoddftout,idftinp)
      IMPLICIT NONE
      INTEGER N
      COMPLEX*16 decoddftout( * )
      DOUBLE PRECISION idftinp( * )
c     ------------------------------------- end interface ----------------
      DOUBLE PRECISION RE(1000),IM(1000)
c     loop indexes
      INTEGER i
c     ====================================================================

      DO i = 1,N
         RE(i)=DBLE(decoddftout(i))
         IM(i)=DIMAG(decoddftout(i))
      ENDDO
      
      idftinp(1) = RE(1)
      DO i = 2,N/2
         idftinp(i) = RE(i)
         idftinp(N/2+i) = IM(N/2-i+2)  
      ENDDO
      idftinp(N/2+1) = RE(N/2+1)
          
      RETURN
      END
