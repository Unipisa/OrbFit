c     ********************************************************************
c     GIVEN A VALUE OF t COMPUTES THE CORRESPONDING VALUE OF s THAT
c     SATISFIES THE SYSTEM 
c     ********************************************************************
c     ************ written by GIOVANNI F. GRONCHI (2001) *****************
c     ********** Department of Mathematics, UNIVERSITY of PISA ***********
c     ====================================================================
      SUBROUTINE solvesystem(nroots,tzero,szero,sflag)
      IMPLICIT NONE
      INTEGER nroots 
      DOUBLE PRECISION tzero(20),szero(20)
c     solving system flag
c     sflag(1) = 1: OK! ;  
c     sflag(1) = 0: WARNING! They are both good solutions!;
c     sflag(2) = 1: OK! ;  
c     sflag(2) = 0: WARNING! Neither of two is close to 0!
      INTEGER sflag(2)
c     ------------------------------ end interface -----------------------
      DOUBLE PRECISION stmp(2)
      DOUBLE PRECISION su,cu,u,supl(2),cupl(2),upl(2)
      DOUBLE PRECISION suplchk(20),cuplchk(20)
      DOUBLE PRECISION suchk(20),cuchk(20)
      DOUBLE PRECISION uu(20),v 
      DOUBLE PRECISION deg2s,deg1s,deg0s
c     
      DOUBLE PRECISION A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      DOUBLE PRECISION A11,A12,A13,A14,A15
      COMMON/Aj1to15/ A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15
c     
      DOUBLE PRECISION t,s,evalst(2)
c     choice flag ( choice = 0 ------> init value
c                 ( choice = 1 ------> stmp(1)
c                 ( choice = 2 ------> stmp(2)
      INTEGER choice
c     smallness parameter
      DOUBLE PRECISION eps
c     loop indexes
      INTEGER i,j,k
c     ====================================================================

c     flag initialization
c      sflag(1) = 1
c      sflag(2) = 1
      choice = 0
      eps = 1.d-5

      DO i = 1,nroots
         t=tzero(i)

c     ==============================================
c     deg2s(t)*(s**2) + deg1s(t)*(s) + deg0s(t) = 0 
         deg2s = -A8 + A11 + 4.d0*A1*t - 4.d0*A1*t**3 - 4.d0*A3*t +
     *        4.d0*A3*t**3 + A8*t**4 + 2.d0*A10*t + 2.d0*A10*t**3 -
     *        A11*t**4 - 2.d0*A12*t - 2.d0*A12*t**3
         deg1s = 2.d0*A7 - 2.d0*A7*t**4 - 4.d0*A9*t**3 - 4.d0*A9*t
         deg0s = 4.d0*A1*t - 4.d0*A3*t - 2.d0*A10*t + A8 + A11 - 
     *        4.d0*A1*t**3 + 4.d0*A3*t**3 - A8*t**4 - 2.d0*A10*t**3 -
     *        A11*t**4 - 2.d0*A12*t**3 - 2.d0*A12*t
         
c     NORMALIZING COEFFICIENTS (the order is IMPORTANT!)
         IF(abs(deg2s).lt.1.d-18)THEN
            WRITE(*,*)'BETTER NOT TO NORMALIZE!'
            GOTO 124
         ENDIF
         deg0s = deg0s/deg2s
         deg1s = deg1s/deg2s
         deg2s = deg2s/deg2s
 124     CONTINUE
         
c     solutions
         stmp(1) = (-deg1s + dsqrt(deg1s**2 - 4.d0*deg0s*deg2s))/
     *        (2.d0*deg2s) 
         stmp(2) = (-deg1s - dsqrt(deg1s**2 - 4.d0*deg0s*deg2s))/
     *        (2.d0*deg2s) 

c     ======================================================
c     deg4s(t)*(s**4) + deg3s(t)*(s**3) + deg2s(t)*(s**2) + 
c     deg1s(t)*(s) + deg0s(t) = 0 

c     conversion into eccentric anomalies
         su = 2.d0*tzero(i)/(1.d0+(tzero(i))**2) 
         cu = (1.d0-(tzero(i))**2)/(1.d0+(tzero(i))**2) 
         u = datan2(su,cu)
         uu(i) = u
         suchk(i) = su
         cuchk(i) = cu
c     
         DO j = 1,2
            supl(j) =  2.d0*stmp(j)/(1.d0+(stmp(j))**2) 
            cupl(j) = (1.d0-(stmp(j))**2)/(1.d0+(stmp(j))**2) 
            upl(j) = datan2(supl(j),cupl(j))
            v = upl(j)

            evalst(j) = 2.d0*A4*sin(v)*cos(v)-2.d0*A6*cos(v)*sin(v)+
     *           A7*sin(u)*cos(v) - A8*sin(u)*sin(v) + A9*cos(u)*cos(v)
     *           - A10*cos(u)*sin(v) + A13*cos(v) - A14*sin(v)            
         ENDDO
         
c     SELECTING the CORRESPONDING SOLUTION
c     ************ checking smallness of the evaluations *************
         IF((abs(evalst(1)).lt.eps).or.(abs(evalst(2)).lt.eps))THEN
c     choosing the one that gives the minimum value
            IF (abs(evalst(1)).lt.abs(evalst(2))) THEN
               choice = 1
               szero(i) = stmp(1)
               suplchk(i) =  supl(1) 
               cuplchk(i) =  cupl(1) 
            ELSEIF (abs(evalst(1)).gt.abs(evalst(2))) THEN
               choice = 2
               szero(i) = stmp(2)
               suplchk(i) =  supl(2) 
               cuplchk(i) =  cupl(2) 
            ENDIF            
            
            
c     if they are BOTH CLOSE TO ZERO
            IF((abs(evalst(1)).lt.eps).and.
     *           (abs(evalst(2)).lt.eps))THEN
c     ***********************
               sflag(1) = 0
c     ***********************
c     writing on screen
c               WRITE(*,*)'++++++++++++++++++++++++++++++++++++++++++++'
c               WRITE(*,*)'SOLVING SYSTEM WARNING!     '
c               WRITE(*,*)'THEY ARE BOTH GOOD SOLUTIONS'
c               WRITE(*,*)'EVALST(1)=',evalst(1),'  EVALST(2)=',evalst(2)
c               WRITE(*,*)'++++++++++++++++++++++++++++++++++++++++++++'
            ENDIF
            
c     ================== DOUBLE CORRESPONDENCE CHECK ====================
c     if sflag(1) has become 0, then this check has to be done for
c     all the following values of i

c     -------------------------------------------------------------
c     CASE 1 (i=1)
            IF((sflag(1).eq.0).and.(i.eq.1)) THEN               
c     definitively choosing previously chosen value of s 
c     -------------------------------------------------------------
c     CASE 2 (i>1)
            ELSEIF((sflag(1).eq.0).and.(i.gt.1))THEN
c     at first choosing previously chosen value of s
               DO k = 1,i-1
                  
                  IF( (abs(suchk(i)-suchk(k)).lt.1.d-1).and.
     *                 (abs(cuchk(i)-cuchk(k)).lt.1.d-1).and.
     *                 (abs(suplchk(i)-suplchk(k)).lt.1.d-1).and.
     *                 (abs(cuplchk(i)-cuplchk(k)).lt.1.d-1) )THEN
c     write(*,*)'OK, MODIFYING szero to',stmp(2)
                     
                     IF(choice.eq.1) THEN
                        szero(i) = stmp(2)  
                        suplchk(i) = supl(2) 
                        cuplchk(i) = cupl(2)
                     ELSEIF(choice.eq.2) THEN
                        szero(i) = stmp(1)  
                        suplchk(i) = supl(1) 
                        cuplchk(i) = cupl(1)
                     ENDIF
                     
                  ELSE
c     leaving previous value of s
                     
                  ENDIF
               ENDDO
            ENDIF    
            
c     neither of the evaluation is close to zero
         ELSE
c     ***********************
            sflag(2) = 0
c     ***********************
c            WRITE(*,*)'++++++++++++++++++++++++++++++++++++++++++++'
c            WRITE(*,*)'        SOLVING SYSTEM WARNING!        '
c            WRITE(*,*)' NEITHER OF THE EVALUATIONS IS CLOSE TO ZERO'
c            WRITE(*,*)' SELECTING THE ONE WITH THE SMALLEST VALUE '
c            WRITE(*,*)'++++++++++++++++++++++++++++++++++++++++++++'
c            WRITE(*,*)'EVALST(1)=',evalst(1),'EVALST(2)=',evalst(2)
c            WRITE(*,*)'THE SMALLEST IS',min(abs(evalst(1)),
c     *           abs(evalst(2)))
c     choosing the one that gives the minimum value
            IF (abs(evalst(1)).lt.abs(evalst(2))) THEN
               szero(i) = stmp(1)
            ELSEIF (abs(evalst(1)).gt.abs(evalst(2))) THEN
               szero(i) = stmp(2)
            ENDIF

         ENDIF
         
      ENDDO
      
      RETURN 
      END
