c     ************************************************************
c     ** C O M P U T I N G   S T A T I O N A R Y    P O I N T S **
c     ***  and   W R I T I N G   on   O U T P U T    F I L E S ***
c     ************************************************************
c     ********** written by GIOVANNI F. GRONCHI (2001) ***********
c     ******* Department of Mathematics, UNIVERSITY of PISA ******
c     ============================================================
      SUBROUTINE writestat(name,lnam,iunit,iun1,iun2,iun3,iun4,
     *    ekpl,elkep)    
      IMPLICIT NONE
c     asteroid name and name lenght
      CHARACTER*9  name
      INTEGER lnam
c     file identifiers
      INTEGER iunit,iun1,iun2,iun3,iun4
c     elements of the Earth and of the asteroid
      DOUBLE PRECISION ekpl(6),elkep(6)
c     ----------------------------------- end interface -------------------
c     eccentric anomalies 
      DOUBLE PRECISION u(20),upl(20)
c     SQUARED DISTANCE function 
      DOUBLE PRECISION D2(20)
c     
      INTEGER nstat,nummin,nummax,answer(20)
c     ========================= check flags ==============================
c     morse:       1 = "OK!" ;   0 = "CHECK with MORSE FAILED!"
c     weier:       1 = "OK!" ;   0 = "CHECK with WEIERSTRASS FAILED!"
c     warnflag(1): 1 = "OK!" ;   0 = "lead. coef. of resultant = 0!"
c     warnflag(2): 1 = "OK!" ;   0 = "higher deg. terms not small!"
c     warnflag(3): 1 = "OK!" ;   0 = "low precision in res. coef. (fft)!"
c     sflag(1):    1 = "OK!" ;   0 = "WARNING! They are both good solutions!"
c     sflag(2):    1 = "OK!" ;   0 = "WARNING! Neither of two is close to 0!"
      INTEGER morse,weier,warnflag(3),sflag(2)
c
      CHARACTER*7 singtype
c     loop index
      INTEGER j
c     =====================================================================

c     flags initialization
c      weier = 1
c      morse = 1
c      sflag(1) = 1
c      sflag(2) = 1
c      warnflag(1) = 1
c      warnflag(2) = 1
c      warnflag(2) = 1

      CALL compute_statpts(ekpl,elkep,u,upl,D2,
     *   nstat,nummin,nummax,answer,
     *   morse,weier,warnflag,sflag)

c     *******************************************************************
c     ******************** WRITING OUTPUT ON SCREEN *********************
c     *******************************************************************
      WRITE(*,*)'#######################################################
     *########'
      WRITE(*,*)'######### S T A T I O N A R Y   P O I N T S ###########
     *########'
      WRITE(*,*)'#######################################################
     *########'
      WRITE(*,*)'       u              upl             DIST          TYP
     *E'
      WRITE(*,*)'=======================================================
     *========'
c     loop on number of stationary points 
      DO j = 1,nstat
         write(*,*)
         IF (answer(j).eq.1) THEN
            singtype='MAXIMUM'
         ELSEIF (answer(j).eq.-1) THEN
            singtype='MINIMUM'
         ELSEIF (answer(j).eq.0) THEN
            singtype='SADDLE'
         ELSEIF (answer(j).eq.2) THEN
            singtype='ERROR'
         ENDIF
         
c     writing on screen
         WRITE(*,105)u(j),upl(j),dsqrt(D2(j)),singtype
c     writing on file iunit.statpts
         WRITE(iunit,105)u(j),upl(j),dsqrt(D2(j)),singtype
 105     FORMAT(2x,f13.8,3x,f13.8,3x,f13.8,5x,a7)
         
c     END MAIN DO LOOP
      ENDDO

c     ******************************************************************
c         *******  W R I T I N G  on  C H E C K   F I L E S  *******
c     ******************************************************************

c     ============== SOLVING SYSTEM WARNINGS ==============
      IF(sflag(1).eq.0) THEN
c         WRITE(iun1,*)'++++++++++++++++++++++++++++++++++++++++++++'
c         WRITE(iun1,*)'SOLVING SYSTEM WARNING!     '
c         WRITE(iun1,*)'THERE ARE TWO EQUALLY GOOD SOLUTIONS'
c         WRITE(iun1,*)'NEGATIVE CHECK FOR ',name(1:lnam)
c         WRITE(iun1,*)'++++++++++++++++++++++++++++++++++++++++++++'
      ENDIF
      IF(sflag(2).eq.0) THEN
c         WRITE(iun1,*)'++++++++++++++++++++++++++++++++++++++++++++'
c         WRITE(iun1,*)'        SOLVING SYSTEM WARNING!        '
c         WRITE(iun1,*)' NEITHER OF THE EVALUATIONS IS CLOSE TO ZERO'
c         WRITE(iun1,*)' SELECTING THE ONE WITH THE SMALLEST VALUE '
c         WRITE(iun1,*)'NEGATIVE CHECK FOR ',name(1:lnam)
c         WRITE(iun1,*)'++++++++++++++++++++++++++++++++++++++++++++'
      ENDIF

c     ============= ADDITIONAL WARNING FLAGS ==============

      IF(warnflag(1).eq.0) THEN
c         WRITE(iun2,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
c         WRITE(iun2,*)'!!! LEADING COEFFICIENT = 0 !!!'
c         WRITE(iun2,*)'NEGATIVE CHECK FOR ',name(1:lnam)
c         WRITE(iun2,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      ENDIF
      IF(warnflag(2).eq.0) THEN
c         WRITE(iun2,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
c         WRITE(iun2,*)'!!!!!!!! WARNING !!!!!!!!!!'
c         WRITE(iun2,*)'!!! HIGHER DEGREE TERMS !!!'
c         WRITE(iun2,*)'!!!!!!!! NOT SMALL !!!!!!!!'
c         WRITE(iun2,*)'NEGATIVE CHECK FOR ',name(1:lnam)
c         WRITE(iun2,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      ENDIF
      IF(warnflag(3).eq.0) THEN
c         WRITE(iun2,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
c         WRITE(iun2,*)'FFT WARNING: DIFFERENCE IN THE EVALUATIONS' 
c         WRITE(iun2,*)'OF THE RESULTANT IS NOT SMALL!!!:'
c         WRITE(iun2,*)'NEGATIVE CHECK FOR ',name(1:lnam)
c         WRITE(iun2,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      ENDIF

c     ============= CHECK WITH MORSE THEORY ===============
      IF (morse.eq.0) THEN
c     writing on screen
         WRITE(*,*)'ERROR: MORSE THEORY VIOLATION!!!'
         WRITE(*,*)'FAILED COMPUTATION FOR ',name(1:lnam)
c     writing on CHECK.morse_weier file
         WRITE(iun3,*)'ERROR: MORSE THEORY VIOLATION!!!'                
         WRITE(iun3,*)'FAILED COMPUTATION FOR ',name(1:lnam)
         WRITE(iun3,*)'NUMSTAT=',nstat
         WRITE(iun3,*)'NUMMIN=',nummin
         WRITE(iun3,*)'NUMMAX=',nummax
      ENDIF
c     =====================================================

c     ============ CHECK WITH WEIERSTRASS THEOREM =========
      IF(weier.eq.0)THEN
c     writing on screen
         WRITE(*,*)'ERROR: WEIERSTRASS THEOREM VIOLATION!!!'
         WRITE(*,*)'FAILED COMPUTATION FOR ',name(1:lnam)
c     writing on CHECK.morse_weier file
         WRITE(iun3,*)'ERROR: WEIERSTRASS THEOREM VIOLATION!!!'
         WRITE(iun3,*)'FAILED COMPUTATION FOR ',name(1:lnam)
      ENDIF
c     ===============================================================

c     ===============================================================
c     IF PREVIOUS CHECK ARE SUCCESSFULL, THEN CHECK FOR HIGH NUMBER
c     OF MINIMA (>2) AND STATIONARY POINTS (>=8)
c     ===============================================================
      IF((weier.ne.0).and.(morse.ne.0))THEN
c     write on file CHECK.hi-st
         IF (nummin.gt.2) THEN
            write(iun4,*)'THERE ARE',nummin,' MINIMA for ',name(1:lnam)
         ENDIF
         IF ((nummin+nummax)*2.ge.8) THEN
            write(iun4,*)'THERE ARE',(nummin+nummax)*2,' STATPTS for ',
     *           name(1:lnam)
            write(iun4,*)'===================================='
         ENDIF
      ENDIF
c     ================================================================
      
      WRITE(*,*)
      WRITE(*,*)'=======================================================
     *========'
      WRITE(*,*)'#######################################################
     *########'

      
      RETURN
      END
