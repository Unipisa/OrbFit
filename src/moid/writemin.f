c     ************************************************************
c     **** C O M P U T I N G     M I N I M U M    P O I N T S ****
c     ***  and   W R I T I N G   on   O U T P U T    F I L E S ***
c     ************************************************************
c     ********** written by GIOVANNI F. GRONCHI (2001) ***********
c     ******* Department of Mathematics, UNIVERSITY of PISA ******
c     ============================================================
      SUBROUTINE writemin(name,lnam,iunit,iun1,iun2,iun3,iun4,
     *     ekpl,elkep)    
      IMPLICIT NONE
c     asteroid name and name lenght
      CHARACTER*9  name
      INTEGER lnam
c     file identifiers
      INTEGER iunit,iun1,iun2,iun3,iun4
c     elements of the Earth and of the asteroid
      DOUBLE PRECISION ekpl(6),elkep(6)
c     in cartesian coordinates
      DOUBLE PRECISION car(6),carpl(6)
c     cartesian coordinates at minima
      DOUBLE PRECISION cmin(6,20),cplmin(6,20)
c     check variables
      DOUBLE PRECISION contrmin(6),contrplmin(6)
      DOUBLE PRECISION contrelkep(6),contrekpl(6)
c     
      DOUBLE PRECISION enne
c     ----------------------------------- end interface -------------------
c     eccentric anomalies 
      DOUBLE PRECISION umin(20),uplmin(20)
c     SQUARED DISTANCE function 
      DOUBLE PRECISION D2(20)
c     
      INTEGER nstat,nummin,nummax
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
c     INCLUDE 'pldata.h'
c     data for masses
      INCLUDE 'jplhdr.h'
c     masses
      INCLUDE 'sunmass.h'
      INCLUDE 'parbep.h'
      INCLUDE 'masses.h'
c     =====================================================================
      
c     flags initialization
c      weier = 1
c      morse = 1
c      sflag(1) = 1
c      sflag(2) = 1
c      warnflag(1) = 1
c      warnflag(2) = 1
c      warnflag(2) = 1
      
c     duplicate computation of gmse, in case masjpl has not been called yet
      gmse=gms*(1.d0+cval(11)/cval(18))
      
c     check
c     write(*,*)'elkep prima:',elkep
c     write(*,*)'ekpl prima:',ekpl
      
c     switching to cartesian coordinates
      CALL coocha(elkep,'KEP',gms,car,'CAR',enne)
      CALL coocha(ekpl,'KEP',gmse,carpl,'CAR',enne)
      
c     check
c     write(*,*)'gmse,gms,enne',gmse,gms,enne
      write(*,*)'car:',car
      write(*,*)'carpl:',carpl
      
c     additional check
c     CALL coocha(car,'CAR',gms,elkep,'KEP',enne)
c     CALL coocha(carpl,'CAR',gmse,ekpl,'KEP',enne)
c     write(*,*)'elkep dopo:',elkep
c     write(*,*)'ekpl dopo:',ekpl
      
c     =======================================================
      CALL compute_minima(car,carpl,cmin,cplmin,D2,
     *     nummin,morse,weier,warnflag,sflag)
c     =======================================================
      
c     *******************************************************************
c     ******************** WRITING OUTPUT ON SCREEN *********************
c     *******************************************************************
      WRITE(*,*)'##################################################'
      WRITE(*,*)'######### M I N I M U M    P O I N T S ###########'
      WRITE(*,*)'##################################################'
      WRITE(*,*)'=======================================================
     *======================='
c     loop on number of stationary points 
      DO j = 1,nummin
         write(*,*)
c     writing on screen
         WRITE(*,*)'minimum label =',j,'      Distance =',dsqrt(D2(j))
         WRITE(*,*)' '
         WRITE(*,*)'Asteroid cartesian coordinates at minimum         '
         WRITE(*,*)'      x             y              z         derx     
     $      dery        derz'
         WRITE(*,107)cmin(1,j),cmin(2,j),cmin(3,j),cmin(4,j),cmin(5,j)
     $        ,cmin(6,j)
         WRITE(*,*)'Earth cartesian coordinates at minimum            '
         WRITE(*,*)'      xpl           ypl            zpl       derxpl     
     $      derypl      derzpl'
         WRITE(*,107)cplmin(1,j),cplmin(2,j),cplmin(3,j),cplmin(4,j)
     $        ,cplmin(5,j),cplmin(6,j)
         WRITE(*,*)'====================================================
     *=========================='
c     writing on file iunit.minpts
         WRITE(iunit,106)j,dsqrt(D2(j))
         WRITE(iunit,107)cmin(1,j),cmin(2,j),cmin(3,j),cmin(4,j),cmin(5
     $        ,j),cmin(6,j)
         WRITE(iunit,107)cplmin(1,j),cplmin(2,j),cplmin(3,j),cplmin(4,j)
     $        ,cplmin(5,j),cplmin(6,j)
 106     FORMAT(2x,i2,15x,f13.8)
 107     FORMAT(2x,f11.6,2x,f11.6,2x,f11.6,2x,f11.6,2x,f11.6,2x,f11.6)

c     check on the distances
c         write(*,*)'-----------------------------------------'
c         WRITE(*,*)'controllo di D2(',j,')=',dsqrt((cmin(1,j)-cplmin(1,j
c     $        ))**2+(cmin(2,j)-cplmin(2,j))**2+(cmin(3,j)-cplmin(3,j))
c     $        **2)
c         write(*,*)'-----------------------------------------'
         
c     END MAIN DO LOOP
      ENDDO
      

c     ========================== CHECK ===========================
c      DO j = 1,6
c         contrmin(j) = cmin(j,1)
c         contrplmin(j) = cplmin(j,1)
c      ENDDO
c     switching to keplerian coordinates
c      CALL coocha(contrmin,'CAR',gms,contrelkep,'KEP',enne)
c      CALL coocha(contrplmin,'CAR',gmse,contrekpl,'KEP',enne)
c      write(*,*)'=========================================='
c      write(*,*)'keplerian coordinates of asteroid at MOID:'
c      write(*,*)contrelkep
c      write(*,*)'keplerian coordinates of Earth at MOID:   '
c      write(*,*)contrekpl
c      write(*,*)'=========================================='
c     ============================================================
      
c     **********************************************************
c     *******  W R I T I N G  on  C H E C K   F I L E S  *******
c     **********************************************************
      
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
c     OF MINIMA (>2) 
c     ===============================================================
      IF((weier.ne.0).and.(morse.ne.0))THEN
c     write on file CHECK.hi-min
         IF (nummin.gt.2) THEN
            write(iun4,*)'THERE ARE',nummin,' MINIMA for ',name(1:lnam)
            write(iun4,*)'===================================='
         ENDIF
      ENDIF
c     ================================================================
      
      WRITE(*,*)
      WRITE(*,*)'==================================================='
      WRITE(*,*)'###################################################'
      
      RETURN
      END
