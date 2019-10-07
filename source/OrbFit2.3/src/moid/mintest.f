c     ********************************************************************
c     *******************  PROGRAM   M I N T E S T  **********************
c     ********************************************************************
c     ************* written by GIOVANNI F. GRONCHI (2001) ****************
c     ********************************************************************
      PROGRAM mintest
      IMPLICIT NONE
c     elements of the Earth and of the asteroid (cartesian coordinates)
      DOUBLE PRECISION car(6),carpl(6)
c     cartesian coordinates at minima
      DOUBLE PRECISION cmin(6,20),cplmin(6,20)
c     check variables
      DOUBLE PRECISION contrmin(6),contrplmin(6)
c     
      DOUBLE PRECISION enne
c     SQUARED DISTANCE function 
      DOUBLE PRECISION D2(20)
c     
      INTEGER nummin
c     ========================= check flags ==============================
c     morse:       1 = "OK!" ;   0 = "CHECK with MORSE FAILED!"
c     weier:       1 = "OK!" ;   0 = "CHECK with WEIERSTRASS FAILED!"
c     warnflag(1): 1 = "OK!" ;   0 = "lead. coef. of resultant = 0!"
c     warnflag(2): 1 = "OK!" ;   0 = "higher deg. terms not small!"
c     warnflag(3): 1 = "OK!" ;   0 = "low precision in res. coef. (fft)!"
c     sflag(1):    1 = "OK!" ;   0 = "WARNING! They are both good solutions!"
c     sflag(2):    1 = "OK!" ;   0 = "WARNING! Neither of two is close to 0!"
      INTEGER morse,weier,warnflag(3),sflag(2)
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

c     data for (433) Eros
      car(1) =  0.75307703 
      car(2) =  1.01818502
      car(3) =  0.228917553 
      car(4) = -0.0143138666 
      car(5) =  0.0070586288
      car(6) = -0.00149599397

      carpl(1) = -0.979829992
      carpl(2) = -0.19614698
      carpl(3) =  1.54932952E-06 
      carpl(4) =  0.00310409504
      carpl(5) = -0.016932854
      carpl(6) = -5.72122365E-07

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
         WRITE(*,*)'minimum label =',j,'    Distance=',dsqrt(D2(j))
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
      
      
c     ============= CHECK WITH MORSE THEORY ===============
      IF (morse.eq.0) THEN
c     writing on screen
         WRITE(*,*)'ERROR: MORSE THEORY VIOLATION!!!'
      ENDIF
c     =====================================================
      
c     ============ CHECK WITH WEIERSTRASS THEOREM =========
      IF(weier.eq.0)THEN
c     writing on screen
         WRITE(*,*)'ERROR: WEIERSTRASS THEOREM VIOLATION!!!'
      ENDIF
c     ===============================================================
            
      WRITE(*,*)
      WRITE(*,*)'==================================================='
      WRITE(*,*)'###################################################'
              
c     STOP
      END
