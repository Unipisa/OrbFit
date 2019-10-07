c     ************************************************************
c     ***** C O M P U T I N G    M I N I M U M   P O I N T S *****
c     ************************************************************
c     ********** written by GIOVANNI F. GRONCHI (2001) ***********
c     ******* Department of Mathematics, UNIVERSITY of PISA ******
c     ============================================================
      SUBROUTINE compute_minima(car,carpl,iplam,cmin,cplmin,D2,nummin)
      IMPLICIT NONE
c     elements of the planet and of the asteroid (cartesian coordinates)
      DOUBLE PRECISION car(6),carpl(6)
      INTEGER iplam
      DOUBLE PRECISION cmin(6,20),cplmin(6,20)
c     SQUARED DISTANCE function 
      DOUBLE PRECISION D2(20)
c     number of relative minima found
      INTEGER nummin
c     type of singular point : answer(j) =  1    MAXIMUM      
c                              answer(j) = -1    MINIMUM
c                              answer(j) =  0    SADDLE
c                              answer(j) =  2    CANNOT DECIDE
      INTEGER answer(20)
c     check with morse theory (morse = 1 OK ;  morse = 0 ERROR)
      INTEGER morse
c     check with weierstrass theory (weier = 1 OK ;  weier = 0 ERROR)
      INTEGER weier
c        warnflag(1) = 1: OK ; 
c        warnflag(1) = 0: leading coefficient of resultant is very small; 
c        warnflag(2) = 1: OK ; 
c        warnflag(2) = 0: higher degree terms in resultant are not small; 
c        warnflag(3) = 1: OK ; 
c        warnflag(3) = 0: low precision for resultant coefficients; 
      INTEGER warnflag(3)
c     solving system flags
      INTEGER sflag(2)
c     --------------------------------- end interface ---------------------
c     number of relative maxima found
      INTEGER nummax
c     number of stationary points (maximum = 20)
      INTEGER nstat
c     temporary elements of the Earth and of the asteroid (cartesian coord.)
      DOUBLE PRECISION cmintmp(6),cplmintmp(6)
c     elements of the Earth and of the asteroid (keplerian elements)
      DOUBLE PRECISION ekpl(6),elkep(6)
c     eccentric anomalies 
      DOUBLE PRECISION u(20),upl(20),umin(20),uplmin(20)
c     mean anomalies (asteroid and planet)
      DOUBLE PRECISION lmin(20),lplmin(20)
c     auxiliary variables
      DOUBLE PRECISION uu,uupl,chu,chupl
      DOUBLE PRECISION DD2 
c     ans=1:MAXIMUM; ans=-1:MINIMUM; ans=0:SADDLE; ans=2:CANNOT DECIDE;
      INTEGER ans
c     mutual reference system
      DOUBLE PRECISION mutI,mutom,mutompl
c     mutual orbital elements
      DOUBLE PRECISION a,e,i,om,apl,epl,ompl
      DOUBLE PRECISION beta,betapl
      COMMON/bbpl/beta,betapl
c     
      DOUBLE PRECISION enne
      DOUBLE PRECISION coe20
c     ordered 3-uples
      DOUBLE PRECISION ordu(20),ordupl(20),ordD2(20)
      INTEGER ordans(20)
c     loop index
      INTEGER k,j
      INCLUDE 'trig.h'
c     INCLUDE 'pldata.h'
c     data for masses
      INCLUDE 'jplhdr.h'
c     masses
      INCLUDE 'sunmass.h'
      INCLUDE 'parbep.h'
      INCLUDE 'masses.h'
      DOUBLE PRECISION gmsp
c     ======================================================================

c     flags initialization
      weier = 1
      morse = 1
      sflag(1) = 1
      sflag(2) = 1
      warnflag(1) = 1
      warnflag(2) = 1
      warnflag(2) = 1
      
c     mass of planet plus mass of the sun
      gmsp=gms+gm(iplam)
c      gmse = 0.000295913108d0

c      write(*,*)'gms,gmse:',gms,gmse
      
c     switching to keplerian coordinates
      CALL coocha(car,'CAR',gms,elkep,'KEP',enne)
      CALL coocha(carpl,'CAR',gmse,ekpl,'KEP',enne)

c     check
c      WRITE(*,*)'asteroid, cartesian:',car
c      WRITE(*,*)'Earth, cartesian:',carpl
c      WRITE(*,*)'asteroid, Keplerian:',elkep
c      WRITE(*,*)'Earth, Keplerian:',ekpl
      
c     mutual reference system 
      CALL mutualrefcha(ekpl,elkep,mutI,mutom,mutompl)   
      
c     rename variables in mutual reference frame
      a    = elkep(1)
      e    = elkep(2)
      apl  = ekpl(1)
      epl  = ekpl(2)
      i    = mutI
      om   = mutom
      ompl = mutompl
      
      beta=dsqrt(1-e**2)
      betapl=dsqrt(1-epl**2)

c     ====================================================================
c     COMPUTING FUNCTIONS of the CONSTANT ORBITAL ELEMENTS 
c     ====================================================================
      CALL aical(a,e,i,om,apl,epl,ompl,coe20)
c     write(*,*)'elements',a,e,i,om,apl,epl,ompl
      
c     ===================================================================
c     COMPUTING STATIONARY POINTS FOR D2 AND DECIDING WHICH ARE 
c     RELATIVE MAXIMA/MINIMA BETWEEN THEM
c     ===================================================================
      CALL comp_heart(u,upl,nstat,nummin,nummax,answer,
     *     warnflag,sflag)
      
c     loop on number of stationary points 
      DO j = 1,nstat
         
c     conversion in radians for D2eval
         uu=u(j)*radeg
         uupl=upl(j)*radeg   
         
c     ==================================================================
c     COMPUTE SQUARED DISTANCE in the points u,upl 
c     ==================================================================     
         CALL D2eval(uu,uupl,DD2)
c     ==================================================================
         
         D2(j)=DD2
         
c     angles between 0 and 360 degrees
         CALL choosedeg(u(j),chu)
         CALL choosedeg(upl(j),chupl)
         
         u(j) = chu
         upl(j) = chupl
         
c     END MAIN DO LOOP
      ENDDO
      
c     =====================================================
c     sorting 3-uples (u,upl,D2) according to value of D2 
      CALL sortd2(nstat,u,upl,D2,answer)
c     =====================================================
      
      DO k = 1,nummin
         umin(k) = u(k)
         uplmin(k) = upl(k)
c     check
c     write(*,*)'umin(',k,')=',umin(k)
c     write(*,*)'uplmin(',k,')=',uplmin(k)
         lmin(k) = u(k)*radeg - elkep(2)*dsin(u(k)*radeg)
         lplmin(k) = upl(k)*radeg - ekpl(2)*dsin(upl(k)*radeg)
         elkep(6) = lmin(k)
         ekpl(6) = lplmin(k)
c     switching to cartesian coordinates
         CALL coocha(ekpl,'KEP',gmse,cplmintmp,'CAR',enne)
         CALL coocha(elkep,'KEP',gms,cmintmp,'CAR',enne)
         DO j = 1,6
            cmin(j,k) = cmintmp(j)
            cplmin(j,k) = cplmintmp(j)
c     check
c            write(*,*)'j,k,cmin,cplmin',j,k,cmin(j,k),cplmin(j,k)
         ENDDO
      ENDDO
      
c     *****************************************************
c     ******************** C H E C K S ********************
c     *****************************************************
c     ============= CHECK WITH MORSE THEORY ===============
      IF (2*(nummax+nummin).ne.nstat) THEN
         morse = 0
      ENDIF
c     WRITE(*,*)'MORSE=',morse
c     =====================================================
c     ============ CHECK WITH WEIERSTRASS THEOREM =========
      IF((nummax.lt.1).or.(nummin.lt.1))THEN
         weier = 0
      ENDIF
c     WRITE(*,*)'WEIER=',morse
c     ===============================================================
      
      RETURN
      END
