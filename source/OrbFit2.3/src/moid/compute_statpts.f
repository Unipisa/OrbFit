c     ************************************************************
c     ** C O M P U T I N G   S T A T I O N A R Y    P O I N T S **
c     ************************************************************
c     ********** written by GIOVANNI F. GRONCHI (2001) ***********
c     ******* Department of Mathematics, UNIVERSITY of PISA ******
c     ============================================================
      SUBROUTINE compute_statpts(ekpl,elkep,u,upl,D2,
     *   nstat,nummin,nummax,answer,
     *   morse,weier,warnflag,sflag)
      IMPLICIT NONE
c     elements of the Earth and of the asteroid
      DOUBLE PRECISION ekpl(6),elkep(6)
c     eccentric anomalies 
      DOUBLE PRECISION u(20),upl(20)
c     SQUARED DISTANCE function 
      DOUBLE PRECISION D2(20)
c     number of stationary points (maximum = 20)
      INTEGER nstat
c     number of relative minima/maxima found
      INTEGER nummin,nummax
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
c     auxiliary variables
      DOUBLE PRECISION uu,uupl,chu,chupl
      DOUBLE PRECISION DD2 
c     ans=1:MAXIMUM; ans=-1:MINIMUM; ans=0:SADDLE; ans=2:CANNOT DECIDE;
      INTEGER ans
c
c     mutual reference system
      DOUBLE PRECISION mutI,mutom,mutompl
c     mutual orbital elements
      DOUBLE PRECISION a,e,i,om,apl,epl,ompl
      DOUBLE PRECISION beta,betapl
      COMMON/bbpl/beta,betapl
c
      DOUBLE PRECISION coe20
c     ordered 3-uples
      DOUBLE PRECISION ordu(20),ordupl(20),ordD2(20)
      INTEGER ordans(20)
c     loop index
      INTEGER j
      INCLUDE 'trig.h'
c     ======================================================================

c     flags initialization
      weier = 1
      morse = 1
      sflag(1) = 1
      sflag(2) = 1
      warnflag(1) = 1
      warnflag(2) = 1
      warnflag(2) = 1

c     mutual reference system 
         CALL mutualrefcha(ekpl,elkep,
     *        mutI,mutom,mutompl)   

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

c     ******************************************************************
c           ******************** C H E C K S ********************
c     ******************************************************************
c     ============= CHECK WITH MORSE THEORY ===============
      IF (2*(nummax+nummin).ne.nstat) THEN
         morse = 0
      ENDIF
c      WRITE(*,*)'MORSE=',morse
c     =====================================================
c     ============ CHECK WITH WEIERSTRASS THEOREM =========
      IF((nummax.lt.1).or.(nummin.lt.1))THEN
         weier = 0
      ENDIF
c      WRITE(*,*)'WEIER=',morse
c     ===============================================================

      RETURN
      END
