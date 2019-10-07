c =============mult_sol==============
c
c PUBLIC ROUTINES
c               mmulti 
c  abolished:
c               [fmulti]
c 
c MODULE CONTAINS
c  ROUTINES
c                propsig
c                intstep
c                difvin
c                graha1
c                weakdi
c
c
c  HEADERS
c
c mult_sol.o: comdif.h phase.h mag.h
c
c            parmul.h internal, used by calling routines 
c            trig.h  public
c            parobx.h only to dimension weights and residuals
c
c
c =======================================
c  MMULTI
c =======================================
c multiple solutions for a differential correction problem
c version with adaptive stepsize
c ===============INTERFACE========================
      SUBROUTINE mmulti(obsc,inic,ok,covc,tc,
     +     eqc,gc,cc,csinor,delnor,
     +     mc,iobs,tauc,istc,alc,dec,rmsa,rmsd,
     +     rmsmag,smag,hc,gmag,sel,
     +     iun20,eqm,gm,cm,hmu,csinom,delnom,
     +           sigma,imult,imim,imip,imi0)
      IMPLICIT NONE
c ================INPUT===========================
c ======observations====
c number of obs., time, station code, alpha, delta
      INTEGER mc,istc(mc)
      DOUBLE PRECISION tauc(mc),alc(mc),dec(mc)
c rms of alpha, delta, rms magnitudes, magnitudes (string)
      DOUBLE PRECISION rmsa(mc),rmsd(mc),rmsmag(mc)
      CHARACTER*6 smag(mc)
c selection flags, obs type
      INTEGER sel(mc),iobs(mc)
c initial conditions: epoch, elements, abs. magnitude, gmag
      DOUBLE PRECISION tc,eqc(6),hc,gmag
c normal and covariance matrices, norms of residuals, of last correction
      DOUBLE PRECISION gc(6,6),cc(6,6),csinor,delnor
c output unit
      INTEGER iun20
c logical flags
      LOGICAL obsc,inic,covc
c ===============OUTPUT=========================
c number of alternate solutions, maximum
      INTEGER imult
      INCLUDE 'parmul.h'
c multiple solutions elements
      DOUBLE PRECISION eqm(6,mulx)
c normal and covariance matrices  
      DOUBLE PRECISION gm(6,6,mulx),cm(6,6,mulx)
c norms of residuals, of last correction
      DOUBLE PRECISION csinom(mulx),delnom(mulx),sigq(mulx),rescov
c magnitude for each alternate solution
      DOUBLE PRECISION hmu(mulx)
c no. of sigma along the line  
      DOUBLE PRECISION sigma
c first, last, reference solution
      INTEGER imim,imip,imi0
c success flas
      LOGICAL ok
c ==============END INTERFACE=====================
c interactive?
      INCLUDE 'comdif.h'
c weak direction, rms along it
      DOUBLE PRECISION wdir(6),wdir0(6),sdir,direc,prscag
c weighing to be computed
      INCLUDE 'parobx.h'
      DOUBLE PRECISION wc(nob2x)
c residuals
      DOUBLE PRECISION csirc(nob2x)
c no obs, used 
      INTEGER nused
c ======== differential correction flags and controls ======
c scalar temporaries
      DOUBLE PRECISION dn,sigmam
c loop indexes 
      INTEGER i, j, imi,jj
c success flag for diff. correction
      LOGICAL succ
c rms magnitudes, residuals
      DOUBLE PRECISION rmsh,resmag(nobx)
c ======== output moid =====================
      DOUBLE PRECISION moid(mulx), dnp(mulx), dnm(mulx)
      INTEGER iconv(mulx)
c hack for catalog
      INTEGER iunctc,le
      CHARACTER*9 astna0
c control of integration method
      INTEGER imint,k
      DOUBLE PRECISION hh,ratio,eqf(6)
      LOGICAL fail
c control of bizarre orbits
      LOGICAl bizarre
c verbosity control
      INCLUDE 'verbosity.h'
      INTEGER iun,iundump,iunint
      iun=iabs(iun20)
      IF(.not.batch.and.verb_mul.ge.20)THEN
         iundump=abs(iun20)
      ELSE
         iundump=-abs(iun20)
      ENDIF
      IF(.not.batch.and.verb_mul.ge.5)THEN
         iunint=abs(iun20)
      ELSE
         iunint=-abs(iun20)
      ENDIF
c =====================================================================
c check availability of observations and initial condition
      IF(.not.obsc)THEN
         WRITE(*,*)'mmulti: NO OBSERVATIONS'
         ok=.false.
         RETURN
      ENDIF
      CALL chereq(2,inic,covc,tc,iun,-1,ok)
      IF(.not.ok)THEN 
         WRITE(*,*)' mmulti: no initial data ',inic,covc,tc
         RETURN
      ENDIF
c =====================================================================
c check availability of JPL ephemerides and ET-UT table
      CALL chetim(tauc(1),tauc(mc),ok)
      IF(.not.ok) THEN
         WRITE(*,*)' mmulti: no JPL ephem/ET-UT table ',tauc(1),tauc(mc)
         RETURN
      ENDIF
c =====================================================================
c weights (removed set to zero!)
      call fitwgt(rmsa,rmsd,dec,sel,iobs,wc,mc,.true.)
c =====================================================================
c compute line of variations
      CALL weakdi(gc,wdir,sdir,iunint)
      CALL vcopy(6,wdir,wdir0)
c input parameters of segment on the variations line
      IF(batch)THEN
c sigma and imult have to be passed in the call
         IF(2*imult+1.gt.mulx)THEN
            WRITE(*,*)' too many; max is ',mulx
            STOP ' mmulti: too many mutsol'
         ENDIF
      ELSE
         WRITE(*,*)' how many sigma?'
         READ(*,*) sigma
 259     WRITE(*,*)' how many steps on each side?'
         READ(*,*) imult
         WRITE(iun20,*)' sigma=',sigma,'  in ',imult,' steps'
         IF(2*imult+1.gt.mulx)THEN
            WRITE(*,*)' too many; max is ',mulx
            GOTO 259
         ENDIF
      ENDIF
c =====================================================================
c nominal stepsize (in the sigma space)
      dn=1.d0/float(imult)
      IF(batch)THEN
c preselect a  stepsize which could not result in hyperbolic orbit
c hypothetical final point
         ratio=sqrt(2.d0)
         DO k=1,10
            DO j=1,6
               eqf(j)=eqc(j)+wdir(j)*sdir*dn*sigma*imult
            ENDDO
            IF(bizarre(eqf))THEN
               dn=dn/ratio
               IF(iun20.gt.0)WRITE(iun20,*)k,dn,eqf
            ELSE
               GOTO 9
            ENDIF
         ENDDO
 9       CONTINUE
      ENDIF
c =====================================================================
c nominal solution at the center of the list
      imi0=imult+1
      CALL vcopy(6,eqc,eqm(1,imi0))
      CALL mcopy(6,6,gc,gm(1,1,imi0))
      CALL mcopy(6,6,cc,cm(1,1,imi0))
      csinom(imi0)=csinor
      delnom(imi0)=delnor
      sigq(imi0)=0.d0
      hmu(imi0)=hc
c orbital distance
      CALL nomoid(tc,eqm(1,imi0),moid(imi0),
     +              iconv(imi0),dnp(imi0),dnm(imi0))
c =====================================================================
c main loop on number of steps (positive side)
c =====================================================================
       DO 5 i=1,imult
        imi=imult+i+1
        IF(.not.batch)WRITE(*,*)' alternate solution no. ',imi
        IF(.not.batch)WRITE(iun,*)' alternate solution no. ',imi
        CALL propsig(eqm(1,imi-1),eqm(1,imi),dn,sigma,
     +       mc,wc,sel,iobs,tauc,istc,tc,
     +           alc,dec,iundump,wdir,sdir,fail)
c check for hyperbolic
        IF(fail)THEN
           IF(.not.batch)WRITE(*,*)'step ',imi,' hyperbolic'
           IF(.not.batch)WRITE(*,*)(eqm(j,imi),j=1,6)
           imi=imi-1
           GOTO 6
        ELSEIF(bizarre(eqm(1,imi)))THEN
           IF(.not.batch)WRITE(*,*)'step ',imi,' byzarre'
           IF(.not.batch)WRITE(*,*)(eqm(j,imi),j=1,6)
           imi=imi-1
           GOTO 6
        ELSE
c differential corrections: 
           CALL difvin(mc,wc,sel,iobs,tauc,istc,tc,eqm(1,imi),
     +           alc,dec,wdir,gmag,
     +           iundump,eqm(1,imi),gm(1,1,imi),cm(1,1,imi),
     +           csinom(imi),delnom(imi),csirc,nused,succ)
c exit if not convergent
           IF(.not.succ) THEN
              imi=imi-1
              GOTO 6
           ENDIF
c estimate magnitude here
           hmu(imi)=0.d0
           CALL magest(smag,rmsmag,sel,mc,hmu(imi),resmag,rmsh)
c orbital distance
           CALL nomoid(tc,eqm(1,imi),moid(imi),
     +              iconv(imi),dnp(imi),dnm(imi))
c check for sigQ.le.sigma
           sigq(imi)=sqrt(abs(csinom(imi)**2-csinom(imi0)**2)*nused)/
     +             rescov(6,nused,csinor)
           IF(batch.and.sigq(imi).gt.sigma)THEN
              GOTO 6
           ENDIF
c compute line of variations
           CALL weakdi(gm(1,1,imi),wdir,sdir,iundump)
           direc=prscag(6,wdir,wdir0)
           IF(direc.lt.0.d0)THEN
              DO jj=1,6
                 wdir(jj)=-wdir(jj)
              ENDDO
           ENDIF
           CALL vcopy(6,wdir,wdir0)
        ENDIF
 5    CONTINUE
c =====================================================================
c recompute line of variations
 6    imip=imi
      CALL weakdi(gc,wdir,sdir,iundump)
      direc=prscag(6,wdir,wdir0)
      IF(direc.lt.0.d0)THEN
         DO jj=1,6
            wdir(jj)=-wdir(jj)
         ENDDO
      ENDIF
      CALL vcopy(6,wdir,wdir0)
      DO 7 i=1,imult
c =====================================================================
c main loop on number of steps (negative side)
c =====================================================================
        imi=imult-i+1
        IF(.not.batch)WRITE(*,*)' alternate solution no. ',imi
        IF(.not.batch)WRITE(iun20,*)' alternate solution no. ',imi
        sigmam=-sigma
        CALL propsig(eqm(1,imi+1),eqm(1,imi),dn,sigmam,
     +       mc,wc,sel,iobs,tauc,istc,tc,
     +           alc,dec,iundump,wdir,sdir,fail)
c check for hyperbolic
        IF(fail)THEN
           IF(.not.batch)WRITE(*,*)'step ',imi,' hyperbolic'
           IF(.not.batch)WRITE(*,*)(eqm(j,imi),j=1,6)
           imi=imi+1
           GOTO 8
        ELSEIF(bizarre(eqm(1,imi)))THEN
           IF(.not.batch)WRITE(*,*)'step ',imi,' byzarre'
           IF(.not.batch)WRITE(*,*)(eqm(j,imi),j=1,6)
           imi=imi+1
           GOTO 8
        ELSE
c differential corrections: 
           iun=-iun20
c for verbose        iun=iun20
           CALL difvin(mc,wc,sel,iobs,tauc,istc,tc,eqm(1,imi),
     +           alc,dec,wdir,gmag,
     +           iundump,eqm(1,imi),gm(1,1,imi),cm(1,1,imi),
     +           csinom(imi),delnom(imi),csirc,nused,succ)
c exit if not convergent
           IF(.not.succ)THEN
               imi=imi+1
               GOTO 8
           ENDIF
c estimate magnitude here
           hmu(imi)=0.d0
           CALL magest(smag,rmsmag,sel,mc,hmu(imi),resmag,rmsh)
c orbital distance
           CALL nomoid(tc,eqm(1,imi),moid(imi),
     +              iconv(imi),dnp(imi),dnm(imi))
c check for sigQ.le.sigma
           sigq(imi)=sqrt(abs(csinom(imi)**2-csinom(imi0)**2)*nused)/
     +             rescov(6,nused,csinor)
           IF(batch.and.sigq(imi).gt.sigma)THEN
              GOTO 8
           ENDIF
c compute line of variations
           CALL weakdi(gm(1,1,imi),wdir,sdir,iundump)
           direc=prscag(6,wdir,wdir0)
           IF(direc.lt.0.d0)THEN
              DO jj=1,6
                wdir(jj)=-wdir(jj)
              ENDDO
           ENDIF
           CALL vcopy(6,wdir,wdir0)
        ENDIF
 7    ENDDO
c =====================================================================
c summary table
c =====================================================================
 8    imim=imi
      IF(batch)RETURN
      CALL tee(iun20,'SUMMARY OF MULTIPLE SOLUTIONS=')
      CALL tee(iun20,
     +  'no       a      h      k      p      q      lambda=')
      CALL filopn(iunctc,'mult.ctc','unknown')
      CALL wromlh(iunctc,'ECLM','J2000')
      DO i=imim,imip
        WRITE(*,144)i,(eqm(j,i),j=1,6)
        WRITE(iun20,144)i,(eqm(j,i),j=1,6)
 144    FORMAT(i3,6f12.8)
      ENDDO
      CALL tee(iun20,'no  RMS ,lastcor,  magn,  MOID ,nod+,nod-, sigQ=')
      DO i=imim,imip
        WRITE(*,145)i,csinom(i),delnom(i),hmu(i)
     +                    ,moid(i),dnp(i),dnm(i),iconv(i),sigq(i)
        WRITE(iun20,145)i,csinom(i),delnom(i),hmu(i)
     +                    ,moid(i),dnp(i),dnm(i),iconv(i),sigq(i)
 145    FORMAT(i5,1x,1p,e13.5,e11.3,2x,0p,f5.2,1x,     
     +              f8.5,1x,f8.5,1x,f8.5,1x,i2,1x,f6.3)
        WRITE(astna0,111)i
 111    FORMAT('mult_',i4)
        CALL rmsp(astna0,le)
        CALL wromlr(iunctc,astna0(1:le),eqm(1,i),'EQU',tc,gm(1,1,i),
     +    .true.,cm(1,1,i),.true.,hmu(i),gmag,0.d0)
      ENDDO
      CALL filclo(iunctc,' ')
      RETURN
      END
c =========================================
c PROPSIG
c propagator in sigma space
      SUBROUTINE propsig(eq1,eq2,dn,sigma,
     +       mc,wc,sel,iobs,tauc,istc,tc,
     +           alc,dec,iun20,wdir,sdir,fail)
      IMPLICIT NONE
c ==============input=================
c current elements and epoch; stepsize factor, target sigma
      DOUBLE PRECISION eq1(6),tc,dn,sigma
c output unit
      INTEGER iun20
c ======observations====
c number of obs., time, station code, alpha, delta
      INTEGER mc,istc(mc)
      DOUBLE PRECISION tauc(mc),alc(mc),dec(mc)
c rms of alpha, delta, rms magnitudes, magnitudes (string)
c      DOUBLE PRECISION rmsa(mc),rmsd(mc),rmsmag(mc)
c      CHARACTER*6 smag(mc)
c selection flags, obs type
      INTEGER sel(mc),iobs(mc)
c weighing to be computed
      INCLUDE 'parobx.h'
      DOUBLE PRECISION wc(nob2x)
c ==============output================
c elements, failure flag
      DOUBLE PRECISION eq2(6)
      LOGICAl fail
c ===========end interface============
      DOUBLE PRECISION hh,ch,h,ridmax
      INTEGER iun,imint
c weak direction, rms along it
      DOUBLE PRECISION wdir(6),sdir
c intermediate point: state, covariance, normal matr., RMS, norm corr
      DOUBLE PRECISION eqt(6),g(6,6),c(6,6),csinor,delnor
c chi**2 rms for each obs,residuals, control, flag (for difcor)
      DOUBLE PRECISION x2(nobx),csi(nob2x),delcr
      INTEGER inew, icor(6),ncor,itmaxold,j
c success control is dummy
      LOGICAL succ
      INCLUDE 'comdif.h'
      INCLUDE 'verbosity.h'
c use weak direction to find initial conditions 
      hh=dn*sigma
      IF(.not.batch.and.verb_mul.ge.9)THEN
         iun=abs(iun20)
      ELSE
         iun=-abs(iun20)
      ENDIF
c 1=Euler 2= RK2
      imint=2
c iteration on stepsize reductions
c     ridmax=30.d0
      ridmax=130.d0
      h=hh
      ch=0.d0   
      CALL vcopy(6,eq1,eqt)
c try to do the step at once...
 1    CONTINUE
c      WRITE(*,*)'mmulti: ch,h, eqt ', ch,h,(eqt(j),j=1,3)
      CALL intstep(eqt,eq2,h,imint,
     +        mc,wc,sel,iobs,tauc,istc,tc,
     +        alc,dec,iun,wdir,sdir,fail)
c     write(*,*)fail,eq2
c if failed, half step
      IF(fail)THEN
         IF(abs(h).ge.abs(hh)/ridmax)THEN
            h=h/2.d0
            IF(verb_mul.ge.20)THEN
                WRITE(iun,*)'propsig: halving', h,hh
            ENDIF
            GOTO 1
         ELSE
            IF(verb_mul.ge.5)THEN
                 WRITE(iun,*)'propsig: too many halvings ',h,hh,ridmax
            ENDIF
            RETURN
         ENDIF
      ELSE
         ch=ch+h
      ENDIF
      IF(abs(ch).lt.abs(hh))THEN
         CALL vcopy(6,eq2,eqt)
c compute line of variations
c compute vectorfield at intermediate point
         itmaxold=itmax
         itmax=0
         CALL whicor(0,icor,ncor,inew)
c compute covariance and residual norm
         CALL difcor(mc,wc,sel,tc,iobs,tauc,istc,eqt,alc,dec,
     +        icor,inew,iun,delcr,
     +        eqt,g,c,csinor,delnor,csi,x2,succ)
         itmax=itmaxold
c compute weak direction and length
         CALL weakdi(g,wdir,sdir,iun) 
         GOTO 1
      ENDIF 
      IF(iun20.gt.0)WRITE(iun20,*)' propsig: stepsize ',h
      RETURN
      END 
c=================================================
c INTSTEP integration step for propagation along LOV
      SUBROUTINE intstep(eq0,eq1,hh,imint,
     +       mc,wc,sel,iobs,tauc,ioco,tc,
     +           alc,dec,iun,wdir,sdir,fail)
      IMPLICIT NONE
c stepsize, integration method
      DOUBLE PRECISION hh
      INTEGER imint
c ======observations====
c number of obs., time, station code, alpha, delta
      INTEGER mc,ioco(mc)
      DOUBLE PRECISION tauc(mc),alc(mc),dec(mc)
c rms of alpha, delta
c     DOUBLE PRECISION rmsa(mc),rmsd(mc)
c selection flags, obs type
      INTEGER sel(mc),iobs(mc)
c weighing to be computed
c weighing to be computed
      INCLUDE 'parobx.h'
      DOUBLE PRECISION wc(nob2x)
c initial conditions: epoch, elements, new elements after step
      DOUBLE PRECISION tc,eq0(6),eq1(6)
c normal and covariance matrices, norms of residuals, computed at eq0
c     DOUBLE PRECISION gc(6,6),csinoc
c output unit
      INTEGER iun
c failure flag
      LOGICAl fail
c =========end interface================
c intermediate point: state, covariance, normal matr., RMS, norm corr
      DOUBLE PRECISION eq12(6),g12(6,6),c12(6,6),csino12,delno12
      INCLUDE 'comdif.h'
c weak direction (in input computed at eq0), length of axis
      DOUBLE PRECISION wdir(6),sdir,wdir12(6),sdir12,wdirst(6),sdirst
c chi**2 rms for each obs,residuals, control, flag (for difcor)
      DOUBLE PRECISION x2(nobx),csi(nob2x),delcr
      INTEGER inew, icor(6),ncor,itmaxold
c success control is dummy,bizarre is impor
      LOGICAL succ,bizarre
c loop indexes
      INTEGER j
c controls of fixed point iteration
      INTEGER itx,it
      DOUBLE PRECISION eps,cosa,prscag,ang,dsi,e1,e12,direc
c ================================================
      IF(imint.eq.1)THEN
c Euler method
        DO  j=1,6
           eq1(j)=eq0(j)+wdir(j)*sdir*hh
        ENDDO
        fail=.false.
      ELSEIF(imint.eq.2)THEN
c Runge-Kutta-Gauss of order 2:
        itx=5
        eps=1.d-4
c first store vectorfield as it is on eq0
        CALL vcopy(6,wdir,wdirst)
        sdirst=sdir
c compute intermediate point (first guess)        
        DO j=1,6
           eq12(j)=eq0(j)+wdir(j)*sdir*hh*0.5d0
        ENDDO 
c before running differential corrections, tests that orbit is still elliptic
        e12=sqrt(eq12(2)**2+eq12(3)**2)
        IF(e12.gt.0.99d0.or.eq12(1).le.0.d0)THEN
           fail=.true.
           IF(iun.gt.0)WRITE(iun,*)'intstep: bizarre ',eq12
           call vcopy(6,eq12,eq1)
           RETURN
        ELSEIF(bizarre(eq12))THEN
           fail=.true.
           IF(iun.gt.0)WRITE(iun,*)'intstep: bizarre ',eq12
           call vcopy(6,eq12,eq1)
           RETURN
        ENDIF
c fixed point iteration
        DO it=1,itx
c compute vectorfield at intermediate point
           itmaxold=itmax
           itmax=0
           CALL whicor(0,icor,ncor,inew)
c compute covariance and residual norm
c          write(*,*)mc,iobs(1),ioco(1)
           CALL difcor(mc,wc,sel,tc,iobs,tauc,ioco,eq12,alc,dec,
     +     icor,inew,iun,delcr,
     +     eq12,g12,c12,csino12,delno12,csi,x2,succ)
           itmax=itmaxold
c compute weak direction and length
           CALL weakdi(g12,wdir12,sdir12,iun)
           direc=prscag(6,wdir12,wdir)
           IF(direc.lt.0.d0)THEN
              DO j=1,6
                 wdir12(j)=-wdir12(j)
              ENDDO
           ENDIF
c recompute intermediate point
           DO j=1,6
              eq12(j)=eq0(j)+wdir12(j)*sdir12*hh*0.5d0
           ENDDO 
c          write(*,*)'weak sigma, iter, eq12',sdir12,it,(eq12(j),j=1,3)
c before a new iteration, tests that orbit is still elliptic
           e12=sqrt(eq12(2)**2+eq12(3)**2)
           IF(e12.gt.0.99d0.or.eq12(1).le.0.d0)THEN
              IF(iun.gt.0)WRITE(iun,*)'intstep: bizarre ',eq12
              fail=.true.
              call vcopy(6,eq12,eq1)
c               write(*,*)'intstep: fail at iteration ',it,fail
              RETURN
           ELSEIF(bizarre(eq12))THEN
              IF(iun.gt.0)WRITE(iun,*)'intstep: bizarre ',eq12
              fail=.true.
              call vcopy(6,eq12,eq1)
              RETURN
           ENDIF
c convergence control
           cosa=prscag(6,wdir12,wdirst)
           IF(cosa.le.1.d0)THEN
              ang=acos(cosa)

           ELSE
              ang=0.d0
           ENDIF
           dsi=abs(sdirst-sdir12)/sdirst
c          WRITE(*,*)' intstep: it, ang, dsi, sdir12,sdirst ',
c    +          it,ang,dsi,sdir12,sdirst
           IF(ang.lt.eps.and.dsi.lt.eps)THEN
c             WRITE(*,*)' intstep: it, ang, dsir ',it,ang,dsi
c convergence achieved
              GOTO 2
           ENDIF
           CALL vcopy(6,wdir12,wdirst)
           sdirst=sdir12           
        ENDDO
c convergence failed
c       WRITE(*,*)' intstep: failed convergence, it, ang, dsig'
c       WRITE(*,*)it-1,ang,(sdir12-sdirst)/sdirst
        fail=.true.
        IF(iun.gt.0)WRITE(iun,*)'intstep: non convergent ',eq12
        call vcopy(6,eq12,eq1)
        RETURN
c convergence succeeded    
  2     CONTINUE
c compute final point         
        DO j=1,6
           eq1(j)=eq0(j)+wdir12(j)*sdir12*hh
        ENDDO 
c before a new iteration, tests that orbit is still elliptic
        e1=sqrt(eq1(2)**2+eq1(3)**2)
        IF(e1.gt.0.99d0.or.eq1(1).le.0.d0)THEN
           fail=.true.
           IF(iun.gt.0)WRITE(iun,*)'intstep: bizarre ',eq12
        ELSEIF(bizarre(eq1))THEN
           fail=.true.
           IF(iun.gt.0)WRITE(iun,*)'intstep: bizarre ',eq12
        ELSE
           fail=.false.
        ENDIF
c but we go on
      ELSE
         WRITE(*,*)' intstep: not invented yet, imint=',imint
         STOP
      ENDIF
c     write(*,*)'intstep: at exit ',fail, (eq1(j),j=1,3)
      RETURN
      END
c Copyright 1998, The Orbfit Consortium 
c Modified Jan 25 199 by A. Milani to estimate magnitudes
c To correct for new calls to: sortoss,unsort,fitwgt
c ===================================================================
c DIFVIN  linearly constrained differential corrector (single step)
c ===================================================================
c version 1.4.12, A. Milani, November 24, 1997
c        works on equinoctal elements: 
c        a,h=e sin(long.peri),k=e cos(long.peri),
c        p=tg I/2 sin(long.node),q=tg I/2 cos(long.node), lambda=mean long.
c Input: m observations number
c        w weights (only diagonal matrix)
c        sel selection flag
c        tau observations time vector
c        iobs observation type
c        ioco station codes vector
c        t0 epoch time for asteroid elements
c        eq0 asteroid equinoctal elements at time t0 (first guess)
c        al,de real observations vectors
c        peq vector to which deq must be orthogonal   
c        gmag G magnitude par for magnitude estimate     
c        iunf = unit file for output; if <0, no residuals output
c Output eq corrected orbital elements at time t0
c        gamma covariance matrix
c        gtwg inverse of covariance matrix
c           warning: if only some elements are corrected, only the 
c              corresponding entries in gamma and gtwg are nonzero 
c        csinor residuals norm
c        delnor differential corrections norm
c        csir residuals (in radians) 
c                r.a. first, then declination, for each obs
c        succ success flag; the solution is meaningful if succ=.true.
c ============= REMARK ===============================================
c The weights are constant in this routine
c =============INTERFACE===== =========================================
      SUBROUTINE difvin(m,w,sel,iobs,tau,ioco,t0,eq0,al,de,peq,gmag,
     +        iunf,eq,gamma,gtwg,csinor,delnor,csir,nused,succ)
c =====================================================================
      IMPLICIT NONE
c include files
      INCLUDE 'parobx.h'
      INCLUDE 'trig.h'
c elongation,distance to Earth, distance to Sun (to compute magnitude)
      INCLUDE 'phase.h'
c magnitude data
      INCLUDE 'mag.h'
      DOUBLE PRECISION appmag
c ================input data==========================
c no. observations, observatory codes, obs. type, selection flags
      INTEGER m, ioco(m), iobs(m),sel(m)
c  times, alpha, delta, 
      DOUBLE PRECISION tau(m),al(m),de(m)
c weights
      DOUBLE PRECISION w(2*m)
c unit file to write output
      INTEGER iunf
c epoch time, initial equinoctal elements, orthogonal vector, G 
      DOUBLE PRECISION t0, eq0(6),peq(6),gmag
c ================output ==========================
c corrected equinoctal elements 
      DOUBLE PRECISION eq(6)
c normal and covar. matrix
      DOUBLE PRECISION gtwg(6,6),gamma(6,6)
c  corr. and residuals norm
      DOUBLE PRECISION delnor,csinor
c no obs, used 
      INTEGER nused
c residuals
      DOUBLE PRECISION csir(nob2x)
c success flag
      LOGICAL succ
c =============END INTERFACE============================================
c input data sorted, in scalar form 
      DOUBLE PRECISION tsort(nobx),als(nobx),des(nobx)
      integer iocos(nobx),iposs(nobx),iocj,no
      DOUBLE PRECISION tauj,alj,dej
c weights
      DOUBLE PRECISION ws(nob2x),wsec(nob2x)
c for sorting routine 
      INTEGER selsrt(nobx),iobsrt(nobx)
c residuals, condition number,  scalar residuals
      DOUBLE PRECISION csi(nob2x),cond,ra,rd
c ====================================================================
c first and second derivatives of alpha, delta w.r. to elements
      DOUBLE PRECISION dade(6),ddde(6),ddade(6,6),dddde(6,6)
c first and second derivatives of observations
      DOUBLE PRECISION g(nob2x,6),h(nob2x,6,6),gr(nob2x,6)
c  corr. and residuals norm, their controls
      integer inew,ider,icor(6),icor6(6),nsolv,itmax,itgmax
      DOUBLE PRECISION delcr,divrat,csino1,mu,rescov,gam(6,6),c(6,6)
c differential correction, eccentricity
      DOUBLE PRECISION deq(6),ecc
c iteration indexes
      integer it,itg
c matrices for coordinate change
      DOUBLE PRECISION v(6,6),deqv(6)
c loop indexes: j=1,m; i,k=1,6
      integer j,k,i
c DOUBLE PRECISION functions
      DOUBLE PRECISION snorm,snormd
c control of two body approximation (must be false)
      logical twobo
c function for control of bizarre orbit, deserving to give up
      LOGICAL bizarre
c unit for output
      integer iun
      INCLUDE 'verbosity.h'
c ====================================================================
c number of solve-for variables
      icor(1)=0
      do  i=2,6
        icor(i)=1
      enddo
      do  i=1,6
        icor6(i)=1
      enddo
c ====================================================================
c assign ider depending from inew
      inew=2
      ider=1
c full n-body
      twobo=.false.
c sort of times and reordering of obs. and weights
      call srtoss(t0,iobs,tau,al,de,ioco,sel,w,m,iobsrt,tsort,
     +   iposs,als,des,iocos,selsrt,ws)
c Initialisation with starting value for elements
      do  k=1,6
        eq(k)=eq0(k)
      enddo
c ====================================================================
      iun=abs(iunf)
      if(iunf.gt.0)then
         write(*,220) eq
         write(iun,220) eq
 220     format(' starting values'/6f13.7)
      endif
c ================== main loop ==============================
c iteration control parameters 
      itmax=80
      itgmax=75
      delcr=1.d-3
      divrat=0.99d0 
c Loop on iterations of NEWTON's METHOD
      itg=0
      do 60 it=1,itmax
c ================== iteration ==============================
c Compute observations and derivatives
         CALL set_restart(.true.)
         do 61 j=1,m
            tauj=tsort(j)
            iocj=iocos(j)
            twobo=.false.
            IF(iobsrt(j)/1000.eq.1)THEN
               call alfdel(eq,t0,tauj,iocj,alj,dej,dade,ddde,ider,twobo,
     +              ddade,dddde)
c  compute magnitude difference (apparent minus absolute)
               dmagns(j)=appmag(0.d0,gmag,dsun,dis,pha)
            ELSEIF(iobsrt(j)/1000.eq.2)THEN
               CALL rrdot(eq,iobsrt(j),t0,tauj,iocj,alj,dej,dade,ddde,
     +              ider,twobo)
            ELSE
               WRITE(*,*)'difvin: iobs= ',iobsrt(j), ' not known'
               STOP 'difvin: iobs'
            ENDIF
            CALL set_restart(.false.)
c Compute residuals, form matrix g=-d(csi)/d(eq)
            csi(2*j-1)=als(j)-alj
            csi(2*j)=des(j)-dej
            DO i=1,6
               g(2*j-1,i)=dade(i)
               g(2*j,i)=ddde(i)
            ENDDO
 61      continue 
         CALL set_restart(.true.)
c ===========linear constraint============================
c find orthonormal basis with peq as first vector
         call graha1(peq,6,v)
c convert derivatives
         no=2*m
         call mulmat_cut(g,nob2x,no,6,v,6,6,gr)
c ===========one differential corrections step=================
c Compute solution of linear least squares
         call minsol(csi,no,ws,gr,h,inew,icor,iun,
     +        c,deqv,gam,cond)
c convert elements
         call mulmav(v,6,6,deqv,6,deq)
c Update solution (only solve-for variables) 
         do  k=1,6
            if(icor(k).ne.0)then
               eq(k)=eq(k)+deq(k)
            endif
         enddo
c ================ end iteration =========================
c Convergence control: norm of the residuals
         csinor=snormd(csi,ws,no,nused)
c reordering the residuals for output
         call unsort(iposs,m,no,csi,csir,selsrt,sel,tsort,tau)   
c unsort dmagn
         do i=1,m
            j=iposs(i)
            dmagn(j)=dmagns(i)
         enddo
c norm of the correction: the zeros do not matter!
         delnor=snorm(deqv,c,6,6)
         if(iunf.gt.0)then
            write(*,200)it,csinor,delnor,eq
            write(iun,200)it,csinor,delnor,eq
 200        format(' *** iteration ',i3,' RMS residuals =',1p,d12.4,
     +        '   norm corr =',d12.4,'  new elem values:'/0p,6f13.7/)
         endif
c control against hyperbolic and bizarre orbits
         IF(bizarre(eq))THEN
            ecc=sqrt(eq(2)**2+eq(3)**2)
            IF(verb_mul.ge.5)write(*,*)' bizarre; e=',ecc,' a=',eq(1)
            write(iun,*)' bizarre; e=',ecc,' a=',eq(1)
            succ=.false.
            return
         endif
c Check if we need another iteration
         if(delnor.lt.delcr)then
            if(iunf.gt.0)then
               IF(verb_mul.ge.9)write(*,*)' corrections small'
               write(iun,*)' convergence corrections small'
            endif
            succ=.true.
            goto 70
         endif
         if(it.gt.1)then
            if(csinor.gt.csino1*1.1d0)then
               itg=itg+1
               if(iunf.gt.0)then
                  IF(verb_mul.ge.9)
     +                 write(*,*)' target function increasing '
                  write(iun,*)' target function increasing '
               endif
               succ=.false.
            elseif(csinor.gt.csino1*divrat)then
               succ=.true.
               if(iunf.gt.0)then
                  IF(verb_mul.ge.9)write(*,*)' target funct. paralyzed '
                  write(iun,*)' target function paralyzed '
               endif
               itg=itg+1
            endif
            if(itg.gt.itgmax)then
                goto 70
            endif
         endif
         csino1=csinor
 60   continue
      succ=.true.
 70   continue
      IF(.not.succ.and.verb_mul.ge.9)WRITE(*,*)' non convergent ',it,itg
c =====================compute covariance matrix for last iteration
c Compute solution of linear least squares
      no=2*m
      call minsol(csi,no,ws,g,h,inew,icor6,iun,
     +                gtwg,deqv,gamma,cond)
c covariance (and normal matrix) rescaling
      nsolv=6
      mu=rescov(nsolv,nused,csinor)
c apply scaling to both covariance and normal matrix
      DO  i=1,6
         DO  j=1,6
            gamma(i,j)=gamma(i,j)*mu**2
            gtwg(i,j)=gtwg(i,j)/mu**2
         ENDDO
      ENDDO
c Convergence control: norm of the residuals
      IF(verb_mul.ge.9)write(*,201)it,csinor,delnor
      IF(iunf.gt.0)write(iun,201)it,csinor,delnor
 201  format(' iterations ',i3,'  RMS=',1p,d12.4,' last corr. =',d12.4)
c =====================
      return
      end
c ====================================================================
c Graham- Schmidt procedure to generate an orthonormal basis v
c starting from 1  n-vector a
c The new basis must be such that the first  vector is a
      SUBROUTINE graha1(a,n,v)
      implicit none
      integer n,nx
      parameter (nx=10)
      double precision a(n),v(n,n)
      integer j,jok,jj,i
      double precision prscag,cc,cc1,epsi,vl
      double precision ws(nx)
      logical ize
c dimension check
      if(n.gt.nx)then
         write(*,*)'n =',n,' larger than nx=',nx,' in graha'
         stop
      endif 
c selection of the control for "zero" vectors
      cc1=sqrt(prscag(n,a,a))
      epsi=1.d-12*cc1
      if(epsi.eq.0.d0)then
         write(*,*)' a has rank zero'
c        stop
      endif
c
c V1 is the versor of A
      call versor(n,a,epsi,v(1,1),vl,ize)
      if(ize)then
         write(*,*)' first vector of a is too small'
c        stop
      endif 
c we now use the vectors of the canonic basis to supplement the span of A1,A2
      jok=0
      do 1 j=1,n
c remove the components along span(A), that is along V1 
        cc1=-v(j,1)
        do  i=1,n
          ws(i)=cc1*v(i,1)
        enddo
        ws(j)=ws(j)+1.d0
        call versor(n,ws,epsi,v(1,2+jok),vl,ize)
        if(.not.ize)then
c now V(3+jok) is orthogonal to span(A); remove the components along
c the previous ones (unless it is the first)
           if(jok.gt.0)then
              do  jj=1,jok
                cc=-prscag(n,v(1,2+jok),v(1,1+jj))
                call lincog(n,v(1,2+jok),1.d0,v(1,1+jj),cc,v(1,2+jok))
              enddo
              call versor(n,v(1,2+jok),epsi,v(1,2+jok),vl,ize)
              if(ize)then
                 goto 1
              endif
           endif
c the new versor is a good one
           jok=jok+1
           if(jok.eq.n-1)then
              goto 2
           endif
        endif
 1    continue
 2    continue
      if(jok.lt.n-1)then
         write(*,*)' something went wrong, jok=',jok
      endif
      return
      end
c =====================================================================
c WEAKDI
c =====================================================================
c  weak direction and sigma
c   input:       gamma = covariance matrix
c                iun8  = output unit (if positive)
c   output:      wdir  = weak direction vector
c                sdir  = sigma
c ============ INTERFACE====================
      SUBROUTINE weakdi(gamma,wdir,sdir,iun8)
      implicit none
      double precision gamma(6,6), wdir(6),sdir
      integer iun8
c ========END INTERFACE====================
c error flag
      integer ierr
c loop index
      integer i
c eigenvalues, eigenvectors, workspace
      double precision eigvec(6,6),eigval(6),fv1(6),fv2(6)
c eigenvalues
      call rs(6,6,gamma,eigval,1,eigvec,fv1,fv2,ierr)
      call vcopy(6,eigvec(1,6),wdir)
      IF(wdir(1).lt.0.d0)THEN
         DO i=1,6
           wdir(i)=-wdir(i)
         ENDDO
      ENDIF
      sdir=sqrt(eigval(6))
      if(iun8.gt.0)then
          write(iun8,*)
          call tee(iun8,'WEAK DIRECTION =')
          write(iun8,109) (wdir(i),i=1,6)
 109      format(6e24.16)
          write(*,110) (wdir(i),i=1,6)
 110      format(1p,6e12.4)
          write(iun8,*)
          write(iun8,111) sdir
          write(*,111) sdir
 111      format(' WEAK SIGMA ',1p,e12.4)
          write(iun8,*)
      endif
      return 
      end
