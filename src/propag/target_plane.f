c ========MODULE target_plane==================
c CONTAINS
c PUBLIC
c             v_infty  vel at infinity (w.r. to Earth only) public 
c             fclan2   close approach analysis driver
c ROUTINES
c             mtpsel2
c             newton_tp
c             tpslin
c             virimp
c             aftclov
c             preob4 move to semilin
c             linobs4  "
c             ellip4   "
c             elemo4   "
c             slinel4  "
c             [mtpref]
c             [mtprot]
c  
c HEADERS
c target_plane.o: npoint.h comdif.h phase.h 
c               
c             closapl.h with close.app.f
c             cloapp.h        "
c             tpcana.h        "
c             sunmass.h trig.h public
c             parbep.h masses.h names and disk radius, masses
c             parobx.h only to dimension weights
c
c ================================================
c velocity at infinity with repect to Earth
c (circular approx) in AU/day
      DOUBLE PRECISION FUNCTION v_infty(eq0)
      IMPLICIT NONE
c input equinoctal elements
      DOUBLE PRECISION eq0(6)
c end interface
      DOUBLE PRECISION cosi
c v_infinity computation
      cosi=(1.d0-eq0(4)**2-eq0(5)**2)/
     +     (1.d0+eq0(4)**2+eq0(5)**2)
      v_infty=sqrt(3.d0-1.d0/eq0(1) - 
     +     2.d0*sqrt(eq0(1)*(1.d0-eq0(2)**2-eq0(3)**2))*cosi)
c normalization in AU/day by Gauss constant
      v_infty=v_infty*0.01720209895d0
      RETURN
      END 
c =======================================
c  FCLAN2
c vers. 2.2.5, A. Milani, 5 april 2001
c =======================================
c close approach analysis
c ===============INTERFACE========================
      SUBROUTINE fclan2(batchcl,t1,iun20,ok,siglim,
     +     eq0,t0,g0,c0,csinor,delnor,hmag,gmag,astna0,
     +     obs0,m,objid,iobs,tau,tutm,idsta,aln,den,rmsa,rmsd,sel)
      IMPLICIT NONE
c ================INPUT===========================
c control for batch mode
      LOGICAL batchcl
c time to monitor up to, passed only in batch mode
      DOUBLE PRECISION t1
c max value of sigma change allowed for newton step
      DOUBLE PRECISION siglim
c ======observations====
c data availability
      LOGICAL obs0
c number of observations, station codes, obs. type, object id
      INTEGER m,idsta(m),iobs(m)
      CHARACTER*(9) objid(m)
c selection flags
      INTEGER sel(m)
c observation times (ET, UT), alpha, delta, a priori rms
      DOUBLE PRECISION tau(m),tutm(m),aln(m),den(m),rmsa(m),rmsd(m)
c =======orbit==============
c magnitude
      DOUBLE PRECISION hmag,gmag
c asteroid names (18 CHARACTERs)
      CHARACTER*18 astna0
c nominal initial conditions: epoch, elements
      DOUBLE PRECISION t0,eq0(6)
c normal and covariance matrices
      DOUBLE PRECISION g0(6,6),c0(6,6)
c norms of residuals, of last correction
      DOUBLE PRECISION csinor,delnor
c output unit, verbosity
      INTEGER iun20,iverb
c ===============HIDDEN INPUT========================
c arrays of close approach data and controls
      INCLUDE 'cloapp.h'
      INCLUDE 'closapl.h'
c ===============OUTPUT=========================
c success flag
      LOGICAL ok
c ==============HIDDEN OUTPUT===================
c covariance matrix passed through hidden common
c ==============END INTERFACE=====================
c target plane type, nonlinearity handling, flag for ephemerides
      INTEGER iclan,inl,ikey
c logical controls
      LOGICAL error
c menu empty items
      CHARACTER*20 menunam
      CHARACTER*70 s3,s4,s5,s6,s7,s8,s9,s10
c final elements after propagation
      DOUBLE PRECISION eq1(6),tcl,c1(6,6),g1(6,6)
c time of close approach passed to subroutines
      DOUBLE PRECISION tclo
c initialisation
      INTEGER lflag
c memory model
      SAVE
      DATA lflag /0/
      IF(lflag.eq.0)THEN
         ipla0=0
         lflag=1
      ENDIF
c ================================================================
c t1 has been passed as dummy variable, nothing to do
c initialize close approach count
      eprdot=1.d-10
      njc=0
      iplam=0
c check availability of JPL ephemerides
      CALL chetim(t0,t1,ok)
      IF(.not.ok)THEN
         WRITE(*,*)' JPL ephemerides not available for t1=',t1
         ok=.false.
         RETURN
      ENDIF
c availability of covariance
      CALL cov_avai(g0)
c explorative propagation, incorporates linear target plane analysis
      CALL proelc('EQU',t0,eq0,g0,c0,t1,eq1,g1,c1)
c select target plane
      IF(batchcl)THEN
         iverb=1
      ELSE
         iverb=20
      ENDIF
      CALL mtpsel2(batchcl,iverb,error)
      IF(error)THEN
         ok=.false.
         RETURN
      ENDIF
c =====================================================================
c selection of the procedure to use the target plane information
c ===================================================================
c restart from menu
 58   CONTINUE
c choose handling of nonlinearity
      IF(batchcl)THEN
c batch use is for Newton's method (from clonew.f)
         inl=4
      ELSE
         menunam='closnonl'
         CALL menu(inl,menunam,4,'How to handle nonlinearity?=',
     +        'linear map=',
     +        'semilinear n-body=',
     +        'Newtons method on LOV, 1 step=',
     +        'Newtons method on LOV, auto=',
     +        s5,s6,s7,s8,s9,s10)
c exit
         IF(inl.eq.0)RETURN
      ENDIF
c ===================================================================
c operation mode depending upon inl
      IF(inl.eq.1.or.inl.eq.2)THEN
c ============target plane semilinear (and linear) analysis:
         CALL tpslin(batchcl,inl,iun20,astna0,t0,eq0,g0,c0,ok)
      ELSEIF(inl.eq.3.or.inl.eq.4)THEN
c ============Newton's method to find virtual impactors:
         IF(.not.batchcl)THEN
            WRITE(*,*)' assign limit to change in sigma, 0=default'
            READ(*,*)siglim
            IF(siglim.eq.0.d0)siglim=3.d0
         ENDIF
         tclo=tcla(jcsel)
         CALL newton_tp(batchcl,inl,iun20,siglim,ok,
     +        tclo,astna0,t0,eq0,hmag,gmag,g0,c0,csinor,delnor,
     +     m,objid,iobs,tau,tutm,idsta,aln,den,rmsa,rmsd,sel)
c
         IF(batchcl)THEN
            RETURN
         ENDIF
      ENDIF
c ===========================================================
c back to the menu
      GOTO 58
      END
c ====================================
c  MTPSEL2 
c selects target plane (in case there are many)
c also handles output in interactive case
c in batch, it selects the lowest minimum distance
c ===================================== 
      SUBROUTINE mtpsel2(batchcl,verbose,error)
      IMPLICIT NONE
c interface
c input: batch? write?
      LOGICAL batchcl
      INTEGER verbose
c output: something wrong
      LOGICAL error
c end interface
c hidden input:
c arrays of close approach data and controls
      INCLUDE 'cloapp.h'
c basis adapted to MTP output from mtprot
      INCLUDE 'tpcana.h'
c planetary data (names)
      INCLUDE 'parbep.h'
      INCLUDE 'masses.h'
c planet names
      CHARACTER*30 planam
      INTEGER lpla,lench
c index of target planes, indexes of coordinates on tp
      INTEGER jc,i,j,jj,k
c minimum distance, angle on TP
      DOUBLE PRECISION dist_min,tpth
c are there close approaches?
      IF(iplam.eq.0.or.njc.eq.0)THEN
         WRITE(*,*)' no close approaches found', iplam, njc
c error flag
         error=.true.
         RETURN
      ELSE
         error=.false.
         planam=ordnam(iplam)
         lpla=lench(planam)
         IF(.not.batchcl)THEN
            WRITE(*,101)planam(1:lpla)
 101        FORMAT(' close approach to planet ',a)
         ENDIF
      ENDIF
c selection of minimum
      IF(.not.batchcl.or. verbose.gt.19)THEN
         WRITE(*,102)(jc,tcla(jc),rmin(jc),jc=1,njc)
 102     FORMAT(' time(MJD), min.dis.(AU), normal to MTP'/
     +        (i3,1x,f12.5,1x,f11.8/))
      ENDIF
      IF(njc.eq.1)THEN
         IF(batchcl)THEN
            jcsel=1
         ELSE
 198        WRITE(*,*)' select target plane 1=yes, 0=quit'
            READ(*,*)jcsel
            IF(jcsel.eq.0)THEN
               WRITE(*,*)' no target plane analysis'
               error=.true.
               RETURN            
            ELSEIF(jcsel.lt.0.or.jcsel.gt.njc)THEN
               WRITE(*,*) 'should be 0< ', jcsel,' < ',njc+1
               GOTO 198
            ENDIF
         ENDIF
      ELSE
         IF(batchcl)THEN
c selects the minimum; is this right???
            dist_min=1.d2
            DO jc=1,njc
              IF(tpr(jc).lt.dist_min)THEN
                 dist_min=tpr(jc)
                 jcsel=jc
              ENDIF
            ENDDO
         ELSE
 199        WRITE(*,*)' select target plane among the above,0=quit'
            READ(*,*)jcsel
            IF(jcsel.eq.0)THEN
               WRITE(*,*)' no target plane analysis'
               error=.true.
               RETURN
            ELSEIF(jcsel.lt.0.or.jcsel.gt.njc)THEN
               WRITE(*,*) 'should be 0< ', jcsel,' < ',njc+1
               GOTO 199
            ENDIF
         ENDIF
      ENDIF
c output (only in interactive mode)
      IF(verbose.lt.9) RETURN
      tpth=atan2(tpc(2,jcsel),tpc(1,jcsel))
      WRITE(*,170) (tpc(j,jcsel),j=1,3),tpr(jcsel),tpth
 170  FORMAT('MTP coordinates ',3f12.8,' dist ',f12.8,
     +     ' angle ',f8.4)
      WRITE(*,*)' partial derivatives'
      DO k=1,3
        WRITE(*,171) (dtpdet(jj,k,jcsel),jj=1,6)
      ENDDO
 171  FORMAT(6(f11.3,1x))
      WRITE(*,173)(sig(i,jcsel),i=1,2)
 173  FORMAT(' semiaxes of target ellipse: ',1p,2d10.3)
      WRITE(*,172)((axes(i,j,jcsel),i=1,2),j=1,2)
 172  FORMAT(' directions ',2f10.6/12x,2f10.6)
      RETURN
      END
c =================================
c NEWTON
c called by fclan2
c target plane Newton's method to find minimum distance
c on long axis of TP ellipse
c ==================================
      SUBROUTINE newton_tp(batchcl,inl,iun20,siglim,ok,
     +     tclo,astna0,t0,eq0,hmag,gmag,g0,c0,csinor,delnor,
     +     m,objid,iobs,tau,tutm,idsta,aln,den,rmsa,rmsd,sel)
      IMPLICIT NONE
c INPUT
c batch control
      LOGICAL batchcl
c inl=1 linear inl=2 semilinear
      INTEGER inl
c log unit
      INTEGER iun20
c max value of sigma change allowed for newton step
      DOUBLE PRECISION siglim
c time of close approach being explored
      DOUBLE PRECISION tclo
c =========elements=================
c asteroid names (18 CHARACTERs)
      CHARACTER*18 astna0
c nominal initial conditions: epoch, elements, magnitude
      DOUBLE PRECISION t0,eq0(6),hmag,gmag
c normal and covariance matrices
      DOUBLE PRECISION g0(6,6),c0(6,6)
c residuals norm (for nominal orbit), convergence of last diff. corr.
      DOUBLE PRECISION csinor,delnor
c =========observations=================
c number of observations, station codes, obs. type, object id
      INTEGER m,idsta(m),iobs(m)
      CHARACTER*(*) objid(m)
c selection flags
      INTEGER sel(m)
c observation times (ET, UT), alpha, delta, a priori rms
      DOUBLE PRECISION tau(m),tutm(m),aln(m),den(m),rmsa(m),rmsd(m)
c HIDDEN INPUT
c max distance (from nominal point on target plane) at which we can 
c reasonably apply a linear approx
      DOUBLE PRECISION linlim
c     PARAMETER (linlim=0.03d0)
      PARAMETER (linlim=0.02d0)
c arrays of close approach data 
      INCLUDE 'cloapp.h'
c basis adapted to MTP output from mtprot
      INCLUDE 'tpcana.h'
c OUTPUT
c success?
      LOGICAL ok
c ==================END INTERFACE=====================
c ====== new starting point=================
c correspondent ellipse in orbital elements space
      DOUBLE PRECISION b(2,2),ceicel(4,2),v(6,6)
c variables for Newton's method
      DOUBLE PRECISION s1,tpc1(2),dtpc(2),s1max
      LOGICAL s1_change
      DOUBLE PRECISION del2(2),del4(4),dels(6),delem(6),elem(6)
c new expected min.distance, normal and for shortened step,change in distance
      DOUBLE PRECISION dminew,dminew2,deldis
c ================for difcor=================
      DOUBLE PRECISION csinew,delnew,cnew(6,6),gnew(6,6)
      INTEGER icor(6),inew,inter,ncor,itsav,iun20m,jcselold
      DOUBLE PRECISION delcr
      LOGICAL succ,error
      INCLUDE 'comdif.h'
      INCLUDE 'parobx.h'
      DOUBLE PRECISION w(nob2x),csir(nob2x),x2(nobx)
c =============== for new propagation===============
c vel at infinity
      DOUBLE PRECISION v_inf,v_infty
c time of the final elements, new elements, covariance
c time too early for close app. in the same shower
      DOUBLE PRECISION tafter,eq1(6),c1(6,6),g1(6,6),tbefore
c newton iterations
      INTEGER itnewma,it,it_too_long, itlongmax,iverb
      DOUBLE PRECISION newcont, disnow
c  moid routine
      DOUBLE PRECISION moid,dnp,dnm
      INTEGER iconv
c ephem of VI
      INTEGER ikey
      CHARACTER*18 astnavi
c =================================
c loop indexes
      INTEGER jj
c =============================================
c Newton's method setup
      IF(.not.batchcl)WRITE(*,*)' Newton method'
      ok=.false.
      IF(inl.eq.3)THEN
         itnewma=1
      ELSEIf(inl.eq.4)THEN
         itnewma=10
      ENDIF
      itlongmax=5
      newcont=1.d-6
      it_too_long=0
c iteration loop
      DO 1 it=1,itnewma
c         IF(.not.batchcl)
         WRITE(*,*)' ===== iteration ',it,' ========'
         disnow=sqrt(tpc(1,jcsel)**2+tpc(2,jcsel)**2)
c find closest point on MTP according to linear approx
         s1=-tpc(1,jcsel)*axes(1,2,jcsel)-tpc(2,jcsel)*axes(2,2,jcsel)
c new target point on MTP
         DO jj=1,2
            dtpc(jj)=s1*axes(jj,2,jcsel)
            tpc1(jj)=tpc(jj,jcsel)+dtpc(jj)
         ENDDO
         dminew=sqrt(tpc1(1)**2+tpc1(2)**2)
c         IF(.not.batchcl) 
         WRITE(*,120) s1,s1/sig(2,jcsel),dminew
 120     FORMAT(' displ.=',1p,d9.2,' deltasig=',d9.2,' dmin=',d9.2)
         WRITE(iun20,*) s1,s1/sig(2,jcsel),dminew
c control on the size of displacement s1 (in AU)
         s1_change=.false.
         s1max=linlim
         IF(abs(s1).gt.s1max)THEN
c            IF(.not.batchcl)
            WRITE(*,121)s1,s1max
 121        FORMAT(' Newton step too long s1=',1p,d9.2,' max=',d9.2)
            s1=s1*s1max/abs(s1)
            s1_change=.true.
         ENDIF
c control on the size of change in the sigma space
         IF(abs(s1/sig(2,jcsel)).gt.siglim)THEN
            WRITE(*,122)s1/sig(2,jcsel),siglim
 122        FORMAT(' too large deltasig=',1p,d9.2,' max=',d10.2)
c           ok=.false.
c           RETURN
c experiment!!!
            it_too_long=it_too_long+1
            s1=s1*siglim/abs(s1/sig(2,jcsel))
            s1_change=.true. 
        ENDIF
c new target point on MTP
        IF(s1_change)THEN
           DO jj=1,2
              dtpc(jj)=s1*axes(jj,2,jcsel)
              tpc1(jj)=tpc(jj,jcsel)+dtpc(jj)
           ENDDO
           dminew2=sqrt(tpc1(1)**2+tpc1(2)**2)
c           IF(.not.batchcl) 
           WRITE(*,120) s1,s1/sig(2,jcsel),dminew2
           WRITE(iun20,120) s1,s1/sig(2,jcsel),dminew
        ENDIF
c compute ellipse in the elements space 
        CALL slinel(dtpdet(1,1,jcsel),g0,c0,ceicel,b,v)
c find corresponding orbital elements at epoch
        CALL mulmav(b,2,2,dtpc,2,del2)
        CALL mulmav(ceicel,4,2,del2,2,del4)
        CALL vcopy(2,del2,dels)
        DO jj=1,4
           dels(2+jj)=-del4(jj)
        ENDDO
        CALL mulmav(v,6,6,dels,6,delem)
        CALL vsumg(6,eq0,delem,elem)
c setup of weights (non-zero!)
        CALL fitwgt(rmsa,rmsd,den,sel,iobs,w,m,.false.)
c solve for all elements
        inter=0
        CALL whicor(inter,icor,ncor,inew)
c zero iterations differential corrections, batch mode
        itsav=itmax
        itmax=0
        iun20m=-iun20
c        batch=.true.
        CALL difcor(m,w,sel,t0,iobs,tau,idsta,elem,aln,den,icor,inew,
     +       iun20m,delcr,elem,gnew,cnew,csinew,delnew,csir,x2,succ)
        itmax=itsav
c time up to which to search for close approaches
        v_inf=v_infty(elem)
        CALL aftclov(iplam,t0,tclo,v_inf,tbefore,tafter)
c make new covariance matrix available for target plane analysis
        CALL cov_avai(gnew)
c reset storage of close encounters
        njc=0
        iplam=0
c propagation to search for new target plane point
        CALL proelc('EQU',t0,elem,gnew,cnew,tafter,eq1,g1,c1)
c select target plane: the same local minimum!!!????
        jcselold=jcsel
c       IF(batchcl)THEN
c          iverb=1
c       ELSE
           iverb=20
c       ENDIF
        CALL mtpsel2(.true.,iverb,error)
c different causes of failure
        IF(error)THEN
           ok=.false.
           WRITE(*,*)' No target plane, Newton failed iter ',it
           RETURN
        ELSEIF(jcselold.ne.jcsel)THEN
           WRITE(*,*)' change of jcsel, was ',jcselold,' now ',jcsel
c should we stop????
c          WRITE(*,*)' Newton failed at iteration ',it
c          ok=.false.
c          RETURN
        ELSEIF(tcla(jcsel).gt.max(tbefore,tafter).or.
     +          tcla(jcsel).lt.min(tbefore,tafter))THEN
           WRITE(*,130)tcla(jcsel), tbefore, tafter
 130       FORMAT(' time_cla=',f11.2,' out of range, must be in '
     +          ,2f11.2)
           WRITE(*,*)' Newton failed iter ',it
           ok=.false.
           RETURN
        ENDIF
c if there is an acceptable close approach, adopt the newton orbit
c as nominal
        CALL vcopy(6,elem,eq0)
        CALL mcopy(6,6,gnew,g0)
        CALL mcopy(6,6,cnew,c0)
        csinor=csinew
        delnor=delnew
c assess where we stand after this iteration
        deldis=tpr(jcsel)-disnow
        disnow=tpr(jcsel)
        write(*,*)' Newton controls=',dminew-disnow,deldis
        write(iun20,*)' Newton controls=',dminew-disnow,deldis
        IF(abs(dminew-disnow).lt.newcont.or.abs(deldis).lt.newcont)THEN  
           ok=.true.  
           IF(batchcl)THEN
c Output the data for resret.pl:
c($time,$mindist,$tpc[0],$tpc[1],$semiwidth,$stretch,$rms,$nconv)=
c            write (iun20,145) tcla(1),rmin(1),tpc(1),tpc(2),
c     +           sig(1),sig(2), csinor, dminew-tpr
c            write (*,145) tcla(1),rmin(1),tpc(1),tpc(2),
c     +            sig(1),sig(2),csinor,dminew-tpr
c assuming reference system on MTP has first axis along the  
c direction of the last close approach found
              write (iun20,145) tcla(jcsel),disnow,disnow,0.d0,
     +             sig(1,jcsel),sig(2,jcsel), csinor, dminew-disnow
              write (*,145) tcla(jcsel),disnow,disnow,0.d0,
     +             sig(1,jcsel),sig(2,jcsel),csinor,dminew-disnow
 145          format(f14.7,7(1x,g18.10))
           ELSE  
c convergence already achieved
               write(*,*)' newton converged'
               astnavi='virimp'
               CALL wromlr (iun20,astnavi,eq0,'EQU',t0,g0,.true.,
     +                  c0,.true.,hmag,gmag,0.d0)
               CALL nomoid(t0,eq0,moid,iconv,dnp,dnm)
               write(iun20,198)moid,iconv,dnp,dnm
 198           format('!MOID ',f8.5,1x,i4/'!NODES ',f8.5,1x,f8.5)
               write(*,*)' ephemeris of virtual impactor? 1=yes 0=no'
               read(*,*)ikey
               IF(ikey.eq.1)THEN
                  CALL virimp(tpc(1,jcsel),dtpdet(1,1,jcsel),
     +                 axes(1,1,jcsel),sig(1,jcsel),ceicel,b,v,
     +                 t0,eq0,g0,c0,hmag,gmag,iun20)
               ENDIF
           ENDIF
           GOTO 2
        ELSEIF(it_too_long.ge.itlongmax)THEN
           GOTO 3
        ENDIF
c end iteration loop
 1    ENDDO
      IF(.not.batchcl.and.inl.eq.4)THEN
         WRITE(*,*)' Done ',itmax,' Newton iterations'
      ENDIF
c     controversial choice: if Newton has not converged, but not gone
c     out of target plane, then take where it stops as output
      ok=.true.
      RETURN
 3    IF(.not.batchcl.and.inl.eq.4)THEN
         WRITE(*,*)' Done ',it,' Newton iterations ', it_too_long,
     +        ' too long'
      ENDIF
c     controversial choice: if Newton has stopped because of too many
c     short steps, but not gone out of target plane, then take where it
c     stops as output 
      ok=.true.
      RETURN
 2    IF(.not.batchcl)WRITE(*,*)' Newton convergence at iter. ',it
      RETURN
      END
c =================================
c TPSLIN
c called by fclan2
c target plane semilinear (and linear) analysis:
c ==================================
      SUBROUTINE tpslin(batchcl,inl,iun20,astna0,t0,eq0,gc,cc,ok)
      IMPLICIT NONE
c INPUT
c batch control
      LOGICAL batchcl
c inl=1 linear inl=2 semilinear
      INTEGER inl
c log unit
      INTEGER iun20
c asteroid names (18 CHARACTERs)
      CHARACTER*18 astna0
c nominal initial conditions: epoch, elements
      DOUBLE PRECISION t0,eq0(6)
c normal and covariance matrices
      DOUBLE PRECISION gc(6,6),cc(6,6)
c HIDDEN INPUT
c arrays of close approach data 
      INCLUDE 'cloapp.h'
c basis adapted to MTP output from mtprot
      INCLUDE 'tpcana.h'
c OUTPUT
c success?
      LOGICAL ok
c END INTERFACE
c propagations: v_inf, close approach times, elems after 
      DOUBLE PRECISION v_inf,v_infty,tbefore,tafter,tcla0,eq1(6)
c ====== multiple target plane intersections============
c string of virtual asteroids: elements, target plane points
      INTEGER npox
      PARAMETER (npox=4000)
      DOUBLE PRECISION elm(6,npox),xcl(npox),ycl(npox)
c confidence boundary, line of max variation
      INTEGER ibv,npo,npo1,npoc,nc
      DOUBLE PRECISION sigma,maxsig,minsig
c correspondent ellipse in orbital elements space
      DOUBLE PRECISION b(2,2),ceicel(4,2),v(6,6)
c target plane reference system
c      DOUBLE PRECISION v1(3),v2(3),vv3(3,3),vvt3(3,3),xc(3),vc(3)
c    +  ,vc0,vsize,prscag,dxl,dyl
      DOUBLE PRECISION xc(3),vc(3),prscag,dxl,dyl
c index of planet (for comparison)
      INTEGER iplamold
c output units
      INTEGER iun7,iun8,iun9
c device flag,labels for plots
      INTEGER idev
      CHARACTER*60 ylabel,xlabel,title
c loop indexes
      INTEGER i,n
c==========================================
c input specification of set of points
      CALL asscbd(iun20,npox,npo,sigma,ibv)
c If ibv=0 then use automatic selection method
      IF(ibv.eq.0)THEN
         maxsig=max(sig(1,jcsel),sig(2,jcsel))
         minsig=min(sig(1,jcsel),sig(2,jcsel))
         if(maxsig/minsig.le.200.d0)then
            ibv=1
         else
            ibv=2
         endif
      endif    
c compute ellipse in the elements space 
      CALL slinel(dtpdet(1,1,jcsel),gc,cc,ceicel,b,v)
c compute line of orbital elements
      CALL linobs(ibv,npo,eq0,axes(1,1,jcsel),sig(1,jcsel),
     +        b,v,sigma,ceicel,elm,npo1)
c open output files for graphics
      IF(inl.eq.1)THEN
c linear analysis is by definition on a fixed TP
         tcla0=tcla(jcsel)
         CALL filopn(iun7,'mtp.fla','unknown')
      ELSEIF(inl.eq.2)THEN
c semilinear analysys is on a variable TP, the close approach manifold
         CALL filopn(iun8,'clo.fla','unknown')
c also resonant return analysys could be done by looking at the elements
         CALL filopn(iun9,'a.fla','unknown')
c need to find tolerances on target plane time
         v_inf=v_infty(eq0)
         tcla0=tcla(jcsel)
         CALL aftclov(iplam,t0,tcla0,v_inf,tbefore,tafter)
      ELSE
         WRITE(*,*)'tpslin: option inl=',inl,' unknown'
         RETURN
      ENDIF
c ===========================================================
c main loop on the number of output points
      nc=0
      iplamold=iplam
      DO 7 n=1,npo1
        IF(inl.eq.1)THEN
c linear map from ellipse
           dxl=prscag(6,elm(1,n),dtpdet(1,1,jcsel))
           dyl=prscag(6,elm(1,n),dtpdet(1,2,jcsel))
           nc=nc+1
           xcl(nc)=dxl+tpc(1,jcsel)
           ycl(nc)=dyl+tpc(2,jcsel)
c file output: target plane linear map
           WRITE(iun7,107)xcl(nc),ycl(nc)
 107       FORMAT(7e20.12)
           CALL vsumg(6,eq0,elm(1,n),elm(1,n))
        ELSEIF(inl.eq.2)THEN
c full n-body propagation from ellipse
           CALL vsumg(6,eq0,elm(1,n),elm(1,n))
           njc=0
           iplam=0
           CALL proele('EQU',t0,elm(1,n),tafter,eq1)
           WRITE(iun9,109)eq1
 109       FORMAT(6f20.15)
c check for existence of data
c PROBLEM: jcsel is always the same???
           IF(njc.eq.0.or.iplam.ne.iplamold)THEN
              WRITE(*,*)' no close approach to planet ',iplamold
           ELSEIF(njc.lt.jcsel)THEN
              WRITE(*,*)n,' multiple close approach',njc,jcsel
           ELSEIF(tcla(jcsel).gt.max(tbefore,tafter).or.
     +          tcla(jcsel).lt.min(tbefore,tafter))THEN
              WRITE(*,130)n,tcla(jcsel), tbefore, tafter
 130          FORMAT(i4,' time_cla=',f11.2,' out of range, must be in '
     +             ,2f11.2)
           ELSE
              nc=nc+1
c reference system with third axis normal to ecliptic, first=velocity
c is available from tpcana.h
c             vc0=vsize(vcla(1,jcsel))
c             DO i=1,3
c                v1(i)=vcla(i,jcsel)/vc0
c             ENDDO
c             v2(1)=0.d0
c             v2(2)=0.d0
c             v2(3)=1.d0
c             CALL mtpref(v1,v2,vv3,vvt3) 
              CALL mulmav (vt3,3,3,xcla(1,jcsel),3,xc)
              CALL mulmav (vt3,3,3,vcla(1,jcsel),3,vc)
              xcl(nc)=xc(3)
              ycl(nc)=xc(2)
c file output: close approach manifold
              WRITE(*,131)xc(3),xc(2),tcla(jcsel),n,nc
 131          FORMAT(f7.4,1x,f7.4,1x,f11.4,1x,i4,1x,i4)
              WRITE(iun8,107)xc,vc,tcla(jcsel)
           ENDIF
        ENDIF
 7      ENDDO
      IF(inl.eq.1)THEN
         CALL filclo(iun7,' ')
      ELSEIF(inl.eq.2)THEN
         CALL filclo(iun8,' ')
         CALL filclo(iun9,' ')
      ENDIF
      iplam=iplamold
c ======================================================================
c graphics: 
c close approach manifold number of points
      IF(ibv.eq.1.and.nc.gt.0)THEN
         npoc=nc+1
         xcl(npoc)=xcl(1)
         ycl(npoc)=ycl(1)
      ELSE
         npoc=nc
      ENDIF
      IF(npoc.eq.0)THEN
         WRITE(*,*)' nothing to plot!'
         RETURN
      ENDIF    
c labels
      IF(inl.eq.1)THEN
         xlabel=' target plane xi (AU), linear approximation '
      ELSE
         xlabel=' target plane xi (AU), semilinear approximation '
      ENDIF
      ylabel=' target plane zeta (AU) (proj, ecl. normal)'
c selection of graphics device
      WRITE(title,200)astna0,tcla0
 200  FORMAT(a18,' encounter at MJD ',f9.2)
 2    CALL getdev(idev)
      IF(idev.eq.0)RETURN
c If we are making .ps files then the MTP will be overwritten!
      if(idev.eq.5.or.idev.eq.6)then
         write(*,*)'The PostScript file ''giffv.ps'' is '//
     +        'about to be overwritten. If you wish to save it'//
     +        ' you should rename it before proceeding.'
         pause
      endif
c plot, mark Earth
      CALL plotob(xcl,ycl,0.d0,0.d0,npoc,xlabel,ylabel,title,idev,1)
      GOTO 2
      END
c ==============================================
c  VIRtual IMPactor ephemerides
      SUBROUTINE virimp(tpc,dtpdet,axes,sig,ceicel,b,v,tc,eqc,gc,cc,
     +     hmag,gmag,iun20)
      IMPLICIT NONE
c ===============input=================
c magnitude
      DOUBLE PRECISION hmag,gmag
c output unit
      INTEGER iun20
c target plane coordinates: cartesian
      DOUBLE PRECISION tpc(2)
c computation of target ellipse
      DOUBLE PRECISION axes(2,2),sig(2),dtpdet(6,2)
c correspondent ellipse in orbital lements space
      DOUBLE PRECISION b(2,2),ceicel(4,2),v(6,6)
c initial conditions: epoch, elements
      DOUBLE PRECISION tc,eqc(6)
c normal and covariance matrices
      DOUBLE PRECISION gc(6,6),cc(6,6)
c ==============end interface=================
c option flag
      INTEGER ivir
c sigma level across the LOV
      DOUBLE PRECISION sigimp
c safe distance from the Earth
      DOUBLE PRECISION dsafe
c rectangle, and its copy in the elements space
      DOUBLE PRECISION tpstr(2),tpwea(2)
      DOUBLE PRECISION elems(6,4),del2(2),del4(4),dels(6)
      DOUBLE PRECISION delems(6),delemw(6),selems(6),selemw(6)
c observation data (proper motion, elongation, distance)
      INCLUDE 'phase.h'
c times
      DOUBLE PRECISION tmjd,tut,sec1,sec2
      INTEGER mjd1,mjd2
c call to seleph
      DOUBLE PRECISION  tut1,tdt1,tut2,tdt2,dt     
      INTEGER idsta
      CHARACTER*3 scale
c for preobs
      DOUBLE PRECISION alpha,delta,appmag
c station code, observation type
      INTEGER ids,iob1
c integer indexes, lengths, functions, units
      INTEGER j,jj,intlo,ln,iuneph,iunvir
c for outobc: covariance of the observations
      DOUBLE PRECISION gamad(2,2),axesky(2,2),sigsky(2)
c ephemeris options
      CHARACTER*80 fields
      DOUBLE PRECISION mass
      CHARACTER*45 file
      INCLUDE 'trig.h'
c file name
      CHARACTER*80 titnam
c      DOUBLE PRECISION av(5),dv(5)
c     DOUBLE PRECISION dxl,dyl,prscag
c tests
c     DOUBLE PRECISION aa(2,6),id(2,2),dtpcde(2,6)
c acll to preob4
c confidence boundary, line of max variation (alpha, delta, app. magnitude)
      INTEGER npo, npox, ibv, npo1, inl
      PARAMETER (npox=4000)
      DOUBLE PRECISION sigma, al(npox),de(npox),appmagv(npox)
c line of elements
      DOUBLE PRECISION elm(6,npox)
c menu empty items
      CHARACTER*20 menunam
      CHARACTER*50 s4,s5,s6,s7,s8,s9,s10
c multiple data for confidence boundary
      INCLUDE 'npoint.h'
c =====================================================
c chase is open for the virtual impactor; what to do?
 10   menunam='null'
      CALL menu(ivir,menunam,3,'how to catch it?=',
     +     'exploratory ephemerides=',
     +     'select observation time=',
     +     'output impact orbital elements=',
     +     s4,s5,s6,s7,s8,s9,s10)
      IF(ivir.eq.0)RETURN
      IF(ivir.eq.1)THEN
c ======= GENERATE EPHEMERIS =========
c select time interval, step
         CALL seleph(tut1,tdt1,tut2,tdt2,dt,idsta)
         CALL filnam('.','virimp','eph',file,ln)
         CALL filopn(iuneph,file(1:ln),'unknown')
         fields='cal,mjd,coord,mag,elong,glat,r,delta,appmot,skyerr'
         scale='UTC'
         IF(nint(abs(tdt2-tdt1)/dt).gt.500)THEN
            write(*,*)'Too many ephemeris points: ',
     +           nint(abs(tdt2-tdt1)/dt)
            write(*,*)'Select a time interval and time span to ',
     +        'ensure that there are fewer than 500 points.'
            goto 10
         ELSE
            CALL ephemc(iuneph,'EQU',tc,eqc,gc,.true.,tdt1,tdt2,
     +           dt,mass,hmag,gmag,idsta,scale,fields)
         ENDIF
         CALL filclo(iuneph,' ')
         WRITE(*,*)' look at the ephemerides in file ./virimp.eph '
      ELSEIF(ivir.eq.3)THEN
         CALL filopn(iunvir,'virimp.eq0','unknown')
c output header 
         CALL wromlh (iunvir,'ECLM','J2000')
         CALL wromlr (iunvir,'virimp',eqc,'EQU',tc,gc,.true.,
     +        cc,.true.,hmag,gmag,0.d0)
         CALL filclo(iunvir,' ')
      ELSEIF(ivir.eq.2)THEN
c ============PREDICT OBSERVATION==================
c ======two dimensional preimage===============
c compute rectangle on the MTP  enclosing all the impact points;
c note that the target plane point tpc is used as origin
         sigimp=5.d0
         dsafe=5.d0*4.2e-5
         DO j=1,2
            tpstr(j)=axes(j,1)*sig(1)*sigimp
            tpwea(j)=axes(j,2)*dsafe
         ENDDO
c find corresponding orbital elements at epoch
         CALL mulmav(b,2,2,tpstr,2,del2)
         CALL mulmav(ceicel,4,2,del2,2,del4)
         CALL vcopy(2,del2,dels)
         DO jj=1,4
            dels(2+jj)=-del4(jj)
         ENDDO
         CALL mulmav(v,6,6,dels,6,delems)
c        WRITE(*,*)tpstr
c        WRITE(*,*)delems
c linear map from ellipse
c        dxl=prscag(6,delems,dtpdet(1,1))
c        dyl=prscag(6,delems,dtpdet(1,2))
c        WRITE(*,*)dxl,dyl
c find corresponding orbital elements at epoch
         CALL mulmav(b,2,2,tpwea,2,del2)
         CALL mulmav(ceicel,4,2,del2,2,del4)
         CALL vcopy(2,del2,dels)
         DO jj=1,4
            dels(2+jj)=-del4(jj)
         ENDDO
         CALL mulmav(v,6,6,dels,6,delemw)
c        WRITE(*,*)tpwea
c        WRITE(*,*)delemw
c linear map from ellipse
c        dxl=prscag(6,delemw,dtpdet(1,1))
c        dyl=prscag(6,delemw,dtpdet(1,2))
c        WRITE(*,*)dxl,dyl
c check 
c        CALL transp(dtpdet,6,2,dtpcde)
c        CALL mulmat(dtpcde,2,6,v,6,6,aa)
c        WRITE(*,*)' aa=',aa
c        CALL mulmat(aa,2,2,b,2,2,id)
c        WRITE(*,*)'id=',id
c compute 4 corners
         CALL vsumg(6,eqc,delems,elems(1,1))
         CALL vsumg(6,elems(1,1),delemw,elems(1,2))
         DO jj=1,6
            selemw(jj)=-delemw(jj)
            selems(jj)=-delems(jj)
         ENDDO
         CALL vsumg(6,elems(1,1),selemw,elems(1,1))
         CALL vsumg(6,eqc,selems,elems(1,3))
         CALL vsumg(6,elems(1,3),delemw,elems(1,4))
         CALL vsumg(6,elems(1,3),selemw,elems(1,3))
c give time
         WRITE(*,*)' give time for prediction (MJD)'
         READ(*,*)tmjd
c universal time of the required observation 
         mjd1=intlo(tmjd)
         sec1=(tmjd-float(mjd1))*86400.d0
         CALL cnvtim(mjd1,sec1,'TDT',mjd2,sec2,'UTC')
         tut=sec2/86400.d0+float(mjd2)
         ids=500
         iob1=1001
c predict
         CALL preobc('EQU',tc,ids,tmjd,eqc,hmag,gmag,gc,
     +        iob1,alpha,delta,appmag,gamad,sigsky,axesky)
         CALL outobc(iun20,iob1,ids,tut,alpha,delta,appmag,adot,ddot,
     +        elo,dis,2,gamad,sigsky,axesky)
         DO jj=1,4
           CALL preobs('EQU',tc,ids,tmjd,elems(1,jj),iob1,al(jj),de(jj),
     +           hmag,gmag,appmag)
            al(jj)=al(jj)-alpha
            de(jj)=de(jj)-delta
c           write(*,*)jj,al(jj)*degrad,de(jj)*degrad
         ENDDO
c ======four dimensional preimage===============
         menunam='prednonl'
         CALL menu(inl,menunam,3,'How to handle nonlinearity?=',
     +        'linear map=',
     +        '2-body nonlinearity=',
     +        'full n-body nonlinearity=',
     +        s4,s5,s6,s7,s8,s9,s10)
         IF(inl.eq.0)RETURN
c input specification of set of points
         CALL asscbd(iun20,npox,npo,sigma,ibv)   
c        sigma=5.d0
c        ibv=0
c        inl=1
c        npo=20
         CALL preob4(tc,ids,tmjd,eqc,hmag,gmag,gc,
     +        cc,v,sigma,npo,ibv,inl,al(5),de(5),appmagv,elm,
     +        alpha,delta,appmag,gamad,sigsky,axesky,npo1)
         al(5+npo1)=al(5)
         de(5+npo1)=de(5)
c        DO j=1,npo1
c           DO jj=1,6
c              delemw(jj)=elm(jj,j)-eqc(jj)
c           ENDDO
c linear map from 4-d axis
c           dxl=prscag(6,delemw,dtpdet(1,1))
c           dyl=prscag(6,delemw,dtpdet(1,2))
c           WRITE(*,*)' alternate no. ',j,dxl,dyl,delemw
c        ENDDO
         DO j=1,npo1+5
            write(*,*)'diff. observation ', j,al(j)*degrad,de(j)*degrad
         ENDDO
         titnam='virtual impactor' 
         CALL plocbd(titnam,alpha,delta,5.d0,tut,al,de,5+npo1,iob1)
      ENDIF
      GOTO 10
      END
c ================================================================
c AFTCLOV
c select time interval to get after the close approach
c taking into account the relative velocity w.r. to Earth
      SUBROUTINE aftclov(iplam,t0,tcla,v_inf,tbefore,tafter)
      IMPLICIT NONE
c input: planet number, time of initial conditions, 
c of known close approach, velocity
      INTEGER iplam
      DOUBLE PRECISION t0,tcla,v_inf
c output: time "after", time "before" (depending upon sense of propagation
      DOUBLE PRECISION tafter,tbefore
c end interface
c time interval to exit from TP disk
      DOUBLE PRECISION delt_tp
c target plane disk radius
      INCLUDE 'parbep.h'
      INCLUDE 'masses.h'
c ===================================================================
c warning: really done only for Earth
      IF(ordnam(iplam).ne.'EARTH') THEN
         WRITE(*,*)' aftclov: not to be used for planet ',ordnam(iplam)
      ENDIF
c time interval to exit from TP disk
      delt_tp=2*dmin(iplam)/v_inf
c forced to avoid short intervals for fast encounters
      IF(delt_tp.lt.50.d0)delt_tp=50.d0
c time interval to be clear out of the TP disk is given as deltat
      IF(t0.lt.tcla)THEN
c future close approaches
         tafter=tcla+delt_tp
         tbefore=tcla-delt_tp
      ELSE
c past close approaches
         tafter=tcla-delt_tp
         tbefore=tcla+delt_tp
      ENDIF
      RETURN
      END
c =====================================================================
c PREOB4- virtual impactor negative observation skyhole
c ============INTERFACE===================================================
      SUBROUTINE preob4(t0,idsta,t1,eq,h,g,gameq,
     +    ceq,v6,sigma,npo,ibv,inl,al,de,hmagv,elm,
     +    alpha,delta,hmagn,gamad,sig,axes,npo1)
      IMPLICIT NONE      
c ============= input ====================================================
c elements and epoch times, covariance and normal matrices at t0,
c sigmas for the boundary
      DOUBLE PRECISION eq(6),t0,t1,gameq(6,6),ceq(6,6),sigma
c number of points, flag for confidence bd/line of variation, nonlinearity
      INTEGER npo,ibv,inl
c magnitude
      DOUBLE PRECISION h,g,hmagn
c station code
      INTEGER idsta
c 6x6 rotation matrix giving the reference system adapted to the target plane
      DOUBLE PRECISION v6(6,6)
c ============= output ===================================================
c points on the confidence boundary (difference w.r. to alpha,delta)
c WARNING! the output number of points is npo1.le.npo; 
c this beacuse hyperbolic points are discarded
      INCLUDE 'npoint.h'
      INTEGER npo1
      DOUBLE PRECISION al(npo),de(npo),hmagv(npo)
c line of elements
      DOUBLE PRECISION elm(6,npo)
c best fit observations
      DOUBLE PRECISION alpha,delta
c covariance
      DOUBLE PRECISION gamad(2,2),axes(2,2),sig(2)
c ============END INTERFACE===============================================
c workspace
      DOUBLE PRECISION tmp(6,6),daddelt(6,2),v6t(6,6)
c partials in new reference system
      DOUBLE PRECISION daddelt4(4,2),g4(4,4),gamv(6,6),c4(4,4),cv(6,6)
      DOUBLE PRECISION ws(4)
      INTEGER ii,jj,indp
c partial derivatives of alpha, delta, w.r. to elements (by columns)
      DOUBLE PRECISION daddet(6,2),dummy(6)
c second derivatives of alpha, delta, w.r. to elements (not used)
      DOUBLE PRECISION ddade(6,6),dddde(6,6)
c ===================================================================
c orthonormal basis, matrix defining the plane of the ellipse
      DOUBLE PRECISION v(4,4),ceicel(2,2)
c transformation matrix between the two planes
      DOUBLE PRECISION b(2,2)
c number of full revolutions around the sky
      INTEGER ng,nrev
c functions
      DOUBLE PRECISION appmag,prscag
c elongation,distance to Earth, distance to Sun (to compute magnitude)
      INCLUDE 'phase.h'
      DOUBLE PRECISION adot0,ddot0 
c ===================================================================
c constant of gravitation, trigonometric constants 
      INCLUDE 'sunmass.h'
      INCLUDE 'trig.h'
c temporaries, indexes
      DOUBLE PRECISION dal,ddl,maxsig,minsig
      INTEGER n
c flag for 2-body approximation; must be .false. for full n-body computation
      LOGICAL twobo
      twobo=.false.
c****************
c   static memory not required
c****************
c =====================================================================
c compute observation; derivatives (of order 1) required            
      CALL alfdel (eq,t0,t1,idsta,alpha,delta,daddet(1,1),daddet(1,2),
     +        1,twobo,ddade,dddde)
c store proper motion
      adot0=adot
      ddot0=ddot
c compute derivatives in the reference system adapted to the target plane
      CALL transp(v6,6,6,v6t)
      CALL mulmat(v6t,6,6,daddet,6,2,daddelt)
c     WRITE(1,*)'daddet ', daddet
c     WRITE(1,*)'daddelt ',daddelt
c     WRITE(1,*)'v6 ',v6
c compute covariance matrix in the new reference system
      CALL mulmat(v6t,6,6,gameq,6,6,tmp)
      CALL mulmat(tmp,6,6,v6,6,6,gamv)
      CALL mulmat(v6t,6,6,ceq,6,6,tmp)
      CALL mulmat(tmp,6,6,v6,6,6,cv)
c reduce all to 4-d
      DO jj=1,4
        DO ii=1,4
c          g4(jj,ii)=gamv(jj+2,ii+2)
          c4(jj,ii)=cv(jj+2,ii+2)
        ENDDO
        DO ii=1,2
          daddelt4(jj,ii)=daddelt(jj+2,ii)
        ENDDO
      ENDDO
      CALL tchinv(c4,4,g4,ws,indp)
c     WRITE(1,*)'daddelt4 ',daddelt4
c     WRITE(1,*)'c4 ',c4
c     WRITE(1,*)'g4 ',g4
c compute ellipse of covariance of alpha,delta
      CALL ellip4(daddelt4,g4,sig,axes,gamad)
c     write(1,*)' sig,axes,gamad ',sig,axes,gamad
c use nonlinear method
      IF(ibv.eq.0)THEN
         maxsig=max(sig(1),sig(2))
         minsig=min(sig(1),sig(2))
         if(maxsig/minsig.le.200.d0)then
            ibv=1
         else
            ibv=2
         endif
      endif    
c =====================================================================
c compute ellipse in the elements space 
      CALL slinel4(daddelt4,g4,c4,ceicel,b,v)
c     WRITE(1,*)'ceicel,b,v4 ',ceicel,b,v
c ===========================================================
c compute line of orbital elements
      CALL linobs4(ibv,npo,eq,axes,sig,b,v,sigma,ceicel,elm,v6,npo1)
c     WRITE(1,*)' npo1',npo1
c     DO jj=1,npo1
c        WRITE(1,*)'elm,jj=',jj,(elm(ii,jj),ii=1,6)
c     ENDDO
c ===========================================================
      ng=0
      DO 7 n=1,npo1
c chose method to handle nonlinearity
        IF(inl.eq.1)THEN
c linear map from ellipse
           dal=prscag(6,elm(1,n),daddet(1,1))
           ddl=prscag(6,elm(1,n),daddet(1,2))
           al(n)=dal
           de(n)=ddl
           CALL vsumg(6,eq,elm(1,n),elm(1,n))
c apparent magnitude is the one of the nominal orbit
           hmagv(n)=hmagn
        ELSEIF(inl.eq.2)THEN
           write(*,*)' inl=2 is forbidden'
           npo1=0
           RETURN
        ELSEIF(inl.eq.3)THEN
c full n-body propagation from ellipse 
           CALL vsumg(6,eq,elm(1,n),elm(1,n))
           CALL alfdel (elm(1,n),t0,t1,idsta,al(n),de(n),
     +          dummy,dummy,0,twobo,ddade,dddde)
           al(n)=al(n)-alpha
           de(n)=de(n)-delta
c          write(*,*)n,al(n),de(n)
c other prediction data stored in common
           phav(n)=pha
           disv(n)=dis
           dsunv(n)=dsun
           elov(n)=elo
           gallav(n)=gallat
           adotv(n)=adot
           ddotv(n)=ddot
c compute apparent magnitude at time of observation
           hmagv(n)=appmag(h,g,dsun,dis,pha)
        ELSE
           WRITE(*,*)' preobn: this we have not invented yet ', inl
           RETURN           
        ENDIF
c keep count of lost revolutions
        IF(n.eq.1)THEN
           IF(al(n).gt.pig)al(n)=al(n)-dpig
        ELSE
           CALL angupd(al(n),al(n-1),ng)
        ENDIF
c temporary output
c       write(*,*)'Solution ',n,', RA/DEC (deg)',
c    +       al(n)*degrad,de(n)*degrad,ng
 7    continue
c =====================================================================
c ensure that LOV is consistent with nominal point
c first find midpoint of LOV, assume npo is even
      if(ibv.eq.2)then
         nrev=nint((al(npo/2)+al(npo/2+1))/2.d0/dpig)
c        write(*,*)'debug: nrev:',nrev
         if(nrev.ne.0)then
            do n=1,npo1
               al(n)=al(n)-nrev*dpig
            enddo
         endif
      endif
c restore original apparent motion
      adot=adot0
      ddot=ddot0
      RETURN
      END
c ===========================================================
c common subroutines for preobn and fclan
c patch 1.6.1, A. Milani, May 2, 1998
c ===========================================================
c LINOBS defines line of changes in orbital elements to be used for
c confidence boundary/variations line 
c ===========================================================
      SUBROUTINE linobs4(ibv,npo,eq,axes,sig,b,v,sigma,ceicel,elm,
     +  v6,npo1)
      IMPLICIT NONE
c ====================INPUT==================================
      INTEGER ibv,npo
      DOUBLE PRECISION eq(6)
      DOUBLE PRECISION axes(2,2),sig(2),b(2,2),sigma
c matrix defining the plane of the ellipse,new orthonormal reference
      DOUBLE PRECISION ceicel(2,2),v(4,4)
c matrix defining the target plane adapted referemce system
      DOUBLE PRECISION v6(6,6)
c ===================OUTPUT==================================
      INTEGER npo1
      DOUBLE PRECISION elm(6,npo)
c ==================END INTERFACE============================
      INTEGER nn,n,i
      DOUBLE PRECISION elm4(4)
      DOUBLE PRECISION s,x,y,vad(2),xv,yv,dn,dth,theta,xa,yd
      DOUBLE PRECISION eqnew(6)
      DOUBLE PRECISION alde(2),ecc
      INCLUDE 'trig.h'
c =====================================================================
c line of maximum variation: in the alpha-delta plane
      DO i=1,2
        vad(i)=axes(i,2)*sig(2)
      ENDDO
c in the elements space
      xv=(b(1,1)*vad(1)+b(1,2)*vad(2))
      yv=(b(2,1)*vad(1)+b(2,2)*vad(2))
c     WRITE(*,*)xv,yv
c direction not used any more 
c     theta0=atan2(yv,xv)
c linear step for variation axis parametrisation
      dn=2.d0/float(npo-1)
c angular step for ellipse parametrisation
      dth=dpig/float(npo)
c ===========================================================
c main loop on the number of output points
      nn=0
      DO 7 n=1,npo
c ===========================================================
c choice between two output options
        IF(ibv.eq.2)THEN
c ===========================================================
c line of maximum variation in the elements space
           s=(n-1)*dn-1.d0
           x=sigma*s*xv
           y=sigma*s*yv
        ELSEIF(ibv.eq.1)THEN
c =====================================================================
c parametrisation of the ellipse in the subspace of elements, based upon the
c parametrisation of the ellipse in the alpha-delta plane
c WARNING: npo must be divisible by 2, otherwise one tip of the
c banana would be missed
           theta=(n-1)*dth
           xa=sig(1)*cos(theta)*sigma
           yd=sig(2)*sin(theta)*sigma
           CALL lincog(2,axes(1,1),xa,axes(1,2),yd,alde)
c transfer of parametrisation in the V1,V2 plane
           x=(b(1,1)*alde(1)+b(1,2)*alde(2))
           y=(b(2,1)*alde(1)+b(2,2)*alde(2))
        ELSE
           write(*,*)' linobs: this should not happen,ibv=',ibv
        ENDIF
c compute displacement on the confidence ellipsoid corresponding to x,y
        nn=nn+1    
        CALL elemo4(x,y,v,ceicel,elm4)
        CALL mulmav(v6(1,3),6,4,elm4,4,elm(1,nn))
c       write(*,*)(elm(i,nn),i=1,6)
c add to the original center of the ellipsoid of confidence
        CALL vsumg(6,eq,elm(1,nn),eqnew)
        ecc=sqrt(eqnew(2)**2+eqnew(3)**2)
        IF(ecc.ge.1.d0.or.eqnew(1).le.0.d0)THEN
           write(*,*)' Hyperbolic, ecc=',ecc,' a=',eqnew(1)
           nn=nn-1
        ELSEIF(ecc.ge.0.99d0)THEN
           write(*,*)' Almost Hyperbolic, ecc=',ecc,' a=',eqnew(1)
           nn=nn-1
        ENDIF
 7    continue
c final count of non hyperbolic orbits
      npo1=nn  
      RETURN
      END
c =====================================================================
c ELLIPS
c compute covariance ellipse of two observables
c =====================================================================
      SUBROUTINE ellip4(daddet,gamm0,sig,axes,gamad)
      IMPLICIT NONE
c input covariance matrix
      DOUBLE PRECISION gamm0(4,4)
c input partial derivatives of alpha, delta, w.r. to elements (by columns)
      DOUBLE PRECISION daddet(4,2)
c output covariance
      DOUBLE PRECISION gamad(2,2),axes(2,2),sig(2)
c ==============END INTERFACE==========================================
c eigenvalues, workspace, transposed
      DOUBLE PRECISION eigval(2),tmp24(2,4),fv1(2),fv2(2),dadde(2,4)
c loop indexes
      INTEGER i
c error flag
      INTEGER ierr
c =====================================================================
      CALL transp(daddet,4,2,dadde)
      CALL mulmat(dadde,2,4,gamm0,4,4,tmp24)
      CALL mulmat(tmp24,2,4,daddet,4,2,gamad)
c =====================================================================
c compute ellipse of confidence
c eigenvalues
      CALL rs(2,2,gamad,eigval,1,axes,fv1,fv2,ierr)
      DO  i=1,2
        IF(eigval(i).gt.0.d0)THEN
           sig(i)=sqrt(eigval(i))
        ELSE
           write(*,*) 'non positive eigenvalue'
           sig(i)=0.d0
        ENDIF
      ENDDO
      RETURN
      END
c =====================================================================
c ELEMOV
c compute displacement on the confidence ellipsoid corresponding to x,y
c on the plane of the gradients of alpha-delta
c =====================================================================
      SUBROUTINE elemo4(x,y,v,ceicel,del)
      IMPLICIT NONE
c inout/output
      DOUBLE PRECISION x,y,v(4,4),del(4)
      DOUBLE PRECISION ceicel(4,2)
c workspace
      DOUBLE PRECISION dee(2),deel(4)
c ===================
      CALL lincog(4,v(1,1),x,v(1,2),y,del)
      CALL lincog(2,ceicel(1,1),-x,ceicel(1,2),-y,dee)
      CALL mulmav(v(1,3),4,2,dee,2,deel)
      CALL vsumg(4,del,deel,del)
      RETURN
      END
c ===========================================================
c SLINEL
c semilinear boundary ellipse computation
c =========================================================== 
      SUBROUTINE slinel4(dtpdet,gc,cc,ceicel,b,v)
      IMPLICIT NONE
      INTEGER ndim,ndimm2
      PARAMETER (ndim=4,ndimm2=ndim-2)
c 6 by 2 matrix with columns= gradients
      DOUBLE PRECISION dtpdet(ndim,ndimm2)
c normal and covariance matrices
      DOUBLE PRECISION gc(ndim,ndim),cc(ndim,ndim)
c orthonormal basis
      DOUBLE PRECISION v(ndim,ndim),vt(ndim,ndim)
c       ,gamv(ndim,ndim)
      DOUBLE PRECISION cv(ndim,ndim),tmp(ndim,ndim)
c partial matrices
      DOUBLE PRECISION c4(ndimm2,ndimm2),cinv(ndimm2,ndimm2)
      DOUBLE PRECISION c42(ndimm2,2),ceicel(ndimm2,2)
c line of maximum variation
      DOUBLE PRECISION a(2,2),b(2,2),deta
c loop indexes ii=1,2, ij,ijj=1,ndim
      INTEGER ii, ij, ijj
c for inversion with tcholevski: workspace, error flag
      DOUBLE PRECISION ws(ndimm2)
      INTEGER ierr
      DOUBLE PRECISION prscag
c =====================================================================
c adapted orthonormal basis, covariance and normal matrix in the new basis
      CALL graha(dtpdet,ndim,v)
      CALL transp(v,ndim,ndim,vt)
c     CALL mulmat(vt,ndim,ndim,gc,ndim,ndim,tmp)
c     CALL mulmat(tmp,ndim,ndim,v,ndim,ndim,gamv)
      CALL mulmat(vt,ndim,ndim,cc,ndim,ndim,tmp)
      CALL mulmat(tmp,ndim,ndim,v,ndim,ndim,cv)
c =====================================================================
c 4x4 and 4x2 submatrices of normal matrix
      do 15 ijj=1,ndimm2
        DO ij=1,ndimm2
          c4(ijj,ij)=cv(ijj+2,ij+2)
        ENDDO
        DO  ii=1,2
          c42(ijj,ii)=cv(ijj+2,ii)
        ENDDO
 15   continue
c ===========================================================
c Cholewski method for inversion
      CALL tchinv(c4,ndimm2,cinv,ws,ierr)
      IF(ierr.ne.0)THEN
         write(*,*)' decide what to do, ierr=',ierr
      ENDIF
c ===========================================================
c matrix to be used for out of plane component
      CALL mulmat(cinv,ndimm2,ndimm2,c42,ndimm2,2,ceicel)
c ===========================================================
c linear map from the elements space (with base V) and the alpha-delta plane 
      a(1,1)=prscag(ndim,dtpdet(1,1),v(1,1))
      a(1,2)=prscag(ndim,dtpdet(1,1),v(1,2))
      a(2,1)=prscag(ndim,dtpdet(1,2),v(1,1))
      a(2,2)=prscag(ndim,dtpdet(1,2),v(1,2))
      CALL inv22(a,b,deta)
      RETURN
      END









