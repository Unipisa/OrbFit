c ===========MODULE close_app=============
c PUBLIC ROUTINES
c            cloapp
c            clocms
c            strclan3
c            cov_avai/cov_not_av
c            clotest
c
c MODULE CONTAINS
c ROUTINES
c               stepcon
c               falsi
c               [strmtp]
c               strclo  
c               rkg
c               fctcl
c               wriclan3
c               reaclan3
c               reaclorec
c               mtpref
c               mtprot3
c               marg2
c
c HEADERS     
c close_app.o: proout.h comdif.h 
c
c             covariance.h private to close_app and target_plane
c             closapl.h with target_plane.f
c             cloapp.h        "
c             tpcana.h        "
c             dx0de.h  with propag_state.f
c             nvarx.h        "
c             parint.h       "
c             comint.h       "
c             rkcoef.h       "
c             closta.h       "
c             sunmass.h trig.h public
c             parbep.h masses.h only number, names, disk radius 
c             iclap.h  only iclap
c
c Note: close_app USES target_plane
c 
c ==============================================================
c CLOAPP 
c close approach control; driver for falsi
c vers. 1.9.2, A.Milani, August 25, 1999
c with decreasing stepsize near closest approach
c also with end of close approach flag
c ==============================================================
      SUBROUTINE cloapp(tcur,th,xa,va,nv,idc,xpla,xldir,dir,nes,cloend)
      IMPLICIT NONE
c close approach control: which planet is close (0=none)
      INTEGER idc
c stepsize control: flag for fixed step, stepsize, direction
      LOGICAL nes
      DOUBLE PRECISION xldir,dir
c current time, previous stepsize
      DOUBLE PRECISION tcur,th
c positon, velocity of asteroid, of close encounter planet
      INTEGER nv
      DOUBLE PRECISION xa(nv),va(nv)
      DOUBLE PRECISION xpla(6)
c logical flag for close approach termination
c to allow for storage of state transition matrix
      LOGICAL cloend
c =========END INTERFACE========================================
c stepsize tricks
      LOGICAL nesold
      DOUBLE PRECISION xlold,tmin,dt,thfirst
c               ,tleft
c which planet is being approached, vector difference
      INTEGER ic
      DOUBLE PRECISION x(3),v(3)
c
      DOUBLE PRECISION vsize
      INCLUDE 'parbep.h'
      INCLUDE 'masses.h'
      INCLUDE 'proout.h'
c common data on close appr. 
      INCLUDE 'closapl.h'
      INCLUDE 'cloapp.h'
c counters for arrays of multiple minima, sections
      INTEGER jc,jt
c loop indexes i=1,3
      INTEGER i
c logical flags for regulae falsi being initiaited
      LOGICAL first
c startup from a non-close approaching state at each integration restart
      INCLUDE 'closta.h'
c static memory model required
      SAVE 
c ===========================================
      IF(clost)THEN
         ic=0
         clost=.false.
      ENDIF
c close approach is not ending, in most cases
      cloend=.false.
c ========================================================
c   control on close approaches
      IF(ic.eq.0)THEN
c close approach not going on; check if it is beginning
         IF(idc.ne.0)THEN
c close approach detected (first step inside influence region)
c relative position and velocity
            DO  i=1,3
              x(i)=xa(i)-xpla(i)
              v(i)=va(i)-xpla(i+3)
            ENDDO
c planet with which close approach is taking place
            ic=idc
            iplam=idc
c           WRITE(iuncla,*)' approach to planet', iplam
c close approach minima counter
            jc=0
c target plane crossing counter
            jt=0
c first call to falsi
            first=.true.
            CALL falsi(tcur,xa,va,nv,xpla,jc,jt,first,iplam)
c setup of fixed stepsize
            nesold=nes
            nes=.true.
            xlold=xldir
c  stepsize based upon angular velocity
            CALL stepcon(x,v,npoint,dmin(idc),dir,dt)
c not allowed to be longer than the one used before the close approach
            xldir=min(th,dt)
c which is stored for alter use
            thfirst=th
         ENDIF
      ELSE
c close approach taking place
c    control of inconsistences
         IF(idc.ne.0)THEN
            IF(idc.ne.ic.or.idc.gt.nmass.or.idc.lt.0)THEN
               WRITE(*,*)' closap: this should not happen',idc,ic
            ENDIF
c relative position and velocity
            DO  i=1,3
              x(i)=xa(i)-xpla(i)
              v(i)=va(i)-xpla(i+3)
            ENDDO
            first=.false.
            CALL falsi(tcur,xa,va,nv,xpla,jc,jt,first,iplam)
c  stepsize based upon angular velocity
            CALL stepcon(x,v,npoint,dmin(idc),dir,dt)
c not allowed to be longer than the one used before the close approach
            xldir=min(thfirst,dt)
         ELSE
c   close approach ended; reset flags etc.
            cloend=.true.
            xldir=xlold
            nes=nesold
            ic=0
            njc=jc
            jc=0
         ENDIF
      ENDIF
c =================================================
      RETURN
      END
c===============================
c CLOCMS
c close approach control
      SUBROUTINE clocms(idc,tt,xxpla)
      IMPLICIT NONE
      INTEGER idc
      DOUBLE PRECISION tt,xxpla(6)
      include 'proout.h'
      include 'parbep.h'
      include 'masses.h'
      include 'iclap.h'
      IF(iorb.eq.11)THEN
         if(iclap.ne.0.and.idc.ne.0)then
* to be improved with a real safety feature
            write(*,*)'t =',tt,' close approach to=',ordnam(idc)
            write(iuncla,*)'t =',tt,' close approach to =',ordnam(idc)
         endif 
      ELSEIF(iorb.eq.9)THEN
         if(idc.ne.0)then
            write(*,*)'t =',tt,' close approach code=',idc
            write(iuncla,*)'t =',tt,' close approach code =',idc
         endif
      ENDIF
      RETURN
      END
c ==================================================
c STRCLAN3 3-d version 30 Jan 2002
c compute close approach MTP analysys, and store in .clolin record
c ==================================================
      SUBROUTINE strclan3(dx1dx0)
      IMPLICIT NONE
c ============INPUT========================
c accumulated state transition matrix
      DOUBLE PRECISION dx1dx0(6,6)
c ========HIDDEN INPUT======================
c target plane coordinates and derivatives, moid
      INCLUDE 'cloapp.h'
c derivatives with respect to elements
      INCLUDE 'dx0de.h'
c covariance matrix
      INCLUDE 'covariance.h'
c planetary data (names)
      INCLUDE 'parbep.h'
      INCLUDE 'masses.h'
c output unit
      INCLUDE 'proout.h'
c =========HIDDEN OUTPUT====================
c in the .clo file
c csi, zeta, s,w,alpha
      DOUBLE PRECISION csi, zeta, stretch, width, alpha
c via common
c basis adapted to MTP 
c output from mtprot
      INCLUDE 'tpcana.h'
c =========END INTERFACE===============
c partial derivatives
      INCLUDE 'nvarx.h'
      DOUBLE PRECISION dxdx1(6,6),dxde(6,6),dxdx0(6,6),y2(nvarx)
      INTEGER nv
c basis adapted to MTP 
      INTEGER jc
      DOUBLE PRECISION v2(3)
c minimum distance
      DOUBLE PRECISION r,rdot,vsize,prscal
c angle alpha
      DOUBLE PRECISION cosa,sina
c weak direction
      DOUBLE PRECISION wdir(6),sdir,prscag
c planet names
      CHARACTER*30 planam
      INTEGER lpla,lench
c calendar date variables
      INTEGER iyear,imonth,iday
      DOUBLE PRECISION hour
      CHARACTER*16 date
c mtp normal
      DOUBLE PRECISION tpno(3),vc0
c loop indexes
      INTEGER i
c =====================================================================
c partials are always available if we are here
      nv =21
c derivative are required: select first minimum always
      IF(njc.eq.0)THEN
c close approach ended, but closest approach not recorded; initial close app?
c       WRITE(*,*)' initial close app.? iplam=',iplam
        RETURN
      ELSEIF(njc.gt.1)THEN
         WRITE(*,*)' multiple minima, njc=',njc,tcla(1),tcla(njc),iplam
      ENDIF
c if not, do everything for each minimum
      DO 1 jc=1,njc

c ...but the covariance matrix is not always available
         IF(.not.covava)THEN
c propagation with derivatives, but covariance matrix not available;
c do we need to output the close approach record? Yes, in the old format
c for compatibility reasons
            planam=ordnam(iplam)
            lpla=lench(planam)
            r=vsize(xcla(1,jc))
            rdot=prscal(xcla(1,jc),vcla(1,jc))/r
            call mjddat(tcla(jc),iday,imonth,iyear,hour)
            write(date,'(i4,a1,i2.2,a1,i2.2,f6.5)')
     +           iyear,'/',imonth,'/',iday,hour/24d0
            WRITE(iuncla,100) planam(1:lpla),date,tcla(jc),r,rdot,
     +           (xcla(i,jc),i=1,3),(vcla(i,jc),i=1,3)
 100        FORMAT(a,1x,a16,f12.5,1x,f11.8,e11.3,1x,6(1x,f11.8))
         ELSE
c First parzial derivatives: rewrap vector into 6x6 matrix
            DO i=1,nv
               y2(i)=xcla(i,jc)
               y2(i+nv)=vcla(i,jc)
            ENDDO
            CALL varwra(y2,dxdx1,nv*2,nv)
c Chain rule: we compute dxde=dxdx1*dx1dx0*dx0de
            call mulmat(dxdx1,6,6,dx1dx0,6,6,dxdx0)
            call mulmat(dxdx0,6,6,dx0de,6,6,dxde)
c reference system with tpno as first axis, 
c and the projection of the normal to the ecliptic as third axis
            vc0=vsize(vcla(1,jc))
            DO i=1,3
               tpno(i)=vcla(i,jc)/vc0
            ENDDO
            v2(1)=0.d0
            v2(2)=0.d0
            v2(3)=1.d0
            CALL mtpref(tpno,v2,v3(1,1,jc),vt3(1,1,jc))
c obsolete:
c projection of planet velocity on normal plane to tpno as second axis
c           CALL mtpref(tpno,xplaj(4,jc),v3,vt3)
c rotation to v3 reference system
            CALL mtprot3(.true.,vt3(1,1,jc),xcla(1,jc),vcla(1,jc),dxde,
     +           gc,tpc(1,jc),dtpdet(1,1,jc),sig(1,jc),
     +           axes(1,1,jc),tpr(jc),svv(jc),cxv(jc),czv(jc))
c output variables
            csi=tpc(1,jc)
            zeta=tpc(2,jc)
            stretch=sig(2,jc)
            width=sig(1,jc)
            IF(axes(1,2,jc)*dtpdet(1,1,jc)+axes(2,2,jc)*dtpdet(1,2,jc)
     +           .lt.0.d0)THEN
               axes(1,2,jc)=-axes(1,2,jc)
               axes(2,2,jc)=-axes(2,2,jc)
            ENDIF
            sina=-axes(1,2,jc)
            cosa=axes(2,2,jc)
            alpha=atan2(sina,cosa)
c weak direction
            CALL weakdi(gc,wdir,sdir,-1)
            IF(wdir(1).lt.0.d0)sdir=-sdir
            wtp(1,jc)=prscag(6,wdir,dtpdet(1,1,jc))*sdir
            wtp(2,jc)=prscag(6,wdir,dtpdet(1,2,jc))*sdir
            wtpv(jc)=prscag(6,wdir,dtpdet(1,3,jc))*sdir
            wtpr(jc)=sqrt(wtp(1,jc)**2+wtp(2,jc)**2)
            wtpal(jc)=atan2(-wtp(1,jc),wtp(2,jc))
c date and planet name
            planam=ordnam(iplam)
            lpla=lench(planam)
            numcla=numcla+1
            CALL wriclan3(iuncla,planam(1:lpla),tcla(jc),xcla(1,jc),
     +           vcla(1,jc),csi,zeta,stretch,width,alpha,
     +           moid0,angmoid,wtpr(jc),wtpal(jc),wtpv(jc)
     +           ,svv(jc),cxv(jc),czv(jc))
         ENDIF
 1    ENDDO
      RETURN
      END
c ==============================================================
c COV_AVAI/COV_NOT_AV
c routines to manipulate the covariance common used
c for online target plane analysys
c =============================================================
      SUBROUTINE cov_avai(g)
      IMPLICIT NONE
c INPUT: covariance matrix
      DOUBLE PRECISION g(6,6) 
c HIDDEN OUTPUT: covariance matrix and flag of availability
      INCLUDE 'covariance.h'
      CALL mcopy(6,6,g,gc)
      covava=.true.
      RETURN
      END
      SUBROUTINE cov_not_av
      IMPLICIT NONE
c HIDDEN OUTPUT: flag of availability
      INCLUDE 'covariance.h'
      covava=.false.
      RETURN
      END
c ====================================================
c CLOTEST
c tests for presence of a close approach at the time t0 for 
c asteroid with equinoctal elements east. On output iplnet is the number
c of the approaching planet, 0 if none. 
c texit is delay time advisable to get out of the close approach sphere (in days).
c ********************************************************
c WARNING: this version actually checks only for close approaches to the Earth
c should be fixed to test for all planets.
c=====================================================
      SUBROUTINE clotest(t0,east,iplanet,texit)
      IMPLICIT NONE
c input: time, asteroid elements
      DOUBLE PRECISION t0,east(6)
c output: planet being approached (0 if none), time to get clear
      INTEGER iplanet
      DOUBLE PRECISION texit
c end interface
      DOUBLE PRECISION xea(6),xast(6),dist,enne
c mass of the sun
      INCLUDE 'sunmass.h'
c radius of the close appooach sphere for the planets
      INCLUDE 'parbep.h'
      INCLUDE 'masses.h'
c JPL Earth vector at observation time
      CALL earcar(t0,xea,1)
c cartesian coordinates of the asteroid
      CALL coocha(east,'EQU',gms,xast,'CAR',enne)
c distance 
      dist=sqrt((xea(1)-xast(1))**2+(xea(2)-xast(2))**2+
     +      (xea(3)-xast(3))**2)
c test
      IF(dist.le.dmin(3))THEN
c        WRITE(*,*)' clotest: warning! close approach to Earth'
c        WRITE(*,*)' at the initial epoch, dist=',dist
         iplanet=3
         texit=20.d0
      ELSE
         iplanet=0
         texit=0.d0
      ENDIF
      RETURN
      END
c =================================================
c STEPCON
c numerical stepsize control based on true anomaly
c =================================================
      SUBROUTINE stepcon(x,v,npoin,disc,dir,dt)
      IMPLICIT NONE
c input: geocentric state
      DOUBLE PRECISION x(3),v(3)
c input: min no points, size of MTP disc, direction of time
      INTEGER npoin
      DOUBLE PRECISION disc,dir
c outpput: suggested stepsize
      DOUBLE PRECISION dt
c end interface
      INCLUDe 'trig.h'
      DOUBLE PRECISION prscal,vsize
      DOUBLE PRECISION dtheta,cosdt1,dtheta1,tmin,dtold
      DOUBLE PRECISION x1(3),cgeo(3),jgeo
      dtheta=pig/npoin
      CALL prvec(x,v,cgeo)
      jgeo=vsize(cgeo)
      dt=prscal(x,x)*dtheta/jgeo*dir
      CALL lincom(x,1.d0,v,dt,x1)
      cosdt1=prscal(x,x1)/(vsize(x)*vsize(x1))
      dtheta1=acos(cosdt1)
      IF(dtheta1.gt.dtheta)THEN
         dt=(dtheta/dtheta1)*dt
      ENDIF
      tmin=2*disc/vsize(v)
      dtold=tmin/npoin*dir
c     WRITE(*,*)dt, dtold,vsize(x),vsize(v),dtheta1,jgeo
      RETURN
      END
c =====================================
c FALSI  regula falsi
c =====================================
      SUBROUTINE falsi(t,xa,va,nv,xpla,jc,jt,first,iplam)
      IMPLICIT NONE
c INPUT 
c current position and velocity
      INTEGER nv
      DOUBLE PRECISION t,xa(nv),va(nv),xpla(6)
c is this the beginning of an integration?
      LOGICAL first
c planet being approached
      INTEGER iplam
c INPUT AND OUTPUT
c close approach and target plane counters
      INTEGER jc,jt
c =====================================
c fixed stepsize time interval
      DOUBLE PRECISION tt1,tt2
c data for regula falsi
      DOUBLE PRECISION r1,r2,r0,rdot1,rdot2,rdot0,t1,t2,t0
      DOUBLE PRECISION z1,z2,z0,zdot1,zdot2,zdot0
c functions
      DOUBLE PRECISION vsize,prscal
c state vectors: without partials
      DOUBLE PRECISION x(3),v(3),xt(3),vt(3),xplat(6)
c with partials
      INCLUDE 'nvarx.h'
      DOUBLE PRECISION xat(nvar2x),vat(nvar2x)
c time steps 
      DOUBLE PRECISION dt,tt,di,hh
c iterations
      INTEGER it,i
      INTEGER itmax
      PARAMETER (itmax=50)
c close approach data
      INCLUDE 'closapl.h'
c planetary data
      INCLUDE 'parbep.h'
      INCLUDE 'masses.h'
c ======for MOID========
      DOUBLE PRECISION moid0
c get gms (GM_sun)
      INCLUDE 'sunmass.h'
      include 'trig.h'
c variables for call to compute_minima
      DOUBLE PRECISION x6(6)
c cartesian coordinates asteroid and planet at minimum
      DOUBLE PRECISION cmin(6,20),cplmin(6,20)
c     SQUARED DISTANCE function 
      DOUBLE PRECISION d2(20)
c     number of relative minima found
      INTEGER nummin
c =====================
c output control
      INCLUDE 'proout.h'
c memory model is static
      SAVE
c planetocentric position
      DO i=1,3
        x(i)=xa(i)-xpla(i)
        v(i)=va(i)-xpla(i+3)
      ENDDO
c initialisation at the beginning of each close approach
      IF(first)THEN
         tt1=t
         r1=vsize(x)
         rdot1=prscal(x,v)/r1    
c MOID at the beginning of each encounter
         DO i=1,3
            x6(i)=xa(i)
            x6(i+3)=va(i)
         ENDDO
         CALL compute_minima(x6,xpla,iplam,cmin,cplmin,d2,nummin)
c end initialisation of close approach
         RETURN
      ENDIF
c compute r, rdot
      r2=vsize(x)
      rdot2=prscal(x,v)/r2
      tt2=t
c direction of time propagation
      IF(tt2.gt.tt1)THEN
         di=1.d0
      ELSEIF(tt2.lt.tt1)THEN
         di=-1.d0
      ELSE
         write(*,*)' falsi: zero step'
         write(*,*)tt1,tt2,x,v,xpla
      ENDIF
c check for minimum distance
      IF(abs(rdot2).lt.eprdot)THEN
c already at stationary point; WARNING: is it a minimum?
         CALL strclo(ordnam(iplam),t,xpla,xa,va,nv,jc,r2,rdot2,
     +        d2,cplmin,nummin)
         rdot2=0.d0
      ELSEIF(rdot1*di.lt.0.d0.and.rdot2*di.gt.0.d0)THEN
c rdot changes sign, in a way appropriate for a minimum
         t1=tt1
         t2=tt2
*        WRITE(*,*) ' r. f. begins, t,r,rdot=',t1,r1,rdot1,t2,r2,rdot2
*        WRITE(iuncla,*) ' r. f. begins ',t1,r1,rdot1,t2,r2,rdot2
         dt=-rdot2*(t1-t2)/(rdot1-rdot2)
         tt=dt+t2
         hh=tt-t
c iterate regula falsi
         DO it=1,itmax
           CALL rkg(t,xa,va,nv,hh,xat,vat,xplat)
c planetocentric position
           DO i=1,3
             xt(i)=xat(i)-xplat(i)
             vt(i)=vat(i)-xplat(i+3)
           ENDDO
           r0=vsize(xt)
           rdot0=prscal(xt,vt)/r0
           t0=tt
c selection of next couple of points with opposite sign
           IF(abs(rdot0*dt).lt.eprdot)THEN
c at minimum already
              CALL strclo(ordnam(iplam),tt,xplat,xat,vat,nv,jc,r0,rdot0,
     +             d2,cplmin,nummin)
              GOTO 2
           ELSEIF(rdot0*rdot1.lt.0.d0)THEN
              r2=r0
              rdot2=rdot0
              t2=t0
           ELSEIF(rdot0*rdot2.lt.0.d0)THEN
              r1=r0
              rdot1=rdot0
              t1=t0 
           ENDIF
           dt=-rdot2*(t1-t2)/(rdot1-rdot2)
           tt=dt+t2
           hh=tt-t
         ENDDO
c failed convergence
         CALL strclo(ordnam(iplam),t0,xplat,xat,vat,nv,jc,r0,rdot0,
     +        d2,cplmin,nummin)
         WRITE(ierrou,*)' falsi: failed convergence for ',ordnam(iplam)
         WRITE(ierrou,*)'t1,r1,rdot1,t2,r2,rdot2,dt,tt'
         WRITE(ierrou,111)t1,r1,rdot1,t2,r2,rdot2,dt,tt
c        WRITE(*,*)' falsi: failed convergence for ',ordnam(iplam)
c        WRITE(*,*)'t1,r1,rdot1,t2,r2,rdot2,dt,tt'
c        WRITE(*,111)t1,r1,rdot1,t2,r2,rdot2,dt,tt
  111      format(8(1x,f16.9))
         numerr=numerr+1
      ENDIF
c found zero
 2    CONTINUE
c save current state to be previous state next time
      tt1=tt2
      r1=r2
      rdot1=rdot2
      RETURN
      END
c ===================================================
c STRCLO
c store close approach 
c ===================================================
      SUBROUTINE strclo(planam,tcur,xpla,xa,va,nv,jc,r,rdot,
     +             d2,cplmin,nummin)
      IMPLICIT NONE
c input 
c planet name
      CHARACTER*30 planam
c time current, distance, radial velocity (should be small)
      DOUBLE PRECISION tcur,r,rdot
c cartesian coordinates of planet at minimum
      DOUBLE PRECISION cplmin(6,20)
c     SQUARED DISTANCE function 
      DOUBLE PRECISION d2(20)
c     number of relative minima found
      INTEGER nummin
c  planet position and velocity
      DOUBLE PRECISION xpla(6)
c heliocentrci position and velocity
c  with partial derivatives if nv=21, without if nv=3
      INTEGER nv
      DOUBLE PRECISION xa(nv),va(nv)
c  input/output : close approach minima conter 
      INTEGER jc   
c =========end interface====================
      DOUBLE PRECISION vl,vsize
      INTEGER i
      INTEGER lpla,lench
c calendar date variables
      INTEGER iyear,imonth,iday
      DOUBLE PRECISION hour
      CHARACTER*16 date
c multiple minima analysis
      DOUBLE PRECISION cosa,angle,prscal
      INTEGER j
c common data on close appr. 
      INCLUDE 'cloapp.h'
      INCLUDE 'proout.h'
c trig constants
      INCLUDE 'trig.h'
c store current local moid (just before encounter)
      moid0=sqrt(d2(1))
      angmoid=4.d2
      DO j=1,nummin
         cosa=prscal(cplmin(1,j),xpla)/(vsize(cplmin(1,j))*vsize(xpla))
         angle=acos(cosa)
         IF(angle.lt.angx*radeg)THEN
            moid0=sqrt(d2(j))
            angmoid=angle*degrad
         ENDIF
      ENDDO
c store close approach point
      jc=jc+1
      IF(jc.gt.njcx)THEN
          WRITE(*,*)' strclo: jc>njcx ',jc,njcx
          STOP
      ENDIF
      tcla(jc)=tcur
      rmin(jc)=r
      DO i=1,3
        xcla(i,jc)=xa(i)-xpla(i)
        vcla(i,jc)=va(i)-xpla(3+i)
        xplaj(i,jc)=xpla(i)
        xplaj(i+3,jc)=xpla(3+i)
      ENDDO
      lpla=lench(planam)
      numcla=numcla+1
      call mjddat(tcur,iday,imonth,iyear,hour)
      IF(nv.eq.3)THEN
         write(date,'(i4,a1,i2.2,a1,i2.2,f6.5)')
     +        iyear,'/',imonth,'/',iday,hour/24d0
         WRITE(iuncla,100) planam(1:lpla),date,tcur,r,rdot,
     +      (xcla(i,jc),i=1,3),(vcla(i,jc),i=1,3)
 100     FORMAT(a,1x,a16,f12.5,1x,f11.8,e11.3,1x,6(1x,f11.8))
         WRITE(*,97)planam(1:lpla),date,tcur,r
 97      FORMAT(' Close approach to ',a,' on ',a16,f12.5,' MJD at ',
     +        f10.8,' AU.')
      ELSEIF(nv.eq.21)THEN
c store partials at close approach time
         DO i=4,21
            xcla(i,jc)=xa(i)
            vcla(i,jc)=va(i)
         ENDDO
c planet velocity is needed to define reference system
         DO i=1,6
           xplaj(i,jc)=xpla(i)
         ENDDO
      ELSE
         write(*,*)' strclo: nv=',nv
         STOP
      ENDIF
      RETURN
      END
c =================================================
c RKG
c  Runge-Kutta-Gauss, used as a pure single step
c =================================================
      SUBROUTINE rkg(t1,xa,va,nv,h,xat,vat,xplat)
      IMPLICIT NONE
c INPUT
c time, stepsize
      DOUBLE PRECISION t1,h
c state: position, velocities
      INTEGER nv
      DOUBLE PRECISION xa(nv),va(nv)
c OUTPUT
c state: position, velocities, approached planet position
      DOUBLE PRECISION xat(nv),vat(nv),xplat(6)
c END INTERFACE
      INCLUDE 'nvarx.h'
      DOUBLE PRECISION y1(nvarx),dery(nvarx),yat(nvarx),de
      INCLUDE 'parint.h'
      INCLUDE 'comint.h'
      INCLUDE 'rkcoef.h'
      INTEGER nvar
c     PARAMETER (nvar=6)
      DOUBLE PRECISION ep(itmax),ck(ismax,nvarx),t(ismax)
      INTEGER i,j,it,jj
      INTEGER imem,ips
c control of derivatives: left as it is
c     INTEGER ide,ideold
c     COMMON/deriv/ide
c     ideold=ide
c     ide=0
c state vector length
      nvar=2*nv
c state vector before the step
      DO i=1,nv
        y1(i)=xa(i)
        y1(i+nv)=va(i)
      ENDDO
c set intermediate times; array ck initialized to zero
      DO j=1,isrk
        t(j)=t1+h*c(j)
        DO i=1,nvar
          ck(j,i)=0.d0
        ENDDO
      ENDDO
c controls on convergence inititalized to zero
      DO  it=1,itmax
        ep(it)=0.d0
      ENDDO
c  gauss-seidel iterations for array ck
      it=1
      ips=0
c  main loop on intermediate points
 1    DO 11 j=1,isrk
        DO 12 i=1,nvar
          de=0.d0
          DO  jj=1,isrk
            de=de+a(j,jj)*ck(jj,i)
          ENDDO
          yat(i)=de*h+y1(i)
 12     continue
c memory location not used by ra15v
        imem=j+10  
        CALL fctcl(t(j),yat,dery,nvar,xplat,ips,imem)
        DO i=1,nvar
          ep(it)=ep(it)+dabs(dery(i)-ck(j,i))
        ENDDO
        DO i=1,nvar
          ck(j,i)=dery(i)
        ENDDO
 11   continue
c  control on end of gauss-seidel iterations
      ep(it)=ep(it)/isrk
      IF(it.gt.1.and.ep(it).gt.ep(it-1)*1.1d0)THEN
         WRITE(*,*)' rkg: stop at iteration ',it, ' before too late'
         WRITE(*,*)t1,(ep(jj),jj=it-1,it)
         GOTO 77
      ENDIF
c  new state vector y3
      DO i=1,nvar
        de=0.d0
        DO  j=1,isrk
          de=de+b(j)*ck(j,i)
          yat(i)=y1(i)+h*de
        ENDDO
      ENDDO
      IF(ep(it).gt.eprk)THEN
         IF(it.ge.lit1)THEN
c  too many gauss-seidel iterations
            WRITE(*,*)' rkg: non convergent after ',it,' iterations'
            WRITE(*,*)t1,ep(it)
         ENDIF
         it=it+1
         ips=-1
         GOTO 1
      ENDIF
c
c right hand side at new state
 77   CALL fctcl(t1+h,yat,dery,nvar,xplat,0,1)
c copy pos. vel.
      DO i=1,nv
        xat(i)=yat(i)
        vat(i)=yat(i+nv)
      ENDDO
c control of derivatives reset: NO
c     ide=ideold
      RETURN
      END
c ==========================================================
c   FCTCL
c reduction to first order, to use with RKG
c ==========================================================
c   subroutine secondo membro
c   equazioni ridotte all'ordine 1
c   a partire dall'eq del secondo ordine
c   calcolata da force
      SUBROUTINE fctcl(t,y,dery,nvar,xxpla,ips,imem)
      implicit none
      integer nvar
      double precision y(nvar),dery(nvar)
      double precision xxpla(6)
      double precision t
      INTEGER ips,imem
c end interface
      integer nvar2,idc,i
c****************
c   static memory not required
c****************
      nvar2=nvar/2
      call force(y,y(nvar2+1),t,dery(nvar2+1),nvar2,idc,xxpla,ips,imem)
      do  i=1,nvar2
        dery(i)=y(nvar2+i)
      enddo
      return
      end

c ================================================
c WRICLAN
c writes close approach record with linear MTP analysis
c ====================================================
      SUBROUTINE wriclan3(iuncla,planam,tcla,xcla,vcla,
     +     csi,zeta,stretch,width,alpha,moid0,angmoid,
     +     wtpr,wtpal,wtpv,svv,cxv,czv)
      IMPLICIT NONE
c INPUT
c unit for output
      INTEGER iuncla
c planetocentric position and velocity, MJD time of closest approach
      DOUBLE PRECISION xcla(3),vcla(3),tcla
c target plane coordinates, 1-sigma axes of target plane ellipse and
c angle between zeta axis and long axis (sign chosen with increasing a)
      DOUBLE PRECISION csi,zeta,stretch,width,alpha
c name of approached planet
      CHARACTER*(*) planam
c LOCAL MOID at beginning of encounter, angle between planet  positions
      DOUBLE PRECISION moid0,angmoid
c projection of LOV on target plane, polar coordinates, velocity component
      DOUBLE PRECISION wtpr,wtpal,wtpv
c third row of covariance on extended MTP
      DOUBLE PRECISION svv,cxv,czv
c hidden input: dimensions of array
      INCLUDE 'maxclo.h'
c hidden output: array written on .clo file
      DOUBLE PRECISION arrline(ncol)
c date, planet name
      CHARACTER*16 date
c end interface
c distance, radial velocity (should be small,just for check) 
      DOUBLE PRECISION r,rdot,vsize,prscal
c calendar date variables
      INTEGER iyear,imonth,iday
      DOUBLE PRECISION hour
      INTEGER j
c batch operations?
      INCLUDE 'comdif.h'
c trig constants
      INCLUDE 'trig.h'
c ====================================================
      r=vsize(xcla)
      rdot=prscal(xcla,vcla)/r
      call mjddat(tcla,iday,imonth,iyear,hour)
      write(date,'(i4,a1,i2.2,a1,i2.2,f6.5)')
     +     iyear,'/',imonth,'/',iday,hour/24d0
c copy in array
      arrline(1)=tcla
      arrline(2)=r
      arrline(3)=rdot
      DO j=1,3
        arrline(j+3)=xcla(j)
        arrline(j+6)=vcla(j)
      ENDDO
      arrline(10)=csi
      arrline(11)=zeta
      arrline(12)=stretch
      arrline(13)=width
      arrline(14)=alpha*degrad
      arrline(15)=moid0
      arrline(16)=angmoid
      arrline(17)=wtpr
      arrline(18)=wtpal*degrad
      arrline(19)=wtpv
      arrline(20)=svv
      arrline(21)=cxv
      arrline(22)=czv
      arrline(23)=vsize(vcla)
c file output
      WRITE(iuncla,100) planam,date,(arrline(j),j=1,23)
 100  FORMAT(a15,1x,a16,f12.5,1x,f11.8,e11.3,1x,6(1x,f11.8)
     +     ,2(1x,f11.8),2(1x,1p,e12.4),1x,0p,f10.5,1x,f10.6,1x,f7.2,
     +     1x,1p,e12.4,0p,1x,f10.5,1x,1p,e12.4,1x,e12.4,1x,e12.4,
     +     1x,e12.4,1x,e12.4)
c standard output
      IF(.not.batch)THEN
      WRITE(*,97)planam,date,tcla,r,csi,zeta,stretch,width,alpha*degrad
     +        ,moid0,angmoid,
     +      wtpr,wtpal*degrad
 97   FORMAT(' Close approach to ',a,' on ',a16,f12.5,' MJD at ',
     +     f10.8,' AU.'/2(1x,f11.8),2(1x,1p,e12.4),1x,0p,f10.5,1x,f10.6,
     +     1x,f7.2,1x,1p,e12.4,0p,1x,f10.5)
      ENDIF
      RETURN
      END
c =================================================================
c REACLOREC reads the next close approach record of the required planet
c reqpla; if reqpla=' ', then all planets are required
c ==================================================================
      SUBROUTINE reaclorec(iunclo,reqpla,imulcur,arrline,
     +     planam,error,eof)
      IMPLICIT NONE
c INPUT
c required planet
      CHARACTER*15 reqpla
c input unit
      INTEGER iunclo
c OUTPUT
c planet found
      CHARACTER*15 planam
c arrays with close approach data
      INTEGER imulcur
      INCLUDE 'maxclo.h'
      DOUBLE PRECISION arrline(ncol)
c error flag, end of file
      LOGICAl error,eof
c end interface
      CHARACTER*512 record
      INTEGER ii
      INTEGER le,le1
c ---------------------------
      eof=.false.
      error=.false.
 1    READ(iunclo,101,end=2)record
 101  FORMAT(a)
      CALL rmsp(reqpla,le)
      ii=index(record,'mult_')
      IF(ii.ne.0)THEN
         IF(imulcur.gt.0)THEN
            READ(record(ii+5:),102)imulcur
 102        FORMAT(i4)
            GOTO 1
         ELSE
            WRITE(*,*)' reaclorec: multiple solution clo file',record
            GOTO 1
         ENDIF
      ELSE
         ii=index(record,reqpla(1:le))
c         WRITE(*,*)ii,reqpla, record
         IF(ii.ne.0.and.le.ne.0)THEN
            planam=reqpla
            CALL reaclan3(record,arrline,error)
         ELSEIF(le.eq.0)THEN
            READ(record,103)planam
            CALL rmsp(planam,le1)
 103        FORMAT(a15)
            CALL reaclan3(record,arrline,error)
         ELSE
            GOTO 1
         ENDIF
      ENDIF
      RETURN  
 2    eof=.true.
      RETURN
      END
c ================================================
c REACLAN3 reads record written by wriclan3, only data part, places in array
      SUBROUTINE reaclan3(record,arrline,error)
      IMPLICIT NONE
c INPUT
      CHARACTER*512 record
      INCLUDE 'maxclo.h'
c OUTPUT
      LOGICAL error     
      DOUBLE PRECISION arrline(ncol)
c
      INTEGER j
      error=.false.
      READ(record, 100, err=2) (arrline(j),j=1,23)
 100  FORMAT(32x,f12.5,1x,f11.8,e11.3,1x,6(1x,f11.8)
     +     ,2(1x,f11.8),2(1x,1p,e12.4),1x,0p,f10.5,1x,f10.6,1x,f7.2,
     +     1x,1p,e12.4,0p,1x,f10.5,1x,
     +     1p,e12.4,1x,e12.4,1x,e12.4,1x,e12.4,1x,e12.4)
      RETURN
 2    error=.true.
      RETURN
      END
c =================================================================
c reference system with tpno as first axis, +dx as third axis
c to be used with dx= normal to ecliptic
      SUBROUTINE mtpref(tpno,dx,v3,vt3)
      IMPLICIT NONE
      DOUBLE PRECISION tpno(3),dx(3),v3(3,3),vt3(3,3)
c end interface
      DOUBLE PRECISION vl,vsize,vv,prscal
      INTEGER i
c the first vetor is normal to the target plane
      CALL vcopy(3,tpno,v3(1,1))
c the third vector is the projection of -dx on the plane orthogonal to tpno       
      vv=prscal(tpno,dx)
      DO i=1,3
         v3(i,3)=-(dx(i)-vv*tpno(i))
      ENDDO
c reduced to unit length
      vl=vsize(v3(1,3))
      DO i=1,3
         v3(i,3)=v3(i,3)/vl
      ENDDO
c the second vector must be orthogonal to both
      CALL prvec(tpno,v3(1,3),v3(1,2))
c paranoia check of unit length
      vl=vsize(v3(1,2))
      DO i=1,3
         v3(i,2)=v3(i,2)/vl
      ENDDO
c the inverse is the transposed matrix
      CALL transp(v3,3,3,vt3)
      RETURN
      END
c =================================================================
c MTPROT3
c rotation to v3 reference system, but with information
c enough to convert to Opik Target plane
      SUBROUTINE mtprot3(batchcl,vt3,dx,dv,dxde,gc,
     +        tpc,dtpdet,sig,axes,tpr,svv,cxv,czv)
c    + tcla,xpla not needed unless we compute yddot
      IMPLICIT NONE
c input
      LOGICAL batchcl
c vt3= matrix to change coord to MTP reference (V_1=velocity unit vector;
c      V_3= projection of some axis, e.g. orthogonal to ecliptic, on MTP
c      V_2=V_1 x V_3; note reversal of orientation)
c dx,dv= geocentric position and velocity, in ecliptic reference frame
c dxde= jacobian of dx,dv with respect to elements
c gc  = covariance matrix of orbital elements
      DOUBLE PRECISION vt3(3,3),dx(3),dv(3),dxde(6,6),gc(6,6)
c approached planet: coordinates, time of close approach
c     DOUBLE PRECISION xpla(6),tcla
c output 
c tpc= 3-d MTP vector: x,z, V
c dtpdet= jacobian of tpc w.r. to elements, trasnposed
c sig= semiaxes, axes= unit vector of axes of MTP confidence ellipse (sigma=1)
c tpr= minimum distance 
c svv= rms of velocity V
c cxv,czv= correlation of V with x,z
      DOUBLE PRECISION tpc(3),dtpdet(6,3),sig(2),axes(2,2),tpr
      DOUBLE PRECISION svv,cxv,czv
c end interface
c xx,vv= position and velocity in MTP reference system
c dxde3,dvde3= reduced jacobians (only for multiplication)
c dtpcde= jacobian of tpc w.r. to elements
      DOUBLE PRECISION xx(3),vv(3),dxde3(3,6),dvde3(3,6),dtpcde(3,6)
c tpth= angle on MTP
      DOUBLE PRECISION tpth
c dxxde,dvvde= jacobians w.r. to elements but in the MTP reference frame
      DOUBLE PRECISION dxxde(3,6),dvvde(3,6)
c xhel,vhel= helicoentric position and velocity
c f=heliocentric acceleration
c yddot= component of geocentric acceleration in the velocity direction
c xxpla= position of the planet
c     DOUBLE PRECISION xhel(3),vhel(3),f(3),xxpla(6),yddot
c     INTEGER idc
c gxz= 2x2 covariance matrix on MTP
c gmtp=3x3 covariance matrix on MTP
c tmp1= workspace for multiplications
      DOUBLE PRECISION gxz(2,2),gmtp(3,3),tmp1(3,6)
c eigenvalues, workspace, transposed
      DOUBLE PRECISION eigval(2),tmp26(2,6),fv1(2),fv2(2),dadde(2,6)
c error flag
      INTEGER ierr
c loop index
      INTEGER i,ii,jj
      DOUBLE PRECISION vl,vsize,bsd,d,v,prscal,tmp
c computation of MTP coordinates and derivatives
      CALL mulmav (vt3,3,3,dx,3,xx)
      IF(.not.batchcl)WRITE(*,*)' distance from target plane=',xx(1)
      CALL mulmav (vt3,3,3,dv,3,vv)
      DO i=1,3
        DO ii=1,6
          dxde3(i,ii)=dxde(i,ii)
          dvde3(i,ii)=dxde(i+3,ii)
        ENDDO
      ENDDO
      CALL mulmat (vt3,3,3,dxde3,3,6,dxxde)
      CALL mulmat (vt3,3,3,dvde3,3,6,dvvde)
c MTP x,z coordinates and derivatives with respect to elements
      DO ii=1,2
        tpc(ii)=xx(ii+1)
        DO jj=1,6
          dtpcde(ii,jj)=dxxde(ii+1,jj)-(vv(ii+1)/vv(1))*dxxde(1,jj)
        ENDDO
      ENDDO
c MTP polar coordinates
      tpr=sqrt(tpc(1)**2+tpc(2)**2)
c velocity and derivatives with respect to elements
      v=vsize(vv)
c      DO jj=1,3
c       xhel(jj)=dx(jj)+xpla(jj)
c       vhel(jj)=dv(jj)+xpla(3+jj)
c     ENDDO
c     CALL force(xhel,vhel,tcla,f,3,idc,xxpla,0,10)
c      yddot=prscal(f,vv)/v
c add velocity as third coordinate
      tpc(3)=vv(1)
      DO jj=1,6
          dtpcde(3,jj)=dvvde(1,jj)
c two-body approx: the acceleration of the Earth does not matter
c for the other bodies, only differential accel. which is neglected
c                                  -(yddot/vv(1))*dxxde(1,jj)
      ENDDO
      CALL transp(dtpcde,3,6,dtpdet)
c output for interactive mode
      IF(.not.batchcl)THEN
         tpth=atan2(tpc(2),tpc(1))
         WRITE(*,170) tpc,tpr,tpth
 170     FORMAT('MTP coordinates ',2f12.8,' dist ',f12.8,' angle ',f8.4)
         WRITE(*,*)' partial derivatives'
         DO ii=1,3
            WRITE(*,171) (dtpcde(ii,jj),jj=1,6)
         ENDDO
 171     FORMAT(6(f12.5,2x))
      ENDIF
c ========================================================
c compute target ellipse of confidence
      CALL mulmat(dtpcde,3,6,gc,6,6,tmp1)
      CALL mulmat(tmp1,3,6,dtpdet,6,3,gmtp)
      CALL marg2(gmtp,gxz,svv,cxv,czv)
c compute ellipse of confidence
c eigenvalues
      CALL rs(2,2,gxz,eigval,1,axes,fv1,fv2,ierr)
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
      SUBROUTINE marg2(gmtp,gxz,svv,cxv,czv)
      IMPLICIT NONE
c INPUT 
c gmtp= 3x3 covariance matrix on extended MTP
      DOUBLE PRECISION gmtp(3,3)
c OUPUT
c gxz= 2x2 covarioance matrix on MTP
c svv= rms of velocity V
c cxv,czv= correlation of V with x,z
      DOUBLE PRECISION gxz(2,2),svv,cxv,czv
      INTEGER i,j
      DO i=1,2
        DO j=1,2
          gxz(i,j)=gmtp(i,j)
        ENDDO
      ENDDO
      svv=sqrt(gmtp(3,3))
      cxv=gmtp(1,3)/sqrt(gmtp(3,3)*gmtp(1,1))
      czv=gmtp(2,3)/sqrt(gmtp(3,3)*gmtp(2,2))
      RETURN
      END
