c ===========MODULE FORCE_MODEL===============
c PUBLIC ROUTINES
c               rmodel    read options and initialize
c               selmet    adjust force model for orbit
c               radar_ob  to adjust model if there is radar
c MODULE CONTAINS
c ROUTINES
c               force     right hand side
c                 planast  planetary and asteroid positions
c                 j2sun    sun oblateness perturbation
c                 genrel   general relativity, sun's contribution only
c                 eihrel   general relativity, n-body
c                 yarkdi   yarkovsky diurnal
c                 yarkse   yarkovsky seasonal
c               yarkinit   initialisation of yarkovsky
c               selpert    avoiding self perturb. of massive asteroids
c
c HEADERS
c force_model.o: model.h iclap.h proout.h selast.h \
c               comint.h parbep.h masses.h jplhdr.h \
c               vlight.h sysdep.h
c
c               model.h ????
c               bifina.h, combep.h with jpl_ephem.f
c               yarkov.h, yarkom.h internal
c               parbep.h (jpl_ephem.f)
c               combep.h (jpl_ephem.f)
c               sunmass.h, trig.h public
c               closapl.h 
c               radar.h only in the 3 public routines
c               
c
c ======================================================
c  RADAR_OB
c  sets the radar flag if there is radar data, to be used
c  by selmet to adjust the force model
c ======================================================
      SUBROUTINE radar_ob(iobs,m)
      IMPLICIT NONE
c INPUT
      INTEGER m
      INTEGER iobs(m)
c HIDDEN OUTPUT
      INCLUDE 'radar.h'
c END INTERFACE
      INTEGER j
      radar=.false.
      DO j=1,m
         IF(iobs(j).ge.2000.and.iobs(j).le.2999)radar=.true.
      ENDDO
      RETURN
      END
* Copyright (C) 1998 by OrbFit Consortium
* Version: December 14, 1998 Mario Carpino
c ================================================================
c RMODEL
c ================================================================
c this subroutine reads from file 'propag.def' all the propagator options;
c the only input argument is the string "run", used to open 
c output files with appropriate names
c  WARNING: we advise the user against changing the file propag.def;
c                do it at your risk !....
c COMMON STORED:
c    model.h : for dinamical model, numerical integrators parameters
c                and controls.
c input options for phisical model! version 1.2 of propagator
c       ilun=0!  0=no moon  1= yes 
c       imerc=1! 0=no mercury 1=yes (recommended =1 for asteroids)
c       iplut=0! 0=no pluto 1=yes (recommended =0 for asteroids, =1 for 
c                                   transneptunian)
c       irel=0!  0=newtonian 1=gen. relativity (recommended =0 for main belt,
c                                   =1 for Earth crossing)
c       iast=0!  0=no asteroids with mass n=no. of massive asteroids 
c       filbe='CPV'! name of the asteroid ephemerides file (written by BINEPH)
c       iclap=1! 0=no close approach control 1=yes (recommended =1)
c       iaber=1! aberration 0=no 1=yes (recommended =1 always)
c       istat=1! 0=no topocentric corr 1= yes (recommended =1 always)
c
c       iclap = 1 means for Radau alg. to calculate time pos. and vel.
c                at closest approach; for multistep alg. the propagator just
c                detects the time of close-appr.
c      iclap =0 the subroutine force does not check for possible close
c                approaches to major planets andor massive asteroids. 
c       iyark=1 ! Yarkovski force 
c       iyarpt=0! partials of yarkovski are generally not computed
c
c   cloapp.h: npoint = minimal number of close-approach records
c              ndeg  = degree of polimomial interpolation 
c                        (ndeg < 2*npoint-1);
c   The minimal distance to define a close-approch is controlled by:
c              dmea (for the Earth only), dter (for Mercury, Venus and Mars)
c              dmoon (for the moon), dmjup (for giant planets), 
c              dmast (for massive asteroids). 
c ================================================================
c ========INTERFACE============================
      SUBROUTINE rmodel(run)
      implicit none
      character*80 run
c ========END INTERFACE========================
      include 'model.h'
      include 'iclap.h'
      include 'bifina.h'
      include 'proout.h'
      include 'parbep.h'
      include 'combep.h'
      include 'selast.h'
      INCLUDE 'closapl.h'
      INCLUDE 'yarkov.h'
      INCLUDE 'yarkom.h'
      INCLUDE 'radar.h'
      INCLUDE 'verbosity.h'
      character*80 filbe
c     character*80 file
      logical fail,fail1,found
      integer ll,ia
c controls for bizarre orbits
      DOUBLE PRECISION ecclim, samin,samax,phmin,ahmax
c****************
c   static memory not required (used only once)
c****************
* read time-range
      call trange
c restart flag set to default value
      CALL set_restart(.true.)
* read options
* force model flags
      fail=.false.
      call rdnint('propag.','ilun',ilun,.true.,found,fail1,fail)
      call rdnint('propag.','imerc',imerc,.true.,found,fail1,fail)
      call rdnint('propag.','iplut',iplut,.true.,found,fail1,fail)
      call rdnint('propag.','irel',irel,.true.,found,fail1,fail)
c in case selmet is not called
      icrel=irel
      call rdnint('propag.','iast',iast,.true.,found,fail1,fail)
      call rdncha('propag.','filbe',filbe,.true.,found,fail1,fail)
      call rmsp(filbe,ll)
      if(ll.le.0)STOP '**name of asteroid ephem file is wrong**'
      filbep=filbe(1:ll)//'.bep'
      filbec=filbe(1:ll)//'.bai'
      DO ia=1,iast
        astid(ia)=ia
      ENDDO
      iatrue=iast
c close approach control
      call rdnint('propag.','iclap',iclap,.true.,found,fail1,fail)
      iorb=11
      call rdnint('propag.','iaber',iaber,.true.,found,fail1,fail)
      call rdnint('propag.','istat',istat,.true.,found,fail1,fail)
      eprdot=1.d-10
c non gravitational perturbations
      call rdnint('propag.','iyark',iyark,.true.,found,fail1,fail)
      call rdnint('propag.','iyarpt',iyarpt,.true.,found,fail1,fail)
      call rdncha('propag.','yardir',yardir,.true.,found,fail1,fail)
      yarfil=.false.
      yarini=.false.
c radar flag defaults to false
      radar=.false.
* numerical integrator options
      call rmsp(run,ll)
      call inipro
*
      if(iclap.ne.0) then
* close approach control if requeststed
        call rdnint('propag.','npoint',npoint,.true.,found,fail1,fail)
        call rdnrea('propag.','dmea',dmea,.true.,found,fail1,fail)
        call rdnrea('propag.','dmoon',dmoon,.true.,found,fail1,fail)
        call rdnrea('propag.','dmjup',dmjup,.true.,found,fail1,fail)
        call rdnrea('propag.','dmast',dmast,.true.,found,fail1,fail)
        call rdnrea('propag.','dter',dter,.true.,found,fail1,fail)
      endif
* Options for difcor (including outlier rejection)
      call difini
      CALL rejini
* Options for stopping difcor at bizarre orbits; use default
      ecclim=0.d0
      samin=0.d0
      samax=0.d0
      phmin=0.0d0
      ahmax=0.d0
      CALL bizset(ecclim,samin,samax,phmin,ahmax)
c availability of covarinace matrix is false at start
      CALL cov_not_av
c verbosity is set at the minimum level by default
      verb_clo=1
      verb_pro=1
      verb_dif=1
      verb_mul=1
      verb_rej=1
c
      if(fail) stop '**** rmodel: abnormal end ****'
      return
      end
* version 1.8.3 A. Milani Feb 1999
* *** SELMET ***
* Choice of numerical integration method
*
* This routine is called from propag only if  imet=0
* in 'namerun'.top.
*            INPUT: equinoctal elements of asteroid
*            OUTPUT : icmet etc. stored in model.h
      SUBROUTINE selmet(eq)
      implicit none
      include 'model.h'
      include 'iclap.h'
      include 'comint.h'
      include 'proout.h'
      INCLUDE 'radar.h'
      INCLUDE 'verbosity.h'
* asteroid elements, eccentricity, perielion, aphelion
      double precision eq(6),ecc,q,qg
      integer iord,iork
* controls to define a main belt asteroid
      double precision qmin,qgmax,eccmax 
**********************************
c ==========modification 28/10/2000========
c trick to save the ilun of the option file, 
c being free to change it for the current case
      integer lflag,ilunold
*  static memory only for:
      save qmin,qgmax,eccmax,ilunold,lflag
      data lflag /0/
c ========================================
**********************************
* selection of minimum perihelion: all NEO 
      qmin=1.3d0
* selection of maximum aphelion: almost Jupiter crossing
      qgmax=4.3d0
* selection of maximum eccentricity for multistep
      eccmax=0.25d0
* current values of the elements
      ecc=sqrt(eq(2)**2+eq(3)**2)
      q=eq(1)*(1.d0-ecc)
      qg=eq(1)*(1.d0+ecc)
c =========================3/9/2001
c a priori setup        
c relativity flag set to the choice done in the option file
      icrel=irel
c lunar flag set to  the choice done in the option fi
      if(lflag.eq.0)then
         ilunold=ilun
         lflag=1
      endif
      ilun=ilunold
c if there are radar data, the dynamical model is always the same 
c and the integrator must be radau anyway
      IF(radar)THEN
         icrel=2
         ilun=1
         imerc=1
         llev=12
         iast=3
         icmet=3
c ==========================end modif.
c asteroids
      ELSEIF(eq(1).lt.5.d0)THEN
c pluto is irrelevant anyway
         iplut=0
c mercury is required
         imerc=1
* main belt case
         if(q.gt.qmin.and.qg.lt.qgmax.and.ecc.lt.eccmax)then
c physical model
c numerical integration method
            icmet=1
            iord=8
            mms=iord-2
            iork=8
            isrk=iork/2
            epms= 1.0d-12
            lit1= 10
            lit2= 4
c non mainbelt: NEA, Mars crosser, high eccentricity, Jupiter crossing
         else
c physical model
            if(q.lt.qmin.and.irel.lt.1)then
               icrel=1
            endif
c ==================28/10/00
            if(q.lt.qmin.and.ilun.eq.0)then
               ilunold=ilun
               ilun=1
            endif
c ================end modif.
c numerical integration method
            icmet=3
            llev=12
         endif
c EKO
      ELSEIF(qg.gt.32.d0.and.ecc.lt.eccmax)THEN
c physical model
         iplut=1
         icrel=0
         imerc=0
c numerical integration method
         icmet=1
         iord=8
         mms=iord-2
         iork=8
         isrk=iork/2
         epms= 1.0d-13
         lit1= 10
         lit2= 4
c Trojans
      elseif(q.gt.4.d0.and.qg.lt.7.d0)then
c physical model
         iplut=0
         icrel=0
         imerc=0
c numerical integration method
         icmet=1
         iord=8
         mms=iord-2
         iork=8
         isrk=iork/2
         epms= 1.0d-12
         lit1= 10
         lit2= 4
c centaurs
      else
         icmet=3
         icrel=0
         iplut=1
c         llev=9
      endif
c write option selected
      IF(verb_pro.gt.10)THEN
         write(ipirip,*) 'a,e= ',eq(1),ecc,' icmet=',icmet
         write(ipirip,*) 'hms,epms,eprk,deltos,error, mms,isrk',
     +        'lit1,lit2,iusci,icha,llev,hev'
         write(ipirip,*) hms,epms,eprk,deltos,error, mms,isrk
     +        ,lit1,lit2,iusci,icha,llev,hev
         write(ipirip,*)'Force: ilun,imerc,iplut,irel,iast,iaber,istat'
         write(ipirip,*) ilun,imerc,iplut,irel,iast,iaber,istat
         write(ipirip,*) 'Close-approaches param; iclap=',iclap
         if(iclap.ne.0)then
            write(ipirip,*) 'dmea,dmoon,dmjup,dmast', 
     +      dmea,dmoon,dmjup,dmast
         endif
      ENDIF
      return
      end
c
c ===========================================================
c FORCE : accelerations acting on a massless asteroid
c ===========================================================
c version 2.2.9; A. Milani, January 9, 2002
c with no recomputation of JPL ephemerides
c if ips=0 reconsult JPL ephem (and asteroid .bep file)
c       and store in location imem
c if ips>0 use JPL ephem in location imem
c           and recompute all 
c if ips<0 use JPL ephem in location imem 
c           partial recomputation (not implemented yet).
c Accelerations on the asteroid, computed by using
c JPL ephemerides as source for the position of the planets
c (variational equation second part)
c close approach control
c asteroid perturbations with trick to avoid self-perturbations
c ===========================================================
      SUBROUTINE force(x,v,t0,f,nd,idc,xxpla,ips,imem)
      implicit none
c model controls
      include 'model.h'
      include 'iclap.h'
c asteroid maximum number
      include 'parbep.h'
c planetary masses and other model parameters, including asteroid masses
      include 'masses.h'
c controls and physical parameters for nongravitational forces
      INCLUDE 'yarkom.h'
      INCLUDE 'yarkov.h' 
c ======INPUT===================
c dimension of position vector
      integer nd
c Position, velocity,  time
      double precision x(nd), v(nd) ,t0
c flag for recomputation, memory location
      INTEGER ips,imem
c ======OUTPUT===================
c acceleration
      double precision f(nd)
c Positions and vel of the planet involved in  close-app
c stored only if idc.ne.0
      integer idc
      double precision xxpla(6)
c ======END INTERFACE==============
c flag for interpolation of position only or pos and vel.
      integer istate
      common/cstate/istate
c control of derivatives
      integer ide
      common/deriv/ide
c ================================================
c Distances 
      double precision r(nmassx), d(nmassx),rast
      double precision derf(3,3),dfb(3,3),dfc(3,3)
c Positions of the planets also vel.; space also for asteroids
      double precision xpla(6,nmassx)
c Relativistic perturbations
      double precision drgr(3,7),frel(3)
c Yarkovsky force
      double precision yarkv(21)
c loop indexes i=1,npla j=1,3
      integer i,j,k,ir,ic,icqfw
c scalar temporaries
      double precision sum,var1
c ===========================================================
c JPL ephemerides, and asteroid ephemerides if required
      CALL planast(t0,ips,imem,istate,xpla)
c ===========================================================
c Computation of planet vector lengths
      do 20 i=1,nmass
         r(i)=sqrt(xpla(1,i)**2+xpla(2,i)**2+xpla(3,i)**2)
 20   continue
      rast=sqrt(x(1)**2+x(2)**2+x(3)**2)
c ===========================================================
c Computation of planets-asteroid distances
      idc=0
      do 30 i=1,nmass
c       d(i)=vsize(x-xpla(1:3,i))
        d(i)=sqrt((x(1)-xpla(1,i))**2+(x(2)-xpla(2,i))**2+
     +             (x(3)-xpla(3,i))**2)
        if(iclap.gt.0)then
           if(d(i).lt.dmin(i))then
              if(idc.eq.0)then
                 idc=i
                 DO  icqfw=1,3
                   xxpla(icqfw)=xpla(icqfw,i)
                   IF(istate.eq.2)xxpla(icqfw+3)=xpla(icqfw+3,i)
                 ENDDO
              else
                 write(*,*)' force: this should not happen',t0, idc,i
                 write(*,*)nmass,(d(j),dmin(j),j=1,nmass)
                 stop
              endif
           endif
        endif
30    continue
c ===========================================================
c initialize force
      DO j=1,3
        f(j)=0.d0
      ENDDO
c ===========================================================
c general relativistic correction
      if(icrel.eq.1)then
         call genrel(x,v,drgr)
         do j=1,3
            f(j)=f(j)+drgr(j,1)
         enddo
      elseif(icrel.eq.2)then
         call eihrel(x,v,xpla,d,r,rast,frel)
         do j=1,3
            f(j)=f(j)+frel(j)
         enddo
c if we need the refined relativity then also include J_2 for the sun
         call j2sun(x,frel)
         do j=1,3
            f(j)=f(j)+frel(j)
         enddo
      endif
c ===========================================================
c Sitarski force (Acta Astronomica vol. 48 (1998), pp. 547-561)     
c         do j=1,3
c            f(j)=f(j)+(-0.15987d-10/2d0/2.5119760)*v(j)
c         enddo
c ===========================================================
c yarkovsky effect, if required and data are avilable
      if(iyark.ge.1.and.yarfil)then
         IF(.not.yarini)THEN
            WRITE(*,*)' contradiction in non gravitational parameters'
            WRITE(*,*)'iyark=',iyark,' yarfil=',yarfil,' yarini=',yarini
            STOP 
         ENDIF
c diurnal
         call yarkdi(x,yarkv,iyarpt)
         do j=1,3
            f(j)=f(j)+yarkv(j)
         enddo
c seasonal
         if(iyark.gt.1)then
            call yarkse(x,v,yarkv,iyarpt)
            do j=1,3
               f(j)=f(j)+yarkv(j)
            enddo
         endif
      endif
c ===========================================================
c Computation of indirect force FI
      DO i=1,nmass
        DO j=1,3
          f(j)=f(j)-gm(i)/r(i)**3*xpla(j,i)
        ENDDO
      ENDDO
c ===========================================================
c Adding planets-asteroid attractions
      DO i=1,nmass
        DO j=1,3
            f(j)=f(j)+gm(i)/d(i)**3*(xpla(j,i)-x(j))
        ENDDO
      ENDDO
c ===========================================================
c Adding solar attraction
      DO j=1,3
        f(j)=f(j)-gm0*x(j)/rast**3
      ENDDO
c ===========================================================
      if(ide.lt.1)return
      IF(nd.eq.3)RETURN
c Computation of partial derivatives matrix
      do 1 ir=1,3
         do 2 ic=ir,3
           var1=3.d0*gm0*x(ir)*x(ic)/rast**5
           if(ir.eq.ic)var1=var1-(gm0/(rast**3))
           sum=0.d0
           do  3 i=1,nmass
             sum=sum+3.d0*gm(i)*
     +             (xpla(ir,i)-x(ir))*(xpla(ic,i)-x(ic))/(d(i)**5)
             if(ir.eq.ic)sum=sum-gm(i)/(d(i)**3)
 3         continue
           derf(ir,ic)=var1+sum
           if(ir.ne.ic) then
              derf(ic,ir)=derf(ir,ic)
           endif
 2       continue
 1     continue
c ===========================================================
c Computation of variational equations second part
c           ndf=3+1
c           call vetmat(x(ndf),9,b,3,3)
c           call vetmat(x(ndf+9),9,c,3,3)
c           call prodmm(dfb,derf,b)
c           call prodmm(dfc,derf,c)
            do 24 j=1,3
              do 25 k=1,3
                dfb(j,k)=0.d0
                dfc(j,k)=0.d0
                do 26 i=1,3
                  dfb(j,k)=dfb(j,k)+derf(j,i)*x(i+3*k)
                  dfc(j,k)=dfc(j,k)+derf(j,i)*x(i+3*k+9)
 26             continue
 25           continue
 24         continue
c           call matvet(dfb,3,3,vecdfb)
c           call matvet(dfc,3,3,vecdfc)
c           do  4 ind=1,9
c               nind=ndf-1+ind
c               f(nind)=vecdfb(ind)
c               f(nind+9)=vecdfc(ind)
c4          continue
            do 27 j=1,3
              do 28 k=1,3
                f(j+3*k)=dfb(j,k)
                f(j+3*k+9)=dfc(j,k)
 28           continue
 27         continue
      return
      end
c ======================================
c PLANAST
c subroutine providing planets and asteroids 
c in ecliptic coordinates
c at given time t0
c ======================================
      SUBROUTINE planast(t0,ips,imem,istate,xpla)
      IMPLICIT NONE
c input
c   time MJD 
      DOUBLE PRECISION t0      
c   flag for recomputation, memory location, flag for velocities, no. bodies
      INTEGER ips,imem,istate
c output
c asteroid maximum number
      include 'parbep.h'
c planetary masses and other model parameters, including asteroid masses
      include 'masses.h'
      DOUBLE PRECISION xpla(6,nmassx)
c hidden input/output
c flags to control units and center
      double precision pvsun(6)
      logical km,bary
      common/stcomx/km,bary,pvsun
c JPL header
      include 'jplhdr.h'
c asteroid requests to rdbep
c list of asteroid indexes in the asteroid ephemerides arrays
      INCLUDE 'selast.h'
c end interface
c Workspace for state call
      double precision et0(2),rot(3,3),pv(6,12),pnut(4),xx(3),xxp(3)
c stored planets array
      INTEGER memx
      PARAMETER (memx=30)
      DOUBLE PRECISION xp(6,nmassx,memx)
c positions of the massive asteroids
      double precision xast(3,nbepx),vast(3,nbepx)
c loop indexes; ia=1,nast
      INTEGER i,j,k,ia
c initial call flag
      integer lflag
c****************
c static memory required
      SAVE
c****************
      data lflag/0/
c ===========================================================
c reference system rotation matrix
      if(lflag.eq.0)THEN 
         call rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0)
         lflag=1
      ENDIF
c computations to be done only if ips>0
      IF(ips.eq.0)THEN
c Read planetary positions from JPL files
         et0(1)=2400000.5d0
         et0(2)=t0
         bary=.false.
         call state(et0,listpl,pv,pnut,istate)
c reorder data (see pleph)
         DO 9 i=1,npla
            IF(itarg(i).ne.3.and.itarg(i).ne.10.and.itarg(i).ne.13)THEN 
               DO j=1,3
                  if(istate.eq.2)xxp(j)=pv(j+3,itarg(i))
                  xx(j)=pv(j,itarg(i))
               ENDDO
            ELSEIF(itarg(i).eq.13)THEN
               DO j=1,3
                  if(istate.eq.2)xxp(j)=pv(j+3,3)
                  xx(j)=pv(j,3)
               ENDDO
            ELSEIF(itarg(i).eq.10)THEN
               DO j=1,3
                  if(istate.eq.2)
     +                 xxp(j)=pv(j+3,10)*emrat/(1.d0+emrat)+pv(j+3,3)
                  xx(j)=pv(j,10)*emrat/(1.d0+emrat)+pv(j,3)
               ENDDO
            ELSEIF(itarg(i).eq.3)THEN
               DO j=1,3
                  if(istate.eq.2)xxp(j)=pv(j+3,3)-pv(j+3,10)/
     +                 (1.d0+emrat)
                  xx(j)=pv(j,3)-pv(j,10)/(1.d0+emrat)
               ENDDO
            ENDIF
c Change of reference system EQUM00 ---> ECLM00
c        call prodmv(xpla(1,i),rot,xx)
c        if(istate.eq.2)call prodmv(xpla(4,i),rot,xxp)
            DO j=1,3
               xp(j,i,imem)=0.d0
               xp(j+3,i,imem)=0.d0
               DO k=1,3
                  xp(j,i,imem)=xp(j,i,imem)+rot(j,k)*xx(k)
                  IF(istate.eq.2)
     +                 xp(j+3,i,imem)=xp(j+3,i,imem)+rot(j,k)*xxp(k)
               ENDDO
            ENDDO
 9       ENDDO
c ===========================================================
         IF(iatrue.gt.0)THEN
c read asteroid positions from binary ephemerides
            CALL rdbep(t0,iatrue,astid,xast,vast)
c stacking selected asteroids in the planets array
            DO ia=1,iatrue
               DO j=1,3
                  xp(j,npla+ia,imem)=xast(j,ia)
                  IF(istate.eq.2)xp(j+3,npla+ia,imem)=vast(j,ia)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
c ========================================================
c copy into output array
      DO i=1,nmass
         DO j=1,3
            xpla(j,i)=xp(j,i,imem)
            if(istate.eq.2)xpla(j+3,i)=xp(j+3,i,imem)
         ENDDO
      ENDDO
      RETURN
      END

c ***************************************************************
c J2SUN
c ***************************************************************
c This subroutine computes the J_2 acceleration on an asteroid
c due to the Sun; J_2(Sun) = 2 x 10^-7 according to JPL DE405.
c Neglects tilt of the solar spin axis to the ecliptic.
c ***************************************************************
      SUBROUTINE j2sun(x,accj2)
      implicit none
c planetary masses and other model parameters
      double precision solj2,radsun
      parameter (solj2=2.d-7,radsun=4.6527174d-3)
      include 'parbep.h'
      include 'masses.h'
c scalars
      double precision xsun2,xsun5,ratio,brac1,brac2
c vectors
      double precision x(3),accj2(3)
c function
      double precision prscal
c ---------------------------------------------------------------
      xsun2=prscal(x,x)
      xsun5=xsun2*xsun2*dsqrt(xsun2)
      ratio=1.5d0*gm0*solj2*(radsun*radsun/xsun2)/xsun5
      brac1=5.d0*x(3)*x(3)-xsun2
      brac2=5.d0*x(3)*x(3)-3.d0*xsun2
      accj2(1)=ratio*x(1)*brac1
      accj2(2)=ratio*x(2)*brac1
      accj2(3)=ratio*x(3)*brac2
c      write(*,*)accj2(1),accj2(2)
      return
      end
c genrel
c ******************************************************************* c
c this routine treates PN terms in the kinematically                  c
c nonrotating frame according to the scheme given by Damour,          c
c Soffel and Xu (DSX IV, Phys Rev D, Jan 1994).                       c
c                                                                     c
c only the Sun monopole term is accepted 
c output format: 1) drgr(i,1) ... components of the acceleration      c
c                2) drgr(i,j+1) ... derivatives with respect to the   c
c                                 initial conditions                  c
c                                                                     c
c                                    D. Vokrouhlick\'y, 1/3/94        c
c  Modified for ORBFIT, Milani & Baccili 1997                         c
c     computation of derivatives disabled                             c
c ******************************************************************* c
      SUBROUTINE genrel(x,vs,drgr)
      implicit none
c position, velocity, relativistic effects
      double precision x(3),vs(3),drgr(3,7)
c intermediate for partials
c     double precision dmgrx(3,6)
c asteroid maximum number
      include 'parbep.h'
c planetary masses and other model parameters, vlight is required
      include 'masses.h'
      include 'vlight.h'
c scalar temporaries
      double precision rsate2,rsate,xv,v2,c2,eafac,fac,brac
c loop indexes
      integer i
c     integer j
c ---------------------------------------------------------------------
      c2=vlight*vlight
c
      rsate2=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
      rsate=dsqrt(rsate2)
      xv=x(1)*vs(1)+x(2)*vs(2)+x(3)*vs(3)
      v2=vs(1)*vs(1)+vs(2)*vs(2)+vs(3)*vs(3)
c  this refer to the Sun, not the earth
      eafac=gm0/rsate
c
c  Internal terms:
c     --  monopole (DSX IV, 3.13)
c         (rem. also accepted as the IERS Standard)
      brac=4.d0*eafac-v2
      fac=eafac/rsate2/c2
c  acceleration and partial derivatives with respect to the X and V
      do  i=1,3
        drgr(i,1)=fac*(brac*x(i)+4.d0*xv*vs(i))
      enddo
c ---------------------------------------------------------------------
c  partial derivatives available for later use
c     do 21 i=1,3
c       do 22 j=1,3
c         dmgrx(i,j)=fac*(4.d0*vs(i)*(vs(j)-3.d0*xv*x(j)/rsate2)-
c    .       (3.d0*brac+4.d0*eafac)*x(i)*x(j)/rsate2)
c         dmgrx(i,j+3)=2.d0*fac*(2.d0*vs(i)*x(j)-x(i)*vs(j))
c         if (i.eq.j) then 
c            dmgrx(i,j)=dmgrx(i,j)+brac*fac
c            dmgrx(i,j+3)=dmgrx(i,j+3)+4.d0*xv*fac
c         endif
c22     enddo
c21   enddo
c
c     do 31 i=1,3 
c       drgr(i,2)=dmgrx(i,1)
c       drgr(i,3)=dmgrx(i,2)
c       drgr(i,4)=dmgrx(i,3)
c       do 30 j=4,6
c         drgr(i,j+1)=dmgrx(i,j)
c30     enddo
c31   enddo
c ---------------------------------------------------------------------
c
      return
      end
c ***************************************************************
c EIHREL
c ***************************************************************
c This subroutine computes heliocentric acceleration of the
c asteroid on the (1/c^2)(= PN) level. Both the direct solar term
c (~Schwarzschild contribution) and the planetary contributions
c are included. The monopole terms considered only (= EIH
c approximation).
c
c Written by D Vokrouhlicky and S Chesley, Nov 3, 1999
c ***************************************************************
      SUBROUTINE eihrel(x,vs,xpla,d,rpla,xsun,drgr)
      implicit none
c planetary masses and other model parameters, vlight is required
      include 'parbep.h'
      include 'masses.h'
      include 'vlight.h'
c position, velocity, relativistic effects
      double precision x(3),vs(3),xpla(6,nmassx),d(nmassx),drgr(3)
      double precision rpla(nmassx)
c scalar temporaries
      double precision temp,rplapla,strqua,potall,potalls,gmall,tempp
      double precision xsun,xsun3,xdsun2,vast2,vsun2,rplapla2,scals
      double precision vlight2
c loop indexes
      integer i,ii,j,ib
      double precision g(3),vsun(3),vast(3),vpla(3,nmassx),tbtvec(3)
      double precision xastpla(3),nastpla(3),vastpla(3)
      double precision rplavec(3),tmpvec(3),rpla3(nmassx)
c ---------------------------------------------------------------------
c initialization
      do i=1,3
         drgr(i)=0.d0
      enddo
c compute potential of all bodies on asteroid, planets on Sun
c and heliocentric distances of planets and the total mass of the
c system
      potall=gm0/xsun
      gmall=gm0
      potalls=0.d0
      do i=1,npla
         potall=potall+gm(i)/d(i)
         rpla3(i)=rpla(i)**3
         potalls=potalls+gm(i)/rpla(i)
         gmall=gmall+gm(i)
      enddo
c compute barycentric velocities of sun, asteroid, and planets
      do j=1,3
         vsun(j)=0.d0
         do i=1,npla
            vsun(j)=vsun(j)-gm(i)*xpla(j+3,i)/gmall
         enddo
         vast(j)=vsun(j)+vs(j)
      enddo
      vsun2=vsun(1)*vsun(1)+vsun(2)*vsun(2)+vsun(3)*vsun(3)
      vast2=vast(1)*vast(1)+vast(2)*vast(2)+vast(3)*vast(3)
      do i=1,npla
         do j=1,3
            vpla(j,i)=vsun(j)+xpla(3+j,i)
         enddo
      enddo
c ................................................................
c compute indirect (1/c^2) planetary accelerations:
      do ib=1,npla
         do j=1,3
            g(j)=-gm(ib)*xpla(j,ib)/rpla3(ib)
         enddo
c        strange quantity
         strqua=1.5d0*gm0/rpla(ib)
         do ii=1,npla
            if(ii.ne.ib) then
               do j=1,3
                  rplavec(j)=xpla(j,ii)-xpla(j,ib)
               enddo
               rplapla2=rplavec(1)**2+rplavec(2)**2+rplavec(3)**2
               rplapla=dsqrt(rplapla2)
               scals=xpla(1,ib)*rplavec(1)+xpla(2,ib)*rplavec(2)+
     +               xpla(3,ib)*rplavec(3)
               strqua=strqua+gm(ii)*(1.d0-0.5d0*scals/rplapla2)/rplapla
            endif
         enddo
         tempp=(xpla(1,ib)*vpla(1,ib)+xpla(2,ib)*vpla(2,ib)+
     +          xpla(3,ib)*vpla(3,ib))/rpla(ib)
         temp=2.d0*(xpla(4,ib)*xpla(4,ib)+xpla(5,ib)*xpla(5,ib)+
     +        xpla(6,ib)*xpla(6,ib))-vsun2-1.5d0*tempp*tempp-
     +        4.d0*potalls-strqua
c add "g*temp" terms
         do i=1,3
            drgr(i)=drgr(i)+temp*g(i)
         enddo
c add "velocity-dependent" terms
         do i=1,3
            tmpvec(i)=-4.d0*xpla(3+i,ib)+vpla(i,ib)
         enddo
         temp=g(1)*tmpvec(1)+g(2)*tmpvec(2)+g(3)*tmpvec(3)
         do i=1,3
            drgr(i)=drgr(i)+temp*xpla(3+i,ib)
         enddo
c "third-body" term (tbt)
         do i=1,3
            tbtvec(i)=gm0*xpla(i,ib)/rpla3(ib)
         enddo
         do ii=1,npla
            if(ii.ne.ib) then
               do j=1,3
                  rplavec(j)=xpla(j,ib)-xpla(j,ii)
               enddo
               rplapla=dsqrt(rplavec(1)**2+rplavec(2)**2+rplavec(3)**2)
               temp=gm(ii)/(rplapla**3)
               do i=1,3
                  tbtvec(i)=tbtvec(i)+temp*rplavec(i)
               enddo
            endif
         enddo
         temp=3.5d0*gm(ib)/rpla(ib)
         do i=1,3
            drgr(i)=drgr(i)+temp*tbtvec(i)
         enddo
      enddo
c compute direct (1/c^2) accelerations:
c -- planetary terms
      do ib=1,npla
         do j=1,3
            g(j)=gm(ib)*(x(j)-xpla(j,ib))/(d(ib)**3)
         enddo
c        pos rel to planet (unit vector)
         do j=1,3
            xastpla(j)=x(j)-xpla(j,ib)
            nastpla(j)=xastpla(j)/d(ib)
         enddo
c        vast rel to planet
         do j=1,3
            vastpla(j)=vast(j)-vpla(j,ib)
         enddo
c        strange quantity
         scals=xastpla(1)*xpla(1,ib)+xastpla(2)*xpla(2,ib)+
     +         xastpla(3)*xpla(3,ib)
         strqua=gm0*(1.d0-0.5d0*scals/(rpla(ib)**2))/rpla(ib)
         do ii=1,npla
            if(ii.ne.ib)then
               do j=1,3
                  rplavec(j)=xpla(j,ii)-xpla(j,ib)
               enddo
               rplapla2=rplavec(1)**2+rplavec(2)**2+rplavec(3)**2
               rplapla=dsqrt(rplapla2)
               scals=xastpla(1)*rplavec(1)+xastpla(2)*rplavec(2)+
     +               xastpla(3)*rplavec(3)
               strqua=strqua+gm(ii)/rplapla*(1.d0+0.5d0*scals/rplapla2)
            endif
         enddo
         tempp=nastpla(1)*vpla(1,ib)+nastpla(2)*vpla(2,ib)+
     +         nastpla(3)*vpla(3,ib)
         temp=2.d0*(vastpla(1)**2+vastpla(2)**2+vastpla(3)**2)-
     +        vast2-1.5d0*tempp*tempp-4.d0*potall-strqua
c        add g*temp to force
         do i=1,3
            drgr(i)=drgr(i)-g(i)*temp
         enddo
c add "velocity-dependent" terms
         do i=1,3
            tmpvec(i)=4.d0*vastpla(i)+vpla(i,ib)
         enddo
         temp=g(1)*tmpvec(1)+g(2)*tmpvec(2)+g(3)*tmpvec(3)
         do i=1,3
            drgr(i)=drgr(i)+vastpla(i)*temp
         enddo
c "third-body" term (tbt)
         do i=1,3
            tbtvec(i)=gm0*xpla(i,ib)/rpla3(ib)
         enddo
         do ii=1,npla
            if(ii.ne.ib)then
               do j=1,3
                  rplavec(j)=xpla(j,ib)-xpla(j,ii)
               enddo
               rplapla=dsqrt(rplavec(1)**2+rplavec(2)**2+rplavec(3)**2)
               temp=gm(ii)/(rplapla**3)
               do i=1,3
                  tbtvec(i)=tbtvec(i)+temp*rplavec(i)
               enddo
            endif
         enddo
         temp=-3.5d0*gm(ib)/d(ib)
         do i=1,3
            drgr(i)=drgr(i)+temp*tbtvec(i)
         enddo
      enddo
c compute direct (1/c^2) accelerations:
c -- solar (~ Schwarzschild) term
      xsun3=xsun**3
      do i=1,3
         g(i)=gm0*x(i)/xsun3
      enddo
      xdsun2=vs(1)*vs(1)+vs(2)*vs(2)+vs(3)*vs(3)
c     strange quantity
      strqua=0.d0
      do i=1,npla
         scals=x(1)*xpla(1,i)+x(2)*xpla(2,i)+x(3)*xpla(3,i)
         strqua=strqua+gm(i)*(1.d0+0.5d0*scals/(rpla(i)**2))/rpla(i)
      enddo
      tempp=(x(1)*vsun(1)+x(2)*vsun(2)+x(3)*vsun(3))/xsun
      temp=2.d0*xdsun2-vast2-1.5d0*tempp*tempp-4.d0*potall-strqua
c add g*temp to force
      do i=1,3
        drgr(i)=drgr(i)-g(i)*temp
      enddo
c add velocity terms ...
      do i=1,3
         tmpvec(i)=4.d0*vs(i)+vsun(i)
      enddo
      temp=g(1)*tmpvec(1)+g(2)*tmpvec(2)+g(3)*tmpvec(3)
      do i=1,3
         drgr(i)=drgr(i)+temp*vs(i)
      enddo
c "third-body" term (tbt)
      do i=1,3
         tbtvec(i)=0.d0
      enddo
      do ii=1,npla
         temp=gm(ii)/rpla3(ii)
         do i=1,3
            tbtvec(i)=tbtvec(i)+temp*xpla(i,ii)
         enddo
      enddo
      temp=3.5d0*gm0/xsun
      do i=1,3
         drgr(i)=drgr(i)+temp*tbtvec(i)
      enddo
c ..................................................................
c scale by 1/c^2
      vlight2=vlight**2
      do i=1,3
         drgr(i)=drgr(i)/vlight2
      enddo
      return
      end
c ******************************************************************
      SUBROUTINE yarkdi(xast,a,iparti)
c ******************************************************************
c
c This subroutine computes the heliocentric components of the
c Yarkovsky thermal acceleration -- the diurnal variant only.
c If the flag (iparti) in the common block is set to 1, one gets at
c the output also partials wrt to some parameters of the thermal
c model (if these are desired to be adjusted).
c
c Input parameters:
c -----------------
c
c - via header   xast(3) ... heliocentric coordinates of the body (in AU)
c - via common   yarkp(1-3) ... sx, sy, sz (unit vector of the
c                               body's spin axis orientation)
c                yarkp(4-5) ... k_0 and k_1 parameters of the
c                               surface thermal conductivity
c                               [K(T) = k_0 + k_1 T_av^3]
c                yarkp(6) ... density of the surface layer
c                yarkp(7) ... radius of the body
c                yarkp(8) ... rotation frequency
c                yarkp(9) ... surface absorptivity
c                iparti   ... partials (yes=1/no=0)
c                             [presently only partials listed below,
c                              a(4) - a(21) are available]
c
c
c Output parameters: a(1-3) ... diurnal acceleration
c ------------------ a(4-6) ... partials wrt the radius of the body
c                    a(7-9) ... partials wrt the thermal conductivity
c                               parameter k_0
c                    a(10-12) ... partials wrt the thermal conductivity
c                                 parameter k_1
c                    a(13-15) ... partials wrt the x-component of the
c                                 spin axis unit vector
c                    a(16-18) ... partials wrt the y-component of the
c                                 spin axis unit vector
c                    a(19-21) ... partials wrt the z-component of the
c                                 spin axis unit vector
c
c SI units are assumed internally in the subroutine, but the results
c (e.g. accelerations) are given in AU and days.
c
c Written by: D. Vokrouhlicky, Oct 99
c (queries to vokrouhl@mbox.cesnet.cz)
c ..................................................................
      implicit double precision (a-h,o-z)
c here we specify two more parameters that eventually might be changed:
c -- the average bulk density of the body (densityb) which is now set
c    to 2 g/cm^3
c -- the heat capacity of the surface layer (capacity) which is now set
c    to 680 J/kg/K
      parameter (densityb=2.d3,capacity=680.d0,solcon=1371.d0)
      parameter (emiss=0.9d0,stefboltz=5.66962d-8,clight3=8.99377374d8)
      parameter (dsqrt2=1.414213562373d0,dsqrt23=1414.213562373d0)
      parameter (aceuni=0.049900176d0)
c input: asteroid position, flag for partials 
      double precision xast(3)
      integer iparti
c output: acceleration and partials
      double precision a(21)
c internal variables
      double precision vprod1(3),vprod2(3)
c physical data on the current asteroid
      INCLUDE 'yarkov.h'    
c -----------------------------------------------------------------------
      rau2=xast(1)*xast(1)+xast(2)*xast(2)+xast(3)*xast(3)
      rau=dsqrt(rau2)
      xn=xast(1)/rau
      yn=xast(2)/rau
      zn=xast(3)/rau
c initializations & constants
      radflu=solcon/rau2
c - subsolar temperature
      tstar=(yarkp(9)*radflu/emiss/stefboltz)**0.25d0
      tav1000=tstar/dsqrt23
c - surface conductivity
      surcon=yarkp(4)+yarkp(5)*(tav1000**3)
c - thermal inertia & diurnal thermal parameter
      bgama=dsqrt(surcon*yarkp(6)*capacity)
      theta=bgama*dsqrt(yarkp(8))/emiss/stefboltz/(tstar**3)
      diudepth=dsqrt(surcon/yarkp(6)/capacity/yarkp(8))
c - radius of the body scaled by the depth of the diurnal wave
      rp=yarkp(7)/diudepth
      al=dsqrt2*rp
      tau=theta/al
      tau1=1.d0+tau
c - the auxiliary functions A-D, a,b
      cal=dcos(al)
      sal=dsin(al)
      if (al.lt.90.d0) then
       ealm=dexp(-al)
      else
       ealm=0.d0
      endif
      af=3.d0*(al+2.d0)*ealm+(3.d0*(al-2.d0)*cal+al*(al-3.d0)*sal)
      bf=al*(al+3.d0)*ealm+(-al*(al-3.d0)*cal+3.d0*(al-2.d0)*sal)
      caf=-(al+2.d0)*ealm+(-(al-2.d0)*cal+al*sal)
      cbf=-al*ealm-(al*cal+(al-2.d0)*sal)
      ccf=caf+tau*af/tau1
      cdf=cbf+tau*bf/tau1
c - G exp(i delta) & amplitude computed
      facp=aceuni*yarkp(9)*radflu/yarkp(7)/densityb/clight3
      deno=ccf*ccf+cdf*cdf
      deno1=deno*tau1
      gcosd=(caf*ccf+cbf*cdf)/deno1
      gsind=(cbf*ccf-caf*cdf)/deno1
c geometric products
c - r x s
      vprod1(1)=yn*yarkp(3)-zn*yarkp(2)
      vprod1(2)=zn*yarkp(1)-xn*yarkp(3)
      vprod1(3)=xn*yarkp(2)-yn*yarkp(1)
c - s x (r x s) = r - (r.s) s
      scalar=xn*yarkp(1)+yn*yarkp(2)+zn*yarkp(3)
      vprod2(1)=xn-scalar*yarkp(1)
      vprod2(2)=yn-scalar*yarkp(2)
      vprod2(3)=zn-scalar*yarkp(3)
c diurnal acceleration
      a(1)=facp*(gsind*vprod1(1)+gcosd*vprod2(1))
      a(2)=facp*(gsind*vprod1(2)+gcosd*vprod2(2))
      a(3)=facp*(gsind*vprod1(3)+gcosd*vprod2(3))
c Partials?
      if (iparti.eq.0) return
c - general
      cafp=-ealm+cal+(2.d0*al-1.d0)*sal
      cbfp=-ealm-(2.d0*al-1.d0)*cal+sal
      afp=3.d0*ealm+(al*al-3.d0)*cal+(al*(al-4.d0)+3.d0)*sal
      bfp=(2.d0*al+3.d0)*ealm-(al*(al-4.d0)+3.d0)*cal
     .     +(al*al-3.d0)*sal
c - thermal conductivity parameters (k_0,k_1)
      xi1r=caf*ccf-cbf*cdf
      xi1i=cbf*ccf+caf*cdf
      xi2r=cafp*af-cbfp*bf
      xi2i=cbfp*af+cafp*bf
      xi2r=xi2r-caf*afp+cbf*bfp
      xi2i=xi2i-cbf*afp-caf*bfp
      deno=xi1r*xi1r+xi1i*xi1i
      facr=1.d0+0.5d0*al*(xi2r*xi1r+xi2i*xi1i)/deno
      faci=     0.5d0*al*(xi2i*xi1r-xi2r*xi1i)/deno
      derikr=-tau*(gcosd*facr-gsind*faci)/tau1
      deriki=-tau*(gsind*facr+gcosd*faci)/tau1
      a(7)=facp*(deriki*vprod1(1)+derikr*vprod2(1))
      a(8)=facp*(deriki*vprod1(2)+derikr*vprod2(2))
      a(9)=facp*(deriki*vprod1(3)+derikr*vprod2(3))
      a(10)=a(7)*(tav1000**3)
      a(11)=a(8)*(tav1000**3)
      a(12)=a(9)*(tav1000**3)
c - radius of the body
      rfac=(tau+tau1)/tau1
      a(4)=-a(1)*rfac-2.d0*a(7)
      a(5)=-a(2)*rfac-2.d0*a(8)
      a(6)=-a(3)*rfac-2.d0*a(9)
c - partials d_K (a), d_R (a) ...
      a(4)=a(4)/yarkp(7)
      a(5)=a(5)/yarkp(7)
      a(6)=a(6)/yarkp(7)
      a(7)=a(7)/surcon
      a(8)=a(8)/surcon
      a(9)=a(9)/surcon
      a(10)=a(10)/surcon
      a(11)=a(11)/surcon
      a(12)=a(12)/surcon
c - spin axis components
c ... sx
      a(13)=-facp*gcosd*(xn*yarkp(1)+scalar)
      a(14)=facp*(gsind*zn-gcosd*xn*yarkp(2))
      a(15)=-facp*(gsind*yn+gcosd*xn*yarkp(3))
c ... sy
      a(16)=-facp*(gsind*zn+gcosd*yn*yarkp(1))
      a(17)=-facp*gcosd*(yn*yarkp(2)+scalar)
      a(18)=facp*(gsind*xn-gcosd*yn*yarkp(3))
c ... sz
      a(19)=facp*(gsind*yn-gcosd*zn*yarkp(1))
      a(20)=-facp*(gsind*xn+gcosd*zn*yarkp(2))
      a(21)=-facp*gcosd*(zn*yarkp(3)+scalar)
      return
      end
c ******************************************************************
      SUBROUTINE yarkse(xast,vast,a,iparti)
c ******************************************************************
c
c This subroutine computes the heliocentric components of the
c Yarkovsky thermal acceleration -- the seasonal variant only.
c If the flag (iparti) is set to 1, one gets at the output also
c partials wrt to some parameters of the thermal model (if these
c are desired to be adjusted).
c
c Input parameters:
c -----------------
c
c - via header   iparti  ... partials (yes=1/no=0)
c                (xast,vast) ... state vector of the asteroid
c - via common   yarkp(1-3) ... sx, sy, sz (unit vector of the
c                               body's spin axis orientation)
c                yarkp(4-5) ... k_0 and k_1 parameters of the
c                               surface thermal conductivity
c                               [K(T) = k_0 + k_1 T_av^3]
c                yarkp(6) ... density of the surface layer
c                yarkp(7) ... radius of the body
c                yarkp(8) ... rotation frequency
c                yarkp(9) ...  surface absorptivity
c                + some more precomputed useful variables
c
c Output parameters: a(1-3) ... seasonal acceleration
c ------------------ a(4-6) ... partials wrt the radius of the body
c                    a(7-9) ... partials wrt the thermal conductivity
c
c REM. PARTIALS ARE DISABLED AT THIS MOMENT
c
c SI units are assumed throughout the subroutine, but the results
c (e.g. accelerations) are given in AU and days.
c
c Written by: D. Vokrouhlicky, Oct 99
c (queries to vokrouhl@mbox.cesnet.cz)
c ..................................................................
      implicit double precision (a-h,o-z)
      parameter (napprox=7)
      parameter (densityb=2.d3,capacity=680.d0,dsqrt2=1.414213562373d0)
      parameter (emiss=0.9d0,clight3=8.99377374d8,aceuni=0.049900176d0)
      dimension xast(3),vast(3)
      dimension brac(7),bras(7),gcosd(7),gsind(7),a(21)
      INCLUDE 'sunmass.h'
      INCLUDE 'trig.h'
      INCLUDE 'yarkov.h'
c -----------------------------------------------------------------------
c - thermal inertia & seasonal thermal parameter
       bgama=dsqrt(yarkp(4)*yarkp(6)*capacity)
       theta=bgama*thfacya/emiss
       seadepth=dsqrt(yarkp(4)/yarkp(6)/capacity/fmeaya)
c - radius of the body scaled by the depth of the seasonal wave
       rp=yarkp(7)/seadepth
       rp2=dsqrt2*rp
       tau=theta*etaya75/rp2
       tau1=1.d0+tau
c - amplitude of the effect
       fac=aceuni*yarkp(9)*radfluya/yarkp(7)/densityb/clight3/tau1
c - G_k cos(d_k) & G_K sin(d_k) functions computed
       do 10 k=1,napprox
        fk=dfloat(k)
        alk=dsqrt(fk)*rp2
c - the auxiliary functions A-D, a,b
        cal=dcos(alk)
        sal=dsin(alk)
        if (alk.lt.90.d0) then
         ealm=dexp(-alk)
        else
         ealm=0.d0
        endif
        af=3.d0*(alk+2.d0)*ealm+(3.d0*(alk-2.d0)*cal+alk*(alk-3.d0)*sal)
        bf=alk*(alk+3.d0)*ealm+(-alk*(alk-3.d0)*cal+3.d0*(alk-2.d0)*sal)
        caf=-(alk+2.d0)*ealm+(-(alk-2.d0)*cal+alk*sal)
        cbf=-alk*ealm-(alk*cal+(alk-2.d0)*sal)
        ccf=caf+tau*af/tau1
        cdf=cbf+tau*bf/tau1
c - G exp(i delta)
        deno=ccf*ccf+cdf*cdf
        gcosd(k)=(caf*ccf+cbf*cdf)/deno
        gsind(k)=(cbf*ccf-caf*cdf)/deno
c compute cos- & sin-related brackets
        brac(k)=spya*alya(k)*gcosd(k)+sqya*beya(k)*gsind(k)
        bras(k)=sqya*beya(k)*gcosd(k)-spya*alya(k)*gsind(k)
10     continue 
c mean anomaly detremined
c - approximated by a linear term only
c      anomaly=ele0(6)+(fmea*t)
c - computed from the state vector
      r2=xast(1)*xast(1)+xast(2)*xast(2)+xast(3)*xast(3)
      v2=vast(1)*vast(1)+vast(2)*vast(2)+vast(3)*vast(3)
      rdot=xast(1)*vast(1)+xast(2)*vast(2)+xast(3)*vast(3)
      r=dsqrt(r2)
      aaxi=1.d0/(2.d0/r-(v2/gms))
      esinu=rdot/dsqrt(aaxi*gms)
      ecosu=(r*v2/gms)-1.d0
      uano=datan2(esinu,ecosu)
      anomaly=uano-esinu
      if (anomaly.lt.0.d0) anomaly=anomaly+dpig
c compute the sum...
      fact=0.d0
      do 100 k=napprox,1,-1
       fk=dfloat(k)
       canomaly=dcos(fk*anomaly)
       sanomaly=dsin(fk*anomaly)
       fact=fact+(brac(k)*canomaly+bras(k)*sanomaly)
100   continue
      fact=fact*fac
c seasonal acceleration (~ factor * {\bf s})
      a(1)=fact*yarkp(1)
      a(2)=fact*yarkp(2)
      a(3)=fact*yarkp(3)
c Partials? -- DISABLED AT THE MOMENT
c      if (iparti.eq.0) return
c - general
c      cafp=-ealm+cal+(2.d0*al-1.d0)*sal
c      cbfp=-ealm-(2.d0*al-1.d0)*cal+sal
c      afp=3.d0*ealm+(al*al-3.d0)*cal+(al*(al-4.d0)+3.d0)*sal
c      bfp=(2.d0*al+3.d0)*ealm-(al*(al-4.d0)+3.d0)*cal
c     .     +(al*al-3.d0)*sal
c - thermal conductivity parameters (k_0,k_1)
c      xi1r=caf*ccf-cbf*cdf
c      xi1i=cbf*ccf+caf*cdf
c      xi2r=cafp*af-cbfp*bf
c      xi2i=cbfp*af+cafp*bf
c      xi2r=xi2r-caf*afp+cbf*bfp
c      xi2i=xi2i-cbf*afp-caf*bfp
c      deno=xi1r*xi1r+xi1i*xi1i
c      facr=1.d0+0.5d0*al*(xi2r*xi1r+xi2i*xi1i)/deno
c      faci=     0.5d0*al*(xi2i*xi1r-xi2r*xi1i)/deno
c      derikr=-tau*(gcosd*facr-gsind*faci)/tau1
c      deriki=-tau*(gsind*facr+gcosd*faci)/tau1
c      a(7)=fac*(deriki*vprod1(1)+derikr*vprod2(1))
c      a(8)=fac*(deriki*vprod1(2)+derikr*vprod2(2))
c      a(9)=fac*(deriki*vprod1(3)+derikr*vprod2(3))
c      a(10)=a(7)*(tav1000**3)
c      a(11)=a(8)*(tav1000**3)
c      a(12)=a(9)*(tav1000**3)
c - radius of the body
c      rfac=(tau+tau1)/tau1
c      a(4)=-a(1)*rfac-2.d0*a(7)
c      a(5)=-a(2)*rfac-2.d0*a(8)
c      a(6)=-a(3)*rfac-2.d0*a(9)
c - partials d_K (a), d_R (a) ...
c      a(4)=a(4)/yarkp(7)
c      a(5)=a(5)/yarkp(7)
c      a(6)=a(6)/yarkp(7)
c      a(7)=a(7)/surcon!!!!!! --> yarkp(4)
c      a(8)=a(8)/surcon
c      a(9)=a(9)/surcon
c      a(10)=a(10)/surcon
c      a(11)=a(11)/surcon
c      a(12)=a(12)/surcon
c - spin axis components
      return
      end
c ==================================================================
c yarkinit: initialisation of the yarkovsky force model for a given asteroid
c           written by A. Milani & D. Vokrouhlicky, Oct 99
      SUBROUTINE yarkinit(astnam,eltype,elem)
      IMPLICIT NONE
      CHARACTER*(*) astnam
      CHARACTER*3 eltype
      DOUBLE PRECISION elem(6)
      CHARACTER*80 file
      DOUBLE PRECISION lat,long,emiss,stefboltz,argu,argu2,tstarya,eta
      DOUBLE PRECISION elkep(6),pvya(3),qvya(3),nvya(3),enne,cgam,obli
      INTEGER unit,le
      INCLUDE 'yarkov.h'  
      INCLUDE 'yarkom.h'
      INCLUDE 'trig.h'
      INCLUDE 'model.h'
      INCLUDE 'sunmass.h'
      INCLUDE 'sysdep.h'
c yar is the logical flag for the existence of the physical data 
c allowing computation of Yarkovsky; otherwise, the non gravitational
c force is set to zero 
c
      IF(iyark.eq.0)RETURN
      yarini=.true.
c convert elements to keplerian
      call coocha(elem,eltype,gms,elkep,'KEP',enne)
c compute the name of the file which could contain the yarkovsky data
      CALL filnam(yardir,astnam,'yar',file,le)
      INQUIRE(file=file(1:le),exist=yarfil)
      IF(yarfil)THEN
         call filopn(unit,file(1:le),'old')
         read(unit,*,end=111)
c ecliptic latitude and longitude of the spin axis
         read(unit,*,end=111)long
         read(unit,*,end=111)lat
c - via common   yarkp(1-3) ... sx, sy, sz (unit vector of the
c                               body's spin axis orientation)
c                yarkp(4-5) ... k_0 and k_1 parameters of the
c                               surface thermal conductivity
c                               [K(T) = k_0 + k_1 T_av^3]
c                yarkp(6) ... density of the surface layer
c                yarkp(7) ... radius of the body
c                yarkp(8) ... rotation frequency
c                yarkp(9) ... surface absorptivity
         yarkp(1)=dcos(lat*radeg)*dcos(long*radeg)
         yarkp(2)=dcos(lat*radeg)*dsin(long*radeg)
         yarkp(3)=dsin(lat*radeg)
         read(unit,*,end=111)yarkp(4)
         read(unit,*,end=111)yarkp(5)
         read(unit,*,end=111)yarkp(6)
         read(unit,*,end=111)yarkp(7)
         read(unit,*,end=111)yarkp(8)
         read(unit,*,end=111)yarkp(9)
c precompute some variables for the seasonal variant of the Yarkovsky
c effect:
c - constants
         emiss=0.9d0
         stefboltz=5.66962d-8
c - mean motion & solar radiation flux at r=a
         fmeaya=(1.9909837d-7)/elkep(1)/dsqrt(elkep(1))
         radfluya=1371.d0/elkep(1)/elkep(1)
c - subsolar temperature
         tstarya=(yarkp(9)*radfluya/emiss/stefboltz)**0.25d0
         thfacya=dsqrt(fmeaya)/stefboltz/(tstarya**3)
c - projections s_P and s_Q of the spin axis computed
         pvya(1)=dcos(elkep(4))*dcos(elkep(5))-
     .           dcos(elkep(3))*dsin(elkep(4))*dsin(elkep(5))
         pvya(2)=dsin(elkep(4))*dcos(elkep(5))+
     .           dcos(elkep(3))*dcos(elkep(4))*dsin(elkep(5))
         pvya(3)=dsin(elkep(3))*dsin(elkep(5))
         qvya(1)=-dcos(elkep(4))*dsin(elkep(5))-
     .           dcos(elkep(3))*dsin(elkep(4))*dcos(elkep(5))
         qvya(2)=-dsin(elkep(4))*dsin(elkep(5))+
     .         dcos(elkep(3))*dcos(elkep(4))*dcos(elkep(5))
         qvya(3)=dsin(elkep(3))*dcos(elkep(5))
         nvya(1)=dsin(elkep(3))*dsin(elkep(4))
         nvya(2)=-dsin(elkep(3))*dcos(elkep(4))
         nvya(3)=dcos(elkep(3))
         spya=yarkp(1)*pvya(1)+yarkp(2)*pvya(2)+yarkp(3)*pvya(3)
         sqya=yarkp(1)*qvya(1)+yarkp(2)*qvya(2)+yarkp(3)*qvya(3)
         cgam=yarkp(1)*nvya(1)+yarkp(2)*nvya(2)+yarkp(3)*nvya(3)
         obli=dacos(cgam)/radeg
         write(*,*)' Obliquity of the spin axis; Yarkovsky: ',obli
c - compute the \alpha(k) and \beta(k) coefficients
         eta=dsqrt(1.d0-elkep(2)*elkep(2))
         etaya75=eta**0.75d0
c -- \beta_1(x) ... \beta_7(x) functions
         argu=elkep(2)
         argu2=argu*argu
         beya(1)=eta*(1.d0+argu2*(-1152.d0+argu2*(48.d0-argu2))/9216.d0)
         argu=2.d0*elkep(2)
         argu2=argu*argu
         beya(2)=eta*argu*(1.d0+argu2*(-1920.d0+argu2*(60.d0-argu2))
     .           /23040.d0)
         argu=3.d0*elkep(2)
         argu2=argu*argu
         beya(3)=3.d0*eta*argu2*(1.d0+argu2*(-40.d0+argu2)/640.d0)/8.d0
         argu=4.d0*elkep(2)
         argu2=argu*argu
         beya(4)=eta*argu2*argu*(1.d0+argu2*(-48.d0+argu2)/960.d0)/12.d0
         argu=5.d0*elkep(2)
         argu2=argu*argu
         beya(5)=5.d0*eta*argu2*argu2*(1.d0-argu2/24.d0)/384.d0
         argu=6.d0*elkep(2)
         argu2=argu*argu
         beya(6)=eta*argu2*argu2*argu*(1.d0-argu2/28.d0)/640.d0
         argu=7.d0*elkep(2)
         argu2=argu*argu
         beya(7)=7.d0*eta*argu2*argu2*argu2/46080.d0
c -- \alpha_1(x) ... \alpha_7(x) functions
         argu=elkep(2)
         argu2=argu*argu
         alya(1)=1.d0+argu2*(-3456.d0+argu2*(240.d0-7.d0*argu2))/9216.d0
         argu=2.d0*elkep(2)
         argu2=argu*argu
         alya(2)=argu*(1.d0+argu2*(-960.d0+argu2*(45.d0-argu2))/5760.d0)
         argu=3.d0*elkep(2)
         argu2=argu*argu
         alya(3)=3.d0*argu2*(1.d0+argu2*(-200.d0+7.d0*argu2)/1920.d0)/
     .           8.d0
         argu=4.d0*elkep(2)
         argu2=argu*argu
         alya(4)=argu*argu2*(1.d0+argu2*(-36.d0+argu2)/480.d0)/12.d0
         argu=5.d0*elkep(2)
         argu2=argu*argu
         alya(5)=argu2*argu2*(1.d0-7.d0*argu2/120.d0)/76.8d0
         argu=6.d0*elkep(2)
         argu2=argu*argu
         alya(6)=argu*argu2*argu2*(1.d0-argu2/21.d0)/640.d0
         argu=7.d0*elkep(2)
         argu2=argu*argu
         alya(7)=7.d0*argu2*argu2*argu2/46080.d0
c close the input file
         call filclo(unit,' ')
      ELSE
         WRITE(*,*)' Yarkovsky datafile not found:',file(1:le)
         stop     
      ENDIF
      WRITE(*,*)' Yarkovsky data loaded for asteroid ', astnam
      RETURN
 111  yarfil=.false.
      WRITE(*,*)' incomplete yarkovsky file for asteroid ', astnam
      RETURN
      END
c
      SUBROUTINE selpert(name,found)
      IMPLICIT NONE
      LOGICAL found
      CHARACTER*(*) name
      CHARACTER*9 nam1
      CHARACTER*30 string
      INCLUDE 'parbep.h'
      INCLUDE 'combep.h'
      INCLUDE 'selast.h'
      INCLUDE 'bifina.h'
c controls of the force model
      include 'model.h'
c
      INTEGER iabe,ln,ls,ia,iat
c
      found=.false.
      IF(iast.eq.0)RETURN
c
      nam1=name
      CALL rmsp(nam1,ln)
      call filopn(iabe,filbec,'old')
c
      read(iabe,*)
      read(iabe,*)
      iat=0
      do  ia=1,iast
         read(iabe,201,err=202,end=202)masbep(ia),string
         call rmsp(string,ls)
         IF(ls.eq.ln.and.nam1(1:ln).eq.string(1:ls))THEN
            WRITE(*,*)' self perturbation of ',nam1(1:ln),' avoided'         
            found=.true.
         ELSE
            iat=iat+1
            astid(iat)=ia
         ENDIF
      enddo
      iatrue=iat
      goto 203
 201  FORMAT(1P,E18.10,1X,A)
 202  WRITE(*,*)'selpert: too many asteroids requested, iast=',iast
      iast=ia-1
      WRITE(*,*)'selpert: asteroids available ',iast
 203  call filclo(iabe,' ')
      RETURN
      END
      SUBROUTINE selpert2(nam0,namp,nfound)
      IMPLICIT NONE 
      CHARACTER*(*) nam0,namp
      LOGICAL found0,foundp
      INTEGER nfound
      INCLUDE 'model.h'
      CALL selpert(nam0,found0)
      CALL selpert(namp,foundp)
      IF(found0.and.foundp)THEN
         WRITE(*,*)' please do not try to identify ',nam0,' with ', namp
         WRITE(*,*)' All perturbations by massive asteroids disabled'
         nfound=2
         iast=0
      ELSEIF(found0)THEN
         WRITE(*,*)' you should not do identification with an'
         WRITE(*,*)' asteroid with mass, such as ',nam0
         WRITE(*,*)' perturbations by ',nam0,' disabled'
         nfound=1
      ELSEIF(foundp)THEN
         WRITE(*,*)' you should not do identification with an'
         WRITE(*,*)' asteroid with mass, such as ',namp
         WRITE(*,*)' perturbations by ',namp,' disabled'
         nfound=1
      ELSE
         nfound=0
      ENDIF
      RETURN
      END

