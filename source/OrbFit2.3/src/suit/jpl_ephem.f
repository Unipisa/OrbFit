c ================MODULE jpl_ephem.f======================
c JPL AND BINARY EPHEMERIS Standish & Newhall, but fixed by Orbfit
c version for asteroids by Carpino
c PUBLIC ROUTINES 
c earth  orbital elelemts 
c earcar cartesian 
c rdbep		read and interpolate binary ephemeris files 
c masjpl        setup of masses, list of planets, closapp distances
c trange        max and min time in ephem, sets some constants
c jpllis	get JPL masses and IDs for a list of planets
c
c MODULE CONTAINS
c SUBROUTINES
c state           access to JPL ephem (clumsy...)
c fszer2          internal routines for accessing JPL DE ephemerides
c dpleph          id.
c interp          id.
c const           id.
c split           id.
c
c HEADERS
c    jplhdr.h
c    timespan.h
c    vlight.h
c
c  note: consider directory ../jpl for insertion
c
c  note: use jpllis in masjpl (otherwise unstable for changes in DExxx)
c
c ==========================================
c EARTH
c get Earth elements
      SUBROUTINE earth(t0,eqp)
      IMPLICIT NONE
c input: epoch time
      DOUBLE PRECISION t0
c output: elements of Earth (equinoctal, ecliptic)
      DOUBLE PRECISION eqp(6)
c =============JPL EPHEM===============
c data for masses
      INCLUDE 'jplhdr.h'
c output of JPL routine, Julian date, rotation matrix
      double precision et(2),rot(3,3),rrd(6),xea(6),enne
c integers for call to JPl routines
      integer ntarg,ncent,istate
c masses
      INCLUDE 'sunmass.h'
      INCLUDE 'parbep.h'
      INCLUDE 'masses.h'
c ====================================
c JPL Earth vector at observation time
      et(1)=2400000.5d0
      et(2)=t0
      ntarg=3
      ncent=11
c duplicate computation of gmse, in case masjpl has not been called yet
      gmse=gms*(1.d0+cval(11)/cval(18))
* ****** added on Sat Jun 14 1997 ******
* first istate need to be=2  (dpleph calculates also vel.)
      istate=2
      call dpleph(et,ntarg,ncent,rrd,istate)
* Change of reference system EQUM00 ---> ECLM00
      call rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0)
      call prodmv(xea,rot,rrd)
      call prodmv(xea(4),rot,rrd(4))
c elements of Earth
      call coocha(xea,'CAR',gmse,eqp,'EQU',enne)
      RETURN
      END
c ======================================================================
c EARCAR - get Earth cartesian coordinates (ecliptic J2000)
c ======================================================================
      SUBROUTINE earcar(t0,xea,ifla)
      IMPLICIT NONE
c input: epoch time, flag for getting Earth (heliocentric; ifla=1)
c        or Sun (barycentric; ifla=2)
      DOUBLE PRECISION t0
      INTEGER ifla
c output: heliocentric state vector of Earth (equinoctal, ecliptic)
c      or barycentric state vector of Sun (equinoctal, ecliptic)
      DOUBLE PRECISION xea(6)
c =============JPL EPHEM===============
c data for masses
      INCLUDE 'jplhdr.h'
c output of JPL routine, Julian date, rotation matrix
      double precision et(2),rot(3,3),rrd(6)
c integers for call to JPl routines
      integer ntarg,ncent,istate
c ====================================
c JPL Earth vector at observation time
      et(1)=2400000.5d0
      et(2)=t0
      if (ifla.eq.1) then
       ntarg=3
       ncent=11
      else
       ntarg=11
       ncent=12
      endif
* ****** added on Sat Jun 14 1997 ******
* first istate need to be=2  (dpleph calculates also vel.)
      istate=2
      call dpleph(et,ntarg,ncent,rrd,istate)
* Change of reference system EQUM00 ---> ECLM00
      call rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0)
      call prodmv(xea,rot,rrd)
      call prodmv(xea(4),rot,rrd(4))
      RETURN
      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 20, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          R D B E P                            *
*  *                                                               *
*  *         Reads and interpolates binary ephemeris files         *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    T         -  Time (MJD, TDT)
*           NB        -  Number of bodies
*           ID        -  Identity (order number) of requested bodies
*
* OUTPUT:   X,V       -  Position and velocity vectors of requested
*                        bodies (refsys: ECLM J2000)
*
* EXAMPLE:
*     NB=2, ID(1)=3, ID(2)=1 means that the third and first bodies
*           contained in the binary ephemeris files are requested
*     on output, X(i,1) will contain the position vector of body
*     number 3 and X(i,2) the position vector of body number 1
*
      SUBROUTINE rdbep(t,nb,id,x,v)
      IMPLICIT NONE

      INCLUDE 'parbep.h'
      INCLUDE 'combep.h'

      CHARACTER*4 kind
      PARAMETER (kind='AM  ')

      INTEGER nb
      INTEGER id(nb)
      DOUBLE PRECISION t,x(3,nb),v(3,nb)

      INTEGER unit,nhl,vers,n,recl,i,k,ib,lf
      DOUBLE PRECISION t1,t2,dt,gma(nbepx),gma1(nbepx)
      DOUBLE PRECISION c1,c2
      DOUBLE PRECISION xv1(6),xv2(6),enne
      LOGICAL first

      SAVE first,unit,nhl,n,t1,t2,dt,gma,gma1

* Buffer
      INTEGER ipt1,ipt2
      DOUBLE PRECISION tb1,tb2,elem1(6,nbepx),elem2(6,nbepx)
      DOUBLE PRECISION enne1(nbepx),enne2(nbepx),m01(nbepx),m02(nbepx)

      SAVE ipt1,ipt2,tb1,tb2,elem1,elem2,enne1,enne2,m01,m02

      INTEGER lench
      EXTERNAL lench
      
      INTEGER ivel

      DATA first/.true./

      IF(first) THEN
          CALL filass(unit,filbep)

* Open with a small record length (just to read the number of bodies)
          OPEN(unit,FILE=filbep,ACCESS='DIRECT',RECL=24,STATUS='OLD')
          READ(unit,REC=1) vers,nbep,n
          IF(vers.NE.102) THEN
              lf=lench(filbep)
              WRITE(*,120) filbep(1:lf),vers
              STOP '**** rdbep: abnormal end ****'
          END IF
          CLOSE(unit)
          CALL chkpdf(nbep,nbepx,'nbepx')
          IF(n.LT.2) STOP '**** rdbep: n < 2 ****'
          nhl=2+2*nbep

* Open again with the correct record length
          recl=8*6*nbep
          OPEN(unit,FILE=filbep,ACCESS='DIRECT',RECL=recl,STATUS='OLD')
          READ(unit,REC=2) t1,t2,dt

* Read masses
          DO 1 i=1,nbep
             READ(unit,REC=2+i) masbep(i),gma(i),gma1(i)
c            WRITE(*,*)i,masbep(i),gma(i),gma1(i)
 1        CONTINUE

* Initialization of buffer
          ipt1=1
          ipt2=2
          READ(unit,REC=nhl+ipt1) ((elem1(i,k),i=1,6),k=1,nbep)
          READ(unit,REC=nhl+ipt2) ((elem2(i,k),i=1,6),k=1,nbep)
          tb1=t1+dt*(ipt1-1)
          tb2=t1+dt*(ipt2-1)
          DO 4 k=1,nbep
          enne1(k)=SQRT(gma1(k)/elem1(1,k)**3)
          enne2(k)=SQRT(gma1(k)/elem2(1,k)**3)
          m01(k)=elem1(6,k)
          m02(k)=elem2(6,k)
 4        CONTINUE

          first=.false.
      END IF
 120  FORMAT(' ERROR: file ',A,' has unsupported version',I4)

      IF(t.LT.tb1 .OR. t.GT.tb2) THEN
          IF(t.LT.t1 .OR. t.GT.t2)
     +        STOP '**** rdbep: T is out of bounds ****'
          ipt1=(t-t1)/dt+1
          ipt2=ipt1+1
          IF(ipt1.LT.1) STOP '**** rdbep: internal error (01) ****'
          IF(ipt2.GT.n) STOP '**** rdbep: internal error (02) ****'
          READ(unit,REC=nhl+ipt1) ((elem1(i,k),i=1,6),k=1,nbep)
          READ(unit,REC=nhl+ipt2) ((elem2(i,k),i=1,6),k=1,nbep)
          tb1=t1+dt*(ipt1-1)
          tb2=t1+dt*(ipt2-1)
          IF(t.LT.tb1 .OR. t.GT.tb2)
     +        STOP '**** rdbep: internal error (03) ****'
          DO 5 k=1,nbep
          enne1(k)=SQRT(gma1(k)/elem1(1,k)**3)
          enne2(k)=SQRT(gma1(k)/elem2(1,k)**3)
          m01(k)=elem1(6,k)
          m02(k)=elem2(6,k)
 5        CONTINUE
      END IF

* Coefficients for interpolation
      c1=(tb2-t)/dt
      c2=(t-tb1)/dt

      DO 2 ib=1,nb
      k=id(ib)
      elem1(6,k)=m01(k)+enne1(k)*(t-tb1)
      elem2(6,k)=m02(k)+enne2(k)*(t-tb2)
      ivel=1
      CALL kepcar(elem1(1,k),gma1(k),ivel,xv1)
      CALL kepcar(elem2(1,k),gma1(k),ivel,xv2)
c     CALL coocha(elem1(1,k),'KEP',gma1(k),xv1,'CAR',enne)
c     CALL coocha(elem2(1,k),'KEP',gma1(k),xv2,'CAR',enne)
      DO 3 i=1,3
      x(i,ib)=c1*xv1(i)+c2*xv2(i)
      v(i,ib)=c1*xv1(i+3)+c2*xv2(i+3)
 3    CONTINUE
 2    CONTINUE

      END
c
      subroutine fszer2(nrecl,ksize,nrfile,namfil)
c
c++++++++++++++++++++++++
c  this subroutine opens the file, 'namfil', with a phony record length, reads 
c  the first record, and uses the info to compute ksize, the number of single 
c  precision words in a record.  
c
c  the subroutine also sets the values of  nrecl, nrfile, and namfil.

* changed by sabrina baccili on Wed Oct 30
*      save

      implicit double precision(a-h,o-z)
      character*6 ttl(14,3),cnam(400)
      character*(*) namfil
      character*150 namtmp
      logical found
c     logical fail,fail1

      dimension ss(3)

      integer ipt(3,12),lpt(3)
      save
* end change

* NEEDED common blocks:
      INCLUDE 'comlib.h'

c  *****************************************************************
c  *****************************************************************
c
c  the parameters nrecl, nrfile, and namfil are to be set by the user
c
c  *****************************************************************

c  nrecl=1 if "recl" in the open statement is the record length in s.p. words
c  nrecl=4 if "recl" in the open statement is the record length in bytes
c  (for unix, it is probably 4)
c
      nrecl=4

      IF(iiclib.NE.36) STOP '**** fszer2: internal error (01) ****'
c changed by A. Milani (October 24, 1998)
c  binary ephemeris file: old method to search for it with a name namfil
c was used in bineph, now incompatibele with fitobs/orbfit
c      fail=.false.
c     CALL rdncha('JPLDE.','file',namfil,.false.,found,fail1,fail)
c      IF(fail1) STOP '**** fszer2: abnormal end ****'
c      IF(.NOT.found) THEN
c  jpleph is alwaysthe external name of the binary ephemeris file
          namfil='jpleph'
          INQUIRE(FILE=namfil,EXIST=found)
          IF(found) GOTO 2
          namtmp=libdir(1:lenld)//namfil
          namfil=namtmp
          INQUIRE(FILE=namfil,EXIST=found)
          IF(found) GOTO 2
          GOTO 10
c      END IF
 2    CONTINUE
c  nrfile is the internal unit number used for the ephemeris file
      CALL filass(nrfile,namfil)
* end of change

c  *****************************************************************
c  *****************************************************************

c  **  open the direct-access file and get the pointers in order to 
c  **  determine the size of the ephemeris record

      mrecl=nrecl*1000

        open(nrfile,
     *       file=namfil,
     *       access='direct',
     *       form='unformatted',
     *       recl=mrecl,
     *       status='old',
     *       err=10)

      read(nrfile,rec=1)ttl,cnam,ss,ncon,au,emrat,ipt,numde,lpt

      close(nrfile)

c  find the number of ephemeris coefficients from the pointers

      kmx = 0
      khi = 0

      do 1 i = 1,12
         if (ipt(1,i) .gt. kmx) then
            kmx = ipt(1,i)
            khi = i
         endif
 1    continue
      if (lpt(1) .gt. kmx) then
          kmx = lpt(1)
          khi = 13
      endif

      nd = 3
      if (khi .eq. 12) nd=2

      if(khi.eq.13) then
          ksize = 2*(lpt(1)+nd*lpt(2)*lpt(3)-1)
      else
          ksize = 2*(ipt(1,khi)+nd*ipt(2,khi)*ipt(3,khi)-1)
      endif

      return

 10   continue
      lf=lench(namfil)
      WRITE(*,100) namfil(1:lf)
 100  FORMAT('ERROR opening file ',A)
      STOP '**** fszer2: abnormal end ****'

      end
c++++++++++++++++++++++++++
c
      subroutine dpleph(et2z,ntarg,ncent,rrd,istate)
c                 pleph ( et, ntarg, ncent, rrd, istate )
c
c  note: (A. Milani, March 11, 1998)
c  we have added a new fifth argument to control the interpolation of 
c  derivatives
c
c     this subroutine reads the jpl planetary ephemeris
c     and gives the position and velocity of the point 'ntarg'
c     with respect to 'ncent'.
c
c     calling sequence parameters:
c
c       et = d.p. julian ephemeris date at which interpolation
c            is wanted.
c
c       ** note the entry dpleph for a doubly-dimensioned time **
c          the reason for this option is discussed in the 
c          subroutine state
c
c     ntarg = integer number of 'target' point.
c
c     ncent = integer number of center point.
c
c            the numbering convention for 'ntarg' and 'ncent' is:
c
c                1 = mercury           8 = neptune
c                2 = venus             9 = pluto
c                3 = earth            10 = moon
c                4 = mars             11 = sun
c                5 = jupiter          12 = solar-system barycenter
c                6 = saturn           13 = earth-moon barycenter
c                7 = uranus           14 = nutations (longitude and obliq)
c                            15 = librations, if on eph file
c
c             (if nutations are wanted, set ntarg = 14. for librations,
c              set ntarg = 15. set ncent=0.)
c
c      rrd = output 6-word d.p. array containing position and velocity
c            of point 'ntarg' relative to 'ncent'. the units are au and
c            au/day. for librations the units are radians and radians
c            per day. in the case of nutations the first four words of
c            rrd will be set to nutations and rates, having units of
c            radians and radians/day.
c
c            the option is available to have the units in km and km/sec.
c            for this, set km=.true. in the stcomx common block.
c
      implicit double precision (a-h,o-z)
      dimension rrd(6),et2z(2),et2(2),pv(6,13)
      dimension pvsun(6)
c header
      logical bsave,km,bary
      common/stcomx/km,bary,pvsun
c standard header  
      include 'jplhdr.h'
      integer list(12)
* memory model static
      save

      et2(1)=et2z(1)
      et2(2)=et2z(2)

  11  ettot=et2(1)+et2(2)

      do i=1,6
        rrd(i)=0.d0
      enddo

  96  if(ntarg .eq. ncent) return

      do i=1,12
        list(i)=0
      enddo
c
c     check for nutation call
c
      if(ntarg.eq.14)then
        if(ipt(2,12).gt.0) then
          list(11)=2
          call state(et2,list,pv,rrd,istate)
          return
        else
          write(6,297)
  297     format(' *****  no nutations on the ephemeris file  *****')
          stop
        endif
      endif
c
c     check for librations
c
      if(ntarg.eq.15)then
        if(lpt(2).gt.0) then
          list(12)=2
          call state(et2,list,pv,rrd,istate)
          do i=1,6
          rrd(i)=pv(i,11)
          enddo
          return
        else
          write(6,298)
  298     format(' *****  no librations on the ephemeris file  *****')
          stop
        endif
      endif
c
c       force barycentric output by 'state'
      bsave=bary
      bary=.true.
c
c       set up proper entries in 'list' array for state call
c
      do i=1,2
        k=ntarg
        if(i .eq. 2) k=ncent
        if(k .le. 10) list(k)=2
        if(k .eq. 10) list(3)=2
        if(k .eq. 3) list(10)=2
        if(k .eq. 13) list(3)=2
      enddo
c
c       make call to state
c
      call state(et2,list,pv,rrd,istate)
c
      if(ntarg .eq. 11 .or. ncent .eq. 11) then
         do i=1,6
           pv(i,11)=pvsun(i)
         enddo
      endif
c
      if(ntarg .eq. 12 .or. ncent .eq. 12) then
         do i=1,6
           pv(i,12)=0.d0
         enddo
      endif

      if(ntarg .eq. 13 .or. ncent .eq. 13) then
         do i=1,6
           pv(i,13)=pv(i,3)
         enddo
      endif

      if(ntarg*ncent .eq. 30 .and. ntarg+ncent .eq. 13) then
         do i=1,6
           pv(i,3)=0.d0
         enddo
         go to 99
      endif

      if(list(3) .eq. 2) then
         do i=1,6
           pv(i,3)=pv(i,3)-pv(i,10)/(1.d0+emrat)
         enddo
      endif

      if(list(10) .eq. 2) then
         do i=1,6
           pv(i,10)=pv(i,3)+pv(i,10)
         enddo
      endif

  99  do i=1,6
        rrd(i)=pv(i,ntarg)-pv(i,ncent)
      enddo

      bary=bsave

      return
      end
c++++++++++++++++++++++++++++++++
c
      subroutine state(et2,list,pv,pnut,istate)
c
c++++++++++++++++++++++++++++++++
c
c this subroutine reads and interpolates the jpl planetary ephemeris file
c
c     calling sequence parameters:
c
c     input:
c
c         et2   dp 2-word julian ephemeris epoch at which interpolation
c               is wanted.  any combination of et2(1)+et2(2) which falls
c               within the time span on the file is a permissible epoch.
c
c                a. for ease in programming, the user may put the
c                   entire epoch in et2(1) and set et2(2)=0.
c
c                b. for maximum interpolation accuracy, set et2(1) =
c                   the most recent midnight at or before interpolation
c                   epoch and set et2(2) = fractional part of a day
c                   elapsed between et2(1) and epoch.
c
c                c. as an alternative, it may prove convenient to set
c                   et2(1) = some fixed epoch, such as start of integration,
c                   and et2(2) = elapsed interval between then and epoch.
c
c        list   12-word integer array specifying what interpolation
c               is wanted for each of the bodies on the file.
c
c                         list(i)=0, no interpolation for body i
c                                =1, position only
c                                =2, position and velocity
c
c               the designation of the astronomical bodies by i is:
c
c                         i = 1: mercury
c                           = 2: venus
c                           = 3: earth-moon barycenter
c                           = 4: mars
c                           = 5: jupiter
c                           = 6: saturn
c                           = 7: uranus
c                           = 8: neptune
c                           = 9: pluto
c                           =10: geocentric moon
c                           =11: nutations in longitude and obliquity
c                           =12: lunar librations (if on file)
c
c
c     output:
c
c          pv   dp 6 x 11 (note: 6 x 12)
c               array that will contain requested interpolated
c               quantities.  the body specified by list(i) will have its
c               state in the array starting at pv(1,i).  (on any given
c               call, only those words in 'pv' which are affected by the
c               first 10 'list' entries (and by list(12) if librations are
c               on the file) are set.  the rest of the 'pv' array
c               is untouched.)  the order of components starting in
c               pv(1,i) is: x,y,z,dx,dy,dz.
c
c               all output vectors are referenced to the earth mean
c               equator and equinox of epoch. the moon state is always
c               geocentric; the other nine states are either heliocentric
c               or solar-system barycentric, depending on the setting of
c               common flags (see below).
c
c               lunar librations, if on file, are put into pv(k,11) if
c               list(12) is 1 or 2.
c
c         nut   dp 4-word array that will contain nutations and rates,
c               depending on the setting of list(11).  the order of
c               quantities in nut is:
c
c                        d psi  (nutation in longitude)
c                        d epsilon (nutation in obliquity)
c                        d psi dot
c                        d epsilon dot
c
c           *   statement # for error return, in case of epoch out of
c               range or i/o errors.
c
c
c     common area stcomx:
c
c          km   logical flag defining physical units of the output
c               states. km = .true., km and km/sec
c                          = .false., au and au/day
c               default value = .false.  (km determines time unit
c               for nutations and librations.  angle unit is always radians.)
c
c        bary   logical flag defining output center.
c               only the 9 planets are affected.
c                        bary = .true. =\ center is solar-system barycenter
c                             = .false. =\ center is sun
c               default value = .false.
c
c       pvsun   dp 6-word array containing the barycentric position and
c               velocity of the sun.
c
c
c change by A.Milani and Z. Knezevic, March 10, 1998, for compatibility with 
c Lahey compiler
      IMPLICIT NONE
      DOUBLE PRECISION  et2(2),pv(6,12),pnut(4),t(2),pjd(4),buf(1500)
      integer list(12)
      logical first
c
      character*6 ttl(14,3),cnam(400)
      common/chrhdr/cnam,ttl
c
      character*150 namfil
c
      DOUBLE PRECISION  pvsun(6)
      logical km,bary
      common/stcomx/km,bary,pvsun
c variable declared for implicit none
      INTEGER istate
      INTEGER nrecl, ksize, nrfile, irecsz, ncoeffs, numde, nrl
      INTEGER i, j, nr, k
      DOUBLE PRECISION s, aufac
c
      include 'jplhdr.h'
      save
      data first/.true./
c end change 1998
c
c       entry point - 1st time in, get pointer data, etc., from eph file
c
      if(first) then
        first=.false.

c ************************************************************************
c ************************************************************************

c the user must select one of the following by deleting the 'c' in column 1

c ************************************************************************

c        call fszer1(nrecl,ksize,nrfile,namfil)
         call fszer2(nrecl,ksize,nrfile,namfil)
c        call fszer3(nrecl,ksize,nrfile,namfil)

      if(nrecl .eq. 0) write(*,*)'  ***** fszer is not working *****'

c ************************************************************************
c ************************************************************************

      irecsz=nrecl*ksize
      ncoeffs=ksize/2

        open(nrfile,
     *       file=namfil,
     *       access='direct',
     *       form='unformatted',
     *       recl=irecsz,
     *       status='old',
     *       err=10)

      read(nrfile,rec=1)ttl,cnam,ss,ncon,au,emrat,
     . ((ipt(i,j),i=1,3),j=1,12),numde,lpt

      read(nrfile,rec=2)cval

      do i=1,3
      ipt(i,13)=lpt(i)
      enddo
      nrl=0

      endif
c
c       ********** main entry point **********
c
      if(et2(1) .eq. 0.d0) return
c
      s=et2(1)-.5d0
      call split(s,pjd(1))
      call split(et2(2),pjd(3))
      pjd(1)=pjd(1)+pjd(3)+.5d0
      pjd(2)=pjd(2)+pjd(4)
      call split(pjd(2),pjd(3))
      pjd(1)=pjd(1)+pjd(3)
c
c       error return for epoch out of range
c
      if(pjd(1)+pjd(4).lt.ss(1) .or. pjd(1)+pjd(4).gt.ss(2)) go to 98
c
c       calculate record # and relative time in interval
c
      nr=idint((pjd(1)-ss(1))/ss(3))+3
      if(pjd(1).eq.ss(2)) nr=nr-1
      t(1)=((pjd(1)-(dble(nr-3)*ss(3)+ss(1)))+pjd(4))/ss(3)
c
c       read correct record if not in core
c
      if(nr.ne.nrl) then
        nrl=nr
        read(nrfile,rec=nr,err=99)(buf(k),k=1,ncoeffs)
      endif

      if(km) then
      t(2)=ss(3)*86400.d0
      aufac=1.d0
      else
      t(2)=ss(3)
      aufac=1.d0/au
      endif
c
c   interpolate ssbary sun
* ****** changed on Sat Jun 14 1997 ******
*      call interp(buf(ipt(1,11)),t,ipt(2,11),3,ipt(3,11),2,pvsun)
      call interp(buf(ipt(1,11)),t,ipt(2,11),3,ipt(3,11),istate,pvsun)
* **************************************

      do i=1,6
      pvsun(i)=pvsun(i)*aufac
      enddo
c
c   check and interpolate whichever bodies are requested
c
      do 4 i=1,10
        if(list(i).eq.0) go to 4
* ****** changed on Sat Jun 14 1997 ******
*      call interp(buf(ipt(1,i)),t,ipt(2,i),3,ipt(3,i),
*     & list(i),pv(1,i))
        call interp(buf(ipt(1,i)),t,ipt(2,i),3,ipt(3,i),
     & istate,pv(1,i))
* **************************************
        do j=1,6
          if(i.le.9 .and. .not.bary) then
             pv(j,i)=pv(j,i)*aufac-pvsun(j)
          else
             pv(j,i)=pv(j,i)*aufac
          endif
        enddo
c
   4  continue
c
c       do nutations if requested (and if on file)
c
      if(list(11).gt.0 .and. ipt(2,12).gt.0)
     * call interp(buf(ipt(1,12)),t,ipt(2,12),2,ipt(3,12),
     * list(11),pnut)
c
c       get librations if requested (and if on file)
c
      if(list(12).gt.0 .and. ipt(2,13).gt.0)
     * call interp(buf(ipt(1,13)),t,ipt(2,13),3,ipt(3,13),
     * list(12),pv(1,11))
c
      return
c
  98  write(*,198)et2(1)+et2(2),ss(1),ss(2)
 198  format(' ***  requested jed,',f12.2,
     * ' not within ephemeris limits,',2f12.2,'  ***')
c
      return
c
   99 write(*,'(2f12.2,a80)')
     & et2,'error return in state'
c
      stop
c
 10   continue
      STOP '**** state: error opening JPL DE file ****'
c
      end
c+++++++++++++++++++++++++++++++++
c
      subroutine interp(buf,t,ncf,ncm,na,ifl,pv)
c
c+++++++++++++++++++++++++++++++++
c
c     this subroutine differentiates and interpolates a
c     set of chebyshev coefficients to give position and velocity
c
c     calling sequence parameters:
c
c       input:
c
c         buf   1st location of array of d.p. chebyshev coefficients of position
c
c           t   t(1) is dp fractional time in interval covered by
c               coefficients at which interpolation is wanted
c               (0 .le. t(1) .le. 1).  t(2) is dp length of whole
c               interval in input time units.
c
c         ncf   # of coefficients per component
c
c         ncm   # of components per set of coefficients
c
c          na   # of sets of coefficients in full array
c               (i.e., # of sub-intervals in full interval)
c
c          ifl  integer flag: =1 for positions only
c                             =2 for pos and vel
c
c
c       output:
c
c         pv   interpolated quantities requested.  dimension
c               expected is pv(ncm,ifl), dp.
c
c
      implicit double precision (a-h,o-z)
c
      save
c CHANGE 27-9-2001 because of out of bounds
c      double precision buf(ncf,ncm,*),t(2),pv(ncm,ifl),pc(18),vc(18)
      double precision buf(ncf,ncm,*),t(2),pv(ncm,*),pc(18),vc(18)
c
      data np/2/
      data nv/3/
      data twot/0.d0/
      data pc(1),pc(2)/1.d0,0.d0/
      data vc(2)/1.d0/
c
c       entry point. get correct sub-interval number for this set
c       of coefficients and then get normalized chebyshev time
c       within that subinterval.
c
      dna=dble(na)
      dt1=dint(t(1))
      temp=dna*t(1)
      l=idint(temp-dt1)+1

c         tc is the normalized chebyshev time (-1 .le. tc .le. 1)

      tc=2.d0*(dmod(temp,1.d0)+dt1)-1.d0

c       check to see whether chebyshev time has changed,
c       and compute new polynomial values if it has.
c       (the element pc(2) is the value of t1(tc) and hence
c       contains the value of tc on the previous call.)

      if(tc.ne.pc(2)) then
        np=2
        nv=3
        pc(2)=tc
        twot=tc+tc
      endif
c
c       be sure that at least 'ncf' polynomials have been evaluated
c       and are stored in the array 'pc'.
c
      if(np.lt.ncf) then
        do  i=np+1,ncf
          pc(i)=twot*pc(i-1)-pc(i-2)
        enddo
        np=ncf
      endif
c
c       interpolate to get position for each component
c
      do 2 i=1,ncm
        pv(i,1)=0.d0
        do 3 j=ncf,1,-1
c =====================================================
c OUT OF BOUNDS 
c PGF90-F-Subscript out of range for array buf (jplsub.f: 297)
c    subscript=2, lower bound=1, upper bound=1, dimension=3

          pv(i,1)=pv(i,1)+pc(j)*buf(j,i,l)
c =====================================================
    3   continue
    2 continue
c modification to avoid computing derivatives when not needed
      IF(ifl.le.1)THEN
         DO  i=1,ncm
c =====================================================
c OUT OF BOUNDS 
c 0: Subscript out of range for array pv (jpl_ephem.f: 912)
c    subscript=2, lower bound=1, upper bound=1, dimension=2
           pv(i,2)=0.d0
c =====================================================
         ENDDO
         RETURN
      ENDIF
c
c       if velocity interpolation is wanted, be sure enough
c       derivative polynomials have been generated and stored.
c
      vfac=(dna+dna)/t(2)
      vc(3)=twot+twot
      if(nv.lt.ncf) then
        do 4 i=nv+1,ncf
          vc(i)=twot*vc(i-1)+pc(i-1)+pc(i-1)-vc(i-2)
    4   continue
        nv=ncf
      endif
c
c       interpolate to get velocity for each component
c
      do 5 i=1,ncm
        pv(i,2)=0.d0
        do 6 j=ncf,2,-1
          pv(i,2)=pv(i,2)+vc(j)*buf(j,i,l)
    6   continue
        pv(i,2)=pv(i,2)*vfac
    5 continue
c
      return
c
      end
c+++++++++++++++++++++++++++++
c
      subroutine const(nam,val,sss,n)
c
c+++++++++++++++++++++++++++++
c
c     this entry obtains the constants from the ephemeris file
c
c     calling seqeunce parameters (all output):
c
c       nam = character*6 array of constant names
c
c       val = d.p. array of values of constants
c
c       sss = d.p. jd start, jd stop, step of ephemeris
c
c         n = integer number of entries in 'nam' and 'val' arrays
c
* changed by sabrina baccili on Wed Oct 30
*      save
c
      implicit double precision (a-h,o-z)
      character*6 nam(*)
      double precision val(*),sss(3)
c ***
      dimension pp1(2),ipp2(12),pp3(6,12),pp4(4)
c ***
      include 'jplhdr.h'

      character*6 ttl(14,3),cnam(400)
      common/chrhdr/cnam,ttl
      save
c  call state to initialize the ephemeris and read in the constants
* end change
c ***
      pp1(1)=0.d0
      pp1(2)=0.d0
      do 321 ijk=1,12
 321     ipp2(ijk)=0
      do 322 ijk=1,6
         do 323 ijj=1,12
 323        pp3(ijk,ijj)=0.d0
 322  continue
      do 324 iki=1,4
 324     pp4(iki)=0.d0
      istate =2
      call state(pp1,ipp2,pp3,pp4,istate)
c      call state(0.d0,0,0,0.d0)
c ***
      n=ncon

      do i=1,3
        sss(i)=ss(i)
      enddo
c
      do i=1,n
        nam(i)=cnam(i)
        val(i)=cval(i)
      enddo
c
      return
      end
c+++++++++++++++++++++++++
c
      subroutine split(tt,fr)
c
c+++++++++++++++++++++++++
c
c     this subroutine breaks a d.p. number into a d.p. integer
c     and a d.p. fractional part.
c
c     calling sequence parameters:
c
c       tt = d.p. input number
c
c       fr = d.p. 2-word output array.
c            fr(1) contains integer part
c            fr(2) contains fractional part
c
c            for negative input numbers, fr(1) contains the next
c            more negative integer; fr(2) contains a positive fraction.
c
c       calling sequence declarations
c
      implicit double precision (a-h,o-z)

      dimension fr(2)
      save

c       main entry -- get integer and fractional parts

      fr(1)=dint(tt)
      fr(2)=tt-fr(1)

      if(tt.ge.0.d0 .or. fr(2).eq.0.d0) return

c       make adjustments for negative input number

      fr(1)=fr(1)-1.d0
      fr(2)=fr(2)+1.d0

      return
      end
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: October 13, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         J P L L I S                           *
*  *                                                               *
*  *      Get the list of masses and IDs from JPL DE header        *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    NAMES     -  Planet names
*           N         -  Number of planets
*
* OUTPUT:   GMP       -  G*M(planet)
*           ID        -  Planet ID number
*           FAIL      -  Error flag
*
      SUBROUTINE jpllis(names,n,gmp,id,fail)
      IMPLICIT NONE

      INTEGER n
      CHARACTER*(*) names(n)
      DOUBLE PRECISION gmp(n)
      INTEGER id(n)
      LOGICAL fail

      DOUBLE PRECISION et2(2),pv(6,12),pnut(4)
      INTEGER list(12),i,ln,emopt
      CHARACTER*12 n1
      LOGICAL emwarn

* JPLDE header
      DOUBLE PRECISION cval(400),ss(3),au,emrat
      INTEGER denum,ncon,ipt(3,13),lpt(3)
      CHARACTER*6 cnam(400),ttl(14,3)
      COMMON/ephhdr/cval,ss,au,emrat,denum,ncon,ipt,lpt
      COMMON/chrhdr/cnam,ttl

      INTEGER lench
      EXTERNAL lench

      DATA et2/2*0.d0/
      DATA list/12*0/

* Dummy call to STATE for reading JPLDE header
      CALL state(et2,list,pv,pnut,1)

* Choice for Earth/Moon: 1=barycenter, 2=distinct bodies
      emopt=0
      emwarn=.false.
      DO 1 i=1,n
      n1=names(i)
      CALL upcase(n1)
      IF(n1.EQ.'MERCURY') THEN
          id(i)=1
          gmp(i)=cval(9)
          IF(cnam(9).NE.'GM1')
     +        STOP '**** jpllis: internal error (C01) ****'
      ELSEIF(n1.EQ.'VENUS') THEN
          id(i)=2
          gmp(i)=cval(10)
          IF(cnam(10).NE.'GM2')
     +        STOP '**** jpllis: internal error (C02) ****'
      ELSEIF(n1.EQ.'EARTH') THEN
          id(i)=3
          gmp(i)=cval(11)*emrat/(1+emrat)
          IF(cnam(11).NE.'GMB')
     +        STOP '**** jpllis: internal error (C03) ****'
          IF(emopt.NE.0 .AND. emopt.NE.2) emwarn=.true.
          emopt=2
      ELSEIF(n1.EQ.'MOON') THEN
          id(i)=10
          gmp(i)=cval(11)/(1+emrat)
          IF(cnam(11).NE.'GMB')
     +        STOP '**** jpllis: internal error (C10) ****'
          IF(emopt.NE.0 .AND. emopt.NE.2) emwarn=.true.
          emopt=2
      ELSEIF(n1.EQ.'EARTH+MOON') THEN
          id(i)=13
          gmp(i)=cval(11)
          IF(cnam(11).NE.'GMB')
     +        STOP '**** jpllis: internal error (C13) ****'
          IF(emopt.NE.0 .AND. emopt.NE.1) emwarn=.true.
          emopt=1
      ELSEIF(n1.EQ.'MARS') THEN
          id(i)=4
          gmp(i)=cval(12)
          IF(cnam(12).NE.'GM4')
     +        STOP '**** jpllis: internal error (C04) ****'
      ELSEIF(n1.EQ.'JUPITER') THEN
          id(i)=5
          gmp(i)=cval(13)
          IF(cnam(13).NE.'GM5')
     +        STOP '**** jpllis: internal error (C05) ****'
      ELSEIF(n1.EQ.'SATURN') THEN
          id(i)=6
          gmp(i)=cval(14)
          IF(cnam(14).NE.'GM6')
     +        STOP '**** jpllis: internal error (C06) ****'
      ELSEIF(n1.EQ.'URANUS') THEN
          id(i)=7
          gmp(i)=cval(15)
          IF(cnam(15).NE.'GM7')
     +        STOP '**** jpllis: internal error (C07) ****'
      ELSEIF(n1.EQ.'NEPTUNE') THEN
          id(i)=8
          gmp(i)=cval(16)
          IF(cnam(16).NE.'GM8')
     +        STOP '**** jpllis: internal error (C08) ****'
      ELSEIF(n1.EQ.'PLUTO') THEN
          id(i)=9
          gmp(i)=cval(17)
          IF(cnam(17).NE.'GM9')
     +        STOP '**** jpllis: internal error (C09) ****'
      ELSE
          ln=lench(names(i))
          WRITE(*,100) names(i)(1:ln)
          fail=.true.
      END IF
 100  FORMAT(' ERROR: ',A,' is unknown among JPL planets')
 1    CONTINUE

      IF(emwarn) THEN
          WRITE(*,101)
          fail=.true.
      END IF
 101  FORMAT(' ERROR in the list of JPL planets: please DO NOT select'/
     +       '       "Earth+Moon" (Earth-Moon barycenter) together'/
     +       '       with "Earth" and/or "Moon"')

      END
c ================================
c MASJPL 
c
c version 1.4, A. Milani, Nov. 10, 1997
c 1.5.2 corrected January 3, 1997
c  input of planetary masses from jpl header
c  it assumes that dpleph has already been called
c ================================
      subroutine masjpl
      implicit none
c controls of the force model
      include 'model.h'
c controls of close approaches monitoring
      include 'iclap.h'
c name of asteroid binary file
      include 'bifina.h'
c asteroid parameter and common
      include 'parbep.h'
      include 'combep.h'
      include 'selast.h'
c planetary masses, constants, etc.
      include 'masses.h'
c ======== constant of gravitation ==============
      include 'sunmass.h'
      double precision gmerc,gm1,gmsun
c ======== JPL EPHEM ============================
c  jpl header data read by dpleph and subr.  
      character*6 cnam(400),ttl(14,3)
      common/chrhdr/cnam,ttl
c standard JPL header
      include 'jplhdr.h'
c ===============================================
c  loop indexes, lengths, output units
      integer m,mm,ia,i,j,ln,iabe
      logical openfl
c functions
      integer lench
c  strings with asteroid names
      character*30 astnam(nbepx)
      character*30 string
      integer mlun,mea,mjup,mnep
c  strings with planet names
      character*30 nomi(11)
      data nomi/'MERCURY','VENUS','EARTH_MOON','MARS',
     + 'JUPITER','SATURN','URANUS','NEPTUNE','PLUTO','MOON','EARTH'/

c      LOGICAL first
c      SAVE first,iabe
c      DATA first/.true./

c  number of planets in our integration
c  also sets the distance for close approaches
      npla=7
      if(imerc.gt.0)then
         npla=npla+1
         m=0
         mlun=9
         mea=3
         mjup=5
         mnep=8
      else
         m=1
         mlun=8
         mea=2
         mjup=4
         mnep=7
      endif
      if(iplut.gt.0)then
         npla=npla+1
         mlun=mlun+1
      endif
c  distance defining a close approach: choise in model.opt
      do 39 j=1,npla
c let dist for mercury, venus and mars be determined by 'dter' (default 0.1 AU)
c while for the earth it is determined by 'dmea' (default 0.1 AU)
c if non-default values are desired, set the propag-options .dter and .dmea
c in the option file (e.g., .dter=0.05d0 and .dmea=0.2d0 for NEAs)
        IF(imerc.gt.0)THEN
           IF(j.eq.1.or.j.eq.2.or.j.eq.4)dmin(j)=dter
        ELSE
           IF(j.eq.1.or.j.eq.3)dmin(j)=dter
        ENDIF
        IF(j.eq.mea)dmin(j)=dmea
        IF(j.ge.mjup.and.j.le.mnep)dmin(j)=dmjup 
 39   continue
c
c  gms in solar masses
      gmsun=cval(18)
c  flags to require data
      do  mm=1,12
        listpl(mm)=0
      enddo
c  setup of the planetary masses
      do 40 mm=1,npla
          m=m+1
          gm(mm)=cval(8+m)/gmsun
          ordnam(mm)=nomi(m)
          if(m.eq.3)then
             itarg(mm)=13
          else
             itarg(mm)=m
          endif
          listpl(m)=2
 40   continue
c store mass of Earth-Moon plus Sun (in solar masses)
      gmse=(gmsun+cval(11))/gmsun
c  store mercury/sun mass ratio; Mercury position is needed anyway, 
c to use Sun-Mercury center of mass
      if(imerc.eq.0)then
          gmerc=cval(9)/gmsun
          listpl(1)=2
      endif
c  the moon as possible additional massive body (then replace E+M center
c  of mass with E).
      if(ilun.gt.0)then
         npla=npla+1
c  close approaches to the Moon are not reported separately 
c  from close approaches to Earth
         dmin(mlun)=dmoon
c  mass of the Moon is computed from Earth/Moon mass ratio
         emrat=cval(8)
         gm(mlun)=cval(11)/((1.d0+emrat)*gmsun)
         ordnam(mlun)=nomi(10)
         itarg(mlun)=10
         listpl(10)=2
c  if the Moon is included, the Earth is included as a single body
c  and not as E-M center of mass
         gm(mea)=cval(11)/((1.d0+1.d0/emrat)*gmsun)
         ordnam(mea)=nomi(11)
         itarg(mea)=3
         listpl(3)=2
      endif
c
      if(iast.ne.0)then
c part relative to massives asteroids
c read names of asteroids and asteroid masses from ascii file 'CPV.abe'
        call filopn(iabe,filbec,'old')
c       inquire(file=filbec,opened=openfl)
c       if(.not.openfl) then
c        IF(first) THEN
c            call filopn(iabe,filbec,'old')
c            write(*,*) 'opening filbec'
c            first=.false.
c        endif
        read(iabe,*)
        read(iabe,*)
        do  ia=1,iast
            read(iabe,201,err=202,end=202)masbep(ia),string
            ln=lench(string)
            astnam(ia)=string(1:ln)
c            astid(ia)=ia
        enddo
c        iatrue=iast
        goto 203
 201    FORMAT(1P,E18.10,1X,A)
 202    WRITE(*,*)'masjpl: too many asteroids requested, iast=',iast
        iast=ia-1
        WRITE(*,*)'masjpl: asteroids available ',iast
c203    close(iabe)
 203    call filclo(iabe,' ')
c  asteroid close approach distance
c  temporary choice: 0.2 AU
        do 13 ia=1,iast
 13        dmin(npla+ia)=dmast
        do 11 ia=1,iatrue
c  asteroid names
           ordnam(npla+ia)=astnam(astid(ia))
c  asteroid masses (unita` masse solari)
 11        gm(npla+ia)=masbep(astid(ia))
      endif
c  gmsun in the chosen units (au, day)
      gm0=gms
      nmass=npla+iatrue
      do 12 i=1,nmass
        gm(i)=gm0*gm(i)
 12   continue
      gmse=gmse*gm0
c  if mercury is not included, its mass is included in the mass of the sun
      if(imerc.eq.0) then
         gmerc=gmerc*gm0
         gm1=gm0+gmerc
         gmu=gmerc
         gm0=gm1
      endif

      return
      end
c
      subroutine trange
      implicit none
*
      include 'timespan.h'
      include 'vlight.h'
*
      double precision et2(2),pv(6,12),pnut(4)
      integer list(12)
      double precision tt,deltt
*
* JPL  header
      include 'jplhdr.h'
      character*6 cnam(400),ttl(14,3)
      common/chrhdr/cnam,ttl
*
      data et2/2*0.d0/
      data list/12*0/
*
* Dummy call to STATE for reading JPLDE header
      CALL state(et2,list,pv,pnut,1)
* store in common timespan time span of jpleph
* transformation from JD to MJD
      tejpl1=ss(1)-2400000.5d0
      tejpl2=ss(2)-2400000.5d0
*      write(*,*) tejpl1,tejpl2
*
c Dummy call to deltt to read the ET-UT data
      tt=deltt(50000.d0)
c speed of light (IAU 1976) in km/s
      ckm=299792.458d0
c id in km/day 
      ckm=ckm*8.64d4
c conversion to AU/day
      vlight=ckm/au
      return
      end


