c ===========MODULE ref_system=====================
c REFERENCE SYSTEMS: Carpino 1996-1998
c PUBLIC ROUTINES
c posobs	heliocentric position of the observer (mean equator J2000)
c pvobs		position of the observer wrt Earth center of mass
c pvobs2        position of the observer wrt Earth center of mass, for radar
c CONTAINS
c SUBROUTINES
c rotpn		transformation between different reference systems
c obleq		mean obliquity of ecliptic
c rnut80	nutation matrix according to Wahr (IAU-1980) theory
c prec		precession matrix
c nutn80	nutation angles according to Wahr (IAU 1980) theory
c equequ	equation of the equinoxes
c gmst		Greenwich Mean Sidereal Time as a function of UT1
c obscoo	body-fixed coordinates of an observatory
c
c chkref	check existence of a reference system indicator
c *ch2ref	translation of a character string into reference system
c
c HEADERS
c
c ref_syst.o:   
c               public trig.h proout.h
c               private t2000.h  
c
c WARNING: pvobs2 uses the rotation routine rotpv belonging to iers_ser.f
c
* Copyright (C) 1997-2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: May 31, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         P O S O B S                           *
*  *                                                               *
*  *   Heliocentric position of the observer (mean equator J2000)  *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    TDT       -  Time (MJD, TDT)
*           OBSCOD    -  Observatory code
*           N         -  Number of observations
*
* OUTPUT:   X         -  Heliocentric position of the observer
*
      SUBROUTINE posobs(tdt,obscod,n,x)
      IMPLICIT NONE

      INTEGER n,obscod(n)
      DOUBLE PRECISION tdt(n),x(3,n)

      INTEGER i,k
      DOUBLE PRECISION et2(2),r6(6),rot(3,3),dxe(3),dve(3),dx(3)
      LOGICAL first
      DATA first/.true./

      SAVE first,rot

      IF(first) THEN
          CALL rotpn(rot,'ECLM','J2000',0.d0,'EQUM','J2000',0.d0)
          first=.false.
      END IF

      et2(1)=2400000.5d0

      DO 1 i=1,n
      et2(2)=tdt(i)
      CALL dpleph(et2,3,11,r6,1)
      CALL pvobs(tdt(i),obscod(i),dxe,dve)
      CALL prodmv(dx,rot,dxe)
      DO 2 k=1,3
      x(k,i)=r6(k)+dx(k)
 2    CONTINUE
 1    CONTINUE

      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 4, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          P V O B S                            *
*  *                                                               *
*  *      Position of the observer with respect to the center      *
*  *         of mass of the Earth (mean ecliptic J2000)            *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    T         -  Time (MJD, TDT)
*           IDSTA     -  Observatory code
*
* OUTPUT:   DX,DV     -  Position and velocity of the observer w.r.
*                        to the center of the Earth (mean ecl. J2000)
*
      SUBROUTINE pvobs(t,idsta,dx,dv)
      IMPLICIT NONE

      DOUBLE PRECISION rot(3,3),dxbf(3),dxtod(3)
      DOUBLE PRECISION dx(3),omega(3),dvbf(3),dvtod(3),dv(3)
c station identifiers
      INTEGER idsta
c times
      INTEGER mjd1,mjd2
      DOUBLE PRECISION t,sec1,sec2,tut,gast
c functions
      DOUBLE PRECISION gmst,equequ
      EXTERNAL gmst,equequ
***************************************
c startup
      CHARACTER*20 name
      LOGICAL first

      INCLUDE 'trig.h'

c static memory allocation only for:
      SAVE first,omega
      DATA first/.true./

      IF(first) THEN
          first=.false.
* Earth angular velocity (rad/d)
          omega(1)=0.d0
          omega(2)=0.d0
          omega(3)=dpig*1.00273790934d0
      END IF
***************************************
* Station name 
      CALL obscoo(idsta,dxbf,name)
* Station position and velocity in the body fixed (equatorial) frame
      CALL prvec(omega,dxbf,dvbf)
* ET decomposed in days plus seconds; no trick (every day is 86400 sec)
      mjd1=t
      sec1=(t-mjd1)*86400
* Computation of UT1
      CALL cnvtim(mjd1,sec1,'ET ',mjd2,sec2,'UT1')
      tut=sec2/86400+mjd2
* 
* Greenwich Apparent Sidereal Time = Greenwich Mean Sidereal Time +
* Equation of the Equinoxes
      gast=gmst(tut)+equequ(t)
* Diurnal rotation matrix (transformation from body-fixed to
* true-of-date frames), neglecting polar motion
      CALL rotmt(-gast,rot,3)
* 
      CALL prodmv(dxtod,rot,dxbf)
      CALL prodmv(dvtod,rot,dvbf)
* Station position and velocity in the J2000 (ecliptic) frame
      CALL rotpn(rot,'EQUT','OFDATE',t,'ECLM','J2000',0.d0)
      CALL prodmv(dx,rot,dxtod)
      CALL prodmv(dv,rot,dvtod)

      END

* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 8, 1997
* revised: end 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          P V O B S 2                          *
*  *                                                               *
*  *      Position of the observer with respect to the center      *
*  *         of mass of the Earth (mean ecliptic J2000)            *
*  *         radar version using IERS series                       *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    T         -  Time (MJD, TDT)
*           IDSTA     -  Observatory code
*
* OUTPUT:   DX,DV     -  Position and velocity of the observer w.r.
*                        to the center of the Earth (mean ecl. J2000)
*
      SUBROUTINE pvobs2(t,idsta,dx,dv)
      IMPLICIT NONE

      DOUBLE PRECISION dxbf(3)
      DOUBLE PRECISION dx(3),dvbf(3),dv(3)
c station identifiers
      INTEGER idsta
c times
      INTEGER mjd1,mjdjpl
      DOUBLE PRECISION t,sec1,secjpl
c indexes
      INTEGER i
***************************************
c station name
      CHARACTER*20 name
***************************************
* Station name 
      CALL obscoo(idsta,dxbf,name)
      DO 1 i=1,3
      dvbf(i)=0.d0
 1    CONTINUE
c date of J2000
      mjdjpl = 51544
      secjpl = 43200.d0 
* ET decomposed in days plus seconds; no trick (every day is 86400 sec)
      mjd1=t
      sec1=(t-mjd1)*86400
c rotation of coordinates, with dragging of velocities
      CALL rotpv('BF  ',.true.,mjd1,sec1,dxbf,dvbf,
     +                 'ECLM',.true.,mjdjpl,secjpl,dx,dv)
* Transformation of velocity from AU/s to AU/d
      DO 2 i=1,3
      dv(i)=dv(i)*86400
 2    CONTINUE

      END
c===================END PUBLIC=============================
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 7, 1997
* ---------------------------------------------------------------------
*
*  ***************************************************************
*  *                                                             *
*  *                          R O T P N                          *
*  *                                                             *
*  *                  General purpose routine                    *
*  *  for transformations between different reference systems    *
*  *                                                             *
*  ***************************************************************
*
*
* INPUT:    RSYS1     -  Starting reference system (EQUM/EQUT/ECLM)
*           EPOCH1    -  Starting epoch (J2000/OFDATE)
*           DATE1     -  Starting date (MJD, TDT): relevant only if
*                        EPOCH1=OFDATE
*           RSYS2     -  Final reference system (EQUM/EQUT/ECLM)
*           EPOCH2    -  Final epoch (J2000/OFDATE)
*           DATE2     -  Final date (MJD, TDT): relevant only if
*                        EPOCH2=OFDATE
*
* OUTPUT:   ROT(3,3)  -  Rotation matrix giving the transformation from
*                        starting to final reference systems:
*                        X2 = ROT X1
*
      SUBROUTINE rotpn(rot,rsys1,epoch1,date1,rsys2,epoch2,date2)
      IMPLICIT NONE

      DOUBLE PRECISION rot(3,3),date1,date2
      CHARACTER*(*) rsys1,rsys2,epoch1,epoch2

      DOUBLE PRECISION eps
      PARAMETER (eps=1.d-6)

      INCLUDE 't2000.h'

      DOUBLE PRECISION r(3,3),date,obl
      INTEGER i,j,nit
      LOGICAL error,epdif
      CHARACTER rsys*4,epoch*10

      DOUBLE PRECISION obleq
      EXTERNAL obleq

      CALL chkref(rsys1,epoch1,error)
      IF(error) THEN
          WRITE(*,110) ' starting',rsys1,epoch1
          STOP '**** rotpn: abnormal end ****'
      END IF
      CALL chkref(rsys2,epoch2,error)
      IF(error) THEN
          WRITE(*,110) ' final',rsys1,epoch1
          STOP '**** rotpn: abnormal end ****'
      END IF
 110  FORMAT(' ERROR: unsupported ',A,' reference system:'/
     +       10X,'RSYS  = ',A/
     +       10X,'EPOCH = ',A)

* Starting point
      rsys=rsys1
      epoch=epoch1
      IF(epoch.EQ.'J2000') THEN
          date=t2000
      ELSEIF(epoch.EQ.'OFDATE') THEN
          date=date1
      ELSE
          STOP '**** rotpn: internal error (01) ****'
      END IF

* Initialization of the rotation matrix (equal to the unit matrix)
      DO 1 i=1,3
      DO 1 j=1,3
 1    rot(i,j)=0.d0
      DO 2 i=1,3
 2    rot(i,i)=1.d0

* Building of the rotation matrix
      nit=0
 3    CONTINUE

* Understand whether the final epoch and the current epoch are the
* same of not
      IF(epoch.EQ.epoch2) THEN
          IF(epoch.EQ.'J2000') THEN
              epdif=.false.
          ELSEIF(epoch.EQ.'OFDATE') THEN
              epdif=(ABS(date-date1).GT.eps)
          ELSE
              STOP '**** rotpn: internal error (02) ****'
          END IF
      ELSE
          epdif=.true.
      END IF

* Different epochs: the transformation have to pass through J2000
* equatorial system
      IF(epdif) THEN
          IF(epoch.NE.'J2000') THEN
              IF(rsys.EQ.'ECLM') THEN
* Transformation ecliptical --> equatorial
                  obl=obleq(date)
                  CALL rotmt(-obl,r,1)
                  CALL mult3(r,rot)
                  rsys='EQUM'
              ELSEIF(rsys.EQ.'EQUT') THEN
* Transformation true equator --> mean equator
                  CALL rnut80(date,r)
                  CALL trsp(r)
                  CALL mult3(r,rot)
                  rsys='EQUM'
              ELSEIF(rsys.EQ.'EQUM') THEN
* Transformation to J2000 (precession)
                  CALL prec(date,r)
                  CALL trsp(r)
                  CALL mult3(r,rot)
                  epoch='J2000'
                  date=t2000
              ELSE
                  STOP '**** rotpn: internal error (03) ****'
              END IF
          ELSE
              IF(rsys.EQ.'ECLM') THEN
* Transformation ecliptical --> equatorial
                  obl=obleq(t2000)
                  CALL rotmt(-obl,r,1)
                  CALL mult3(r,rot)
                  rsys='EQUM'
              ELSEIF(rsys.EQ.'EQUT') THEN
* Transformation true equator --> mean equator
                  CALL rnut80(t2000,r)
                  CALL trsp(r)
                  CALL mult3(r,rot)
                  rsys='EQUM'
              ELSEIF(rsys.EQ.'EQUM') THEN
                  IF(epoch2.EQ.'OFDATE') THEN
                      CALL prec(date2,r)
                      CALL mult3(r,rot)
                      epoch=epoch2
                      date=date2
                  ELSE
                      STOP '**** rotpn: internal error (04) ****'
                  END IF
              ELSE
                  STOP '**** rotpn: internal error (05) ****'
              END IF
          END IF
* Same epoch
      ELSE
          IF(rsys.EQ.rsys2) RETURN
* Transformation of reference system at the same epoch (date)
          IF(rsys.EQ.'EQUT') THEN
* Transformation true equator --> mean equator
              CALL rnut80(date,r)
              CALL trsp(r)
              CALL mult3(r,rot)
              rsys='EQUM'
          ELSEIF(rsys.EQ.'ECLM') THEN
* Transformation ecliptical --> equatorial
              obl=obleq(date)
              CALL rotmt(-obl,r,1)
              CALL mult3(r,rot)
              rsys='EQUM'
          ELSEIF(rsys.EQ.'EQUM') THEN
              IF(rsys2.EQ.'EQUT') THEN
* Transformation mean equator --> true equator
                  CALL rnut80(date,r)
                  CALL mult3(r,rot)
                  rsys='EQUT'
              ELSEIF(rsys2.EQ.'ECLM') THEN
* Transformation equatorial --> ecliptical
                  obl=obleq(date)
                  CALL rotmt(obl,r,1)
                  CALL mult3(r,rot)
                  rsys='ECLM'
              ELSE
                  STOP '**** rotpn: internal error (06) ****'
              END IF
          ELSE
              STOP '**** rotpn: internal error (07) ****'
          END IF
      END IF
      nit=nit+1
      IF(nit.GT.20) STOP '**** rotpn: internal error (08) ****'
      GOTO 3

      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 7, 1997
* ---------------------------------------------------------------------
*
*  ***************************************************************
*  *                                                             *
*  *                          O B L E Q                          *
*  *                                                             *
*  *                 Mean obliquity of ecliptic                  *
*  *            (see Astronomical Almanac 1987, B18)             *
*  *                                                             *
*  ***************************************************************
*
*  INPUT:    TJM   -   Modified Julian Time (TDT)
*
*  OUTPUT:   in radians
*
      double precision function obleq(tjm)
      implicit none

      include 'trig.h'
      include 't2000.h'

      double precision tjm
      double precision ob0,ob1,ob2,ob3,t
      logical first

      save first,ob0,ob1,ob2,ob3

      data first/.true./
      if(first)then
          first=.false.
* IAU value
*         ob0   =  (float(23*3600+26*60)+21.45d0) * radsec
* Improved value
          ob0   =  (float(23*3600+26*60)+21.448d0) * radsec
          ob1   =  -46.815d0  * radsec
          ob2   =  -0.0006d0  * radsec
          ob3   =   0.00181d0 * radsec
      end if
      t     = ( tjm - t2000 ) / 36525.d0
      obleq = (( ob3 * t + ob2 ) * t + ob1 ) * t + ob0
      end
* Author: Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 5, 1996
*
*  ***************************************************************
*  *                                                             *
*  *                         R N U T 8 0                         *
*  *                                                             *
*  *    Nutation matrix according to Wahr (IAU-1980) theory      *
*  *                                                             *
*  ***************************************************************
*
*  INPUT:    TJM         -   Modified Julian Time (TDT)
*
*  OUTPUT:   RNUT(3,3)   -   Nutation matrix, transforming MEAN
*                            coordinates into TRUE coordinates
*                            Xtrue = RNUT Xmean
*
      subroutine rnut80(tjm,rnut)
      implicit none

      include 'trig.h'

      double precision tjm,rnut(3,3)
      double precision r1(3,3),r2(3,3),r3(3,3),rp(3,3)
      double precision epsm,epst,dpsi,deps
      double precision obleq
      external obleq
      integer i,j

      epsm = obleq(tjm)
      call nutn80(tjm,dpsi,deps)
      dpsi = dpsi * radsec
      epst = epsm + deps * radsec
      call rotmt(  epsm , r1 , 1)
      call rotmt( -dpsi , r2 , 3)
      call rotmt( -epst , r3 , 1)
      do 1 i=1,3
      do 1 j=1,3
 1    rp(i,j)=r2(i,1)*r1(1,j)+r2(i,2)*r1(2,j)+r2(i,3)*r1(3,j)
      do 2 i=1,3
      do 2 j=1,3
 2    rnut(i,j)=r3(i,1)*rp(1,j)+r3(i,2)*rp(2,j)+r3(i,3)*rp(3,j)

      end
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 7, 1997
* ---------------------------------------------------------------------
*
*  ***************************************************************
*  *                                                             *
*  *                           P R E C                           *
*  *                                                             *
*  *                      Precession matrix                      *
*  *                                                             *
*  ***************************************************************
*
*  INPUT:    TJM         -   Modified Julian Time (TDT)
*
*  OUTPUT:   RPREC(3,3)  -   Precession matrix, transforming J2000
*                            equatorial coordinates into mean
*                            equatorial coordinates at epoch TJM
*                            Xtjm = RPREC Xj2000
*
      subroutine prec(tjm,rprec)
      implicit none

      include 'trig.h'
      include 't2000.h'

      double precision tjm,rprec(3,3)
      double precision r1(3,3),r2(3,3),r3(3,3),rp(3,3)
      double precision t,zeta,theta,z
      double precision zed,zd,thd,zedd,zdd,thdd,zeddd,zddd,thddd
      integer i,j
      logical first

      save first,zed,zd,thd,zedd,zdd,thdd,zeddd,zddd,thddd
      data first/.true./
*
* Calcolo costanti usate (vedi Astronomical Almanac 1987, B18)
*
      if(first)then
          first=.false.
* Termine lineare
          zed   =   0.6406161d0 * radeg
          zd    =   0.6406161d0 * radeg
          thd   =   0.5567530d0 * radeg
* Termine quadratico
          zedd  =   0.0000839d0 * radeg
          zdd   =   0.0003041d0 * radeg
          thdd  = - 0.0001185d0 * radeg
* Termine cubico
          zeddd =   0.0000050d0 * radeg
          zddd  =   0.0000051d0 * radeg
          thddd = - 0.0000116d0 * radeg
      end if
*
* Calcolo argomenti fondamentali
*
      t     = ( tjm - t2000 ) / 36525.d0
      zeta  = ( ( zeddd * t + zedd ) * t + zed ) * t
      z     = ( (  zddd * t +  zdd ) * t +  zd ) * t
      theta = ( ( thddd * t + thdd ) * t + thd ) * t
*
* Calcolo matrice di rotazione
*
      call rotmt(- zeta , r1 , 3)
      call rotmt( theta , r2 , 2)
      call rotmt(-    z , r3 , 3)
      do 1 i=1,3
      do 1 j=1,3
 1    rp(i,j)=r2(i,1)*r1(1,j)+r2(i,2)*r1(2,j)+r2(i,3)*r1(3,j)
      do 2 i=1,3
      do 2 j=1,3
 2    rprec(i,j)=r3(i,1)*rp(1,j)+r3(i,2)*rp(2,j)+r3(i,3)*rp(3,j)

      end
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 7, 1997
* ---------------------------------------------------------------------
*
*  ***************************************************************
*  *                                                             *
*  *                         N U T N 8 0                         *
*  *                                                             *
*  *    Nutation angles according to Wahr (IAU 1980) theory      *
*  *                                                             *
*  ***************************************************************
*
*  INPUT:    TJM    -   Modified Julian Time (TDT)
*
*  OUTPUT:   DPSI   -   Nutation in longitude (arcsec)
*            DEPS   -   Nutation in obliquity (arcsec)
*
      subroutine nutn80(tjm,dpsi,deps)
      double precision tjm,dpsi,deps,dl,dp,df,dd,dn,rs,p2,t1
      include 't2000.h'
      real l,n
      save
      data rs,p2/4.84813681109536d-6, 6.2831853071795865d0/
      data z1,z2/1.0, 2.0/


      t1=(tjm-t2000)/36525.d0
      t=sngl(t1)
      t2=t*t
      t3=t2*t
*
* Fundamental arguments (IAU 1980)
*
      dl=( 485866.733d0 +1717915922.633d0*t1 +31.310*t2 +0.064*t3)*rs
      dp=(1287099.804d0 + 129596581.224d0*t1 - 0.577*t2 -0.012*t3)*rs
      df=( 335778.877d0 +1739527263.137d0*t1 -13.257*t2 +0.011*t3)*rs
      dd=(1072261.307d0 +1602961601.328d0*t1 - 6.891*t2 +0.019*t3)*rs
      dn=( 450160.280d0 -   6962890.539d0*t1 + 7.455*t2 +0.008*t3)*rs
      l=sngl(dmod(dl,p2))
      p=sngl(dmod(dp,p2))
      x=sngl(dmod(df,p2)*2.d0)
      d=sngl(dmod(dd,p2))
      n=sngl(dmod(dn,p2))
      cl=cos(l)
      sl=sin(l)
      cp=cos(p)
      sp=sin(p)
      cx=cos(x)
      sx=sin(x)
      cd=cos(d)
      sd=sin(d)
      cn=cos(n)
      sn=sin(n)
      cp2=z2*cp*cp -z1
      sp2=z2*sp*cp
      cd2=z2*cd*cd -z1
      sd2=z2*sd*cd
      cn2=z2*cn*cn -z1
      sn2=z2*sn*cn
      cl2=z2*cl*cl -z1
      sl2=z2*sl*cl
      ca=cx*cd2 +sx*sd2
      sa=sx*cd2 -cx*sd2
      cb=ca*cn -sa*sn
      sb=sa*cn +ca*sn
      cc=cb*cn -sb*sn
      sc=sb*cn +cb*sn
      cv=cx*cd2 -sx*sd2
      sv=sx*cd2 +cx*sd2
      ce=cv*cn -sv*sn
      se=sv*cn +cv*sn
      cf=ce*cn -se*sn
      sf=se*cn +ce*sn
      cg=cl*cd2 +sl*sd2
      sg=sl*cd2 -cl*sd2
      ch=cx*cn2 -sx*sn2
      sh=sx*cn2 +cx*sn2
      cj=ch*cl -sh*sl
      sj=sh*cl +ch*sl
      ck=cj*cl -sj*sl
      sk=sj*cl +cj*sl
      cm=cx*cl2 +sx*sl2
      sm=sx*cl2 -cx*sl2
      cq=cl*cd +sl*sd
      sq=sl*cd -cl*sd
      cr=z2*cq*cq -z1
      sr=z2*sq*cq
      cs=cx*cn -sx*sn
      ss=sx*cn +cx*sn
      ct=cs*cl -ss*sl
      st=ss*cl +cs*sl
      cu=cf*cl +sf*sl
      su=sf*cl -cf*sl
      cw=cp*cg -sp*sg
      sw=sp*cg +cp*sg
*
* Series for DPSI
*
      dpsi=-(171996.+174.2*t)*sn +(2062.+0.2*t)*sn2 +46.*(sm*cn+cm*sn)
     1-11.*sm -3.*(sm*cn2+cm*sn2) -3.*(sq*cp-cq*sp) -2.*(sb*cp2-cb*sp2)
     2+(sn*cm-cn*sm) -(13187.+1.6*t)*sc +(1426.-3.4*t)*sp
     3-(517.-1.2*t)*(sc*cp+cc*sp) +(217.-0.5*t)*(sc*cp-cc*sp)
     4+(129.+0.1*t)*sb +48.*sr -22.*sa +(17.-0.1*t)*sp2
     5-15.*(sp*cn+cp*sn) -(16.-0.1*t)*(sc*cp2+cc*sp2) -12.*(sn*cp-cn*sp)
      dpsi=dpsi -6.*(sn*cr-cn*sr) -5.*(sb*cp-cb*sp) +4.*(sr*cn+cr*sn)
     1+4.*(sb*cp+cb*sp) -4.*sq +(sr*cp+cr*sp) +(sn*ca-cn*sa)
     2-(sp*ca-cp*sa) +(sp*cn2+cp*sn2) +(sn*cq-cn*sq) -(sp*ca+cp*sa)
     3-(2274.+0.2*t)*sh +(712.+0.1*t)*sl -(386.+0.4*t)*ss -301.*sj
     4-158.*sg +123.*(sh*cl-ch*sl) +63.*sd2 +(63.+0.1*t)*(sl*cn+cl*sn)
     5-(58.+0.1*t)*(sn*cl-cn*sl) -59.*su -51.*st -38.*sf +29.*sl2
      dpsi=dpsi +29.*(sc*cl+cc*sl) -31.*sk +26.*sx +21.*(ss*cl-cs*sl)
     1+16.*(sn*cg-cn*sg) -13.*(sn*cg+cn*sg) -10.*(se*cl-ce*sl)
     2-7.*(sg*cp+cg*sp) +7.*(sh*cp+ch*sp) -7.*(sh*cp-ch*sp)
     3-8.*(sf*cl+cf*sl) +6.*(sl*cd2+cl*sd2) +6.*(sc*cl2+cc*sl2)
     4-6.*(sn*cd2+cn*sd2) -7.*se +6.*(sb*cl+cb*sl) -5.*(sn*cd2-cn*sd2)
     5+5.*(sl*cp-cl*sp) -5.*(ss*cl2+cs*sl2) -4.*(sp*cd2-cp*sd2)
      dpsi=dpsi +4.*(sl*cx-cl*sx) -4.*sd -3.*(sl*cp+cl*sp)
     1+3.*(sl*cx+cl*sx) -3.*(sj*cp-cj*sp) -3.*(su*cp-cu*sp)
     2-2.*(sn*cl2-cn*sl2) -3.*(sk*cl+ck*sl) -3.*(sf*cp-cf*sp)
     3+2.*(sj*cp+cj*sp) -2.*(sb*cl-cb*sl)
      dpsi=dpsi +2.*(sn*cl2+cn*sl2) -2.*(sl*cn2+cl*sn2)
     1+2.*(sl*cl2+cl*sl2) +2.*(sh*cd+ch*sd) +(sn2*cl-cn2*sl)
     2-(sg*cd2-cg*sd2) +(sf*cl2-cf*sl2) -2.*(su*cd2+cu*sd2)
     3-(sr*cd2-cr*sd2) +(sw*ch+cw*sh) -(sl*ce+cl*se) -(sf*cr-cf*sr)
     4+(su*ca+cu*sa) +(sg*cp-cg*sp) +(sb*cl2+cb*sl2) -(sf*cl2+cf*sl2)
     5-(st*ca-ct*sa) +(sc*cx+cc*sx) +(sj*cr+cj*sr) -(sg*cx+cg*sx)
      dpsi=dpsi +(sp*cs+cp*ss) +(sn*cw-cn*sw) -(sn*cx-cn*sx)
     1-(sh*cd-ch*sd) -(sp*cd2+cp*sd2) -(sl*cv-cl*sv) -(ss*cp-cs*sp)
     2-(sw*cn+cw*sn) -(sl*ca-cl*sa) +(sl2*cd2+cl2*sd2)
     3-(sf*cd2+cf*sd2) +(sp*cd+cp*sd)
*
* Series for DEPS
*
      deps=(92025.+8.9*t)*cn -(895.-0.5*t)*cn2 -24.*(cm*cn-sm*sn)
     1+(cm*cn2-sm*sn2) +(cb*cp2+sb*sp2) +(5736.-3.1*t)*cc
     2+(54.-0.1*t)*cp +(224.-0.6*t)*(cc*cp-sc*sp)
     3-(95.-0.3*t)*(cc*cp+sc*sp) -70.*cb +cr +9.*(cp*cn-sp*sn)
     4+7.*(cc*cp2-sc*sp2) +6.*(cn*cp+sn*sp) +3.*(cn*cr+sn*sr)
     5+3.*(cb*cp+sb*sp) -2.*(cr*cn-sr*sn) -2.*(cb*cp-sb*sp)
      deps=deps +(977.-0.5*t)*ch -7.*cl +200.*cs +(129.-0.1*t)*cj -cg
     1-53.*(ch*cl+sh*sl) -2.*cd2 -33.*(cl*cn-sl*sn) +32.*(cn*cl+sn*sl)
     2+26.*cu +27.*ct +16.*cf -cl2 -12.*(cc*cl-sc*sl) +13.*ck -cx
     3-10.*(cs*cl+ss*sl) -8.*(cn*cg+sn*sg) +7.*(cn*cg-sn*sg)
     4+5.*(ce*cl+se*sl) -3.*(ch*cp-sh*sp) +3.*(ch*cp+sh*sp)
     5+3.*(cf*cl-sf*sl) -3.*(cc*cl2-sc*sl2) +3.*(cn*cd2-sn*sd2) +3.*ce
      deps=deps -3.*(cb*cl-sb*sl) +3.*(cn*cd2+sn*sd2)
     1+3.*(cs*cl2-ss*sl2) +(cj*cp+sj*sp) +(cu*cp+su*sp)
     2+(cn*cl2+sn*sl2) +(ck*cl-sk*sl) +(cf*cp+sf*sp) -(cj*cp-sj*sp)
     3+(cb*cl+sb*sl)
      deps=deps -(cn*cl2-sn*sl2) +(cl*cn2-sl*sn2) -(ch*cd-sh*sd)
     1-(cn2*cl+sn2*sl) -(cf*cl2+sf*sl2) +(cu*cd2-su*sd2) -(cw*ch-sw*sh)
     2+(cl*ce-sl*se) +(cf*cr+sf*sr) -(cb*cl2-sb*sl2)
*
      dpsi=dpsi*1.d-4
      deps=deps*1.d-4
      end
* Author: Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 5, 1996
*
*  ***************************************************************
*  *                                                             *
*  *                         E Q U E Q U                         *
*  *                                                             *
*  *                 Equation of the equinoxes                   *
*  *                                                             *
*  ***************************************************************
*
*
* INPUT:    TJM       -  Modified Julian Time (TDT)
*
* OUTPUT:   EQUEQU    -  Equation of the equinoxes (difference
*                        between apparent sidereal time and mean
*                        sidereal time) in radians
*
      double precision function equequ(tjm)
      implicit none

      include 'trig.h'

      double precision tjm
      double precision oblm,dpsi,deps

      double precision obleq
      external obleq

      oblm = obleq(tjm)
      call nutn80(tjm,dpsi,deps)
      equequ=radsec*dpsi*cos(oblm)

      end
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 7, 1997
* ---------------------------------------------------------------------
*
*  ***************************************************************
*  *                                                             *
*  *                           G M S T                           *
*  *                                                             *
*  *                Greenwich Mean Sidereal Time                 *
*  *                    as a function of UT1                     *
*  *                                                             *
*  ***************************************************************
*
*
* INPUT:    TJM       -  Modified Julian Time (UT1)
*
* OUTPUT:   GMST      -  Greenwich Mean Sidereal Time referred
*                        to the mean equinox of date (rad)
*
      double precision function gmst(tjm)
      implicit none

      include 'trig.h'
      include 't2000.h'

      double precision tjm

      double precision c0,c1,c2,c3,rap
      parameter (c0=24110.54841d0,c1=8640184.812866d0,
     +           c2=9.3104d-2,c3=-6.2d-6)
      parameter (rap=1.00273790934d0)

      integer itjm,i
      double precision t,gmst0,h

* Sidereal time at 0h UT1
      itjm=tjm
      t=(itjm-t2000)/36525.d0
      gmst0=((c3*t+c2)*t+c1)*t+c0
      gmst0=gmst0*dpig/86400.d0

* Increment in GMST from 0h
      h=(tjm-itjm)*dpig
      gmst=gmst0+h*rap
      i=gmst/dpig
      if(gmst.lt.0.d0)i=i-1
      gmst=gmst-i*dpig

      end
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: February 24, 1997
*
*  ***************************************************************
*  *                                                             *
*  *                         O B S C O O                         *
*  *                                                             *
*  *          Body-fixed coordinates of an observatory           *
*  *                                                             *
*  ***************************************************************
*
*
* INPUT:    IDOBS     -  Identifier of the observatory (0-999)
*
* OUTPUT:   XBF(3)    -  Body-fixed position of the observatory (AU)
*           NAME      -  Name of the observatory
*
      subroutine obscoo(idobs,xbf,name)
      implicit none

      include 'trig.h'
c     include 'obscoo.h'
      integer lnobnx
      parameter (lnobnx=47)
c error file
      include 'proout.h'
c INPUT
      integer idobs
c OUTPUT
      double precision xbf(3)
      character*(*) name
c 
      integer ns1,ns2
      parameter (ns1=0,ns2=6200)

      double precision aukm,eradkm,eradau
* Astronomical Unit in km
      parameter (aukm=149.59787066d6)
* Earth radius in km
      parameter (eradkm=6378.140d0)
* Earth radius in AU
      parameter (eradau=eradkm/aukm)

      double precision al1,pxy1,pz1,xbfv(3,ns1:ns2)
      integer unit,i,k
      character*3 ocod
      character*(lnobnx) name1,namev(ns1:ns2)
      character*80 rec
      logical first,loaded(ns1:ns2)
      save first,loaded,xbfv,namev
      data first/.true./

      if(first)then
* Get parallax data from MPC
          do 1 i=ns1,ns2
 1        loaded(i)=.false.

          call filopl(unit,'OBSCODE.dat')

 2        continue
          read(unit,110,end=3)ocod,al1,pxy1,pz1,name1
          call statcode(ocod,k)
          if(k.lt.ns1.or.k.gt.ns2)
     +        stop ' **** obscoo: observatory code out of range ****'
          al1=al1*radeg
          xbfv(1,k)=eradau*pxy1*cos(al1)
          xbfv(2,k)=eradau*pxy1*sin(al1)
          xbfv(3,k)=eradau*pz1
          namev(k)=name1
          loaded(k)=.true.
          goto 2
 3        continue

          call filclo(unit,' ')

************************************
* Get parallax data for radar sites
          call filopl(unit,'RADCODE.dat')
 4        continue
             read(unit,'(a)',end=5)rec
             if(rec(1:1).eq.'!'.or.rec(1:1).eq.' ')goto 4
             read(rec,*)k,al1,pxy1,pz1,name1
             if(k.lt.ns1.or.k.gt.ns2)
     +            stop ' **** obscoo: internal error (11) ****'
             al1=al1*radeg
             xbfv(1,k)=1d-10*pxy1*cos(al1)
             xbfv(2,k)=1d-10*pxy1*sin(al1)
             xbfv(3,k)=1d-10*pz1
             namev(k)=name1
             loaded(k)=.true.
          goto 4
 5        call filclo(unit,' ')
************************************

          first=.false.
      end if
 110  format(a3,f10.5,f8.6,f9.5,a)

      if(loaded(idobs)) then
          do i=1,3
            xbf(i)=xbfv(i,idobs)
          enddo
          name=namev(idobs)
      else
          write(*,101)idobs
 101      format(' obscoo: observatory',i4,
     +       ' is not listed in file "OBSCODE.dat"')
c          stop ' **** obscoo: abnormal end ****'
          write(ierrou,101)idobs
          numerr=numerr+1
          do i=1,3
            xbf(i)=0.d0
          enddo
          name='UNKNOWN'
      end if

      end

* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 7, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         C H K R E F                           *
*  *                                                               *
*  *       Check existence of a reference system indicator         *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    RSYS      -  Reference system
*           EPOCH     -  Epoch
*
* OUTPUT:   ERROR     -  Error flag
*
      SUBROUTINE chkref(rsys,epoch,error)
      IMPLICIT NONE

      CHARACTER*(*) rsys,epoch
      LOGICAL error

      error=.true.

      IF(rsys.EQ.'EQUM') GOTO 1
      IF(rsys.EQ.'EQUT') GOTO 1
      IF(rsys.EQ.'ECLM') GOTO 1
      RETURN

 1    CONTINUE
      IF(epoch.EQ.'J2000') GOTO 2
      IF(epoch.EQ.'OFDATE') GOTO 2
      RETURN

 2    CONTINUE
      error=.false.

      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 7, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         C H 2 R E F                           *
*  *                                                               *
*  *              Translation of a character string                *
*  *              into reference system description                *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    STRING    -  Character string
*
* OUTPUT:   RSYS      -  Reference system type (EQUM/EQUT/ECLM)
*           EPOCH     -  Epoch specification (J2000/OFDATE)
*           ERROR     -  Error flag
*
* WARNING: if EPOCH=OFDATE, the reference system specification must be
*          completed with the date, which is derived from a different
*          source (date of orbital elements or observations)
*
      SUBROUTINE ch2ref(string,rsys,epoch,error)
      IMPLICIT NONE

      CHARACTER*(*) string,rsys,epoch
      LOGICAL error

      CHARACTER cont*20,rest*200,rec*200

      error=.false.

* First item: reference system type
      rec=string
      CALL strcnt(rec,cont,rest,error)
      IF(error) RETURN
      rsys=cont
      CALL upcase(rsys)

* Second item: epoch
      rec=rest
      CALL strcnt(rec,cont,rest,error)
      IF(error) RETURN
      epoch=cont
      CALL upcase(epoch)

      CALL chkref(rsys,epoch,error)

      END
