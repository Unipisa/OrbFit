c Copyright (C) 1998 by Zoran Knezevic (zoran@aob.bg.ac.yu)
c Version: November, 1998
c
c VAISALA
c
c =============================================================
c Subroutine that computes preliminary orbit
c INPUT:
c     TDT     - times of observations (MJD TDT)
c     A       - observed right ascension (RAD)
c     D       - observed declinations (RAD)
c     OBSCOD  - observatory codes
c     IFO     - 1=interactive (fitobs); 2=noninteractive (orbfit)
c
c OUTPUT:
c     ELEK    - elements of the preliminary orbit (KEP)
c     EPOCH   - epoch of the elements 
c     FAIL    - logical control of the succes of computation
c
c =============================================================
      subroutine vaisala(tdt,a,d,obscod,ifo,elek,epoch,fail)
      implicit none
c 
      include 'comvais.h'
      include 'trig.h'
      include 'parobx.h'
c
      integer i,obscod(3),ifo
c 
      double precision a(3),d(3),tdt(3),xobs(3,3)
      double precision epoch,an,ele(6),eleq(6),rot(3,3)
      double precision elek(6),rho,zddif,x2,y2,z2,x2p,y2p,z2p
      double precision eq(6)
c
      logical fail   
c ==============================================================
      fail=.false.
c computation of observer's position
      do i=1,3
        call posobs(tdt(i),obscod(i),1,xobs(1,i))
      enddo
c copy the variables: right ascensions, declinations, MJD(TDT),
c rectangular equatorial topocentric coordinates of the Sun for
c the 3 selected observations
      do i=1,3
        alpha(i)=a(i)
        delta(i)=d(i)
        tjd0(i)=tdt(i)
        xt(i)=-xobs(1,i) 
        yt(i)=-xobs(2,i) 
        zt(i)=-xobs(3,i) 
      enddo
c ==============================================================
c Vaisala
c ==============================================================
      call vaisub(ifo,rho,zddif,x2,y2,z2,x2p,y2p,z2p,fail)
c ==============================================================
c check the succes of computation
      if(fail)then
        return
      endif
c ======================================================
c Orbital elements 
c ======================================================
c copying coord. and veloc. 
      ele(1)=x2
      ele(2)=y2
      ele(3)=z2
      ele(4)=x2p
      ele(5)=y2p
      ele(6)=z2p
      epoch=tjd(2)
      CALL rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0)
      call prodmv(eleq,rot,ele)
      call prodmv(eleq(4),rot,ele(4))
      call coocha(eleq,'CAR',1.d0,eq,'EQU',an)
      IF(an.eq.0.d0)THEN
c orbit as computed is parabolic/hyperbolic
         fail=.true.
         return
      ENDIF
      call coocha(eq,'EQU',1.d0,elek,'KEP',an)

      end
