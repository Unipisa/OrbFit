c Copyright (C) 1998 by Zoran Knezevic (zoran@aob.bg.ac.yu)
c Version: November, 1998
c =============================================================
c
c Subroutine ORBEL
c
c
c Computes preliminary orbital elements from the position
c and velocity at the instant of the second observation
c
c INPUT:
c     x2,y2,z2    -  coordinates at t2 
c     x2p,y2p,z2p -  velocity components at t2
c  
c OUTPUT:
c     ele         -  preliminary orbital elements (KEP)
c     epoch       -  epoch of elements
c     fail        -  logical control on the success of computation    
c ================================================================
c
      subroutine orbel(x2,y2,z2,x2p,y2p,z2p,
     +           ele,epoch,fail)
      implicit none
c
      include 'comvais.h'
      include 'trig.h'
c     
      integer j
c
c auxiliaries used in computations of elements
      double precision x2,y2,z2,x2p,y2p,z2p
      double precision r2,r2s,v2s,rv,ese,ece,p1,p2,an,ecan
      double precision aph1,aph2,beta1,beta2,ses,px,py,pz,qx,qy,qz
      double precision somsi,comsi,sepsi,cepsi,som,com,snod,cnod
      double precision eps,seps,ceps,epoch,ele(6),tmp(6)
c
      logical fail
c    
      double precision obleq
      external obleq
c initialize control of failures in orbit computation
      fail=.false.
c radius and velocity at t2
      r2=sqrt(x2**2+y2**2+z2**2)
      r2s=r2**2
      v2s=x2p**2+y2p**2+z2p**2
      rv=x2*x2p+y2*y2p+z2*z2p
c computation of semimajor axis and eccentricity
      ele(6)=1.d0/(2.d0/r2-v2s)
      ese=rv/sqrt(ele(6))
      ece=1.d0-r2/ele(6)
      ele(5)=sqrt(ese**2+ece**2)
c control against ellipse parameter
      p1=ele(6)*(1-ele(5)**2)
      p2=r2s*v2s-rv**2
      if(abs(p1-p2).gt.1.d-6)then
        write(*,*)'Ellipse parameter problem; dp = ', p1-p2
c        stop
        fail=.true.
        return
      endif
c computation of the mean anomaly for the epoch
      ecan=atan2(ese,ece)
      ele(1)=(ecan-ele(5)*sin(ecan))
      do 4 while(ele(1).ge.2.d0*pig)
 4    ele(1)=ele(1)-2.d0*pig
c computation of vectors P and Q 
      aph1=cos(ecan)/r2
      aph2=-sqrt(ele(6))*sin(ecan)
      beta1=sin(ecan)/r2
      beta2=sqrt(ele(6))*(cos(ecan)-ele(5))
      px=aph1*x2+aph2*x2p
      py=aph1*y2+aph2*y2p
      pz=aph1*z2+aph2*z2p
      ses=sqrt(1.d0-ele(5)**2)
      qx=(beta1*x2+beta2*x2p)/ses
      qy=(beta1*y2+beta2*y2p)/ses
      qz=(beta1*z2+beta2*z2p)/ses
c control for P and Q
      if(abs(px*qx+py*qy+pz*qz).gt.1d-6)then
        write(*,*)'Scalar product (PQ) different from 0'
c        stop 
        fail=.true.
        return
      endif
c mean obliquity of ecliptic
      eps=obleq(tjd(2))
      seps=sin(eps)
      ceps=cos(eps)
      somsi=-py*seps+pz*ceps
      comsi=-qy*seps+qz*ceps
      ele(4)=asin(sqrt(somsi**2+comsi**2))
      ele(2)=atan2(somsi,comsi)
      som=sin(ele(2))
      com=cos(ele(2))
      snod=(py*com-qy*som)/ceps
      cnod=px*com-qx*som
      ele(3)=atan2(snod,cnod)
      do j=1,3
        if(ele(j).lt.0.d0)ele(j)=ele(j)+dpig
      enddo
      do j=1,6
        tmp(7-j)=ele(j)
      enddo
      do j=1,6
        ele(j)=tmp(j)
      enddo       
      epoch=tjd(2)
      return
      end
