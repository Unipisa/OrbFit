c ================(part of) MODULE orb_els===============
c ORBITAL ELEMENTS:
c CONTAINS
c SUBROUTINES 
c   COORDINATE CHANGE
c   coocha      general purpose coordinate/element change (with and w/o der)
c prop2b        2-body propagator, with derivatives (2nd derivatives??)
c *carequ	cartesian to equinoctal transformation
c *equcar	equinoctal to cartesian transformation
c *kepequ	keplerian to equinoctal transformation
c *equkep	equinoctal to keplerian transformation
c *eqpequ	equinoctal polar to equinoctal transformation
c *equeqp	equinoctal to equinoctal polar transformation
c *ekensd	keplerian to equinoctal transformation
c *kepcar       fast conversion directly from keplerian to cartesian
c   INPUT/OUTPUT
c rdoef		input of orbital elements
c *oporbf	open an orbital element file (OEF format)
c *rdorb	read orbital elements from a file (OEF format)
c *clorbf	close an orbital element file (OEF format)
c *fixcnm	generates normal matrix from covariance matrix or viceversa
c rdelem	read orbital els for a list of objects from a list of files
c oefdet	auto-detects format of orbital element files
c rdastb	read orbital elements from Bowell's astorb.dat format
c rdmpca	read orbital elements from MPC format for asteroids
c mpcpds	computes MPC-style packed designation from official IAU code
c iaucod        computes IAU official designations from MPC-style packed(asteroids only)
c outele	verbose output of a set orbital elements to a report file
c wro1lh	writes the header of an orbital element file (1L format)
c wro1lr	writes an orbital element record in an orbital element file (1L format)
c wromlh	writes the header of an orbital element file (ML format)
c wromlr	writes an orbital element record in an orbital element file (ML format)
c
c mpcdat	computes MJD from MPC-style packed dates

c
c
c
c ===================================================================
c COOCHA-COODER
c ===================================================================
c   general purpose coordinate to/from
c   elements change
c
c  specifically designed to handle correctly low
c  eccentricity and low inclination orbits;
c  polar orbits are supported without derivatives
c
c  limitations: only elliptic orbits are handled in this
c  version; singularity for 180 degrees inclination
c
c  coordinate types implememted:
c    'CAR' : cartesian positions and velocities in a vector
c            of length 6
c    'EQU' : equinoctal elements, see carequ, in a vector of
c            length 6; the mean motion enne is also provided
c    'EQP' : equinoctal elements for highly inclined orbit
c    'KEP' : classical keplerian elements,  see kepequ, in a
c             vector of length 6; singular for
c             0 eccentricity, 0 and 180 degrees inclination;
c
c    INPUT:  x(6), coox (string defining the coord. type of x)
c            gm=G x Mass of the Sun (anyway the central mass, e.g.
c                 mass of Sun = mass of planet for massive palnet)
c            cooy (string defining the coord. type)
c
c    OUTPUT: y(6)
c            enne=mean motion; WARNING: enne is not computed
c            for identical transformation
c ====================================================================
c COOCHA
c ====================================================================
c in this simpler version the elements are computed without derivatives
c ===============INTERFACE========================================
      subroutine coocha(x,coox,gm,y,cooy,enne)
      implicit none
c ================= input/output ================
      double precision x(6),y(6),gm,enne
      character*3 coox,cooy
c =============END INTERFACE====================
c ================ workspace ============
      double precision z(6)
c ================ loop indexes ==============
      integer j
c ================ rounding off problem ===========
      double precision roff,rouoff,eps
      integer lflag,nb
**********************************
*  static memory only for:
      save lflag,rouoff
**********************************
c error return with enne=0
      enne=0.d0
c  machine rounding off is computed to decide accuracy
c  required in Kepler equation and controls for singular cases
      data lflag/0/
      if(lflag.eq.0)then
         rouoff=roff(nb)
         lflag=1
      endif
c  dummy case
      if(coox.eq.cooy)then
         do 1 j=1,6
 1         y(j)=x(j)
         return
      endif
c  input is interpreted according to input type coox,
c  and anyway transformed to equinoctal elements z; enne is
c  also computed in any case apart from identical transformation
      if(coox.eq.'EQU')then
         do 2 j=1,6
 2         z(j)=x(j)
         if(x(1).gt.0.d0)then
            enne=sqrt(gm/x(1)**3)
         else
            enne=0.d0
         endif
      elseif(coox.eq.'KEP')then
         call kepequ(x,z)
         if(x(1).gt.0.d0)then
            enne=sqrt(gm/x(1)**3)
         else
            enne=0.d0
         endif
      elseif(coox.eq.'CAR')then
         call carequ(x,gm,z,enne)
         IF(enne.eq.0.d0)RETURN
      elseif(coox.eq.'EQP')then
         call eqpequ(x,z)
         if(x(1).gt.0.d0)then
            enne=sqrt(gm/x(1)**3)
         else
            enne=0.d0
         endif
      else
         write(*,*)'**** coocha: in coord.type ',coox,' unknown ****'
         stop
      endif
c  transformation to the output type cooy
      if(cooy.eq.'EQU')then
         do 3 j=1,6
 3         y(j)=z(j)
      elseif(cooy.eq.'KEP')then
c  this is a potentially singular case; the control for negligible
c  eccentricity and inclination is set to $100 \times$ rounding off
         eps=1.d2*rouoff
         call equkep(z,eps,y)
      elseif(cooy.eq.'EQP')then
c  this case is singular only for zero inclination
         eps=1.d2*rouoff
         call equeqp(z,eps,y)
      elseif(cooy.eq.'CAR')then
c  control for convergence in kepler equation is set to $100 \times$
c  rounding off
         eps=1.d2*rouoff
         call equcar(z,gm,eps,y)
      else
         write(*,*)'**** coocha: out coord.type ',cooy,' unknown ****'
         stop
      endif
      return
      end
c ======================================================
c COODER
c ======================================================
c   general purpose coordinate to/from
c   elements change: version with partial derivatives
c   definitions and input as for coocha (but no equinoctal polar)
c   in output, derpar(6,6) = matrix of partial derivatives
c ==============INTERFACE============================
      subroutine cooder(x,coox,gm,y,cooy,enne,derpar)
      implicit none
c ================= input/output ================
      double precision x(6),y(6),gm,enne,derpar(6,6)
      character*3 coox,cooy
c =============END INTERFACE====================
c ================ workspace ============
      double precision z(6),derws(6,6),derws2(6,6),w(6),ddxde(3,6,6)
c =============== scalars ===============
      double precision t0,det
c ================ loop indexes ==============
      integer j
c ================ derivatives, inversion control =======
      integer ider,ising
c ================ rounding off problem ===========
      double precision roff,rouoff,eps
      integer lflag,nb
**********************************
*  static memory only for:
      save lflag,rouoff
**********************************
c error return with enne=0
      enne=0.d0
c  machine rounding off is computed to decide accuracy
c  required in Kepler equation and controls for singular cases
      data lflag/0/
      if(lflag.eq.0)then
         rouoff=roff(nb)
         lflag=1
      endif
c  dummy case
      if(coox.eq.cooy)then
         do 1 j=1,6
 1         y(j)=x(j)
         call eye(6,derpar)
         return
      endif
c  input is interpreted according to input type coox,
c  and anyway transformed to equinoctal elements z; enne is
c  also computed in any case apart from identical transformation
      if(coox.eq.'EQU')then
         do 2 j=1,6
 2         z(j)=x(j)
         if(x(1).gt.0.d0)then
            enne=sqrt(gm/x(1)**3)
         else
            enne=0.d0
         endif
         call eye(6,derws)
      elseif(coox.eq.'KEP')then
         call ekensd(x,z,derws)
         if(x(1).gt.0.d0)then
            enne=sqrt(gm/x(1)**3)
         else
            enne=0.d0
         endif
      elseif(coox.eq.'CAR')then
         call carequ(x,gm,z,enne)
         IF(enne.eq.0.d0)RETURN
         t0=0.d0
         ider=1
         call prop2b(t0,z,t0,w,gm,ider,derws,ddxde)
         call matin(derws,det,6,0,6,ising,1)
      elseif(coox.eq.'EQP')then
         write(*,*)' partial derivatives for EQP not implemented'
         stop
      else
         write(*,*)'**** coocha: in coord.type ',coox,' unknown ****'
         stop
      endif
c  transformation to the output type cooy
      if(cooy.eq.'EQU')then
         do 3 j=1,6
 3         y(j)=z(j)
         call eye(6,derws2)
      elseif(cooy.eq.'KEP')then
c  this is a potentially singular case; the control for negligible
c  eccentricity and inclination is set to $100 \times$ rounding off
         eps=1.d2*rouoff
         call equkep(z,eps,y)
         call ekensd(y,w,derws2)
* ***    write(*,*)(w(ii)-z(ii),ii=1,6)
         call matin(derws2,det,6,0,6,ising,1)
      elseif(cooy.eq.'EQP')then
         write(*,*)' partial derivatives for EQP not implemented'
         stop
      elseif(cooy.eq.'CAR')then
c  control for convergence in kepler equation is set to $100 \times$
c  rounding off
         t0=0.d0
         ider=1
         call prop2b(t0,z,t0,y,gm,ider,derws2,ddxde)
         eps=1.d2*rouoff
         call equcar(z,gm,eps,w)
* ***    write(*,*)(w(ii)-y(ii),ii=1,6)
      else
         write(*,*)'**** coocha: out coord.type ',cooy,' unknown ****'
         stop
      endif
c  multiplication of jacobian matrices to get the jacobian of the composite
      call mulmat(derws,6,6,derws2,6,6,derpar)
      return
      end
c
c  PROP2B
c
c  Task:       Solves the two body problem in equinoctal coordinates
c
c  Features:   the Newton's method is exploited in the interval
c              [W, 2PI + W] where W is the peri-planet longitude:
c              W = w_p + RAAN
c
c  Input arguments:
c              E         =     equinoctal orbital elements 
c                               e=(a,h,k,p,q,lambda)
c              T0        =     epoch time
c              T1        =     prediction time
c              GM        =     gravitational constant of the system GM
c              IDER       =     flag for derivatives option :
c                                0 = only position and velocites
c                                1 = first derivatives
c                                2 = first and second partial derivatives 
c  Output arguments:
c              X         =     position and velocity components
c                               in absolute cartesian coordinates
c              DXDE      =     first derivatives of position vector x
c                               with respect to elements
c              DDXDE     =     second derivatives of position vector x   
c                               with respect to elements 
c
c****************
c   static memory not required
c****************
      subroutine prop2b(t0,e,t1,x,gm,ider,dxde,ddxde)
      implicit none
      double precision e(6),x(6)
      double precision f(3),g(3),w(3)
      double precision dxde(6,6),ddxde(3,6,6)
c scalars
      double precision t0,t1,gm
      integer ider
c iteration maximum
      integer iter
      parameter (iter=25)
c scalar temporaries
      double precision enne, pml,ecc2,errm,roff,eps,pol,princ
      double precision sinel,cosel,rf,rdf,del,el,beta,r,upq,
     +  xe,ye,coe,xpe,ype,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,
     +  tmp9,tmp10,tmp11,tmp12,dx1de2,dx2de2,dx1de3,dx2de3,dx4de2,
     +  dx5de2,dx4de3,dx5de3
      integer j,n
c  Trgonometric constants 
      include 'trig.h'
c  Mean motion 
      enne=sqrt(gm/e(1)**3)
c  Mean longitude PML at time T1
      pml=e(6)+enne*(t1-t0)
c  Absolute peri-planet longitude POL at epoch T0
      ecc2=e(2)**2+e(3)**2
      errm=roff(n)
      eps=errm*1.d2
      if(ecc2.lt.eps)then
         pol=0.d0
      elseif(ecc2.ge.1.d0)then
         write(*,*)' dcc.ge.1, ecc**2=',ecc2
         write(*,*)e
         stop
      else
         pol=atan2(e(2),e(3))
         pol=princ(pol)
      endif
c  Mean longitude restriction to [POL, POL + 2*PIGR]
c  (Newton's method continuity conditions)
      pml = princ(pml)
      if (pml.lt.pol)then
         pml = pml + dpig
      endif
c  Newton's method for motion equation solution:
c  R(F,lambda) = F - ksinF + hcosF - lambda = 0
c  search for F for a given lambda within [POL, POL + 2PIGR]
      el = pig+pol
      do 70 j=1,iter
             sinel=sin(el)
             cosel=cos(el)
             rf = el - e(3)*sinel + e(2)*cosel - pml
             rdf = 1.d0 - e(3)*cosel - e(2)*sinel
             del = - rf/rdf
             el = el + del
             if (abs(del).lt.eps) goto 100
 70   continue
      write(*,*)' Too many iter. in newton, iter=',iter,' del=',del
      write(*,*) ' eq,eps ',e, eps
      stop
c  Computation of position and velocity on the orbit plane
c  in equinoctal cartesian coordinates (f,g,w)
 100  beta=1.d0/(1.d0+sqrt(1.d0-ecc2))
      xe=e(1)*((1.d0-beta*e(2)**2)*cosel+e(2)*e(3)*beta*sinel-e(3))
      ye=e(1)*((1.d0-beta*e(3)**2)*sinel+e(2)*e(3)*beta*cosel-e(2))
c  Equinoctal reference frame
      upq=1.d0+e(4)**2+e(5)**2
      f(1)=(1.d0-e(4)**2+e(5)**2)/upq
      f(2)=2.d0*e(4)*e(5)/upq
      f(3)=-2.d0*e(4)/upq
      g(1)=2.d0*e(4)*e(5)/upq
      g(2)=(1.d0+e(4)**2-e(5)**2)/upq
      g(3)=2.d0*e(5)/upq
c  Conversion from equinoctal to absolute coordinates
      call lincom(f,xe,g,ye,x)
c  Computation of velocities 
      coe=enne*e(1)**2/sqrt(xe**2+ye**2)
      xpe=coe*(e(2)*e(3)*beta*cosel-(1.d0-beta*e(2)**2)*sinel)
      ype=coe*((1.d0-beta*e(3)**2)*cosel-e(2)*e(3)*beta*sinel)
      call lincom(f,xpe,g,ype,x(4))
c  Computation of partials if required
      if(ider.lt.1)return
c  Equinoctal reference frame, third vector 
      w(1)=2.d0*e(4)/upq
      w(2)=-2.d0*e(5)/upq
      w(3)=(1.d0-e(4)**2-e(5)**2)/upq
c  Computation of same temporary variables
      r=sqrt(xe**2+ye**2)
      tmp1=pml-el
      tmp2=beta+e(2)**2*beta**3/(1.d0-beta)
      tmp3=e(2)*e(3)*beta**3/(1.d0-beta)
      tmp4=beta*e(2)-sinel
      tmp5=beta*e(3)-cosel
      tmp6=beta+e(3)**2*beta**3/(1.d0-beta)
      tmp7=1.d0-r/e(1)
      tmp8=sinel-e(2)
      tmp9=cosel-e(3)
      tmp10=e(1)*cosel/r
      tmp11=e(1)*sinel/r
      tmp12=enne*e(1)**2/r
c  Computation of derivatives of position vector w. r. to the elements 
      dxde(1,1)=(x(1)-3.d0*x(4)*(t1-t0)/2.d0)/e(1)
      dxde(2,1)=(x(2)-3.d0*x(5)*(t1-t0)/2.d0)/e(1)
      dxde(3,1)=(x(3)-3.d0*x(6)*(t1-t0)/2.d0)/e(1)
      dx1de2=-e(1)*(tmp1*tmp2+e(1)*cosel*tmp4/r)          
      dx2de2=e(1)*(tmp1*tmp3-1.d0+e(1)*cosel*tmp5/r)
      call lincom(f,dx1de2,g,dx2de2,dxde(1,2)) 
      dx1de3=-e(1)*(tmp1*tmp3+1.d0-e(1)*sinel*tmp4/r)
      dx2de3=e(1)*(tmp1*tmp6-e(1)*sinel*tmp5/r)
      call lincom(f,dx1de3,g,dx2de3,dxde(1,3))
      dxde(1,4)=2.d0*(e(5)*(ye*f(1)-xe*g(1))-xe*w(1))/upq
      dxde(2,4)=2.d0*(e(5)*(ye*f(2)-xe*g(2))-xe*w(2))/upq
      dxde(3,4)=2.d0*(e(5)*(ye*f(3)-xe*g(3))-xe*w(3))/upq
      dxde(1,5)=2.d0*(e(4)*(-ye*f(1)+xe*g(1))+ye*w(1))/upq          
      dxde(2,5)=2.d0*(e(4)*(-ye*f(2)+xe*g(2))+ye*w(2))/upq
      dxde(3,5)=2.d0*(e(4)*(-ye*f(3)+xe*g(3))+ye*w(3))/upq
      dxde(1,6)=x(4)/enne
      dxde(2,6)=x(5)/enne
      dxde(3,6)=x(6)/enne
c  Computation of derivatives of velocity vector w. r. to the elements
      dxde(4,1)=-(x(4)-3.d0*gm*x(1)*(t1-t0)/r**3)/(2.d0*e(1))
      dxde(5,1)=-(x(5)-3.d0*gm*x(2)*(t1-t0)/r**3)/(2.d0*e(1))
      dxde(6,1)=-(x(6)-3.d0*gm*x(3)*(t1-t0)/r**3)/(2.d0*e(1))
      dx4de2=tmp12*(tmp7*tmp2+e(1)**2*tmp8*tmp4/r**2+tmp10*cosel)          
      dx5de2=-tmp12*(tmp7*tmp3+e(1)**2*tmp8*tmp5/r**2-tmp10*sinel)
      call lincom(f,dx4de2,g,dx5de2,dxde(4,2))
      dx4de3=tmp12*(tmp7*tmp3+e(1)**2*tmp9*tmp4/r**2-tmp11*cosel)
      dx5de3=-tmp12*(tmp7*tmp6+e(1)**2*tmp9*tmp5/r**2+tmp11*sinel)
      call lincom(f,dx4de3,g,dx5de3,dxde(4,3))
      dxde(4,4)=2.d0*(e(5)*(ype*f(1)-xpe*g(1))-xpe*w(1))/upq
      dxde(5,4)=2.d0*(e(5)*(ype*f(2)-xpe*g(2))-xpe*w(2))/upq
      dxde(6,4)=2.d0*(e(5)*(ype*f(3)-xpe*g(3))-xpe*w(3))/upq
      dxde(4,5)=2.d0*(e(4)*(-ype*f(1)+xpe*g(1))+ype*w(1))/upq          
      dxde(5,5)=2.d0*(e(4)*(-ype*f(2)+xpe*g(2))+ype*w(2))/upq
      dxde(6,5)=2.d0*(e(4)*(-ype*f(3)+xpe*g(3))+ype*w(3))/upq
      dxde(4,6)=-enne*e(1)**3*x(1)/r**3
      dxde(5,6)=-enne*e(1)**3*x(2)/r**3
      dxde(6,6)=-enne*e(1)**3*x(3)/r**3
c  Computation of second derivatives if required
      if(ider.lt.2)return
      write(*,*) 'second derivatives not supported in this version'
      stop
      end
c ======================================================
c   {\bf carequ}: coordinate change from cartesian to
c   equinoctal elements
c
c     x= position and velocity, cartesian coordinates
c
c     gm= $G \times Mass$
c
c         eq(1)= a
c
c         eq(2)=ecc*sin(dig)
c
c         eq(3)=ecc*cos(dig) ;  dig=longitude pericentre
c
c         eq(4)=tgim*cos(omega) ; tgim=tg(i/2)
c
c         eq(5)=tgim*sin(omega) ; omega=longitude ascending node
c
c         eq(6)=mean longitude
c
c         enne= mean motion
c
c  if the energy is not negative, eq(6) is set equal to true longitude,
c  eq(1)=negative a and enne=0
c
c  bugs: zero angular momentum and/or energy are handled in such
c  a way that zero divide are avoided, but overflows can occur
c
c  this definition of elements is singular only for a planar
c  retrograde orbit and parabolic orbits
c=======================================================
      subroutine carequ(x,gm,eq,enne)
      implicit double precision (a-h,o-z)
      dimension x(6),eq(6)
      dimension ang(3),f(3),g(3),vlenz(3)
c error return when enne is zero
      enne=0.d0
c  radius and velocity squared
      vel2=prscal(x(4),x(4))
      r=vsize(x)
c  angular momentum
      call prvec(x(1),x(4),ang)
c   zero divide occurs for zero angular momentum
      gei=vsize(ang)
      IF(gei.eq.0.d0)THEN
         write(*,*) '****** carequ: zero angular momentum ******'
         return
      ENDIF
c  angular momentum unit vector, Lenz vector
      call prvec(x(4),ang,vlenz)
      do 1 i=1,3
        vlenz(i)=vlenz(i)/gm-x(i)/r
 1      ang(i)=ang(i)/gei
c   zero divide occurs for inclination of 180 degrees
      d=1.d0+ang(3)
      IF(d.eq.0.d0)THEN
         WRITE(*,*) '****** carequ: 180 deg. inclination ******'
         RETURN
      ENDIF
c  unit vectors of the equinoctal reference system (Broucke and
c  Cefola 1972, CM 5, 303--310) are f, g, ang
      f(1)=1.d0-ang(1)**2/d
      f(2)=-ang(1)*ang(2)/d
      f(3)=-ang(1)
      call prvec(ang,f,g)
c  elements related to eccentricity and inclination
      eq(2)=prscal(vlenz,g)
      eq(3)=prscal(vlenz,f)
      eq(4)=ang(1)/d
      eq(5)=-ang(2)/d
c     tgim1=dsqrt(eq(4)**2+eq(5)**2)
c  test on energy
      ainv=2/r-vel2/gm
      if(ainv.eq.0.d0)then
c eq(1) is q
         eq(1)=(r+prscal(vlenz,x))/(1+vsize(vlenz))
         enne=0
         cosf=prscal(x,f)
         sinf=prscal(x,g)
         eq(6)=atan2(sinf,cosf)
         eq(6)=princ(eq(6))
         return
c         stop '****** carequ: parabolic orbit ******'
      elseif(ainv.lt.0.d0)then
         eq(1)=1.d0/ainv
         enne=0
         cosf=prscal(x,f)
         sinf=prscal(x,g)
         eq(6)=atan2(sinf,cosf)
         eq(6)=princ(eq(6))
         return
      endif
c   semimajor axis and mean motion
      eq(1)=1.d0/ainv
      enne=dsqrt(gm/eq(1)**3)
c   mean longitude from non--singular Kepler equation
      ecc2=eq(2)**2+eq(3)**2
      rad=dsqrt(1.d0-ecc2)
      beta=1.d0/(1.d0+rad)
      chk=eq(2)*eq(3)*beta
      ch=1.d0-eq(2)**2*beta
      ck=1.d0-eq(3)**2*beta
      x2=prscal(x,f)
      y2=prscal(x,g)
      cosf=eq(3)+(ck*x2-chk*y2)/(eq(1)*rad)
      sinf=eq(2)+(ch*y2-chk*x2)/(eq(1)*rad)
c   eccentric longitude
      fe=datan2(sinf,cosf)
      eq(6)=fe+eq(2)*cosf-eq(3)*sinf
c   reduction to principal value
      eq(6)=princ(eq(6))
      return
      end
c ======================================================
c   {\bf equcar}: coordinate change from
c   equinoctal elements to cartesian
c
c     x= position and veloctiy, cartesian coordinates
c
c     gm= $G \times Mass$
c
c     eps= convergence control for non--singular Kepler equation
c
c     eq: equinoctal elements, see carequ
c
c         enne= mean motion
c
c  bugs: if  eq(1) is negative, the hyperbolic case is not handled in
c  this version
c
c  this definition of elements is singular only for a planar
c  retrograde orbit and parabolic orbits
c=======================================================
      subroutine equcar(eq,gm,eps,x)
c  ====================================================================
      implicit double precision (a-h,o-z)
      dimension x(6),eq(6),f(3),g(3)
      include 'trig.h'
c   test for hyperbolic orbit
      if(eq(1).le.0.d0)then
         stop'****** equcar: hyperbolic/parabolic orbit ******'
      endif
c  non--singular intermediate variables
      ecc2=eq(2)**2+eq(3)**2
      rad=dsqrt(1.d0-ecc2)
      beta=1.d0/(1.d0+rad)
      chk=eq(2)*eq(3)*beta
      ch=1.d0-eq(2)**2*beta
      ck=1.d0-eq(3)**2*beta
      tgim2=eq(4)**2+eq(5)**2
c     tgim=dsqrt(tgim2)
      opwz=1.d0+tgim2
c   mean motion
      enne=dsqrt(gm/eq(1)**3)
c  unit vectors of the equinoctal reference system (Broucke and
c  Cefola 1972, CM 5, 303--310) are f, g and the angular momentum
c  unit vector
      f(1)=(1.d0-eq(4)**2+eq(5)**2)/opwz
      f(2)=2*eq(4)*eq(5)/opwz
      f(3)=-2*eq(4)/opwz
      g(1)=2*eq(4)*eq(5)/opwz
      g(2)=(1.d0+eq(4)**2-eq(5)**2)/opwz
      g(3)=2*eq(5)/opwz
c  Non singular Kepler equation
      ecc=dsqrt(ecc2)
      if(ecc.lt.eps)then
c  for negligible eccentricity, the eccentric longitude is
c  set equal to the mean longitude
         fe=eq(6)
         cosf=cos(fe)
         sinf=sin(fe)
      else
c  mean longitude is reduced to the interval (dig, dig+2*pig)
c  with dig=longitude of pericentre
         tlong=eq(6)
         dig=atan2(eq(2),eq(3))
         tlong=princ(tlong-dig)+dig
c  initial condition for Newton's method is always
c  the longitude of apocentre; this ensures the right
c  convexity for secure convergence
         fe=dig+pig
         do 16 i=1,100
           cosf=cos(fe)
           sinf=sin(fe)
           df=(fe-tlong+eq(2)*cosf-eq(3)*sinf)/
     +        (1.d0-eq(2)*sinf-eq(3)*cosf)
           if(dabs(df).lt.eps)goto 17
 16        fe=fe-df
c  convergence problems -- this should happen only for
c  extremely high eccentricity
         stop'****** equcar: 100 iterations of Newton ******'
      endif
c  cartesian coordinates in the equinoctal frame x2, y2, 0
 17   x2=eq(1)*(ch*cosf+chk*sinf-eq(3))
      y2=eq(1)*(ck*sinf+chk*cosf-eq(2))
      do 18 i=1,3
 18     x(i)=x2*f(i)+y2*g(i)
c  cartesian velocities in the equinoctal frame xp2, yp2, 0
      de=enne*eq(1)**2/dsqrt(x2**2+y2**2)
      xp2=de*(chk*cosf-ch*sinf)
      yp2=de*(ck*cosf-chk*sinf)
      do 19 i=1,3
 19     x(i+3)=xp2*f(i)+yp2*g(i)
      return
      end
c=======================================================
c   {\bf kepequ}: coordinate change from keplerian to
c   equinoctal elements
c     eq: equinoctal elements, see carequ
c
c     el(1)=a
c
c     el(2)=ecc
c
c     el(3)=inclination (radians)
c
c     el(4)=longitude asc. node (radians)
c
c     el(5)=argument of pericentre (radians)
c
c     el(6)=mean anomaly (radians)
c
c   bugs: no test for hyperbolic/parabolic case; element a
c   is copied whatever its meaning
c=======================================================
      subroutine kepequ(el,eq)
      implicit double precision (a-h,o-z)
      dimension el(6),eq(6)
      eq(1)=el(1)
      dig=el(4)+el(5)
      ecc=el(2)
      eq(2)=ecc*sin(dig)
      eq(3)=ecc*cos(dig)
      tgim=tan(el(3)/2.d0)
      eq(4)=tgim*dsin(el(4))
      eq(5)=tgim*dcos(el(4))
      eq(6)=dig+el(6)
      eq(6)=princ(eq(6))
      return
      end
c=======================================================
c   {\bf equkep}: coordinate change from equinoctal elements
c   to keplerian
c
c     eq: equinoctal elements, see equcar/carequ
c
c     el: keplerian elements, see kepequ
c
c     eps=control on eccentricity/inclination; for smaller values,
c     the angles eq(4), eq(5) are set to arbitrary value 0
c
c   bugs: no test for hyperbolic/parabolic case; element a
c   is copied whatever its meaning
c=======================================================
      subroutine equkep(eq,eps,el)
      implicit double precision (a-h,o-z)
      dimension el(6),eq(6)
      el(1)=eq(1)
c  test on eccentricity
      ecc=sqrt(eq(2)**2+eq(3)**2)
      if(ecc.lt.eps)then
         dig=0.d0
      else
         dig=atan2(eq(2),eq(3))
      endif
      el(2)=ecc
c   test on tangent of half inclination
      tgi2=sqrt(eq(4)**2+eq(5)**2)
      if(tgi2.lt.eps)then
         el(4)=0.d0
      else
         el(4)=atan2(eq(4),eq(5))
      endif
      el(3)=2.d0*atan(tgi2)
c   angular variables
      el(5)=dig-el(4)
      el(6)=eq(6)-dig
      el(4)=princ(el(4))
      el(5)=princ(el(5))
      el(6)=princ(el(6))
      return
      end
c=======================================================
c   {\bf eqpequ}: coordinate change from equinoctal 
c   polar to equinoctal elements
c     eq: equinoctal elements, see carequ
c
c     el(1)=a
c
c     el(2)=e sin (omega)
c
c     el(3)=e cos (omega)
c
c     el(4)=tg (I/2) sin (Omega)
c
c     el(5)=tg (I/2) cos (Omega)
c
c     el(6)=mean argument of latitude (=mean anom+ omega)
c
c       omega= argument of pericentre
c       Omega= longitude of asc. node
c
c   limitations: no test for hyperbolic/parabolic case; element a
c   is copied whatever its meaning
c=======================================================
      subroutine eqpequ(el,eq)
      implicit double precision (a-h,o-z)
      dimension el(6),eq(6)
      eq(1)=el(1)
      omnod=atan2(el(4),el(5))
      co=cos(omnod)
      so=sin(omnod)
      eq(2)=el(2)*so+el(3)*co
      eq(3)=el(3)*co-el(2)*so
      eq(4)=el(4)
      eq(5)=el(5)
      eq(6)=princ(el(6)+omnod)
      return
      end
c=======================================================
c   {\bf equeqp}: coordinate change from equinoctal elements
c   to equinoctal polar
c
c     eq: equinoctal elements, see equcar/carequ
c
c     el: equinoctal polar elements, see eqpequ
c
c     eps=control on inclination; for smaller values,
c       the computation stops 
c   limitation: no test for hyperbolic/parabolic case; element a
c   is copied whatever its meaning
c=======================================================
      subroutine equeqp(eq,eps,el)
      implicit double precision (a-h,o-z)
      dimension el(6),eq(6)
      el(1)=eq(1)
c  test on inclination
      tgi2=sqrt(eq(4)**2+eq(5)**2)
      if(tgi2.lt.eps)then
         write(*,*)' inclination zero cannot be handled in eqp'
         stop
      else
         omnod=atan2(eq(4),eq(5))
      endif
      co=cos(omnod)
      so=sin(omnod)
      el(2)=eq(2)*co-eq(3)*so
      el(3)=eq(3)*co+eq(2)*so
      el(4)=eq(4)
      el(5)=eq(5)
      el(6)=princ(eq(6)-omnod)
      return
      end
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: February 24, 1997
*
* Transformation between keplerian and non-singular (equinoctal) elements,
* with jacobian matrix
* (all angles in radians)
*
* INPUT:    elem      -  keplerian elements
*                            elem(1) = semimajor axis
*                            elem(2) = eccentricity
*                            elem(3) = inclination
*                            elem(4) = longitude of ascending node
*                            elem(5) = argument of pericenter
*                            elem(6) = mean anomaly
*
* OUTPUT:   ekns      -  non-singular elements
*                            ekns(1) = semimajor axis
*                            ekns(2) = h = e*sin(varpi)
*                            ekns(3) = k = e*cos(varpi)
*                            ekns(4) = P = tan(i/2)*sin(Omega)
*                            ekns(5) = Q = tan(i/2)*cos(Omega)
*                            ekns(6) = mean longitude
*           dkdns     -  jacobian matrix d(ekns)/d(elem)
*
      subroutine ekensd(elem,ekns,dkdns)
      implicit none

      double precision elem(6),ekns(6),dkdns(6,6)
      double precision alp,ti2,sinalp,cosalp,dti2,sinom,cosom
      integer i,k

      include 'trig.h'

      ekns(1)=elem(1)
* Longitude of pericenter
      alp=elem(5)+elem(4)
      sinalp=sin(alp)
      cosalp=cos(alp)
      ekns(2)=elem(2)*sinalp
      ekns(3)=elem(2)*cosalp
      ti2=tan(elem(3)/2)
      sinom=sin(elem(4))
      cosom=cos(elem(4))
      ekns(4)=ti2*sinom
      ekns(5)=ti2*cosom
      ekns(6)=alp+elem(6)
      ekns(6)=mod(ekns(6),dpig)

      do 1 i=1,6
      do 1 k=1,6
 1    dkdns(i,k)=0.d0

      dkdns(1,1)=1.d0
      dkdns(2,2)=sinalp
      dkdns(2,5)=ekns(3)
      dkdns(2,4)=ekns(3)
      dkdns(3,2)=cosalp
      dkdns(3,5)=-ekns(2)
      dkdns(3,4)=-ekns(2)
      dti2=1.d0/(2*cos(elem(3)/2)**2)
      dkdns(4,3)=dti2*sinom
      dkdns(4,4)=ekns(5)
      dkdns(5,3)=dti2*cosom
      dkdns(5,4)=-ekns(4)
      dkdns(6,4)=1
      dkdns(6,5)=1
      dkdns(6,6)=1

      end
c  ====================================================================
c  KEPCAR - fast coordinate change
c  ====================================================================
      SUBROUTINE kepcar(ek,gm,ivel,x)
c  ====================================================================
      IMPLICIT NONE
c INPUT
c keplerian elements, mass of central body
      DOUBLE PRECISION ek(6),gm
c ivel=1 also velocities, ivel=0 no
      INTEGER ivel
c OUTPUT
      DOUBLE PRECISION x(6)
c END INTERFACE
      DOUBLE PRECISION f(3),g(3),ecc,h,k,rad,beta,dig,chk,ch,ck,p,q,ecc2
      DOUBLE PRECISION tgim2,tlong,fe,sinf,cosf,princ,x2,y2,xp2,yp2,de
      DOUBLE PRECISION opwz,enne,df
c controls: convergence of kepler's equation
      DOUBLE PRECISION eps, roff
      INTEGER lflag,nb,i
      DATA lflag/0/
      SAVE lflag, eps
      include 'trig.h'
c ======================================================================
c initialize control of rounding off
      if(lflag.eq.0)then
         eps=roff(nb)*1.d2
         lflag=1
      endif
c   test for hyperbolic orbit
      if(ek(1).le.0.d0)then
         WRITE(*,*)'****** kepcar: hyperbolic/parabolic orbit ******'
      endif
c  non--singular intermediate variables
      ecc=ek(2)
      ecc2=ecc**2
      rad=dsqrt(1.d0-ecc2)
      beta=1.d0/(1.d0+rad)
      dig=ek(4)+ek(5)
      k=ecc*cos(dig)
      h=ecc*sin(dig)
      chk=h*k*beta
      ch=1.d0-h**2*beta
      ck=1.d0-k**2*beta
      tgim2=tan(ek(3)/2.d0)
c     tgim=dsqrt(tgim2)
      opwz=1.d0+tgim2
      p=tgim2*sin(ek(4))
      q=tgim2*cos(ek(4))
c   mean motion
      enne=dsqrt(gm/ek(1)**3)
c  unit vectors of the equinoctal reference system (Broucke and
c  Cefola 1972, CM 5, 303--310) are f, g and the angular momentum
c  unit vector
      f(1)=(1.d0-p**2+q**2)/opwz
      f(2)=2*p*q/opwz
      f(3)=-2*p/opwz
      g(1)=2*p*q/opwz
      g(2)=(1.d0+p**2-q**2)/opwz
      g(3)=2*q/opwz
c  Non singular Kepler equation
      IF(ecc.lt.eps)THEN
c  for negligible eccentricity, the eccentric longitude is
c  set equal to the mean longitude
         fe=princ(dig+ek(6))
         cosf=cos(fe)
         sinf=sin(fe)
      ELSE
c  mean longitude is reduced to the interval (dig, dig+2*pig)
c  with dig=longitude of pericentre
         tlong=dig+ek(6)
         tlong=princ(tlong-dig)+dig
c  initial condition for Newton's method is always
c  the longitude of apocentre; this ensures the right
c  convexity for secure convergence
         fe=dig+pig
         DO i=1,10
            cosf=cos(fe)
            sinf=sin(fe)
            df=(fe-tlong+h*cosf-k*sinf)/
     +           (1.d0-h*sinf-k*cosf)
            fe=fe-df
            IF(dabs(df).lt.eps)GOTO 17
         ENDDO
c  convergence problems -- this should happen only for
c  extremely high eccentricity
         WRITE(*,*)'****** kepcar: 10 iterations of Newton ******'
      ENDIF
c  cartesian coordinates in the equinoctal frame x2, y2, 0
 17   x2=ek(1)*(ch*cosf+chk*sinf-k)
      y2=ek(1)*(ck*sinf+chk*cosf-h)
      DO i=1,3
        x(i)=x2*f(i)+y2*g(i)
      ENDDO
      IF(ivel.eq.0)RETURN
c  cartesian velocities in the equinoctal frame xp2, yp2, 0
      de=enne*ek(1)**2/dsqrt(x2**2+y2**2)
      xp2=de*(chk*cosf-ch*sinf)
      yp2=de*(ck*cosf-chk*sinf)
      DO  i=1,3
        x(i+3)=xp2*f(i)+yp2*g(i)
      ENDDO
      return
      end
* Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it),
*                            Andrea Milani (milani@dm.unipi.it),
*                            Zoran Knezevic (zoran@aob.aob.bg.ac.yu)
* Version: February 12, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D E L E M                           *
*  *                                                               *
*  *          Read orbital elements for a list of objects          *
*  *                     from a list of files                      *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Output FORTRAN unit (report file)
*           OBJNAM    -  Object names
*           NOBJ      -  Number of objects
*           INFILS    -  Input files
*           NFIL      -  Number of input files
*           DEFORB    -  Is the orbit already defined?
*
* OUTPUT:   DEFORB    -  Was the orbit found?
*           DEFCN     -  Tells whether covariance/normal matrices
*                            are defined
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
*           ELEM      -  Orbital elements (ECLM J2000)
*           COVE      -  Covariance matrix of orbital elements
*           NORE      -  Normal matrix of orbital elements
*           MASS      -  Mass (solar masses)
*           HMAG      -  H absolute magnitude (if <-100, missing)
*           GMAG      -  G slope parameter
*           COMELE    -  Comment on orbital elements
*
* WARNING: the routine is designed to allow multiple calls on
*          different lists of files/objects: for this reason,
*          orbital elements of objects which have already an
*          orbit defined (DEFORB=.true.) ARE NOT modified.
*          For this reason, the user must define suitably DEFORB
*          BEFORE calling RDELEM
*
      SUBROUTINE rdelem(unit,objnam,nobj,infils,nfil,deforb,defcn,
     +                  eltype,telem,elem,cove,nore,
     +                  mass,h,g,comele)
      IMPLICIT NONE
      INTEGER nfilx
      PARAMETER(nfilx=20)
      INTEGER unit,nobj,nfil
      DOUBLE PRECISION telem(nobj),elem(6,nobj),cove(6,6,nobj)
      DOUBLE PRECISION nore(6,6,nobj),mass(nobj),h(nobj),g(nobj)
      CHARACTER*(*) objnam(nobj),infils(nfil),eltype(nobj),comele(nobj)
      LOGICAL deforb(nobj),defcn(nobj)

      INTEGER i,lf,is1,lfo,uniin,n
      CHARACTER form*10,infil1*100
      LOGICAL opened

      INTEGER lench
      EXTERNAL lench

* Number of input files
      IF(nfil.GT.nfilx) STOP '**** rdelem: nfil > nfilx ****'

* Nothing to do
      IF(nfil.LE.0) RETURN

* Loop on files
      DO 2 i=1,nfil
      lf=lench(infils(i))
      opened=.false.

* Understand file format
      is1=index(infils(i)(1:lf),'[')
      IF(is1.GT.0 .AND. infils(i)(lf:lf).EQ.']') THEN
          form=infils(i)(is1+1:lf-1)
          lf=is1-1
          infil1=infils(i)(1:lf)
      ELSE
          infil1=infils(i)
          form=' '
          CALL filopn(uniin,infil1,'OLD')
          opened=.true.
          CALL oefdet(uniin,infil1,form)
          REWIND(uniin)
      END IF
      lfo=lench(form)
      IF(form.NE.' ') WRITE(*,102) infil1(1:lf),form(1:lfo)
 102  FORMAT('Scanning file "',A,'" (format: ',A,')')

* Reading file
      IF(form.EQ.'OEF') THEN
          IF(opened) CALL filclo(uniin,' ')
          CALL rdoef(infil1,objnam,nobj,deforb,defcn,eltype,telem,
     +               elem,cove,nore,mass,h,g,comele)
      ELSEIF(form.EQ.'BA1') THEN
          IF(.NOT.opened) CALL filopn(uniin,infil1,'OLD')
          CALL rdast1(uniin,infil1,objnam,nobj,deforb,defcn,
     +                eltype,telem,elem,cove,nore,mass,h,g,comele)
          CALL filclo(uniin,' ')
      ELSEIF(form.EQ.'BA2') THEN
          IF(.NOT.opened) CALL filopn(uniin,infil1,'OLD')
          CALL rdast2(uniin,infil1,objnam,nobj,deforb,defcn,
     +                eltype,telem,elem,cove,nore,mass,h,g,comele)
          CALL filclo(uniin,' ')
      ELSEIF(form.EQ.'MPC-A') THEN
          IF(.NOT.opened) CALL filopn(uniin,infil1,'OLD')
          CALL rdmpca(uniin,infil1,objnam,nobj,deforb,defcn,
     +                eltype,telem,elem,cove,nore,mass,h,g,comele)
          CALL filclo(uniin,' ')
      ELSEIF(form.EQ.' ') THEN
          WRITE(*,120) infil1(1:lf)
          STOP '**** rdelem: abnormal end ****'
      ELSE
          WRITE(*,121) form(1:lfo),infil1(1:lf)
          STOP '**** rdelem: abnormal end ****'
      END IF
 120  FORMAT('ERROR: unknown format type for file "',A,'"')
 121  FORMAT('ERROR: unsupported format "',A,'" for file "',A,'"')
c set default for magnitude data
      DO n=1,nobj
        IF(h(n).lt.-1.d6)THEN
c leave it            
        ENDIF
        IF(g(n).lt.-1.d6)THEN
           g(n)=0.15d0
        ENDIF
      ENDDO
* End of loop on files
 2    CONTINUE

 30   CONTINUE

      END



* Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 21, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         W R O 1 L H                           *
*  *                                                               *
*  *       Writes the header of an orbital element file            *
*  *   (single-line format, different epochs, keplerian elements)  *
*  *                                                               *
*  *****************************************************************
*
* OUTPUT:   UNIT      -  Output FORTRAN unit
*           RSYS      -  Reference system type (EQUM/EQUT/ECLM)
*           EPOCH     -  Reference system epoch (J2000/OFDATE)
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*
      SUBROUTINE wro1lh(unit,rsys,epoch,eltype)
      IMPLICIT NONE

      INTEGER unit
      CHARACTER*(*) rsys,epoch,eltype

      INCLUDE 'parcmc.h'

      INTEGER l1,l2,l3,nb
      CHARACTER*100 bl

      INTEGER lench
      EXTERNAL lench

      l1=lench(rsys)
      l2=lench(epoch)
      l3=lench(eltype)
      nb=MAX(14-l1-l2,1)
      bl=' '

      WRITE(unit,100) comcha,comcha,eltype(1:l3),comcha,
     +                rsys(1:l1),epoch(1:l2),
     +                bl(1:nb),comcha
 100  FORMAT('format  = ''OEF1.1''       ',A,' file format'/
     +       'rectype = ''1L''           ',A,' record type (1L/ML)'/
     +       'elem    = ''',A,'''          ',A,
     +                   ' type of orbital elements'/
     +       'refsys  = ',A,1X,A,A,A,' default reference system'/
     +       'END_OF_HEADER')

      IF(eltype.EQ.'KEP') THEN
          WRITE(unit,201) comcha
      ELSEIF(eltype.EQ.'CAR') THEN
          WRITE(unit,202) comcha
      ELSEIF(eltype.EQ.'EQU') THEN
          WRITE(unit,203) comcha
      END IF

 201  FORMAT(A,' Name, Epoch(MJD), a, e, i, long. node,',
     +         ' arg. peric., mean anomaly')
 202  FORMAT(A,' Name, Epoch(MJD), cartesian position and velocity',
     +         ' vectors')
 203  FORMAT(A,' Name, Epoch(MJD), a, e*sin(LP), e*cos(LP),',
     +         ' tan(i/2)*sin(LN), tan(i/2)*cos(LN), mean long.')

      END
* Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 21, 1998
* Version December 18, 1998 Steven Chesley (chesley@dm.unipi.it)
* Added ELTYPE, Changed time field to F13.4 to handle jpleph406
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         W R O 1 L R                           *
*  *                                                               *
*  *  Writes an orbital element record in an orbital element file  *
*  *   (single-line format, different epochs, keplerian elements)  *
*  *                     (multi-line format)                       *
*  *                                                               *
*  *****************************************************************
*
* OUTPUT:   UNIT      -  Output FORTRAN unit
*           NAME      -  Name of planet/asteroid/comet
*           ELEM(6)   -  Orbital element vector (keplerian elements)
*           ELTYPE    -  Type of orbital elements (KEP/EQU/CAR)
*           T0        -  Epoch of orbital elements (MJD, TDT)
*           H         -  H absolute magnitude (if <-100, missing)
*           G         -  G slope parameter
*
* WARNING: the routine does not write the header of the file: this
*          must be generated by calling subroutine wro1lh
*
      SUBROUTINE wro1lr(unit,name,elem,eltype,t0,h,g)
      IMPLICIT NONE

      INTEGER unit
      DOUBLE PRECISION elem(6),t0,h,g
      CHARACTER*(*) name,eltype

      INCLUDE 'trig.h'

* Expected max length of name
      INTEGER namtl
      PARAMETER (namtl=12)
      DOUBLE PRECISION cnvele(6)

      INTEGER ln,nb,i
      CHARACTER*(namtl) blanks

      INTEGER lench
      EXTERNAL lench
      CHARACTER*80 elefmt
* Name
      ln=MAX(1,lench(name))
      nb=MAX(1,namtl-ln)
      blanks=' '
* Convert to degrees
      DO i=1,6
         cnvele(i)=elem(i)
      ENDDO
      IF(eltype .eq. 'CAR')THEN
         elefmt='('''''''',A,'''''''',A,1x,F13.6,'//
     +        '6(1x,F17.13),'//
     +        '2F6.2)'
      ELSEIF(eltype .eq. 'KEP')THEN
         DO i=3,6
            cnvele(i)=cnvele(i)*degrad
         ENDDO
         elefmt='('''''''',A,'''''''',A,F13.6,'//
     +        '1p,6(e25.16),'//
     +        '0p,2F6.2)'
c        elefmt='('''''''',A,'''''''',A,F13.6,'//
c    +        'F16.12,1x,F12.10,4(1x,F12.7),'//
c    +        '2F6.2)'
      ELSEIF(eltype .eq. 'EQU')THEN
         cnvele(6)=cnvele(6)*degrad
         elefmt='('''''''',A,'''''''',A,F13.6,'//
     +        'F16.12,4(1x,f12.9),1x,F12.7,'//
     +        '2F6.2)'
      ELSE
         STOP '*** wro1lr internal error' 
      ENDIF
* Write
      IF(h.GT.-100.d0) THEN
          WRITE(unit,elefmt) name(1:ln),blanks(1:nb),t0,cnvele,h,g
      ELSE
          WRITE(unit,elefmt) name(1:ln),blanks(1:nb),t0,cnvele
      ENDIF

      END






* Copyright (C) 1997-2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 7, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         O U T E L E                           *
*  *                                                               *
*  *           Verbose output of a set orbital elements            *
*  *                       to a report file                        *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Output FORTRAN unit (report file)
*           ELEM      -  Orbital elements (ECLM J2000)
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           T0        -  Epoch of orbital elements (MJD, TDT)
*           LABEL     -  Label
*           MULTI     -  Multi-line output
*           STDOUT    -  Standard output
*
      SUBROUTINE outele(unit,elem,eltype,t0,label,multi,stdout)
      IMPLICIT NONE

      INTEGER unit
      DOUBLE PRECISION elem(6),t0
      CHARACTER*(*) eltype,label
      LOGICAL multi,stdout

      INCLUDE 'trig.h'

      INTEGER lt,i,day,month,year,ll
      CHARACTER cm*3
      DOUBLE PRECISION hour

      INTEGER lench
      CHARACTER*3 chmon
      EXTERNAL lench,chmon

      ll=lench(label)

      IF(multi .AND. (ll.GT.0)) THEN
          IF(unit.GT.0) WRITE(unit,133) label(1:ll)
          IF(stdout) WRITE(*,133) label(1:ll)
      END IF
      IF(eltype.EQ.'KEP') THEN
          IF(multi) THEN
              IF(unit.GT.0) WRITE(unit,100) elem(1),elem(2),
     +                                      (elem(i)*degrad,i=3,6)
              IF(stdout) WRITE(*,100) elem(1),elem(2),
     +                                (elem(i)*degrad,i=3,6)
          ELSE
              IF(ll.GT.0) THEN
                  IF(unit.GT.0) WRITE(unit,120) label(1:ll),elem(1),
     +                                          elem(2),
     +                                         (elem(i)*degrad,i=3,6),t0
                  IF(stdout) WRITE(*,120) label(1:ll),elem(1),elem(2),
     +                                    (elem(i)*degrad,i=3,6),t0
              ELSE
                  IF(unit.GT.0) WRITE(unit,130) elem(1),elem(2),
     +                                         (elem(i)*degrad,i=3,6),t0
                  IF(stdout) WRITE(unit,130) elem(1),elem(2),
     +                                       (elem(i)*degrad,i=3,6),t0
              END IF
          END IF
      ELSEIF(eltype.EQ.'EQU') THEN
          IF(multi) THEN
              IF(unit.GT.0) WRITE(unit,101) (elem(i),i=1,5),
     +                                      elem(6)*degrad
              IF(stdout) WRITE(*,101) (elem(i),i=1,5),elem(6)*degrad
          ELSE
              IF(ll.GT.0) THEN
                  IF(unit.GT.0) WRITE(unit,121) label(1:ll),
     +                                          (elem(i),i=1,5),
     +                                          elem(6)*degrad,t0
                  IF(stdout) WRITE(*,121) label(1:ll),(elem(i),i=1,5),
     +                                    elem(6)*degrad,t0
              ELSE
                  IF(unit.GT.0) WRITE(unit,131) (elem(i),i=1,5),
     +                                          elem(6)*degrad,t0
                  IF(stdout) WRITE(*,131) (elem(i),i=1,5),
     +                                    elem(6)*degrad,t0
              END IF
          END IF
      ELSEIF(eltype.EQ.'CAR') THEN
          IF(multi) THEN
              IF(unit.GT.0) WRITE(unit,102) elem
              IF(stdout) WRITE(*,102) elem
          ELSE
              IF(ll.GT.0) THEN
                  IF(unit.GT.0) WRITE(unit,122) label(1:ll),elem,t0
                  IF(stdout) WRITE(*,122) label(1:ll),elem,t0
              ELSE
                  IF(unit.GT.0) WRITE(unit,132) elem,t0
                  IF(stdout) WRITE(*,132) elem,t0
              END IF
          END IF
      ELSE
          lt=lench(eltype)
          WRITE(*,200) eltype(1:lt)
          STOP '**** outele: abnormal end ****'
      END IF
 120  FORMAT(8X,'KepElem(',A,'):',1P,E15.7,0P,F13.8,4F10.5,
     +           ' (T=',F10.3,')')
 130  FORMAT(8X,'KepElem:',1P,E15.7,0P,F13.8,4F10.5,
     +           ' (T=',F10.3,')')
 121  FORMAT(8X,'EQUElem(',A,'):',1P,E15.7,0P,4F13.8,F10.5,
     +           ' (T=',F10.3,')')
 131  FORMAT(8X,'EQUElem:',1P,E15.7,0P,4F13.8,F10.5,' (T=',F10.3,')')
 122  FORMAT(8X,'PosVel(',A,'):',1P,6E15.7,' (T=',F10.3,')')
 132  FORMAT(8X,'PosVel:',1P,6E15.7,' (T=',F10.3,')')
 100  FORMAT(8X,'Semimajor axis     =',1P,E23.14,0P,' AU'/
     +       8X,'Eccentricity       =',F20.15/
     +       8X,'Inclination        =',F18.13,' deg'/
     +       8X,'Long. of node      =',F18.13,' deg'/
     +       8X,'Arg. of pericenter =',F18.13,' deg'/
     +       8X,'Mean anomaly       =',F18.13,' deg')
 101  FORMAT(8X,'Semimajor axis     =',1P,E23.14,0P,' AU'/
     +       8X,'h [e*sin(w)]       =',F20.15/
     +       8X,'k [e*cos(w)]       =',F20.15/
     +       8X,'P [tg(i/2)*sin(N)] =',F20.15/
     +       8X,'Q [tg(i/2)*cos(N)] =',F20.15/
     +       8X,'Mean longitude     =',F18.13,' deg')
 102  FORMAT(8X,'Position vector   =',1X,1P,3E22.14,' AU'/
     +       8X,'Velocity vector   =',1X,3E22.14,' AU/d')
 200  FORMAT('ERROR: unknown type "',A,'" of orbital elements')
 133  FORMAT(8X,'Orbital elements for ',A,':')

      IF(multi) THEN
          CALL mjddat(t0,day,month,year,hour)
          cm=chmon(month)
          WRITE(unit,110) t0,cm,day,year,hour
          IF(stdout) WRITE(*,110) t0,cm,day,year,hour
      END IF
 110  FORMAT(8X,'Epoch of elements  : MJD',F17.8,' TDT (',A,I3,',',
     +       I5,',',F10.6,' h)')

      END
* Copyright (C) 1997-2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 13, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         O U T E L 1                           *
*  *                                                               *
*  *           Verbose output of a set orbital elements            *
*  *                       to a report file                        *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Output FORTRAN unit (report file)
*           ELEM      -  Orbital elements (ECLM J2000)
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           T0        -  Epoch of orbital elements (MJD, TDT)
*           LABEL     -  Label
*           MULTI     -  Multi-line output
*           STDOUT    -  Standard output
*           W         -  Indentation width
*
      SUBROUTINE outel1(unit,elem,eltype,t0,label,multi,stdout,w)
      IMPLICIT NONE

      INTEGER unit,w
      DOUBLE PRECISION elem(6),t0
      CHARACTER*(*) eltype,label
      LOGICAL multi,stdout

      INCLUDE 'trig.h'

      INTEGER lt,i,day,month,year,ll
      CHARACTER cm*3
      DOUBLE PRECISION hour

      CHARACTER*50 b
      DATA b/'                                                  '/
      SAVE b

      INTEGER lench
      CHARACTER*3 chmon
      EXTERNAL lench,chmon

      IF(w.LE.0) STOP '**** outel1: internal error (01) ****'
      IF(w.GT.LEN(b)) STOP '**** outel1: internal error (02) ****'

      ll=lench(label)

      IF(multi .AND. (ll.GT.0)) THEN
          IF(unit.GT.0) WRITE(unit,133) b(1:w),label(1:ll)
          IF(stdout) WRITE(*,133) b(1:w),label(1:ll)
      END IF
      IF(eltype.EQ.'KEP') THEN
          IF(multi) THEN
              IF(unit.GT.0) WRITE(unit,100)
     +            b(1:w),elem(1),b(1:w),elem(2),
     +            (b(1:w),elem(i)*degrad,i=3,6)
              IF(stdout) WRITE(*,100)
     +            b(1:w),elem(1),b(1:w),elem(2),
     +            (b(1:w),elem(i)*degrad,i=3,6)
          ELSE
              IF(ll.GT.0) THEN
                  IF(unit.GT.0) WRITE(unit,120)
     +                b(1:w),label(1:ll),elem(1),elem(2),
     +                (elem(i)*degrad,i=3,6),t0
                  IF(stdout) WRITE(*,120)
     +                 b(1:w),label(1:ll),elem(1),elem(2),
     +                (elem(i)*degrad,i=3,6),t0
              ELSE
                  IF(unit.GT.0) WRITE(unit,130)
     +                b(1:w),elem(1),elem(2),(elem(i)*degrad,i=3,6),t0
                  IF(stdout) WRITE(unit,130)
     +                b(1:w),elem(1),elem(2),(elem(i)*degrad,i=3,6),t0
              END IF
          END IF
      ELSEIF(eltype.EQ.'EQU') THEN
          IF(multi) THEN
              IF(unit.GT.0) WRITE(unit,101)
     +            (b(1:w),elem(i),i=1,5),b(1:w),elem(6)*degrad
              IF(stdout) WRITE(*,101)
     +            (b(1:w),elem(i),i=1,5),b(1:w),elem(6)*degrad
          ELSE
              IF(ll.GT.0) THEN
                  IF(unit.GT.0) WRITE(unit,121)
     +                b(1:w),label(1:ll),(elem(i),i=1,5),
     +                elem(6)*degrad,t0
                  IF(stdout) WRITE(*,121)
     +                b(1:w),label(1:ll),(elem(i),i=1,5),
     +                elem(6)*degrad,t0
              ELSE
                  IF(unit.GT.0) WRITE(unit,131)
     +                b(1:w),(elem(i),i=1,5),elem(6)*degrad,t0
                  IF(stdout) WRITE(*,131)
     +                b(1:w),(elem(i),i=1,5),elem(6)*degrad,t0
              END IF
          END IF
      ELSEIF(eltype.EQ.'CAR') THEN
          IF(multi) THEN
              IF(unit.GT.0) WRITE(unit,102)
     +            b(1:w),(elem(i),i=1,3),b(1:w),(elem(i),i=4,6)
              IF(stdout) WRITE(*,102)
     +            b(1:w),(elem(i),i=1,3),b(1:w),(elem(i),i=4,6)
          ELSE
              IF(ll.GT.0) THEN
                  IF(unit.GT.0) WRITE(unit,122)
     +                b(1:w),label(1:ll),elem,t0
                  IF(stdout) WRITE(*,122)
     +                b(1:w),label(1:ll),elem,t0
              ELSE
                  IF(unit.GT.0) WRITE(unit,132) b(1:w),elem,t0
                  IF(stdout) WRITE(*,132) b(1:w),elem,t0
              END IF
          END IF
      ELSE
          lt=lench(eltype)
          WRITE(*,200) eltype(1:lt)
          STOP '**** outel1: abnormal end ****'
      END IF
 120  FORMAT(A,'KepElem(',A,'):',1P,E15.7,0P,F13.8,4F10.5,
     +           ' (T=',F10.3,')')
 130  FORMAT(A,'KepElem:',1P,E15.7,0P,F13.8,4F10.5,
     +           ' (T=',F10.3,')')
 121  FORMAT(A,'EQUElem(',A,'):',1P,E15.7,0P,4F13.8,F10.5,
     +           ' (T=',F10.3,')')
 131  FORMAT(A,'EQUElem:',1P,E15.7,0P,4F13.8,F10.5,' (T=',F10.3,')')
 122  FORMAT(A,'PosVel(',A,'):',1P,6E15.7,' (T=',F10.3,')')
 132  FORMAT(A,'PosVel:',1P,6E15.7,' (T=',F10.3,')')
 100  FORMAT(A,'Semimajor axis     =',1P,E23.14,0P,' AU'/
     +       A,'Eccentricity       =',F20.15/
     +       A,'Inclination        =',F18.13,' deg'/
     +       A,'Long. of node      =',F18.13,' deg'/
     +       A,'Arg. of pericenter =',F18.13,' deg'/
     +       A,'Mean anomaly       =',F18.13,' deg')
 101  FORMAT(A,'Semimajor axis     =',1P,E23.14,0P,' AU'/
     +       A,'h [e*sin(w)]       =',F20.15/
     +       A,'k [e*cos(w)]       =',F20.15/
     +       A,'P [tg(i/2)*sin(N)] =',F20.15/
     +       A,'Q [tg(i/2)*cos(N)] =',F20.15/
     +       A,'Mean longitude     =',F18.13,' deg')
 102  FORMAT(A,'Position vector   =',1X,1P,3E22.14,' AU'/
     +       A,'Velocity vector   =',1X,3E22.14,' AU/d')
 200  FORMAT('ERROR: unknown type "',A,'" of orbital elements')
 133  FORMAT(A,'Orbital elements for ',A,':')

      IF(multi) THEN
          CALL mjddat(t0,day,month,year,hour)
          cm=chmon(month)
          WRITE(unit,110) b(1:w),t0,cm,day,year,hour
          IF(stdout) WRITE(*,110) b(1:w),t0,cm,day,year,hour
      END IF
 110  FORMAT(A,'Epoch of elements  : MJD',F17.8,' TDT (',A,I3,',',
     +       I5,',',F10.6,' h)')

      END
* Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 21, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         W R O M L H                           *
*  *                                                               *
*  *       Writes the header of an orbital element file            *
*  *                     (multi-line format)                       *
*  *                                                               *
*  *****************************************************************
*
* OUTPUT:   UNIT      -  Output FORTRAN unit
*           RSYS      -  Reference system type (EQUM/EQUT/ECLM)
*           EPOCH     -  Reference system epoch (J2000/OFDATE)
*
      SUBROUTINE wromlh(unit,rsys,epoch)
      IMPLICIT NONE

      INTEGER unit
      CHARACTER*(*) rsys,epoch

      INCLUDE 'parcmc.h'

      INTEGER l1,l2,nb
      CHARACTER*100 bl

      INTEGER lench
      EXTERNAL lench

      l1=lench(rsys)
      l2=lench(epoch)
      nb=MAX(14-l1-l2,1)
      bl=' '

      WRITE(unit,100) comcha,comcha,rsys(1:l1),epoch(1:l2),
     +                bl(1:nb),comcha
 100  FORMAT('format  = ''OEF1.1''       ',A,' file format'/
     +       'rectype = ''ML''           ',A,' record type (1L/ML)'/
     +       'refsys  = ',A,1X,A,A,A,' default reference system'/
     +       'END_OF_HEADER')

      END
* Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Modified by Andrea Milani, vers. 1.8.3, January 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         W R O M L R                           *
*  *                                                               *
*  *  Writes an orbital element record in an orbital element file  *
*  *                     (multi-line format)                       *
*  *                                                               *
*  *****************************************************************
*
* OUTPUT:   UNIT      -  Output FORTRAN unit
*           NAME      -  Name of planet/asteroid/comet
*           ELEM(6)   -  Orbital element vector
*           ELTYPE    -  Type of orbital elements (KEP/EQU/CAR)
*           T0        -  Epoch of orbital elements (MJD, TDT)
*           COVE      -  Covariance matrix of orbital elements
*           DEFCOV    -  Tells whether the covariance matrix is defined
*           NORE      -  Normal matrix of orbital elements
*           DEFNOR    -  Tells whether the normal matrix is defined
*           H         -  H absolute magnitude (if <-100, missing)
*           G         -  G slope parameter
*           MASS      -  Mass (solar masses)
*
* WARNING: the routine does not write the header of the file: this
*          must be generated by calling subroutine wromlh
*
      SUBROUTINE wromlr(unit,name,elem,eltype,t0,cove,defcov,
     +                  nore,defnor,h,g,mass)
      IMPLICIT NONE

      INTEGER unit
      DOUBLE PRECISION elem(6),t0,h,g,cove(6,6),nore(6,6),mass
      CHARACTER*(*) name,eltype
      LOGICAL defcov,defnor

      INCLUDE 'trig.h'
      INCLUDE 'parcmc.h'

      INTEGER l1,ln,i,k
      DOUBLE PRECISION cnv(6),std(6)

      INTEGER lench
      EXTERNAL lench
c eigenvalues, eigenvectors
      DOUBLE PRECISION eigvec(6,6),eigval(6),fv1(6),fv2(6)
      INTEGER ierr
* Name
      ln=lench(name)
      IF(ln.LE.0) THEN
          name='????'
          ln=4
      END IF
      l1=1
      IF(name(l1:l1).EQ.' ') THEN
          DO 1 l1=1,ln
          IF(name(l1:l1).NE.' ') GOTO 2
 1        CONTINUE
 2        CONTINUE
          IF(l1.GT.ln) THEN
              name='????'
              ln=4
          END IF
      END IF
      WRITE(unit,100) name(l1:ln)
 100  FORMAT(A)

* Orbital elements
      DO 3 i=1,6
      cnv(i)=1
 3    CONTINUE
      IF(eltype.EQ.'KEP') THEN
          DO 4 i=3,6
          cnv(i)=degrad
 4        CONTINUE
          if(elem(6).lt.0.d0)elem(6)=elem(6)+dpig
          if(elem(5).lt.0.d0)elem(5)=elem(5)+dpig
          if(elem(4).lt.0.d0)elem(4)=elem(4)+dpig
          WRITE(unit,201) comcha
          WRITE(unit,101) (elem(i)*cnv(i),i=1,6)
      ELSEIF(eltype.EQ.'CAR') THEN
          WRITE(unit,202) comcha
          WRITE(unit,102) elem
      ELSEIF(eltype.EQ.'EQU') THEN
          cnv(6)=degrad
          IF(elem(6).lt.0.d0)THEN
             WRITE(*,*) ' wromlr: negative mean anomaly', elem(6)*cnv(6)
             elem(6)=elem(6)+dpig
          ENDIF
          WRITE(unit,203) comcha
          WRITE(unit,103) (elem(i)*cnv(i),i=1,6)
      ELSE
          STOP '**** wromlr: unsupported orbital element type ****'
      END IF
 101  FORMAT(' KEP ',1P,E22.14,0P,F18.15,4F18.13)
 102  FORMAT(' CAR ',1P,6E22.14)
 103  FORMAT(' EQU ',1P,E22.14,0P,4F19.15,F18.13)
 201  FORMAT(A,' Keplerian elements: a, e, i, long. node,',
     +         ' arg. peric., mean anomaly')
 202  FORMAT(A,' Cartesian position and velocity vectors')
 203  FORMAT(A,' Equinoctial elements: a, e*sin(LP), e*cos(LP),',
     +         ' tan(i/2)*sin(LN), tan(i/2)*cos(LN), mean long.')

* Epoch
      WRITE(unit,104) t0
 104  FORMAT(' MJD ',F18.8,' TDT')

* Mass
      IF(mass.NE.0.d0) WRITE(unit,105) mass
 105  FORMAT(' MAS ',1P,E20.12)

* Magnitudes
      IF(h.GT.-100.d0) WRITE(unit,106) h,g
 106  FORMAT(' MAG ',2F7.3)

* Covariance matrix
      IF(defcov) THEN
          DO 5 i=1,6
          std(i)=SQRT(cove(i,i))
 5        CONTINUE
c eigenvalues
          CALL rs(6,6,cove,eigval,1,eigvec,fv1,fv2,ierr)

          DO i=1,6
            IF(eigval(i).gt.0.d0)THEN
               eigval(i)=sqrt(eigval(i))
            ELSE
               WRITE(*,*)'wromlr: zero/negative eigenvalue', eigval(i)
               WRITE(*,*) '  for asteroid ',name(l1:ln)
               eigval(i)=-sqrt(-eigval(i))
            ENDIF
          ENDDO
c RMS, eigenvalues and weak direction are so far commented
          WRITE(unit,107) comcha,(std(i)*cnv(i),i=1,6)
 107      FORMAT(A1,' RMS ',1P,6E14.5)
          WRITE(unit,111) comcha,eigval
 111      FORMAT(A1,' EIG',1P,6E14.5)
          WRITE(unit,110) comcha,(eigvec(i,6),i=1,6)
 110      FORMAT(A1,' WEA',6F10.5)
c covariance matrix is given uncommented, to be readable
          WRITE(unit,108) ((cove(i,k)*cnv(i)*cnv(k),k=i,6),i=1,6)
 108      FORMAT(' COV ',1P,3E23.15)
      END IF

* Normal matrix
      IF(defnor) THEN
          WRITE(unit,109) ((nore(i,k)/(cnv(i)*cnv(k)),k=i,6),i=1,6)
      END IF
 109  FORMAT(' NOR ',1P,3E23.15)

      END



* Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: February 12, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D A S T 1                           *
*  *                                                               *
*  *          Read orbital elements for a list of objects          *
*  *       from a file written in Bowell's astorb.dat format       *
*  *          (version reading OLD, pre-1999 format)               *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Input FORTRAN unit (must be already opened)
*           FILNAM    -  Input file name (for error messages)
*           OBJNAM    -  Object names
*           NOBJ      -  Number of objects
*
* OUTPUT:   DEFORB    -  Tells whether orbital elements are defined
*           DEFCN     -  Tells whether covariance/normal matrices
*                            are defined
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
*           ELEM      -  Orbital elements (ECLM J2000)
*           COVE      -  Covariance matrix of orbital elements
*           NORE      -  Normal matrix of orbital elements
*           MASS      -  Mass (solar masses)
*           H         -  H absolute magnitude (if <-100, missing)
*           G         -  G slope parameter
*           COMELE    -  Comment on orbital elements
*
* WARNING: the routine assumes that objects having DEFORB=.true.
*          have already orbital elements defined (possibly from another
*          input file) and does not overwrite them
*
* OBJECT NAME TRANSLATION: all names of objects appearing in the
* file are modified (removing all blanks) before comparison with
* the name requested by the calling module
*
      SUBROUTINE rdast1(unit,filnam,objnam,nobj,deforb,defcn,
     +                  eltype,telem,elem,cove,nore,mass,h,g,comele)
      IMPLICIT NONE

      INTEGER unit,nobj
      DOUBLE PRECISION telem(nobj),elem(6,nobj),cove(6,6,nobj)
      DOUBLE PRECISION nore(6,6,nobj),mass(nobj),h(nobj),g(nobj)
      CHARACTER*(*) filnam,objnam(nobj),eltype(nobj),comele(nobj)
      LOGICAL deforb(nobj),defcn(nobj)

      INCLUDE 'trig.h'

      INTEGER ln,nr,lf,nrem,k,flags(6),year,month,day,lc
      DOUBLE PRECISION el1(6)
      CHARACTER n1*4,n2*18,name*18,hc*5,gc*5,krc*10

      INTEGER lench
      DOUBLE PRECISION tjm1
      EXTERNAL lench,tjm1

* Number of remaining object (orbit not yet found)
      nrem=0
      DO 10 k=1,nobj
      IF(deforb(k)) GOTO 10
      nrem=nrem+1
 10   CONTINUE
      IF(nrem.LE.0) RETURN
      lf=lench(filnam)

      nr=0
 1    CONTINUE
      READ(unit,100,END=2) n1,n2
      nr=nr+1
      IF(n1.EQ.'    ') THEN
          name=n2
      ELSE
          name=n1
      END IF
      CALL rmsp(name,ln)
      IF(ln.LE.0) THEN
          WRITE(*,200) filnam(1:lf),nr
          GOTO 1
      END IF
 200  FORMAT('ERROR in reading file "',A,'": no object name at record',
     +       I6)

      DO 3 k=1,nobj
      IF(deforb(k)) GOTO 3
      IF(name.EQ.objnam(k)) THEN
          BACKSPACE(unit)
          READ(unit,100) n1,n2,hc,gc,flags,year,month,day,el1

          deforb(k)=.true.
          defcn(k)=.false.
          eltype(k)='KEP'
          telem(k)=tjm1(day,month,year,0.D0)
          elem(1,k)=el1(6)
          elem(2,k)=el1(5)
          elem(3,k)=el1(4)*radeg
          elem(4,k)=el1(3)*radeg
          elem(5,k)=el1(2)*radeg
          elem(6,k)=el1(1)*radeg
          mass(k)=0.d0
          IF(hc.EQ.'     ') THEN
              h(k)=-1.D9
          ELSE
              READ(hc,101) h(k)
          END IF
          IF(gc.EQ.'     ') THEN
              g(k)=0.15D0
          ELSE
              READ(gc,101) g(k)
          END IF
          WRITE(krc,107) nr
          CALL rmsp(krc,lc)
          comele(k)='read from file "'//filnam(1:lf)//
     +              '" at record '//krc(1:lc)

          nrem=nrem-1
          IF(nrem.LE.0) RETURN
      END IF
 100  FORMAT(A4,1X,A18,17X,A5,1X,A5,17X,6I4,12X,I4,2I2,1X,
     +       3(F10.6,1X),F9.6,1X,F10.8,1X,F12.8)
 101  FORMAT(F5.2)
 107  FORMAT(I6)

 3    CONTINUE

      GOTO 1
 2    CONTINUE

      END
* Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: February 12, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D A S T 2                           *
*  *                                                               *
*  *          Read orbital elements for a list of objects          *
*  *       from a file written in Bowell's astorb.dat format       *
*  *          (version reading NEW, post-1999 format)              *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Input FORTRAN unit (must be already opened)
*           FILNAM    -  Input file name (for error messages)
*           OBJNAM    -  Object names
*           NOBJ      -  Number of objects
*
* OUTPUT:   DEFORB    -  Tells whether orbital elements are defined
*           DEFCN     -  Tells whether covariance/normal matrices
*                            are defined
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
*           ELEM      -  Orbital elements (ECLM J2000)
*           COVE      -  Covariance matrix of orbital elements
*           NORE      -  Normal matrix of orbital elements
*           MASS      -  Mass (solar masses)
*           H         -  H absolute magnitude (if <-100, missing)
*           G         -  G slope parameter
*           COMELE    -  Comment on orbital elements
*
* WARNING: the routine assumes that objects having DEFORB=.true.
*          have already orbital elements defined (possibly from another
*          input file) and does not overwrite them
*
* OBJECT NAME TRANSLATION: all names of objects appearing in the
* file are modified (removing all blanks) before comparison with
* the name requested by the calling module
*
      SUBROUTINE rdast2(unit,filnam,objnam,nobj,deforb,defcn,
     +                  eltype,telem,elem,cove,nore,mass,h,g,comele)
      IMPLICIT NONE

      INTEGER unit,nobj
      DOUBLE PRECISION telem(nobj),elem(6,nobj),cove(6,6,nobj)
      DOUBLE PRECISION nore(6,6,nobj),mass(nobj),h(nobj),g(nobj)
      CHARACTER*(*) filnam,objnam(nobj),eltype(nobj),comele(nobj)
      LOGICAL deforb(nobj),defcn(nobj)

      INCLUDE 'trig.h'

      INTEGER ln,nr,lf,nrem,k,flags(6),year,month,day,lc
      DOUBLE PRECISION el1(6)
      CHARACTER n1*5,n2*18,name*18,hc*5,gc*5,krc*10

      INTEGER lench
      DOUBLE PRECISION tjm1
      EXTERNAL lench,tjm1

* Number of remaining object (orbit not yet found)
      nrem=0
      DO 10 k=1,nobj
      IF(deforb(k)) GOTO 10
      nrem=nrem+1
 10   CONTINUE
      IF(nrem.LE.0) RETURN
      lf=lench(filnam)

      nr=0
 1    CONTINUE
      READ(unit,100,END=2) n1,n2
      nr=nr+1
      IF(n1.EQ.'     ') THEN
          name=n2
      ELSE
          name=n1
      END IF
      CALL rmsp(name,ln)
      IF(ln.LE.0) THEN
          WRITE(*,200) filnam(1:lf),nr
          GOTO 1
      END IF
 200  FORMAT('ERROR in reading file "',A,'": no object name at record',
     +       I6)

      DO 3 k=1,nobj
      IF(deforb(k)) GOTO 3
      IF(name.EQ.objnam(k)) THEN
          BACKSPACE(unit)
          READ(unit,100) n1,n2,hc,gc,flags,year,month,day,el1

          deforb(k)=.true.
          defcn(k)=.false.
          eltype(k)='KEP'
          telem(k)=tjm1(day,month,year,0.D0)
          elem(1,k)=el1(6)
          elem(2,k)=el1(5)
          elem(3,k)=el1(4)*radeg
          elem(4,k)=el1(3)*radeg
          elem(5,k)=el1(2)*radeg
          elem(6,k)=el1(1)*radeg
          mass(k)=0.d0
          IF(hc.EQ.'     ') THEN
              h(k)=-1.D9
          ELSE
              READ(hc,101) h(k)
          END IF
          IF(gc.EQ.'     ') THEN
              g(k)=0.15D0
          ELSE
              READ(gc,101) g(k)
          END IF
          WRITE(krc,107) nr
          CALL rmsp(krc,lc)
          comele(k)='read from file "'//filnam(1:lf)//
     +              '" at record '//krc(1:lc)

          nrem=nrem-1
          IF(nrem.LE.0) RETURN
      END IF
 100  FORMAT(A5,1X,A18,17X,A5,1X,A5,17X,6I4,12X,I4,2I2,1X,
     +       3(F10.6,1X),F9.6,1X,F10.8,1X,F12.8)
 101  FORMAT(F5.2)
 107  FORMAT(I6)

 3    CONTINUE

      GOTO 1
 2    CONTINUE

      END
* Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
*                            Andrea Milani (milani@dm.unipi.it)
*
* Version: February 12, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         O E F D E T                           *
*  *                                                               *
*  *         Auto-detects format of orbital element files          *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Input FORTRAN unit (must be already opened)
*           FILNAM    -  Input file name (for error messages)
*
* OUTPUT:   FORM      -  Format:
*                          1) OEF:   ORBFIT orbital element file
*                          2) BA1:   Bowell's astrob.dat (pre-1999)
*                          3) BA2:   Bowell's astrob.dat (post-1999)
*                          4) MPC-A: MPC (asteroids)
*                          5) BAC:   Bowell private format with C
* in the future:
*                          6) MPC-C: MPC (comets)
*
      SUBROUTINE oefdet(unit,filnam,form)
      IMPLICIT NONE

      INTEGER unit
      CHARACTER*(*) filnam,form

      INCLUDE 'parcmc.h'

* Number of supported formats
      INTEGER nfx
      PARAMETER (nfx=5)

      INTEGER lf,lr,i,k,npf,ipf,lrwb
      LOGICAL poss(nfx),neoh,error
      CHARACTER rec*300,recwb*200,b4*4,p4*4
      CHARACTER*50 tmp,tmp1,tmp2

      INTEGER lench
      EXTERNAL lench

* Names of supported formats
      CHARACTER*10 fname(nfx)
      DATA fname/'OEF','BA1','BA2','MPC-A','BAC'/

      form=' '

* METHOD: at the beginning, we flag as possible all the supported
* formats. Then we scan the first records of the file and discard
* formats as we find lines which are not compatible with them.
* At the end, detection is successful if one and only one format
* remains.
      DO 1 i=1,nfx
      poss(i)=.true.
 1    CONTINUE

      neoh=.true.
* Scan only first 100 records of the file
      DO 2 k=1,100
      READ(unit,100,END=3) rec
 100  FORMAT(A)
      lr=lench(rec)

* Format 1 (OEF)
      IF(poss(1).AND.neoh) THEN
          recwb=rec
          CALL rmsp(recwb,lrwb)
          IF(lrwb.GT.0) THEN
              IF(recwb.EQ.'END_OF_HEADER') THEN
                  neoh=.false.
              ELSEIF(recwb(1:1).NE.comcha) THEN
                  IF(index(recwb,'=').EQ.0) THEN
                      poss(1)=.false.
                  ELSEIF(recwb(1:7).EQ.'format=') THEN
                      tmp=recwb(8:)
                      i=index(tmp,comcha)
                      IF(i.GT.0) THEN
                          tmp1=tmp(1:i-1)
                          tmp=tmp1
                      END IF
                      CALL strcnt(tmp,tmp1,tmp2,error)
                      IF(error) THEN
                          poss(1)=.false.
                      ELSEIF(tmp1.NE.'OEF1.1') THEN
                          poss(1)=.false.
                      END IF
                  END IF
              END IF
          END IF
      END IF

* Format 2 (BA1)
      IF(poss(2)) THEN
          IF(lr.NE.265) poss(2)=.false.
          b4(1:1)=rec(5:5)
          b4(2:2)=rec(24:24)
          b4(3:3)=rec(40:40)
          b4(4:4)=rec(46:46)
          IF(b4.NE.'    ') poss(2)=.false.
          p4(1:1)=rec(117:117)
          p4(2:2)=rec(128:128)
          p4(3:3)=rec(139:139)
          p4(4:4)=rec(149:149)
          IF(p4.NE.'....') poss(2)=.false.
      END IF

* Format 3 (BA2)
      IF(poss(3)) THEN
          IF(lr.NE.266) poss(3)=.false.
          b4(1:1)=rec(6:6)
          b4(2:2)=rec(25:25)
          b4(3:3)=rec(41:41)
          b4(4:4)=rec(47:47)
          IF(b4.NE.'    ') poss(3)=.false.
          p4(1:1)=rec(118:118)
          p4(2:2)=rec(129:129)
          p4(3:3)=rec(140:140)
          p4(4:4)=rec(150:150)
          IF(p4.NE.'....') poss(3)=.false.
      END IF

* Format 4 (MPC-A)
      IF(poss(4)) THEN
          b4(1:1)=rec(8:8)
          b4(2:2)=rec(26:26)
          b4(3:3)=rec(58:58)
          b4(4:4)=rec(80:80)
          IF(b4.NE.'    ') poss(4)=.false.
          p4(1:1)=rec(30:30)
          p4(2:2)=rec(41:41)
          p4(3:3)=rec(52:52)
          p4(4:4)=rec(63:63)
          IF(p4.NE.'....') poss(4)=.false.
      END IF

c Format 5 (BAC)
      IF(poss(5)) THEN
          IF(mod(k,9).eq.2) THEN
              IF(rec(1:1).eq.'B')THEN
                  poss(5)=.true.
              ELSEIF(rec(1:1).eq.'M')THEN
                  poss(5)=.true.
              ELSE
                  poss(5)=.false.
              END IF
          ELSE
              poss(5)=.false.
          END IF
      END IF

 2    CONTINUE
 3    CONTINUE
      IF(neoh) poss(1)=.false.

* Final check
      ipf=0
      npf=0
      DO 10 i=1,nfx
      IF(poss(i)) THEN
          npf=npf+1
          ipf=i
      END IF
 10   CONTINUE

      IF(npf.EQ.1) form=fname(ipf)

      IF(form.EQ.' ') THEN
          lf=lench(filnam)
          WRITE(*,200) filnam(1:lf)
      END IF
 200  FORMAT('ERROR: format auto-detection failed for file "',A,'":'/
     +       '       please specify format explicitly')

      END
* Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 16, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D M P C A                           *
*  *                                                               *
*  *          Read orbital elements for a list of objects          *
*  *        from a file written in MPC format for asteroids        *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Input FORTRAN unit (must be already opened)
*           FILNAM    -  Input file name (for error messages)
*           OBJNAM    -  Object names
*           NOBJ      -  Number of objects
*
* OUTPUT:   DEFORB    -  Tells whether orbital elements are defined
*           DEFCN     -  Tells whether covariance/normal matrices
*                            are defined
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
*           ELEM      -  Orbital elements (ECLM J2000)
*           COVE      -  Covariance matrix of orbital elements
*           NORE      -  Normal matrix of orbital elements
*           MASS      -  Mass (solar masses)
*           H         -  H absolute magnitude (if <-100, missing)
*           G         -  G slope parameter
*           COMELE    -  Comment on orbital elements
*
* WARNING: the routine assumes that objects having DEFORB=.true.
*          have already orbital elements defined (possibly from another
*          input file) and does not overwrite them
*
      SUBROUTINE rdmpca(unit,filnam,objnam,nobj,deforb,defcn,
     +                  eltype,telem,elem,cove,nore,mass,h,g,comele)
      IMPLICIT NONE

      INTEGER unit,nobj
      DOUBLE PRECISION telem(nobj),elem(6,nobj),cove(6,6,nobj)
      DOUBLE PRECISION nore(6,6,nobj),mass(nobj),h(nobj),g(nobj)
      CHARACTER*(*) filnam,objnam(nobj),eltype(nobj),comele(nobj)
      LOGICAL deforb(nobj),defcn(nobj)

      INCLUDE 'trig.h'

      INTEGER nobjx
      PARAMETER (nobjx=10)

      INTEGER nrem,k,lf,nr,ln,lc
      DOUBLE PRECISION el1(6)
      CHARACTER*7 nmpc(nobjx),nmpc1
      CHARACTER hc*5,ep5*5,krc*10,gc*5
      LOGICAL error

      INTEGER lench
      EXTERNAL lench

      IF(nobj.GT.nobjx) STOP '**** rdmpca: nobj > nobjx ****'

* Number of remaining object (orbit not yet found)
      nrem=0
      DO 10 k=1,nobj
      IF(deforb(k)) GOTO 10
      nrem=nrem+1
      CALL mpcpds(objnam(k),nmpc(k),error)
      IF(error) THEN
          ln=lench(objnam(k))
          WRITE(*,110) objnam(k)(1:ln)
      END IF
 110  FORMAT('rdmpca: cannot understand asteroid code "',A,'"')
 10   CONTINUE
      IF(nrem.LE.0) RETURN
      lf=lench(filnam)

      nr=0
 1    CONTINUE
      READ(unit,100,END=2) nmpc1
 100  FORMAT(A7,A5,2X,A5,1X,A5,1X,F9.5,2X,F9.5,2X,F9.5,2X,F9.5,
     +       2X,F9.7,13X,F11.7)
      nr=nr+1

      DO 3 k=1,nobj
      IF(deforb(k)) GOTO 3
      IF(nmpc1.EQ.nmpc(k)) THEN
          BACKSPACE(unit)
          READ(unit,100) nmpc1,hc,gc,ep5,el1
          deforb(k)=.true.
          defcn(k)=.false.
          eltype(k)='KEP'
          CALL mpcdat(ep5,telem(k),error)
          IF(error) THEN
              WRITE(*,111) nr,filnam(1:lf)
              STOP '**** rdmpca: abnormal end ****'
          END IF
          elem(1,k)=el1(6)
          elem(2,k)=el1(5)
          elem(3,k)=el1(4)*radeg
          elem(4,k)=el1(3)*radeg
          elem(5,k)=el1(2)*radeg
          elem(6,k)=el1(1)*radeg
          mass(k)=0.d0
          IF(hc.EQ.'     ') THEN
              h(k)=-1.D9
          ELSE
              READ(hc,101) h(k)
          END IF
          IF(gc.EQ.'     ') THEN
              g(k)=-1.D9
          ELSE
              READ(gc,101) g(k)
          END IF
          WRITE(krc,107) nr
          CALL rmsp(krc,lc)
          comele(k)='read from file "'//filnam(1:lf)//
     +              '" at record '//krc(1:lc)

          nrem=nrem-1
          IF(nrem.LE.0) RETURN
      END IF
 101  FORMAT(F5.2)
 107  FORMAT(I6)
 111  FORMAT('INPUT ERROR: illegal date code at record',I6,
     +       ' of file "',A,'"')

 3    CONTINUE

      GOTO 1
 2    CONTINUE

      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 15, 1997
* Modified November 9, 1998 by Steven Chesley (chesley@dm.unipi.it)
* in order to handle three digit subscripts in IAU codes.
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         M P C P D S                           *
*  *                                                               *
*  * Computes MPC-style packed designation from official IAU code  *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    IAUCO    -  IAU code
*
* OUTPUT:   MPCCOD    -  MPC-style packed code
*           ERROR     -  Error flag (cannot understand input code)
*
      SUBROUTINE mpcpds(iauco,mpccod,error)
      IMPLICIT NONE

      CHARACTER*(*) iauco,mpccod
      LOGICAL error

      INTEGER ln,nd,i,head
      CHARACTER tsn*1

      INTEGER lench
      LOGICAL isnum,islett
      EXTERNAL lench,isnum,islett

      error=.false.
      mpccod=' '

      ln=lench(iauco)
      IF(ln.LE.0) GOTO 10

* Numbered asteroids
      IF(isnum(iauco(1:ln)) .AND. ln.LE.5) THEN
           DO 1 i=1,5-ln
           mpccod(i:i)='0'
 1         CONTINUE
           mpccod(5-ln+1:5)=iauco(1:ln)
           RETURN
      END IF

* Asteroid provisional designations (e.g., 1982QB1)
      IF(ln.LT.6 .OR. ln.GT.9) GOTO 2
      IF(.NOT.isnum(iauco(1:4))) GOTO 2
      IF(.NOT.islett(iauco(5:6))) GOTO 2
      nd=ln-6
      IF(nd.GT.0) THEN
          IF(.NOT.isnum(iauco(7:ln))) GOTO 2
      END IF

      IF(iauco(1:2).EQ.'19') THEN
          mpccod(1:1)='J'
      ELSEIF(iauco(1:2).EQ.'20') THEN
          mpccod(1:1)='K'
      ELSEIF(iauco(1:2).EQ.'18') THEN
          mpccod(1:1)='I'
      ELSE
          GOTO 2
      END IF

      mpccod(2:4)=iauco(3:5)
      IF(nd.EQ.0) THEN
          mpccod(5:6)='00'
      ELSEIF(nd.EQ.1) THEN
          mpccod(5:5)='0'
          mpccod(6:6)=iauco(7:7)
      ELSEIF(nd.EQ.2) THEN
          mpccod(5:6)=iauco(7:8)
      ELSEIF(nd.EQ.3) THEN
          read(iauco,103) head
 103      format(6x,i2)
          mpccod(5:5)=char(head+55)
          mpccod(6:6)=iauco(9:9)
      ELSE
          GOTO 2
      END IF
      mpccod(7:7)=iauco(6:6)
      RETURN

 2    CONTINUE

* Palomar-Leiden survey
      i=index(iauco(1:ln),'P-L')
      IF(i.GT.0) THEN
          IF(ln.NE.i+2) GOTO 3
          nd=i-1
          IF(nd.LT.1) GOTO 3
          mpccod(1:3)='PLS'
          DO 4 i=1,4-nd
          mpccod(3+i:3+i)='0'
 4        CONTINUE
          mpccod(8-nd:7)=iauco(1:nd)
          RETURN
      END IF
 3    CONTINUE

* Trojan surveys
      i=index(iauco(1:ln),'T-')
      IF(i.GT.0) THEN
          IF(ln.NE.i+2) GOTO 5
          tsn=iauco(i+2:i+2)
          IF(tsn.NE.'1' .AND. tsn.NE.'2' .AND. tsn.NE.'3') GOTO 5
          mpccod(1:1)='T'
          mpccod(2:2)=tsn
          mpccod(3:3)='S'
          nd=i-1
          DO 6 i=1,4-nd
          mpccod(3+i:3+i)='0'
 6        CONTINUE
          mpccod(8-nd:7)=iauco(1:nd)
          RETURN
      END IF
 5    CONTINUE

* Cannot understand input code
 10   CONTINUE
      mpccod=iauco
      error=.true.

      END
* Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it),
* Version: June 19, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          R D O E F                            *
*  *                                                               *
*  *          Read orbital elements for a list of objects          *
*  *         from a file written in internal ORBFIT format         *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  Input file name
*           OBJNAM    -  Object names (without embedded blanks)
*           NOBJ      -  Number of objects
*
* OUTPUT:   DEFORB    -  Tells whether orbital elements are defined
*           DEFCN     -  Tells whether covariance/normal matrices
*                            are defined
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
*           ELEM      -  Orbital elements (ECLM J2000)
*           COVE      -  Covariance matrix of orbital elements
*           NORE      -  Normal matrix of orbital elements
*           MASS      -  Mass (solar masses)
*           H         -  H absolute magnitude (if <-100, missing)
*           G         -  G slope parameter
*           COMELE    -  Comment on orbital elements
*
* WARNING: the routine assumes that objects having DEFORB=.true.
*          have already orbital elements defined (possibly from another
*          input file) and does not overwrite them
*
      SUBROUTINE rdoef(file,objnam,nobj,deforb,defcn,eltype,telem,
     +                 elem,cove,nore,mass,h,g,comele)
      IMPLICIT NONE

      INTEGER nobj
      DOUBLE PRECISION telem(nobj),elem(6,nobj),cove(6,6,nobj)
      DOUBLE PRECISION nore(6,6,nobj),mass(nobj),h(nobj),g(nobj)
      CHARACTER*(*) file,objnam(nobj),eltype(nobj),comele(nobj)
      LOGICAL deforb(nobj),defcn(nobj)

      INTEGER kr,k,j1,j2,lf,ln,lc,lc1,nrem
      DOUBLE PRECISION t1,gmsun,gma,gma1,enne,h1,g1,m1
      DOUBLE PRECISION rot(3,3)
      DOUBLE PRECISION elem1(6),xv(6),cove1(6,6),cove2(6,6)
      DOUBLE PRECISION de(6,6),nore1(6,6),nore2(6,6)
      CHARACTER name1*80,eltyp1*3,rsys*10,epoch*10,krc*10,nc1*80
      LOGICAL defcov,defnor,end,error

      INTEGER lench
      EXTERNAL lench

      gmsun=0.01720209895d0**2

* Number of remaining object (orbit not yet found)
      nrem=0
      DO 1 k=1,nobj
      IF(.NOT.deforb(k)) nrem=nrem+1
 1    CONTINUE
      IF(nrem.LE.0) RETURN

      CALL oporbf(file)
      lf=lench(file)

 3    CONTINUE
      CALL rdorb(name1,elem1,eltyp1,t1,cove1,defcov,nore1,defnor,
     +           h1,g1,m1,rsys,epoch,kr,end)
      IF(end) GOTO 20
* Name match is performed disregarding embedded blanks
      nc1=name1
      CALL rmsp(nc1,lc1)
      IF(lc1.LE.0) GOTO 3

      DO 4 k=1,nobj
      IF(deforb(k)) GOTO 4
      IF(nc1(1:lc1).EQ.objnam(k)) THEN
          deforb(k)=.true.
           IF(rsys.EQ.'ECLM' .AND. epoch.EQ.'J2000') THEN
              DO 14 j1=1,6
              elem(j1,k)=elem1(j1)
              DO 14 j2=1,6
              cove(j1,j2,k)=cove1(j1,j2)
              nore(j1,j2,k)=nore1(j1,j2)
 14           CONTINUE
              eltype(k)=eltyp1
          ELSE
              gma=gmsun*m1
              gma1=gma+gmsun
              IF(defcov.OR.defnor) THEN
* Transformation in cartesian coordinates
                  CALL cooder(elem1,eltyp1,gma1,xv,'CAR',enne,de)
                  IF(defcov) CALL covprs(cove1,de,6,cove2)
                  IF(defnor) THEN
                      CALL norprs(nore1,de,6,nore2,error)
                      IF(error) THEN
                          ln=lench(name1)
                          WRITE(*,120) file(1:lf),name1(1:ln)
                          defnor=.false.
                      END IF
                  END IF
* Transformation of reference system
                  CALL rotpn(rot,rsys,epoch,t1,'ECLM','J2000',0.d0)
                  CALL prodmv(elem(1,k),rot,xv(1))
                  CALL prodmv(elem(4,k),rot,xv(4))
                  DO 15 j1=1,3
                  DO 15 j2=1,3
                  de(j1,j2)=rot(j1,j2)
                  de(j1+3,j2)=0
                  de(j1,j2+3)=0
                  de(j1+3,j2+3)=rot(j1,j2)
 15               CONTINUE
                  IF(defcov) CALL covprs(cove2,de,6,cove(1,1,k))
                  IF(defnor) THEN
                      CALL norprs(nore2,de,6,nore(1,1,k),error)
                      IF(error) THEN
                          ln=lench(name1)
                          WRITE(*,120) file(1:lf),name1(1:ln)
                          defnor=.false.
                      END IF
                  END IF
              ELSE
                  CALL coocha(elem1,eltyp1,gma1,xv,'CAR',enne)
                  CALL rotpn(rot,rsys,epoch,t1,'ECLM','J2000',0.d0)
                  CALL prodmv(elem(1,k),rot,xv(1))
                  CALL prodmv(elem(4,k),rot,xv(4))
              END IF
              eltype(k)='CAR'
          END IF
          telem(k)=t1
          CALL fixcnm(defcov,defnor,defcn(k),cove(1,1,k),nore(1,1,k))
          mass(k)=m1
          h(k)=h1
          g(k)=g1
          WRITE(krc,101) kr
          CALL rmsp(krc,lc)
          comele(k)='read from file "'//file(1:lf)//
     +              '" at record '//krc(1:lc)
          nrem=nrem-1
          IF(nrem.LE.0) GOTO 20
      END IF
 101  FORMAT(I6)
 120  FORMAT(' rdoef: error in transforming normal matrix'/
     +       '        (file "',A,'", object "',A,'")')
 4    CONTINUE
      GOTO 3
 20   CONTINUE

      CALL clorbf

      END
* Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it),
*
* Version: June 19, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         F I X C N M                           *
*  *                                                               *
*  *                Fix covariance/normal matrices                 *
*  *                                                               *
*  *****************************************************************
*
* IN/OUT:   DEFCOV    -  Tells whether the covariance matrix is defined
*           DEFNOR    -  Tells whether the normal matrix is defined
*           COVE      -  Covariance matrix of orbital elements
*           NORE      -  Normal matrix of orbital elements
*
* OUTPUT:   DEFCN     -  Tells whether covariance/normal matrices
*                            are defined
*
* The purpose of this routine is to compute the covariance or
* normal matrix of orbital elements when only one of the two is
* available
*
      SUBROUTINE fixcnm(defcov,defnor,defcn,cove,nore)
      IMPLICIT NONE

      DOUBLE PRECISION cove(6,6),nore(6,6)
      LOGICAL defcov,defnor,defcn
      DOUBLE PRECISION err,roff
      INTEGER nb
      DOUBLE PRECISION tmp(6)
      INTEGER i,k,indp

      err=roff(nb)*100
      defcn=(defcov.AND.defnor)
      IF(defcn) RETURN

      IF(defcov) THEN
          DO 1 i=1,6
          DO 1 k=1,6
          nore(i,k)=cove(i,k)
 1        CONTINUE
          CALL tchol(nore,6,6,indp,err)
          defnor=(indp.EQ.0)
          IF(defnor) CALL inver(nore,tmp,6,6)
      END IF
      IF(defnor) THEN
          DO 2 i=1,6
          DO 2 k=1,6
          cove(i,k)=nore(i,k)
 2        CONTINUE
          CALL tchol(cove,6,6,indp,err)
          defcov=(indp.EQ.0)
          IF(defcov) CALL inver(cove,tmp,6,6)
      END IF

      defcn=(defcov.AND.defnor)
      IF(defcn) RETURN

      DO 3 i=1,6
      DO 3 k=1,6
      cove(i,k)=0
      nore(i,k)=0
 3    CONTINUE
      defcov=.false.
      defnor=.false.

      END
* Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: 2.1.1 April 29, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D M P C A 2                         *
*  *                                                               *
*  *          Read orbital elements of all objects                 *
*  *        as given in a file written in MPC format for asteroids *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Input FORTRAN unit (must be already opened)
*           FILNAM    -  Input file name (for error messages)
*           NOBJX      -  Number of objects (maximum)
*
* OUTPUT:   OBJNAM    -  Object names
*           DEFCN     -  Tells whether covariance/normal matrices
*                            are defined
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
*           ELEM      -  Orbital elements (ECLM J2000)
*           COVE      -  Covariance matrix of orbital elements
*           NORE      -  Normal matrix of orbital elements
*           MASS      -  Mass (solar masses)
*           H         -  H absolute magnitude (if <-100, missing)
*           G         -  G slope parameter
*           COMELE    -  Comment on orbital elements
*           NOBJ      -  Number of objects, actual number found
*
*
*
      SUBROUTINE rdmpca2(unit,filnam,objnam,nobjx,nobj,defcn,
     +                  eltype,telem,elem,cove,nore,mass,h,g,comele)
      IMPLICIT NONE

      INTEGER unit,nobj,nobjx
      DOUBLE PRECISION telem(nobjx),elem(6,nobjx),cove(6,6,nobjx)
      DOUBLE PRECISION nore(6,6,nobjx),mass(nobjx),h(nobjx),g(nobjx)
      CHARACTER*(*) filnam,objnam(nobjx),eltype(nobjx),comele(nobjx)
      LOGICAL defcn(nobjx)

      INCLUDE 'trig.h'

      INTEGER nrem,k,lf,nr,ln,lc
      DOUBLE PRECISION el1(6)
      CHARACTER*7 nmpc1
      CHARACTER hc*5,ep5*5,krc*10,gc*5
      LOGICAL error

      INTEGER lench
      EXTERNAL lench
      CALL rmsp(filnam,lf)

      nr=0
      DO 3 k=1,nobjx
         READ(unit,100,END=2) nmpc1,hc,gc,ep5,el1
 100     FORMAT(A7,A5,2X,A5,1X,A5,1X,F9.5,2X,F9.5,2X,F9.5,2X,F9.5,
     +        2X,F9.7,13X,F11.7)
         nr=nr+1
         nobj=nr
c name conversion
         CALL iaucod2(nmpc1,objnam(k),error)
         IF(error) THEN
            WRITE(*,112) nr,filnam(1:lf)
 112        FORMAT('INPUT ERROR: illegal name code at record',I6,
     +       ' of file "',A,'"')
            STOP '**** rdmpca2: name conversion error ****'
         END IF
         defcn(k)=.false.
         eltype(k)='KEP'
         CALL mpcdat(ep5,telem(k),error)
         IF(error) THEN
            WRITE(*,111) nr,filnam(1:lf)
 111        FORMAT('INPUT ERROR: illegal date code at record',I6,
     +       ' of file "',A,'"')
            STOP '**** rdmpca2: date conversion error ****'
         END IF
         elem(1,k)=el1(6)
         elem(2,k)=el1(5)
         elem(3,k)=el1(4)*radeg
         elem(4,k)=el1(3)*radeg
         elem(5,k)=el1(2)*radeg
         elem(6,k)=el1(1)*radeg
         mass(k)=0.d0
         IF(hc.EQ.'     ') THEN
            h(k)=-1.D9
         ELSE
            READ(hc,101) h(k)
 101        FORMAT(F5.2)
         END IF
         IF(gc.EQ.'     ') THEN
            g(k)=-1.D9
         ELSE
            READ(gc,101) g(k)
         END IF
         WRITE(krc,107) nr
 107     FORMAT(I6)
         CALL rmsp(krc,lc)
         comele(k)='read from file "'//filnam(1:lf)//
     +              '" at record '//krc(1:lc)
 3    ENDDO
      WRITE(*,*)' file not completely read, record ',nr
 2    CONTINUE
      RETURN
      END
* Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 21, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         O P O R B F                           *
*  *                                                               *
*  *                Open an orbital element file                   *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  File name
*
      SUBROUTINE oporbf(file)
      IMPLICIT NONE

      CHARACTER*(*) file

      INTEGER kr,lf,mjd,mjde,nn
      LOGICAL first,found,end
      CHARACTER form*20,scale*3,rec*200
      DOUBLE PRECISION sec,sece

      SAVE first

      INCLUDE 'comorb.h'

      INTEGER lench,nitchs
      EXTERNAL lench,nitchs

      DATA first/.true./

      IF(first) THEN
          orbunt=0
          orbfn=' '
          first=.false.
      END IF

      IF(orbunt.NE.0) STOP '**** oporbf: internal error (01) ****'

      orbfn=file
      CALL filopn(orbunt,orbfn,'OLD')
      CALL rdfnam(orbunt,orbfn,orbnr)
      lf=lench(orbfn)

* Format
      CALL rdfcha(orbunt,'format',.true.,form,found,kr)
      IF(form.NE.'OEF1.1') THEN
          WRITE(*,100) orbfn(1:lf)
          STOP '**** oporbf: abnormal end ****'
      END IF
 100  FORMAT('ERROR: unsupported format in file ',A)

* Record type and default orbital element type
      CALL rdfcha(orbunt,'rectype',.false.,rectyp,found,kr)
      IF(.NOT.found) THEN
 1        CONTINUE
          CALL getrsc(orbunt,rec,orbnr,end)
          IF(end) THEN
              WRITE(*,104) orbfn(1:lf)
              STOP '**** oporbf: abnormal end ****'
          END IF
          orbnr=orbnr-1
          BACKSPACE(orbunt)
          nn=nitchs(rec)
          IF(nn.EQ.1) THEN
              rectyp='ML'
          ELSEIF(nn.GE.7) THEN
              rectyp='1L'
          ELSE
              orbnr=orbnr+1
              GOTO 10
          END IF
      END IF
 104  FORMAT(' ERROR: file ',A,' is empty')
      IF(rectyp.EQ.'1L') THEN
          CALL rdfcha(orbunt,'elem',.true.,deltyp,found,kr)
      ELSEIF(rectyp.EQ.'ML') THEN
          deltyp=' '
      ELSE
          WRITE(*,101) orbfn(1:lf)
          STOP '**** oporbf: abnormal end ****'
      END IF
 101  FORMAT('ERROR: unsupported record type in file ',A)

* Default reference system
      CALL rdfref(orbunt,'refsys',.false.,dfrsty,dfrsep,found,kr)
      IF(.NOT.found) THEN
          IF(rectyp.EQ.'1L') THEN
              WRITE(*,105) orbfn(1:lf)
              STOP '**** oporbf: abnormal end ****'
          END IF
          dfrsty=' '
          dfrsep=' '
      END IF
 105  FORMAT(' ERROR: missing keyword "refsys" in file ',A)

* Default epoch for orbital elements
      CALL rdftim(orbunt,'epoch',.false.,depstr,mjd,sec,scale,deft0,kr)
      IF(deft0) THEN
          CALL cnvtim(mjd,sec,scale,mjde,sece,'TDT')
          dept0=mjde+sece/86400.d0
      END IF

      nxtend=.false.
      iicorb=36
      RETURN

 10   CONTINUE
      WRITE(*,102) orbfn(1:lf),orbnr
 102  FORMAT(' FORMAT ERROR in file ',A,' at line',I5)
      STOP '**** oporbf: abnormal end ****'

      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 7, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         C L O R B F                           *
*  *                                                               *
*  *                Close an orbital element file                  *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  File name
*
      SUBROUTINE clorbf
      IMPLICIT NONE

      INCLUDE 'comorb.h'

      IF(orbunt.LE.0) STOP '**** clorbf: internal error (01) ****'

      CALL filclo(orbunt,' ')

      orbunt=0
      orbfn=' '
      iicorb=0

      END
* Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 21, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          R D O R B                            *
*  *                                                               *
*  *     Read orbital elements from a file opened with OPORBF      *
*  *                                                               *
*  *****************************************************************
*
* The routine operates in a sequential way, returning each time
* the orbital elements of the next object contained in the file
*
* OUTPUT:   NAME      -  Name of planet/asteroid/comet
*           ELEM(6)   -  Orbital element vector
*           ELTYPE    -  Type of orbital elements (KEP/EQU/CAR)
*           T0        -  Epoch of orbital elements (MJD, TDT)
*           COVE      -  Covariance matrix of orbital elements
*           DEFCOV    -  Tells whether the covariance matrix is defined
*           NORE      -  Normal matrix of orbital elements
*           DEFNOR    -  Tells whether the normal matrix is defined
*           H         -  H absolute magnitude (if <-100, missing)
*           G         -  G slope parameter
*           MASS      -  Mass (solar masses)
*           RSYS      -  Reference system type (EQUM/EQUT/ECLM)
*           EPOCH     -  Epoch specification (J2000/OFDATE)
*           KR        -  Record number at which object is found
*           EOF       -  End-of-file flag
*
      SUBROUTINE rdorb(name,elem,eltype,t0,cove,defcov,nore,defnor,
     +                 h,g,mass,rsys,epoch,kr,eof)
      IMPLICIT NONE

      DOUBLE PRECISION elem(6),t0,h,g,cove(6,6),nore(6,6),mass
      CHARACTER*(*) name,eltype,rsys,epoch
      LOGICAL defcov,defnor,eof
      INTEGER kr

      INCLUDE 'trig.h'

      CHARACTER rec*200,rest*200,scale*3
      INTEGER lf,nit,i,k,mjd,mjde,ik
      DOUBLE PRECISION sec,sece,tmp(21),cnv(6)
      LOGICAL error,end1,noep


      INCLUDE 'comorb.h'

      INTEGER lench,nitchs
      EXTERNAL lench,nitchs

      IF(iicorb.NE.36) STOP '**** rdorb: internal error (01) ****'
      IF(orbunt.LE.0) STOP '**** rdorb: internal error (02) ****'

      rsys=dfrsty
      epoch=dfrsep
      mass=0.d0
      DO 1 i=1,6
      DO 2 k=1,6
      cove(i,k)=0
      nore(i,k)=0
 2    CONTINUE
 1    CONTINUE
      defcov=.false.
      defnor=.false.
      h=-1.d9
      g=-1.d9

      IF(rectyp.EQ.'1L') THEN
          CALL getrsc(orbunt,rec,orbnr,eof)
          IF(eof) RETURN
          CALL strcnt(rec,name,rest,error)
          IF(error) GOTO 20
          kr=orbnr
          nit=nitchs(rest)
          IF(deft0) THEN
              t0=dept0
              IF(nit.LT.8) THEN
                  READ(rest,*,ERR=20) elem
              ELSE
                  READ(rest,*,ERR=20) elem,h,g
              END IF
          ELSE
              IF(nit.LT.9) THEN
                  READ(rest,*,ERR=20) t0,elem
              ELSE
                  READ(rest,*,ERR=20) t0,elem,h,g
              END IF
          END IF
          eltype=deltyp
      ELSEIF(rectyp.EQ.'ML') THEN
          IF(nxtend) THEN
              eof=.true.
              RETURN
          END IF
          noep=.true.
* Name
          CALL getrsc(orbunt,name,orbnr,eof)
          IF(eof) RETURN
          kr=orbnr
* Orbital elements (mandatory, immediately after the name)
          CALL getrsc(orbunt,rec,orbnr,end1)
          IF(end1) GOTO 20
          IF(rec(1:4).EQ.' KEP') THEN
              READ(rec(5:),*,ERR=20) elem
              eltype='KEP'
          ELSEIF(rec(1:4).EQ.' EQU') THEN
              READ(rec(5:),*,ERR=20) elem
              eltype='EQU'
          ELSEIF(rec(1:4).EQ.' CAR') THEN
              READ(rec(5:),*,ERR=20) elem
              eltype='CAR'
          ELSE
              GOTO 20
          END IF

* Other keywords
 3        CONTINUE
          CALL getrsc(orbunt,rec,orbnr,end1)
          IF(end1) THEN
              nxtend=.true.
              GOTO 4
          END IF
          IF(rec(1:1).NE.' ') THEN
              BACKSPACE(orbunt)
              orbnr=orbnr-1
              GOTO 4
          END IF
* Epoch of elements
          IF(rec(1:4).EQ.' MJD' .OR. rec(1:4).EQ.' JD ' .OR.
     +       rec(1:4).EQ.' CAL') THEN
              CALL ch2tim(rec,mjd,sec,scale,error)
              IF(error) GOTO 20
              CALL cnvtim(mjd,sec,scale,mjde,sece,'TDT')
              t0=mjde+sece/86400.d0
              noep=.false.
          ELSEIF(rec(1:4).EQ.' MAG') THEN
              READ(rec(5:),*,ERR=20) h,g
          ELSEIF(rec(1:4).EQ.' MAS') THEN
              READ(rec(5:),*,ERR=20) mass
          ELSEIF(rec(1:4).EQ.' COV') THEN
              READ(rec(5:),*,ERR=20) (tmp(i),i=1,3)
              DO 17 k=1,6
              CALL getrsc(orbunt,rec,orbnr,end1)
              IF(end1) GOTO 20
              IF(rec(1:4).NE.' COV') GOTO 20
              READ(rec(5:),*,ERR=20) (tmp(i),i=3*k+1,3*k+3)
 17           CONTINUE
              ik=0
              DO 8 i=1,6
              DO 7 k=i,6
              ik=ik+1
              cove(i,k)=tmp(ik)
 7            CONTINUE
 8            CONTINUE
              IF(ik.NE.21) STOP '**** rdorb: internal error (03) ****'
              defcov=.true.
          ELSEIF(rec(1:4).EQ.' NOR') THEN
              READ(rec(5:),*,ERR=20) (tmp(i),i=1,3)
              DO 27 k=1,6
              CALL getrsc(orbunt,rec,orbnr,end1)
              IF(end1) GOTO 20
              IF(rec(1:4).NE.' NOR') GOTO 20
              READ(rec(5:),*,ERR=20) (tmp(i),i=3*k+1,3*k+3)
 27           CONTINUE
              ik=0
              DO 38 i=1,6
              DO 37 k=i,6
              ik=ik+1
              nore(i,k)=tmp(ik)
 37           CONTINUE
 38           CONTINUE
              IF(ik.NE.21) STOP '**** rdorb: internal error (04) ****'
              defnor=.true.
          ELSE
              GOTO 20
          END IF
          GOTO 3
 4        CONTINUE
          IF(noep) THEN
              IF(deft0) THEN
                  t0=dept0
              ELSE
                  GOTO 20
              END IF
          END IF
      ELSE
          STOP '**** rdorb: internal error (05) ****'
      END IF

* Transformation of angles in orbital elements
* and covariance matrix
      DO 9 i=1,6
      cnv(i)=1
 9    CONTINUE
      IF(eltype.EQ.'KEP') THEN
          DO 10 i=3,6
          cnv(i)=radeg
 10       CONTINUE
      ELSEIF(eltype.EQ.'EQU') THEN
          cnv(6)=radeg
      ELSEIF(eltype.EQ.'CAR') THEN
          CONTINUE
      ELSE
          STOP '**** rdorb: internal error (06) ****'
      END IF
      DO 11 i=1,6
      elem(i)=elem(i)*cnv(i)
 11   CONTINUE
      IF(defcov) THEN
          DO 13 i=1,6
          cove(i,i)=cove(i,i)*(cnv(i)**2)
          DO 12 k=i+1,6
          cove(i,k)=cove(i,k)*cnv(i)*cnv(k)
          cove(k,i)=cove(i,k)
 12       CONTINUE
 13       CONTINUE
      END IF
      IF(defnor) THEN
          DO 15 i=1,6
          nore(i,i)=nore(i,i)/(cnv(i)**2)
          DO 14 k=i+1,6
          nore(i,k)=nore(i,k)/(cnv(i)*cnv(k))
          nore(k,i)=nore(i,k)
 14       CONTINUE
 15       CONTINUE
      END IF

      eof=.false.
      RETURN

 20   CONTINUE
      lf=lench(orbfn)
      WRITE(*,200) orbfn(1:lf),orbnr
 200  FORMAT(' ERROR in file ',A,' at line',I6)
      STOP '**** rdorb: abnormal end ****'

      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 10, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         M P C D A T                           *
*  *                                                               *
*  *         Computes MJD from MPC-style packed dates              *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    PDATE     -  MPC-style packed date
*
* OUTPUT:   TJM       -  Modified Julian Date (TDT)
*           ERROR     -  Error flag (cannot understand input)
*
      SUBROUTINE mpcdat(pdate,tjm,error)
      IMPLICIT NONE

      CHARACTER*(*) pdate
      DOUBLE PRECISION tjm
      LOGICAL error

      INTEGER year,yy,month,day

      INTEGER lench
      LOGICAL isnum
      DOUBLE PRECISION tjm1
      EXTERNAL lench,isnum,tjm1

      error=.true.
      tjm=0.d0
      IF(lench(pdate).NE.5) RETURN

* Year
      IF(pdate(1:1).EQ.'I') THEN
          year=1800
      ELSEIF(pdate(1:1).EQ.'J') THEN
          year=1900
      ELSEIF(pdate(1:1).EQ.'K') THEN
          year=2000
      ELSE
          RETURN
      END IF
      READ(pdate(2:3),100,ERR=10) yy
 100  FORMAT(I2)
      year=year+yy

* Month
      IF(isnum(pdate(4:4))) THEN
          READ(pdate(4:4),101,ERR=10) month
      ELSEIF(pdate(4:4).EQ.'A') THEN
          month=10
      ELSEIF(pdate(4:4).EQ.'B') THEN
          month=11
      ELSEIF(pdate(4:4).EQ.'C') THEN
          month=12
      ELSE
          RETURN
      END IF
 101  FORMAT(I1)

* Day
      IF(isnum(pdate(5:5))) THEN
          READ(pdate(5:5),101,ERR=10) day
      ELSE
          day=ichar(pdate(5:5))-55
          IF(day.LT.10 .OR. day.GT.31) GOTO 10
      END IF
      tjm=tjm1(day,month,year,0.d0)
      error=.false.

 10   CONTINUE

      END



