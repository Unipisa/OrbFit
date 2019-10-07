c ===========LIBRARY math_lib=====================
c note: many routines to be removed in fortran90
c MATHEMATICAL ROUTINES:
c   ANGLES
c princ		principal value of an angle
c primea	mean of two angle variables 
c pridif        difference of two angles
c intlo		lower integer
c inthi         higher integer
c polar         transformation from cartesian to polar coordinates
c   3D LINEAR ALGEBRA
c prscal	scalar product of two 3D vectors
c prvec		vector product of two 3D vectors
c vsize		size of a 3D vector
c prodmv	product of a 3x3 matrix by a 3D vector
c prodmm	product of two 3x3 matrices
c lincom	linear combination of two 3D vectors
c vdiff          difference of 2-vectors
c   GENERAL LINEAR ALGEBRA
c prscag 	scalar product of two general vectors
c mulmat	product between generic matrices
c mulmat_cut    product between generic matrices, when one is not full
c lincog        linear combination of two general vectors
c matvet        generic matrix to vector copy
c vetmat	generic vector to matrix copy
c mulmav        general matrix by vector product
c transp  	transposition of a general matrix
c mcopy		matrix copy
c vcopy		vector copy
c vsumg           vector sum
c eye		identity matrix
c bilin           bilinear form of two general vectors
c   LINEAR ALGEBRA, but used only in iers_ser.f
c assmat
c pdmat
c pd1mat
c pd2mat
c prodvs
c rotmt1
c rotmt2
c summat
c sumv
c trsp3
c   LEAST SQUARES
c snorm		weighed RMS norm sqrt(xTWx/n)
c snormd        RMS norm sqrt(xTx/n)
c tchol		Tcholesky factorization of a positive-defined matrix
c inver         Tcholesky inversion
c matin		inversion of a matrix (Gauss method)
c tchinv        standard interface for tchol, inver
c inv22         inversion of a 2x2 matrix
c covprs	covariance propagation for a NxN (square) matrix
c norprs        normal matrix propagation for a NxN (square) matrix
c   SORTING
c heapsort      sort by real, leaving in place and storing indexes
c heapsorti     sort by integer, leaving in place and storing indexes
c   ROUNDING OFF
c roff		computation of rounding-off
c truncat       truncation accounting for rounding off
c
c
c HEADERS math_lib.o: trig.h (public)
c
c==============================================
*
* Computes the principal value of an angle
*
      double precision function princ (a)
      double precision a
*
      include 'trig.h'
*
      princ=dmod(a,dpig)
      if(princ.lt.0.d0)princ=princ+dpig
      return
      end
c
      double precision function prscal(a,b)
c
c    Computes the scalar product of two 3d vectors
c
      double precision a(3),b(3)
      prscal=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end
*
*  ***************************************************************
*  *                                                             *
*  *                          P R V E C                          *
*  *                                                             *
*  *                 Vector product:   C = A x B                 *
*  *                                                             *
*  ***************************************************************
*
      subroutine prvec(a,b,c)
      implicit double precision (a-h,o-z)
      dimension a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end
* Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: August 16, 1996
* ---------------------------------------------------------------------
*
*  ***************************************************************
*  *                                                             *
*  *                          V S I Z E                          *
*  *                                                             *
*  *                     Size of a 3-D vector                    *
*  *                                                             *
*  ***************************************************************
*
      DOUBLE PRECISION FUNCTION vsize(x)
      IMPLICIT NONE

      DOUBLE PRECISION x(3),s

      s=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
      vsize=SQRT(s)

      END
* Author: Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 5, 1996
*
*  ***************************************************************
*  *                                                             *
*  *                         P R O D M V                         *
*  *                                                             *
*  *            Product of a 3x3 matrix by a 3-D vector          *
*  *                           y = Ax                            *
*  *                                                             *
*  ***************************************************************
*
      SUBROUTINE prodmv(y,a,x)
      IMPLICIT NONE
      DOUBLE PRECISION x(3),y(3),a(3,3),s
      INTEGER j,l

      DO 2 j=1,3
      s=0.d0
      DO 1 l=1,3
 1    s=s+a(j,l)*x(l)
 2    y(j)=s

      END
      subroutine mulmat(a,na,ma,b,nb,mb,c)
      implicit double precision (a-h,o-z)
      dimension a(na,ma),b(nb,mb),c(na,mb)
      if (ma.ne.nb) then
         write (*,*) 'dimensioning error in mulmat'
         stop
      endif
      do 3 i=1,na
        do 2 j=1,mb
          c(i,j)=0.d0
          do 1 k=1,ma
 1           c(i,j)=c(i,j)+a(i,k)*b(k,j)
 2      continue
 3    continue
      return
      end
      subroutine mulmat_cut(a,nax,na,ma,b,nb,mb,c)
      implicit double precision (a-h,o-z)
      dimension a(nax,ma),b(nb,mb),c(nax,mb)
      if (ma.ne.nb) then
         write (*,*) 'dimensioning error in mulmat'
         stop
      endif
      do 3 i=1,na
        do 2 j=1,mb
          c(i,j)=0.d0
          do 1 k=1,ma
 1           c(i,j)=c(i,j)+a(i,k)*b(k,j)
 2      continue
 3    continue
      return
      end
c
      subroutine lincom(v1,a,v2,b,v3)
c
c    Computes the linear combination vector v3 of the 
c    three-elements vectors v1,v2 with coefficients a,b
c
      implicit double precision (a-h,o-z)
      dimension v1(3),v2(3),v3(3)
      do 40 i=1,3
        v3(i)=a*v1(i)+b*v2(i)
40    continue
      return
      end

c ===========================================================
c                                                               
c                            t c h o l                          
c                                                               
c      Tcholesky factorization of a positive-defined matrix     
c                                                               
c             originally written by prof. Luigi Mussio                
c         Istituto di Topografia del Politecnico di Milano     
c                                                               
c ===========================================================
c
c    warning: only the upper triangle of the matrix is processed
c
c
c input:    a(nmax,nmax)   - matrix to be processed
c           nmax           - first dimension of a as declared in the dim
c                            statement of the calling program
c           n              - actual dimension of a
c           err            - control on pivot (if <0, control is automatic)
c
c
c output:   a              - contains the tcholesky factorization of the
c                            matrix (to be supplied to subroutine inver)
c           indp           - if it is different from zero, the input mat
c                            was not positive defined and the algorithm
c                            does not work
c
      subroutine tchol(a,nmax,n,indp,err)
      implicit none
      integer n,nmax,indp
      double precision a(nmax,n),err
c ========END INTERFACE=========================
      double precision as,errl,s
      integer i,k,ii,l,j
      errl=err
      as=0.d0
      do 1 i=1,n
      do 1 k=1,i
        as=max(as,abs(a(k,i)))
 1    continue
c in this version we are not checking the roundoff value
c of the machine, beacuse the routine roff is too easily machine dependent
c thus err must be positive 
      errl=errl*as 
      do 40 i=1,n
        l=i-1
        s=0.d0
        if(l.eq.0) goto 15
        do  k=1,l
          s=s+a(k,i)*a(k,i)
        enddo
   15   a(i,i)=a(i,i)-s
        if(a(i,i).le.errl)then
           indp=i
           return
        endif
        a(i,i)=sqrt(a(i,i))
        if(i.eq.n) goto 40
        ii=i+1
        do 30 j=ii,n
          s=0.d0
          if(l.eq.0)goto 25
          do  k=1,l
            s=s+a(k,i)*a(k,j)
          enddo
   25     a(i,j)=(a(i,j)-s)/a(i,i)
   30   continue
   40 continue
      indp=0
      return
      end
c ===========================================================
c                                                               
c                            i n v e r                          
c                                                               
c              inversion of a positive-defined matrix           
c                                                               
c                  written by prof. luigi mussio                
c         istituto di topografia del politecnico di milano      
c                                                               
c
c input:    a(nmax,nmax)   - input matrix, previously factorized by
c                            subroutine tchol
c           v(nmax)        - work area
c           nmax           - first dimension of a as declared in the dim
c                            statement of the calling program
c           n              - actual dimension of a
c
c
c output:   a              - inverse of input matrix
c 
c
      subroutine inver(a,v,nmax,n)
      implicit none
      integer nmax,n
      double precision a(nmax,n),v(n)
c ==========END INTERFACE=========
      integer i,l,k,j
      double precision s
c ================================
      do 8 i=n,1,-1
        l=i+1
        if(i.eq.n)goto 4
        do 3 j=n,l,-1
          s=0.d0
          do 2 k=n,l,-1
            if(j.lt.k)goto 1
            s=s+a(i,k)*a(k,j)
            goto 2
    1       s=s+a(i,k)*a(j,k)
 2        continue
          v(j)=-s/a(i,i)
    3   continue
    4   s=0.d0
        if(i.eq.n)goto 7
        do  k=n,l,-1
          s=s+a(i,k)*v(k)
        enddo
        do  j=n,l,-1
          a(i,j)=v(j)
        enddo
    7   a(i,i)=(1.d0/a(i,i)-s)/a(i,i)
    8 continue
      do 9 i=2,n
        do 9 j=1,i-1
    9     a(i,j)=a(j,i)
      return
      end
* Author: Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 5, 1996
*
*  ***************************************************************
*  *                                                             *
*  *                          M U L T 3                          *
*  *                                                             *
*  *         Multiplication (to the left) of a 3x3 matrix        *
*  *                                                             *
*  ***************************************************************
*
* INPUT:    R,ROT     -  3x3 matrices
*
* OUTPUT:   ROT = R * ROT
*
      subroutine mult3(r,rot)
      implicit none

      double precision r(3,3),rot(3,3),rot1(3,3)
      integer i,j

      do 1 i=1,3
      do 1 j=1,3
 1    rot1(i,j)=r(i,1)*rot(1,j)+r(i,2)*rot(2,j)+r(i,3)*rot(3,j)
      do 2 i=1,3
      do 2 j=1,3
 2    rot(i,j)=rot1(i,j)
      end
* Author: Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 5, 1996
*
*  ***************************************************************
*  *                                                             *
*  *                          R O T M T                          *
*  *                                                             *
*  *              Rotation matrix around k axis	                 *
*  *                                                             *
*  ***************************************************************
*
*  If X are "old" coordinates and X' are "new" coordinates (referred
*  to a frame which is rotated by an angle alpha around k axis in
*  direct sense), then X' = R X
*
      subroutine rotmt(alpha,r,k)
      implicit none

      double precision alpha,r(3,3)
      double precision cosa,sina
      integer k,i1,i2,i3

      if(k.lt.1.or.k.gt.3)stop' **** ROTMT: k = ??? ****'
      cosa=cos(alpha)
      sina=sin(alpha)
      i1=k
      if(i1.gt.3)i1=i1-3
      i2=i1+1
      if(i2.gt.3)i2=i2-3
      i3=i2+1
      if(i3.gt.3)i3=i3-3

      r(i1,i1)=1.d0
      r(i1,i2)=0.d0
      r(i1,i3)=0.d0
      r(i2,i1)=0.d0
      r(i2,i2)=cosa
      r(i2,i3)=sina
      r(i3,i1)=0.d0
      r(i3,i2)=-sina
      r(i3,i3)=cosa

      end
* Author: Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 5, 1996
*
*  ***************************************************************
*  *                                                             *
*  *                           T R S P                           *
*  *                                                             *
*  *                  Transpose of a 3x3 matrix                  *
*  *                                                             *
*  ***************************************************************
*
*
* INPUT:    R(3,3)   -  3x3 matrix
*
* OUTPUT:   R(3,3)   -  Transpose of input matrix
*
      subroutine trsp(r)
      implicit none

      double precision r(3,3),rt
      integer i,j

      do 1 i=1,2
      do 1 j=i+1,3
      rt=r(i,j)
      r(i,j)=r(j,i)
 1    r(j,i)=rt
      end
      subroutine matvet(b,nr,nc,vec)
      implicit double precision (a-h,o-z)
      dimension b(nr,nc),vec(nr*nc)
      do 21 j=1,nc
        do 22 i=1,nr
           vec(nr*(j-1)+i)=b(i,j)
 22    continue
 21   continue
      return
      end
* Author: Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 5, 1996
*
*  ***************************************************************
*  *                                                             *
*  *                         P R O D M M                         *
*  *                                                             *
*  *                 Product of two 3x3 matrices                 *
*  *                           A = BC                            *
*  *                                                             *
*  ***************************************************************
*
      subroutine prodmm(a,b,c)
      implicit none

      double precision a(3,3),b(3,3),c(3,3),s
      integer j,k,l

      do 2 j=1,3
      do 2 k=1,3
      s=0.d0
      do 1 l=1,3
 1    s=s+b(j,l)*c(l,k)
 2    a(j,k)=s

      end
      subroutine vetmat(vec,nv,b,nr,nc)
      implicit double precision (a-h,o-z)      
      dimension vec(nv),b(nr,nc)
      if ((nr*nc).ne.nv) then
        write(*,*)'error: nr*nc must be = nv'
        stop
      endif
      do 12 j=1,nc
        do 13 i=1,nr
          b(i,j)=vec(nr*(j-1)+i)
 13     continue
 12   continue
      return
      end
c ==============================================
c SNORM norm of vector v according to metric a
      DOUBLE PRECISION FUNCTION snorm(v,a,n,nnx)
      IMPLICIT NONE
      INTEGER nnx,n,i,k
      DOUBLE PRECISION  v(nnx),a(nnx,nnx)
      snorm=0.d0
      DO  i=1,n
        DO  k=1,n
          snorm=snorm+v(i)*v(k)*a(i,k)
        ENDDO
      ENDDO
      snorm=sqrt(snorm/n)
      RETURN
      END
c ===========================================
c SNORMD RMS norm of vector v with diagonal metric a
c warning: if a weight is zero, the component is not
c counted in the vector length
      DOUBLE PRECISION FUNCTION snormd(v,a,n,nused)
      IMPLICIT NONE
      INTEGER n,i,nused
      DOUBLE PRECISION v(n),a(n)
      snormd=0.d0
      nused=0
      DO  i=1,n
        snormd=snormd+(v(i)**2)*a(i)
        IF(a(i).gt.0.d0)nused=nused+1
      ENDDO
      snormd=sqrt(snormd/nused)
      RETURN
      END
c
      subroutine vdiff(v,w,d)
c
c  Computes the difference vector d between the vectors v,w 
c
      implicit double precision (a-h,o-z)
      dimension v(3),w(3),d(3)
      do 3 i=1,3
 3       d(i)=v(i)-w(i)
      return
      end
      subroutine transp(a,n,m,b)
      double precision a(n,m),b(m,n)
      do 50 i=1,n
        do 60 j=1,m
          b(j,i)=a(i,j)
 60     continue  
 50   continue  
      return
      end 
c ================================================
*
*  ***************************************************************
*  *                                                             *
*  *                          M A T I N                          *
*  *                                                             *
*  *              Gauss' method for matrix inversion             *
*  *         and solution of systems of linear equations         *
*  *                                                             *
*  ***************************************************************
*
* INPUT:    A(i,j)    -   Matrix of coefficients of the system
*           N         -   Order of A(i,j)
*           L         -   Number of linear systems to be solved (the
*                         right hand  sides are stored in A(i,N+j),
*                         j=1,L)
*           NL        -   First dimension of A(i,j)
*           INVOP     -   If =1 the inverse of A(i,j) is computed
*                         explicitly,
*                         if =0 only the linear systems are solved
*
* OUTPUT:   A(i,j)    -   If INVOP=1, contains the inverse of the input matrix;
*                         if INVOP=0, contains the triangular
*                         factorization of the input matrix.
*                         In both cases A(i,N+j), j=1,L contain the
*                         solutions of the input systems.
*           DET       -   Determinant of the input matrix A(i,j)
*           ISING     -   If =0 the algorithm was completed
*                         successfully,
*                         if =1 the input matrix was singular
*
      subroutine matin(a,det,n,l,nl,ising,invop)
      implicit double precision (a-h,o-z)
      parameter (nmax=1000)
      dimension a(nl,n+l),vet(nmax)
      integer p(nmax)
      ising=1
      det=0.d0
      if(n.gt.nmax)stop ' **** matin: n > nmax ****'
      do 5 j=1,n
 5    p(j)=j
      do 4 j=1,n-1
      amax=0.d0
      do 1 i=j,n
      if(amax.ge.dabs(a(i,j)))goto 1
      npiv=i
      amax=dabs(a(i,j))
 1    continue
      if(amax.eq.0.d0)return
      if(npiv.eq.j)goto 6
      nsc=p(npiv)
      p(npiv)=p(j)
      p(j)=nsc
      do 2 i=1,n+l
      sc=a(npiv,i)
      a(npiv,i)=a(j,i)
 2    a(j,i)=sc
 6    do 3 i=j+1,n
 3    a(i,j)=a(i,j)/a(j,j)
      do 4 i=j+1,n
      do 4 k=j+1,n+l
 4    a(i,k)=a(i,k)-a(i,j)*a(j,k)
      if(l.eq.0)goto 7
      do 9 i=1,l
      do 9 j=n,1,-1
      if(j.eq.n)goto 9
      do 8 k=j+1,n
 8    a(j,n+i)=a(j,n+i)-a(j,k)*a(k,n+i)
 9    a(j,n+i)=a(j,n+i)/a(j,j)
 7    det=1.d0
      do 10 j=1,n
 10   det=det*a(j,j)
      ising=0
      if(invop.eq.0)return
      if(n.ne.1)goto 20
      a(1,1)=1.d0/a(1,1)
      return
 20   do 12 i=2,n
      do 12 j=1,i-1
      a(i,j)=-a(i,j)
      if(j+2.gt.i)goto 12
      do 11 k=j+1,i-1
 11   a(i,j)=a(i,j)-a(i,k)*a(k,j)
 12   continue
      do 15 k=n,1,-1
      do 14 i=1,n
      vet(i)=0.d0
      if(i.eq.k)vet(i)=1.d0
      if(k.gt.i)vet(i)=a(k,i)
      if(k.eq.n)goto 14
      do 13 j=k+1,n
 13   vet(i)=vet(i)-a(k,j)*a(j,i)
 14   continue
      sc=a(k,k)
      do 15 i=1,n
 15   a(k,i)=vet(i)/sc
 18   nc=0
      do 16 j=1,n
      if(j.eq.p(j))goto 16
      nsc=p(j)
      p(j)=p(nsc)
      p(nsc)=nsc
      do 17 i=1,n
      sc=a(i,nsc)
      a(i,nsc)=a(i,j)
 17   a(i,j)=sc
      nc=1
 16   continue
      if(nc.eq.1)goto 18
      return
      end
c==============================================
c identitiy matrix ndimxndim in output a
         subroutine eye(ndim,a)
         implicit none
         integer ndim,i,j
         double precision a(ndim,ndim)
         do  i=1,ndim
           do  j=1,ndim
             if(i.eq.j)then
                a(i,j)=1.d0
             else
                a(i,j)=0.d0
             endif
           enddo
         enddo
         return
         end
c ======================================
c VCOPY - vector copy for n-vector
c
      subroutine vcopy(n,a,b)
      implicit none
      integer n,i
      double precision a(n),b(n)
      do  i=1,n
        b(i)=a(i)
      enddo
      return
      end
c ======================================
c MCOPY - matrix copy for n x m matrix
c
      subroutine mcopy(n,m,a,b)
      implicit none
      integer n,m,i,j
      double precision a(n,m),b(n,m)
      do   i=1,n
        do   j=1,m
          b(i,j)=a(i,j)
        enddo
      enddo
      return
      end
c========================
c  primea
c  arithmetic mean of two angles, 
c  expressed as principal values
       double precision function primea(a,b)
       implicit none
       include 'trig.h'
       double precision a,b,princ
       a=princ(a)
       b=princ(b)
       if(abs(b-a).gt.pig)then
           primea=(a+b)/2.d0+pig
           primea=princ(primea)
       else
           primea=(a+b)/2.d0
       endif
       return
       end

c========================
c  pridif
c  difference of two angles, 
c  expressed as principal values
       double precision function pridif(a,b)
       implicit none
       include 'trig.h'
       double precision a,b,princ

       a=princ(a)
       b=princ(b)
       pridif=a-b
       if(pridif.gt.pig)then
           pridif=pridif-dpig
       elseif(pridif.lt.-pig)then
           pridif=pridif+dpig
       endif
       return
       end

* Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: April 2, 1996
* ---------------------------------------------------------------------
*
*  ***************************************************************
*  *                                                             *
*  *                          I N T L O                          *
*  *                                                             *
*  *                        Lower integer                        *
*  *                                                             *
*  ***************************************************************
*
* INPUT:    A         -  Value
*
* OUTPUT:   INTLO     -  Lower integer INTLO <= A < INTLO+1
*
      INTEGER FUNCTION intlo(a)
      IMPLICIT NONE

      DOUBLE PRECISION a

      intlo=a
      IF(intlo.GT.a) intlo=intlo-1

      END
* Copyright (C) 1999 by A. Milani
* Version: February 4, 1999
* ---------------------------------------------------------------------
*
*  ***************************************************************
*  *                                                             *
*  *                          I N T H I                          *
*  *                                                             *
*  *                        Higher integer                       *
*  *                                                             *
*  ***************************************************************
*
* INPUT:    A         -  Value
*
* OUTPUT:   INTLO     -  Higher integer INTLO-1 <= A < INTLO
*
      INTEGER FUNCTION inthi(a)
      IMPLICIT NONE

      DOUBLE PRECISION a

      inthi=a
      IF(inthi.LT.a) inthi=inthi+1

      END
c =======================================================
c computation of a bilinear form
      double precision function bilin(v1,v2,nv,w,nx,n)
      implicit none
      integer nv,nx,n
      double precision v1(nv),v2(nv),w(nx,n)
      integer i,j
c  control on dimensions
      if(nv.ne.n)write(*,*)' dimension quadr, nv=',nv,' ne n=',n
      if(n.gt.nx)write(*,*)' dimension quadr, n=',n,' gt nx=',nx
      bilin=0.d0
      do 10 i=1,n
      do 10 j=1,n
 10   bilin=bilin+v1(i)*w(i,j)*v2(j)
      return
      end
c =======================================================
c linear combination of general n vectors
      subroutine lincog(n,v1,a,v2,b,v3)
c
c    Computes the linear combination vector v3 of the 
c    n-elements vectors v1,v2 with coefficients a,b
c
      implicit none
      integer n,i
      double precision v1(n),v2(n),v3(n),a,b
      do 40 i=1,n
        v3(i)=a*v1(i)+b*v2(i)
40    continue
      return
      end
c =======================================================
c  Computes the sum vector d between the n-vectors v,w 
c
      subroutine vsumg(n,v,w,d)
      implicit none
      integer n,i
      double precision v(n),w(n),d(n)
      do 3 i=1,n
 3       d(i)=v(i)+w(i)
      return
      end
c =======================================================
c  genneral scalar product of n-vectors
      double precision function prscag(n,a,b)
      implicit none
      integer n,i
      double precision a(n),b(n)
      prscag=0.d0
      do  i=1,n
        prscag=prscag+a(i)*b(i)
      enddo
      return
      end
c =====================================================================
c inverse of a 2x2 matrix
      subroutine inv22(a,b,deta)
      implicit none
      double precision a(2,2),b(2,2),deta
      deta=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      if(deta.eq.0.d0)then
         write(*,*)' deta is zero '
      endif
      b(1,1)=a(2,2)/deta
      b(2,2)=a(1,1)/deta
      b(1,2)=-a(1,2)/deta
      b(2,1)=-a(2,1)/deta
      return
      end
c =====================================================================
c matrix multiplied by column vector
      subroutine mulmav(a,na,ma,b,nb,c)
      implicit none
      integer na,ma,nb,i,k
      double precision a(na,ma),b(nb),c(na)
      if (ma.ne.nb) then
         write (*,*) 'dimensioning error in mulmav'
         stop
      endif
      do 3 i=1,na
        c(i)=0.d0
        do 1 k=1,ma
          c(i)=c(i)+a(i,k)*b(k)
 1      continue
 3    continue
      return
      end      
c ===========================================================
c Cholewski method for inversion
c requires a workspace vector, but does not have dimension limits
c ======================================================
      subroutine tchinv(c,n,cinv,ws,indp)
      implicit none
c input: dimension of input matrix c
      integer n
      double precision c(n,n)
c output: possible degeneracy (0=no, else row where pivot is too small)
      integer indp
      double precision cinv(n,n)
c workspace
      double precision ws(n)
c end interface
      integer i,nb
      double precision err,roff,omax,omin,cond
c =========================================
      call mcopy(n,n,c,cinv)
      err=roff(nb)*100
c matrix factorized by Tcholesky method: 
      call tchol(cinv,n,n,indp,err)
      if(indp.ne.0)then
         write(*,*)' indp=',indp,' in tchol'
      endif
c ===========================================================
c Control of conditioning number and inversion of matrix
      omax=cinv(1,1)
      omin=cinv(1,1)
      do 19 i=2,n
        if (cinv(i,i).gt.omax) then
           omax=cinv(i,i)
        endif
        if (cinv(i,i).lt.omin) then
           omin=cinv(i,i)
        endif
 19   continue
      cond=(omax/omin)**2
      write(*,111)n,n,cond
 111  format(' Conditioning of ',i2,'x',i2,'matrix=',d12.4)
c ===========================================================
c inversion
      call inver(cinv,ws,n,n)
      return
      end

* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 26, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         C O V P R S                           *
*  *                                                               *
*  *      Covariance propagation for a NxN (square) matrix         *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    C1        -  Covariance matrix E<dX1 dX1'>
*           JAC       -  Jacobian matrix dX2/dX1
*           N         -  Dimension (logical AND physical) of matrices
*
* OUTPUT:   C2        -  Covariance matrix E<dX2 dX2'> = JAC C1 JAC'
*
      SUBROUTINE covprs(c1,jac,n,c2)
      IMPLICIT NONE

      INTEGER n
      DOUBLE PRECISION c1(n,n),c2(n,n),jac(n,n)

      INTEGER nx
      PARAMETER (nx=10)

      INTEGER i,j,k
      DOUBLE PRECISION s,tmp(nx,nx)

      IF(n.GT.nx) STOP '**** covprs: n > nx ****'

      DO 2 i=1,n
      DO 2 k=1,n
      s=0
      DO 1 j=1,n
      s=s+jac(i,j)*c1(j,k)
 1    CONTINUE
      tmp(i,k)=s
 2    CONTINUE

      DO 4 i=1,n
      DO 4 k=1,n
      s=0
      DO 3 j=1,n
      s=s+tmp(i,j)*jac(k,j)
 3    CONTINUE
      c2(i,k)=s
 4    CONTINUE

      END
* Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 19, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         N O R P R S                           *
*  *                                                               *
*  *              Propagation of a normal matrix (NxN)             *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    NOR1      -  Normal matrix (X1 coord.)
*           JAC       -  Jacobian matrix dX2/dX1
*           N         -  Dimension (logical AND physical) of matrices
*
* OUTPUT:   NOR2      -  Normal matrix (X2 coord.) =
*                           = JAC'^(-1) NOR1 JAC^(-1)
*           ERROR     -  Error flag (singular Jacobian matrix)
*
      SUBROUTINE norprs(nor1,jac,n,nor2,error)
      IMPLICIT NONE

      INTEGER n
      DOUBLE PRECISION nor1(n,n),nor2(n,n),jac(n,n)
      LOGICAL error

      INTEGER nx
      PARAMETER (nx=10)

      INTEGER i,j,k,ising
      DOUBLE PRECISION s,tmp(nx,nx),jacinv(nx,nx),det

      IF(n.GT.nx) STOP '**** norprs: n > nx ****'

* Inversion of Jacobian matrix
      DO 5 i=1,n
      DO 5 k=1,n
      jacinv(i,k)=jac(i,k)
 5    CONTINUE
      CALL matin(jacinv,det,n,0,nx,ising,1)
      error=(ising.NE.0)
      IF(error) THEN
          DO 6 i=1,n
          DO 6 k=1,n
          nor2(i,k)=0
 6        CONTINUE
          RETURN
      END IF

      DO 2 i=1,n
      DO 2 k=1,n
      s=0
      DO 1 j=1,n
      s=s+jacinv(j,i)*nor1(j,k)
 1    CONTINUE
      tmp(i,k)=s
 2    CONTINUE

      DO 4 i=1,n
      DO 4 k=1,n
      s=0
      DO 3 j=1,n
      s=s+tmp(i,j)*jacinv(j,k)
 3    CONTINUE
      nor2(i,k)=s
 4    CONTINUE

      END
* Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 11, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          P O L A R                            *
*  *                                                               *
*  *    Transformation from cartesian to polar coordinates         *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    X         -  Cartesian coordinates
*
* OUTPUT:   ALPHA     -  Right ascension (rad)
*           DELTA     -  Declination (rad)
*           R         -  Distance
*
      SUBROUTINE polar(x,alpha,delta,r)
      IMPLICIT NONE

      DOUBLE PRECISION x(3),alpha,delta,r

      INCLUDE 'trig.h'

      DOUBLE PRECISION send,cosd,sina,cosa

      r=SQRT(x(1)**2+x(2)**2+x(3)**2)
      IF(r.EQ.0.d0) THEN
          alpha=0.d0
          delta=0.d0
          RETURN
      END IF

      send=x(3)/r
      delta=ASIN(send)
      cosd=COS(delta)
      IF(cosd.EQ.0.d0) THEN
          alpha=0.d0
          RETURN
      END IF

      cosa=x(1)/(r*cosd)
      sina=x(2)/(r*cosd)
      alpha=ATAN2(sina,cosa)
      IF(alpha.LT.0.d0) alpha=alpha+dpig

      END
* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 10, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         A S S M A T                           *
*  *                                                               *
*  *                   Matrix assignment (copy)                    *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    R2        -  Input matrix
*
* OUTPUT:   R1        -  R1 + R2
*
      SUBROUTINE assmat(r1,r2)
      IMPLICIT NONE

      DOUBLE PRECISION r1(3,3),r2(3,3)

      INTEGER i,j

      DO 1 i=1,3
      DO 1 j=1,3
      r1(i,j)=r2(i,j)
 1    CONTINUE

      END
* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 9, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          P D M A T                            *
*  *                                                               *
*  *                  Product of two matrices                      *
*  *  (on output, the second input matrix contains the product)    *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    R1        -  R1 matrix
*           R2        -  R2 matrix
*
* OUTPUT:   R2        -  R1 R2
*
      SUBROUTINE pdmat(r1,r2)
      IMPLICIT NONE

      DOUBLE PRECISION r1(3,3),r2(3,3)

      DOUBLE PRECISION r(3,3)

      CALL prodmm(r,r1,r2)
      CALL assmat(r2,r)

      END
* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 10, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         P D 1 M A T                           *
*  *                                                               *
*  *      First time derivative of the product of two matrices     *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    R1        -  R1 matrix
*           R1D       -  dR1/dt
*           R2        -  R2 matrix
*           R2D       -  dR2/dt
*
* OUTPUT:   R2D       -  d(R1 R2)/dt
*
* The time derivative of the product R1 R2 is computed according
* to the formula
* 
* d(R1 R2)/dt = dR1/dt R2 + R1 dR2/dt
*
      SUBROUTINE pd1mat(r1,r1d,r2,r2d)
      IMPLICIT NONE

      DOUBLE PRECISION r1(3,3),r1d(3,3),r2(3,3),r2d(3,3)

      DOUBLE PRECISION p1(3,3),p2(3,3)

      CALL prodmm(p1,r1d,r2)
      CALL prodmm(p2,r1,r2d)
      CALL summat(r2d,p1,p2)

      END
* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 10, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         P D 2 M A T                           *
*  *                                                               *
*  *     Second time derivative of the product of two matrices     *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    R1        -  R1 matrix
*           R1D       -  dR1/dt
*           R1DD      -  d^2 R1/dt^2
*           R2        -  R2 matrix
*           R2D       -  dR2/dt
*           R2DD      -  d^2 R2/dt^2
*
* OUTPUT:   R2DD      -  d^2 (R1 R2)/dt^2
*
* The time derivative of the product R1 R2 is computed according
* to the formula
* 
* d^2 (R1 R2)/dt^2 =
*       = d^2 R1/dt^2 R2 + 2 (dR1/dt)*(dR2/dt) + R1 d^2 R2/dt^2
*
      SUBROUTINE pd2mat(r1,r1d,r1dd,r2,r2d,r2dd)
      IMPLICIT NONE

      DOUBLE PRECISION r1(3,3),r1d(3,3),r1dd(3,3)
      DOUBLE PRECISION r2(3,3),r2d(3,3),r2dd(3,3)

      INTEGER i,j
      DOUBLE PRECISION p1(3,3),p2(3,3),p3(3,3)

      CALL prodmm(p1,r1dd,r2)
      CALL prodmm(p2,r1d,r2d)
      CALL prodmm(p3,r1,r2dd)

      DO 1 i=1,3
      DO 1 j=1,3
      r2dd(i,j)=p1(i,j)+2.d0*p2(i,j)+p3(i,j)
 1    CONTINUE

      END
* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 10, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         P R O D V S                           *
*  *                                                               *
*  *            Scalar by vector product (Y = A X)                 *
*  *                                                               *
*  *****************************************************************
*
      SUBROUTINE prodvs(y,a,x)
      IMPLICIT NONE

      DOUBLE PRECISION y(3),x(3),a

      INTEGER i

      DO 1 i=1,3
      y(i)=a*x(i)
 1    CONTINUE

      END
* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 10, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R O T M T 1                           *
*  *                                                               *
*  *           First time derivative of a rotation matrix          *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    A         -  Rotation angle alpha (rad)
*           K         -  (=1,2,3) Rotation axis (1=x, 2=y, 3=z)
*           ADOT      -  First time derivative of the rotation angle
*                        dalpha/dt (rad/arbitrary time unit)
*
* OUTPUT:   R         -  First time derivative of the rotation matrix
*                        dR_k/dt (rad/arbitrary time unit)
*
* The time derivative of the rotation matrix R_k(alpha) (for its
* definition, see SUBROUTINE ROTMT) is computed according to the formula
* 
* dR_k/dt = dR_k/dalpha * dalpha/dt
* 
      SUBROUTINE rotmt1(a,r,k,adot)
      IMPLICIT NONE

      INTEGER k
      DOUBLE PRECISION a,r(3,3),adot

      INTEGER i1,i2,i3
      DOUBLE PRECISION cosa,sina

      IF(k.LT.1 .OR. k.GT.3) STOP '**** rotmt1: k = ??? ****'

      cosa=COS(a)
      sina=SIN(a)
      i1=k
      IF(i1.GT.3) i1=i1-3
      i2=i1+1
      IF(i2.GT.3) i2=i2-3
      i3=i2+1
      IF(i3.GT.3) i3=i3-3

      r(i1,i1) =  0.d0
      r(i1,i2) =  0.d0
      r(i1,i3) =  0.d0
      r(i2,i1) =  0.d0
      r(i2,i2) = -sina*adot
      r(i2,i3) =  cosa*adot
      r(i3,i1) =  0.d0
      r(i3,i2) = -cosa*adot
      r(i3,i3) = -sina*adot

      END
* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 10, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R O T M T 2                           *
*  *                                                               *
*  *           Second time derivative of a rotation matrix         *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    A         -  Rotation angle alpha (rad)
*           K         -  (=1,2,3) Rotation axis (1=x, 2=y, 3=z)
*           ADOT      -  First time derivative of the rotation angle
*                        dalpha/dt (rad/arbitrary time unit)
*           ADDOT     -  Second time derivative of the rotation angle
*                        d^2 alpha/dt^2 (rad/(arbitrary time unit)^2)
*
* OUTPUT:   R         -  Second time derivative of the rotation matrix
*                        d^2 R_k/dt^2 (rad/(arbitrary time unit)^2)
*
* The time derivative of the rotation matrix R_k(alpha) (for its
* definition, see SUBROUTINE ROTMT) is computed according to the formula
* 
* dR^2_k/dt^2 =
*  = d^2 R_k/dalpha^2 * (dalpha/dt)^2 + (dR_k/dalpha) * (d^2 alpha/dt^2)
* 
      SUBROUTINE rotmt2(a,r,k,adot,addot)
      IMPLICIT NONE

      INTEGER k
      DOUBLE PRECISION a,r(3,3),adot,addot

      INTEGER i1,i2,i3
      DOUBLE PRECISION cosa,sina,adot2

      IF(k.LT.1 .OR. k.GT.3) STOP '**** rotmt2: k = ??? ****'

      cosa=COS(a)
      sina=SIN(a)
      adot2=adot**2

      i1=k
      IF(i1.GT.3) i1=i1-3
      i2=i1+1
      IF(i2.GT.3) i2=i2-3
      i3=i2+1
      IF(i3.GT.3) i3=i3-3

      r(i1,i1) =  0.d0
      r(i1,i2) =  0.d0
      r(i1,i3) =  0.d0
      r(i2,i1) =  0.d0
      r(i2,i2) = -cosa*adot2 -sina*addot
      r(i2,i3) = -sina*adot2 +cosa*addot
      r(i3,i1) =  0.d0
      r(i3,i2) =  sina*adot2 -cosa*addot
      r(i3,i3) = -cosa*adot2 -sina*addot

      END

* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 10, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         S U M M A T                           *
*  *                                                               *
*  *                     Sum of two matrices                       *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    B         -  B matrix
*           C         -  C matrix
*
* OUTPUT:   A         -  A = B + C
*
      SUBROUTINE summat(a,b,c)
      IMPLICIT NONE

      DOUBLE PRECISION a(3,3),b(3,3),c(3,3)

      INTEGER i,j

      DO 1 i=1,3
      DO 1 j=1,3
      a(i,j)=b(i,j)+c(i,j)
 1    CONTINUE

      END
* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 10, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                           S U M V                             *
*  *                                                               *
*  *             Sum of two vectors (Y = X1 + X2)                  *
*  *                                                               *
*  *****************************************************************
*
      SUBROUTINE sumv(y,x1,x2)
      IMPLICIT NONE

      DOUBLE PRECISION y(3),x1(3),x2(3)

      INTEGER i

      DO 1 i=1,3
      y(i)=x1(i)+x2(i)
 1    CONTINUE

      END
* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 10, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          T R S P 3                            *
*  *                                                               *
*  *                     Transpose of a matrix                     *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    R         -  Input matrix
*
* OUTPUT:   R         -  R' (transpose of input matrix)
*
      SUBROUTINE trsp3(r)
      IMPLICIT NONE

      DOUBLE PRECISION r(3,3)

      INTEGER i,j
      DOUBLE PRECISION rt

      DO 1 i=1,2
      DO 1 j=i+1,3
      rt=r(i,j)
      r(i,j)=r(j,i)
      r(j,i)=rt
 1    CONTINUE

      END
c ===================================
c HEAPSORT
c sorting by heapsort of a double precision array a of lenght n
c the input array a is not changed, in output
c ind is the indirect address of the sorted array
c ===================================
      SUBROUTINE heapsort(a,n,ind)
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION a(n)
      INTEGER ind(n)
c =====end interface========
      integer j,l,ir,indt,i
      double precision q
c initialise indexes to natural order
      DO j=1,n
        ind(j)=j
      ENDDO
      IF(n.eq.1)RETURN
c counters for recursion, length of heap
      l=n/2+1
      ir=n
c recursive loop
 1    CONTINUE
      IF(l.gt.1)THEN
         l=l-1
         indt=ind(l)
         q=a(indt)
      ELSE
         indt=ind(ir)
         q=a(indt)
         ind(ir)=ind(1)
         ir=ir-1
         IF(ir.eq.1)THEN
            ind(1)=indt
            RETURN
         ENDIF
      ENDIF
      i=l
      j=l+l
 2    IF(j.le.ir)THEN
         IF(j.lt.ir)THEN
            IF(a(ind(j)).lt.a(ind(j+1)))j=j+1
         ENDIF
         IF(q.lt.a(ind(j)))THEN
            ind(i)=ind(j)
            i=j
            j=j+j
         ELSE
            j=ir+1
         ENDIF
         GOTO 2
      ENDIF
      ind(i)=indt
      GOTO 1
      END
c ===================================
c HEAPSORTI
c sorting by heapsort of an integer array a of lenght n
c the input array a is not changed, in output
c ind is the indirect address of the sorted array
c ===================================
      SUBROUTINE heapsorti(a,n,ind)
      IMPLICIT NONE
      INTEGER n
      INTEGER a(n)
      INTEGER ind(n)
c =====end interface========
      integer j,l,ir,indt,i
      double precision q
c initialise indexes to natural order
      DO j=1,n
        ind(j)=j
      ENDDO
      IF(n.eq.1)RETURN
c counters for recursion, length of heap
      l=n/2+1
      ir=n
c recursive loop
 1    CONTINUE
      IF(l.gt.1)THEN
         l=l-1
         indt=ind(l)
         q=a(indt)
      ELSE
         indt=ind(ir)
         q=a(indt)
         ind(ir)=ind(1)
         ir=ir-1
         IF(ir.eq.1)THEN
            ind(1)=indt
            RETURN
         ENDIF
      ENDIF
      i=l
      j=l+l
 2    IF(j.le.ir)THEN
         IF(j.lt.ir)THEN
            IF(a(ind(j)).lt.a(ind(j+1)))j=j+1
         ENDIF
         IF(q.lt.a(ind(j)))THEN
            ind(i)=ind(j)
            i=j
            j=j+j
         ELSE
            j=ir+1
         ENDIF
         GOTO 2
      ENDIF
      ind(i)=indt
      GOTO 1
      END
c ======================================================
c   {\bf roff}: machine accuracy tested in a machine independent way
c
c    warning: if multiple floating point operations use temporary storage
c    in the floating point unit, the result will indicate an 
c    illusory precision. This happens e.g. wtih f2c and g77.
c    fix: roff less than 2e-16 is refused
c
c    warning: roff is only the 1/2 the value of the last stored bit;
c    rounding off in multiple floating point operations should
c    average to a much smaller value in a computer with IEEE 754
c    floating point unit (see Knuth, chap. 4).
c   
c    in output nb-1 is the number of bits in the mantissa
c     (in base 2; in base 16 this is not that clear)
c ======================================================
      double precision function roff(nb)
      double precision roffp
      integer nbp
      logical first
      save first,roffp,nbp
      data first/.true./

      if(first) then
          roff=1.d0
          do 1 nb=1,200
            roff=roff/2.d0
            if(1.d0+roff.eq.1.d0)goto 2
 1        continue
 2        if(roff.lt.2.2d-16)roff=2.2d-16
          roffp=roff
          nbp=nb
          first=.false.
      else
          nb=nbp
          roff=roffp
      endif
      end
c This function returns the truncated integer of the input. 
c Importantly it rounds to the nearest integer if the input 
c is within machine precision of an integer value
      integer function truncat(flt,eps)
      implicit none

      double precision flt,one,eps

      one=1.d0-eps
      truncat=flt
      if(abs(truncat-flt).ge.one) truncat=nint(flt)
      end
