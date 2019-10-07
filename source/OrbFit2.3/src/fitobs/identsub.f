c ===================================================================== 
c recomputation of the suggested identification,
c given the covariance and normal matrices for both arcs at the same time
c ===================================================================== 
      subroutine idno5(eq0,eqp,g1,c1,g2,c2,
     +       d2,da2,dq,dqalt,dista,detc2,detc5,eqf,ff,aval5)
      implicit none
c ==============INPUT==========================
c observation times, orbital elements, their covariance and normal matrices
      double precision eq0(6),eqp(6),aval5(5)
      double precision g1(6,6),g2(6,6),c1(6,6),c2(6,6)
c =============OUTPUT==========================
      double precision eqf(6)
      double precision d2,da2,dista
      double precision dq,dqalt,ddel,detc2,detc5
c logical failure flag
      logical ff
c =============END INTERFACE===================
c 2x2 reduced matrices (for inclination)
      double precision ci1(2,2),ci2(2,2)
      double precision gami(2,2),detgami
      double precision ceq0red(2),pq(2),fact
      double precision pl,primea
c 5x5 matrices
      double precision gamma55(5,5),c15(5,5),c25(5,5)
      double precision aval15(5),aval25(5),det15,det25
      integer ier15,ier25
      integer icode
c loop indexes
      integer i,j
c =============================================
      fact=1.d6
      call mesura(eq0,eqp,fact,6,dista)
c =====================================================================
c Arc 1: computation of reduced normal matrices
c Definition of the covariance matrix Gamma_I
      DO  i=1,2
        DO j=1,2
          gami(i,j)=g1(i+3,j+3)
        ENDDO
      ENDDO
c Computation of the normal matrix C^I
      CALL inv22(gami,ci1,detgami)
c computation of C0red*X0red 
      call mulmav(ci1,2,2,eq0(4),2,ceq0red)
c =====================================================================
c Arc 2: computation of reduced normal matrices
c Definition of the covariance matrix Gamma_I
      DO  i=1,2
        DO j=1,2
          gami(i,j)=g2(i+3,j+3)
        ENDDO
      ENDDO
c Computation of the normal matrix C^I
      CALL inv22(gami,ci2,detgami)
c ===================================================================== 
c  identification according to inclination (only)
      CALL astid2(ci1,eq0(4),ceq0red,ci2,eqp(4),
     +       pq,d2,da2,detc2,icode)
      if(icode.eq.0)then
         write(*,*)' icode=',icode
         d2=0.d0
         da2=0.d0
      endif
c ======================================================================
c compute C1^55,C2^55
       do i=1,5
         do j=1,5
             gamma55(i,j)=g1(i,j)
         enddo
       enddo
       call qrinv5(gamma55,c15,ier15,det15,aval15)
       write(*,*)aval15(1),aval15(5)/aval15(1)
       do i=1,5
         do j=1,5
             gamma55(i,j)=g2(i,j)
         enddo
       enddo
       call qrinv5(gamma55,c25,ier25,det25,aval25)
       write(*,*)aval25(1),aval25(5)/aval25(1)
c now take mean 
        pl=primea(eq0(6),eqp(6))
c now tests difference
        ff=.false.
        call astid5(c15,eq0,c25,eqp,eqf,dq,dqalt,detc5,icode,ddel,aval5)
        eqf(6)=pl
c check if the inversion of c0 failed
        if (icode.eq.0)then
           DO i=1,5
             eqf(i)=0.d0
           ENDDO
           ff=.true.
c check for hyperbolic orbits
        elseif(eqf(1).lt.0.d0.or.eqf(2)**2+eqf(3)**2.gt.1.d0)then
           ff=.true.
        endif
c check for retrograde orbits
        if(eqf(4)**2+eqf(5)**2.gt.1.d0)then
           ff=.true.
        endif

      RETURN
      END
c ===================================================================== 
c recomputation of the suggested identification,
c given the covariance and normal matrices for both arcs at the same time
c ===================================================================== 
      subroutine idno6(eq0,eqp,g1,c1,g2,c2,
     +             d2,da2,dq,dqalt,dista,detc2,detc6,eqf,ff,aval6)
      implicit none
c ==============INPUT==========================
c observation times, orbital elements, their covariance and normal matrices
      double precision eq0(6),eqp(6),aval6(6)
      double precision g1(6,6),g2(6,6),c1(6,6),c2(6,6)
c =============OUTPUT==========================
      double precision eqf(6)
      double precision d2,da2,dista
      double precision dq,dqalt,ddel,detc6
c logical failure flag
      logical ff
c =============END INTERFACE===================
c 2x2 reduced matrices (for inclination)
      double precision ci1(2,2),ci2(2,2)
      double precision gami(2,2),detgami
      double precision ceq0red(2),pq(2),detc2,fact
      integer icode
c loop indexes
      integer i,j
c =============================================
      fact=1.d6
      call mesura(eq0,eqp,fact,6,dista)
c =====================================================================
c Arc 1: computation of reduced normal matrices
c Definition of the covariance matrix Gamma_I
      DO  i=1,2
        DO j=1,2
          gami(i,j)=g1(i+3,j+3)
        ENDDO
      ENDDO
c Computation of the normal matrix C^I
      CALL inv22(gami,ci1,detgami)
c computation of C0red*X0red 
      call mulmav(ci1,2,2,eq0(4),2,ceq0red)
c =====================================================================
c Arc 2: computation of reduced normal matrices
c Definition of the covariance matrix Gamma_I
      DO  i=1,2
        DO j=1,2
          gami(i,j)=g2(i+3,j+3)
        ENDDO
      ENDDO
c Computation of the normal matrix C^I
      CALL inv22(gami,ci2,detgami)
c ===================================================================== 
c  identification according to inclination (only)
c     f=2.0d0/(2.d0*mall)
      CALL astid2(ci1,eq0(4),ceq0red,ci2,eqp(4),
     +       pq,d2,da2,detc2,icode)
      if(icode.eq.0)then
         write(*,*)' icode=',icode
         d2=0.d0
         da2=0.d0
      endif
c now tests difference
        ff=.false.
        call astid6(c1,eq0,c2,eqp,eqf,dq,dqalt,detc6,icode,
     +      ddel,aval6)
c check if the inversion of c0 failed
        if (icode.eq.0)then
           ff=.true.
           dq=0.d0
           dqalt=0.d0
           DO i=1,6
             eqf(i)=0.d0
           ENDDO
c check for hyperbolic orbits
        elseif(eqf(1).lt.0.d0.or.eqf(2)**2+eqf(3)**2.gt.1.d0)then
           ff=.true.
        endif
c check for retrograde orbits
        if(eqf(4)**2+eqf(5)**2.gt.1.d0)then
           ff=.true.
        endif
      RETURN
      END
c ===================================================
c   subroutine astid (linear identification)
c   input:
c         c1 normal matrix 1
c         eq1 equinoctal elements 1
c         ceq1 c1*eq1
c         c2 normal matrix 2
c         eq2 equinoctal elements 2
c         num matrix dimension
c   output:
c         eq0 new first guess
c         deltaq increase in target function (to be mult by 2/m)
c         deltalt the same, alternative computation
c         detc determinant of C
c===============================================================
c   subroutine astid2: to perform identification paying attention 
c   to inclination (only)
c===========================
      subroutine astid2(c1,eq1,ceq1,c2,eq2,eq0,deltaq,
     +    deltalt,detc,icode)
      implicit double precision (a-h,o-z)
      parameter (num=2)
      dimension c1(num,num),c2(num,num),c(num,num),c0(num,num),
     +         g0(num,num),c02(num,num),
     +         calt(num,num),c01(num,num)
      dimension eq1(num),eq2(num),eq0(num),
     +         delta(num),ceq1(num),ceq2(num),cn(num)
c  icode is 1 if we can compute the inverse c0 and is 0, otherwise
      icode=1
c  identification
c  c0=c1+c2
          do 30 i1=1,num
            do 31 i2=1,num
              c0(i1,i2)=c1(i1,i2)+c2(i1,i2)
 31         continue
 30       continue
c normalisation of columns of matrix to be inverted
           do 87 i=1,num
             cn(i)=sqrt(abs(c0(i,i)))
 87        continue
           do 88 i=1,num
             do 89 j=1,num
               c0(i,j)=c0(i,j)/(cn(i)*cn(j))
 89          continue
 88        continue
c  computation of the determinant of c0
          detc0=c0(1,1)*c0(2,2)-c0(1,2)*c0(2,1)
c  if invertible, computation of the inverse of c0
          c0max=max(c0(1,1),c0(2,2),c0(1,2))
          eps=1.d-15*c0max
          if(abs(detc0).lt.eps)then
             icode=0
             return
          endif
          if(detc0.lt.0.d0)then
             write(9,109) detc0
 109         format('detc0 non positive ',f12.3)
          endif
          g0(1,1)=c0(2,2)/detc0
          g0(2,2)=c0(1,1)/detc0         
          g0(1,2)=-c0(2,1)/detc0
          g0(2,1)=g0(1,2)
c unnormalize the matrix by norm of columns
              do 83 i=1,num
                do 84 j=1,num
                  g0(i,j)=g0(i,j)/(cn(i)*cn(j))
 84             continue
 83           continue
c  computation of  C2*X2
          do 1 i=1,2
             ceq2(i)=c2(i,1)*eq2(1)+c2(i,2)*eq2(2)
 1        continue
c  computation of X0
          do 2 i=1,2
            eq0(i)=g0(i,1)*(ceq1(1)+ceq2(1))+g0(i,2)*(ceq1(2)+ceq2(2))
 2        continue
c  matrix c computed by first formula
          do 3 i=1,2
            do 4 j=1,2
              c02(i,j)=g0(i,1)*c2(1,j)+g0(i,2)*c2(2,j)
 4         continue
 3      continue
          do 41 i1=1,num
            do 42 i2=1,num
              c(i1,i2)=c2(i1,i2)-c2(i1,1)*c02(1,i2)-c2(i1,2)*c02(2,i2)
 42         continue
 41       continue
c  matrix c computed by second formula
          do 13 i=1,2
            do 14 j=1,2
              c01(i,j)=g0(i,1)*c1(1,j)+g0(i,2)*c1(2,j)
 14        continue
 13     continue
          do 43 i1=1,num
            do 44 i2=1,num
              calt(i1,i2)=c1(i1,i2)-c1(i1,1)*c01(1,i2)
     +                  -c1(i1,2)*c01(2,i2)
 44        continue
 43     continue
          do 40 i1=1,num
            delta(i1)=eq2(i1)-eq1(i1)
 40       continue
c  determinant
          detc=c(1,1)*c(2,2)-c(2,1)*c(1,2)
c  increase of target function
          deltaq=bilin(delta,delta,num,c,num,num)
          deltalt=bilin(delta,delta,num,calt,num,num)
          ddel=2.d0*(deltaq-deltalt)/(deltaq+deltalt)
          delcr=1.d-1
c          if(abs(ddel).gt.delcr)then
c             write(*,*)'deltaq not consistent',ddel,nam1,nam2
c          endif
      return
      end
c ===================================================
c   subroutine astid (linear identification)
c   input:
c         c1 normal matrix 1
c         eq1 equinoctal elements 1
c         ceq1 c1*eq1
c         c2 normal matrix 2
c         eq2 equinoctal elements 2
c         num matrix dimension
c   output:
c         eq0 new first guess
c         deltaq increase in target function (to be mult by 2/m)
c         deltalt the same, alternative computation
c         detc determinant of C 
c================================================================
c  subroutine astid5: to perform identifications taking into account
c  all the elements
c===================================
      subroutine astid5(c1,eq1,c2,eq2,eq0,deltaq,
     +    deltalt,detc,icode,ddel,aval)
      implicit none
      double precision deltaq,deltalt,ddel
      double precision bilin
      integer num
      parameter (num=5)
      double precision c1(num,num),c2(num,num),c(num,num),
     +         g0(num,num),c02(num,num),cg0(num,num),
     +         calt(num,num),c01(num,num)
      double precision eq1(num),eq2(num),eq0(num),
     +      aval(num), delta(num),ceq1(num),ceq2(num),ceqt(num)
      double precision q(num,num),eigen(num),w1(num),w2(num),detc,detc0
      integer ierr,i,icode,i1,i2
      save
c  icode is 1 if we can compute the inverse c0 and is 0, otherwise
      icode=1
c  identification
      do i1=1,num
        do i2=1,num
           g0(i1,i2)=c1(i1,i2)+c2(i1,i2)
        enddo
      enddo
c ===============================================================
c Covariance matrix through qrinv
      call mcopy(num,num,g0,cg0)
      call qrinv5(cg0,g0,ierr,detc0,aval)
c Non invertible cases are not handle
      if(ierr.ne.0)then
        icode=0
        return
      endif
c  matrix c computed by first formula
      call mulmat(g0,num,num,c2,num,num,c02)
      call mulmat(c2,num,num,c02,num,num,c)
      do  i1=1,num
        do  i2=1,num
          c(i1,i2)=c2(i1,i2)-c(i1,i2)
        enddo
      enddo
c  matrix c computed by second formula
      call mulmat(g0,num,num,c1,num,num,c01)
      call mulmat(c1,num,num,c01,num,num,calt)
      do  i1=1,num
        do i2=1,num
          calt(i1,i2)=(c1(i1,i2)-calt(i1,i2))
        enddo
      enddo
c  computation of  C1*X1
      call mulmav(c1,num,num,eq1,num,ceq1)
c  computation of  C2*X2
      call mulmav(c2,num,num,eq2,num,ceq2)
c  computation of X0
      call vsumg(num,ceq1,ceq2,ceqt)
      call mulmav(g0,num,num,ceqt,num,eq0)
c  compute amount of increase of target function
      do  i1=1,num
        delta(i1)=eq2(i1)-eq1(i1)
      enddo
c  increase of target function
      deltaq=bilin(delta,delta,num,c,num,num)
      deltalt=bilin(delta,delta,num,calt,num,num)
      ddel=2.d0*(deltaq-deltalt)/(deltaq+deltalt)
c determinant of c
      CALL rs(num,num,c,eigen,1,q,w1,w2,ierr)
      detc=1.d0
      DO i=1,num
        detc=detc*eigen(i)
      ENDDO 
      return
      end
c ===================================================
c   subroutine astid (linear identification)
c   input:
c         c1 normal matrix 1
c         eq1 equinoctal elements 1
c         ceq1 c1*eq1
c         c2 normal matrix 2
c         eq2 equinoctal elements 2
c         num matrix dimension
c   output:
c         eq0 new first guess
c         deltaq increase in target function (to be mult by 2/m)
c         deltalt the same, alternative computation
c         detc determinant of C 
c================================================================
c  subroutine astid6: to perform identifications taking into account
c  all the elements
c===================================
      subroutine astid6(c1,eq1,c2,eq2,eq0,deltaq,
     +    deltalt,detc,icode,ddel,aval)
      implicit none
      double precision deltaq,deltalt,ddel
      double precision bilin
      integer num
      parameter (num=6)
      double precision c1(num,num),c2(num,num),c(num,num),
     +         g0(num,num),c02(num,num),cg0(num,num),
     +         calt(num,num),c01(num,num)
      double precision eq1(num),eq2(num),eq0(num),
     +      aval(num), delta(num),ceq1(num),ceq2(num),ceqt(num)
      double precision q(num,num),eigen(num),w1(num),w2(num),detc,detc0
      integer ierr,i,icode,i1,i2
      save
c  icode is 1 if we can compute the inverse c0 and is 0, otherwise
      icode=1
c  identification
      do i1=1,num
        do i2=1,num
           g0(i1,i2)=c1(i1,i2)+c2(i1,i2)
        enddo
      enddo
c ===============================================================
c Covariance matrix through qrinv
      call mcopy(num,num,g0,cg0)
      call qrinv(cg0,g0,ierr,detc0,aval)
c Non invertible cases are not handle
      if(ierr.ne.0)then
        icode=0
        return
      endif
c  matrix c computed by first formula
      call mulmat(g0,num,num,c2,num,num,c02)
      call mulmat(c2,num,num,c02,num,num,c)
      do  i1=1,num
        do  i2=1,num
          c(i1,i2)=c2(i1,i2)-c(i1,i2)
        enddo
      enddo
c  matrix c computed by second formula
      call mulmat(g0,num,num,c1,num,num,c01)
      call mulmat(c1,num,num,c01,num,num,calt)
      do  i1=1,num
        do i2=1,num
          calt(i1,i2)=(c1(i1,i2)-calt(i1,i2))
        enddo
      enddo
c  computation of  C1*X1
      call mulmav(c1,num,num,eq1,num,ceq1)
c  computation of  C2*X2
      call mulmav(c2,num,num,eq2,num,ceq2)
c  computation of X0
      call vsumg(num,ceq1,ceq2,ceqt)
      call mulmav(g0,num,num,ceqt,num,eq0)
c  compute amount of increase of target function
      do  i1=1,num
        delta(i1)=eq2(i1)-eq1(i1)
      enddo
c  increase of target function
      deltaq=bilin(delta,delta,num,c,num,num)
      deltalt=bilin(delta,delta,num,calt,num,num)
      ddel=2.d0*(deltaq-deltalt)/(deltaq+deltalt)
c determinant of c
      CALL rs(num,num,c,eigen,1,q,w1,w2,ierr)
      detc=1.d0
      DO i=1,num
        detc=detc*eigen(i)
      ENDDO 
      return
      end
c ***********************************
      subroutine mesura(eq1,eq2,fact,num,dista)
      implicit double precision (a-h,o-z)
      double precision eq1(num),eq2(num)
      if(num.eq.6)then
         dista = (eq1(2)-eq2(2))**2+(eq1(3)-eq2(3))**2+
     +           (eq1(4)-eq2(4))**2+(eq1(5)-eq2(5))**2+
     +           ((eq1(1)-eq2(1))/(eq1(1)+eq2(1)))**2
         diff = abs(eq1(6)-eq2(6))
         pig=4.d0*atan(1.d0)
         dpig=pig*2.d0
         if (diff.gt.pig)then
           diff = dpig - diff
         endif
         dista = dista + (diff/fact)**2
         dista = sqrt(dista)
      elseif(num.eq.2)then
         dista = (eq1(1)-eq2(1))**2+(eq1(2)-eq2(2))**2
         dista = sqrt(dista)
      elseif(num.eq.5)then
         dista = (eq1(2)-eq2(2))**2+(eq1(3)-eq2(3))**2+
     +          (eq1(4)-eq2(4))**2+(eq1(5)-eq2(5))**2+
     +          ((eq1(1)-eq2(1))/(eq1(1)+eq2(1)))**2
      else
         write(*,*) 'Error in dimension num'
         stop
      endif
      return
      end
c====================================================
      subroutine qrinv(a,b,ierr,det,aval)
      implicit none
c*****************************************************
c a = input matrix
c b = output matrix
c aval= eigenvalues vector
c q = eigenvectors matrix
c*****************************************************
      double precision a(6,6),b(6,6),q(6,6),dstar(6,6)
      double precision qt(6,6),tmp(6,6)
      double precision aval(6),w1(6),w2(6)
      double precision det
      integer i,j,ierr,izer
      call rs(6,6,a,aval,1,q,w1,w2,ierr)
      do 50 i=1,6
        do 60 j=1,6
           dstar(i,j)=0.d0
 60     continue
 50   continue
      det=1.d0
      izer=0
      do 70 i=1,6
         det=det*aval(i)  
         if (aval(i).gt.0.d0) then
           dstar(i,i)=1.d0/aval(i)
         else
           izer=izer+1
         endif
 70   continue
      if(izer.eq.6)then
         ierr=66
         return
      endif
      call transp(q,6,6,qt)   
      call mulmat(dstar,6,6,qt,6,6,tmp)
      call mulmat(q,6,6,tmp,6,6,b)
      return
      end
c ===================================================
      subroutine qrinv5(a,b,ierr,det,aval)
      implicit none
c*****************************************************
c a = input matrix
c b = output matrix
c aval= eigenvalues vector
c q = eigenvectors matrix
c*****************************************************
      double precision a(5,5),b(5,5),q(5,5),dstar(5,5)
      double precision qt(5,5),tmp(5,5)
      double precision aval(5),w1(5),w2(5)
      double precision det
      integer i,j,ierr
      call rs(5,5,a,aval,1,q,w1,w2,ierr)
      do 50 i=1,5
        do 60 j=1,5
           dstar(i,j)=0.d0
 60     continue
 50   continue
      det=1.d0
      do 70 i=1,5
         det=det*aval(i)  
         if (aval(i).gt.0.d0) then
           dstar(i,i)=1.d0/aval(i)
         endif
 70   continue
      call transp(q,5,5,qt)   
      call mulmat(dstar,5,5,qt,5,5,tmp)
      call mulmat(q,5,5,tmp,5,5,b)
      return
      end
c ====================================================
      subroutine qrinv4(a,b,ierr,det,aval)
      implicit none
c*****************************************************
c a = input matrix
c b = output matrix
c aval= eigenvalues vector
c q = eigenvectors matrix
c*****************************************************
      double precision a(4,4),b(4,4),q(4,4),dstar(4,4)
      double precision qt(4,4),tmp(4,4)
      double precision aval(4),w1(4),w2(4)
      double precision det
      integer i,j,ierr
      call rs(4,4,a,aval,1,q,w1,w2,ierr)
      do 50 i=1,4
        do 60 j=1,4
           dstar(i,j)=0.d0
 60     continue
 50   continue
      det=1.d0
      do 70 i=1,4
         det=det*aval(i)  
         if (aval(i).gt.0.d0) then
           dstar(i,i)=1.d0/aval(i)
         endif
 70   continue
      call transp(q,4,4,qt)   
      call mulmat(dstar,4,4,qt,4,4,tmp)
      call mulmat(q,4,4,tmp,4,4,b)
      return
      end

