c------------------------------------------------------------c
c  A real-valued, in-place, split-radix FFT program          c
c Real input and output data in array x                      c
c Length is N=2**M                                           c
c Decimation-in-time                                         c
c Output in order                                            c
c   [Re(0),Re(1),...,Re(N/2),Im(N/2-1),...,Im(1)]            c
c H.V. Sorensen, Rice University, 1985                       c
c Modified by D. Bini 1993                                   c
c------------------------------------------------------------c
c------------------------------------------------------------c
c For convolution (polynomial multiplication)                c
c X must contain the coefficients of the first polynomial    c
c Y the coefficients of the second polynomial                c
c Then compute u=irvfft(x), v=irfft(y) and set               c
c     w(0)=u(0)*v(0), w(i)=u(i)*v(i)-u(n-i+1)*v(n-i+1)       c
c     w(n-i+1)=u(i)*v(n-i+1)+u(n-i+1)*v(i)                   c
c Finally output rvfft(w)/n                                  c
c------------------------------------------------------------c


      subroutine rvfft(x,n,m)
      implicit real*8 (a-h,o-z)
      dimension x(n)
c----DIGIT REVERSE COUNTER NOT REMOVED -----------------c
100    j=1
       n1=n-1
       do 104 i=1,n1
           if(i.ge.j)goto 101
           xt=x(j)
           x(j)=x(i)
           x(i)=xt
101        k=n/2
102        if (k.ge.j)goto 103
               j=j-k
               k=k/2
               goto102
103        j=j+k
104    continue
c----LENGTH TWO BUTTERFLIES-------------------------c
      is=1
      id=4
70    do 60 i0=is,n,id
        i1=i0+1
        r1=x(i0)
        x(i0)=r1+x(i1)
        x(i1)=r1-x(i1)
60    continue
      is=2*id-1
      id=4*id
      if(is.lt.n) goto 70
c-----L SHAPED BUTTERFLIES--------------------------c
      n2=2
      do 10 k=2,m
        n2=n2*2
        n4=n2/4
        n8=n2/8
        e=6.283185307179586D0/n2
        is=0
        id=n2*2
40      do 38 i=is,n-1,id
          i1=i+1
          i2=i1+n4
          i3=i2+n4
          i4=i3+n4
          t1=x(i4)+x(i3)
          x(i4)=x(i4)-x(i3)
          x(i3)=x(i1)-t1
          x(i1)=x(i1)+t1
          if(n4.eq.1)goto 38
          i1=i1+n8
          i2=i2+n8
          i3=i3+n8
          i4=i4+n8
          t1=(x(i3)+x(i4))/dsqrt(2.0d0)
          t2=(x(i3)-x(i4))/dsqrt(2.0d0)
          x(i4)=x(i2)-t1
          x(i3)=-x(i2)-t1
          x(i2)=x(i1)-t2
          x(i1)=x(i1)+t2
38    continue
         is=2*id-n2
         id=4*id
      if(is.lt.n)goto 40
      a=e
      do 32 j=2,n8
         a3=3*a         
         cc1=dcos(a)
         ss1=dsin(a)
         cc3=dcos(a3)
         ss3=dsin(a3)
         a=j*e
         is=0
         id=2*n2
36       do 30 i=is,n-1,id
            i1=i+j
            i2=i1+n4
            i3=i2+n4
            i4=i3+n4
            i5=i+n4-j+2
            i6=i5+n4
            i7=i6+n4
            i8=i7+n4
            t1=x(i3)*cc1+x(i7)*ss1
            t2=x(i7)*cc1-x(i3)*ss1
            t3=x(i4)*cc3+x(i8)*ss3
            t4=x(i8)*cc3-x(i4)*ss3
            t5=t1+t3
            t6=t2+t4
            t3=t1-t3
            t4=t2-t4
            t2=x(i6)+t6
            x(i3)=t6-x(i6)
            x(i8)=t2
            t2=x(i2)-t3
            x(i7)=-x(i2)-t3
            x(i4)=t2
            t1=x(i1)+t5
            x(i6)=x(i1)-t5
            x(i1)=t1
            t1=x(i5)+t4
            x(i5)=x(i5)-t4
            x(i2)=t1
30       continue
            is=2*id-n2
            id=4*id
         if(is.lt.n)goto 36
32       continue
10    continue
      return
      end

