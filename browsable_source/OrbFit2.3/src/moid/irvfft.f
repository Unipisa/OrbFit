c-------------------------------------------------------c
c  A real-valued, in-place, split-radix IFFT program    c
c Hermitian symmetric input and real output in array X  c
c Length is N=2**M                                      c
c Decimation-in-frequency                               c
c Input  order                                          c
c     [Re(0),Re(1),...,Re(N/2),Im(N/2-1),...,Im(1)]     c
c H.V. Sorensen, Rice University, 1985                  c
c Modified by D. Bini 1993                              c
c-------------------------------------------------------c
      subroutine irvfft(x,n,m)
      implicit real*8 (a-h,o-z)
      dimension x(n)
c----L SHAPED BUTTERFLIES-------------------------------c
      n2=2*n
      do 10 k=1,m-1
         is=0
         id=n2
         n2=n2/2
         n4=n2/4
         n8=n4/2
        e=6.283185307179586D0/n2
17       do 15 i=is,n-1,id
            i1=i+1
            i2=i1+n4
            i3=i2+n4
            i4=i3+n4
            t1=x(i1)-x(i3)
            x(i1)=x(i1)+x(i3)
            x(i2)=2*x(i2)
            x(i3)=t1-2*x(i4)
            x(i4)=t1+2*x(i4)
            if (n4.eq.1)goto 15
            i1=i1+n8
            i2=i2+n8
            i3=i3+n8
            i4=i4+n8
            t1=(x(i2)-x(i1))/dsqrt(2.0d0)
            t2=(x(i4)+x(i3))/dsqrt(2.0d0)
            x(i1)=x(i1)+x(i2)
            x(i2)=x(i4)-x(i3)
            x(i3)=2*(-t2-t1)
            x(i4)=2*(-t2+t1)
15       continue
            is=2*id-n2
            id=4*id
         if(is.lt.n-1)goto 17
         a=e
         do 20 j=2,n8
            a3=3*a
            cc1=dcos(a)
            ss1=dsin(a)
            cc3=dcos(a3)
            ss3=dsin(a3)
            a=j*e
            is=0
            id=2*n2
40          do 30 i=is,n-1,id
               i1=i+j
               i2=i1+n4
               i3=i2+n4
               i4=i3+n4
               i5=i+n4-j+2
               i6=i5+n4
               i7=i6+n4
               i8=i7+n4
               t1=x(i1)-x(i6)
               x(i1)=x(i1)+x(i6)
               t2=x(i5)-x(i2)
               x(i5)=x(i2)+x(i5)
               t3=x(i8)+x(i3)
               x(i6)=x(i8)-x(i3)
               t4=x(i4)+x(i7)
               x(i2)=x(i4)-x(i7)
               t5=t1-t4
               t1=t1+t4
               t4=t2-t3
               t2=t2+t3
               x(i3)=t5*cc1+t4*ss1
               x(i7)=-t4*cc1+t5*ss1
               x(i4)=t1*cc3-t2*ss3
               x(i8)=t2*cc3+t1*ss3
30         continue
               is=2*id-n2
               id=4*id
           if(is.lt.n-1)goto 40
20       continue
10    continue
c-----LENGTH TWO BUTTERFLIES--------------------------------c
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
      if(is.lt.n)goto 70
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

c-----DIVISION BY N REMOVED---------------------------------c
       do 99 i=1,n
          x(i)=x(i)/n
 99    continue
       return
       end

