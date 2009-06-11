! SUBROUTINES: rvfft, irvfft, choosedeg, code_input, decode_out

!------------------------------------------------------------c          
!  A real-valued, in-place, split-radix FFT program          c          
! Real input and output data in array x                      c          
! Length is N=2**M                                           c          
! Decimation-in-time                                         c          
! Output in order                                            c          
!   [Re(0),Re(1),...,Re(N/2),Im(N/2-1),...,Im(1)]            c          
! H.V. Sorensen, Rice University, 1985                       c          
! Modified by D. Bini 1993                                   c          
!------------------------------------------------------------c          
!------------------------------------------------------------c          
! For convolution (polynomial multiplication)                c          
! X must contain the coefficients of the first polynomial    c          
! Y the coefficients of the second polynomial                c          
! Then compute u=irvfft(x), v=irfft(y) and set               c          
!     w(0)=u(0)*v(0), w(i)=u(i)*v(i)-u(n-i+1)*v(n-i+1)       c          
!     w(n-i+1)=u(i)*v(n-i+1)+u(n-i+1)*v(i)                   c          
! Finally output rvfft(w)/n                                  c          
!------------------------------------------------------------c          
                                                                        
                                                                        
subroutine rvfft_QP(x,n,m)
  IMPLICIT NONE
  INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
  REAL(KIND=qkind), INTENT(INOUT) :: x(n)
  INTEGER, INTENT(IN):: n,m
!      implicit double precision (a-h,o-z)
!      dimension x(n) 

  INTEGER :: j,n1,k,is,id,i0,i1,n2,n4,n8,i,i2,i3,i4,i5,i6,i7,i8
  REAL(KIND=qkind) :: xt,r1,e,t1,t2,a,a3,cc1,ss1,cc3,ss3,t3,t4,t5,t6

!----DIGIT REVERSE COUNTER NOT REMOVED -----------------c               
  100  j=1 
       n1=n-1 
       do 104 i=1,n1 
           if(i.ge.j)goto 101 
           xt=x(j) 
           x(j)=x(i) 
           x(i)=xt 
  101      k=n/2 
  102      if (k.ge.j)goto 103 
               j=j-k 
               k=k/2 
               goto 102 
  103      j=j+k 
  104  continue 
!----LENGTH TWO BUTTERFLIES-------------------------c                   
      is=1 
      id=4 
   70 do 60 i0=is,n,id 
        i1=i0+1 
        r1=x(i0) 
        x(i0)=r1+x(i1) 
        x(i1)=r1-x(i1) 
   60 continue 
      is=2*id-1 
      id=4*id 
      if(is.lt.n) goto 70 
!-----L SHAPED BUTTERFLIES--------------------------c                   
      n2=2 
      do 10 k=2,m 
        n2=n2*2 
        n4=n2/4 
        n8=n2/8
! *********************************
!        e=6.283185307179586D0/n2 
        e=8.q0*qatan(1.q0)/n2
! *********************************
        is=0 
        id=n2*2 
   40   do 38 i=is,n-1,id 
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
          t1=(x(i3)+x(i4))/qsqrt(2.q0) 
          t2=(x(i3)-x(i4))/qsqrt(2.q0) 
          x(i4)=x(i2)-t1 
          x(i3)=-x(i2)-t1 
          x(i2)=x(i1)-t2 
          x(i1)=x(i1)+t2 
   38 continue 
         is=2*id-n2 
         id=4*id 
      if(is.lt.n)goto 40 
      a=e 
      do 32 j=2,n8 
         a3=3*a 
         cc1=qcos(a) 
         ss1=qsin(a) 
         cc3=qcos(a3) 
         ss3=qsin(a3) 
         a=j*e 
         is=0 
         id=2*n2 
   36    do 30 i=is,n-1,id 
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
   30    continue 
            is=2*id-n2 
            id=4*id 
         if(is.lt.n)goto 36 
   32    continue 
   10 continue 
      return 
    END SUBROUTINE rvfft_QP

! *********************************************
! DOUBLE PRECISION VERSION
! *********************************************
subroutine rvfft_DP(x,n,m)
  IMPLICIT NONE
  INTEGER, PARAMETER :: dkind=kind(1.d0) !for quadruple precision
  REAL(KIND=dkind), INTENT(INOUT) :: x(n)
  INTEGER, INTENT(IN):: n,m
!      implicit double precision (a-h,o-z)
!      dimension x(n) 

  INTEGER :: j,n1,k,is,id,i0,i1,n2,n4,n8,i,i2,i3,i4,i5,i6,i7,i8
  REAL(KIND=dkind) :: xt,r1,e,t1,t2,a,a3,cc1,ss1,cc3,ss3,t3,t4,t5,t6

!----DIGIT REVERSE COUNTER NOT REMOVED -----------------c               
  100  j=1 
       n1=n-1 
       do 104 i=1,n1 
           if(i.ge.j)goto 101 
           xt=x(j) 
           x(j)=x(i) 
           x(i)=xt 
  101      k=n/2 
  102      if (k.ge.j)goto 103 
               j=j-k 
               k=k/2 
               goto 102 
  103      j=j+k 
  104  continue 
!----LENGTH TWO BUTTERFLIES-------------------------c                   
      is=1 
      id=4 
   70 do 60 i0=is,n,id 
        i1=i0+1 
        r1=x(i0) 
        x(i0)=r1+x(i1) 
        x(i1)=r1-x(i1) 
   60 continue 
      is=2*id-1 
      id=4*id 
      if(is.lt.n) goto 70 
!-----L SHAPED BUTTERFLIES--------------------------c                   
      n2=2 
      do 10 k=2,m 
        n2=n2*2 
        n4=n2/4 
        n8=n2/8
! *********************************
!        e=6.283185307179586D0/n2 
        e=8.d0*atan(1.d0)/n2
! *********************************
        is=0 
        id=n2*2 
   40   do 38 i=is,n-1,id 
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
          t1=(x(i3)+x(i4))/sqrt(2.d0) 
          t2=(x(i3)-x(i4))/sqrt(2.d0) 
          x(i4)=x(i2)-t1 
          x(i3)=-x(i2)-t1 
          x(i2)=x(i1)-t2 
          x(i1)=x(i1)+t2 
   38 continue 
         is=2*id-n2 
         id=4*id 
      if(is.lt.n)goto 40 
      a=e 
      do 32 j=2,n8 
         a3=3*a 
         cc1=cos(a) 
         ss1=sin(a) 
         cc3=cos(a3) 
         ss3=sin(a3) 
         a=j*e 
         is=0 
         id=2*n2 
   36    do 30 i=is,n-1,id 
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
   30    continue 
            is=2*id-n2 
            id=4*id 
         if(is.lt.n)goto 36 
   32    continue 
   10 continue 
      return 
    END SUBROUTINE rvfft_DP

!-------------------------------------------------------c               
!  A real-valued, in-place, split-radix IFFT program    c               
! Hermitian symmetric input and real output in array X  c               
! Length is N=2**M                                      c               
! Decimation-in-frequency                               c               
! Input  order                                          c               
!     [Re(0),Re(1),...,Re(N/2),Im(N/2-1),...,Im(1)]     c               
! H.V. Sorensen, Rice University, 1985                  c               
! Modified by D. Bini 1993                              c               
!-------------------------------------------------------c               
      subroutine irvfft_QP(x,n,m) 
        IMPLICIT NONE
        INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
        REAL(KIND=qkind), INTENT(INOUT) :: x(n)
        INTEGER, INTENT(IN):: n,m
!      implicit double precision (a-h,o-z) 
!      dimension x(n) 
        INTEGER :: n2,k,is,id,n4,n8,i,i1,i2,i3,i4,j,i5,i6,i7,i8,i0,n1
        real(KIND=qkind) :: e,t1,t2,a,a3,cc1,ss1,cc3,ss3,t3,t4,t5,r1,xt
!----L SHAPED BUTTERFLIES-------------------------------c               
      n2=2*n 
      do 10 k=1,m-1 
         is=0 
         id=n2 
         n2=n2/2 
         n4=n2/4 
         n8=n4/2 
! *********************************
!        e=6.283185307179586D0/n2 
        e=8.q0*qatan(1.q0)/n2
! *********************************
   17    do 15 i=is,n-1,id 
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
            t1=(x(i2)-x(i1))/qsqrt(2.q0) 
            t2=(x(i4)+x(i3))/qsqrt(2.q0) 
            x(i1)=x(i1)+x(i2) 
            x(i2)=x(i4)-x(i3) 
            x(i3)=2*(-t2-t1) 
            x(i4)=2*(-t2+t1) 
   15    continue 
            is=2*id-n2 
            id=4*id 
         if(is.lt.n-1)goto 17 
         a=e 
         do 20 j=2,n8 
            a3=3*a 
            cc1=qcos(a) 
            ss1=qsin(a) 
            cc3=qcos(a3) 
            ss3=qsin(a3) 
            a=j*e 
            is=0 
            id=2*n2 
   40       do 30 i=is,n-1,id 
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
   30      continue 
               is=2*id-n2 
               id=4*id 
           if(is.lt.n-1)goto 40 
   20    continue 
   10 continue 
!-----LENGTH TWO BUTTERFLIES--------------------------------c           
      is=1 
      id=4 
   70 do 60 i0=is,n,id 
      i1=i0+1 
      r1=x(i0) 
      x(i0)=r1+x(i1) 
      x(i1)=r1-x(i1) 
   60 continue 
      is=2*id-1 
      id=4*id 
      if(is.lt.n)goto 70 
!----DIGIT REVERSE COUNTER NOT REMOVED -----------------c               
  100  j=1 
       n1=n-1 
       do 104 i=1,n1 
           if(i.ge.j)goto 101 
           xt=x(j) 
           x(j)=x(i) 
           x(i)=xt 
  101      k=n/2 
  102      if (k.ge.j)goto 103 
               j=j-k 
               k=k/2 
               goto 102 
  103      j=j+k 
  104  continue 
                                                                        
!-----DIVISION BY N REMOVED---------------------------------c           
       do 99 i=1,n 
          x(i)=x(i)/n 
   99  continue 
       return 
      END SUBROUTINE irvfft_QP


! *********************************************
! DOUBLE PRECISION VERSION
! ********************************************* 
      subroutine irvfft_DP(x,n,m) 
        IMPLICIT NONE
        INTEGER, PARAMETER :: dkind=kind(1.d0) !for quadruple precision
        REAL(KIND=dkind), INTENT(INOUT) :: x(n)
        INTEGER, INTENT(IN):: n,m
!      implicit double precision (a-h,o-z) 
!      dimension x(n) 
        INTEGER :: n2,k,is,id,n4,n8,i,i1,i2,i3,i4,j,i5,i6,i7,i8,i0,n1
        real(KIND=dkind) :: e,t1,t2,a,a3,cc1,ss1,cc3,ss3,t3,t4,t5,r1,xt
!----L SHAPED BUTTERFLIES-------------------------------c               
      n2=2*n 
      do 10 k=1,m-1 
         is=0 
         id=n2 
         n2=n2/2 
         n4=n2/4 
         n8=n4/2 
! *********************************
!        e=6.283185307179586D0/n2 
        e=8.d0*atan(1.d0)/n2
! *********************************
   17    do 15 i=is,n-1,id 
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
            t1=(x(i2)-x(i1))/sqrt(2.d0) 
            t2=(x(i4)+x(i3))/sqrt(2.d0) 
            x(i1)=x(i1)+x(i2) 
            x(i2)=x(i4)-x(i3) 
            x(i3)=2*(-t2-t1) 
            x(i4)=2*(-t2+t1) 
   15    continue 
            is=2*id-n2 
            id=4*id 
         if(is.lt.n-1)goto 17 
         a=e 
         do 20 j=2,n8 
            a3=3*a 
            cc1=cos(a) 
            ss1=sin(a) 
            cc3=cos(a3) 
            ss3=sin(a3) 
            a=j*e 
            is=0 
            id=2*n2 
   40       do 30 i=is,n-1,id 
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
   30      continue 
               is=2*id-n2 
               id=4*id 
           if(is.lt.n-1)goto 40 
   20    continue 
   10 continue 
!-----LENGTH TWO BUTTERFLIES--------------------------------c           
      is=1 
      id=4 
   70 do 60 i0=is,n,id 
      i1=i0+1 
      r1=x(i0) 
      x(i0)=r1+x(i1) 
      x(i1)=r1-x(i1) 
   60 continue 
      is=2*id-1 
      id=4*id 
      if(is.lt.n)goto 70 
!----DIGIT REVERSE COUNTER NOT REMOVED -----------------c               
  100  j=1 
       n1=n-1 
       do 104 i=1,n1 
           if(i.ge.j)goto 101 
           xt=x(j) 
           x(j)=x(i) 
           x(i)=xt 
  101      k=n/2 
  102      if (k.ge.j)goto 103 
               j=j-k 
               k=k/2 
               goto 102 
  103      j=j+k 
  104  continue 
                                                                        
!-----DIVISION BY N REMOVED---------------------------------c           
       do 99 i=1,n 
          x(i)=x(i)/n 
   99  continue 
       return 
     END SUBROUTINE irvfft_DP

! ============================================================
! ******* CODING of the n evaluations of a polynomial ********
! ********* p: p(om(1)), p(om(n)) to apply irvfft ************
! ******** written by GIOVANNI F. GRONCHI (2004) *************
! ============================================================
 SUBROUTINE code_input_QP(N,eval,codeval)
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: N
   INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
! evaluation of poly p(x) in the N-th roots of unity
   COMPLEX(KIND=qkind),DIMENSION(N),INTENT(IN) :: eval
! -----------------------------------------------------------
! if omega is a primitive complex root of unity then we have
! Re(p(omega^s)) = Re(p(omeg^{-s}))
! Im(p(omega^s)) = -Im(p(omeg^{-s}))   for s=1,..,N/2
! -----------------------------------------------------------
   REAL(KIND=qkind),DIMENSION(N),INTENT(OUT) :: codeval ! coded evaluations
! ------ end interface ----------------
   INTEGER :: i ! loop indexes
! ============================================================

   DO i = 1,N/2+1
      codeval(i) = QREAL(eval(i)) ! from 1 to N/2+1 take the real part
   ENDDO
! (roots of unity are taken clockwise)
   DO i = 1,N/2-1 
      codeval(N/2+1+i) = AIMAG(eval(N/2-i+1))
   ENDDO
 END SUBROUTINE code_input_QP

! *********************************************
! DOUBLE PRECISION VERSION
! ********************************************* 
 SUBROUTINE code_input_DP(N,eval,codeval)
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: N
   INTEGER, PARAMETER :: dkind=kind(1.d0) !for quadruple precision
! evaluation of poly p(x) in the N-th roots of unity
   COMPLEX(KIND=dkind),DIMENSION(N),INTENT(IN) :: eval
! -----------------------------------------------------------
! if omega is a primitive complex root of unity then we have
! Re(p(omega^s)) = Re(p(omeg^{-s}))
! Im(p(omega^s)) = -Im(p(omeg^{-s}))   for s=1,..,N/2
! -----------------------------------------------------------
   REAL(KIND=dkind),DIMENSION(N),INTENT(OUT) :: codeval ! coded evaluations
! ------ end interface ----------------
   INTEGER :: i ! loop indexes
! ============================================================

   DO i = 1,N/2+1
      codeval(i) = REAL(eval(i)) ! from 1 to N/2+1 take the real part
   ENDDO
! (roots of unity are taken clockwise)
   DO i = 1,N/2-1 
      codeval(N/2+1+i) = AIMAG(eval(N/2-i+1))
   ENDDO
 END SUBROUTINE code_input_DP

!     ******************************************************************
!     OUTPUT CONVERSION FROM THE STANDARD FORM of DTF TO THE USUAL ONE  
!     ******************************************************************
!     *********** written by GIOVANNI F. GRONCHI (2001) ****************
!     ********** Department of Mathematics, UNIVERSITY of PISA *********
!     ==================================================================
 SUBROUTINE decode_out_QP(N,dftout1,dftout2,dfteval) 
   IMPLICIT NONE 
   INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
   INTEGER N 
   REAL(KIND=qkind) :: dftout1( * ),dftout2( * ) 
   COMPLEX(KIND=qkind) :: dfteval( * ) 
!     ----------------------------------- end interface ----------------
!     loop indexes                                                      
   INTEGER j 
!     ==================================================================
   dfteval(1) = DCMPLX(dftout1(1),0.d0) 
   DO j = 2,N/2 
      dfteval(j) = DCMPLX(dftout1(j),dftout1(N+2-j)) 
      dfteval(N/2+j) = DCMPLX(dftout2(N/2+j),dftout2(N/2+2-j)) 
   ENDDO
   dfteval(N/2+1) = DCMPLX(dftout1(N/2+1),0.d0) 
   
 END SUBROUTINE decode_out_QP
 
! *********************************************
! DOUBLE PRECISION VERSION
! ********************************************* 
 SUBROUTINE decode_out_DP(N,dftout1,dftout2,dfteval) 
   IMPLICIT NONE 
   INTEGER, PARAMETER :: dkind=kind(1.d0) !for quadruple precision
   INTEGER N 
   REAL(KIND=dkind) :: dftout1( * ),dftout2( * ) 
   COMPLEX(KIND=dkind) :: dfteval( * ) 
!     ----------------------------------- end interface ----------------
!     loop indexes                                                      
   INTEGER j 
!     ==================================================================
   dfteval(1) = DCMPLX(dftout1(1),0.d0) 
   DO j = 2,N/2 
      dfteval(j) = DCMPLX(dftout1(j),dftout1(N+2-j)) 
      dfteval(N/2+j) = DCMPLX(dftout2(N/2+j),dftout2(N/2+2-j)) 
   ENDDO
   dfteval(N/2+1) = DCMPLX(dftout1(N/2+1),0.d0) 
   
 END SUBROUTINE decode_out_DP
 
