      integer function izamax(n,zx,incx) 
!                                                                       
!     finds the index of element having max. absolute value.            
!     jack dongarra, 1/15/85.                                           
!     modified 3/93 to return if incx .le. 0.                           
!     modified 12/3/93, array(1) declarations changed to array(*)       
!                                                                       
      complex*16 zx(*) 
      double precision smax 
      integer i,incx,ix,n 
      double precision dcabs1 
!                                                                       
      izamax = 0 
      if( n.lt.1 .or. incx.le.0 )return 
      izamax = 1 
      if(n.eq.1)return 
      if(incx.eq.1)go to 20 
!                                                                       
!        code for increment not equal to 1                              
!                                                                       
      ix = 1 
      smax = dcabs1(zx(1)) 
      ix = ix + incx 
      do 10 i = 2,n 
         if(dcabs1(zx(ix)).le.smax) go to 5 
         izamax = i 
         smax = dcabs1(zx(ix)) 
    5    ix = ix + incx 
   10 continue 
      return 
!                                                                       
!        code for increment equal to 1                                  
!                                                                       
   20 smax = dcabs1(zx(1)) 
      do 30 i = 2,n 
         if(dcabs1(zx(i)).le.smax) go to 30 
         izamax = i 
         smax = dcabs1(zx(i)) 
   30 continue 
      return 
      END                                           
