      subroutine  zscal(n,za,zx,incx) 
!                                                                       
!     scales a vector by a constant.                                    
!     jack dongarra, 3/11/78.                                           
!     modified 3/93 to return if incx .le. 0.                           
!     modified 12/3/93, array(1) declarations changed to array(*)       
!                                                                       
      complex*16 za,zx(*) 
      integer i,incx,ix,n 
!                                                                       
      if( n.le.0 .or. incx.le.0 )return 
      if(incx.eq.1)go to 20 
!                                                                       
!        code for increment not equal to 1                              
!                                                                       
      ix = 1 
      do 10 i = 1,n 
        zx(ix) = za*zx(ix) 
        ix = ix + incx 
   10 continue 
      return 
!                                                                       
!        code for increment equal to 1                                  
!                                                                       
   20 do 30 i = 1,n 
        zx(i) = za*zx(i) 
   30 continue 
      return 
      END                                           
