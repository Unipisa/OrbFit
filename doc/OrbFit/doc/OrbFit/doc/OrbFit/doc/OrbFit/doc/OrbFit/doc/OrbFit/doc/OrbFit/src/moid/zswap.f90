      subroutine  zswap (n,zx,incx,zy,incy) 
!                                                                       
!     interchanges two vectors.                                         
!     jack dongarra, 3/11/78.                                           
!     modified 12/3/93, array(1) declarations changed to array(*)       
!                                                                       
      complex*16 zx(*),zy(*),ztemp 
      integer i,incx,incy,ix,iy,n 
!                                                                       
      if(n.le.0)return 
      if(incx.eq.1.and.incy.eq.1)go to 20 
!                                                                       
!       code for unequal increments or equal increments not equal       
!         to 1                                                          
!                                                                       
      ix = 1 
      iy = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      if(incy.lt.0)iy = (-n+1)*incy + 1 
      do 10 i = 1,n 
        ztemp = zx(ix) 
        zx(ix) = zy(iy) 
        zy(iy) = ztemp 
        ix = ix + incx 
        iy = iy + incy 
   10 continue 
      return 
!                                                                       
!       code for both increments equal to 1                             
   20 do 30 i = 1,n 
        ztemp = zx(i) 
        zx(i) = zy(i) 
        zy(i) = ztemp 
   30 continue 
      return 
      END                                           
