      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)                               
c                                                                               
      integer n,nm,ierr,matz                                                    
      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)                                 
c                                                                               
c     this subroutine calls the recommended sequence of                         
c     subroutines from the eigensystem subroutine package (eispack)             
c     to find the eigenvalues and eigenvectors (if desired)                     
c     of a real symmetric matrix.                                               
c                                                                               
c     on input:                                                                 
c                                                                               
c        nm  must be set to the row dimension of the two-dimensional            
c        array parameters as declared in the calling program                    
c        dimension statement;                                                   
c                                                                               
c        n  is the order of the matrix  a;                                      
c                                                                               
c        a  contains the real symmetric matrix;                                 
c                                                                               
c        matz  is an integer variable set equal to zero if                      
c        only eigenvalues are desired;  otherwise it is set to                  
c        any non-zero integer for both eigenvalues and eigenvectors.            
c                                                                               
c     on output:                                                                
c                                                                               
c        w  contains the eigenvalues in ascending order;                        
c                                                                               
c        z  contains the eigenvectors if matz is not zero;                      
c                                                                               
c        ierr  is an integer output variable set equal to an                    
c        error completion code described in section 2b of the                   
c        documentation.  the normal completion code is zero;                    
c                                                                               
c        fv1  and  fv2  are temporary storage arrays.                           
c                                                                               
c     questions and comments should be directed to b. s. garbow,                
c     applied mathematics division, argonne national laboratory                 
c                                                                               
c     ------------------------------------------------------------------        
c                                                                               
      if (n .le. nm) go to 10                                                   
      ierr = 10 * n                                                             
      go to 50                                                                  
c                                                                               
   10 if (matz .ne. 0) go to 20                                                 
c     :::::::::: find eigenvalues only ::::::::::                               
      call  tred1(nm,n,a,w,fv1,fv2)                                             
      call  tqlrat(n,w,fv2,ierr)                                                
      go to 50                                                                  
c     :::::::::: find both eigenvalues and eigenvectors ::::::::::              
   20 call  tred2(nm,n,a,w,fv1,z)                                               
      call  tql2(nm,n,w,fv1,z,ierr)                                             
   50 return                                                                    
c     :::::::::: last card of rs ::::::::::                                     
      end                                                                       
      subroutine tred1(nm,n,a,d,e,e2)                                           
c                                                                               
      integer i,j,k,l,n,ii,nm,jp1                                               
      double precision a(nm,n),d(n),e(n),e2(n)                                            
      double precision f,g,h,scale                                                        
c                                                                               
c     this subroutine is a translation of the algol procedure tred1,            
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.           
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).           
c                                                                               
c     this subroutine reduces a real symmetric matrix                           
c     to a symmetric tridiagonal matrix using                                   
c     orthogonal similarity transformations.                                    
c                                                                               
c     on input:                                                                 
c                                                                               
c        nm must be set to the row dimension of two-dimensional                 
c          array parameters as declared in the calling program                  
c          dimension statement;                                                 
c                                                                               
c        n is the order of the matrix;                                          
c                                                                               
c        a contains the real symmetric input matrix.  only the                  
c          lower triangle of the matrix need be supplied.                       
c                                                                               
c     on output:                                                                
c                                                                               
c        a contains information about the orthogonal trans-                     
c          formations used in the reduction in its strict lower                 
c          triangle.  the full upper triangle of a is unaltered;                
c                                                                               
c        d contains the diagonal elements of the tridiagonal matrix;            
c                                                                               
c        e contains the subdiagonal elements of the tridiagonal                 
c          matrix in its last n-1 positions.  e(1) is set to zero;              
c                                                                               
c        e2 contains the squares of the corresponding elements of e.            
c          e2 may coincide with e if the squares are not needed.                
c                                                                               
c     questions and comments should be directed to b. s. garbow,                
c     applied mathematics division, argonne national laboratory                 
c                                                                               
c     ------------------------------------------------------------------        
c                                                                               
      do 100 i = 1, n                                                           
  100 d(i) = a(i,i)                                                             
c     :::::::::: for i=n step -1 until 1 do -- ::::::::::                       
      do 300 ii = 1, n                                                          
         i = n + 1 - ii                                                         
         l = i - 1                                                              
         h = 0.0d0                                                              
         scale = 0.0d0                                                          
         if (l .lt. 1) go to 130                                                
c     :::::::::: scale row (algol tol then not needed) ::::::::::               
         do 120 k = 1, l                                                        
  120    scale = scale + abs(a(i,k))                                           
c                                                                               
         if (scale .ne. 0.0d0) go to 140                                        
  130    e(i) = 0.0d0                                                           
         e2(i) = 0.0d0                                                          
         go to 290                                                              
c                                                                               
  140    do 150 k = 1, l                                                        
            a(i,k) = a(i,k) / scale                                             
            h = h + a(i,k) * a(i,k)                                             
  150    continue                                                               
c                                                                               
         e2(i) = scale * scale * h                                              
         f = a(i,l)                                                             
         g = -sign(sqrt(h),f)                                                 
         e(i) = scale * g                                                       
         h = h - f * g                                                          
         a(i,l) = f - g                                                         
         if (l .eq. 1) go to 270                                                
         f = 0.0d0                                                              
c                                                                               
         do 240 j = 1, l                                                        
            g = 0.0d0                                                           
c     :::::::::: form element of a*u ::::::::::                                 
            do 180 k = 1, j                                                     
  180       g = g + a(j,k) * a(i,k)                                             
c                                                                               
            jp1 = j + 1                                                         
            if (l .lt. jp1) go to 220                                           
c                                                                               
            do 200 k = jp1, l                                                   
  200       g = g + a(k,j) * a(i,k)                                             
c     :::::::::: form element of p ::::::::::                                   
  220       e(j) = g / h                                                        
            f = f + e(j) * a(i,j)                                               
  240    continue                                                               
c                                                                               
         h = f / (h + h)                                                        
c     :::::::::: form reduced a ::::::::::                                      
         do 260 j = 1, l                                                        
            f = a(i,j)                                                          
            g = e(j) - h * f                                                    
            e(j) = g                                                            
c                                                                               
            do 260 k = 1, j                                                     
               a(j,k) = a(j,k) - f * e(k) - g * a(i,k)                          
  260    continue                                                               
c                                                                               
  270    do 280 k = 1, l                                                        
  280    a(i,k) = scale * a(i,k)                                                
c                                                                               
  290    h = d(i)                                                               
         d(i) = a(i,i)                                                          
         a(i,i) = h                                                             
  300 continue                                                                  
c                                                                               
      return                                                                    
c     :::::::::: last card of tred1 ::::::::::                                  
      end                                                                       
      subroutine tqlrat(n,d,e2,ierr)                                            
c                                                                               
      integer i,j,l,m,n,ii,l1,mml,ierr                                          
      double precision d(n),e2(n)                                                         
      double precision b,c,f,g,h,p,r,s,machep                                             
*      real*8 dsqrt,dabs,dsign                                                   
c                                                                               
c     this subroutine is a translation of the algol procedure tqlrat,           
c     algorithm 464, comm. acm 16, 689(1973) by reinsch.                        
c                                                                               
c     this subroutine finds the eigenvalues of a symmetric                      
c     tridiagonal matrix by the rational ql method.                             
c                                                                               
c     on input:                                                                 
c                                                                               
c        n is the order of the matrix;                                          
c                                                                               
c        d contains the diagonal elements of the input matrix;                  
c                                                                               
c        e2 contains the squares of the subdiagonal elements of the             
c          input matrix in its last n-1 positions.  e2(1) is arbitrary.         
c                                                                               
c      on output:                                                               
c                                                                               
c        d contains the eigenvalues in ascending order.  if an                  
c          error exit is made, the eigenvalues are correct and                  
c          ordered for indices 1,2,...ierr-1, but may not be                    
c          the smallest eigenvalues;                                            
c                                                                               
c        e2 has been destroyed;                                                 
c                                                                               
c        ierr is set to                                                         
c          zero       for normal return,                                        
c          j          if the j-th eigenvalue has not been                       
c                     determined after 30 iterations.                           
c                                                                               
c     questions and comments should be directed to b. s. garbow,                
c     applied mathematics division, argonne national laboratory                 
c                                                                               
c     ------------------------------------------------------------------        
c                                                                               
c     :::::::::: machep is a machine dependent parameter specifying             
c                the relative precision of floating point arithmetic.           
c                machep = 16.0d0**(-13) for long form arithmetic                
c                on s360 ::::::::::                                             
*      data machep/z'3410000000000000'/
      data machep/1.d-16/
c                                                                               
      ierr = 0                                                                  
      if (n .eq. 1) go to 1001                                                  
c                                                                               
      do 100 i = 2, n                                                           
  100 e2(i-1) = e2(i)                                                           
c                                                                               
      f = 0.0d0                                                                 
      b = 0.0d0                                                                 
      e2(n) = 0.0d0                                                             
c                                                                               
      do 290 l = 1, n                                                           
         j = 0                                                                  
         h = machep * (abs(d(l)) + sqrt(e2(l)))                               
         if (b .gt. h) go to 105                                                
         b = h                                                                  
         c = b * b                                                              
c     :::::::::: look for small squared sub-diagonal element ::::::::::         
  105    do 110 m = l, n                                                        
            if (e2(m) .le. c) go to 120                                         
c     :::::::::: e2(n) is always zero, so there is no exit                      
c                through the bottom of the loop ::::::::::                      
  110    continue                                                               
c                                                                               
  120    if (m .eq. l) go to 210                                                
  130    if (j .eq. 30) go to 1000                                              
         j = j + 1                                                              
c     :::::::::: form shift ::::::::::                                          
         l1 = l + 1                                                             
         s = sqrt(e2(l))                                                       
         g = d(l)                                                               
         p = (d(l1) - g) / (2.0d0 * s)                                          
         r = sqrt(p*p+1.0d0)                                                   
         d(l) = s / (p + sign(r,p))                                            
         h = g - d(l)                                                           
c                                                                               
         do 140 i = l1, n                                                       
  140    d(i) = d(i) - h                                                        
c                                                                               
         f = f + h                                                              
c     :::::::::: rational ql transformation ::::::::::                          
         g = d(m)                                                               
         if (g .eq. 0.0d0) g = b                                                
         h = g                                                                  
         s = 0.0d0                                                              
         mml = m - l                                                            
c     :::::::::: for i=m-1 step -1 until l do -- ::::::::::                     
         do 200 ii = 1, mml                                                     
            i = m - ii                                                          
            p = g * h                                                           
            r = p + e2(i)                                                       
            e2(i+1) = s * r                                                     
            s = e2(i) / r                                                       
            d(i+1) = h + s * (h + d(i))                                         
            g = d(i) - e2(i) / g                                                
            if (g .eq. 0.0d0) g = b                                             
            h = g * p / r                                                       
  200    continue                                                               
c                                                                               
         e2(l) = s * g                                                          
         d(l) = h                                                               
c     :::::::::: guard against underflow in convergence test ::::::::::         
         if (h .eq. 0.0d0) go to 210                                            
         if (abs(e2(l)) .le. abs(c/h)) go to 210                              
         e2(l) = h * e2(l)                                                      
         if (e2(l) .ne. 0.0d0) go to 130                                        
  210    p = d(l) + f                                                           
c     :::::::::: order eigenvalues ::::::::::                                   
         if (l .eq. 1) go to 250                                                
c     :::::::::: for i=l step -1 until 2 do -- ::::::::::                       
         do 230 ii = 2, l                                                       
            i = l + 2 - ii                                                      
            if (p .ge. d(i-1)) go to 270                                        
            d(i) = d(i-1)                                                       
  230    continue                                                               
c                                                                               
  250    i = 1                                                                  
  270    d(i) = p                                                               
  290 continue                                                                  
c                                                                               
      go to 1001                                                                
c     :::::::::: set error -- no convergence to an                              
c                eigenvalue after 30 iterations ::::::::::                      
 1000 ierr = l                                                                  
 1001 return                                                                    
c     :::::::::: last card of tqlrat ::::::::::                                 
      end                                                                       
      subroutine tred2(nm,n,a,d,e,z)                                            
c                                                                               
      integer i,j,k,l,n,ii,nm,jp1                                               
      double precision a(nm,n),d(n),e(n),z(nm,n)                                          
      double precision f,g,h,hh,scale                                                     
*      real*8 dsqrt,dabs,dsign                                                   
c                                                                               
c     this subroutine is a translation of the algol procedure tred2,            
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.           
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).           
c                                                                               
c     this subroutine reduces a real symmetric matrix to a                      
c     symmetric tridiagonal matrix using and accumulating                       
c     orthogonal similarity transformations.                                    
c                                                                               
c     on input:                                                                 
c                                                                               
c        nm must be set to the row dimension of two-dimensional                 
c          array parameters as declared in the calling program                  
c          dimension statement;                                                 
c                                                                               
c        n is the order of the matrix;                                          
c                                                                               
c        a contains the real symmetric input matrix.  only the                  
c          lower triangle of the matrix need be supplied.                       
c                                                                               
c     on output:                                                                
c                                                                               
c        d contains the diagonal elements of the tridiagonal matrix;            
c                                                                               
c        e contains the subdiagonal elements of the tridiagonal                 
c          matrix in its last n-1 positions.  e(1) is set to zero;              
c                                                                               
c        z contains the orthogonal transformation matrix                        
c          produced in the reduction;                                           
c                                                                               
c        a and z may coincide.  if distinct, a is unaltered.                    
c                                                                               
c     questions and comments should be directed to b. s. garbow,                
c     applied mathematics division, argonne national laboratory                 
c                                                                               
c     ------------------------------------------------------------------        
c                                                                               
      do 100 i = 1, n                                                           
c                                                                               
         do 100 j = 1, i                                                        
            z(i,j) = a(i,j)                                                     
  100 continue                                                                  
c                                                                               
      if (n .eq. 1) go to 320                                                   
c     :::::::::: for i=n step -1 until 2 do -- ::::::::::                       
      do 300 ii = 2, n                                                          
         i = n + 2 - ii                                                         
         l = i - 1                                                              
         h = 0.0d0                                                              
         scale = 0.0d0                                                          
         if (l .lt. 2) go to 130                                                
c     :::::::::: scale row (algol tol then not needed) ::::::::::               
         do 120 k = 1, l                                                        
  120    scale = scale + abs(z(i,k))                                           
c                                                                               
         if (scale .ne. 0.0d0) go to 140                                        
  130    e(i) = z(i,l)                                                          
         go to 290                                                              
c                                                                               
  140    do 150 k = 1, l                                                        
            z(i,k) = z(i,k) / scale                                             
            h = h + z(i,k) * z(i,k)                                             
  150    continue                                                               
c                                                                               
         f = z(i,l)                                                             
         g = -sign(sqrt(h),f)                                                 
         e(i) = scale * g                                                       
         h = h - f * g                                                          
         z(i,l) = f - g                                                         
         f = 0.0d0                                                              
c                                                                               
         do 240 j = 1, l                                                        
            z(j,i) = z(i,j) / h                                                 
            g = 0.0d0                                                           
c     :::::::::: form element of a*u ::::::::::                                 
            do 180 k = 1, j                                                     
  180       g = g + z(j,k) * z(i,k)                                             
c                                                                               
            jp1 = j + 1                                                         
            if (l .lt. jp1) go to 220                                           
c                                                                               
            do 200 k = jp1, l                                                   
  200       g = g + z(k,j) * z(i,k)                                             
c     :::::::::: form element of p ::::::::::                                   
  220       e(j) = g / h                                                        
            f = f + e(j) * z(i,j)                                               
  240    continue                                                               
c                                                                               
         hh = f / (h + h)                                                       
c     :::::::::: form reduced a ::::::::::                                      
         do 260 j = 1, l                                                        
            f = z(i,j)                                                          
            g = e(j) - hh * f                                                   
            e(j) = g                                                            
c                                                                               
            do 260 k = 1, j                                                     
               z(j,k) = z(j,k) - f * e(k) - g * z(i,k)                          
  260    continue                                                               
c                                                                               
  290    d(i) = h                                                               
  300 continue                                                                  
c                                                                               
  320 d(1) = 0.0d0                                                              
      e(1) = 0.0d0                                                              
c     :::::::::: accumulation of transformation matrices ::::::::::             
      do 500 i = 1, n                                                           
         l = i - 1                                                              
         if (d(i) .eq. 0.0d0) go to 380                                         
c                                                                               
         do 360 j = 1, l                                                        
            g = 0.0d0                                                           
c                                                                               
            do 340 k = 1, l                                                     
  340       g = g + z(i,k) * z(k,j)                                             
c                                                                               
            do 360 k = 1, l                                                     
               z(k,j) = z(k,j) - g * z(k,i)                                     
  360    continue                                                               
c                                                                               
  380    d(i) = z(i,i)                                                          
         z(i,i) = 1.0d0                                                         
         if (l .lt. 1) go to 500                                                
c                                                                               
         do 400 j = 1, l                                                        
            z(i,j) = 0.0d0                                                      
            z(j,i) = 0.0d0                                                      
  400    continue                                                               
c                                                                               
  500 continue                                                                  
c                                                                               
      return                                                                    
c     :::::::::: last card of tred2 ::::::::::                                  
      end                                                                       
      subroutine tql2(nm,n,d,e,z,ierr)                                          
c                                                                               
      integer i,j,k,l,m,n,ii,l1,nm,mml,ierr                                     
      double precision d(n),e(n),z(nm,n)                                                  
      double precision b,c,f,g,h,p,r,s,machep                                             
*      real*8 dsqrt,dabs,dsign                                                   
c                                                                               
c     this subroutine is a translation of the algol procedure tql2,             
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and             
c     wilkinson.                                                                
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).           
c                                                                               
c     this subroutine finds the eigenvalues and eigenvectors                    
c     of a symmetric tridiagonal matrix by the ql method.                       
c     the eigenvectors of a full symmetric matrix can also                      
c     be found if  tred2  has been used to reduce this                          
c     full matrix to tridiagonal form.                                          
c                                                                               
c     on input:                                                                 
c                                                                               
c        nm must be set to the row dimension of two-dimensional                 
c          array parameters as declared in the calling program                  
c          dimension statement;                                                 
c                                                                               
c        n is the order of the matrix;                                          
c                                                                               
c        d contains the diagonal elements of the input matrix;                  
c                                                                               
c        e contains the subdiagonal elements of the input matrix                
c          in its last n-1 positions.  e(1) is arbitrary;                       
c                                                                               
c        z contains the transformation matrix produced in the                   
c          reduction by  tred2, if performed.  if the eigenvectors              
c          of the tridiagonal matrix are desired, z must contain                
c          the identity matrix.                                                 
c                                                                               
c      on output:                                                               
c                                                                               
c        d contains the eigenvalues in ascending order.  if an                  
c          error exit is made, the eigenvalues are correct but                  
c          unordered for indices 1,2,...,ierr-1;                                
c                                                                               
c        e has been destroyed;                                                  
c                                                                               
c        z contains orthonormal eigenvectors of the symmetric                   
c          tridiagonal (or full) matrix.  if an error exit is made,             
c          z contains the eigenvectors associated with the stored               
c          eigenvalues;                                                         
c                                                                               
c        ierr is set to                                                         
c          zero       for normal return,                                        
c          j          if the j-th eigenvalue has not been                       
c                     determined after 30 iterations.                           
c                                                                               
c     questions and comments should be directed to b. s. garbow,                
c     applied mathematics division, argonne national laboratory                 
c                                                                               
c     ------------------------------------------------------------------        
c                                                                               
c     :::::::::: machep is a machine dependent parameter specifying             
c                the relative precision of floating point arithmetic.           
c                machep = 16.0d0**(-13) for long form arithmetic                
c                on s360 ::::::::::                                             
*      data machep/z'3410000000000000'/
      data machep/1.d-16/
c                                                                               
      ierr = 0                                                                  
      if (n .eq. 1) go to 1001                                                  
c                                                                               
      do 100 i = 2, n                                                           
  100 e(i-1) = e(i)                                                             
c                                                                               
      f = 0.0d0                                                                 
      b = 0.0d0                                                                 
      e(n) = 0.0d0                                                              
c                                                                               
      do 240 l = 1, n                                                           
         j = 0                                                                  
         h = machep * (dabs(d(l)) + dabs(e(l)))                                 
         if (b .lt. h) b = h                                                    
c     :::::::::: look for small sub-diagonal element ::::::::::                 
         do 110 m = l, n                                                        
            if (dabs(e(m)) .le. b) go to 120                                    
c     :::::::::: e(n) is always zero, so there is no exit                       
c                through the bottom of the loop ::::::::::                      
  110    continue                                                               
c                                                                               
  120    if (m .eq. l) go to 220                                                
  130    if (j .eq. 30) go to 1000                                              
         j = j + 1                                                              
c     :::::::::: form shift ::::::::::                                          
         l1 = l + 1                                                             
         g = d(l)                                                               
         p = (d(l1) - g) / (2.0d0 * e(l))                                       
         r = dsqrt(p*p+1.0d0)                                                   
         d(l) = e(l) / (p + dsign(r,p))                                         
         h = g - d(l)                                                           
c                                                                               
         do 140 i = l1, n                                                       
  140    d(i) = d(i) - h                                                        
c                                                                               
         f = f + h                                                              
c     :::::::::: ql transformation ::::::::::                                   
         p = d(m)                                                               
         c = 1.0d0                                                              
         s = 0.0d0                                                              
         mml = m - l                                                            
c     :::::::::: for i=m-1 step -1 until l do -- ::::::::::                     
         do 200 ii = 1, mml                                                     
            i = m - ii                                                          
            g = c * e(i)                                                        
            h = c * p                                                           
            if (dabs(p) .lt. dabs(e(i))) go to 150                              
            c = e(i) / p                                                        
            r = dsqrt(c*c+1.0d0)                                                
            e(i+1) = s * p * r                                                  
            s = c / r                                                           
            c = 1.0d0 / r                                                       
            go to 160                                                           
  150       c = p / e(i)                                                        
            r = dsqrt(c*c+1.0d0)                                                
            e(i+1) = s * e(i) * r                                               
            s = 1.0d0 / r                                                       
            c = c * s                                                           
  160       p = c * d(i) - s * g                                                
            d(i+1) = h + s * (c * g + s * d(i))                                 
c     :::::::::: form vector ::::::::::                                         
            do 180 k = 1, n                                                     
               h = z(k,i+1)                                                     
               z(k,i+1) = s * z(k,i) + c * h                                    
               z(k,i) = c * z(k,i) - s * h                                      
  180       continue                                                            
c                                                                               
  200    continue                                                               
c                                                                               
         e(l) = s * p                                                           
         d(l) = c * p                                                           
         if (dabs(e(l)) .gt. b) go to 130                                       
  220    d(l) = d(l) + f                                                        
  240 continue                                                                  
c     :::::::::: order eigenvalues and eigenvectors ::::::::::                  
      do 300 ii = 2, n                                                          
         i = ii - 1                                                             
         k = i                                                                  
         p = d(i)                                                               
c                                                                               
         do 260 j = ii, n                                                       
            if (d(j) .ge. p) go to 260                                          
            k = j                                                               
            p = d(j)                                                            
  260    continue                                                               
c                                                                               
         if (k .eq. i) go to 300                                                
         d(k) = d(i)                                                            
         d(i) = p                                                               
c                                                                               
         do 280 j = 1, n                                                        
            p = z(j,i)                                                          
            z(j,i) = z(j,k)                                                     
            z(j,k) = p                                                          
  280    continue                                                               
c                                                                               
  300 continue                                                                  
c                                                                               
      go to 1001                                                                
c     :::::::::: set error -- no convergence to an                              
c                eigenvalue after 30 iterations ::::::::::                      
 1000 ierr = l                                                                  
 1001 return                                                                    
c     :::::::::: last card of tql2 ::::::::::                                   
      end                                                                       
