      double precision function dcabs1(z) 
       complex*16 z 
!      complex*16 z,zz                                                  
!      double precision t(2)                                            
!      equivalence (zz,t(1))                                            
!      zz = z                                                           
!      dcabs1 = dabs(t(1)) + dabs(t(2))                                 
!       dcabs1=abs(real(z))+abs(real(z*(0,1)))                          
       dcabs1=abs(DBLE(z))+abs(DBLE(z*(0,1))) 
                                                                        
       return 
      END                                           
