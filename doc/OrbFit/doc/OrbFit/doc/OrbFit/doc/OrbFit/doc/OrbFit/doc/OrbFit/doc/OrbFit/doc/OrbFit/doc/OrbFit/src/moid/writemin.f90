!     ************************************************************      
!     **** C O M P U T I N G     M I N I M U M    P O I N T S ****      
!     ***  and   W R I T I N G   on   O U T P U T    F I L E S ***      
!     ************************************************************      
!     ********** written by GIOVANNI F. GRONCHI (2001) ***********      
!     ******* Department of Mathematics, UNIVERSITY of PISA ******      
!     ******************** modified July 2002 ********************      
!     ============================================================      
      SUBROUTINE writemin(name,lnam,iunit,ekpl,elkep) 
      USE fund_const
      USE planet_masses
      IMPLICIT NONE 
!     asteroid name and name lenght                                     
      CHARACTER*9  name 
      INTEGER lnam 
!     file identifier                                                   
      INTEGER iunit 
!     elements of the Earth and of the asteroid                         
      DOUBLE PRECISION ekpl(6),elkep(6) 
!     in cartesian coordinates                                          
      DOUBLE PRECISION car(6),carpl(6) 
!     cartesian coordinates at minima                                   
      DOUBLE PRECISION cmin(6,20),cplmin(6,20) 
!     check variables                                                   
      DOUBLE PRECISION contrmin(6),contrplmin(6) 
      DOUBLE PRECISION contrelkep(6),contrekpl(6) 
!                                                                       
      DOUBLE PRECISION enne 
!     ------------------ end interface -------------------              
!     eccentric anomalies                                               
      DOUBLE PRECISION umin(20),uplmin(20) 
!     SQUARED DISTANCE function                                         
      DOUBLE PRECISION D2(20) 
!     number of minimum points                                          
      INTEGER nummin 
!                                                                       
      CHARACTER*7 singtype 
!     loop index                                                        
      INTEGER j 
      INTEGER iplam 
!     data for masses                                                   
      INCLUDE 'jplhdr.h90' 
!     ==================================================================
                                                                        
!     duplicate computation of gmse, in case masjpl has not been called 
      gmse=gms*(1.d0+cval(11)/cval(18)) 
                                                                        
!     switching to cartesian coordinates                                
      CALL coocha(elkep,'KEP',gms,car,'CAR',enne) 
      CALL coocha(ekpl,'KEP',gmse,carpl,'CAR',enne) 
                                                                        
!     select Earth as planet                                            
      iplam = 3 
                                                                        
!     =======================================================           
      CALL compute_minima(car,carpl,iplam,cmin,cplmin,D2,nummin) 
!     =======================================================           
      write(*,*)'number of minimum points for ',name(1:lnam), '=',nummin 
                                                                        
!     ******************************************************************
!     ******************** WRITING OUTPUT ON SCREEN ********************
!     ******************************************************************
      WRITE(*,*)'##################################################' 
      WRITE(*,*)'######### M I N I M U M    P O I N T S ###########' 
      WRITE(*,*)'##################################################' 
      WRITE(*,*)'=======================================================&
     &======================='                                          
!     loop on number of stationary points                               
      DO j = 1,nummin 
         write(*,*) 
!     writing on screen                                                 
         WRITE(*,*)'minimum point label =',j,'  local MOID =',dsqrt(D2(j&
     &        ))                                                        
         WRITE(*,*)' ' 
         WRITE(*,*                                                      &
     &        )'Asteroid cartesian coordinates at this minimum point'   
         WRITE(*,*)'      x             y              z         derx   &
     &      dery        derz'                                           
         WRITE(*,107)cmin(1,j),cmin(2,j),cmin(3,j),cmin(4,j),cmin(5,j)  &
     &        ,cmin(6,j)                                                
         WRITE(*,*)'Earth cartesian coordinates at this minimum point ' 
         WRITE(*,*)'      xpl           ypl            zpl       derxpl &
     &      derypl      derzpl'                                         
         WRITE(*,107)cplmin(1,j),cplmin(2,j),cplmin(3,j),cplmin(4,j)    &
     &        ,cplmin(5,j),cplmin(6,j)                                  
         WRITE(*,*)'====================================================&
     &=========================='                                       
!     writing on file iunit.minpts                                      
         WRITE(iunit,106)j,dsqrt(D2(j)) 
         WRITE(iunit,107)cmin(1,j),cmin(2,j),cmin(3,j),cmin(4,j),cmin(5 &
     &        ,j),cmin(6,j)                                             
         WRITE(iunit,107)cplmin(1,j),cplmin(2,j),cplmin(3,j),cplmin(4,j)&
     &        ,cplmin(5,j),cplmin(6,j)                                  
  106    FORMAT(2x,i2,15x,f13.8) 
  107    FORMAT(2x,f11.6,2x,f11.6,2x,f11.6,2x,f11.6,2x,f11.6,2x,f11.6) 
                                                                        
!     check on the distances                                            
!         write(*,*)'-----------------------------------------'         
!         WRITE(*,*)'controllo di D2(',j,')=',dsqrt((cmin(1,j)-cplmin(1,
!     $        ))**2+(cmin(2,j)-cplmin(2,j))**2+(cmin(3,j)-cplmin(3,j)) 
!     $        **2)                                                     
!         write(*,*)'-----------------------------------------'         
                                                                        
!     END MAIN DO LOOP                                                  
      ENDDO 
                                                                        
!     ========================== CHECK ===========================      
!      DO j = 1,6                                                       
!         contrmin(j) = cmin(j,1)                                       
!         contrplmin(j) = cplmin(j,1)                                   
!      ENDDO                                                            
!     switching to keplerian coordinates                                
!      CALL coocha(contrmin,'CAR',gms,contrelkep,'KEP',enne)            
!      CALL coocha(contrplmin,'CAR',gmse,contrekpl,'KEP',enne)          
!      write(*,*)'=========================================='           
!      write(*,*)'keplerian coordinates of asteroid at MOID:'           
!      write(*,*)contrelkep                                             
!      write(*,*)'keplerian coordinates of Earth at MOID:   '           
!      write(*,*)contrekpl                                              
!      write(*,*)'=========================================='           
!     ============================================================      
                                                                        
      WRITE(*,*) 
      WRITE(*,*)'===================================================' 
      WRITE(*,*)'###################################################' 
                                                                        
      RETURN 
      END                                           
