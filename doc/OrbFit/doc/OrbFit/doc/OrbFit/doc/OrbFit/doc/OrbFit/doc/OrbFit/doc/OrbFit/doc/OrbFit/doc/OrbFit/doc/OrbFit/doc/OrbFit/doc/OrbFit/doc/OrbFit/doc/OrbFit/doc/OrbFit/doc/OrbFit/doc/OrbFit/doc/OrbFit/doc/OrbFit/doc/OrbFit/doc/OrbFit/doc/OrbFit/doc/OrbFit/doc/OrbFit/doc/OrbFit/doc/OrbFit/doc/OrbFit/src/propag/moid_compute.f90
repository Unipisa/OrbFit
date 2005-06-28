! =========MODULE moid_compute==============                            
! CONTAINS                                                              
! ROUTINES                                                              
!             nomoid                                                    
!             nodedi                                                    
!                                                                       
!  HEADERS  MODULES                                                            
! moid_compute.o:                                                       
!            parbep.h masses.h masses Earth and Moon                    
!            fund_const.mod for mass of the Sun

! =========================================                             
! NOMOID - compute moid using Gronchi's routine                         
! nodal distances and Minimum Orbital Intersection Distance             
! with respect to the Earth                                             
! for an orbit, given in equinoctal elements, at a given time           
! =========================================                             
SUBROUTINE nomoid(t0,el0,moid,dnp,dnm) 
  USE fund_const
  USE planet_masses
  USE orbit_elements
  IMPLICIT NONE 
! ===========INPUT=====================                                 
! epoch (MJD)                       
  DOUBLE PRECISION, INTENT(IN) :: t0 
  TYPE(orbit_elem), INTENT(IN) :: el0
! ===========OUTPUT=================================                    
! MOID, iteration count and succ. flag                                  
! ascending node distance,  descending node distance                    
  DOUBLE PRECISION moid,dnp,dnm 
! ==========END INTERFACE================================               
! elements of Earth (equinoctal*, cartesian coordinates   
  TYPE(orbit_elem) :: eleq, elcar
  DOUBLE PRECISION eqp(6),xast(6),xea(6) 
! cartesian coordinates asteroid and planet at minimum                  
  DOUBLE PRECISION cmin(3,16),cplmin(3,16) 
!     SQUARED DISTANCE function                                         
  DOUBLE PRECISION d2(16) 
!     number of relative minima found                                   
  INTEGER nummin 
! error in coordinate change
  INTEGER fail_flag, fail_flag1
! ======================================================                
! get Earth elements                                                    
  CALL earth(t0,eqp)
! get equinoctal elements
  CALL coo_cha(el0,'EQU',eleq, fail_flag)
  IF(fail_flag.ge.4)THEN
     WRITE(*,*)' nomoid: failed change to EQU ',el0
     moid=-1.d0
     dnp=-1.d0
     dnm=-1.d0
     RETURN
  ENDIF 
! compute nodal distances                                               
  CALL nodedi(eleq%coord,eqp,dnp,dnm) 
! transform to cartesian coordinates                                    
  CALL coo_cha(el0,'CAR',elcar, fail_flag1)
  xast=elcar%coord
  CALL earcar(t0,xea,1)
! compute moid                                                          
  CALL compute_minima_ta(xast,xea,3,cmin,cplmin,d2,nummin) 
  moid=sqrt(d2(1)) 
END SUBROUTINE nomoid
! ==========================================                            
! NODEDI                                                                
! nodal distances of two elliptic orbits                                
! ==========================================                            
SUBROUTINE nodedi(eq,eqp,dnp,dnm) 
  USE fund_const
  USE planet_masses
  IMPLICIT NONE 
! ===========INPUT=====================                                 
! elements of asteroid, of Earth (equinoctal)                           
  DOUBLE PRECISION eq(6),eqp(6) 
! ===========OUTPUT=================================                    
! output ascending node distance,  descending node distance             
  DOUBLE PRECISION dnp,dnm 
! ==========END INTERFACE================================               
  DOUBLE PRECISION c(3),cp(3),x(6),xp(6),vlenz(3),vlenzp(3) 
  DOUBLE PRECISION enne,ennep,vnod(3),vnl,ome,omep 
  DOUBLE PRECISION ecc,eccp,cosf,cosfp,rp,rm,rpe,rme 
  INTEGER i 
  DOUBLE PRECISION prscal,vsize 
! cartesian coordinates                                                 
  CALL coocha(eq,'EQU',gms,x,'CAR',enne) 
  CALL coocha(eqp,'EQU',gmse,xp,'CAR',ennep) 
!  angular momentum                                                     
  CALL prvec(x,x(4),c) 
  CALL prvec(xp,xp(4),cp) 
! ascending node                                                        
  CALL prvec(cp,c,vnod) 
  vnl=vsize(vnod) 
!  angular momentum unit vector, Lenz vector                            
  CALL  prvec(x(4),c,vlenz) 
  DO i=1,3 
     vlenz(i)=vlenz(i)/gms-x(i)/vsize(x) 
  ENDDO
  ecc=vsize(vlenz) 
  CALL  prvec(xp(4),cp,vlenzp) 
  DO i=1,3 
     vlenzp(i)=vlenzp(i)/gmse-xp(i)/vsize(xp) 
  ENDDO
  eccp=vsize(vlenzp) 
! true anomaly at mutual node= - arg. of perihelion                     
  cosf=prscal(vnod,vlenz)/(vnl*ecc) 
  ome=acos(cosf) 
  cosfp=prscal(vnod,vlenzp)/(vnl*eccp) 
  omep=acos(cosfp) 
! nodal points and distances                                            
  rp=eq(1)*(1.d0-eq(2)**2-eq(3)**2)/(1.d0+ecc*cosf) 
  rm=eq(1)*(1.d0-eq(2)**2-eq(3)**2)/(1.d0-ecc*cosf) 
  rpe=eqp(1)*(1.d0-eqp(2)**2-eqp(3)**2)/(1.d0+eccp*cosfp) 
  rme=eqp(1)*(1.d0-eqp(2)**2-eqp(3)**2)/(1.d0-eccp*cosfp) 
  dnp=rp-rpe 
  dnm=rm-rme 
!      WRITE(*,*)rp,rpe,dnp,rm,rme,dnm,cosf,cosfp,ome,omep              
END SUBROUTINE nodedi
