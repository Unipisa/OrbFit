! Copyright (C) 1998-2000 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version 3.3.2, A. Milani 2006                                               
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         G A U S S D E G 8                     *    
!  *                                                               *    
!  *       Initial orbit determination with Gauss' method          *
!  *       Limited to deg 8 dynamical equation                     *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    TOBS      -  Time of observations (MJD, TDT)                
!           ALPHA     -  Right ascension (rad)                          
!           DELTA     -  Declination (rad)                              
!           OBSCOD    -  Observatory code                               
!           DEBUG     -  Print debug information                        
!                                                                       
! OUTPUT:   EL        -  Orbital elements (mean ecliptic                
!                                          and equinox J2000)           
!           NROOTS    -  Number of positive roots of 8th degree pol.    
!           NSOL      -  Number of solutions (some roots may be         
!                        discarded if they lead to solution with        
!                        eccentricity > ECCMAX)                         
!           FAIL      -  Error flag                                     
!           MSG       -  Error message                                  
!                                                                       
SUBROUTINE gaussdeg8(tobs,alpha,delta,obscod,el,nroots,nsol,fail,msg,debug)
  USE reference_systems 
  USE fund_const
  USE orbit_elements
  IMPLICIT NONE  
  INTEGER, INTENT(IN) ::  obscod(3) 
  DOUBLE PRECISION, INTENT(IN):: tobs(3),alpha(3),delta(3) 
  TYPE(orbit_elem), INTENT(OUT) :: el(3)
  LOGICAL, INTENT(OUT) :: fail
  LOGICAL, INTENT(IN) :: debug
  CHARACTER*(*), INTENT(OUT) :: msg
  INTEGER, INTENT(OUT) :: nsol,nroots
! end interface
  DOUBLE PRECISION elem(6,8),t0(8) 
  INTEGER i,j,k,ising,ir,it 
  DOUBLE PRECISION xt(3,3),sinv0(3,3),a(3),b(3),c(3) 
  DOUBLE PRECISION ra(3),rb(3),coef(0:8),a2star,b2star,r22,s2r2 
  DOUBLE PRECISION esse0(3,3),cosd,det,tau1,tau3,tau13 
  DOUBLE PRECISION roots(8),esse(3,3),esse1(3,3),sinv(3,3),r2m3 
  DOUBLE PRECISION gcap(3),crhom(3),rho(3),xp(3,3),vp(3),xv(6) 
  DOUBLE PRECISION vekp(6),v1(3),v2(3),vs,err,sca,vaber(3),xv1(6) 
  DOUBLE PRECISION xp1(3,3),tis2   !,rot(3,3) 
  DOUBLE PRECISION xve(6),vekpe(6),fs1,gs1,fs3,gs3,fggf 
      
                                                                        
  DOUBLE PRECISION vsize,princ 
  EXTERNAL vsize,princ 
                   
  INTEGER fail_flag
                                                     
  fail=.true. 
  msg=' ' 
  nroots=0 
  nsol=0 
                                                                        
! COMPUTATION OF PRELIMINARY ORBIT                                      
!                                                                       
! ESSE = unit vector pointing in the direction of observations          
  DO 10 k=1,3 
     cosd=COS(delta(k)) 
     esse0(1,k)=cosd*COS(alpha(k)) 
     esse0(2,k)=cosd*SIN(alpha(k)) 
     esse0(3,k)=SIN(delta(k)) 
! Position of the observer                                              
     DO 11 i=1,3 
        CALL posobs(tobs(i),obscod(i),1,xt(1,i)) 
11   END DO
10 END DO
                                                                        
! Inverse of ESSE matrix                                                
  sinv0=esse0 
  CALL matin(sinv0,det,3,0,3,ising,1) 
  IF(ising.NE.0) THEN 
     msg='Singular S matrix (coplanar orbits?)' 
     RETURN 
  END IF
  tis2=tobs(2) 
                                                                        
! A and B vectors                                                       
  tau1=gk*(tobs(1)-tis2) 
  tau3=gk*(tobs(3)-tis2) 
  tau13=tau3-tau1 
  a(1)= tau3/tau13 
  a(2)=-1.d0 
  a(3)=-(tau1/tau13) 
  b(1)=a(1)*(tau13**2-tau3**2)/6.d0 
  b(2)=0.d0 
  b(3)=a(3)*(tau13**2-tau1**2)/6.d0 
  
! Coefficients of 8th degree equation for r2                            
  CALL prodmv(ra,xt,a) 
  CALL prodmv(rb,xt,b) 
  a2star=sinv0(2,1)*ra(1)+sinv0(2,2)*ra(2)+sinv0(2,3)*ra(3) 
  b2star=sinv0(2,1)*rb(1)+sinv0(2,2)*rb(2)+sinv0(2,3)*rb(3) 
  r22=xt(1,2)**2+xt(2,2)**2+xt(3,2)**2 
  s2r2=esse0(1,2)*xt(1,2)+esse0(2,2)*xt(2,2)+esse0(3,2)*xt(3,2) 
  coef(8)=1.d0 
  coef(7)=0.d0 
  coef(6)=-(a2star**2)-r22-(2.d0*a2star*s2r2) 
  coef(5)=0.d0 
  coef(4)=0.d0 
  coef(3)=-(2.d0*b2star*(a2star+s2r2)) 
  coef(2)=0.d0 
  coef(1)=0.d0 
  coef(0)=-(b2star**2) 
  IF(coef(0).EQ.0) THEN 
     msg='coef(0)=0' 
     RETURN 
  END IF
  CALL solv8(coef,roots,nroots) 
  IF(debug) THEN 
!     WRITE(*,520) 
!     WRITE(*,521) (i,coef(i),i=0,8) 
     WRITE(*,513) 2 
     IF(nroots.LE.0) THEN 
        WRITE(*,524) 
     ELSE 
        WRITE(*,514)(i,roots(i),i=1,nroots) 
     END IF
  END IF
!520 FORMAT(12X,'Coefficients of 8-th degree polynomial:') 
!521 FORMAT(16X,'coef(',i1,')  =',F18.10) 
513 FORMAT(12X,'Possible roots for Sun-asteroid distance at ',        &
         &       'observation',I2,':')                                      
514 FORMAT(16X,'r(',i1,')  =',F10.6) 
524 FORMAT(16X,'none') 
                                                                        
  IF(nroots.LE.0) THEN 
     msg='8th degree polynomial has no real roots' 
     RETURN 
  ELSEIF(nroots.gt.3)THEN
     WRITE(*,*)' gaussdeg8: more than 3 positive roots ',nroots
     STOP
  END IF
                                                                        
! Orbital elements of preliminary solution
  nsol=0                              
  DO 20 ir=1,nroots 
     esse=esse0 
     esse1=esse0 
     sinv=sinv0 
     r2m3=1.d0/(roots(ir)**3) 
     c(1)=a(1)+b(1)*r2m3 
     c(2)=-1.d0 
     c(3)=a(3)+b(3)*r2m3 
     CALL prodmv(gcap,xt,c) 
     CALL prodmv(crhom,sinv,gcap) 
     DO 13 k=1,3 
        rho(k)=-(crhom(k)/c(k)) 
13   END DO
                                                                        
! Position of the asteroid at the time of observations                  
     DO 14 k=1,3 
        xp(1:3,k)=xt(1:3,k)+rho(k)*esse(1:3,k) 
14   ENDDO

! remove spurious root
     IF(rho(2).lt.0.01d0)THEN
        IF(debug)WRITE(*,*) ' spurious root , r, rho ',roots(ir),rho(2)
        CYCLE
     ELSE
        IF(debug)WRITE(*,*) ' accepted root , r, rho ',roots(ir),rho(2)
        nsol=nsol+1
     ENDIF
     
! Gibbs' transformation, giving the velocity of the planet at the       
! time of second observation                                            
     CALL gibbs(xp,tau1,tau3,vp,gk) 
! Orbital elements of preliminary orbit                                 
     DO 15 i=1,3 
        xv(i)=xp(i,2) 
        xv(i+3)=vp(i) 
15   END DO
     el(nsol)=undefined_orbit_elem
     el(nsol)%coord=xv
     el(nsol)%coo='CAR'
! planetary aberration 
     el(nsol)%t=tobs(2)-rho(2)/vlight
     CALL coo_cha(el(nsol),'KEP',el(nsol),fail_flag)
     IF(debug) THEN 
        WRITE(*,525) ir, fail_flag 
        IF(el(nsol)%coo.EQ.'KEP') THEN 
           WRITE(*,535) 'a',el(nsol)%coord(1),el(nsol)%coord(2) 
        ELSEIF(el(nsol)%coo.EQ.'COM') THEN 
           WRITE(*,535) 'q',el(nsol)%coord(1),el(nsol)%coord(2)
        ELSE 
           STOP '**** gaussdeg8: unknown elem type ****' 
        END IF
     END IF
525  FORMAT(12X,'ROOT NO.',I2,1x,I2) 
535  FORMAT(16X,'Preliminary orbit: ',A,' =',F10.5,';  ecc =',F10.5) 
                                                                        
                                                                        
20 END DO
                                                                        
  IF(nsol.LE.0) THEN 
     msg='No acceptable solution' 
     RETURN 
  END IF
                                                                        
  fail=.false. 
                                                                        
END SUBROUTINE gaussdeg8
