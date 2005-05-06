
! ************************************************************      
! ***** C O M P U T I N G    M I N I M U M   P O I N T S *****      
! ************************************************************      
! ********** written by GIOVANNI F. GRONCHI (2001) ***********      
! ******* Department of Mathematics, UNIVERSITY of PISA ******      
! last modified June 2004 (GFG)
! ==================================================================      
  SUBROUTINE compute_minima(car,carpl,iplam,cmin,cplmin,D2,nummin)
    USE fund_const 
    USE planet_masses
    IMPLICIT NONE 
! cartesian coordinates of the asteroid and of the planet:
    DOUBLE PRECISION,INTENT(IN) :: car(6),carpl(6) 
    INTEGER,INTENT(IN) :: iplam ! planet number (3 = Earth) 
! cartesian coordinates of the asteroid and of the planet 
! at minimum points:
    DOUBLE PRECISION,INTENT(OUT) :: cmin(6,20),cplmin(6,20) 
    DOUBLE PRECISION,INTENT(OUT) :: D2(20) ! squared distance function 
                                             ! at minimum points 
    INTEGER,INTENT(OUT) :: nummin ! number of relative minima found 
! ========= end interface ==========================================
    INTEGER nstat ! number of stationary points (maximum = 20)   
    INTEGER nummax ! number of relative minima/maxima found
    INTEGER :: answer(20) ! type of singular point : 
!                             answer(j) =  1    MAXIMUM
!                             answer(j) = -1    MINIMUM                
!                             answer(j) =  0    SADDLE                 
!                             answer(j) =  2    CANNOT DECIDE          
    CHARACTER*7 :: singtype 

    LOGICAL :: morse ! check with Morse theory 
!                        morse = true      OK   
!                        morse = false     ERROR  
    LOGICAL :: weier ! check with Weierstrass theory:
!                        weier = true      OK 
!                        weier = false     ERROR
    LOGICAL :: warnflag(3) ! program warning flags:
!      warnflag(1) = true    OK
!      warnflag(1) = false   leading coefficient of resultant is very small
!      warnflag(2) = true    OK                                          
!      warnflag(2) = false   higher degree terms in resultant are not small
!      warnflag(3) = true    OK                                           
!      warnflag(3) = false   low precision for resultant coefficients     
    LOGICAL :: sflag(6) ! solving system messages:
!         sflag(1) = true    OK
!         sflag(1) = false   there are two good solutions             
!         sflag(2) = true    OK                                               
!         sflag(2) = false   neither of the two evaluation is close to 0
!         sflag(3) = true    OK                                               
!         sflag(3) = false   the s-component of the solution is complex
!         sflag(4) = true    OK    
!         sflag(4) = false   leading coeff. of the 2nd degree pol. is small
!         sflag(5) = true    OK    
!         sflag(5) = false   1st and 2nd coeffs except of 2nd degree poly
!                            are small                       
!         sflag(6) = true    OK    
!         sflag(6) = false   2nd degree poly have all coefficients small
    LOGICAL :: hzflag ! hzflag = .true.  OK!
                      ! hzflag = .false. abs(root)>10^5
                      ! (large root of the resultant poly)
    LOGICAL :: hwflag ! hwflag = .true.  OK!
                      ! hwflag = .false. abs(root)>10^5
                      ! (large root of 2nd deg poly)
    LOGICAL :: multfl ! multfl = .true.  OK!
                      ! multfl = .false. 0 has multiplicity > 4

    DOUBLE PRECISION :: ekpl(6),elkep(6) ! keplerian elements
    DOUBLE PRECISION :: u(20),upl(20),umin(20),uplmin(20)!eccentric anomalies
    DOUBLE PRECISION lmin(20),lplmin(20)! mean anomalies
! auxiliary                                               
    DOUBLE PRECISION :: cmintmp(6),cplmintmp(6) 
    DOUBLE PRECISION :: uu,uupl,chu,chupl 
    DOUBLE PRECISION :: DD2 
! ------------------------------------------
! mutual orbital elements 
    DOUBLE PRECISION :: mutI,mutom,mutompl
    DOUBLE PRECISION :: a,e,i,om,apl,epl,ompl 
    DOUBLE PRECISION :: beta,betapl 
    COMMON/bbpl/beta,betapl 
! ------------------------------------------
    DOUBLE PRECISION :: enne ! for coocha
    DOUBLE PRECISION gmsp ! mass of the planet + mass of the sun  
! -------------------------------------------------
! ordered 3-uples, for sorting 
    DOUBLE PRECISION :: ordu(20),ordupl(20),ordD2(20) 
    INTEGER :: ordans(20)                            
! -------------------------------------------------                      
    INTEGER :: k,j ! loop index 
! ==================================================================

!   flags initialization                                              
    weier = .true.
!    weier = .false. ! dummy
    morse = .true.
    sflag(1) = .true. 
    sflag(2) = .true. 
    sflag(3) = .true.
    sflag(4) = .true.
    sflag(5) = .true.
    sflag(6) = .true.
    warnflag(1) = .true.
    warnflag(2) = .true.
    warnflag(3) = .true.
    hzflag = .true.
    hwflag = .true.
    multfl = .true.

    gmsp=gms+gm(iplam) ! mass of planet plus mass of the sun
!   gmse = 0.000295913108d0
    cmin=0.d0
    cplmin=0.d0
!   switch to keplerian coordinates                                
    CALL coocha(car,'CAR',gms,elkep,'KEP',enne) 
    CALL coocha(carpl,'CAR',gmse,ekpl,'KEP',enne) 
    IF(enne.le.0.d0)THEN
       nummin=0
       RETURN              
    ENDIF
! *************************************************
! switch to mutual reference system 
    CALL mutualrefcha(ekpl,elkep,mutI,mutom,mutompl) 
    a    =  elkep(1) 
    e    =  elkep(2) 
    apl  =  ekpl(1) 
    epl  =  ekpl(2) 
    i    =  mutI 
    om   =  mutom 
    ompl =  mutompl 
! *************************************************
                                                                        
    beta=dsqrt(1-e**2) 
    betapl=dsqrt(1-epl**2) 
                                                                        
! compute coefficients of the trigonometric polynomials
    CALL aical(a,e,i,om,apl,epl,ompl) 
                                                                        
! ******************************************************************
! compute stationary points of D2 
    CALL comp_heart_rot(u,upl,nstat,nummin,nummax, &
         & answer,warnflag,sflag,morse,weier,hzflag,hwflag,multfl)
! ******************************************************************
                                                                        
!   check
    if(.not.sflag(3)) then
       write(*,*)'compute_minima: negative discriminant! ecc=',e
    endif

    DO j = 1,nstat 
                                                                        
!   conversion in radians for D2eval                                  
       uu=u(j)*radeg 
       uupl=upl(j)*radeg 

!   compute squared distance at critical points
       CALL D2eval(uu,uupl,DD2) 
       D2(j)=DD2 

!   angles between 0 and 360 degrees                                  
       CALL choosedeg(u(j),chu) 
       CALL choosedeg(upl(j),chupl) 
       u(j) = chu 
       upl(j) = chupl 
                                                                        
    ENDDO
                                                                        
! sorting 3-uples (u,upl,D2) according to the value of D2               
    CALL sortd2(nstat,u,upl,D2,answer) 
                                                                        
    DO k = 1,nummin 
       umin(k) = u(k) 
       uplmin(k) = upl(k) 
       
       lmin(k) = u(k)*radeg - elkep(2)*dsin(u(k)*radeg) 
       lplmin(k) = upl(k)*radeg - ekpl(2)*dsin(upl(k)*radeg) 
       elkep(6) = lmin(k) 
       ekpl(6) = lplmin(k) 

! switch to cartesian coordinates                                
       CALL coocha(ekpl,'KEP',gmse,cplmintmp,'CAR',enne) 
       CALL coocha(elkep,'KEP',gms,cmintmp,'CAR',enne) 
       DO j = 1,6 
          cmin(j,k) = cmintmp(j) 
          cplmin(j,k) = cplmintmp(j) 
       ENDDO
    ENDDO
                                                                        
!    write(*,*)'u,upl,D2',umin(1)*radeg,uplmin(1)*radeg,D2(1)

    RETURN 
  END SUBROUTINE compute_minima
