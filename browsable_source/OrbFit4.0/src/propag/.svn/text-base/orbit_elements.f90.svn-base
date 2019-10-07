! ORBIT_ELEMENTS
! contains:
!
! OUT OF MOD:
!      appmag     correction to absolute magnitude to get apparent


MODULE orbit_elements
USE output_control
USE fund_const
IMPLICIT NONE
PRIVATE

! public routines
PUBLIC coo_cha, equal_orbels, write_elems, read_elems, ecc_peri, convertunc 
PUBLIC cootyp, attelements, eldiff, tpcar_mtpcar, opik_tpcar
PUBLIC centerbody_mass, att_prelim2 ! ??? used by bizarre
PUBLIC find_elems, wrikep

INTEGER, PARAMETER :: ndimx = 6     ! maximum dimension of elements vector

DOUBLE PRECISION :: eccbou ! boundary in e between use of Kepler equation
                           ! and use of f-g series

! Generic orbital elements set
TYPE orbit_elem
CHARACTER(LEN=3)       :: coo     ! elements type; known types are:
!    'CAR' : cartesian positions and velocities 
!    'EQU' : equinoctal elements
!    'KEP' : classical keplerian elements, singular for
!             0 eccentricity, 0 and 180 degrees inclination;
!    'COM' : cometary elements, valid for ecc.ge.1 
!    'ATT' : attributable plus r rdot
!    'OPI' : Opik type elements not yet implemented

INTEGER :: ndim  ! actual dimension of elements vector
DOUBLE PRECISION, DIMENSION(ndimx) :: coord ! elements vector 
              ! DANGER: some of the coordinates may actually be angles...   
DOUBLE PRECISION :: t  ! epoch time, MJD, TDT

LOGICAL :: mag_set ! true if some absolute magnitude is available
DOUBLE PRECISION :: h_mag, g_mag ! Absolute magnitude H, opposition effect G

INTEGER :: center ! 0 for Sun, code for planet as in JPL ephemerides

INTEGER :: obscode ! observatory code, for ATT type only
 
END TYPE orbit_elem

! default value when undefined
DOUBLE PRECISION, DIMENSION(6), PARAMETER :: zero_6d_vect = &
&    (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)

! undefined orbital elements
TYPE(orbit_elem), PARAMETER :: undefined_orbit_elem = ORBIT_ELEM( &
&  'UNK',           & ! no coordynate defined
&   6,              & ! default is 6, for now
& zero_6d_vect,     & ! null elements
&   1.d99,          & ! null epoch 
&   .false.,        & ! default no magnitudes
&   -9.99d0,        & ! null magnitude
&    0.15d0,        & ! default opposition effect
&    0,             & ! default is heliocentric
&    500            & ! default is not topocentric correction 
! &  ,gms            & ! if needed
&     )  ! 

TYPE orb_uncert
! normal and covariance matrix (warning: covariance not full if not 
! all solved)
DOUBLE PRECISION, DIMENSION(ndimx,ndimx):: c,g 
! norm of residuals, of last correction, RMS of magnitudes
LOGICAL :: succ ! logical for success (convergent differential corrections)
INTEGER :: ndim  ! actual dimension of elements vector
END TYPE orb_uncert

! default value when undefined
DOUBLE PRECISION, DIMENSION(6,6), PARAMETER :: zero_6x6_matrix = &
& RESHAPE((/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
& 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
& 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
& 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /),(/ 6, 6 /))

TYPE(orb_uncert), PARAMETER :: undefined_orb_uncert= ORB_UNCERT( &
& zero_6x6_matrix , zero_6x6_matrix, & ! no matrices defined
.false., 6 ) 

TYPE orb_with_uncert
TYPE(orbit_elem) :: orbit ! nominal orbit, solution of least squares
TYPE(orb_uncert) :: uncertainty ! covariance matrix etc.
END TYPE orb_with_uncert

! public entities
PUBLIC orbit_elem,undefined_orbit_elem,orb_uncert,undefined_orb_uncert
PUBLIC ndimx, eccbou

CONTAINS

! ===================================================================   
! COO_CHA                                                         
! ===================================================================   
!   general purpose elements change 
!    
!  fail _flag = 0 conversion OK
!             = 1 zero eccentricity (conv. to KEP/COM fails)
!             = 2 zero inclination               "
!             = 3 zero ecc and zero inclination   "
!             = 5 eccentricity 1 and more (conv to EQU/KEP fails) 
!             = 9 180 deg inclination (conv to KEP/COM/EQU fails)
!             = 10 zero angular momentum (conv to EQU/KEP/COM fails)
!             = 20 inconsistent calling sequence 
!  if fail_flag =1,2,3 the jacobian is undefined for kep_equ, replaced
!                by an arbitrary matrix (to avoid undefined on exit)
!  if fail_flag =5 conversion to KEP, EQU is refused, something else
!                  is provided in output (should be COM)
!  if fail_flag > 8 there is no conversion provided....
!               we should set up the messages in a suitable way,
!               providing for a smooth failure
! =================================================================== 
SUBROUTINE coo_cha(x,cooy,y,fail_flag,del,iplanet,obscode)
TYPE(orbit_elem), INTENT(IN) :: x
CHARACTER*3, INTENT(IN) :: cooy
TYPE(orbit_elem), INTENT(INOUT) :: y ! needed to allow for x=y
INTEGER, INTENT(OUT) :: fail_flag ! see above
! partial derivatives \partial y/\partial x
DOUBLE PRECISION, INTENT(OUT), OPTIONAL, DIMENSION(x%ndim,x%ndim) :: del 
INTEGER, INTENT(IN), OPTIONAL :: obscode ! output observer
INTEGER, INTENT(IN), OPTIONAL :: iplanet ! output center
! END INTERFACE
TYPE(orbit_elem) z
LOGICAL kepgr_in,kepgr_ou,cargr_in,cargr_ou,opposite, done
CHARACTER*3 coox
! target_center is the output gravitational center
INTEGER ndim,obscod,planet,target_center
DOUBLE PRECISION, DIMENSION(x%ndim,x%ndim) :: del1,del2 

done=.false.
fail_flag=0
coox=x%coo
ndim=x%ndim
IF(PRESENT(iplanet))THEN
   planet=iplanet
   target_center=planet
ELSEIF(coox.eq.'ATT')THEN
   IF(rhs.eq.1.or.rhs.eq.3)THEN
      target_center=0
   ELSEIF(rhs.eq.2)THEN
      target_center=3
   ELSE
      WRITE(*,*)' coo_cha: rhs not initialized, rhs=', rhs
      STOP
   ENDIF
ELSE
   target_center=x%center
ENDIF

!IF(coox.ne.'ATT')THEN
!   target_center=x%center
!ELSEIF(cooy.eq.'ATT')THEN
!   planet=3
!   target_center=0
!ELSE
!   planet=0
!   target_center=0
!ENDIF
!IF(x%center.ne.planet)THEN
!   WRITE(*,*)'coo_cha: change of center from ',x%center,' to ',planet
!ENDIF
!NO!!! if change of center???
IF(coox.eq.cooy)THEN
   y=x
   IF(PRESENT(del)) CALL eye(x%ndim,del)
   done=.true.
   RETURN
ENDIF
IF(.not.PRESENT(obscode).and.(cooy.eq.'ATT'))THEN
   obscod=500
ELSEIF(PRESENT(obscode))THEN
   obscod=obscode
ENDIF
kepgr_in=(coox.eq.'KEP'.or.coox.eq.'EQU'.or.coox.eq.'COM'.or.coox.eq.'COT')
kepgr_ou=(cooy.eq.'KEP'.or.cooy.eq.'EQU'.or.cooy.eq.'COM'.or.cooy.eq.'COT')
cargr_in=(coox.eq.'CAR'.or.coox.eq.'ATT')
cargr_ou=(cooy.eq.'CAR'.or.cooy.eq.'ATT')
!opposite=((kepgr_in.and.cargr_ou).or.(cargr_in.and.kepgr_ou))

IF(cargr_in)THEN
   IF(coox.eq.'ATT')THEN
!      IF(planet.eq.0)THEN
! attributable case must be converted to cartesian heliocentric
      IF(PRESENT(del))THEN
         z=car_att(x,target_center,del1) ! convert to CAR
      ELSE
         z=car_att(x,target_center)      ! convert to CAR
      ENDIF
   ELSEIF(coox.eq.'CAR')THEN
      z=x                                 ! already cartesian
      IF(PRESENT(del))THEN
         CALl eye(ndim,del1)     ! initialize by identity
      ENDIF
   !   ELSEIF(coox.eq.'OPI')THEN
   !      WRITE(*,*)' coo_cha: this conversion not ready ', coox, ' ',cooy
   !      STOP
   ENDIF
   IF(cooy.eq.'CAR')THEN
      y=z ! done
      IF(PRESENT(del))del=del1
   ELSEIF(cooy.eq.'ATT')THEN
      IF(PRESENT(del))THEN
         y=att_car(z,obscod,target_center,del2)
         del=MATMUL(del2,del1)
      ELSE
         y=att_car(z,obscod,target_center)
      ENDIF
   ELSEIF(kepgr_ou)THEN
! convert to elements  
      IF(PRESENT(del))THEN
         y=to_elems(z,cooy,fail_flag,del2)
         del=MATMUL(del2,del1)
      ELSE
         y=to_elems(z,cooy,fail_flag)
      ENDIF
! handle fail_flag
      IF(fail_flag.gt.0.and.fail_flag.le.3)THEN
! zero eccentricity/inclination
! EQU are differentiable
         IF(cooy.eq.'EQU')fail_flag=0 
! KEP, COM are defined but not diff.
         IF(.not.PRESENT(del))fail_flag=0 
      ELSEIF(fail_flag.eq.5)THEN 
! force to COT, leaving fail_flag=5 as warning
         IF(PRESENT(del))THEN
 ! should not be necessary
            IF(y%coo.ne.'COT')THEN
               WRITE(ierrou,*)' coo_cha: fail_flag, outcoo ', fail_flag, y%coo
               numerr=numerr+1
               y=to_elems(z,'COT',fail_flag,del2)
            ENDIF
            del=MATMUL(del2,del1)
         ELSE
            IF(y%coo.ne.'COT')THEN ! should not be necessary
               WRITE(ierrou,*)' coo_cha: fail_flag, outcoo ', fail_flag, y%coo
               numerr=numerr+1
               y=to_elems(z,'COT',fail_flag)
            ENDIF
         ENDIF
      ELSE
! fail_flag=9,10,20: nothing to do
      ENDIF
      done=.true.
      RETURN
   ELSE
      PRINT *, ' coo_cha: this should not happen 1, coox=',coox,' cooy=', cooy
      WRITE(*,*) x
      STOP
   ENDIF
ENDIF
IF(kepgr_in.and.cargr_ou)THEN
   IF(PRESENT(del))THEN
      z=car_elems(x,fail_flag,del1)
   ELSE
      z=car_elems(x,fail_flag)
   ENDIF
   IF(cooy.eq.'ATT')THEN
      IF(PRESENT(del))THEN
         y=att_car(z,obscod,target_center,del2)
         del=MATMUL(del2,del1)
      ELSE
         y=att_car(z,obscod,target_center)
      ENDIF
   ELSEIF(cooy.eq.'CAR')THEN
      y=z ! done
      IF(PRESENT(del))del=del1
   ENDIF
! handle fail_flag
   IF(fail_flag.gt.0.and.fail_flag.le.3)THEN
      IF(coox.eq.'EQU')fail_flag=0 ! EQU are differentiable
      IF(.not.PRESENT(del))fail_flag=0 ! KEP, COM are defined but not diff.
   ELSE
! fail_flag=5 cannot happen, 9,10,20 cannot be fixed
   ENDIF
   done=.true.
   RETURN
ENDIF
IF(kepgr_in.and.kepgr_ou)THEN ! twelfe cases
! starting from KEP
   IF(coox.eq.'KEP')THEN
      IF(cooy.eq.'EQU')THEN
         IF(PRESENT(del))THEN
            y=equ_kep(x,del)
         ELSE
            y=equ_kep(x)
         ENDIF
      ELSEIF(cooy.eq.'COM')THEN
         IF(PRESENT(del))THEN
            y=com_kep(x,fail_flag,del)
         ELSE
            y=com_kep(x,fail_flag)
         ENDIF
      ELSEIF(cooy.eq.'COT')THEN
         IF(PRESENT(del))THEN
            y=cot_kep(x,fail_flag,del)
         ELSE
            y=cot_kep(x,fail_flag)
         ENDIF
      ENDIF
      done=.true.
      RETURN
   ENDIF
! starting from EQU
   IF(coox.eq.'EQU')THEN
      IF(cooy.eq.'KEP')THEN
         IF(PRESENT(del))THEN
            y=kep_equ(x,fail_flag,del) 
            IF(fail_flag.gt.0) CALl eye(ndim,del)! del non exi, set to identity
         ELSE
            y=kep_equ(x,fail_flag)
         ENDIF
      ELSEIF(cooy.eq.'COM')THEN
         IF(PRESENT(del))THEN
            z=kep_equ(x,fail_flag,del1)
            IF(fail_flag.eq.0)THEN
               y=com_kep(z,fail_flag,del2)
               del=MATMUL(del2,del1)
            ELSE
               y=com_kep(z,fail_flag)
               del=del1 ! del2 not existent, set to identity 
            ENDIF
         ELSE
            z=kep_equ(x,fail_flag) ! always existent, maybe not C1
            y=com_kep(z,fail_flag) 
         ENDIF
      ELSEIF(cooy.eq.'COT')THEN
         IF(PRESENT(del))THEN
            z=kep_equ(x,fail_flag,del1)
            IF(fail_flag.eq.0)THEN
               y=cot_kep(z,fail_flag,del2)
               del=MATMUL(del2,del1)
            ELSE
               y=cot_kep(z,fail_flag)
               del=del1 ! del2 not existent, set to identity 
            ENDIF
         ELSE
            z=kep_equ(x,fail_flag) ! always existent, maybe not C1
            y=cot_kep(z,fail_flag) 
         ENDIF
      ENDIF
      done=.true.
      RETURN
   ENDIF
! starting from COM
   IF(coox.eq.'COM')THEN
      IF(cooy.eq.'KEP')THEN
         IF(PRESENT(del))THEN
            y=kep_com(x,fail_flag,del)
         ELSE
            y=kep_com(x,fail_flag)
         ENDIF
      ELSEIF(cooy.eq.'EQU')THEN
         IF(PRESENT(del))THEN
            z=kep_com(x,fail_flag,del1)
            IF(fail_flag.lt.4)THEN
               y=equ_kep(z,del2)
               del=MATMUL(del2,del1)
            ELSE
               y=z
               del=del1 ! del 2 non existent, set to identity
            ENDIF
         ELSE
            z=kep_com(x,fail_flag)
            IF(fail_flag.lt.4)THEN
               y=equ_kep(z)
            ELSE
               y=z
            ENDIF
         ENDIF
      ELSEIF(cooy.eq.'COT')THEN
!         WRITE(*,*)coox,cooy,x
         IF(PRESENT(del))THEN
            y=cot_com(x,fail_flag,del)
         ELSE
            y=cot_com(x,fail_flag)
         ENDIF
!      WRITE(*,*)coox,cooy, y
      ENDIF
      done=.true.
!      WRITE(*,*)coox,cooy,done, y
      RETURN
   ENDIF
! starting from COT
   IF(coox.eq.'COT')THEN
      IF(cooy.eq.'KEP')THEN
         IF(PRESENT(del))THEN
            y=kep_cot(x,fail_flag,del)
         ELSE
            y=kep_cot(x,fail_flag)
         ENDIF
      ELSEIF(cooy.eq.'COM')THEN
         IF(PRESENT(del))THEN
            y=com_cot(x,fail_flag,del) 
            IF(fail_flag.gt.0) CALl eye(ndim,del)! del non exi, set to identity
         ELSE
            y=com_cot(x,fail_flag)
         ENDIF
      ELSEIF(cooy.eq.'EQU')THEN
         IF(PRESENT(del))THEN
            z=kep_cot(x,fail_flag,del1)
            IF(fail_flag.lt.4)THEN
               y=equ_kep(z,del2)
               del=MATMUL(del2,del1)
            ELSE
               y=z
               del=del1 ! del 2 non existent, set to identity
            ENDIF
         ELSE
            z=kep_cot(x,fail_flag)
            IF(fail_flag.lt.4)THEN
               y=equ_kep(z)
            ELSE
               y=z
            ENDIF
         ENDIF
      ENDIF
      done=.true.
      RETURN
   ENDIF
   IF(.not.done)THEN
      PRINT *, ' coo_cha: this should not happen 2, coox=',coox,' cooy=', cooy
      WRITE(*,*) x
      STOP
   ENDIF
ENDIF
END SUBROUTINE coo_cha
! ====================================================================
! UTILITIES
! ====================================================================
LOGICAL FUNCTION equal_orbels(el1,el2)
  TYPE(orbit_elem), INTENT(IN) :: el1,el2
  INTEGER i
  equal_orbels=.false.
  IF(el1%coo.ne.el2%coo.or.el1%ndim.ne.el2%ndim.or.el1%center.ne.el2%center)RETURN
  DO i=1,el1%ndim
     IF(el1%coord(i).ne.el2%coord(i))RETURN
  ENDDO
  equal_orbels=.true.
!magnitude does not matter!
END FUNCTION equal_orbels

SUBROUTINE eldiff(el1,el0,del)
  TYPE(orbit_elem), INTENT(IN) :: el1,el0
  DOUBLE PRECISION, DIMENSION(el0%ndim) :: del
  DOUBLE PRECISION pridif
  IF(el1%coo.ne.el0%coo)THEN
     WRITE(ierrou,*)'eldiff: different coordinates ',el1%coo, el0%coo
     numerr=numerr+1
     del=0.d0
     RETURN
  ENDIF
  IF(el0%coo.eq.'CAR')THEN
     del=el1%coord-el0%coord
  ELSEIF(el0%coo.eq.'EQU')THEN
     del(1:5)=el1%coord(1:5)-el0%coord(1:5)
     del(6)=pridif(el1%coord(6),el0%coord(6))
  ELSEIF(el0%coo.eq.'KEP')THEN
     del(1:3)=el1%coord(1:3)-el0%coord(1:3)
     del(4)=pridif(el1%coord(4),el0%coord(4))
     del(5)=pridif(el1%coord(5),el0%coord(5))
     del(6)=pridif(el1%coord(6),el0%coord(6))
  ELSEIF(el0%coo.eq.'COM')THEN
     del(1:3)=el1%coord(1:3)-el0%coord(1:3)
     del(4)=pridif(el1%coord(4),el0%coord(4))
     del(5)=pridif(el1%coord(5),el0%coord(5))
     del(6)=el1%coord(6)-el0%coord(6)
  ELSEIF(el0%coo.eq.'ATT')THEN
     del(2:6)=el1%coord(2:6)-el0%coord(2:6)
     del(1)=pridif(el1%coord(1),el0%coord(1))
  ELSE
     WRITE(*,*)'eldiff: unknown coordinates ',el1%coo, el0%coo
     STOP
  ENDIF
END SUBROUTINE eldiff

! ====================================================================
! ECC_PERI
! computes eccentricty, pericenter, apocenter, mean motion
! if energy>0 then enne=0, qg=1.d+50
SUBROUTINE ecc_peri(el,ecc,q,qg,enne)
TYPE(orbit_elem), INTENT(IN) :: el
DOUBLE PRECISION, INTENT(OUT) :: ecc,q,qg,enne
! ===================================================================
DOUBLE PRECISION gm,eps,a
DOUBLE precision x(3),y(3),ang(3),vlenz(3),gei,gei2,r0,vel2,alpha,prscal,vsize
TYPE(orbit_elem) el1
INTEGER fail_flag
!
eps=100*epsilon(1.d0)
IF(rhs.eq.1)THEN
   gm=centerbody_mass(0) ! we want always eccentricity heliocentric
ELSEIF(rhs.eq.2)THEN
   gm=centerbody_mass(3) ! we want always eccentricity geocentric
ELSE
! WARNING: automatic choice of method/dynamical model not supported in ORBIT9
   WRITE(*,*)'ecc_peri: ORBIT9 not supported here '
   STOP
ENDIF 
IF(el%coo.eq.'EQU')THEN
   ecc=sqrt(el%coord(2)**2+el%coord(3)**2)
   q=el%coord(1)*(1.d0-ecc)
   qg=el%coord(1)*(1.d0+ecc)
   enne=sqrt(gm/el%coord(1)**3)
ELSEIF(el%coo.eq.'KEP')THEN
   ecc=el%coord(2)
   q=el%coord(1)*(1.d0-ecc)
   qg=el%coord(1)*(1.d0+ecc)
   enne=sqrt(gm/el%coord(1)**3)
ELSEIF(el%coo.eq.'COM'.or.el%coo.eq.'COT')THEN
   ecc=el%coord(2)
   q=el%coord(1)
   IF(ecc.gt.1.d0-eps)THEN
      qg=1.d+50
      enne=0.d0
   ELSE
      a=el%coord(1)/(1.d0-ecc)
      qg=a*(1.d0+ecc)
      enne=sqrt(gm/a**3)
   ENDIF
ELSE
   IF(el%coo.eq.'ATT')THEN
! convert to cartesian if it is ATT
      CALL coo_cha(el,'CAR',el1,fail_flag)
   ELSEIF(el%coo.eq.'CAR')THEN
      el1=el
   ELSE
      WRITE(*,*)' ecc_peri: unknown coord type, coo= ',el%coo
      WRITE(*,*) el
      STOP
   ENDIF
! compute eccentricity
   x=el1%coord(1:3)
   y=el1%coord(4:6)
!  radius and velocity squared                                          
   vel2=prscal(y,y)          ! velocity
   r0=vsize(x)               ! distance from center
   call prvec(x,y,ang)       !  angular momentum 
! non singular first element
   gei2=prscal(ang,ang)
   gei=sqrt(gei2)
   IF(gei.eq.0.d0)THEN
      WRITE(*,*) ' ecc_peri: zero angular momentum ',el
   ENDIF
   call prvec(y,ang,vlenz) ! Lenz vector
   vlenz=vlenz*(1.d0/gm)-x(1:3)*(1.d0/r0)
   ecc=vsize(vlenz)
   q=gei2/(gm*(1.d0+ecc))
   alpha=vel2-2*gm/r0 ! 2* energy 
   IF(ecc.gt.1.d0-eps)THEN
      qg=1.d+50
      enne=0.d0
   ELSE
      qg=gei2/(gm*(1.d0-ecc))
      a=-gm/alpha
      enne=sqrt(gm/a**3)
   ENDIF
ENDIF
END SUBROUTINE ecc_peri

DOUBLE PRECISION FUNCTION centerbody_mass(iplanet)
USE planet_masses
INTEGER,INTENT(IN) :: iplanet
IF(iplanet.eq.0)THEN
! WARNING: does not work in ORBIT9
   IF(rhs.gt.2)THEN
      WRITE(*,*)' centerbody_mass: not ready for rhs=', rhs
      STOP
   ENDIF
   centerbody_mass=gms
ELSEIF(iplanet.eq.3)THEN
   centerbody_mass=gmearth
ELSEIF(iplanet.le.npla+iatrue.and.iplanet.ne.3)THEN
   centerbody_mass=gm(iplanet)
ELSE
   WRITE(*,*)' centerbody_mass: center not available ', iplanet
   STOP
ENDIF 
END FUNCTION centerbody_mass
!
TYPE(orbit_elem) FUNCTION coord_translate(el,iplanet)
TYPE(orbit_elem), INTENT(IN) :: el
INTEGER, INTENT(IN) :: iplanet ! planet; JPL ephem codes
DOUBLE PRECISION :: xpla(6)
IF(el%coo.ne.'CAR')THEN
   WRITE(*,*)' coord_translate: wrong input coord ',el
   STOP
ENDIF
coord_translate=el
! check planet, observers are on Earth! 
IF(iplanet.eq.3.and.el%center.eq.0)THEN
   CALL earcar(el%t,xpla,1)
   coord_translate%coord=el%coord-xpla
   coord_translate%center=3
ELSEIF(iplanet.eq.0.and.el%center.eq.3)THEN
   CALL earcar(el%t,xpla,1)
   coord_translate%coord=el%coord+xpla
   coord_translate%center=0
ELSE
   WRITE(*,*)' coord_translate: not yet available for ', iplanet, el%center
   STOP
ENDIF
END FUNCTION coord_translate
! ====================================================================
! END UTILITIES
! ====================================================================
! CARTESIAN TO/FROM ELEMENTS
! ====================================================================
! TO_ELEMS
! ====================================================================
TYPE(orbit_elem) FUNCTION to_elems(el,cooy,fail_flag,del)
USE ever_pitkin
TYPE(orbit_elem),INTENT(IN):: el ! must be cartesian
CHARACTER*3, INTENT(IN) :: cooy ! can be 'KEP', 'EQU', 'COM', 'COT'
INTEGER, INTENT(OUT) :: fail_flag ! error flag
! optional jacobian matrix. in output: derivatives d(to_elems)/d(el)
DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim),INTENT(OUT),OPTIONAL :: del
DOUBLE PRECISION gm
DOUBLE PRECISION, DIMENSION(3) :: x0, y0
DOUBLE PRECISION, DIMENSION(3) :: x, y, ang ! position, velocity, ang. momentum
DOUBLE PRECISION, DIMENSION(3) :: vlenz, fb, gb, ww ! Lenz vector, broucke f and g, pos x vlenz
DOUBLE PRECISION :: gei2, gei, ecc2, ecc, varpi, div, tgi2, omnod
DOUBLE PRECISION :: prscal, vsize, princ, eps
DOUBLE PRECISION :: r0, vel2, alpha, sig0, peri, psi, dt
DOUBLE PRECISION :: enne, rad, beta, chk, ch, ck, fe, x2, y2, cosf, sinf
DOUBLE PRECISION :: eq(6), cosv,sinv,dd
TYPE(orbit_elem) :: elback ! for inversion
DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim) :: del1
INTEGER :: fail_flag1, ising
DOUBLE PRECISION :: det
! paranoia check of the calling arguments
IF(el%coo.ne.'CAR')THEN
   WRITE(*,*)' to_elems: wrong coordinates in input',x
   fail_flag=20
   RETURN
ENDIF
IF(cooy.ne.'COM'.and.cooy.ne.'KEP'.and.cooy.ne.'EQU'.and.cooy.ne.'COT')THEN
   WRITE(*,*)' to_elems: wrong output coordinates', cooy
   fail_flag=20
   RETURN
ENDIF
! initializations
fail_flag=0
gm=centerbody_mass(el%center) ! mass of central body
eps=100*epsilon(1.d0)         ! small number for rounding off
to_elems=el                   ! copy defaults
! =====================================================================
! portion of the algorithm used for all output coordinates
! =====================================================================
x=el%coord(1:3)               ! cartesian position
y=el%coord(4:6)               ! cartesian velocity
vel2=prscal(y,y)              ! velocity squared
r0=vsize(x)                   ! distance from center
CALL prvec(x,y,ang)           ! angular momentum vector 
gei2=prscal(ang,ang)
gei=sqrt(gei2)                ! angular momentum scalar
IF(gei.eq.0.d0)THEN 
   WRITE(*,*) ' to_elems: zero angular momentum ',el 
   fail_flag=10
   RETURN
ENDIF
CALL prvec(y,ang,vlenz)       ! Lenz vector 
vlenz=vlenz*(1.d0/gm)-x(1:3)*(1.d0/r0) 
ang=ang*(1.d0)/gei            ! orbit normal unit vector
!   zero divide occurs for inclination of 180 degrees                   
div=1.d0+ang(3) 
IF(div.eq.0.d0)THEN 
   WRITE(*,*) ' to_elems: 180 deg. inclination ',el
   fail_flag=9
   RETURN
ENDIF   
!  unit vectors of the equinoctal reference system (Broucke and         
!  Cefola 1972, CM 5, 303--310) are fb, gb, ang                           
fb(1)=1.d0-ang(1)**2/div 
fb(2)=-ang(1)*ang(2)/div 
fb(3)=-ang(1) 
CALL prvec(ang,fb,gb) 
!  elements related to eccentricity and inclination                     
eq(2)=prscal(vlenz,gb) 
eq(3)=prscal(vlenz,fb) 
eq(4)=ang(1)/div 
eq(5)=-ang(2)/div 
ecc2=eq(2)**2+eq(3)**2 
ecc=sqrt(ecc2)                ! eccentricity
IF(ecc.lt.eps)THEN
   varpi=0.d0
   fail_flag=1
ELSE
   varpi=atan2(eq(2),eq(3))   ! longitude of pericenter
ENDIF                                          
tgi2=sqrt(eq(4)**2+eq(5)**2)  ! tanget of inclination/2
if(tgi2.lt.eps)then 
   omnod=0.d0 
   fail_flag=fail_flag+2
else 
   omnod=princ(atan2(eq(4),eq(5))) ! Omega (longitude of node) 
endif
alpha=vel2-2*gm/r0            ! 2* energy
IF(alpha.gt.-eps)THEN
  IF(cooy.eq.'EQU'.or.cooy.eq.'KEP') fail_flag=5 ! hyperbolic orbit; is a failure
! if requested coordinates in output are KEP, EQU
ENDIF
sig0=prscal(x,y)              ! scalar product as used in ever_pitkin
! =====================================================================
! part depending upon output coordinates
! =====================================================================
IF(cooy.eq.'COM'.or.cooy.eq.'COT'.or.fail_flag.eq.5)THEN
!         coord(1)=q=a(1-e)
!         coord(2)=e
!         coord(3)=I
!         coord(4)=Omega 
!         coord(5)=omega 
! if cooy=COM
!         coord(6)=t_0 (time at pericenter)
! elseif cooy=COT
!         coord(6)=v (true anomaly)
   eq(1)=gei2/(gm*(1.d0+ecc))
   eq(2)=ecc
   eq(5)=princ(varpi-omnod)
   eq(3)=2*atan(tgi2)
   eq(4)=omnod
   peri=eq(1)
   IF(cooy.eq.'COM')THEN
      IF(ecc.gt.eps)THEN
! WARNING: use of solve_peri may fail, add error flag which affects
!          the output fail_flag
         CALL solve_peri(r0,sig0,peri,gm,alpha,psi,dt)
         eq(6)=el%t+dt
         fail_flag=0
      ELSE
         fail_flag=1  ! numerically circular orbit
         eq(6)=el%t   ! always at perihelion
      ENDIF
      to_elems%coo='COM' 
   ELSEIF(cooy.eq.'COT'.or.fail_flag.eq.5)THEN
      dd=r0*ecc
      cosv=prscal(x,vlenz)/dd  ! cosine of true anomaly
      CALL prvec(vlenz,x,ww)
      sinv=prscal(ang,ww)/dd   ! sine of true anomaly   
      eq(6)=atan2(sinv,cosv)
      IF(eq(6).lt.0.d0) eq(6)=eq(6)+dpig
      to_elems%coo='COT'
   ENDIF
ELSE
   eq(1)=-gm/alpha
! selection of different methods
   IF(ecc.gt.eccbou)THEN
      peri=eq(1)*(1-ecc)
      enne=dsqrt(gm/eq(1)**3)
      CALL solve_peri(r0,sig0,peri,gm,alpha,psi,dt)
      eq(6)=princ(-dt*enne+varpi)
   ELSE
! mean longitude from non--singular Kepler equation                   
      rad=dsqrt(1.d0-ecc2) 
      beta=1.d0/(1.d0+rad) 
      chk=eq(2)*eq(3)*beta 
      ch=1.d0-eq(2)**2*beta 
      ck=1.d0-eq(3)**2*beta 
      x2=prscal(x,fb) 
      y2=prscal(x,gb) 
      cosf=eq(3)+(ck*x2-chk*y2)/(eq(1)*rad) 
      sinf=eq(2)+(ch*y2-chk*x2)/(eq(1)*rad) 
! eccentric longitude                                                 
      fe=datan2(sinf,cosf) 
      eq(6)=fe+eq(2)*cosf-eq(3)*sinf 
! reduction to principal value  
      eq(6)=princ(eq(6))
   ENDIF
   IF(cooy.eq.'KEP'.and.fail_flag.eq.0)THEN
      eq(2)=ecc
      eq(4)=omnod
      eq(3)=2*atan(tgi2)
      eq(5)=princ(varpi-omnod)  
      eq(6)=princ(eq(6)-varpi)
      to_elems%coo='KEP'
   ELSEIF(cooy.eq.'EQU'.or.(fail_flag.ge.1.and.fail_flag.le.3))THEN
!          coord(2)= h = e sin(varpi) ; varpi = Omega + omega
!          coord(3)= k = e cos(varpi)
!          coord(4)= p = tg(I/2) sin(Omega)
!          coord(5)= q = tg(I/2) cos(Omega)
!          coord(6)= mean longitude, measured from fb vector (mean an + varpi)
      fail_flag=0
      to_elems%coo='EQU'
   ELSE
      PRINT *, ' to_elems: this should not happen ', cooy, fail_flag
   ENDIF
ENDIF
to_elems%coord=eq
IF(PRESENT(DEL))THEN
! compute jacobian matrix by using the inverse transformation
! this is dirty, but the derivatives of elems w.r. to cartesian
! are not as important as the inverse
elback=car_elems(to_elems,fail_flag1,del1)
! WARNING: how to use fail_flag1??????
! inversion
CALL matin(del1,det,6,0,6,ising,1)
IF(ising.eq.0)THEN
   del=del1
ELSE
   del=0.d0
! WARNING: how to set fail_flag??? as above????
ENDIF
ENDIF
END FUNCTION to_elems
! ====================================================================
! CAR_ELEMS
! ====================================================================
TYPE(orbit_elem) FUNCTION car_elems(el,fail_flag,del)
USE ever_pitkin
TYPE(orbit_elem),INTENT(IN):: el ! can be 'KEP', 'EQU', 'COM', 'COT'
INTEGER, INTENT(OUT) :: fail_flag
! optional jacobian matrix. in output: derivatives d(car)/d(el)
DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim),INTENT(OUT),OPTIONAL :: del
! END INTERFACE
DOUBLE PRECISION :: gm ! center body mass
CHARACTER*(3) :: coox
TYPE(orbit_elem) :: equ ! used for equinoctal, to remove the KEP case 
DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim) :: del1,del0,del2,del3
DOUBLE PRECISION :: eq(6) ! coordinate 6d-vectors
DOUBLE PRECISION, DIMENSION(3) :: x, y, x0, y0, acc
DOUBLE PRECISION :: princ, vsize, eps ! functions, small no
DOUBLE PRECISION :: r0, v0, dt
DOUBLE PRECISION, DIMENSION(3) :: fb,gb,dfdp,dfdq,dgdp,dgdq ! Brouke's vectors
DOUBLE PRECISION :: tgim2, tgim, ecc, upq, ecc2, varpi, enne, mean_an
! vectors used in computing derivatives (not the complete derivatives)
DOUBLE PRECISION :: dedh, dedk, ddigdh, ddigdk, fact, rsvperi,q, upecv2
DOUBLE PRECISION cosn,sinn,coso,sino,cosi,sini, cosv,sinv,rsq,fgfact,eye3(3,3)
! paranoia check of input arguments
coox=el%coo
IF(coox.ne.'COM'.and.coox.ne.'KEP'.and.coox.ne.'EQU'.and.coox.ne.'COT')THEN
   WRITE(*,*)' car_elems: wrong input coordinates', coox
   fail_flag=20
   RETURN
ENDIF
! initializations
fail_flag=0
gm=centerbody_mass(el%center)
eps=100*epsilon(1.d0) ! small number
car_elems=el ! copy other data (e.g., magnitude)
! ================================================================
! the position and velocity at pericenter are computed for all
!  output coordiantes
! ================================================================ 
IF(coox.eq.'COM'.or.coox.eq.'COT')THEN
   ecc=el%coord(2)
   ecc2=ecc**2
   IF(ecc.lt.eps)THEN
      varpi=0.d0
      fail_flag=1
   ELSE
      varpi=el%coord(4)+el%coord(5)
   ENDIF   
   IF(el%coord(3).lt.eps) fail_flag=fail_flag+2
! pericenter  
   r0=el%coord(1)
   v0=sqrt(gm*(1.d0+ecc)/el%coord(1))
!   WRITE(*,*) ' r0,v0 ',r0,v0
   cosn=cos(el%coord(4))
   sinn=sin(el%coord(4))
   coso=cos(el%coord(5))
   sino=sin(el%coord(5))
   cosi=cos(el%coord(3))
   sini=sin(el%coord(3))
   x0=(/coso*cosn-sino*sinn*cosi, coso*sinn+sino*cosn*cosi, sino*sini/)*r0
   y0=(/-sino*cosn-coso*sinn*cosi, -sino*sinn+coso*cosn*cosi, coso*sini/)*v0
   IF(coox.eq.'COM')THEN
      dt=el%t-el%coord(6) ! time between peihelion and epoch, for 2-body propagation
   ELSEIF(coox.eq.'COT')THEN
      cosv=cos(el%coord(6)) ! cosine and sine of true anomaly, for 2-body formulae
      sinv=sin(el%coord(6))
   ENDIF
ELSE
   IF(el%coo.eq.'KEP')THEN
      IF(PRESENT(del))THEN
         equ=equ_kep(el,del1)
      ELSE
         equ=equ_kep(el)
      ENDIF
   ELSEIF(el%coo.eq.'EQU')THEN
      equ=el
      IF(PRESENT(del))CALL eye(el%ndim,del1)
   ENDIF
   eq=equ%coord
   ecc2=eq(2)**2+eq(3)**2
   ecc=sqrt(ecc2)
! pericenter
   r0=eq(1)*(1.d0-ecc)
   v0=sqrt(gm*(1.d0+ecc)/(eq(1)*(1.d0-ecc)))
   enne=sqrt(gm/eq(1)**3)
!   WRITE(*,*) ' r0,v0, enne ',r0,v0,enne
! varpi= Omega+ omega
   IF(ecc.lt.eps)THEN
      varpi=0.d0
      fail_flag=1
   ELSE
      varpi=atan2(eq(2),eq(3))
   ENDIF
   mean_an=princ(eq(6)-varpi)
   IF(mean_an.gt.pig)mean_an=mean_an-dpig
   dt=mean_an/enne
   tgim2=eq(4)**2+eq(5)**2 
   tgim=dsqrt(tgim2)
   IF(tgim.lt.eps) fail_flag=fail_flag+2                                     
   upq=1.d0+tgim2 
!  unit vectors of the equinoctal reference system (Broucke and         
!  Cefola 1972, CM 5, 303--310) are f, g and the angular momentum       
!  unit vector                                                          
   fb(1)=(1.d0-eq(4)**2+eq(5)**2)/upq 
   fb(2)=2*eq(4)*eq(5)/upq 
   fb(3)=-2*eq(4)/upq 
   gb(1)=2*eq(4)*eq(5)/upq 
   gb(2)=(1.d0+eq(4)**2-eq(5)**2)/upq 
   gb(3)=2*eq(5)/upq    
! pericenter vector, parallel to Lenz vector
   x0=r0*(fb*cos(varpi)+gb*sin(varpi))
! velocity at pericenter
   y0=v0*(-fb*sin(varpi)+gb*cos(varpi))
ENDIF
! WARNING: for now the use of f-g series is avoided only for coox='COT'
IF(coox.eq.'COT')THEN
   q=el%coord(1)
   rsq=(1.d0+ecc)/(1.d0+ecc*cosv) ! r/q
   rsvperi=sqrt((1.d0+ecc)*q**3/gm)/(1+ecc*cosv)
   x=x0*(rsq*cosv)+y0*(rsvperi*sinv)
   fgfact=sqrt(gm/(q**3*(1.d0+ecc)))
   y=(-sinv*fgfact)*x0+((cosv+ecc)/(1.d0+ecc))*y0
ELSE
   IF(PRESENT(del))THEN
      CALL fser_propag_der(x0,y0,0.d0,dt,gm,x,y,del3)
   ELSE
      CALL fser_propag(x0,y0,0.d0,dt,gm,x,y)
   ENDIF
ENDIF
car_elems%coord(1:3)=x
car_elems%coord(4:6)=y
car_elems%coo='CAR'
IF(PRESENT(del))THEN

! compute partials of x0, y0 w.r. to elements eq(1:5)

   IF(coox.eq.'COM'.or.coox.eq.'COT')THEN
      del0=0.d0
      del0(1:3,1)=  x0*(1/el%coord(1))
      del0(4:6,1)=  y0*(-0.5/el%coord(1))
      del0(1:3,2)=  0.d0
      del0(4:6,2)=  y0*0.5/(1.d0+ecc)
      del0(1:3,3)=  (/ sinn*sino*sini, -cosn*sino*sini, sino*cosi/)*r0
      del0(4:6,3)=  (/ sinn*coso*sini, -cosn*coso*sini, coso*cosi/)*v0
      del0(1:3,4)=  (/-sinn*coso-cosn*sino*cosi, cosn*coso-sinn*sino*cosi, 0.d0/)*r0
      del0(4:6,4)=  (/sinn*sino-cosn*coso*cosi, -cosn*sino-sinn*coso*cosi, 0.d0/)*v0
      fact= r0/v0
      del0(1:3,5)=  fact*y0
      del0(4:6,5)= -(1.d0/fact)*x0  
      IF(coox.eq.'COM')THEN    
! now use Goodyear to compute jacobian of propagation from t=eq(6) to t=eqn%t
         del=MATMUL(del3,del0)
! then fix derivatives with respect to time of pericenter
         acc=-x*gm/vsize(x)**3
         del(1:3,6)=  -y
         del(4:6,6)=  -acc
      ELSEIF(coox.eq.'COT')THEN
         del3=0.d0 ! del3=\partial (x,y)/\partial (x0,y0)
         CALL eye(3,eye3)
         del3(1:3,1:3)=rsq*cosv*eye3
         del3(1:3,4:6)=rsvperi*sinv*eye3
         del3(4:6,1:3)=-fgfact*sinv*eye3
         del3(4:6,4:6)=((cosv+ecc)/(1.d0+ecc))*eye3
         del2=0.d0 ! del2= \partial (x,y)/\partial el%coord (not through x0,y0)
         upecv2=(1.d0+ecc*cosv)**(-2)
         del2(1:3,1)=1.5d0*sqrt(q*(1.d0+ecc)/gm)*sinv/(1+ecc*cosv) * y0 ! \partial x/\partial q
         del2(1:3,2)=cosv*(1.d0-cosv)*upecv2 * x0    &! \partial x/\partial ecc
&            +0.5d0*sqrt(q**3/((1.d0+ecc)*gm))*sinv*(1.d0-ecc*cosv-2.d0*cosv)*upecv2 * y0
         del2(1:3,6)=-sinv*(1.d0+ecc)*upecv2 * x0  & ! \partial x/\partial v
&             +sqrt(q**3*(1.d0+ecc)/gm)*(cosv+ecc)*upecv2* y0
         del2(4:6,1)=1.5d0*sqrt(gm/(q**5*(1.d0+ecc)))*sinv * x0 ! \partial y/\partial q
         del2(4:6,2)=0.5d0*sqrt(gm/(q**3*(1.d0+ecc)**3))*sinv * x0  & ! \partial y/\partial e
&             + (1.d0-cosv)/(1.d0+ecc)**2 * y0 
         del2(4:6,6)= -sqrt(gm/(q**3*(1.d0+ecc)))*cosv * x0 - sinv/(1.d0+ecc) * y0 ! \partial y/\partial v
         del=MATMUL(del3,del0)+del2 ! chain rule
      ENDIF
   ELSEIF(coox.eq.'EQU'.or.coox.eq.'KEP')THEN
      del2=0.d0
      del2(1:3,1)=  x0*(1/eq(1))
      del2(4:6,1)=  y0*(-0.5/eq(1))
      dedh=eq(2)/ecc
      dedk=eq(3)/ecc
      ddigdh=eq(3)/ecc2
      ddigdk=-eq(2)/ecc2
      fact=r0/v0
      del2(1:3,2)= -dedh/(1.d0-ecc)*x0+fact*ddigdh*y0
      del2(1:3,3)= -dedk/(1.d0-ecc)*x0+fact*ddigdk*y0
!      del2(1:3,2)=  eq(1)*(-gb+ddigdh*(-sin(varpi)*fb+cos(varpi)*gb))
      del2(1:3,3)= -dedk/(1.d0-ecc)*x0+fact*ddigdk*y0
!      del2(1:3,3)=  eq(1)*(-fb+ddigdk*(-sin(varpi)*fb+cos(varpi)*gb))
      del2(4:6,2)=  dedh/(1.d0-ecc2)*y0-(1.d0/fact)*ddigdh*x0
      del2(4:6,3)=  dedk/(1.d0-ecc2)*y0-(1.d0/fact)*ddigdk*x0
      dfdp(1)=-2*eq(4)
      dfdp(2)=2*eq(5)
      dfdp(3)=-2.d0
      dfdp=(dfdp-fb*(2*eq(4)))/upq
      dfdq(1)=2*eq(5)
      dfdq(2)=2*eq(4)
      dfdq(3)=0.d0
      dfdq=(dfdq-fb*(2*eq(5)))/upq
      dgdp(1)=2*eq(5)
      dgdp(2)=2*eq(4)
      dgdp(3)=0.d0
      dgdp=(dgdp-gb*(2*eq(4)))/upq
      dgdq(1)=2*eq(4)
      dgdq(2)=-2*eq(5)
      dgdq(3)=2.d0
      dgdq=(dgdq-gb*(2*eq(5)))/upq
      del2(1:3,4)=  r0*(cos(varpi)*dfdp+sin(varpi)*dgdp)
      del2(1:3,5)=  r0*(cos(varpi)*dfdq+sin(varpi)*dgdq)
      del2(4:6,4)=  v0*(-sin(varpi)*dfdp+cos(varpi)*dgdp)
      del2(4:6,5)=  v0*(-sin(varpi)*dfdq+cos(varpi)*dgdq)
      del2(1:6,6)=  0.d0
      del0=MATMUL(del2,del1)
! now use Goodyear to compute jacobian of propagation from t=eq(6) to t=eqn%t
      del=MATMUL(del3,del0)
! then fix derivatives with respect to mean anomaly
      acc=-x*gm/vsize(x)**3
      del(1:3,6)=y/enne
      del(4:6,6)=acc/enne
! then correct derivatives with respect to a to account for enne in dt
      del(1:3,1)= del(1:3,1)+(1.5d0*mean_an/(enne*eq(1)))*y
      del(4:6,1)= del(4:6,1)+(1.5d0*mean_an/(enne*eq(1)))*acc
! then correct derivatives with respedct to h,k to account for varpi in dt
      IF(coox.eq.'EQU')THEN
         del(1:3,2)=del(1:3,2)-(ddigdh/enne)*y
         del(1:3,3)=del(1:3,3)-(ddigdk/enne)*y
         del(4:6,2)=del(4:6,2)-(ddigdh/enne)*acc
         del(4:6,3)=del(4:6,3)-(ddigdk/enne)*acc
      ENDIF
   ENDIF
ENDIF

END FUNCTION car_elems
! ====================================================================
! ELEMENTS INTERNAL CONVERSIONS
! ====================================================================
! KEP_EQU
! keplerian elements from equinoctal
! =======================================================
TYPE(orbit_elem) FUNCTION kep_equ(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be equinoctal
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION dig,ecc,ecc2,tgi2,tgi2s,princ,eps
DOUBLE PRECISION,DIMENSION(6) :: eq
!
IF(el%coo.ne.'EQU')THEN
   WRITE(*,*)' kep_equ: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
eps=100*epsilon(1.d0) ! negligible number
eq=el%coord
kep_equ=el  !default
kep_equ%coord(1)=eq(1) !semimajor axis
!  test on eccentricity                                                 
ecc2=eq(2)**2+eq(3)**2
ecc=sqrt(ecc2) 
if(ecc.lt.eps)then 
   fail_flag=1 ! but computation is completed
   dig=0.d0 
else 
   dig=atan2(eq(2),eq(3)) 
endif
kep_equ%coord(2)=ecc 
!   test on tangent of half inclination                                 
tgi2s=eq(4)**2+eq(5)**2
tgi2=sqrt(tgi2s) 
if(tgi2.lt.eps)then 
   fail_flag=fail_flag+2 ! but computation is completed
   kep_equ%coord(4)=0.d0 
else 
   kep_equ%coord(4)=atan2(eq(4),eq(5)) 
endif
kep_equ%coord(3)=2.d0*atan(tgi2) 
!   angular variables                                                   
kep_equ%coord(5)=dig-kep_equ%coord(4) 
kep_equ%coord(6)=eq(6)-dig 
kep_equ%coord(4)=princ(kep_equ%coord(4)) 
kep_equ%coord(5)=princ(kep_equ%coord(5)) 
kep_equ%coord(6)=princ(kep_equ%coord(6)) 
kep_equ%coo='KEP'
IF(PRESENT(del).and.fail_flag.eq.0)THEN
   del=0.d0 ! jacobian is not usable in the nearly singular case
   del(1,1)=  1.d0 
   del(2,2)=  eq(2)/ecc ! ecc=h^2+k^2
   del(2,3)=  eq(3)/ecc
   del(3,4)=  2*eq(4)/((1+tgi2**2)*tgi2) 
   del(3,5)=  2*eq(5)/((1+tgi2**2)*tgi2)
   del(4,4)=  eq(5)/tgi2s ! Omega=atan2(p,q)
   del(4,5)= -eq(4)/tgi2s
   del(5,2)=  eq(3)/ecc2  ! dig=atan2(h,k)
   del(5,3)= -eq(2)/ecc2  ! omega=dig-Omega 
   del(5,4)= -eq(5)/tgi2s
   del(5,5)=  eq(4)/tgi2s
   del(6,2)= -eq(3)/ecc2  ! ell=lambda-dig
   del(6,3)=  eq(2)/ecc2
   del(6,6)=  1.d0
! else the jacobian is not defined!!! 
ENDIF

RETURN
END FUNCTION kep_equ

! =======================================================
! EQU_KEP
! equinoctal elements from keplerian
! =======================================================
TYPE(orbit_elem) FUNCTION equ_kep(el,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be keplerian
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION dig,ecc,tgim,princ
DOUBLE PRECISION sindig,cosdig,sinom,cosom,dti2
!
IF(el%coo.ne.'KEP')THEN
   WRITE(*,*)' equ_kep: wrong input coordinates',el
   STOP
ENDIF 
equ_kep=el  !default
equ_kep%coord(1)=el%coord(1) 
dig=el%coord(4)+el%coord(5) 
ecc=el%coord(2) 
sindig=sin(dig)
cosdig=cos(dig)
equ_kep%coord(2)=ecc*sindig 
equ_kep%coord(3)=ecc*cosdig 
tgim=tan(el%coord(3)/2.d0) 
sinom=sin(el%coord(4))
cosom=cos(el%coord(4))
equ_kep%coord(4)=tgim*sinom
equ_kep%coord(5)=tgim*cosom 
equ_kep%coord(6)=dig+el%coord(6) 
equ_kep%coord(6)=princ(equ_kep%coord(6))
equ_kep%coo='EQU' 
IF(PRESENT(del))THEN
   del=0.d0
   del(1,1)=1.d0 
   del(2,2)=sindig
   del(2,5)=equ_kep%coord(3)
   del(2,4)=equ_kep%coord(3) 
   del(3,2)=cosdig
   del(3,5)=-equ_kep%coord(2)
   del(3,4)=-equ_kep%coord(2)
   dti2=1.d0/(2*cos(el%coord(3)/2)**2) 
   del(4,3)=dti2*sinom 
   del(4,4)= equ_kep%coord(5)
   del(5,3)=dti2*cosom 
   del(5,4)=- equ_kep%coord(4)
   del(6,4)=1.d0 
   del(6,5)=1.d0 
   del(6,6)=1.d0 
ENDIF
RETURN
END FUNCTION equ_kep
! =======================================================
! COM_KEP
! cometary elelements from keplerian
! =======================================================
TYPE(orbit_elem) FUNCTION com_kep(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be keplerian
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION ecc, eps, gm, enne, princ
DOUBLE PRECISION,DIMENSION(6) :: eq
!
IF(el%coo.ne.'KEP')THEN
   WRITE(*,*)' com_kep: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
gm=centerbody_mass(el%center)
eps=100*epsilon(1.d0) ! negligible number
eq=el%coord
com_kep=el  !default
ecc=eq(2)
com_kep%coord(1)=eq(1)*(1.d0-ecc)
!  test on eccentricity                                                 
if(ecc.lt.eps) fail_flag= 1 ! but computation is completed
!   test on tangent of half inclination                                 
if(eq(3).lt.eps) fail_flag=fail_flag+2 ! but computation is completed
! time of pericenter from mean anomaly
enne=sqrt(gm/eq(1)**3)
eq(6)=princ(eq(6))
com_kep%coord(6)=el%t-eq(6)/enne 
com_kep%coo='COM'
IF(PRESENT(del).and.fail_flag.eq.0)THEN
   del=0.d0 ! jacobian is not usable in the nearly singular case
   del(1,1)=  1.d0-ecc
   del(1,2)=  -eq(1) 
   del(2,2)=  1.d0
   del(3,3)=  1.d0
   del(4,4)=  1.d0
   del(5,5)=  1.d0
   del(6,1)= -1.5d0*eq(6)/(enne*eq(1))
   del(6,6)=  -1.d0/enne
ENDIF
END FUNCTION com_kep

! =======================================================
! COT_KEP
! cometary elements with true anomaly from keplerian
! =======================================================
TYPE(orbit_elem) FUNCTION cot_kep(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be keplerian
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION ecc, eps, gm, enne, princ
DOUBLE PRECISION, DIMENSION(6) :: eq
DOUBLE PRECISION ell,u,du,cosu,sinu,cosv,sinv,dvdu,dvde,dudell,dude 
INTEGER, PARAMETER :: itmax=10
INTEGER i
!
IF(el%coo.ne.'KEP')THEN
   WRITE(*,*)' cot_kep: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
gm=centerbody_mass(el%center)
eps=100*epsilon(1.d0) ! negligible number
eq=el%coord
cot_kep=el  !default
ecc=eq(2)
cot_kep%coord(1)=eq(1)*(1.d0-ecc)
!  test on eccentricity                                                 
if(ecc.lt.eps) fail_flag= 1 ! but computation is completed
!   test on tangent of half inclination                                 
if(eq(3).lt.eps) fail_flag=fail_flag+2 ! but computation is completed
! time of pericenter from mean anomaly
!enne=sqrt(gm/eq(1)**3)
! Kepler's equation to find eccentric anomaly
u=pig
DO i=1,itmax
   cosu=cos(u);sinu=sin(u)
   ell=u-ecc*sinu
   du=(eq(6)-ell)/(1.d0-ecc*cosu)
   u=u+du
   IF(abs(du).lt.eps*10) GOTO 3
ENDDO
WRITE(ierrou,*)' cot_kep: non convergent Kepler equation, u,e,du:', u,ecc,du
3 CONTINUE
cosv=(cosu-ecc)/(1.d0-ecc*cosu)
sinv=sqrt(1.d0-ecc**2)*sinu/(1.d0-ecc*cosu)
cot_kep%coord(6)=atan2(sinv,cosv)
cot_kep%coord(6)=princ(cot_kep%coord(6))
cot_kep%coo='COT'
IF(PRESENT(del).and.fail_flag.eq.0)THEN
   del=0.d0 ! jacobian is not usable in the nearly singular case
   del(1,1)=  1.d0-ecc
   del(1,2)=  -eq(1) 
   del(2,2)=  1.d0
   del(3,3)=  1.d0
   del(4,4)=  1.d0
   del(5,5)=  1.d0
   dvdu=sqrt((1.d0+ecc)/(1.d0-ecc))*(cosv+1.d0)/(cosu+1.d0)
   dvde=tan(u/2.d0)*(cosv+1.d0)/sqrt((1.d0+ecc)*(1.d0-ecc)**3)
   dudell=1.d0/(1.d0-ecc*cosu)
   dude=sinu/(1.d0-ecc*cosu)
   del(6,6)=dvdu*dudell
   del(6,2)=dvde+dvdu*dude
ENDIF
END FUNCTION cot_kep
! =======================================================
! KEP_COT
! keplerian elelements from cometary with true anomaly
! =======================================================
TYPE(orbit_elem) FUNCTION kep_cot(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be cometary
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION ecc, eps, gm, enne, a, dldn
DOUBLE PRECISION,DIMENSION(6) :: eq
DOUBLE PRECISION r,x1,y1,sinu,cosu,u,delldu,dellde,dudv,dude,v2
!
IF(el%coo.ne.'COT')THEN
   WRITE(*,*)' kep_cot: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
gm=centerbody_mass(el%center)
eps=100*epsilon(1.d0) ! negligible number
eq=el%coord
kep_cot=el  !default
ecc=eq(2)
if(ecc.ge.1.d0-eps)THEN
! this conversion is impossible
   fail_flag=5
! in output the elements are left as 'COT', del= identity
   IF(PRESENT(del)) CALL eye(el%ndim,del)
   RETURN
ENDIF
a=eq(1)/(1.d0-ecc)
kep_cot%coord(1)=a
!  test on eccentricity                                                 
if(ecc.lt.eps) fail_flag= 1 ! but computation is completed
!   test on tangent of half inclination                                 
if(eq(3).lt.eps) fail_flag=fail_flag+2 ! but computation is completed
! eccentric anomaly from true anomaly
r=a*(1.d0-ecc**2)/(1.d0+ecc*cos(eq(6)))
x1=r*cos(eq(6)); y1=r*sin(eq(6))
sinu=y1/(a*sqrt(1.d0-ecc**2)); cosu=x1/a+ecc
u=atan2(sinu,cosu)
! kepler's equation (easy way)
kep_cot%coord(6)=u-ecc*sinu
kep_cot%coo='KEP' ! successful conversion
IF(PRESENT(del).and.fail_flag.eq.0)THEN
   del=0.d0 ! jacobian is not usable in the nearly singular case
   del(1,1)=  1.d0/(1.d0-ecc)
   del(1,2)=  +eq(1)/(1.d0-ecc)**2 
   del(2,2)=  1.d0
   del(3,3)=  1.d0
   del(4,4)=  1.d0
   del(5,5)=  1.d0
! derivatives of mean anomaly (composite through eccentric anomaly)
   delldu=1.d0-ecc*cosu
   dellde=-sinu
   v2=eq(6)/2.d0
   dudv=sqrt((1.d0-ecc)/(1.d0+ecc))*cos(u/2.d0)**2/cos(v2)**2
   dude=(-2*tan(v2)*cos(u/2.d0)**2)/sqrt((1.d0-ecc)*(1.d0+ecc)**3)
   del(6,6)=delldu*dudv
   del(6,2)=dellde + delldu*dude
ENDIF
END FUNCTION kep_cot
! =======================================================
! COT_COM
! cometary elements with true anomaly from cometary
! =======================================================
TYPE(orbit_elem) FUNCTION cot_com(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be keplerian
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION ecc, eps , gm
DOUBLE PRECISION delcar(6,6),dcarel(6,6)
TYPE(orbit_elem) elcar,elcot
DOUBLE PRECISION dt, alpha,v,dvdq,dvde,dvdt0, vel2
INTEGER fail_flag1, fail_flag2
! ======================================================
IF(el%coo.ne.'COM')THEN
   WRITE(*,*)' cot_com: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
gm=centerbody_mass(el%center)
eps=100*epsilon(1.d0) ! negligible number
! eq=el%coord
cot_com=el  !default
ecc=el%coord(2)
!  test on eccentricity                                                 
if(ecc.lt.eps) fail_flag= 1 ! but computation is completed
!   test on tangent of half inclination                                 
if(el%coord(3).lt.eps) fail_flag=fail_flag+2 ! but computation is completed
! find true anomaly from dt
dt=el%t-el%coord(6)
alpha=gm*(el%coord(2)-1.d0)/el%coord(1) ! 2* energy
! CALL ta_from_t0(el%coord(1),ecc,dt,gm,alpha,v,dvdq,dvde,dvdt0)
! cot_com%coord(6)=v
cot_com%coo='COT'
! use brutal composition of two available changes
IF(PRESENT(del))THEN
   CALL coo_cha(el,'CAR',elcar,fail_flag1,delcar)
   CALL coo_cha(elcar,'COT',elcot, fail_flag2,dcarel)
   del=MATMUL(dcarel,delcar)
ELSE
   CALL coo_cha(el,'CAR',elcar,fail_flag1)
   CALL coo_cha(elcar,'COT',elcot, fail_flag2)
ENDIF 
cot_com%coord(6)=elcot%coord(6)
! cot_com%coord(6)=v
IF(PRESENT(del))THEN
  CALL eye(5,del(1:5,1:5))
!  del(6,1)=dvdq ! NOT READY
!  del(6,2)=dvde
  del(6,3:5)=0.d0
!  del(6,6)=dvdt0
  del(1:5,6)=0.d0
ENDIF
RETURN
! WARNING: what to do with fail_flag1 and fail_flag2 ???
END FUNCTION cot_com
! ======================================================
! TA_FROM_T0 true anomaly from time of prerihelion
! ======================================================
SUBROUTINE ta_from_t0(q,ecc,dt,mu,alpha,v,dvdq,dvde,dvdt0)
  USE ever_pitkin
  DOUBLE PRECISION,INTENT(IN) :: q,ecc ! perihelion distance, eccentricity
  DOUBLE PRECISION,INTENT(IN) :: dt ! dt=t-t0 difference between 
                       ! present time and time of passage at perihelion
  DOUBLE PRECISION,INTENT(IN) :: mu, alpha ! grav. constant, 2*energy
  DOUBLE PRECISION,INTENT(OUT) :: v ! true anomaly
  DOUBLE PRECISION,INTENT(OUT) :: dvdq,dvde,dvdt0 ! derivatives of true anomaly
! end interface                                   ! w.r.t. q,e,t0
  DOUBLE PRECISION :: sig0, rdot, psi, s0,s1,s2,s3, r
  DOUBLE PRECISION :: cosv,dvdr, dpsidt0, drdt0,drdpsi, ds2da,ds0da, drda
  DOUBLE PRECISION s0v,s1v,s2v,s3v,dd, psiv, eccv,qv, vv,cosvv
! *******************************
! compute psi at time t (psi0=0)
! *******************************
  sig0 = 0.d0
! solve universal Kepler's equation and compute Stumpff's functions
  CALL solve_kepuniv2(dt,q,sig0,mu,alpha,ecc,psi,s0,s1,s2,s3,DS2DA=ds2da,DS0DA=ds0da)
  dd=1.d-12
!  CALL s_funct(psi,alpha+dd,s0v,s1v,s2v,s3v)
!  CALL solve_kepuniv2(dt,q,sig0,mu,alpha+dd,ecc,psiv,s0v,s1v,s2v,s3v)
!  WRITE(*,*)' ds2da=',(s2v-s2)/dd, ds2da, s2
!  WRITE(*,*)' ds0da=',(s0v-s0)/dd, ds0da, s0
!  WRITE(*,*)' psi=', psi, psiv,psi-psiv
!  ds2da=(s2v-s2)/dd
!  ds0da=(s0v-s0)/dd
! *******************************
! compute r at time t 
! *******************************
  r = q*s0 + sig0*s1 + mu*s2
  drdpsi = sig0*s0 + (mu+alpha*q)*s1
  cosv =(-1.d0 + q/r*(1.d0+ecc))/ecc 
  v=acos(cosv)*sign(1.d0,drdpsi)

  qv=q+dd
  cosvv =(-1.d0 + qv/r*(1.d0+ecc))/ecc 
  vv=acos(cosvv)*sign(1.d0,drdpsi)
  dvdq=(vv-v)/dd
  WRITE(*,*)' dvdq direct inc=', dvdq
  eccv=ecc+dd
  cosvv =(-1.d0 + q/r*(1.d0+eccv))/eccv 
  vv=acos(cosvv)*sign(1.d0,drdpsi)
  dvde=(vv-v)/dd
  WRITE(*,*)' dvde direct inc=', dvde

! derivatives of v w.r. to q,e 
! ***********WARNING: WRONG***************
! missing dependence of r from q,e 
! dr/dq=s0+q*d(s0)/dq+ mu*d(s2)/dq 
! dr/d(alpha)=q*d(s0)/d(alpha)+ mu*d(s2)/d(alpha)
! d(sj)/dq=d(sj)/d(alpha) d(alpha)/dq ; d(sj)/de=d(sj)/d(alpha) d(alpha)/de
! d(alpha)/dq=-alpha/q ; d(alpha)/de= mu/q
  dvdq = -1.d0/sqrt(1.d0-cosv**2)*(1.d0+ecc)/(r*ecc) *sign(1.d0,drdpsi)
  WRITE(*,*)' dvdq direct=', dvdq
  dvdr =  1.d0/sqrt(1.d0-cosv**2)*q*(1.d0+ecc)/(r**2*ecc)*sign(1.d0,drdpsi)! OK
  drda= q*ds0da + mu*ds2da
!  drdq=s0+(q*ds0da+mu*ds2da)*(-alpha/q)
  dvdq=dvdq + dvdr*(drda*(-alpha/q)+s0)
  dvde = -1.d0/sqrt(1.d0-cosv**2)*(1.d0-q/r)/ecc**2*sign(1.d0,drdpsi)
  WRITE(*,*)' dvde direct=', dvde
  dvde= dvde + dvdr*drda*mu/q
! ************END WARNING*****************
! derivative of v w.r. to t0
  dpsidt0 = -1.d0/(q*s0 + sig0*s1 + mu*s2)
  drdt0 = drdpsi*dpsidt0
  dvdt0 = dvdr*drdt0
END SUBROUTINE ta_from_t0

! =======================================================
! COM_COT
! cometary elelements from cometary with true anomaly
! =======================================================
TYPE(orbit_elem) FUNCTION com_cot(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be cometary
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION ecc, eps, gm !, enne, a, dldn, dt
!DOUBLE PRECISION,DIMENSION(6) :: eq
!DOUBLE PRECISION r,x1,y1,sinu,cosu,u,delldu,dellde,dudv,dude,v2,ell,delldv
DOUBLE PRECISION delcar(6,6),dcarel(6,6)
TYPE(orbit_elem) elcar,elcom
INTEGER fail_flag1, fail_flag2
! =======================================================
IF(el%coo.ne.'COT')THEN
   WRITE(*,*)' com_cot: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
!gm=centerbody_mass(el%center)
eps=100*epsilon(1.d0) ! negligible number
!eq=el%coord
com_cot=el  !default
ecc=el%coord(2) ! eccentricity
!  test on eccentricity                                                 
IF(ecc.lt.eps) fail_flag= 1 ! but computation is completed
!   test on inclination                                 
IF(el%coord(3).lt.eps) fail_flag=fail_flag+2 ! but computation is completed
! use brutal composition of two available changes
IF(PRESENT(del))THEN
   CALL coo_cha(el,'CAR',elcar,fail_flag1,delcar)
   CALL coo_cha(elcar,'COM',elcom, fail_flag2,dcarel)
   del=MATMUL(dcarel,delcar)
ELSE
   CALL coo_cha(el,'CAR',elcar,fail_flag1)
   CALL coo_cha(elcar,'COM',elcom, fail_flag2)
ENDIF 
! WARNING: what to do with fail_flag1 and fail_flag2 ???
com_cot=elcom
END FUNCTION com_cot
! =======================================================
! KEP_COM
! keplerian elelements from cometary
! =======================================================
TYPE(orbit_elem) FUNCTION kep_com(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be cometary
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION ecc, eps, gm, enne, a, dldn
DOUBLE PRECISION,DIMENSION(6) :: eq
!
IF(el%coo.ne.'COM')THEN
   WRITE(*,*)' kep_com: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
gm=centerbody_mass(el%center)
eps=100*epsilon(1.d0) ! negligible number
eq=el%coord
kep_com=el  !default
ecc=eq(2)
if(ecc.ge.1.d0-eps)THEN
! this conversion is impossible
   fail_flag=5
! in output the elements are left as 'COM', del= identity
   IF(PRESENT(del)) CALL eye(el%ndim,del)
   RETURN
ENDIF
a=eq(1)/(1.d0-ecc)
kep_com%coord(1)=a
!  test on eccentricity                                                 
if(ecc.lt.eps) fail_flag= 1 ! but computation is completed
!   test on tangent of half inclination                                 
if(eq(3).lt.eps) fail_flag=fail_flag+2 ! but computation is completed
! time of pericenter from mean anomaly
enne=sqrt(gm/a**3)
kep_com%coord(6)=(el%t-eq(6))*enne 
kep_com%coo='KEP'
IF(PRESENT(del).and.fail_flag.eq.0)THEN
   del=0.d0 ! jacobian is not usable in the nearly singular case
   del(1,1)=  1.d0/(1.d0-ecc)
   del(1,2)=  +eq(1)/(1.d0-ecc)**2 
   del(2,2)=  1.d0
   del(3,3)=  1.d0
   del(4,4)=  1.d0
   del(5,5)=  1.d0
   dldn=      el%t-eq(6)
   del(6,1)=  -1.5d0*enne/eq(1)*dldn
   del(6,2)=  -1.5d0*enne/(1.d0-ecc)*dldn
   del(6,6)=  -enne
ENDIF

RETURN
END FUNCTION kep_com

! ====================================================================
! END ELEMENTS INTERNAL CONVERSIONS
! ======================================================================
! ATT submodule
! ====================================================================
! =======================================================
! ATT_CAR
! attributable elements from cartesian position and velocity, 
! planetocentric  
! attributable is:
!     %coord(1) = alpha
!     %coord(2) = delta
!     %coord(3) = d(alpha)/dt
!     %coord(4) = d(delta)/dt
!     %coord(5) = range
!     %coord(6) = range rate
! WARNING: att in equatorial, car in ecliptic coordinates 
! if the center is different from the Earth
! =======================================================
TYPE(orbit_elem) FUNCTION att_car(el,obscode,target_center,del)
  USE reference_systems
  TYPE(orbit_elem),INTENT(IN):: el
  INTEGER, INTENT(IN) :: obscode
  INTEGER, INTENT(in) ::target_center
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
  DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim),INTENT(OUT),OPTIONAL :: del 
! end interface
  DOUBLE PRECISION,DIMENSION(6) :: att ! attributable
  DOUBLE PRECISION, DIMENSION(3) :: dadx,dddx,d,dv !first derivatives, 
  DOUBLE PRECISION, DIMENSION(3,3) :: ddadx,ddddx ! second derivatives
  DOUBLE PRECISION :: dz,dis0,vsize,prscal,tobs,princ ! distances, functions
  DOUBLE PRECISION, DIMENSION(6) :: xtop,xobs
! aux. var for computation of second derivatives
  DOUBLE PRECISION den,x2,y2,z2,x4,y4 
  DOUBLE PRECISION x2y2,x2z2,y2z2, ddd(3)
! =======================================================               
  IF(el%coo.ne.'CAR')THEN
     WRITE(*,*)' att_car: wrong source coordinates ', el%coo
     STOP
  ENDIF
  if(target_center.ne.0.and.target_center.ne.3)then
     write(*,*)'att_car: Center different from Earth or Sun'
     stop
  end if
! compute observer position, taking into account light time
  CALL iter_obs(el%t,el%coord,obscode,xtop,tobs,dis0,target_center)
  d=xtop(1:3)
  dv=xtop(4:6)
! Computation of observation: right ascension (radians)                 
  dz=d(1)**2+d(2)**2 
  if (dz.le.100*epsilon(1.d0)) then
! remove singularity at poles 
     att(1)=0.d0 
  else 
     att(1)=atan2(d(2),d(1)) 
  endif
! Computation of observation: declination (radians)                     
  att(2)=asin(d(3)/dis0) 
! Computation of first derivatives of $\alpha$ and $\delta$ w.r. to posi
! (if required): we derived eq. (2.20)                                  
  dadx(1)=-d(2)/dz 
  dadx(2)=d(1)/dz 
  dadx(3)=0.d0 
  dddx(1)=-d(3)*(d(1)/(sqrt(dz)*dis0**2)) 
  dddx(2)=-d(3)*(d(2)/(sqrt(dz)*dis0**2)) 
  dddx(3)=sqrt(dz)/dis0**2 
! Apparent motion:                                                      
  att(3)=prscal(dadx,dv) 
  att(4)=prscal(dddx,dv)
! range, range rate
  att(5)=dis0
  att(6)=prscal(d,dv)/dis0
! assign attributable type elements
  att_car=el !defaults
  att_car%center=3
  att_car%coord=att
  att_car%coo='ATT'
! aberration correction NO!!!
  att_car%obscode=obscode
! h magnitude is already the absolute one
! derivatives, if required
  IF(PRESENT(del))THEN
     if(target_center.eq.0)then
        del=0.d0
        del(1,1:3)= MATMUL(roteqec,dadx) ! = MATMUL(dadx,roteceq)
        del(2,1:3)=MATMUL(roteqec,dddx)
        del(3,4:6)=del(1,1:3) ! commutation of d/dt
        del(4,4:6)=del(2,1:3) ! commutation of d/dt
        del(5,1:3)=MATMUL(roteqec,d*(1.d0/dis0))
        del(6,1:3)=MATMUL(roteqec,d*(-att(6)/dis0**2)+dv*(1.d0/dis0))
        del(6,4:6)=del(5,1:3) ! commutation of d/dt
! ===================================================================== 
! partials of adot,ddot with respect to positions require the           
! second derivatives of alpha, delta with respect to the                
! equatorial reference system                                           
! ===================================================================== 
! Computation of second derivatives of $\alpha$ w.r. to positions       
        ddadx(1,1)=2.d0*d(1)*d(2)/dz**2 
        ddadx(1,2)=(d(2)**2-d(1)**2)/dz**2 
        ddadx(2,1)=ddadx(1,2) 
        ddadx(2,2)=-ddadx(1,1) 
        ddadx(3,1)=0.d0 
        ddadx(3,2)=0.d0 
        ddadx(3,3)=0.d0 
        ddadx(2,3)=0.d0 
        ddadx(1,3)=0.d0 
! chain rule for derivatives of adot                                    
        ddd=MATMUL(ddadx,dv) 
! rotation to the equatorial reference system   
        del(3,1:3)=MATMUL(roteqec,ddd)
! =======================================================               
! Computation of second derivatives of $\delta$ w.r. to positions       
        den=1.d0/(dis0**4*dz*sqrt(dz)) 
        x2=d(1)**2
        y2=d(2)**2 
        z2=d(3)**2
        x4=x2*x2 
        y4=y2*y2                                    
        x2y2=x2*y2 
        x2z2=x2*z2 
        y2z2=y2*z2 
        ddddx(1,1)=d(3)*(2.d0*x4+x2y2-y2z2-y4)*den 
        ddddx(2,2)=d(3)*(2.d0*y4+x2y2-x2z2-x4)*den 
        ddddx(1,2)=d(1)*d(2)*d(3)*(z2+3.d0*x2+3.d0*y2)*den 
        ddddx(2,1)=ddddx(1,2) 
        ddddx(3,3)=-2.d0*d(3)*dz**2*den 
        ddddx(1,3)=d(1)*dz*(z2-x2-y2)*den 
        ddddx(3,1)=ddddx(1,3) 
        ddddx(2,3)=d(2)*dz*(z2-x2-y2)*den 
        ddddx(3,2)=ddddx(2,3) 
! chain rule for derivatives of ddot                                    
        ddd=MATMUL(ddddx,dv)
! rotation to the equatorial reference system 
        del(4,1:3)=MATMUL(roteqec,ddd)
! the same for the earth, but without rotation
     elseif(target_center.eq.3)then
        del=0.d0
        del(1,1:3)=dadx
        del(2,1:3)=dddx
        del(3,4:6)=del(1,1:3)
        del(4,4:6)=del(2,1:3)
        del(5,1:3)=d*(1.d0/dis0)
        del(6,1:3)=d*(-att(6)/dis0**2)+dv*(1.d0/dis0)
        del(6,4:6)=del(5,1:3) ! commutation of d/dt
        ddadx(1,1)=2.d0*d(1)*d(2)/dz**2 
        ddadx(1,2)=(d(2)**2-d(1)**2)/dz**2 
        ddadx(2,1)=ddadx(1,2) 
        ddadx(2,2)=-ddadx(1,1) 
        ddadx(3,1)=0.d0 
        ddadx(3,2)=0.d0 
        ddadx(3,3)=0.d0 
        ddadx(2,3)=0.d0 
        ddadx(1,3)=0.d0 
        ddd=MATMUL(ddadx,dv) 
        del(3,1:3)=ddd
        den=1.d0/(dis0**4*dz*sqrt(dz)) 
        x2=d(1)**2
        y2=d(2)**2 
        z2=d(3)**2
        x4=x2*x2 
        y4=y2*y2                                    
        x2y2=x2*y2 
        x2z2=x2*z2 
        y2z2=y2*z2 
        ddddx(1,1)=d(3)*(2.d0*x4+x2y2-y2z2-y4)*den 
        ddddx(2,2)=d(3)*(2.d0*y4+x2y2-x2z2-x4)*den 
        ddddx(1,2)=d(1)*d(2)*d(3)*(z2+3.d0*x2+3.d0*y2)*den 
        ddddx(2,1)=ddddx(1,2) 
        ddddx(3,3)=-2.d0*d(3)*dz**2*den 
        ddddx(1,3)=d(1)*dz*(z2-x2-y2)*den 
        ddddx(3,1)=ddddx(1,3) 
        ddddx(2,3)=d(2)*dz*(z2-x2-y2)*den 
        ddddx(3,2)=ddddx(2,3) 
        ddd=MATMUL(ddddx,dv)
        del(4,1:3)=ddd
     else
        write(*,*)'att_car: Center different from Earth or Sun'
        stop
     end if
  ENDIF
END FUNCTION att_car

! =======================================================
! ITER-OBS
! =======================================================
SUBROUTINE iter_obs(t,x,obscode,xtop,tobs,r,target_center)
USE reference_systems, ONLY: pvobs
! x is heliocentric ecliptic (tc=0), geocentric equatorial (tc=3)
DOUBLE PRECISION, INTENT(IN) :: t, x(6)
INTEGER, INTENT(IN) :: obscode, target_center
DOUBLE PRECISION, DIMENSION(6), INTENT(OUT):: xtop ! topocentric
! pos/vel with respect to observer obscode at appropriate time
! also ecliptic heliocentric position of observer
DOUBLE PRECISION, INTENT(OUT) :: tobs, r !obs time, range
DOUBLE PRECISION, DIMENSION(6) :: xpla
DOUBLE PRECISION vsize,tc, xobs(6)
INTEGER, PARAMETER :: jmax=5
DOUBLE PRECISION, PARAMETER :: tcont=1.d-8
INTEGER j
! change center: get ecliptic coordinates of the observer
tc=t
DO j=1,jmax
   CALL pvobs(tc,obscode,xobs(1:3),xobs(4:6)) ! xobs = geocentric ecliptic
   IF(target_center.eq.0)THEN 
      CALL earcar(tc,xpla,1) ! xpla=heliocentric ecliptic of geocenter
      xobs=xpla+xobs ! now xobs= heliocentric ecliptic of topos
      xtop=x-xobs   
      xtop(1:3)=MATMUL(roteceq,xtop(1:3)) ! xtop=topocentric equatorial
      xtop(4:6)=MATMUL(roteceq,xtop(4:6))
   ELSEIF(target_center.eq.3)THEN
      xobs(1:3)=MATMUL(roteceq,xobs(1:3)) ! xobs=geocentric equatorial
      xobs(4:6)=MATMUL(roteceq,xobs(4:6))
      xtop=x-xobs ! topocentric equatorial
   ELSE
      STOP
   ENDIF
! ecliptic topocentric vector of the asteroid
   xtop=x-xobs ! xtop=topocentric ecliptic of asteroid
   r=vsize(xtop)
   tobs=t +r/vlight !aberration correction
   IF(abs(tc-tobs).lt.tcont)RETURN
   tc=tobs
ENDDO
WRITE(*,*) 'iter_obs: slow convergence ', j-1, tobs, t, r
END SUBROUTINE iter_obs

! =======================================================
! NONITER_OBS
! =======================================================
SUBROUTINE noniter_obs(tobs,obscode,xobs, target_center)
  USE reference_systems
  DOUBLE PRECISION, INTENT(IN) :: tobs
  INTEGER, INTENT(IN) :: obscode, target_center
  DOUBLE PRECISION, DIMENSION(6), INTENT(OUT):: xobs ! topocentric ecliptic
! pos/vel with respect to observer obscode at appropriate time
  DOUBLE PRECISION, DIMENSION(6) :: xpla
! change center: get ecliptic coordinates of the observer
  CALL pvobs(tobs,obscode,xobs(1:3),xobs(4:6)) ! xobs = geocentric ecliptic
  IF(target_center.eq.0)THEN
     CALL earcar(tobs,xpla,1) ! xpla=heliocentric ecliptic of geocenter
     xobs=xpla+xobs ! now xobs= heliocentric ecliptic of topos
  ELSEIF(target_center.eq.3)THEN
     xobs(1:3)=matmul(roteceq,xobs(1:3))
     xobs(4:6)=matmul(roteceq,xobs(4:6))
  ELSE
     STOP
  ENDIF
END SUBROUTINE noniter_obs

! =======================================================
! CAR_ATT
! cartesian position and velocity, planetocentric, from 
! attributable elements 
! WARNING:  car in ecliptic coordinates, att in equatorial
! =======================================================
TYPE(orbit_elem) FUNCTION car_att(el,target_center,del)
  USE reference_systems
  USE fund_const, ONLY: vlight
  TYPE(orbit_elem),INTENT(IN) :: el
! New variable to know the center
  INTEGER, INTENT(in) ::target_center
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
  DOUBLE PRECISION, DIMENSION(el%ndim,el%ndim),INTENT(OUT),OPTIONAL :: del 
! end interface
  DOUBLE PRECISION, DIMENSION(3) :: rhat,ralphat,rdelhat,ralpalp, &
     & rdeldel,ralpdel
  DOUBLE PRECISION, DIMENSION(6) :: att,x
  DOUBLE PRECISION cosa,cosd,sina,sind,r,rdot,tobs
  DOUBLE PRECISION, DIMENSION(6) :: xobs
!
  IF(el%coo.ne.'ATT')THEN
     WRITE(*,*)' car_att: wrong source coordinates ', el%coo
     STOP
  ENDIF

  if(target_center.ne.3.and.target_center.ne.0)then
     write(*,*)'car_att: The center is different frome Earth and Sun'
     stop
  endif

  att=el%coord
! UNIT vectors of the r alpha delta ref. system
  cosa=cos(att(1))
  cosd=cos(att(2))
  sina=sin(att(1))
  sind=sin(att(2))
  r=att(5)
  rhat(1)=cosa*cosd
  rhat(2)=sina*cosd
  rhat(3)=sind
  ralphat(1)=-sina*cosd ! note: contains cosd and is not an unit vector
  ralphat(2)= cosa*cosd ! that is d(rhat)/d(alpha)=ralphat
  ralphat(3)=0.d0
  rdelhat(1)=-cosa*sind
  rdelhat(2)=-sina*sind
  rdelhat(3)=cosd
! second derivatives
  ralpalp(1)=-cosa*cosd
  ralpalp(2)=-sina*cosd
  ralpalp(3)=0.d0
  ralpdel(1)=sina*sind
  ralpdel(2)=-cosa*sind
  ralpdel(3)=0.d0
  rdeldel(1)=-cosa*cosd
  rdeldel(2)=-sina*cosd
  rdeldel(3)=-sind

  if(target_center.eq.0)then
! planetocentric position and velocity, equatorial, converted to ecliptic
     x(1:3)=MATMUL(roteqec,att(5)*rhat)
     x(4:6)=MATMUL(roteqec,att(6)*rhat+r*(att(3)*ralphat+att(4)*rdelhat))
  else
     x(1:3)=att(5)*rhat
     x(4:6)=att(6)*rhat+r*(att(3)*ralphat+att(4)*rdelhat)
  end if
! assign output
  car_att=el
! aberration correction
  tobs=el%t+r/vlight
  CALL noniter_obs(tobs,el%obscode,xobs,target_center)
! change center: heliocentric, ecliptic/ geocentric, equatorial pos. at epoch
  car_att%center=target_center
  car_att%coord=x+xobs ! 
  car_att%coo='CAR'
! Obscode does not change for the Earth case
  car_att%obscode=500
! h magnitude is already absolute
  IF(PRESENT(del))THEN
! for the earth case
     if(target_center.eq.3)then
        del=0.d0
        del(1:3,1)=ralphat*r
        del(1:3,2)=rdelhat*r
        del(1:3,5)=rhat
        del(4:6,1)=att(6)*ralphat+att(5)*(att(3)*ralpalp+att(4)*ralpdel)
        del(4:6,2)=att(6)*rdelhat+att(5)*(att(3)*ralpdel+att(4)*rdeldel)
        del(4:6,3)=del(1:3,1) ! commutation of d/dt
        del(4:6,4)=del(1:3,2) ! commutation of d/dt
        del(4:6,5)=att(3)*ralphat+att(4)*rdelhat
        del(4:6,6)=del(1:3,5) ! commutation of d/dt
     else
        del=0.d0
        del(1:3,1)=MATMUL(roteqec,ralphat*r)
        del(1:3,2)=MATMUL(roteqec,rdelhat*r)
        del(1:3,5)=MATMUL(roteqec,rhat)
        del(4:6,1)=MATMUL(roteqec,att(6)*ralphat+att(5)*(att(3)*ralpalp+att(4)*ralpdel))
        del(4:6,2)=MATMUL(roteqec,att(6)*rdelhat+att(5)*(att(3)*ralpdel+att(4)*rdeldel))
        del(4:6,3)=del(1:3,1) ! commutation of d/dt
        del(4:6,4)=del(1:3,2) ! commutation of d/dt
        del(4:6,5)=MATMUL(roteqec,att(3)*ralphat+att(4)*rdelhat)
        del(4:6,6)=del(1:3,5) ! commutation of d/dt
     end if
  ENDIF
END FUNCTION car_att
!=============================
! ATTELEMENTS 
! ATT elements from attributable and r,rdot
!=============================
SUBROUTINE attelements(att,r,rdot,elatt,unc)
  USE station_coordinates, ONLY: statcode
  USE attributable
  TYPE(attrib), INTENT(IN):: att ! 4-dim attributable
  DOUBLE PRECISION, INTENT(IN) :: r,rdot ! range, range rate
  TYPE(orbit_elem), INTENT(OUT):: elatt ! elements of type ATT
  TYPE(orb_uncert), INTENT(OUT), OPTIONAL :: unc ! covar/normal matrices
! end interface
  TYPE(orbit_elem) xcar
  INTEGER fail_flag,indp
  DOUBLE PRECISION rsun,phase,cosph, xx(3),xpla(6),vsize,prscal,appmag,ws(4)
!=============================
  elatt=undefined_orbit_elem
  elatt%coord(1:4)=att%angles
  elatt%coord(5)=r
  elatt%coord(6)=rdot
  elatt%coo='ATT'
  elatt%t=att%tdtobs-r/vlight
  CALL statcode(att%obscod, elatt%obscode)
! cartesian coordinates needed to get phase right
  CALL coo_cha(elatt,'CAR',xcar,fail_flag)
  xx=xcar%coord(1:3)
  rsun=vsize(xx)
  CALL earcar(att%tdtobs,xpla,1)
  cosph=prscal(xx,xx-xpla(1:3))/(rsun*r)
  IF(cosph.gt.1.d0)cosph=1.d0
  IF(cosph.lt.-1.d0)cosph=-1.d0
  phase=acos(cosph)
  IF(att%apm.gt.0.d0)THEN
     elatt%mag_set=.true.
     elatt%h_mag=att%apm-appmag(0.d0,elatt%g_mag,rsun,r,phase)
  ELSE
     elatt%mag_set=.false.
  ENDIF
  IF(PRESENT(unc))THEN
     unc=undefined_orb_uncert
     unc%g=0.d0
     unc%c=0.d0
     unc%g(1:4,1:4)=att%g
     unc%succ=.true.
     CALL tchinv(att%g,4,unc%c(1:4,1:4),ws,indp)
  ENDIF
END SUBROUTINE attelements

! ================================================
! ATT_PRELIM
! attributable preliminary orbit
! ================================================
SUBROUTINE att_prelim2(name0,attr,elk,uncatt,fail)
  USE attributable
  CHARACTER*(*), INTENT(IN) :: name0 ! asteroid designation
  TYPE(attrib), INTENT(IN) :: attr  ! attributable as computed 
  TYPE(orbit_elem), INTENT(OUT) :: elk ! new elements
  TYPE(orb_uncert), INTENT(OUT) :: uncatt ! its covariance
  LOGICAL, INTENT(OUT) :: fail
! END INTERFACE
  INTEGER iunmat,nroots,nrootsd
  DOUBLE PRECISION ::  c(0:5)! polynomial coefficients
  DOUBLE PRECISION, DIMENSION(3) :: roots ! real roots
  DOUBLE PRECISION, DIMENSION(8) :: rootsd ! roots of the derivative
  DOUBLE PRECISION :: a_orb,E_bound ! bound for the energy
! matrix of the r eps theta ref. system, topocentric equatorial
  DOUBLE PRECISION :: refs(3,3)
 ! observer pos/vel, asteroid pos/vel, equatorial heliocentric
  DOUBLE PRECISION :: xo(3), vo(3)
  DOUBLE PRECISION rr, rdot,ecc

  LOGICAL bizarre
  a_orb = 100.d0 ! comet boundary
!   a_orb = 2.5d0 ! main belt
  fail=.true.
  iunmat=0
  CALL admis_reg(name0,iunmat,attr%tdtobs,attr%angles,attr%obscod,  &
   &  nroots,roots,c,refs,xo,vo,a_orb,E_bound,nrootsd,rootsd)
! handle rare cases
  IF(nroots.eq.3)THEN
     write(iun_log,*)'att_prelim: ',name0,' two c.c.',roots(1:nroots)
     rr=(roots(3)+roots(2))/2.d0  ! privilege given to farther component
  ELSEIF(nroots.eq.2)THEN
     WRITE(iun_log,*)'att_prelim: ',name0, ' rare, nroots= ',nroots,roots(1:nroots)
     rr=roots(1)*0.8d0 ! somewhat inside the admissible region
  ELSEIF(nroots.eq.1)THEN
     rr=roots(1)*0.8d0 ! somewhat inside the admissible region
  ENDIF
  rdot=-c(1)/2.d0
  CALL attelements(attr,rr,rdot,elk,uncatt)
  
  fail=bizarre(elk,ecc)
  IF(verb_prelim.gt.10)THEN
     WRITE(iun_log,100)name0,rr,ecc
100  FORMAT(A,1X,'rho=',f8.3,1X,'ecc=',f8.4)
  ENDIF
END SUBROUTINE att_prelim2


! ====================================================================
! END ATT submodule
! ====================================================================
! ====================================================================
! OPIK submodule
! ====================================================================
! TPCAR_MTPCAR from perihelion position and velocity, obtain position on the 
! asymptote with rotation by -gamma/2 and rescaling by b/d of position,
! d/b of velocity
! ====================================================================
TYPE(orbit_elem) FUNCTION tpcar_mtpcar(el,fail_flag,bd,del)
  TYPE(orbit_elem), INTENT(IN) :: el ! must be cartesian
  INTEGER, INTENT(OUT) :: fail_flag 
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: bd ! b/d rescaling factor
  DOUBLE PRECISION, DIMENSION(6,6), OPTIONAL :: del ! jacobian matrix
! =============END INTERFACE==========================================
  DOUBLE PRECISION gm, v, d, eps, sig0, vsize, prscal
  DOUBLE PRECISION bsd, bsd2, sinhalfgam,coshalfgam
  DOUBLE PRECISION dbsddd, dbsddv, dgdd, dgdv ! derivatives of b/d, sinhalfgamma w.r. to d,v
  DOUBLE PRECISION, DIMENSION(3) :: x,y, ang, xtp,ytp, xmtp,xhat,yhat, xp,yp, xphat,yphat
  DOUBLE PRECISION, DIMENSION(3,3) :: rottp,drotpdg2 ! R, dr/d(sin(gam/2)
  DOUBLE PRECISION, DIMENSION(6,6) :: dlbd, dxpdx ! derivatives of rescaling,
                          ! of rotation by -halfgamma 
  DOUBLE PRECISION, DIMENSION(3,3,3) :: dtaudJ ! d(rotation)/d(ang)
  DOUBLE PRECISION, DIMENSION(3,3,3) :: ericci ! Ricci tensor
  DOUBLE PRECISION, DIMENSION(3,3)   :: dbJdX
  DOUBLE PRECISION, DIMENSION(3,3)   :: dbJdY
  DOUBLE PRECISION                   :: normJ
  INTEGER i,j,k,h ! =1,3
  IF(el%coo.ne.'CAR')THEN
     WRITE(*,*)' tpcar_mtpcar: wrong coords; center ',el%coo,' ', el%center
     STOP
  ENDIF
  fail_flag=0
  gm=centerbody_mass(el%center)
  tpcar_mtpcar=el ! copy defaults
  x=el%coord(1:3) ! position vector 
  y=el%coord(4:6) ! velocity vector
!  radius and velocity                                          
  v=vsize(y)      ! velocity
  xmtp=x-prscal(x,y)*y/v**2
  d=vsize(xmtp)      ! distance from center
  eps=2.d-7
! check on orthogonality
  sig0=prscal(x,y)/(v*d)
  IF(abs(sig0).gt.eps)THEN
     WRITE(*,199) sig0,eps
 199 FORMAT(' tpcar_mtpcar: pos, vel not orthogonal ',1P,2d10.3)
  ENDIF
  call prvec(x,y,ang)           !  angular momentum 
  bsd2=(d*v**2)/(d*v**2-2*gm)
  IF(bsd2.gt.0.d0)THEN
     bsd=sqrt(bsd2)
  ELSE
     fail_flag=6 ! means the orbit is not hyperbolic
     RETURN
  ENDIF
  sinhalfgam=gm/(v**2*d-gm)
  coshalfgam=sqrt(1.d0-sinhalfgam**2)
! halfgamma=asin(sinhalfgam)
  normJ=vsize(ang)
  ang=ang/normJ
! coordinates rotated by -halfgamma around ang
  CALL rotmhalfgam(sinhalfgam,ang,rottp,drotpdg2,dtaudJ)
  xp=MATMUL(rottp,x)
  yp=MATMUL(rottp,y)
!  xp=x
!  yp=y
! rescaling to asymptote
  xtp=xp*bsd
  ytp=yp*(1.d0/bsd)  
  tpcar_mtpcar%coord(1:3)=xtp
  tpcar_mtpcar%coord(4:6)=ytp
!  tpcar_mtpcar%coord(1:3)=xp
!  tpcar_mtpcar%coord(4:6)=yp
  IF(PRESENT(bd))THEN
    bd=bsd
  ENDIF
  tpcar_mtpcar%coo='TPC'
! partial derivatives
  IF(PRESENT(del))THEN
! derivatives d(xp,yp)/d(x,y) 
! neglecting derivatives with respect to ang
     xhat=x*(1.d0/d)
     yhat=y*(1.d0/v)
     dgdd=-v**2/gm*sinhalfgam**2
     dgdv=-2.d0*v*d/gm*sinhalfgam**2
     dxpdx=0.d0
     dxpdx(1:3,1:3)=rottp
     dxpdx(4:6,4:6)=rottp
     DO i=1,3
       DO j=1,3 
         DO k=1,3
             dxpdx(i,j)=dxpdx(i,j)+drotpdg2(i,k)*x(k)*dgdd*xhat(j)
             dxpdx(i+3,j+3)=dxpdx(i+3,j+3)+drotpdg2(i,k)*y(k)*dgdv*yhat(j)
             dxpdx(i,j+3)=dxpdx(i,j+3)+drotpdg2(i,k)*x(k)*dgdv*yhat(j)
             dxpdx(i+3,j)=dxpdx(i+3,j)+drotpdg2(i,k)*y(k)*dgdd*xhat(j)
         ENDDO
       ENDDO
     ENDDO
! derivatives of rottp w.r. to ang
     ericci(1:3,1:3,1:3)=0.d0
     ericci(1,2,3)=1.d0
     ericci(2,3,1)=1.d0
     ericci(3,1,2)=1.d0
     ericci(1,3,2)=-1.d0
     ericci(3,2,1)=-1.d0
     ericci(2,1,3)=-1.d0
     DO j=1,3
        DO h=1,3
         dbJdX(h,j)=1.d0/normJ*(DOT_PRODUCT(ericci(j,1:3,h),y)-v/d*x(j)*ang(h))
         dbJdY(h,j)=1.d0/normJ*(DOT_PRODUCT(ericci(1:3,j,h),x)-d/v*y(j)*ang(h))
      ENDDO
     ENDDO
! add term depending upon d(rottp0/d(ang)
      DO i=1,3
        DO j=1,3
          DO k=1,3
             dxpdx(i,j)=dxpdx(i,j)+DOT_PRODUCT(dtaudJ(i,k,1:3),dbJdX(1:3,j))*x(k)
             dxpdx(i,j+3)=dxpdx(i,j+3)+DOT_PRODUCT(dtaudJ(i,k,1:3),dbJdY(1:3,j))*x(k)
             dxpdx(i+3,j)=dxpdx(i+3,j)+DOT_PRODUCT(dtaudJ(i,k,1:3),dbJdX(1:3,j))*y(k)
             dxpdx(i+3,j+3)=dxpdx(i+3,j+3)+DOT_PRODUCT(dtaudJ(i,k,1:3),dbJdY(1:3,j))*y(k)
          ENDDO
        ENDDO
     ENDDO
! derivatives d(txp,typ)/d(xp,yp) 
     xphat=xp*(1.d0/d)
     yphat=yp*(1.d0/v)
     dbsddd=1.d0/(2.d0*bsd)*(-2.d0*gm*v**2)/(v**2*d-2.d0*gm)**2
     dbsddv=1.d0/(2.d0*bsd)*(-4.d0*gm*v*d)/(v**2*d-2.d0*gm)**2
     dlbd=0.d0
     DO i=1,3
       dlbd(i,i)=bsd
       dlbd(3+i,3+i)=1.d0/bsd
       DO j=1,3
         dlbd(i,j)=dlbd(i,j)+dbsddd*xphat(j)*xp(i)
         dlbd(i,j+3)=dbsddv*yphat(j)*xp(i)
         dlbd(i+3,j)=(-1.d0/bsd**2)*dbsddd*xphat(j)*yp(i)
         dlbd(i+3,j+3)=dlbd(i+3,j+3)+(-1.d0/bsd**2)*dbsddv*yphat(j)*yp(i)
       ENDDO
     ENDDO
     del=MATMUL(dlbd,dxpdx)
!     del=dlbd
!     del=dxpdx
  ENDIF
CONTAINS
SUBROUTINE rotmhalfgam(sgmez,bJ,tau,dtaudsgm,dtaudJ)
  DOUBLE PRECISION,INTENT(IN) :: sgmez !sinus(gamma/2)
  DOUBLE PRECISION,DIMENSION(3),INTENT(IN) :: bJ !ang. momentum unit vector
! rotation of -gamma/2 around bJ and its derivative w.r.t. sinus(gamma/2)
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(OUT) :: tau,dtaudsgm 
! derivative of tau w.r.t. the components of bJ
! dtaudJ(*,*,k) is the derivative w.r.t. the k-th component  
  DOUBLE PRECISION,DIMENSION(3,3,3),INTENT(OUT) :: dtaudJ 
! ------------- end interface -------------------------------------------
  DOUBLE PRECISION :: normJ,HRoy
  DOUBLE PRECISION :: cI,sI,cO,sO, Incl, Omeg, gmez
  DOUBLE PRECISION :: cgmez,tgmez ! cosinus and tangent of gamma/2
! rotation matrices
  DOUBLE PRECISION, DIMENSION(3,3) :: ROmeg
  DOUBLE PRECISION, DIMENSION(3,3) :: RmOmeg
  DOUBLE PRECISION, DIMENSION(3,3) :: RIncl
  DOUBLE PRECISION, DIMENSION(3,3) :: RmIncl
  DOUBLE PRECISION, DIMENSION(3,3) :: Rmhgam
! for the derivatives
  DOUBLE PRECISION,DIMENSION(3) :: dHRoydJ
  DOUBLE PRECISION, DIMENSION(3,3) :: delta
  DOUBLE PRECISION, DIMENSION(3,3) :: dRmhgam
  DOUBLE PRECISION, DIMENSION(3,3) :: RORI,RmIRmO,tauauxr,tauauxl !auxiliary
  DOUBLE PRECISION, DIMENSION(3,3,3) :: dROmegdJ,dRmOmegdJ
  DOUBLE PRECISION, DIMENSION(3,3,3) :: dRIncldJ,dRmIncldJ
  DOUBLE PRECISION, DIMENSION(3,3,3) :: tensor1,tensor2,tensor3,tensor4
  DOUBLE PRECISION :: vsize !functions
  INTEGER :: h,k ! loop indexes
! ---------------------------------------------
  normJ=vsize(bJ)
  IF(abs(normJ-1.d0).gt.100*epsilon(1.d0))THEN
     WRITE(*,*) 'rotmhalfgam: angular momentum not normalized! normJ=',normJ
  ENDIF
  HRoy=sqrt(bJ(1)**2+bJ(2)**2)
  IF(HRoy/normJ.lt.epsilon(1.d0)*100)THEN
     cI=1.d0
     sI=0.d0
     cO=1.d0
     sO=0.d0
  ELSE
     cI = bJ(3) ! bJ(3)/normJ (but J is unitary)
     sI = HRoy  ! HRoy/normJ
     cO = -bJ(2)/HRoy ! singular for zero I
     sO = bJ(1)/HRoy ! singular for zero I
  ENDIF
  cgmez = sqrt(1.d0-sgmez**2)
  tgmez = sgmez/cgmez

  Incl=atan2(sI,cI)
  Omeg=atan2(sO,cO)
  gmez=atan2(sgmez,cgmez)

! *******************************************************
! rotation matrix tau = ROmeg*RIncl*Rmhgam*RmIncl*RmOmeg
! *******************************************************
  CALL rotmt(-Omeg,ROmeg,3)
  CALL rotmt(Omeg,RmOmeg,3)
  CALL rotmt(-Incl,RIncl,1)
  CALL rotmt(Incl,RmIncl,1)
  CALL rotmt(gmez,Rmhgam,3)
  RORI=matmul(ROmeg,RIncl)
  RmIRmO=matmul(RmIncl,RmOmeg)
  tauauxr = matmul(Rmhgam,RmIRmO)
  tauauxl = matmul(RORI,Rmhgam) ! will be used later
  tau = matmul(RORI,tauauxr)

! *********************************************
! derivative of Rmhgam w.r.t. sgmez 
! dtaudsgm = ROmeg*RIncl*dRmhgam*RmIncl*RmOmeg
! *********************************************
  dRmhgam(1:3,1:3)=0.d0! initialization
  dRmhgam(1,1)=-tgmez
  dRmhgam(1,2)=1.d0
  dRmhgam(2,1)=-1.d0
  dRmhgam(2,2)=-tgmez
  dtaudsgm = matmul(dRmhgam, RmIRmO) ! temporary step
  dtaudsgm = matmul(RORI, dtaudsgm)

! *************************************************
! derivative of Rmhgam w.r.t. the components of bJ
! *************************************************
!  IF(HRoy/normJ.lt.epsilon(1.d0)*100)THEN
     dHRoydJ(1) = sO
     dHRoydJ(2) = -cO
!  ELSE
!     dHRoydJ(1) = bJ(1)/HRoy ! singular for zero I
!     dHRoydJ(2) = bJ(2)/HRoy ! singular for zero I
!  ENDIF
  dHRoydJ(3) = 0.d0

  CALL eye(3,delta) ! Kronecker's delta

  dROmegdJ(1:3,1:3,1:3)=0.d0 ! initialization
  IF(HRoy/normJ.ge.epsilon(1.d0)*100)THEN
     dROmegdJ(1,1,1:3)= -(1.d0/HRoy)*delta(2,1:3) + (bJ(2)/HRoy**2)*dHRoydJ(1:3) ! singular for zero I
     dROmegdJ(1,2,1:3)= -(1.d0/HRoy)*delta(1,1:3) + (bJ(1)/HRoy**2)*dHRoydJ(1:3) ! singular for zero I
  ELSE
!     dROmegdJ(1,1,1:3)= -(1.d0/HRoy)*delta(2,1:3) +  (-cO/HRoy)*dHRoydJ(1:3) ! singular for zero I
!     dROmegdJ(1,2,1:3)= -(1.d0/HRoy)*delta(1,1:3) +  (sO/HRoy)*dHRoydJ(1:3) ! singular for zero I
  ENDIF
  dROmegdJ(2,1,1:3)= -dROmegdJ(1,2,1:3)
  dROmegdJ(2,2,1:3)=  dROmegdJ(1,1,1:3)
  dRmOmegdJ(1:3,1:3,1:3)= dROmegdJ(1:3,1:3,1:3)
  dRmOmegdJ(1,2,1:3)= -dROmegdJ(1,2,1:3)
  dRmOmegdJ(2,1,1:3)= -dRmOmegdJ(1,2,1:3)
  dRIncldJ(1:3,1:3,1:3)=0.d0 ! initialization
  dRIncldJ(2,2,1:3)=  delta(3,1:3)
  dRIncldJ(2,3,1:3)= -dHRoydJ(1:3)
  dRIncldJ(3,2,1:3)= -dRIncldJ(2,3,1:3)
  dRIncldJ(3,3,1:3)=  dRIncldJ(2,2,1:3)

  dRmIncldJ(1:3,1:3,1:3)= dRIncldJ(1:3,1:3,1:3)
  dRmIncldJ(2,3,1:3)= -dRIncldJ(2,3,1:3)
  dRmIncldJ(3,2,1:3)= -dRmIncldJ(2,3,1:3)

! **********************************************
! tensor1 = dROmegdJ*RIncl*Rmhgam*RmIncl*RmOmeg
! **********************************************
  DO k=1,3
     tensor1(1:3,1:3,k) = matmul(RIncl,tauauxr)
     tensor1(1:3,1:3,k) = matmul(dROmegdJ(1:3,1:3,k),tensor1(1:3,1:3,k))
  ENDDO
! **********************************************
! tensor2 = ROmeg*dRIncldJ*Rmhgam*RmIncl*RmOmeg
! **********************************************
  DO k=1,3
     tensor2(1:3,1:3,k) = matmul(dRIncldJ(1:3,1:3,k),tauauxr)
     tensor2(1:3,1:3,k) = matmul(ROmeg,tensor2(1:3,1:3,k))
  ENDDO
! **********************************************
! tensor3 = ROmeg*RIncl*Rmhgam*dRmIncldJ*RmOmeg
! **********************************************
  DO k=1,3
     tensor3(1:3,1:3,k) = matmul(tauauxl,dRmIncldJ(1:3,1:3,k))
     tensor3(1:3,1:3,k) = matmul(tensor3(1:3,1:3,k),RmOmeg)
  ENDDO
! **********************************************
! tensor4 = ROmeg*RIncl*Rmhgam*RmIncl*dRmOmegdJ
! **********************************************
  DO k=1,3
     tensor4(1:3,1:3,k) = matmul(tauauxl,RmIncl)
     tensor4(1:3,1:3,k) = matmul(tensor4(1:3,1:3,k),dRmOmegdJ(1:3,1:3,k))
  ENDDO
  
  dtaudJ = tensor1 + tensor2 + tensor3 + tensor4

END SUBROUTINE rotmhalfgam

END FUNCTION tpcar_mtpcar
! ====================================================================
! OPIK_TPCAR from perihelion position and velocity, obtain Opik elements 
! ====================================================================
! opik= U, long, lat, csi, zeta, eta
! U= unperturbed v_infty
! long, lat  = antiradians (in ecliptical coordinates)
! csi, zeta = TP coordinates with axis zeta along projection of normal to 
!                   ecliptic
! eta = U(t-t0), t0 time of pericenter passage
! WARNING: in .clo file we need also the time dependent theta, phi 
! to be computed using the Earth velocity at t0
! ====================================================================
TYPE(orbit_elem) FUNCTION opik_tpcar(el,fail_flag,del,vt3)
  TYPE(orbit_elem), INTENT(IN) :: el ! must be cartesian
  INTEGER, INTENT(OUT) :: fail_flag 
  DOUBLE PRECISION, DIMENSION(6,6), OPTIONAL :: del ! jacobian matrix
  DOUBLE PRECISION, DIMENSION(3,3), OPTIONAL :: vt3 ! rotation of TP to plane eta=0
! =============END INTERFACE==========================================
  DOUBLE PRECISION opik(6), xytp(6), vsize, b, dz, prscal, cosa, dt
  DOUBLE PRECISION v1(3), v2(3), rv(3,3),xx(3),drv(3,3,3)!,rvt(3,3), rv1(3,3)
  INTEGER i,j
  IF(el%coo.ne.'TPC')THEN
     WRITE(*,*)' opik_tpcar: wrong coordinates ', el
     STOP
  ENDIF
  fail_flag=0
  opik_tpcar=el ! set default; note that epoch remains the one of pericenter
  xytp=el%coord
  opik(1)=vsize(xytp(4:6))
  v1=xytp(4:6)*(1.d0/opik(1))
  b=vsize(xytp(1:3))
  cosa=prscal(xytp(1:3),xytp(4:6))/(b*opik(1))
  dt=-b*cosa/opik(1)
  IF(abs(cosa).gt.2.d-7)THEN
!     xytp(1:3)=xytp(1:3)-b*cosa*v1
!     b=vsize(xytp(1:3))
     WRITE(*,199)cosa, b, opik(1)
 199 FORMAT('opik_tpcar: pos not orthogonal to vel ',1p,3D10.3)
!     fail_flag=7
  ENDIF
! polar coordinates                 
  dz=xytp(4)**2+xytp(5)**2 
  if (dz.le.100*epsilon(1.d0)) then
! remove singularity at poles 
     opik(2)=0.d0 
  else 
! ecliptic longitude of velocity (radians)
     opik(2)=atan2(xytp(5),xytp(4)) 
     if (opik(2).lt.0.d0) then 
        opik(2)=opik(2)+dpig 
     endif
  endif
! ecliptic latitude of velocity (radians)                     
  opik(3)=asin(xytp(6)/opik(1))
! target plane 
  CALL drvdytp(xytp(4:6),drv,rv)
  IF(PRESENT(vt3))THEN
     vt3=rv
  ENDIF
  v1=xytp(4:6)*(1.d0/opik(1))
!  v2=(/0.d0, 0.d0, 1.d0/) 
!  CALL mtp_ref(v1,v2,rvt,rv) 
  xx=MATMUL(rv,xytp(1:3))
  opik(4)=xx(2)
  opik(5)=xx(3)
  IF(abs(xx(1)).gt.2.d-7)THEN
     WRITE(*,*)' opik_tpcar: pos, vel not orthogonal ', xx(1), xx(1)/b
  ENDIF
! last element is the time of pericenter passage 
  opik(6)=opik(1)*dt
! resulting elements
  opik_tpcar%t=el%t
  opik_tpcar%coord=opik
  IF(PRESENT(del))THEN
! partial derivatives
     del(1:3,1:3)=0.d0
     del(1,4:6)=v1
     del(2,4)=-xytp(5)/dz 
     del(2,5)=xytp(4)/dz 
     del(2,6)=0.d0 
     del(3,4)=-xytp(6)*(xytp(4)/(sqrt(dz)*opik(1)**2)) 
     del(3,5)=-xytp(6)*(xytp(5)/(sqrt(dz)*opik(1)**2)) 
     del(3,6)=sqrt(dz)/opik(1)**2
     del(4:5,4:6)=0.d0 ! neglecting derivatives of rv
     del(4:5,1:3)=rv(2:3,1:3) 
     del(6,1:3)=-(1.d0/opik(1))*xytp(4:6)
     del(6,4:6)=-(1.d0/opik(1))*xytp(1:3)+b/opik(1)*cosa*v1
! contribution from derivatives of rv
     DO i=4,5
       DO j=4,5
         del(i,j)=del(i,j)+DOT_PRODUCT(drv(i-2,1:3,j-3),xytp(1:3))
       ENDDO
     ENDDO
  ENDIF
  CONTAINS

SUBROUTINE drvdytp(ytp,drv,rv)
DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: ytp
                  ! Velocity vector in TP reference frame
DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: rv
                  ! rotation matrix        
DOUBLE PRECISION, DIMENSION(3,3,3), INTENT(OUT) :: drv 
                  ! Derivatives of rotation matrix rv
DOUBLE PRECISION :: U,sqrx2y2,trhalf,d1,d2, x,y,z
                  ! Auxiliary quantities
! ==================================================================
 U=sqrt(ytp(1)**2+ytp(2)**2+ytp(3)**2)          ! Module of velocity
 x=ytp(1)
 y=ytp(2)
 z=ytp(3)
! Matrix rv
 rv(1,1)=x/U  
 rv(1,2)=y/U  
 rv(1,3)=z/U
 rv(3,1)=-x*z/(U*sqrt(x**2+y**2))  
 rv(3,2)=-y*z/(U*sqrt(x**2+y**2))
 rv(3,3)=(x**2+y**2)/(U*sqrt(x**2+y**2))
 rv(2,1)=y/sqrt(x**2+y**2) 
 rv(2,2)=-x/sqrt(x**2+y**2)
 rv(2,3)=0.d0

! ==================================================================

sqrx2y2=sqrt(ytp(1)**2+ytp(2)**2)
trhalf=sqrx2y2*(ytp(1)**2+ytp(2)**2)
d1=U**3*trhalf
d2=U**3*sqrx2y2
! ===================================================
! Derivatives of first row of matrix R_v
drv(1,1,1)=(ytp(2)**2+ytp(3)**2)/U**3
drv(1,1,2)=-ytp(1)*ytp(2)/U**3
drv(1,1,3)=-ytp(1)*ytp(3)/U**3

drv(1,2,1)=-ytp(1)*ytp(2)/U**3
drv(1,2,2)=(ytp(1)**2+ytp(3)**2)/U**3
drv(1,2,3)=-ytp(2)*ytp(3)/U**3

drv(1,3,1)=-ytp(1)*ytp(3)/U**3
drv(1,3,2)=-ytp(2)*ytp(3)/U**3
drv(1,3,3)=(ytp(1)**2+ytp(2)**2)/U**3
! ===================================================
! Derivatives of third row of matrix R_v
drv(3,1,1)=ytp(3)*((ytp(1)**2-U**2)*(ytp(1)**2+ytp(2)**2)+U**2*ytp(1)**2)/d1
drv(3,1,2)=ytp(1)*ytp(2)*ytp(3)*(ytp(1)**2+ytp(2)**2+U**2)/d1
drv(3,1,3)=ytp(1)*(ytp(3)**2-U**2)/d2

drv(3,2,1)=ytp(1)*ytp(2)*ytp(3)*(ytp(1)**2+ytp(2)**2+U**2)/d1
drv(3,2,2)=ytp(3)*((ytp(2)**2-U**2)*(ytp(1)**2+ytp(2)**2)+U**2*ytp(2)**2)/d1
drv(3,2,3)=ytp(2)*(ytp(3)**2-U**2)/d2

drv(3,3,1)=ytp(1)*ytp(3)**2/d2
drv(3,3,2)=ytp(2)*ytp(3)**2/d2
drv(3,3,3)=-ytp(3)*sqrx2y2/U**3
! ====================================================
! Derivatives of second row of matrix R_v
drv(2,1,1)=-ytp(1)*ytp(2)/trhalf
drv(2,1,2)=ytp(1)**2/trhalf
drv(2,1,3)=0.d0

drv(2,2,1)=-ytp(2)**2/trhalf
drv(2,2,2)=ytp(1)*ytp(2)/trhalf
drv(2,2,3)=0.d0

drv(2,3,1)=0.d0
drv(2,3,2)=0.d0
drv(2,3,3)=0.d0
! =====================================================
END SUBROUTINE drvdytp



END FUNCTION opik_tpcar
! ====================================================================
! END OPIK submodule
! ====================================================================
! ====================================================================
! I/O submodule
! ====================================================================
SUBROUTINE read_elems(el,name,eof,file,unit,covar)
  USE io_elems, ONLY: obscod
  USE station_coordinates, ONLY: statcode
  TYPE(orbit_elem), INTENT(OUT) :: el
  CHARACTER*(*),INTENT(OUT) :: name
  LOGICAL, INTENT(OUT) :: eof
  CHARACTER*(*), INTENT(IN) :: file ! input file
  INTEGER, INTENT(IN), OPTIONAL :: unit ! input unit, only if file already opened
  TYPE(orb_uncert), INTENT(OUT), OPTIONAL :: covar 
! end interface
  DOUBLE PRECISION elem(6),t0,h,g,cove(6,6),nore(6,6),mass !,cnv(6)
  CHARACTER*3 eltype
  CHARACTER*10 rsys,epoch 
  LOGICAL defcov,defnor,defcn
  INTEGER kr, iobs 

  IF(.not.PRESENT(unit))THEN
     CALL oporbf(file,0)
  ENDIF
  el=undefined_orbit_elem
  CALL rdorb(name,elem,eltype,t0,cove,defcov,nore,defnor,     &
     &                 h,g,mass,rsys,epoch,kr,eof)
  IF(eof)THEN
     IF(.not.PRESENT(unit)) CALL clorbf
     RETURN
  ENDIF
! unit conversions already done by rdorb
  IF(eltype.eq.'ATT')THEN
! input obscod
     CALL statcode(obscod,iobs)
  END IF
! what to do with rsys,epoch????
  IF(rsys.ne.'ECLM'.or.epoch.ne.'J2000')THEN
     WRITE(*,*)' read_elem: other reference systems not handled ', rsys,'  ', epoch
     STOP
  ENDIF
! copy into output elements
  el%coo=eltype
  el%coord=elem
  el%t=t0
  IF(h.lt.-10.d0)THEN
     el%mag_set=.false.
  ELSE
     el%mag_set=.true.
     el%h_mag=h
     el%g_mag=g
  ENDIF
  IF(el%coo.eq.'ATT')el%obscode=iobs

  IF(PRESENT(covar))THEN
     covar=undefined_orb_uncert
     IF(defcov) covar%g=cove
     IF(defnor) covar%c=nore
     CALL fixcnm(defcov,defnor,defcn,cove,nore)
     IF(defcn)THEN
        covar%succ=.true.
     ELSE
        covar%succ=.false.
     ENDIF
  ENDIF

  IF(.not.PRESENT(unit)) CALL clorbf
END SUBROUTINE read_elems

SUBROUTINE write_elems(el,name,form,file,unit,covar,incfit)
  USE name_rules
  USE io_elems, ONLY: obscod
  USE station_coordinates, ONLY: codestat
  TYPE(orbit_elem), INTENT(IN) :: el
  CHARACTER*(*),INTENT(IN) :: name
  CHARACTER*(*), INTENT(IN) :: form ! can be either '1L' or 'ML'
  CHARACTER*(*), INTENT(IN), OPTIONAL :: file ! output file 
                       ! (not used if PRESENT(unit))
  INTEGER, INTENT(IN), OPTIONAL :: unit ! input unit, only if already opened
! WARNING: if file already opened, it must also have the header already written
  TYPE(orb_uncert), INTENT(IN), OPTIONAL :: covar 
  INTEGER, INTENT(IN), OPTIONAL :: incfit ! result from incomplete fit;
    ! warn the users of this!!!!
! =========== end interface================
  INTEGER uniout,le
  CHARACTER*(idnamvir_len) namloc
  CHARACTER*10 rsys,epoch
  LOGICAL defcov,defnor
  DOUBLE PRECISION mass, cove(6,6), nore(6,6) 
! open output file if not opened already
  IF(.not.PRESENT(unit))THEN
     CALL filopn(uniout,file,'unknown')
! write file header
     rsys='ECLM'
     epoch='J2000'
     IF(form.eq.'1L')THEN
! what to do when el%coo=ATT ?????
        CALL codestat(el%obscode,obscod)
        CALL wro1lh(uniout,rsys,epoch,el%coo)
     ELSEIF(form.eq.'ML')THEN
        CALL wromlh(uniout,rsys,epoch)
     ENDIF
  ELSE
     uniout=unit
  ENDIF
  namloc=name
  CALL rmsp(namloc,le)
  IF(PRESENT(incfit))THEN
     IF(incfit.eq.5)THEN
        WRITE(uniout,200)namloc(1:le)
200     FORMAT('! CONSTRAINED SOLUTION for ',a)
     ELSEIF(incfit.eq.4)THEN
        WRITE(uniout,201)namloc(1:le)
201     FORMAT('! 4_PARAMETERS SOLUTION for ',a)
     ENDIF
  ENDIF
! COVARIANCE; HOWEVER NOT WRITTEN IF FORMAT= '1L'
  IF(PRESENT(covar))THEN
     IF(covar%succ)THEN
        defcov=.true.
        defnor=.true.
        cove=covar%g
        nore=covar%c
     ELSE
        defcov=.false.
        defnor=.false.
        cove=0.d0
        nore=0.d0
     ENDIF
  ELSE
     defcov=.false.
     defnor=.false.
     cove=0.d0
     nore=0.d0
  ENDIF
  mass=0.d0
  IF(form.eq.'1L')THEN
     CALL wro1lr(uniout,namloc(1:le),el%coord,el%coo,el%t,el%h_mag,el%g_mag)
  ELSEIF(form.eq.'ML')THEN
     CALL codestat(el%obscode,obscod)
     CALL wromlr(uniout,namloc(1:le),el%coord,el%coo,el%t,cove,defcov,      &
          &                  nore,defnor,el%h_mag,el%g_mag,mass)
  ENDIF

  IF(.not.PRESENT(unit)) CALL filclo(uniout,' ')

END SUBROUTINE write_elems

SUBROUTINE convertunc(unc1,dee,unc2)
! unc1 is input, unc2 is output, but they can be on the same memory location
  TYPE(orb_uncert), INTENT(IN) :: unc1
  TYPE(orb_uncert), INTENT(INOUT) :: unc2
  DOUBLE PRECISION, DIMENSION(6,6), INTENT(IN) :: dee
  LOGICAl error
! covariance matrix is propagated by similarity transformation     
  CALL convertcov(unc1%g,dee,unc2%g)
! normal matrix is propagated by similarity transformation  
  CALL norprs(unc1%c,dee,6,unc2%c,error)
  IF(error)THEN
     unc2%succ=.false.
     unc2%ndim=unc1%ndim
  ELSE
     unc2%succ=.true.
     unc2%ndim=unc1%ndim
  ENDIF
END SUBROUTINE convertunc

CHARACTER*3 FUNCTION cootyp(ityp)
  INTEGER, INTENT(IN) :: ityp
  if(ityp.eq.1)then 
     cootyp='KEP' 
  elseif(ityp.eq.2)then 
     cootyp='EQU' 
  elseif(ityp.eq.3)then 
     cootyp='CAR'
  elseif(ityp.eq.4)then 
     cootyp='COM'
  elseif(ityp.eq.5)then 
     cootyp='COT'
  elseif(ityp.eq.6)then 
     cootyp='ATT' 
  else
     WRITE(*,*)' cootyp: coord code not known ',ityp
     STOP
  endif
END FUNCTION cootyp

! Copyright (C) 1997-1998 by Orbfit Consortium                          
! Version: December 16, 1998 Steven Chesley chesley@dm.unipi.it         
! Version 3.1 September 2003 A.Milani                               
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         W R I K E P                           *    
!  *                                                               *    
!  *  Writes an orbital element record in an orbital element file  *    
!  *                     (multi-line format)                       *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INTPUT:   FILE      -  Output file name                               
!           NAME      -  Name of planet/asteroid/comet                  
!           ELEM      -  Orbital elements
!           COVAR     -  Uncertainty of orbital elements          
!                                                                       
! OUPUT:    MOID      -  MOID                                           
!           PHA       -  0 or 1 PHA flag                                
!                                                                       
! WARNING: the routine does write the header of the file: this          
!          must not be generated by calling subroutine wromlh           
!                                                                       
SUBROUTINE wrikep(file,name,elem,covar,moid,pha)
  CHARACTER*(*), INTENT(INOUT) :: name
  TYPE(orbit_elem), INTENT(IN) :: elem
  TYPE(orb_uncert), INTENT(IN) :: covar
  DOUBLE PRECISION, INTENT(OUT) ::  moid 
  CHARACTER*(*), INTENT(INOUT) :: file 
  LOGICAL, INTENT(OUT) :: pha 
  TYPE(orbit_elem) :: elk
  TYPE(orb_uncert) :: covk,cove
  INCLUDE 'parcmc.h90' 
  INTEGER lnn,lnf,i,j,k,unit, fail_flag
  DOUBLE PRECISION cnv(6),std(6)
  DOUBLE PRECISION enne, dkde(6,6), correl(6,6) 
                                                                        
! "Other useful data"                                                   
  DOUBLE PRECISION perihe,aphe,period,dnode,anode,vinfty
  INTEGER iconv 
                                                                        
! Open FIle    
  call rmsp(file,lnf) 
  call filopn(unit,file(1:lnf),'unknown') 
  CALL wromlh (unit,'ECLM','J2000') 
                                                                        
! Object Name     
  call rmsp(name,lnn) 
  WRITE(unit,100) name(1:lnn) 
100 FORMAT(A) 
! Get Keplerian Elements
  CALL coo_cha(elem,'KEP',elk,fail_flag,dkde) 
! WARNING: not working for cometary elements
! Write Orbital elements
  cnv(1:2)=1
  cnv(3:6)=degrad 
  WRITE(unit,201) comcha 
201 FORMAT(A,' Keplerian elements: a, e, i, long. node,',             &
     &         ' arg. peric., mean anomaly')                            
  WRITE(unit,101) (elk%coord(i)*cnv(i),i=1,6) 
101 FORMAT(' KEP ',2F12.6,4F12.3) 
! Epoch
  WRITE(unit,104) elk%t 
104 FORMAT(' MJD ',F12.4,' TDT') 
! Mass                                                                  
!      IF(mass.NE.0.d0) WRITE(unit,105) mass 
!  105 FORMAT(' MAS ',1P,E20.12) 
                                                                        
! Magnitudes                                                            
  IF(elk%mag_set) WRITE(unit,106) elk%h_mag,elk%g_mag 
106 FORMAT(' MAG ',2F7.3) 
                                                                        
! Other useful data:                                                    
! PERIHELION 1.2345                                                     
! APHELION 2.3456                                                       
! ANODE 1.01234                                                         
! DNODE 1.04321                                                         
! MOID 0.00001                                                          
! PERIOD 1.11111                                                        
! PHA F                                                                 
! VINFTY 14.5678 (in km/s)
  perihe=elk%coord(1)*(1.d0-elk%coord(2)) 
  aphe=elK%coord(1)*(1.d0+elk%coord(2)) 
  enne=sqrt(gms/elk%coord(1)**3)
  period=dpig/enne 
  CALL nomoid(elem%t,elem,moid,anode,dnode) 
  pha = moid.le.0.05d0.and.elk%h_mag.lt.22d0.and.elk%mag_set
  vinfty= v_infty0(elk)
  WRITE(unit,102)comcha,perihe,comcha,aphe,comcha,                  &
     &     anode,comcha,dnode,comcha,moid,comcha,period,comcha,pha, &
     &     comcha,vinfty      
  102 FORMAT(A1,' PERIHELION ',F9.4/                                    &
     &       A1,' APHELION ',F9.4/                                      &
     &       A1,' ANODE ',F10.5/                                        &
     &       A1,' DNODE ',F10.5/                                        &
     &       A1,' MOID ',F10.5/                                         &
     &       A1,' PERIOD ',F16.4/                                       &
     &       A1,' PHA ',L1/                                             &
     &       A1,' VINFTY ',F9.4)
! Get Keplerian Covariance  
  IF(covar%succ) THEN 
     cove=covar
     CALL convertunc(cove,dkde,covk)
! Get Correlation Matrix                                                
     DO i=1,6 
        std(i) = sqrt(covk%g(i,i)) 
     ENDDO
     DO i=1,6 
        DO j=1,6 
           correl(i,j) = covk%g(i,j)/std(i)/std(j) 
        ENDDO
     ENDDO                                                                   
! Covariance matrix                                                     
     WRITE(unit,107) comcha,(std(i)*cnv(i),i=1,6) 
     WRITE(unit,108) ((covk%g(i,k)*cnv(i)*cnv(k),k=i,6),i=1,6) 
107  FORMAT(A1,' RMS ',1P,6E12.3) 
108  FORMAT(' COV ',1P,3E23.8) 
! Correlation matrix                                                    
     WRITE(unit,109)  ((correl(i,k),k=i,6),i=1,6) 
109  FORMAT(' COR ',3F23.8) 
  ENDIF
  
  call filclo(unit,' ') 
                                                                        
END SUBROUTINE wrikep

! ===================================================
! velocity at infinity with repect to Earth
! (circular approx) in AU/day
! this is a duplicate of the one in tp_trace, to be used by wrikep
DOUBLE PRECISION FUNCTION v_infty0(el0)
! input elements
  TYPE(orbit_elem), INTENt(IN) :: el0
! END INTERFACE
  TYPE(orbit_elem) eleq
  DOUBLE PRECISION eq0(6) 
  DOUBLE PRECISION cosi,v2
  CHARACTER*3 coo
  INTEGER fail_flag
  coo=el0%coo
  CALL coo_cha(el0,'EQU',eleq,fail_flag)
  eq0=eleq%coord
  IF(fail_flag.eq.0)THEN 
! v_infinity computation in equinoctal
     cosi=(1.d0-eq0(4)**2-eq0(5)**2)/                                  &
     &     (1.d0+eq0(4)**2+eq0(5)**2)                                   
     v2=3.d0-1.d0/eq0(1) -                                             &
     &     2.d0*sqrt(eq0(1)*(1.d0-eq0(2)**2-eq0(3)**2))*cosi            
  ELSEIF(fail_flag.eq.5)THEN
! v_infinity computation in cometary
     cosi=cos(eq0(3))
     v2=3.d0-(1.d0-eq0(2))/eq0(1) - 2.d0*sqrt(eq0(1)*(1.d0+eq0(2)))*cosi
  ELSE
     WRITE(*,*)' v_infty: coordinate conversion failed ', el0,fail_flag  
  ENDIF
! handle non hyperbolic (w.r. to Earth) case
  IF(v2.gt.0.d0)THEN 
     v_infty0=sqrt(v2) 
  ELSE 
     v_infty0=0.d0 
  ENDIF
! normalization in AU/day by Gauss constant
  v_infty0=v_infty0*gk
! conversion to km/s required by wrikep
   v_infty0=v_infty0*aukm/86400.d0
END FUNCTION v_infty0

! ===================================
! FIND_ELEMS
! find elements in .eq0 file, if necessary convert into
! required coordinate system, and report on availability
!  ===================================
SUBROUTINE find_elems(nam0,eledir,buildpath,coox,iunout,el0,unc0,ini0,cov0)
  CHARACTER*(*),INTENT(IN) :: nam0 ! asteroid name, could be long
  CHARACTER*60, INTENT(IN) ::  eledir ! root directory for .eq0 files
  LOGICAL, INTENT(IN):: buildpath ! to compose path with eledir, if not use 
                                  ! eledir as complete pathname  
  CHARACTER*3, INTENT(IN) :: coox ! in which coords the elements are required
   ! if coox=' ' coords are not changed, but left as they are
  INTEGER,INTENT(IN) :: iunout ! log output unit  
  TYPE(orbit_elem), INTENT(OUT) :: el0 ! epoch time (MJD), elements, abs. mag. 
  TYPE(orb_uncert), INTENT(OUT) :: unc0 ! covariance and normal matrices 
  LOGICAL,INTENT(OUT) :: ini0,cov0 ! successful input flag, covar. available
! === END INTERFACE =====================
  INCLUDE 'sysdep.h90'
  LOGICAL eof, error
  CHARACTER*60 elefi0 ! elements  file name
  INTEGER le, fail_flag
  CHARACTER*19 name
  CHARACTER*120 file
  TYPE(orbit_elem) :: elk
  TYPE(orb_uncert) :: unck 
  DOUBLE PRECISION dee(6,6)
!
  ini0=.false. 
  cov0=.false.
  IF(buildpath)THEN
     CALL fidinam(eledir,nam0,'eq0',elefi0,le)
  ELSE
     file=eledir//dircha//nam0//'.eq0'
     CALL rmsp(file,le)
     elefi0(1:le)=file(1:le)
  ENDIF  
  INQUIRE(file=elefi0(1:le),exist=ini0) 
  IF(ini0)THEN 
     CALL read_elems(elk,name,eof,elefi0(1:le),COVAR=unck) 
     ini0=(.not.eof).and.(name.eq.nam0) 
     IF(ini0)THEN
        cov0=unck%succ
        WRITE(iunout,*)nam0,' OK'
     ELSE
        WRITE(*,*)elefi0(1:le),'not read properly'
         CALL write_err(nam0,iunout,' ELEMENTS NOT FOUND')
     ENDIF
  ELSE 
     WRITE(*,*)elefi0(1:le),'not found'
     CALL write_err(nam0,iunout,' ELEMENT FILE NOT FOUND') 
  ENDIF
  IF(coox.eq.' ')RETURN
  IF(ini0)THEN 
! coordinate change to equinoctal
     IF(cov0)THEN
        CALL coo_cha(elk,coox,el0,fail_flag,dee)
        CALL convertunc(unck,dee,unc0)
        cov0=unc0%succ
     ELSE                                       
        CALL coo_cha(elk,coox,el0,fail_flag)
     ENDIF
     IF(fail_flag.gt.5)THEN
        CALL write_err(nam0,iunout,' FAILED COORD CHANGE ')
        WRITE(ierrou,*) elk%coo,fail_flag,coox
        WRITE(iunout,*) elk%coo,fail_flag,coox
        ini0=.false.
        cov0=.false.
     ELSEIF(fail_flag.eq.5.and.(coox.eq.'EQU'.or.coox.eq.'KEP'))THEN
        CALL write_err(nam0,iunout,' HYPERBOLIC COORD CHANGE ')
        WRITE(ierrou,*) elk%coo,fail_flag,coox
        WRITE(iunout,*) elk%coo,fail_flag,coox
        ini0=.false.
        cov0=.false.
     ENDIF
!  no covariances available 
  ENDIF
END SUBROUTINE find_elems


END MODULE orbit_elements

! ===================================================================== 
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                         A P P M A G                         *      
!  *                                                             *      
!  *     Calcolo della magnitudine apparente di un asteroide     *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
!                                                                       
! INPUT:    H         -  Magnitudine assoluta                           
!           G         -  Parametro di slope                             
!           DS        -  Distanza dal Sole (AU)                         
!           DT        -  Distanza dalla Terra (AU)                      
!           BETA      -  Angolo di fase solare (rad): angolo tra il     
!                        Sole e l'osservatore (Terra), visto dall'astero
!                                                                       
DOUBLE PRECISION FUNCTION appmag(h,g,ds,dt,beta) 
  IMPLICIT NONE
  DOUBLE PRECISION h,g,ds,dt,beta 
  DOUBLE PRECISION a1,a2,b1,b2 
  DOUBLE PRECISION tb2,phi1,phi2 
  SAVE a1,a2,b1,b2
! Costanti per il calcolo della magnitudine                             
  DATA a1,a2,b1,b2/3.33d0,1.87d0,0.63d0,1.22d0/ 
  tb2=TAN(beta/2.d0) 
  phi1=EXP(-a1*tb2**b1) 
  phi2=EXP(-a2*tb2**b2) 
  appmag=5.d0*LOG10(ds*dt)+h-2.5d0*LOG10((1.d0-g)*phi1+g*phi2) 
END FUNCTION appmag
