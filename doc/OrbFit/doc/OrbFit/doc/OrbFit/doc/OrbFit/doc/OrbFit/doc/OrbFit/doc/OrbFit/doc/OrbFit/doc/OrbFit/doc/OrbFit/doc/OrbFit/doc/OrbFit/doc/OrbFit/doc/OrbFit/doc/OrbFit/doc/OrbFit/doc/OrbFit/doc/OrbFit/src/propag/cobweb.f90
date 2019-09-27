!======================================================================!
! MODULE COBWEB                                                        !
!======================================================================!
! Authors: G. Tommei, F. Spoto, A. Del Vigna, A. Milani                !
! Date of creation: Nov. 2004                                          !
! Last update: June 2018 (A. Del Vigna)                                !
!                                                                      !
! References:                                                          !
! [1] Spoto, Del Vigna (2018): "Short arc orbit determination and      !
!                               imminent impactors in the Gaia era"    !
!                                                                      !
! [2] Del Vigna (2018): "On Impact Monitoring of Near-Earth Asteroids" !
!                       (PhD Thesis)                                   !
!======================================================================!
! Subroutines contained:                                               !
!    - sample_web_grid                                                 !
!    - syst_ranging                                                    !
!    - sample_grid                                                     !
!    - spider                                                          !
!    - mov                                                             !
!    - write_output_files                                              !
!    - compute_score                                                   !
!    - immediate_imp                                                   !
!    - mc_res_att                                                      !
!======================================================================!
MODULE cobweb
  USE output_control
  USE fund_const
  USE astrometric_observations
  USE orbit_elements
  USE least_squares
  USE name_rules

  IMPLICIT NONE

  PRIVATE
  !*******************************!
  !  PUBLIC PARAMETERS/VARIABLES  !
  !*******************************!
  INTEGER, PARAMETER, PUBLIC :: nmax  = 250       ! Maximum number of points along one axis
  INTEGER, PARAMETER, PUBLIC :: nmax2 = nmax*nmax ! Maximum number of points in the grid
  
  DOUBLE PRECISION, PUBLIC :: sigx                      ! Max of chi for impact search (now set to 5.d0 in fitobs.def)
  DOUBLE PRECISION, PUBLIC :: hmax                      ! Max magnitude for interesting object (now 34.5)
  INTEGER,          PUBLIC :: ndir                      ! Number of rho sample pts
  INTEGER,          PUBLIC :: np                        ! Number of rho-dot sample pts
  INTEGER,          PUBLIC :: ndir2                     ! Number of rho sample pts (densified grid)
  INTEGER,          PUBLIC :: np2                       ! Number of rho-dot sample pts (densified grid)
  LOGICAL,          PUBLIC :: grid_lev_curve            ! If .TRUE., a 50x50 grid is computed even in cobweb case
  LOGICAL,          PUBLIC :: propag_geoc_orbit         ! If .TRUE. the geocentric orbits are propagated for VI search
  LOGICAL,          PUBLIC :: mov_flag                  ! True if the MOV is available (used in fitobs)
  DOUBLE PRECISION, PUBLIC :: elong                     ! Elongation from the Sun (output of the AR computation)
  DOUBLE PRECISION, PUBLIC :: xogeo(3)                  ! Geocentric observer position (Equatorial)
  DOUBLE PRECISION, PUBLIC :: vogeo(3)                  ! Geocentric observer velocity (Equatorial)
  TYPE(orbit_elem), PUBLIC :: el_mov(0:nmax,0:nmax)     ! Corrected elements on MOV (output of mov)
  DOUBLE PRECISION, PUBLIC :: rms_points(0:nmax,0:nmax) ! RMS of the MOV orbits
  LOGICAL,          PUBLIC :: succ_cob(0:nmax,0:nmax)   ! Success flag for the MOV orbit computation

  !********************!
  !  COMMON VARIABLES  !
  !********************!
  ! Parameters
  DOUBLE PRECISION, PARAMETER :: rho_max_unif = 3.16227766016838d0 ! sqrt(10) Maximum rho value to switch to a unif. step in rho
  ! AR computation
  LOGICAL          :: use_nominal               ! If true, do cobweb; if false, grid is used
  LOGICAL          :: select_risk               ! True when an orbit with cov. is available and geoc_chi > 1
  DOUBLE PRECISION :: E_lim                     ! Limit value for the energy
  LOGICAL          :: log_grid                  ! Logarithmic grid in rho
  DOUBLE PRECISION :: r_rdot(0:nmax,0:nmax,2)   ! Coordinates of the cobweb/grid points
  DOUBLE PRECISION :: r_cob(nmax)               ! Array of distances (cobweb case)
  DOUBLE PRECISION :: theta(nmax)               ! Array of angles (cobweb case)
  DOUBLE PRECISION :: e_earth(0:nmax,0:nmax)    ! Two-body energy w.r. to Earth
  DOUBLE PRECISION :: e_sun(0:nmax,0:nmax)      ! Two-body energy w.r. to Sun
  ! MOV computation
  TYPE(orb_uncert) :: unc6_mov(0:nmax,0:nmax)   ! 6-dim uncertainty matrices (output of mov)
  DOUBLE PRECISION :: chi_cob(0:nmax,0:nmax)    ! Chi after the MOV computation
  DOUBLE PRECISION :: delnor_cob(0:nmax,0:nmax) ! Diff. corrections norm
  DOUBLE PRECISION :: jac_sigma(0:nmax,0:nmax)  ! Jacobian det. for the map from the sampling space to the AR
  DOUBLE PRECISION :: jac_mu(0:nmax,0:nmax)     ! Determinant of the map from the AR to the MOV
  LOGICAL          :: inside_AR(0:nmax,0:nmax)  ! Flag which is TRUE if rho > 0 and E < E_lim
  LOGICAL          :: good_first                ! Flag whith value TRUE if the first grid produces at least one point with chi<5
  INTEGER          :: succ_flag(0:nmax,0:nmax)  ! 0 = diff. corr. fails, 2 = not in AR, 1 = good
  ! Imminent impactors analysis
  DOUBLE PRECISION :: tcla(0:nmax,0:nmax)       ! Time of close approach of the (i,j) MOV orbit
  DOUBLE PRECISION :: dcla(0:nmax,0:nmax)       ! Closest approach distance of the (i,j) MOV orbit
  DOUBLE PRECISION :: imppro(0:nmax,0:nmax)     ! Impact probability 
  DOUBLE PRECISION :: csi(0:nmax,0:nmax)        ! Xi coordinate on the MTP
  DOUBLE PRECISION :: zeta(0:nmax,0:nmax)       ! Zeta coordinate on the MTP

  !*******************!
  !  PUBLIC ROUTINES  !
  !*******************!
  PUBLIC :: spider, mov, sample_web_grid, immediate_imp, mc_res_att, compute_score

CONTAINS

  !============================================================!                  
  ! SAMPLE_WEB_GRID                                            !
  !============================================================!
  ! Sampling of admissible region, either by cobweb or by grid !
  !============================================================!
  SUBROUTINE sample_web_grid(astnam,ini0,cov0,elc,unc, &
       &                     csinor,m,obs,obsw,nd,artyp,geoc_chi,acce_chi,chi)
    USE attributable 
    USE arc_control 
    USE triangles
    !=======================================================================================================
    CHARACTER*(name_len), INTENT(IN)    :: astnam      ! Asteroid name
    LOGICAL,              INTENT(IN)    :: ini0        ! Availability of orbital elements
    LOGICAL,              INTENT(IN)    :: cov0        ! Availability of covariance
    TYPE(orbit_elem),     INTENT(INOUT) :: elc         ! Nominal elements
    TYPE(orb_uncert),     INTENT(IN)    :: unc         ! Uncertainty matrices
    DOUBLE PRECISION,     INTENT(IN)    :: csinor      ! Residuals norm
    INTEGER,              INTENT(IN)    :: m           ! Number of observations
    TYPE(ast_obs),        INTENT(IN)    :: obs(m)      ! Observations
    TYPE(ast_wbsr),       INTENT(INOUT) :: obsw(m)     ! Weights, residuals
    INTEGER,              INTENT(IN)    :: nd          ! Dimension of the par space
    INTEGER,              INTENT(OUT)   :: artyp       ! Arc type
    DOUBLE PRECISION,     INTENT(OUT)   :: geoc_chi    ! Normalized value of the geodetic curvature
    DOUBLE PRECISION,     INTENT(OUT)   :: acce_chi    ! Normalized value for the along-track acceleration
    DOUBLE PRECISION,     INTENT(OUT)   :: chi         ! Chi-value of geodetic curvature and acceleration
    !===== Arc type computation ============================================================================
    TYPE(attrib)                :: attrc                 ! Attributable
    LOGICAL                     :: error                 ! Error in the attributable computation
    DOUBLE PRECISION            :: trou                  ! Rounded time for the attributable (at the integer MJD)
    INTEGER                     :: nigarc                ! Number of nights for the arc type
    INTEGER                     :: fail_arty             ! Arc type computation failure
    CHARACTER(LEN=3)            :: staz                  ! Observatory code (from attributable computation)
    DOUBLE PRECISION, PARAMETER :: sphx=2.d0             ! Maximum arc span (in degrees)
    !===== Admissible Region computation ===================================================================
    DOUBLE PRECISION, PARAMETER :: a_max=100.d0          ! Maximum value for the semimajor axis
    DOUBLE PRECISION            :: refs(3,3)             ! Matrix of the reference system
    DOUBLE PRECISION            :: E_bound               ! Bound for E_Sun
    DOUBLE PRECISION            :: c(0:5)                ! Coefficients 
    INTEGER                     :: nrootsf               ! Number of roots of the 6 deg. poly. (V(r), comments below)
    DOUBLE PRECISION            :: roots(3)              ! Roots of the polynomial V(r)
    INTEGER                     :: nrootsd               ! Number of roots of the pol. derivative V'(r)
    DOUBLE PRECISION            :: rootsd(8)             ! Roots of the pol. derivative V'(r)
    DOUBLE PRECISION            :: xo(3)                 ! Heliocentric observer position (Equatorial)     
    DOUBLE PRECISION            :: vo(3)                 ! Heliocentric observer velocity (Equatorial)
    DOUBLE PRECISION            :: rmin                  ! Minimum value for rho
    INTEGER                     :: ind                   ! Index of a zero of V'(r) between rmin and the first pos. root
    DOUBLE PRECISION            :: rder0                 ! Zero of V'(r) between rmin and the first pos. root
    INTEGER                     :: ntot                  ! Total number of points for the boundary computation
    DOUBLE PRECISION            :: frv(npox)             ! Metric value on the boundary of the AR
    DOUBLE PRECISION            :: rhodot(npox)          ! Values of rho_dot in the boundary sampling
    DOUBLE PRECISION            :: rdw                   ! Scaling factor for range in rhodot
    !===== Cobweb/Grid computation =========================================================================
    DOUBLE PRECISION            :: rdotmax               ! Maximum value for rho_dot
    DOUBLE PRECISION            :: rdotmin               ! Minimum value for rho_dot
    DOUBLE PRECISION            :: rho                   ! Rho values in the grid
    DOUBLE PRECISION            :: rmax                  ! Maximum value for rho
    DOUBLE PRECISION            :: rms                   ! RMS for cobweb/grid computation
    DOUBLE PRECISION            :: jacthr(0:nmax,0:nmax) ! Jacobian of (R,theta) -> (rho,rho_dot)
    DOUBLE PRECISION            :: score_nea             ! Score point: NEA
    DOUBLE PRECISION            :: score_mba             ! Score point: MBA 
    DOUBLE PRECISION            :: score_distant         ! Score point: distant
    DOUBLE PRECISION            :: score_scatt           ! Score point: scattered 
    !===== File variables ==================================================================================
    CHARACTER(LEN=100)          :: file                  ! File name
    INTEGER                     :: le                    ! Lenght of file
    INTEGER                     :: iunfla                ! Unit for the att file
    INTEGER                     :: iunpol                ! Unit for the pol file (AR computation output)
    INTEGER                     :: iuncob                ! Unit for the web file (info on the sampling points)
    INTEGER                     :: iuncom                ! Unit for the com file (info on the sampling points)
    INTEGER                     :: iunele                ! Unit for the kepele file (KEP el. for the sampling points)
    INTEGER                     :: iunscore              ! Unit for the score file
    !===== Functions =======================================================================================
    DOUBLE PRECISION            :: prscal                ! Scalar product
    DOUBLE PRECISION            :: vsize                 ! Norm of a vector
    !===== Other variables =================================================================================
    DOUBLE PRECISION            :: crho(2,2)             ! Submatrix of C for (rho,rho_dot)
    DOUBLE PRECISION            :: grho(2,2)             ! Submatrix of G for (rho,rho_dot)
    DOUBLE PRECISION            :: att(4)                ! Attributable for the energy computation
    TYPE(orbit_elem)            :: elk                   ! Keplerian elements of the sampling points
    TYPE(orbit_elem)            :: elcom                 ! Cometarian elements of the sampling points
    INTEGER                     :: fail_flag             ! Failure of the coordinate change
    CHARACTER(LEN=20)           :: mulname               ! Name of the sampling point
    INTEGER                     :: i, j, k, ii, jj, iii  ! Loop indexes
    !=======================================================================================================
    IF(nd.NE.6)THEN
       WRITE(*,*) 'sample_web_grid: not ready for non-grav, nd = ',nd
       STOP
    END IF
    !******************!
    !  Inizialization  !
    !******************!
    log_grid  = .FALSE.
    r_cob     = 0.d0
    theta     = 0.d0
    ! Energy curve level (orbit with a = a_max)
    E_lim = -gms/(2.d0*a_max) 
    ! Determinants 
    jac_sigma = 0.d0
    jac_mu    = 0.d0
    !**************************!
    ! Compute the attributable !
    !**************************!
    CALL tee(iun_log,'------------------------=')
    CALL tee(iun_log,'| COMPUTE ATTRIBUTABLE |=')
    CALL tee(iun_log,'------------------------=')
    CALL attri_comp(m,obs,obsw,attrc,error)
    IF(error)THEN
       WRITE(*,*) 'sample_web_grid: ERROR! In computing the attributable'
       STOP
    END IF
    trou = NINT(attrc%tdtobs)
    ! Write the attributable computation output
    IF(attrc%sph*radeg.GT.sphx)THEN
       WRITE(*,*) 'sample_web_grid: arc too wide = ', attrc%sph*radeg
    ELSE
       CALL wri_attri(0,0,astnam,attrc,trou)
       CALL wri_attri(iun_log,iun_log,astnam,attrc,trou)
    ENDIF
    !**********************!
    ! Compute the arc type !
    !**********************!
    CALL tee(iun_log,'--------------------=')
    CALL tee(iun_log,'| COMPUTE ARC TYPE |=')
    CALL tee(iun_log,'--------------------=')
    artyp = arc_type(obs,obsw,m,geoc_chi,acce_chi,chi,nigarc,fail_arty)
    WRITE(*,176) artyp
    WRITE(*,177) m,nigarc,fail_arty,geoc_chi,acce_chi,chi
    WRITE(iun_log,176) artyp
    WRITE(iun_log,177) m,nigarc,fail_arty,geoc_chi,acce_chi,chi
176 FORMAT(' The arc type is ',I3)
177 FORMAT(' nobs = ',I4,', nights = ',I3,', errcode = ',I4/ &
         & ' geoc_chi = ',1P,D9.2,', acce_chi = ',D9.2,', chi = ',D9.2)
    !*******************************!
    ! Compute the Admissible Region !
    !*******************************!
    !------------------------------------------------------------------!
    ! The most important condition defining the AR is the              !
    ! non-positivity of the two body energy of the body w.r.t the      !
    ! Sun. This condition gives solutions for r_dot if and only if the !
    ! following inequality holds:                                      !
    !                                                                  !
    !                          V(r) <= 4k^4                            !
    !                                                                  !
    ! where V(r) is a 6 degree polynomial. We can have only three      !
    ! possibilities:                                                   !
    !                                                                  !
    ! (1) V(r) has 4 simple real roots -> the AR has 2 connected       !
    !     components;                                                  !
    !                                                                  !
    ! (2) V(r) has 3 distinct real roots, 2 simple and 1 with even     !
    !     multiplicity -> the second c.c. reduces to a point;          !
    !                                                                  !
    ! (3) V(r) has 2 distinct real roots, 1 simple and 1 with odd      !
    !     multiplicity -> the AR has 1 c.c.                            !
    !                                                                  !
    ! Moreover, V'(r) cannot have more than 3 distinct roots. If they  !
    ! are exactly three, then there cannot be any root with            ! 
    ! multiplicity 2.                                                  !
    !------------------------------------------------------------------!
    !---------------------------------------------------!
    !  Open .att file for output of the AR computation  !
    !---------------------------------------------------!
    file = astnam//'.att'
    CALL rmsp(file,le)
    CALL filopn(iunfla,file(1:le),'unknown')
    staz = attrc%obscod
    CALL tee(iun_log,'-----------------------------=')
    CALL tee(iun_log,'| COMPUTE ADMISSIBLE REGION |=')
    CALL tee(iun_log,'-----------------------------=')
    CALL admis_reg(astnam,iunfla,attrc%tdtobs,attrc%angles,staz,nrootsf,roots,c, &
         &         refs,xo,vo,a_max,E_bound,nrootsd,rootsd,attrc%apm,xogeo,vogeo,elong)
    CALL filclo(iunfla,' ')
    ! Output the roots
    WRITE(*,*) 'Roots = ', roots(1:nrootsf) 
    WRITE(iun_log,*) 'Roots = ', roots(1:nrootsf) 
    !----------------------------------------!
    !  Find the range of values for rho_dot  !
    !----------------------------------------!
    rmin = 1.d-4    ! Arbitrary lower boundary for the AR
    ! Search for a zero (rder0) of the derivative V'(r) between rmin and roots(1) 
    rder0 = -1.d0
    ind   = 0
    IF(nrootsd.GE.1)THEN
       IF(nrootsd.EQ.1)THEN
          ! Do nothing
       ELSEIF(nrootsd.EQ.2)THEN
          ! Do nothing
       ELSEIF(nrootsd.EQ.3)THEN
          DO iii=1,nrootsd
             IF(rootsd(iii).GT.rmin .AND. rootsd(iii).LT.roots(1))THEN
                ind=iii
             ENDIF
          ENDDO
          IF(ind.GT.2)THEN
             IF(rootsd(ind-2).GT.rmin)THEN
                rder0=rootsd(ind-2)
                WRITE(*,*) 'rmin = ',rmin,', rder0 = ',rder0,', roots(1) = ',roots(1)
             ENDIF
          ENDIF
          IF(rder0.GE.roots(1))THEN
             WRITE(*,*) 'sample_web_grid: ERROR! rder0 > roots(1)'
          ENDIF
       ELSE
          ! V'(r) cannot have more than 3 distinct roots
          WRITE(*,*) 'sample_web_grid: roots of the der. of V(r) (nrootsd) = ',nrootsd
          WRITE(*,*) '                 It cannot be possible!'
       ENDIF
    ENDIF
    nfunc = 2 ! Use logaritmic scale for the metric used in the boundary sampling
    CALL sample_bound_ne1(E_bound,attrc%eta**2,rmin,rder0,roots(1),8,5, &
         &                ntot,c,frv,rhodot,npox,rdw)
    rdotmin = MINVAL(rhodot(1:ntot))
    rdotmax = MAXVAL(rhodot(1:ntot))
    !-------------------------------------------------!
    !  Open .pol file, containing the data of the AR  !
    !-------------------------------------------------!
    file=astnam//'.pol'
    CALL rmsp(file,le)
    CALL filopn(iunpol,file(1:le),'unknown')
    ! Write header
    WRITE(iunpol,'(A)') '% N. roots      Roots(1)              Roots(2)              Roots(3)              &
         &c(0)                 c(1)                c(2)                 c(3)                 c(4)          &
         &       c(5)                rdotmin              rdotmax'

    WRITE(iunpol,178) nrootsf, roots, c, rdotmin, rdotmax
178 FORMAT(6X,I1,2X,11(1X,F20.16)) 
    CALL filclo(iunpol,' ')
    !*************************!
    ! Cobweb/grid computation !
    !*************************!
    !------------------------------------------------------------------!
    ! The sampling of the AR depends on the availability of a nominal  !
    ! solution. If there is a nominal solution, with its covariance,   !
    ! and if the SNR of the geodetic curvature is > 3, then we trust   !
    ! the nominal solution and we do a cobweb sampling.                !
    !------------------------------------------------------------------!
    use_nominal = (ini0 .AND. cov0 .AND. ABS(geoc_chi).GT.3.d0)
    select_risk = (ini0 .AND. cov0 .AND. ABS(geoc_chi).GT.1.d0)
    DO k=1,2
       IF(k.EQ.1)THEN
          ! When grid_lev_curve is false, skip the grid computation
          IF(use_nominal .AND. .NOT.grid_lev_curve) CYCLE
          !-----------------------------------------------------------------!
          !                  First iteration (first grid)                   !
          !-----------------------------------------------------------------!
          ! Initialization
          good_first = .TRUE.
          ! Selection of logarithmic grid in rho
          IF(roots(1).LT.rho_max_unif) log_grid = .TRUE.
          ! Selection of other grid parameters
          CALL select_grid_param(k,nrootsf,roots(1:3),rmin,rmax,rdotmin,rdotmax)
          ! Save values of the first iteration
          !--------------------!
          ! Systematic ranging !
          !--------------------!
          CALL syst_ranging(astnam,m,obs,obsw,nd,attrc,elc,unc,csinor,c(0:5),nrootsf,roots(1:3), &
               &            rmin,rmax,rdotmin,rdotmax,k)
          !********************************!
          ! Compute_score (only grid case) !
          !********************************!
          IF(.NOT.use_nominal)THEN
             CALL compute_score(score_distant,score_mba,score_nea,score_scatt)
          END IF
       ELSEIF(k.EQ.2)THEN
          IF(.NOT.use_nominal)THEN
             !-------------------------------------------------------------------!
             !                Second iteration (densified grid)                  !
             !-------------------------------------------------------------------!
             !-------------------------------------------------------------!
             ! Once we obtain a first preliminary grid, we densify for a   !
             ! higher resolution. We select the minimum and the maximum    !
             ! value of rho and rho_dot among all the values of the        !
             ! points for which we had convergence.                        !
             !-------------------------------------------------------------!
             ! Selection of a logarithmic sampling if rmax is too big (only with one component)
             IF(nrootsf.EQ.1)THEN
                log_grid=.FALSE.
                IF(score_nea.GT.0.5d0) log_grid = .TRUE.
             END IF
             ! Selection of other grid parameters
             CALL select_grid_param(k,nrootsf,roots(1:3),rmin,rmax,rdotmin,rdotmax)
             ! No good point in the first iteration
             IF(.NOT.good_first)THEN
                CALL tee(iun_log,'sample_web_grid: WARNING! No good points in the previous iteration=')
                CALL tee(iun_log,'                 Output files .comele and .web equal to _1.comele and _1.web=')
                CALL write_output_files(k,roots(2),astnam,-1.d0,csinor,undefined_orbit_elem)
                EXIT
             END IF
             !------------------------------------------------!
             ! WARNING: this part is not useful anymore (ADV) !
             !------------------------------------------------!
             ! If rmin and rmax are equal (vertical line) give more space
             !IF(rmin.EQ.rmax)THEN
             !   rmin = rmin*0.95d0 ! Not done with sums to avoid negative values
             !   rmax = rmax*1.05d0 ! Symmetric
             !END IF
             ! If rdotmin and rdotmax are equal (horizontal line) give more space
             !IF(rdotmin.EQ.rdotmax)THEN
             !   rdotmin = rdotmin - 0.005
             !   rdotmax = rdotmax + 0.005
             !END IF
             !--------------------!
             ! Systematic ranging !
             !--------------------!
             CALL syst_ranging(astnam,m,obs,obsw,nd,attrc,elc,unc,csinor,c(0:5),nrootsf,roots(1:3), &
                  &            rmin,rmax,rdotmin,rdotmax,k)
          ELSE
             !-------------------------------------------------------------------!
             !                    Second iteration (cobweb)                      !
             !-------------------------------------------------------------------!
             ndir = ndir2   ! Number of directions
             np   = np2     ! Number of points per direction
             ! Inizialization
             r_rdot  = 0.d0
             e_earth = 0.d0
             e_sun   = 0.d0
             CALL syst_ranging(astnam,m,obs,obsw,nd,attrc,elc,unc,csinor,c(0:5),nrootsf,roots(1:3), &
                  &            rmin,rmax,rdotmin,rdotmax,k)
          END IF
       END IF
    END DO
    !***************!
    ! Compute score !
    !***************!
    CALL compute_score(score_distant,score_mba,score_nea,score_scatt)
    ! Open output file
    file = astnam//'.score'
    CALL rmsp(file,le)
    CALL filopn(iunscore,file(1:le),'unknown')
    ! Write header score file
    WRITE(iunscore,'(A)') '% Object '//astnam//' score'
    WRITE(iunscore, 595) 'Score Near-Earth : ', score_nea
    WRITE(iunscore, 595) 'Score Main Belt  : ', score_mba
    WRITE(iunscore, 595) 'Score Distant    : ', score_distant
    WRITE(iunscore, 595) 'Score Scattered  : ', score_scatt
595 FORMAT(A, 3X, F4.2)
    CALL filclo(iunscore,' ')

  END SUBROUTINE sample_web_grid
  
  
  !============================================================!                  
  ! SELECT_GRID_PARAM                                          !
  !============================================================!
  ! Selection of grid end-points.                              !
  !============================================================!
  SUBROUTINE select_grid_param(k_ind,n_roots,roots,rmin,rmax,rdotmin,rdotmax)
    !=======================================================================================================
    INTEGER,          INTENT(IN)     :: k_ind     ! Iteration index
    INTEGER,          INTENT(IN)     :: n_roots   ! Number of roots of the 6 deg. poly.
    DOUBLE PRECISION, INTENT(IN)     :: roots(3)  ! Roots of the polynomial defining the AR
    DOUBLE PRECISION, INTENT(INOUT)  :: rmin      ! Starting point for the rho sampling
    DOUBLE PRECISION, INTENT(INOUT)  :: rmax      ! Ending point for the rho sampling
    DOUBLE PRECISION, INTENT(INOUT)  :: rdotmin   ! Starting point for the rhodot sampling
    DOUBLE PRECISION, INTENT(INOUT)  :: rdotmax   ! Ending point for the rhodot sampling
    !=======================================================================================================
    DOUBLE PRECISION, PARAMETER :: rmin_aux = 1000.d0     ! Auxiliary rmin
    DOUBLE PRECISION, PARAMETER :: rmax_aux = 0.d0        ! Auxiliary rmax
    DOUBLE PRECISION, PARAMETER :: rdotmin_aux = 1000.d0  ! Auxiliary rdotmin
    DOUBLE PRECISION, PARAMETER :: rdotmax_aux = -1000.d0 ! Auxiliary rdotmax
    LOGICAL                     :: same_side              ! TRUE when rmin and rmax are on the same side w.r.t. r_1
    INTEGER                     :: i_min,i_max            ! Index of min and max rho
    INTEGER                     :: j_min,j_max            ! Index of min and max rho_dot
    INTEGER                     :: i,j                    ! Loop indexes
    !=======================================================================================================    
    IF(k_ind.EQ.1)THEN
       !-----------------------------------------------!
       ! rmax, rmin: maximum and minumum value for rho !
       !-----------------------------------------------!
       IF(n_roots.EQ.3)THEN
          rmin = 1.d-3     ! rmin passes from 1.d-4 to 1.d-3
          rmax = roots(3)
       ELSE
          rmax = roots(1)
       END IF
       CALL tee(iun_log,'----------------=')
       CALL tee(iun_log,'| COMPUTE GRID |=')
       CALL tee(iun_log,'----------------=')
    ELSE
       !--------------------------------------------------------------------------!
       ! rmax, rmin: Search for the minimum and the maximum rho for which chi < 5 ! 
       !--------------------------------------------------------------------------!
       ! Initialize AR sampling boundaries (*_aux defined above as parameters)
       rmin    = rmin_aux
       rmax    = rmax_aux
       rdotmin = rdotmin_aux
       rdotmax = rdotmax_aux
       DO i=1,ndir
          DO j=1,np
             IF(chi_cob(i,j).LT.5.d0 .AND. el_mov(i,j)%h_mag.LT.hmax .AND. succ_cob(i,j))THEN
                IF(r_rdot(i,j,1).LT.rmin)THEN
                   rmin  = r_rdot(i,j,1)
                   i_min = i
                END IF
                IF(r_rdot(i,j,1).GT.rmax)THEN
                   rmax  = r_rdot(i,j,1)
                   i_max = i
                END IF
                IF(r_rdot(i,j,2).LT.rdotmin)THEN
                   rdotmin = r_rdot(i,j,2)
                   j_min   = j
                END IF
                IF(r_rdot(i,j,2).GT.rdotmax)THEN
                   rdotmax = r_rdot(i,j,2)
                   j_max   = j
                END IF
             END IF
          END DO
       END DO
       IF(rmin.EQ.rmin_aux .OR. rmax.EQ.rmax_aux .OR. rdotmin.EQ.rdotmin_aux .OR. rdotmax.EQ.rdotmax_aux)THEN
          good_first = .FALSE.
          RETURN
       END IF
       !--------------------------------------------------------------!
       ! To avoid loss of good points due to the lower resolution of  !
       ! the first grid, we enlarge to the nearby columns             !
       !--------------------------------------------------------------!
       IF(i_min.EQ.1) i_min = 2
       IF(j_min.EQ.1) j_min = 2
       IF(i_max.EQ.ndir) i_max = ndir-1
       IF(j_max.EQ.np)   j_max = np-1
       rmin    = r_rdot(i_min-1,1,1)
       rmax    = r_rdot(i_max+1,1,1)
       rdotmin = r_rdot(1,j_min-1,2)
       rdotmax = r_rdot(1,j_max+1,2)
       !----------------!
       ! Densified grid !
       !----------------!
       ndir = ndir2
       np   = np2
       CALL tee(iun_log,'--------------------------=')
       CALL tee(iun_log,'| COMPUTE DENSIFIED GRID |=')
       CALL tee(iun_log,'--------------------------=')
    END IF
    !------------------------------------------------------!
    ! ndir, np: no. of sample points in rho and in rho_dot !
    !------------------------------------------------------!
    ! Position of rmin and rmax w.r.t. r_1
    same_side = same_side_r1(rmin,rmax,roots(1))
    IF(n_roots.EQ.3 .AND. .NOT.same_side)THEN
       ndir = 2*ndir    
    END IF
    !-----------------!
    ! Output messages !
    !-----------------!
    WRITE(*,*) '   Grid kind selection: logarithmic grid = ', log_grid
    WRITE(*,101) ndir,np 
    IF(same_side)THEN
       WRITE(*,102) rmin,rmax 
    ELSE
       WRITE(*,103) rmin,roots(1),roots(2),rmax
    END IF
    WRITE(*,104) rdotmin,rdotmax 
101 FORMAT('    Grid dimensions = ',I4,' x 'I4)
102 FORMAT('    Sampling interval in rho = [',F12.6,','F12.6,']')
103 FORMAT('    Sampling interval in rho = [',F12.6,','F12.6,'] U [',F12.6,','F12.6,']')
104 FORMAT('    Sampling interval in rho_dot = [',F12.6,','F12.6,']')

  END SUBROUTINE select_grid_param


  !====================================================================!                  
  ! SYST_RANGING                                                       !
  !====================================================================!
  ! It performs the systematic ranging on the AR: it samples a grid or !
  ! a spider web, compute the Manifold of Variations and write the     !
  ! output files.                                                      !
  !====================================================================!
  SUBROUTINE syst_ranging(astnam,m,obs,obsw,nd,attrc,el_nom,unc,csinor,coeff,n_roots,roots, &
       &                  rmin,rmax,rhodot_min,rhodot_max,k_ind)
    USE attributable 
    USE station_coordinates, ONLY: statcode
    !=======================================================================================================
    CHARACTER*(name_len), INTENT(IN)    :: astnam      ! Asteroid name
    INTEGER,              INTENT(IN)    :: m           ! Number of observations
    TYPE(ast_obs),        INTENT(IN)    :: obs(m)      ! Observations
    TYPE(ast_wbsr),       INTENT(INOUT) :: obsw(m)     ! Weights, residuals
    INTEGER,              INTENT(IN)    :: nd          ! Dimension of the parameters space
    TYPE(attrib),         INTENT(IN)    :: attrc       ! Attributable
    TYPE(orbit_elem),     INTENT(INOUT) :: el_nom      ! Nominal elements
    TYPE(orb_uncert),     INTENT(IN)    :: unc         ! Uncertainty of the nominal elements
    DOUBLE PRECISION,     INTENT(IN)    :: csinor      ! Residuals norm
    DOUBLE PRECISION,     INTENT(IN)    :: coeff(0:5)  ! Coefficients for the AR computation 
    INTEGER,              INTENT(IN)    :: n_roots     ! Number of roots of the 6 deg. poly.
    DOUBLE PRECISION,     INTENT(IN)    :: roots(3)    ! Roots of the polynomial defining the AR
    DOUBLE PRECISION,     INTENT(IN)    :: rmin        ! Minimum value for rho
    DOUBLE PRECISION,     INTENT(IN)    :: rmax        ! Maximum value for rho
    DOUBLE PRECISION,     INTENT(IN)    :: rhodot_max  ! Maximum value for rho_dot
    DOUBLE PRECISION,     INTENT(IN)    :: rhodot_min  ! Minimum value for rho_dot
    INTEGER,              INTENT(IN)    :: k_ind       ! Loop index, for output
    !=======================================================================================================
    DOUBLE PRECISION  :: rms                   ! RMS for cobweb/grid computation
    DOUBLE PRECISION  :: jacthr(0:nmax,0:nmax) ! Jacobian of (R,theta) -> (rho,rho_dot)
    DOUBLE PRECISION  :: att(4)                ! Attributable for the energy computation
    INTEGER :: i,j
    !=======================================================================================================
    IF(use_nominal .AND. k_ind.EQ.2)THEN
       !-------------------------!
       ! Cobweb case (100 x 100) !
       !-------------------------!
       IF(el_nom%coo.NE.'ATT')THEN
          WRITE(*,*) ' syst_ranging: ERROR! Nominal orbit must be in ATT. It is in ', el_nom%coo
          STOP
       ENDIF
       CALL tee(iun_log,'------------------=')
       CALL tee(iun_log,'| COMPUTE COBWEB |=')
       CALL tee(iun_log,'------------------=')
       r_rdot = 0.d0                       ! Initialization of the prev. computed grid
       r_rdot(0,0,1:2) = el_nom%coord(5:6) ! Initialize r_rdot(0,0)
       rms    = csinor  
       CALL spider(el_nom,unc,csinor,theta,r_cob)
    ELSE
       !-----------!
       ! Grid case !
       !-----------!
       r_rdot =  0.d0    ! Initialization
       rms    = -1.d0    ! RMS fixed value -1 (no nominal)
       ! Sampling in the (rho,rho_dot) space 
       CALL sample_grid(rmin,rmax,rhodot_min,rhodot_max,n_roots,roots)
    END IF
    !-----------------------!
    ! Energy w.r.t. the Sun !
    !-----------------------!
    CALL energy_sun(r_rdot(0:ndir,0:np,1:2),coeff,ndir,np,use_nominal,E_lim,rmin,e_sun(0:ndir,0:np), &
         &          inside_AR(0:ndir,0:np),succ_flag(0:ndir,0:np))
    !************************!
    ! Manifold of Variations !
    !************************!
    ! If a reliable nominal solution does not exist, use only the computed attributable
    IF(.NOT.use_nominal)THEN
       el_nom     = undefined_orbit_elem
       el_nom%coo = 'ATT'                          ! ATT coordinates
       el_nom%t   = attrc%tdtobs                   ! Time of the attributable
       CALL statcode(attrc%obscod,el_nom%obscode)  ! Observatory code
       el_nom%coord(1:4) = attrc%angles            ! Attributable for the first 4 coordinates
       el_nom%coord(5:6) = (/0.d0, 0.d0/)          ! Coordinates (0,0) for (rho,rho_dot)
    END IF
    CALL mov(m,obs,obsw,el_nom,rms,nd,k_ind,el_mov,unc6_mov,delnor_cob,succ_cob,CHI_COB_LOC=chi_cob)
    mov_flag = .TRUE.
    ! Energy w.r.t. the Earth
    att = attrc%angles ! Attributable (angles and derivatives)
    CALL energy_earth(att,xogeo,vogeo,r_rdot(0:ndir,0:np,1:2),ndir,np,use_nominal,e_earth(0:ndir,0:np))
    ! Write output
    CALL write_output_files(k_ind,roots(2),astnam,rms,csinor,el_nom)

  END SUBROUTINE syst_ranging


  !===================================================================!
  ! SAME_SIDE_R1                                                      !
  !===================================================================!
  ! Logical function to assess if rmin and rmax are on the same side  !
  ! with respect to r_1.                                              !
  !===================================================================!
  LOGICAL FUNCTION same_side_r1(rmin,rmax,r_1)
    !=======================================================================================================
    DOUBLE PRECISION, INTENT(IN)  :: rmin      ! Starting point for the rho sampling
    DOUBLE PRECISION, INTENT(IN)  :: rmax      ! Ending point for the rho sampling
    DOUBLE PRECISION, INTENT(IN)  :: r_1       ! First root of the polynomial defining the AR
    !=======================================================================================================
    same_side_r1 = (rmin.LE.r_1 .AND. rmax.LE.r_1*1.05d0) .OR. (rmin.GT.r_1 .AND. rmax.GT.r_1)

  END FUNCTION same_side_r1
  

  !============================================================!                  
  ! SAMPLE_GRID                                                !
  !============================================================!
  ! Grid in the (rho,rho_dot) space.                           !
  !============================================================!
  SUBROUTINE sample_grid(rmin,rmax,rdotmin,rdotmax,n_roots,roots)
    !=======================================================================================================
    DOUBLE PRECISION, INTENT(IN)  :: rmin      ! Starting point for the rho sampling
    DOUBLE PRECISION, INTENT(IN)  :: rmax      ! Ending point for the rho sampling
    DOUBLE PRECISION, INTENT(IN)  :: rdotmin   ! Starting point for the rhodot sampling
    DOUBLE PRECISION, INTENT(IN)  :: rdotmax   ! Ending point for the rhodot sampling
    INTEGER,          INTENT(IN)  :: n_roots   ! Number of roots of the 6 deg. poly.
    DOUBLE PRECISION, INTENT(IN)  :: roots(3)  ! Roots of the polynomial defining the AR
    !=======================================================================================================
    DOUBLE PRECISION :: rms       ! RMS value
    DOUBLE PRECISION :: dr        ! Sampling step in rho
    DOUBLE PRECISION :: dr1       ! Sampling step in rho (first c.c.)
    DOUBLE PRECISION :: dr2       ! Sampling step in rho (second c.c.)
    DOUBLE PRECISION :: drdot     ! Sampling step in rho_dot
    LOGICAL          :: same_side ! TRUE when rmin and rmax are on the same side w.r.t. r_1
    INTEGER          :: i,j       ! Loop indexes
    INTEGER          :: ndir_mid  ! Half of ndir
    !=======================================================================================================
    ! Paranoia check: rmin > 0
    IF(rmin.LT.0)THEN
       WRITE(*,*) 'sample_grid: ERROR! rmin must be > 0'
       STOP
    END IF
    ! Initialization
    r_rdot = 0.d0
    same_side =  same_side_r1(rmin,rmax,roots(1))
    IF(same_side)THEN
       IF(log_grid)THEN
          dr = (LOG10(rmax) - LOG10(rmin))/(ndir - 1)
       ELSE
          dr = (rmax - rmin)/(ndir - 1)
       END IF
       !-----------------!
       ! Sampling of rho !
       !-----------------!
       DO i=1,ndir
          IF(log_grid)THEN
             r_rdot(i,1:np,1) = 10**(dr*(i-1) + LOG10(rmin))
          ELSE
             r_rdot(i,1:np,1) = dr*(i-1) + rmin
          END IF
       ENDDO
    ELSE
       IF(MOD(ndir,2).NE.0)THEN
          CALL tee(iun_log,' The value of ndir must be even!=')
          STOP
       END IF
       ndir_mid = ndir/2
       IF(log_grid)THEN
          dr1 = (LOG10(roots(1)) - LOG10(rmin))/(ndir_mid - 1)    
          dr2 = (LOG10(rmax) - LOG10(roots(2)))/(ndir_mid - 1)
       ELSE
          dr1 = (roots(1) - rmin)/(ndir_mid - 1)    
          dr2 = (rmax - roots(2))/(ndir_mid - 1)
       END IF
       !-----------------!
       ! Sampling of rho !
       !-----------------!
       DO i=1,ndir_mid
          IF(log_grid)THEN
             r_rdot(i,1:np,1)          = 10**(dr1*(i-1) + LOG10(rmin))
             r_rdot(ndir_mid+i,1:np,1) = 10**(dr2*(i-1) + LOG10(roots(2)))
          ELSE
             r_rdot(i,1:np,1)          = dr1*(i-1) + rmin
             r_rdot(ndir_mid+i,1:np,1) = dr2*(i-1) + roots(2)
          END IF
       ENDDO
    END IF
    !---------------------!
    ! Sampling of rho_dot !
    !---------------------!
    drdot = (rdotmax - rdotmin)/(np-1)
    DO j=1,np
       r_rdot(1:ndir,j,2) = rdotmin + drdot*(j-1)
    ENDDO
    !--------------------------------------------!
    ! Determinant of the sampling map (cfr. [1]) !
    !--------------------------------------------!
    IF(log_grid)THEN
       jac_sigma(1:ndir,1:np) = LOG(10.d0)*r_rdot(1:ndir,1:np,1)
    ELSE
       jac_sigma(1:ndir,1:np) = 1.d0
    END IF
  END SUBROUTINE sample_grid

  
  !====================================================================!
  ! SPIDER                                                             !
  !====================================================================!
  ! Computation of the points in the (R,theta) space of polar elliptic !
  ! coordinates                                                        !
  !====================================================================!
  SUBROUTINE spider(el_nom,unc,csinor,theta,r)
    !=======================================================================================================
    TYPE(orbit_elem), INTENT(IN)  :: el_nom                  ! Nominal elements (ATT coord.)
    TYPE(orb_uncert), INTENT(IN)  :: unc                     ! Uncertainty of el_nom
    DOUBLE PRECISION, INTENT(IN)  :: csinor                  ! Residuals norm
    DOUBLE PRECISION, INTENT(OUT) :: theta(nmax)             ! Array of angles
    DOUBLE PRECISION, INTENT(OUT) :: r(nmax)                 ! Array of distances (polar coordinates)
    !=======================================================================================================
    DOUBLE PRECISION :: gamma_rho(2,2) ! Covariance matrix for (rho,rho_dot)
    DOUBLE PRECISION :: eigenval(2)    ! Eigenvalues of gamma_rho
    DOUBLE PRECISION :: eigenvec(2,2)  ! Matrix of eigenvectors of gamma_rho
    DOUBLE PRECISION :: fv1(2),fv2(2)  ! Temporary storage arrays for rs
    DOUBLE PRECISION :: v1(2),v2(2)    ! Eigenvectors of gamma_rho
    DOUBLE PRECISION :: a,b            ! Ellipse's semiaxis
    INTEGER          :: ierr           ! Error flag
    INTEGER          :: i,j            ! Loop indexes
    !=======================================================================================================
    !******************!
    !  Inizialization  !
    !******************!
    r_rdot    = 0.d0
    theta     = 0.d0
    r         = 0.d0
    !***************************************!
    !  Covariance matrix for (rho,rho_dot)  !
    !***************************************!
    gamma_rho(1,1) = unc%g(5,5)
    gamma_rho(1,2) = unc%g(5,6)
    gamma_rho(2,1) = unc%g(6,5)
    gamma_rho(2,2) = unc%g(6,6)
    ! Eigenvalues and eigenvectors of gamma_rho  
    CALL rs(2,2,gamma_rho,eigenval,1,eigenvec,fv1,fv2,ierr)
    IF(ierr.NE.0)THEN 
       WRITE(*,*) 'spider: failed rs, ierr = ',ierr
       STOP
    END IF
    ! Eigenvectors and eigenvalues of gamma_rho
    v1 = eigenvec(:,1)     ! Longest axis direction
    v2 = eigenvec(:,2)
    a  = SQRT(eigenval(1)) ! Semi-major axis
    b  = SQRT(eigenval(2)) ! Semi-minor axis
    !***************************************************************!
    !  Regular grid in (R,theta) plane: elliptic polar coordinates  !
    !***************************************************************!
    !------------------------------------------------------------------!
    ! Create a regular grid of points in the space of polar elliptic   !
    ! coordinates (R,theta), where R is in [0,sigx] (sigx is a         !
    ! parameter defining the maximum value of RMS we consider reliable !
    ! in our analysis), and theta is in [0,2*pi). sigx is currently    !
    ! 5.d0, and it is set into fitobs.def.                             !
    !------------------------------------------------------------------!
    DO i=1,ndir
       theta(i) = (dpig/ndir)*i  
    END DO
    DO j=1,np
       r(j) = (sigx/np)*j 
    END DO
    !***************************************!
    !  Map from (R,theta) to (rho,rho_dot)  !
    !***************************************!
    ! Nominal solution at (i,j)=(0,0)
    r_rdot(0,0,1:2) = el_nom%coord(5:6)
    DO i=1,ndir
       DO j=1,np
          r_rdot(i,j,1)  = r(j)*(a*COS(theta(i))*v1(1) - b*SIN(theta(i))*v1(2)) + el_nom%coord(5)
          r_rdot(i,j,2)  = r(j)*(a*COS(theta(i))*v1(2) + b*SIN(theta(i))*v1(1)) + el_nom%coord(6)
          !--------------------------------------------!
          ! Determinant of the sampling map (cfr. [1]) !
          !--------------------------------------------!
          jac_sigma(i,j) = r(j)*a*b
       END DO
    END DO

  END SUBROUTINE spider
  
  
  !====================================================================!
  ! MOV                                                                !
  !====================================================================!
  ! Computation of the Manifold Of Variations                          !
  !====================================================================!
  SUBROUTINE mov(m,obs,obsw,el_nom,rms,nd,k_ind,el_mov_loc,unc6_mov_loc,delnor_cob_loc,succ_cob_loc,chi_cob_loc,unc4_mov)
    INCLUDE 'parobx.h90' 
    !=======================================================================================================
    INTEGER,          INTENT(IN)            :: m                             ! Number of observations
    TYPE(ast_obs),    INTENT(IN)            :: obs(m)                        ! Observations
    TYPE(ast_wbsr),   INTENT(INOUT)         :: obsw(m)                       ! Weights and observations
    TYPE(orbit_elem), INTENT(IN)            :: el_nom                        ! Nominal elements
    DOUBLE PRECISION, INTENT(IN)            :: rms                           ! Residuals norm
    INTEGER,          INTENT(IN)            :: nd                            ! Dimension of the parameters space
    INTEGER,          INTENT(IN)            :: k_ind                         ! Iteration index
    TYPE(orbit_elem), INTENT(OUT)           :: el_mov_loc(0:nmax,0:nmax)     ! Corrected elements
    TYPE(orb_uncert), INTENT(OUT)           :: unc6_mov_loc(0:nmax,0:nmax)   ! 6-dim uncertainty matrices
    DOUBLE PRECISION, INTENT(OUT)           :: delnor_cob_loc(0:nmax,0:nmax) ! Differential corrections norm
    LOGICAL,          INTENT(OUT)           :: succ_cob_loc(0:nmax,0:nmax)   ! Success flag for each point
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: chi_cob_loc(0:nmax,0:nmax)    ! Chi function
    TYPE(orb_uncert), INTENT(OUT), OPTIONAL :: unc4_mov(0:nmax,0:nmax)       ! Uncertainty mat. of the 4-dim fit
    !===== Computation of the Jacobian =====================================================================
    DOUBLE PRECISION, PARAMETER :: jac_mu_min=1.d-4 ! Minimum value for jac_mu
    TYPE(orb_uncert) :: unc4_mov_loc(0:nmax,0:nmax) ! Uncertainty mat. of the 4-dim fit
    INTEGER          :: scal_obs                    ! Number of scalar observations
    DOUBLE PRECISION :: des_mat(nob2x,nd)           ! Design matrix
    DOUBLE PRECISION :: b_rho(nob2x,2)              ! Submatrix of B relative to (rho,rho_dot)
    DOUBLE PRECISION :: b_A(nob2x,4)                ! Submatrix of B relative to A
    DOUBLE PRECISION :: b_A_rho(4,2)                ! Product of b_A and b_rho
    DOUBLE PRECISION :: dAdrho(4,2)                 ! Derivative of the corrected attributable w.r.t. (rho,rho_dot)
    DOUBLE PRECISION :: M_A(2,2)                    ! Matrix of the map AR -> MOV
    !=======================================================================================================
    TYPE(orbit_elem) :: eltemp(0:nmax,0:nmax) ! Temporary orbit elements
    INTEGER          :: nused(0:nmax,0:nmax)  ! Number of observations used in the differential corrections 
    DOUBLE PRECISION :: rmsh(0:nmax,0:nmax)   ! Weighted RMS of the residuals 
    INTEGER          :: i,j,ii                ! Loop indexes
    LOGICAL          :: succ4                 ! Success of fourdim_fit
    DOUBLE PRECISION :: rms_min               ! Minimum RMS among the succesful sampling points
    INTEGER          :: verb_dif_old          ! Temporary verbosity level 
    INTEGER          :: num_zer               ! Number of null eigenvalues of C_A
    DOUBLE PRECISION :: out_mat(4,4)          ! Output matrix of QR inversion
    DOUBLE PRECISION :: eig_val(4)            ! Eigenvalues of C_A
    DOUBLE PRECISION :: det_aux               ! Auxiliary determinant
    !=======================================================================================================
    !****************!
    ! Inizialization !
    !****************!
    succ_cob_loc = .FALSE. 
    verb_dif_old = verb_dif ! Anti-verbose
    verb_dif     = 1
    !*************************************************************!
    ! Assign the nominal elements to the point with indexes (0,0) !
    !*************************************************************!
    ! Non-negative RMS means that a nominal sol. exists, otherwise the RMS is -1
    IF(rms.GE.0.d0)THEN
       el_mov_loc(0,0)   = el_nom
       rms_points(0,0)   = rms
       chi_cob_loc(0,0)  = 0.d0
       succ_cob_loc(0,0) = .TRUE.    
    ENDIF
    !****************************************************!
    ! 4-dim differential correction along each direction !
    !****************************************************!
    ! Initialization to determine the minimum value of the RMS
    rms_min = 100.d0
    dir_loop: DO i=1,ndir
!    dir_loop: DO i=ndir,1,-1
       eltemp(i,0) = el_nom ! For j=0, give the nominal elements or the attributable
       pt_loop: DO j=1,np
          IF(inside_AR(i,j))THEN
             !---------------------------------------------------------!
             ! If the diff. corr. on the previous points were          !
             ! successful, we give the obtained orbit as a first       !
             ! guess. Otherwise we give the elements el. If j.NE.1 and !
             ! we have a nominal, we use a CYCLE to change direction.  !
             !---------------------------------------------------------!
             IF(j.EQ.1)THEN
                eltemp(i,j) = el_nom
             ELSE
                IF(succ_cob_loc(i,j-1))THEN
                   eltemp(i,j) = eltemp(i,j-1)
                ELSE
                   eltemp(i,j) = el_nom
                   IF(use_nominal .AND. k_ind.EQ.2 .AND. inside_AR(0,0))THEN
                      !---------------------------------------------------------------------!
                      ! Remark: the nominal solution has to be inside the AR, otherwise the !
                      ! ray starting at the nominal intersect the AR only for j>0           !
                      !---------------------------------------------------------------------!
                      WRITE(iun_log,'(A,I3,A,I3,A)') 'mov: point i = ',i,', j = ',j, &
                           &           ' with prev. not conv. in cobweb sampling. Change direction.'
                      WRITE(*,'(A,I3,A,I3,A)') 'mov: point i = ',i,', j = ',j, &
                           &     ' with prev. not conv. in cobweb sampling. Change direction.'
!                      CYCLE pt_loop
                      CYCLE dir_loop
                   END IF
                ENDIF
             END IF
             ! (rho_i,rhodot_j) in the last two coordinates
             eltemp(i,j)%coord(5) = r_rdot(i,j,1)
             eltemp(i,j)%coord(6) = r_rdot(i,j,2)
             CALL fourdim_fit(m,obs,obsw,eltemp(i,j),el_mov_loc(i,j),nd,unc4_mov_loc(i,j),unc6_mov_loc(i,j), &
                  &           rms_points(i,j),delnor_cob_loc(i,j),rmsh(i,j),nused(i,j),succ4,des_mat,scal_obs)
             ! Successful 4-dim diff. cor: MOV orbit reached
             succ_cob_loc(i,j) = unc6_mov_loc(i,j)%succ.AND.succ4
             ! Find the minimum RMS among the succesful points
             IF(succ_cob_loc(i,j))THEN
                rms_min = MIN(rms_points(i,j),rms_min)
             END IF
             !*********************************************!
             ! Second order computation of the determinant !
             !*********************************************!
             IF(succ_cob_loc(i,j))THEN
                !-------------------------------------------------------!
                ! Determinant of the 4x4 matrix C_A (eq. (2.22) in [2]) !
                !-------------------------------------------------------!
                !CALL qr_inv(unc4_mov_loc(i,j)%c(1:4,1:4),out_mat(1:4,1:4),4,num_zer,det_C_A(i,j),eig_val)
!---------------------------------------!
! Computation of det(C) in a stable way !
!---------------------------------------!
!!$                DO r=1,nd
!!$                   norm_fact(r) = NORM2(unc6_mov_loc(i,j)%c(1:nd,r))
!!$                END DO
!!$                DO r=1,nd
!!$                   DO s=1,nd
!!$                      C_norm(r,s) = unc6_mov_loc(i,j)%c(r,s)/(norm_fact(r)*norm_fact(s))
!!$                   END DO
!!$                END DO
!!$                err = 1.d2*EPSILON(1.d0)
!!$                P_mat(1:nd,1:nd) = C_norm(1:nd,1:nd)
!!$                CALL tchol(P_mat(1:nd,1:nd),nd,nd,flag,err)
!!$                IF(flag.NE.0)THEN
!!$                   WRITE(*,*) 'Non-zero flag for tchol, flag,i,j', flag,i,j
!!$                END IF
!!$                DO r=1,nd
!!$                   P_mat(r,1:nd) = P_mat(r,1:nd)*norm_fact(r)
!!$                END DO
!!$                det_P(i,j)= 1.d0
!!$                DO r=1,nd
!!$                   det_P(i,j) = det_P(i,j)*P_mat(r,r)
!!$                END DO
!---------------------------------------------------------------------------------------------------
                !----------------------------------------------------!
                ! Determinant of the matrix M_mu (cf. [1], Appendix) !
                !----------------------------------------------------!
                b_rho(1:scal_obs,:) = des_mat(1:scal_obs,5:6)
                b_A(1:scal_obs,:)   = des_mat(1:scal_obs,1:4)
                b_A_rho             = MATMUL(TRANSPOSE(b_A(1:scal_obs,:)),b_rho(1:scal_obs,:))
                dAdrho              = MATMUL(unc4_mov_loc(i,j)%g(1:4,1:4),b_A_rho)
                M_A                 = MATMUL(TRANSPOSE(dAdrho),dAdrho)
                det_aux             = (M_A(1,1)+1.d0)*(M_A(2,2)+1.d0)-M_A(1,2)*M_A(2,1)
                IF(det_aux.LT.0.d0)THEN
                   WRITE(*,*) 'mov: det(M_mu)^2 is negative at sample point ',i,j
                   det_aux = 0.d0
                END IF
                jac_mu(i,j)         = SQRT(det_aux)
                IF(jac_mu(i,j).LE.jac_mu_min) WRITE(*,*) 'mov: quasi-singular matrix M_mu (AR to MOV)'
             END IF
          ELSE
             succ_cob_loc(i,j)=.FALSE.
          ENDIF
       END DO pt_loop
    END DO dir_loop
    ! Compare the minimum RMS with the RMS of the nominal, if given
    IF(rms.LT.0.d0)THEN
       rms_points(0,0) = rms_min ! No nominal, use minimum
    ELSE
       IF(rms_min.LT.rms)THEN
          WRITE(*,*) 'mov: WARNING! Found RMS lower than the nominal one. RMS = ', rms_min
       ENDIF
    ENDIF
    !***********************************!
    !  Computation of the chi function  !
    !***********************************!
    DO i=1,ndir
       DO j=1,np
          chi_cob_loc(i,j) = SQRT(ABS(rms_points(i,j)**2-rms_points(0,0)**2)*nused(i,j))
       END DO
    END DO
    ! Output of unc4_mov (if present)
    IF(PRESENT(unc4_mov)) unc4_mov = unc4_mov_loc
    ! Restore the beginning verbosity level
    verb_dif = verb_dif_old

  END SUBROUTINE mov


  !============================================================!                  
  ! WRITE_OUTPUT_FILES                                         !
  !============================================================!
  ! Write output files: web and element files (only comele).   !
  !============================================================!
  SUBROUTINE write_output_files(k_ind,root_2,astnam,rms,csinor,el_nom)
    !=======================================================================================================
    INTEGER,              INTENT(IN) :: k_ind               ! Grid number (first or second)
    DOUBLE PRECISION,     INTENT(IN) :: root_2              ! Second positive root of V(r) 
    CHARACTER*(name_len), INTENT(IN) :: astnam              ! Asteroid name
    DOUBLE PRECISION,     INTENT(IN) :: rms                 ! RMS
    DOUBLE PRECISION,     INTENT(IN) :: csinor              ! Residuals norm
    TYPE(orbit_elem),     INTENT(IN) :: el_nom              ! Nominal elements
    !=======================================================================================================
    INTEGER             :: i,j           ! Loop indices
    CHARACTER(LEN=100)  :: file_web      ! Cobweb file
    CHARACTER(LEN=100)  :: file_com      ! Cometarian file
    INTEGER             :: iuncob        ! Unit for web file
    INTEGER             :: iuncom        ! Unit for cometarian file
    INTEGER             :: le            ! File lenght
    TYPE(orbit_elem)    :: elcom         ! Cometarian elements of the sampling points
    INTEGER             :: fail_flag     ! Failure of the coordinate change
    CHARACTER(LEN=20)   :: mulname       ! Name of the sampling point 
    INTEGER             :: region        ! Number of components
    !=======================================================================================================
    !**************************!
    ! Open files .web, .comele !
    !**************************!
    IF(k_ind.EQ.1)THEN
       file_web=astnam//'_1.web'
       file_com=astnam//'_1.comele'
    ELSEIF(k_ind.EQ.2)THEN
       file_web=astnam//'.web'
       file_com=astnam//'.comele'
    END IF
    ! Header of the file .web
    CALL rmsp(file_web,le)
    CALL filopn(iuncob,file_web(1:le),'unknown')
    WRITE(iuncob,"(A)")'%  i    j   succ    rho         rho_dot        r           theta        RMS          &
         &chi      Jac(AR-MOV)    delnor    h_mag    E_Sun        E_earth'
    IF(rms.GE.0.d0)THEN
       ! If the nominal does not exist, the RMS is not defined
       IF(succ_cob(0,0)) succ_flag(0,0) = 1
       WRITE(iuncob,130) 0,0,succ_flag(0,0),r_rdot(0,0,1:2),0.d0,0.d0,csinor,chi_cob(0,0),&
            & jac_mu(0,0),delnor_cob(0,0),el_nom%h_mag,e_sun(0,0),e_earth(0,0)
    ENDIF
    ! Header of the file .comele
    CALL rmsp(file_com,le)
    CALL filopn(iuncom,file_com(1:le),'unknown')
    CALL wro1lh_matlab(iuncom,'ECLM','J2000','COM')

    !************************!
    ! Loop on the MOV orbits !
    !************************!
    DO i=1,ndir
       DO j=1,np
          WRITE(mulname,'(I3,1X,I3)') i,j
          IF(succ_cob(i,j)) succ_flag(i,j) = 1
          ! Write in the .web file
          WRITE(iuncob,130) i,j,succ_flag(i,j),r_rdot(i,j,1:2),r_cob(j),theta(i),rms_points(i,j),chi_cob(i,j),jac_mu(i,j),&
               &delnor_cob(i,j),el_mov(i,j)%h_mag,e_sun(i,j),e_earth(i,j)
130       FORMAT(I4,1X,I4,4X,I1,1X,4(F12.7,1X),1P,D12.5,1X,D12.5,1X,D12.5,0P,2X,1P,D10.3,0P,1X,F6.3,1X,&
               &1P,D12.5,0P,1X,1P,D12.5)
          !------------------------------------------------------------------!
          ! Write the points for which the MOV has been computed and chi < 5 !
          !------------------------------------------------------------------!
          IF(succ_cob(i,j) .AND. chi_cob(i,j).LE.5.d0)THEN
             ! Compute AR component
             IF(root_2.GT.0 .AND. r_rdot(i,j,1).GE.root_2)THEN
                region = 2
             ELSE
                region = 1
             END IF
             ! Write in the .comele file
             CALL coo_cha(el_mov(i,j),'COM',elcom,fail_flag)
             IF(fail_flag.GT.0)THEN
                WRITE(*,*) 'write_output_files: failed conversion to COM. Fail flag = ', fail_flag
                WRITE(*,*) '     i, j, elements ', i,j,el_mov(i,j)
                WRITE(*,*) '     elcom ',elcom
             ELSE
                CALL wro1lr_matlab(iuncom,mulname,elcom%coord,elcom%coo,elcom%t,elcom%h_mag, &
                     &             CHI=chi_cob(i,j),RMS=rms_points(i,j),SUCC=region)
             END IF
          END IF
       END DO
    END DO
    CALL filclo(iuncob,' ')
    CALL filclo(iuncom,' ')

  END SUBROUTINE write_output_files


  !====================================================================!
  ! COMPUTE_SCORE                                                      !
  !====================================================================!
  ! Compute the score of the object to be a NEO, a MBA, a distant      !
  ! object, or everything else.                                        !
  !====================================================================!
  SUBROUTINE compute_score(score_dist,score_mba,score_neo,score_scatt)
    DOUBLE PRECISION, INTENT(OUT) :: score_neo     ! Score of the object: NEO
    DOUBLE PRECISION, INTENT(OUT) :: score_mba     ! Score of the object: MBA / Trojan
    DOUBLE PRECISION, INTENT(OUT) :: score_dist    ! Score of the object: Distant
    DOUBLE PRECISION, INTENT(OUT) :: score_scatt   ! Score of the object: Scattered
    !=======================================================================================================
    DOUBLE PRECISION :: nposs                   ! Number of possible points
    DOUBLE PRECISION :: n_neo,n_mba,n_dis,n_sca ! Number of possible points per class
    LOGICAL          :: mba                     ! Flag for MBA category
    INTEGER          :: i,j                     ! Loop indexes
    TYPE(orbit_elem) :: elcom                   ! Cometary orbital elements
    INTEGER          :: fail_flag               ! Fail flag in changing coordinates
    DOUBLE PRECISION :: q_NEO, q_dist           ! Perihelion for NEO and for distant objects
    DOUBLE PRECISION :: q,ecc,sma               ! Perihelion, eccentricity, semimajor axis
    DOUBLE PRECISION :: sum                     ! Sum of the four scores
    !=======================================================================================================
    nposs = 0
    n_neo = 0
    n_mba = 0
    n_dis = 0
    n_sca = 0
    !------------------------------------------------------------!
    ! CATEGORIES DEFINITION. Each class of objects has a proper  !
    ! definition.                                                !
    !    - NEO: q < 1.3 au                                       !
    !                                                            !
    !    - MBA: non-NEO and ((1.7 < a < 4.5 and e < 0.4) or      !
    !                       (4.5 < a < 5.5 and e < 0.3))         !
    !                                                            !
    !    - DISTANT: q > 28 au                                    ! 
    !                                                            !
    !    - SCATTERED: everything else                            !
    !------------------------------------------------------------!
    q_NEO  = 1.3d0
    q_dist = 28.d0
    ! Loop on MOV points
    DO i=0,ndir
       DO j=0,np
          IF(i.EQ.0 .AND. j.NE.0) CYCLE
          IF(i.NE.0 .AND. j.EQ.0) CYCLE
          IF(i.EQ.0 .AND. j.EQ.0 .AND. (.NOT.use_nominal)) CYCLE
          IF(.NOT.succ_cob(i,j) .OR. chi_cob(i,j).GT.5.d0) CYCLE
          IF(el_mov(i,j)%h_mag.GT.hmax) CYCLE
          nposs = nposs + EXP(-chi_cob(i,j)**2/2)*jac_sigma(i,j)*jac_mu(i,j)
          CALL coo_cha(el_mov(i,j),'COM',elcom,fail_flag)
          ! Perihelion and aphelion
          q   = elcom%coord(1)
          ecc = elcom%coord(2)
          sma = q*(1-ecc)
          mba = (sma.LT.4.5d0 .AND. sma.GT.1.7d0 .AND. ecc.LT.0.4d0).OR.(sma.LT.5.5d0 .AND. sma.GT.4.5d0 .AND. ecc.LT.0.3d0)
          IF(chi_cob(i,j).LT.5.d0)THEN
             IF(q.LT.q_NEO)THEN
                n_neo = n_neo + EXP(-chi_cob(i,j)**2/2)*jac_sigma(i,j)*jac_mu(i,j)
             ELSEIF(mba)THEN
                n_mba = n_mba + EXP(-chi_cob(i,j)**2/2)*jac_sigma(i,j)*jac_mu(i,j)
             ELSEIF(q.GT.q_dist)THEN
                n_dis = n_dis + EXP(-chi_cob(i,j)**2/2)*jac_sigma(i,j)*jac_mu(i,j)
             ELSE
                n_sca = n_sca + EXP(-chi_cob(i,j)**2/2)*jac_sigma(i,j)*jac_mu(i,j)
             END IF
          END IF
       ENDDO
    ENDDO
    IF(nposs.NE.0)THEN
       score_neo   = n_neo/nposs
       score_mba   = n_mba/nposs
       score_dist  = n_dis/nposs
       score_scatt = n_sca/nposs
       !-------------------!
       ! Output adjustment !
       !-------------------!
       IF(score_neo.GT.1.d0)   score_neo   = 1.d0
       IF(score_mba.GT.1.d0)   score_mba   = 1.d0
       IF(score_dist.GT.1.d0)  score_dist  = 1.d0
       IF(score_scatt.GT.1.d0) score_scatt = 1.d0
       sum = score_neo + score_mba + score_dist + score_scatt
       IF(sum.NE.1.d0)THEN
          score_neo   = score_neo/sum
          score_mba   = score_mba/sum
          score_dist  = score_dist/sum
          score_scatt = score_scatt/sum
       END IF
    ELSE
       WRITE(iun_log,*) 'WARNING: skipped score computation because nposs = 0'
       WRITE(*,*) 'WARNING: skipped score computation because nposs = 0'
    END IF

  END SUBROUTINE compute_score


  !====================================================================!
  ! IMMEDIATE_IMP                                                      !
  !====================================================================!
  ! Search for immediate impactors, given cobweb or grid               !
  !====================================================================!
  SUBROUTINE immediate_imp(astnam,tf_search,dtclo,csinor,artyp,geoc_chi,acce_chi,chi,m,obs)
    USE propag_state, ONLY: pro_ele 
    USE tp_trace,     ONLY: wri_tppoint, dtpde, tpplane, mtp_store, njc 
    INCLUDE 'vers.h90'
    !=======================================================================================================
    CHARACTER*(name_len), INTENT(IN) :: astnam    ! Asteroid name 
    DOUBLE PRECISION,     INTENT(IN) :: tf_search ! End time for the scan
    DOUBLE PRECISION,     INTENT(IN) :: dtclo     ! Time (days) for the scan
    DOUBLE PRECISION,     INTENT(IN) :: csinor    ! Residuals norm 
    INTEGER,              INTENT(IN) :: artyp     ! Arc type
    DOUBLE PRECISION,     INTENT(IN) :: geoc_chi  ! Normalized value of the geodetic curvature
    DOUBLE PRECISION,     INTENT(IN) :: acce_chi  ! Normalized value for the along-track acceleration
    DOUBLE PRECISION,     INTENT(IN) :: chi       ! Chi-value of geodetic curvature and acceleration
    INTEGER,              INTENT(IN) :: m         ! Number of observations
    TYPE(ast_obs),        INTENT(IN) :: obs(m)    ! Observations
    !===== Impact Probability computation ==================================================================
    INTEGER            :: nposs       ! Number of points with chi <= 5, H <= H_max
    INTEGER            :: nimp        ! Number of impacting points
    LOGICAL            :: significant ! Significant flag
    TYPE(orbit_elem)   :: el1         ! Propagation of a sampling orbit
    DOUBLE PRECISION   :: pronorm     ! Normalization factor for the Impact Probability
    DOUBLE PRECISION   :: iptot       ! Total IP
    !===== Auxiliary data ==================================================================================
    DOUBLE PRECISION   :: p_earth_bound ! Probability to be an Earth-bounded object
    DOUBLE PRECISION   :: h_mag_min     ! Minimum H among MOV impactors
    DOUBLE PRECISION   :: h_mag_max     ! Maximum H among MOV impactors
    DOUBLE PRECISION   :: t_min_imp     ! Minimum impact time among MOV impactors
    DOUBLE PRECISION   :: t_max_imp     ! Maximum impact among MOV impactors
    DOUBLE PRECISION   :: hour_min      ! Hour of t_min_imp
    INTEGER            :: iday_min      ! Day of t_min_imp
    INTEGER            :: imonth_min    ! Month of t_min_imp
    INTEGER            :: iyear_min     ! Year of t_min_imp
    DOUBLE PRECISION   :: hour_max      ! Hour of t_max_imp
    INTEGER            :: iday_max      ! Day of t_max_imp
    INTEGER            :: imonth_max    ! Month of t_max_imp
    INTEGER            :: iyear_max     ! Year of t_max_imp
    INTEGER            :: imp_flag      ! Flag for impactors
    CHARACTER(LEN=16)  :: t_min_str     ! Minimum closest approach distance among MOV impactors (CAL)
    CHARACTER(LEN=16)  :: t_max_str     ! Maximum closest approach distance among MOV impactors (CAL)
    !===== Other variables =================================================================================
    CHARACTER(LEN=8)   :: date        ! Date of computation
    CHARACTER(LEN=10)  :: time        ! Time of computation
    INTEGER            :: i,j         ! Loop indexes 
    DOUBLE PRECISION   :: dt          ! Arc length (hours)
    CHARACTER(LEN=100) :: file        ! File name
    INTEGER            :: le          ! Length of file
    INTEGER            :: iunimp      ! File unit for the .mtp file (MTP analysis)
    INTEGER            :: iunrisk     ! File unit for the risk file
    CHARACTER(LEN=20)  :: mulname     ! Name of the Virtual Asteroid
    INTEGER            :: dataflag    ! Data policy flag (values from 0 to 4)
    !=======================================================================================================
    ! Select the MTP plane (not the TP)
    tpplane = .FALSE.
    ! Clean up the previous close app. records
    REWIND(iuncla) 
    !**************************************!
    !  File .mtp for the MTP data storage  !
    !**************************************!
    file=astnam//'.mtp'
    CALL rmsp(file,le)
    CALL filopn(iunimp,file(1:le),'unknown')
    WRITE(iunimp,'(A94)') '% i    j        csi          zeta      t_cla           d_cla     H         ip_fac       imp_flag'
    
    !**********************************!
    !  COMPUTE THE IMPACT PROBABILITY  !
    !**********************************!
    !******************!
    !  Inizialization  !
    !******************!
    nposs         = 0      ! Counter of "good" points
    pronorm       = 0.d0   ! Normalization factor
    nimp          = 0
    imppro        = 0.d0
    iptot         = 0.d0
    p_earth_bound = 0.d0
    h_mag_min     = 100.d0
    h_mag_max     = 0.d0
    t_min_imp     = 100000.d0
    t_max_imp     = 0.d0
    DO i=0,ndir
       DO j=0,np
          IF(i.EQ.0 .AND. j.NE.0) CYCLE                          ! Points (0,j) (with j non-zero) does not exist
          IF(i.NE.0 .AND. j.EQ.0) CYCLE                          ! Points (i,0) (with i non-zero) does not exist
          IF(i.EQ.0 .AND. j.EQ.0 .AND. (.NOT.use_nominal)) CYCLE ! Skip the nominal if it is not available
          IF(.NOT.succ_cob(i,j) .OR. chi_cob(i,j).GT.5.d0) CYCLE ! Consider only the successful points with chi <= 5
          IF(el_mov(i,j)%h_mag.LE.hmax)THEN                      ! Consider the points below the shooting star limit
             !---------------------------------------------------!
             ! Normalization factor for IP (cf. [2], eq. (2.22)) !
             !---------------------------------------------------!
             nposs   = nposs + 1
             pronorm = pronorm + EXP(-chi_cob(i,j)**2/2.d0)*jac_mu(i,j)*jac_sigma(i,j)
             !----------------------!
             ! Earth-bounded orbits !
             !----------------------!
             imp_flag = 0
             IF(e_earth(i,j).LE.0.d0)THEN
                ! Take into account the VA with negative two-body energy with the Earth
                p_earth_bound = p_earth_bound + EXP(-chi_cob(i,j)**2/2.d0)*jac_mu(i,j)*jac_sigma(i,j)
                IF(.NOT.propag_geoc_orbit) CYCLE
             END IF
             WRITE(mulname,590) i,j
590          FORMAT('cob_',I2,'_',I2)
             CALL rmsp(mulname,le)
             WRITE(iuncla,*) mulname
             CALL cov_not_av ! Propagate without covariance
             njc = 0         ! Number of minima during a close approach
             WRITE(iun_log,*) 'MOV orbit ',i,'  ',j,' propagation'
             CALL pro_ele(el_mov(i,j),tf_search,el1) 
             !****************!
             !  MTP analysis  !
             !****************!
             IF(njc.GT.0)THEN
                !-------------------------------------------------------!
                ! WARNING: we should handle double minima with a DO on  !
                ! jc and selecting the minimum mtp_store(jc)%r          !
                !-------------------------------------------------------!
                IF(njc.GT.1)THEN
                   WRITE(iun_log,*) '   multiple minima, njc = ',njc
                END IF
                IF(mtp_store(1)%iplam.EQ.3)THEN
                   csi(i,j)  = mtp_store(1)%tp_coord(1)/reau
                   zeta(i,j) = mtp_store(1)%tp_coord(2)/reau
                   tcla(i,j) = mtp_store(1)%tcla
                   dcla(i,j) = mtp_store(1)%rcla
                   IF(csi(i,j)**2+zeta(i,j)**2.LT.1.d0)THEN
                      imp_flag = 1
                      nimp = nimp+1 ! Counter for impacting points
                      !------------------------------------------------!
                      ! Contribution from the VI (cf. [2], eq. (2.22)) !
                      !------------------------------------------------!
                      imppro(i,j) = EXP(-chi_cob(i,j)**2/2.d0)*jac_mu(i,j)*jac_sigma(i,j)
                      iptot       = iptot + imppro(i,j)                   
                      !************************************************!
                      ! Auxiliary data: H and impact time of impactors !
                      !************************************************!
                      h_mag_max = MAX(h_mag_max,el_mov(i,j)%h_mag)
                      h_mag_min = MIN(h_mag_min,el_mov(i,j)%h_mag)
                      t_max_imp = MAX(t_min_imp,tcla(i,j))
                      t_min_imp = MIN(t_min_imp,tcla(i,j))
                   END IF
                   WRITE(iunimp,591) i, j, csi(i,j), zeta(i,j), tcla(i,j), dcla(i,j), el_mov(i,j)%h_mag, imppro(i,j), imp_flag
591                FORMAT(I4,1X,I4,4(1X,F12.6),2X,F5.2,2X,E16.10,5X,I1)         
                END IF
             END IF
          END IF
       END DO
    END DO
    CALL filclo(iunimp,' ')
    ! Normalization of probabilities
    IF((use_nominal .AND. nposs.GE.2) .OR. (.NOT.use_nominal .AND. nposs.GE.1))THEN
       p_earth_bound = p_earth_bound/pronorm
       iptot         = iptot/pronorm
    END IF

    !***********************!
    !  WRITE THE RISK FILE  !
    !***********************!
    !-------------------------------------------------------------!
    ! We consider significant a case that has more than 3         !
    ! observations, and an arc of more than 0.5 hours (i.e., 30   !
    ! minutes), or a case which has been selected through the     !
    ! select_risk flag, which is true when geoc_chi > 1.          !
    !-------------------------------------------------------------!  
    dt          = (obs(m)%time_tdt-obs(1)%time_tdt)*24.d0       ! Arc duration (hours)
    significant = (m.GT.2 .AND. dt.GT.0.5d0) .OR. (select_risk) ! Significant flag
    !------------------------------------------------------------------!
    ! If significant is false we open a .nosig file. Otherwise, if     !
    ! nimp = 0 (there are no impacting points) we open a .norisk file, !
    ! while if nimp > 0 (and thus IP > 0) we open a .risk file.        !
    !------------------------------------------------------------------!
    IF(.NOT.significant)THEN
       file=astnam//'.nosig'
    ELSEIF(nimp.GT.0)THEN
       file=astnam//'.risk'
    ELSE
       file=astnam//'.norisk' 
    ENDIF
    !-----------------------------------------------------------------------!
    ! Data policy flag assignment. We give the following data policy flag:  !
    !                                                                       !
    ! - 0 for IP <= 1E-6                                                    !
    !                                                                       !
    ! - 1 for 1E-6 < IP <= 1E-3                                             !
    !                                                                       !
    ! - 2 for 1E-3 < IP <= 1E-2                                             !
    !                                                                       !
    ! - 3 for IP > 1E-2, and when a nominal does not exists or exists       !
    !   but curvature chi^2 <= 10                                           !
    !                                                                       !
    ! - 4 for IP > 1E-2, when a nominal exists and has curvature chi^2 > 10 !
    !-----------------------------------------------------------------------!
    IF(iptot.LE.1.d-6)THEN
       dataflag = 0
    ELSEIF(iptot.LE.1.d-3)THEN
       dataflag = 1
    ELSEIF(iptot.LE.1.d-2)THEN
       dataflag = 2
    ELSEIF(use_nominal .AND. chi**2.GT.10)THEN
       dataflag = 4
    ELSE
       dataflag = 3
    ENDIF
    !**************!
    !  File .risk  !
    !**************!
    CALL rmsp(file,le)
    CALL filopn(iunrisk,file(1:le),'unknown')
    WRITE(iunrisk,201) 'NEOCP     # obs      dt    Arc        Significance of curvature       ', &
         &             'RMS      Total    MOV    Imp      IP      Data '
    WRITE(iunrisk,201) 'Name                 (h)   type    Geodesic     Accel.     Overall    ', &
         &             '         # pts    pts    pts              flag '
    WRITE(iunrisk,201) '----------------------------------------------------------------------', &
         &             '-----------------------------------------------'
201 FORMAT(A70,A47) 
    IF(use_nominal)THEN
       !***************!
       !  Cobweb case  !
       !***************!
       ! Content (IP in exp. format if < 1E-3, otherwise in floating point format)
       IF(iptot.GT.0.d0 .AND. iptot.LT.1.d-3)THEN
          WRITE(iunrisk,555) astnam,m,dt,artyp,geoc_chi,acce_chi,chi,csinor, &
               &             ndir*np,nposs,nimp,iptot,dataflag
       ELSE
          WRITE(iunrisk,556) astnam,m,dt,artyp,geoc_chi,acce_chi,chi,csinor, &
               &             ndir*np,nposs,nimp,iptot,dataflag
       END IF
555    FORMAT(A9,2X,I3,3X,F8.2,2X,I2,4X,F9.2,2X,F9.2,3X,F9.2,3X,F7.3,2X,I7,2X,I5,2X,I5,2X,1P,E10.3,4X,I1) 
556    FORMAT(A9,2X,I3,3X,F8.2,2X,I2,4X,F9.2,2X,F9.2,3X,F9.2,3X,F7.3,2X,I7,2X,I5,2X,I5,5X,F5.3,5X,I1)   
    ELSE
       !********************************!
       !  Grid case (no RMS avaliable)  !
       !********************************!
       ! Content (IP in exp. format if < 1E-3, otherwise in floating point format)
       IF(iptot.GT.0.d0 .AND. iptot.LT.1.d-3)THEN
          WRITE(iunrisk,557) astnam,m,dt,artyp,geoc_chi,acce_chi,chi,        &
               &             ndir*np,nposs,nimp,iptot,dataflag
       ELSE
          WRITE(iunrisk,558) astnam,m,dt,artyp,geoc_chi,acce_chi,chi,        &
               &             ndir*np,nposs,nimp,iptot,dataflag
       END IF

557    FORMAT(A9,2X,I3,3X,F8.2,2X,I2,4X,F9.2,2X,F9.2,3X,F9.2,6X,'-',5X,I7,2X,I5,2X,I5,2X,1P,E10.3,4X,I1) 
558    FORMAT(A9,2X,I3,3X,F8.2,2X,I2,4X,F9.2,2X,F9.2,3X,F9.2,6X,'-',5X,I7,2X,I5,2X,I5,5X,F5.3,5X,I1)   
    ENDIF
    !***************!
    !  Write notes  !
    !***************!
    !-------------------------!
    ! Elongation from the Sun !
    !-------------------------!
    WRITE(iunrisk,595) elong*degrad
    WRITE(iunrisk,*)
    !--------------------------------------------------------!
    ! Absolute magnitude, closest approach time and distance !
    !--------------------------------------------------------!
    IF(nimp.GT.0)THEN
       CALL convert_mjd_cal(t_min_imp,t_min_str)
       CALL convert_mjd_cal(t_max_imp,t_max_str)
       WRITE(iunrisk,596) h_mag_min,h_mag_max
       WRITE(iunrisk,597) t_min_str,t_max_str
    ENDIF
    !-------------------------------------------!
    ! Probability to be an Earth-bounded object !
    !-------------------------------------------!
    WRITE(iunrisk,599) p_earth_bound
    !---------------------!
    ! Shooting star limit !
    !---------------------!
    WRITE(iunrisk,600) hmax
    WRITE(iunrisk,*)
    !-------------------------------------!
    ! OrbFit version, date of computation !
    !-------------------------------------!
    CALL date_and_time(date,time)
    WRITE(iunrisk,121) pvers, date, time
    ! Reason for which the case is not significant
    IF(.NOT.significant)THEN
       WRITE(iunrisk,*)
       IF(m.LE.2)THEN
          WRITE(iunrisk,122) m
       ELSEIF(dt.LE.0.5d0)THEN
          WRITE(iunrisk,123) dt
       ENDIF
    ENDIF
595 FORMAT(/'Elongation from the Sun = ',F5.1)
596 FORMAT(/'Absolute magnitude of impactors: H_min = ',F5.2,' and H_max = ',F5.2)
597 FORMAT(/'Closest approach time of impactors: t_min = ',A16,' TDB and t_max = ',A16,' TDB')
599 FORMAT(/'Probability to be an Earth-bounded object = ', F5.3)
600 FORMAT(/'Shooting star limit: ',F5.2, ' magnitudes')
121 FORMAT(/'OrbFit software version ',A17,' Date of computation = ',A10,1X,A10,' CET')
122 FORMAT(/'This case is not significant because the number of observations is ',I3)
123 FORMAT(/'This case is not significant because dt is ',F8.2,' hours')
    CALL filclo(iunrisk,' ')

  END SUBROUTINE immediate_imp
  

  !====================================================================!                  
  ! MC_RES_ATT                                                         !
  !====================================================================!
  ! Search for immediate impactors with MC method on residuals         !
  ! assuming a nominal solution is available                           !
  !====================================================================!
  SUBROUTINE mc_res_att(astnac,nmc,tf_search,m,obs,obsw,el0,obscodn,nd,nimp,isucc)
    USE tp_trace,     ONLY: tpplane, mtp_store, njc
    USE close_app,    ONLY: fix_mole, kill_propag
    USE propag_state, ONLY: pro_ele
    INCLUDE 'parobx.h90'
    !=======================================================================================================
    CHARACTER*(name_len), INTENT(IN)  :: astnac     ! Asteroid name
    INTEGER,              INTENT(IN)  :: nmc        ! Number of MC sample points
    DOUBLE PRECISION,     INTENT(IN)  :: tf_search  ! Final time of propagation
    INTEGER,              INTENT(IN)  :: obscodn    ! Observatory code
    INTEGER,              INTENT(IN)  :: m          ! Number of observations
    TYPE(ast_obs),        INTENT(IN)  :: obs(nobx)  ! Observations
    TYPE(ast_wbsr),       INTENT(IN)  :: obsw(nobx) ! Residuals and weigths
    TYPE(orbit_elem),     INTENT(IN)  :: el0        ! Nominal elements
    INTEGER,              INTENT(IN)  :: nd         ! Dimension of the par space
    INTEGER,              INTENT(OUT) :: nimp       ! Number of impacts
    INTEGER,              INTENT(OUT) :: isucc      ! Number of convergent differential corrections
    !=======================================================================================================
    CHARACTER(LEN=30) :: file        ! File name
    !===== MC computation ==================================================================================
    INTEGER           :: k           ! Index for the MC points
    INTEGER           :: kk          ! For one istance of noise
    DOUBLE PRECISION  :: rvnorm      ! Random noise
    TYPE(ast_obs)     :: obs1(nobx)  ! Perturbed observations
    TYPE(ast_wbsr)    :: obsw1(nobx) ! Perturbed weights
    INTEGER           :: icor(6)     ! Parameters to be corrected
    TYPE(orbit_elem)  :: elk         ! Keplerian elements (after diff_cor)
    TYPE(orb_uncert)  :: unck        ! Uncertainty matrices for elk
    TYPE(orbit_elem)  :: elatt       ! Attributable orbital elements
    TYPE(orbit_elem)  :: el1         ! Propagated elements
    LOGICAL           :: succ        ! Successful differential corrections
    DOUBLE PRECISION  :: csinok      ! Norm of the residuals
    DOUBLE PRECISION  :: delnok      ! Norm of the last correction
    DOUBLE PRECISION  :: csi         ! Xi on the TP
    DOUBLE PRECISION  :: zeta        ! Zeta on the TP
    !===== Other variables =================================================================================
    INTEGER           :: fail_flag   ! Failure flag for the coordinate change
    INTEGER           :: le          ! File length
    INTEGER           :: iuncob      ! Output file for the MC points
    !=======================================================================================================
    IF(nd.NE.6)THEN
       WRITE(*,*) 'mc_res_att: not ready for non-grav, nd = ',nd
       STOP
    END IF
    !******************!
    !  Inizialization  !
    !******************!
    verb_dif=1
    tpplane=.FALSE. ! Use MTP plane, not TP
    !***************************!
    !  File .mc for the output  !
    !***************************!
    file=astnac//'.mc'
    CALL rmsp(file,le)
    CALL filopn(iuncob,file(1:le),'unknown')
    
    !***************************************!
    !  LOOP ON INSTANCES OF GAUSSIAN NOISE  !
    !***************************************!
    nimp=0  ! Impact counter
    isucc=0 ! Counter of nominal solutions found
    DO k=1,nmc
       !********************************!
       !  Create one instance of noise  !
       !********************************!
       DO kk=1,m 
          IF(obs(kk)%type.EQ.'O')THEN
             obs1(kk)=obs(kk)
             obs1(kk)%coord(1)=obs1(kk)%coord(1)+rvnorm()*obsw(kk)%rms_coord(1)/COS(obs1(kk)%coord(2))
             obs1(kk)%coord(2)=obs1(kk)%coord(2)+rvnorm()*obsw(kk)%rms_coord(2)
             obsw1(kk)=obsw(kk)
          ENDIF
       ENDDO
       !****************************!
       !  Compute nominal solution  !
       !****************************!
       icor=1    
       CALL diff_cor(m,obs1,obsw1,el0,icor,-1,elk,unck,csinok,delnok,succ,nd)
       ! Check convergence
       IF(succ)THEN
          isucc=isucc+1
          ! First output the data point on the (rho,rho_dot) plane
          CALL coo_cha(elk,'ATT',elatt,fail_flag,OBSCODE=obscodn)
          WRITE(iuncob,130) k,isucc,elatt%coord(5:6),0.d0,0.d0,csinok,0.d0,0.d0,0.d0, &
               &            1,delnok,elatt%h_mag,0.d0
130       FORMAT(I4,1X,I4,1X,4(F12.7,1X),1P,D12.5,1X,D12.5,1X,D12.5,1X,D12.5,0P,1X,I1,1X,1P,D10.3,0P,1X,F6.3,1X,1P,D12.5)
          !******************************!
          !  Propagate until final time  !
          !******************************!
          CALL cov_not_av
          njc=0
          CALL pro_ele(elk,tf_search,el1)
          kill_propag=.FALSE.
          !******************************************!
          !  Read close approach data, if available  !
          !******************************************!
          IF(njc.GT.0)THEN
             IF(mtp_store(1)%iplam.EQ.3)THEN
                csi=mtp_store(1)%tp_coord(1)/reau
                zeta=mtp_store(1)%tp_coord(2)/reau
                ! Impactor?
                IF(csi**2+zeta**2.LT.1.d0)THEN
                   nimp=nimp+1
                ENDIF
             ENDIF
          ENDIF
       ENDIF
       WRITE(*,*) k, isucc, nimp, csinok
    ENDDO
    CALL filclo(iuncob,' ')
    verb_dif=20
  END SUBROUTINE mc_res_att

END MODULE cobweb
