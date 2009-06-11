MODULE cat_debias
  USE healpix_types
  IMPLICIT NONE
  PRIVATE 
! common data: biases per catalog
  DOUBLE PRECISION         :: catrecord(49152,26)
  INTEGER(KIND=i4b), PARAMETER :: ns_max=8192 ! 2^13 : largest nside available

  !initialise array x2pix, y2pix and pix2x, pix2y used in several routines
  integer(KIND=i4b), save, dimension(128) :: x2pix=0,y2pix=0

  integer(KIND=i4b), save, dimension(0:1023) :: pix2x=0, pix2y=0

! public routines
  PUBLIC cat_find,cat_bias
  
  CONTAINS

  SUBROUTINE cat_bias_init
    INTEGER  :: error, i,iun
    CALL filopl(iun,'xcat.bias')
    READ(iun,*)
704 FORMAT(7X,F10.6,25(F16.6))
    DO i=1,49152
       READ(iun,704,IOSTAT=error) catrecord(i,:)
       IF (error/=0)THEN
         WRITE(*,*)' cat_bias_init: error in reading, record=',i 
         WRITE(*,*) catrecord(i,:)
       ENDIF
    ENDDO
   END SUBROUTINE cat_bias_init

  !------------------------------------------------------------------
  ! SUBROUTINE cat_find
  ! The purpose of this routine is to provide the astrometric catalog
  ! code starting from the MPC observations flag.
  ! Special cases are handled for observatory codes 691, 703 and 704
  ! Creation on 27 May 2009 by F. Bernardi
  !-----------------------------------------------------------------
  SUBROUTINE cat_find(mpcflag,t_cat,obscode,catcode)
    IMPLICIT NONE
    CHARACTER(LEN=2),           INTENT(OUT)   :: catcode
    CHARACTER(LEN=2),           INTENT(INOUT) :: mpcflag
    CHARACTER(LEN=3),           INTENT(IN)    :: obscode
    DOUBLE PRECISION,           INTENT(IN)    :: t_cat
    IF (mpcflag.EQ.'  ') THEN
       !----------------
       ! Special cases
       !----------------
       IF (obscode.EQ.'691') THEN
          IF(t_cat.GE.48499.AND.t_cat.LE.51449) THEN
             catcode='G1'
             mpcflag='sh'
          ELSEIF(t_cat.GT.51449.AND.t_cat.LE.51899) THEN
             catcode='A1'
             mpcflag='sa'
          ELSEIF(t_cat.GT.51899.AND.t_cat.LE.54095) THEN 
             catcode='A2'
             mpcflag='sc'
          ENDIF
       ELSEIF (obscode.EQ.'703') THEN
          IF(t_cat.LE.53371) THEN
             catcode='A2'
             mpcflag='sc'
          ENDIF
       ELSEIF (obscode.EQ.'704') THEN
          IF(t_cat.LE.51544) THEN
             catcode='A1'
             mpcflag='sa'
          ELSEIF(t_cat.GT.51544) THEN 
             catcode='A2'
             mpcflag='sc'
          ENDIF
       ELSE
          catcode='UN'                                        ! Unknown
       ENDIF
    ELSEIF (mpcflag.EQ.' a'.OR.mpcflag.EQ.' b') THEN 
       catcode='A1'                                        ! USNO-A1
    ELSEIF (mpcflag.EQ.' c'.OR.mpcflag.EQ.' d') THEN 
       catcode='A2'                                        ! USNO-A2
    ELSEIF (mpcflag.EQ.' e'.OR.mpcflag.EQ.' q'.OR.mpcflag.EQ.' r'.OR.mpcflag.EQ.' t'.OR.mpcflag.EQ.' u') THEN
       catcode='UC'                                        ! UCAC
    ELSEIF (mpcflag.EQ.' f'.OR.mpcflag.EQ.' g') THEN 
       catcode='TY'                                        ! Tycho
    ELSEIF (mpcflag.EQ.' h'.OR.mpcflag.EQ.' i'.OR.mpcflag.EQ.' j'.OR.mpcflag.EQ.' z') THEN 
       catcode='G1'                                        ! GSC-1
    ELSEIF (mpcflag.EQ.' k') THEN 
       catcode='G2'                                        ! GSC-2
    ELSEIF (mpcflag.EQ.' l') THEN 
       catcode='AC'                                        ! ACT
    ELSEIF (mpcflag.EQ.' m') THEN 
       catcode='GA'                                        ! GSC-ACT
    ELSEIF (mpcflag.EQ.' n') THEN 
       catcode='TR'                                        ! TRC
    ELSEIF (mpcflag.EQ.' o'.OR.mpcflag.EQ.' s') THEN 
       catcode='B1'                                        ! USNO-B1
    ELSEIF (mpcflag.EQ.' p') THEN 
       catcode='PP'                                        ! PPM
    ELSEIF (mpcflag.EQ.' v') THEN 
       catcode='NO'                                         ! NOMAD
    ELSEIF (mpcflag.EQ.' w') THEN 
       catcode='CM'                                        ! CMC
    ELSEIF (mpcflag.EQ.' x') THEN 
       catcode='HI'                                        ! Hip-2
    ELSE
       catcode='UN'                                        ! Unknown
    ENDIF
  END SUBROUTINE cat_find

  !------------------------------------------------------------------
  ! SUBROUTINE cat_bias
  ! The purpose of this routine is to provide the astrometric catalog
  ! biases from the xcat.bias file. Biases provided by S. Chesley and
  ! Baer.
  ! Outputs are bias in RA and DEC, given the coordinates of the 
  ! observations
  ! Creation on 21 May 2009 by F. Bernardi
  !-----------------------------------------------------------------

  SUBROUTINE cat_bias(phi,del,catcode,bias_RA,bias_DEC,biasflag)
    USE healpix_types
    DOUBLE PRECISION,      INTENT(OUT) :: bias_DEC
    DOUBLE PRECISION,      INTENT(OUT) :: bias_RA
    LOGICAL,               INTENT(OUT) :: biasflag        
    CHARACTER*2,           INTENT(IN)  :: catcode               ! this is our internal cat coding
!    DOUBLE PRECISION,      INTENT(IN)  :: catrecord(49152,26)   ! This array contains the cat biases
    DOUBLE PRECISION,      INTENT(IN)  :: del           ! DEC of the observations in radians
    DOUBLE PRECISION                   :: dpig=6.28318530717958648d0
    INTEGER(KIND=4)                    :: ipring
    INTEGER(KIND=4),       PARAMETER   :: nside=64      ! This is the SIDE number for the HEALpix pixellization
    DOUBLE PRECISION,      INTENT(IN)  :: phi           ! RA of the observations in radians
    DOUBLE PRECISION                   :: theta         ! colatitude in equatorial coo.
    LOGICAL,                      SAVE :: first=.true.
    IF(first)THEN
       CALL cat_bias_init
       first=.false.
    ENDIF
     theta=dpig/4-del         
    CALL ang2pix_ring(nside, theta, phi, ipring)
    IF(catcode.EQ.'TY') THEN                         ! Case for Tycho catalog
       bias_RA=catrecord((ipring+1),3)
       bias_DEC=catrecord((ipring+1),(3+2))
       biasflag=.true.
    ELSEIF(catcode.EQ.'UC') THEN                     ! Case for UCAC catalog
       bias_RA=catrecord((ipring+1),7)
       bias_DEC=catrecord((ipring+1),(7+2))
       biasflag=.true.
    ELSEIF(catcode.EQ.'B1') THEN                     ! Case for USNO-B1 catalog
       bias_RA=catrecord((ipring+1),11)
       bias_DEC=catrecord((ipring+1),(11+2))
       biasflag=.true.
    ELSEIF(catcode.EQ.'A1') THEN                     ! Case for USNO-A1 catalog
       bias_RA=catrecord((ipring+1),15)
       bias_DEC=catrecord((ipring+1),(15+2))
       biasflag=.true.
    ELSEIF(catcode.EQ.'A2') THEN                     ! Case for USNO-A2 catalog
       bias_RA=catrecord((ipring+1),19)
       bias_DEC=catrecord((ipring+1),(19+2))
       biasflag=.true.
    ELSEIF(catcode.EQ.'G2') THEN                     ! Case for GSC-2 catalog
       bias_RA=catrecord((ipring+1),23)
       bias_DEC=catrecord((ipring+1),(23+2))
       biasflag=.true.
    ELSE
       bias_RA=0.d0
       bias_DEC=0.d0
       biasflag=.false.
    ENDIF
    ! Handle cases where biases are not provided
    IF(biasflag.eq..true.)THEN
       IF(bias_RA.gt.9000.or.bias_DEC.gt.9000)THEN
          bias_RA=0.d0
          bias_DEC=0.d0
          biasflag=.false. 
       ENDIF
    ENDIF
  END SUBROUTINE cat_bias

  !=======================================================================
  subroutine ang2pix_ring(nside, theta, phi, ipix)
    !=======================================================================
    !     renders the pixel number ipix (RING scheme) for a pixel which contains
    !     a point on a sphere at coordinates theta and phi, given the map
    !     resolution parameter nside
    !=======================================================================
    INTEGER(KIND=i4b), INTENT(IN) :: nside
    INTEGER(KIND=i4b), INTENT(OUT) :: ipix
    REAL (KIND=dp), INTENT(IN) ::  theta, phi

    INTEGER(KIND=i4b) ::  nl4, jp, jm
    REAL (KIND=dp)  ::  z, za, tt, tp, tmp, temp1, temp2
    INTEGER(KIND=i4b) ::  ir, ip, kshift

    !-----------------------------------------------------------------------
    interface fatal_error
       module procedure fatal_error_womsg, fatal_error_msg
    end interface
    if (nside<1 .or. nside>ns_max) call fatal_error ("nside out of range")
    if (theta<0.0_dp .or. theta>pi)  then
       print*,"ANG2PIX_RING: theta : ",theta," is out of range [0, Pi]"
       call fatal_error
    endif

    z = COS(theta)
    za = ABS(z)
    tt = MODULO( phi, twopi) / halfpi  ! in [0,4)


    if ( za <= twothird ) then ! Equatorial region ------------------
       temp1 = nside*(.5_dp+tt)
       temp2 = nside*.75_dp*z
       jp = int(temp1-temp2) ! index of  ascending edge line
       jm = int(temp1+temp2) ! index of descending edge line

       ir = nside + 1 + jp - jm ! in {1,2n+1} (ring number counted from z=2/3)
       kshift = 1 - modulo(ir,2) ! kshift=1 if ir even, 0 otherwise

       nl4 = 4*nside
       ip = INT( ( jp+jm - nside + kshift + 1 ) / 2 ) ! in {0,4n-1}
       if (ip >= nl4) ip = ip - nl4

       ipix = 2*nside*(nside-1) + nl4*(ir-1) + ip

    else ! North & South polar caps -----------------------------

       tp = tt - INT(tt)      !MODULO(tt,1.0_dp)
       tmp = nside * SQRT( 3.0_dp*(1.0_dp - za) )

       jp = INT(tp          * tmp ) ! increasing edge line index
       jm = INT((1.0_dp - tp) * tmp ) ! decreasing edge line index

       ir = jp + jm + 1        ! ring number counted from the closest pole
       ip = INT( tt * ir )     ! in {0,4*ir-1}
       if (ip >= 4*ir) ip = ip - 4*ir

       if (z>0._dp) then
          ipix = 2*ir*(ir-1) + ip
       else
          ipix = 12*nside**2 - 2*ir*(ir+1) + ip
       endif

    endif

    return
  end subroutine ang2pix_ring

!-----------------------------------------------------
!  subroutine fatal_error (msg)
!    character(len=*), intent(in), optional :: msg
!
!    if (present(msg)) then
!       print *,'Fatal error: ', trim(msg)
!    else
!       print *,'Fatal error'
!    endif
!    call exit_with_status(1)
!  end subroutine fatal_error

  subroutine fatal_error_msg (msg)
    character(len=*), intent(in) :: msg
       print *,'Fatal error: ', trim(msg)
    call exit_with_status(1)
  end subroutine fatal_error_msg

  subroutine fatal_error_womsg
      print *,'Fatal error'
    call exit_with_status(1)
  end subroutine fatal_error_womsg

  subroutine exit_with_status (code, msg)
    ! ===========================================================
    integer, intent(in) :: code
    character (len=*), intent(in), optional :: msg
    ! ===========================================================
!  USE f90_unix, ONLY : iargc, getarg, exit
  
    if (present(msg)) print *,trim(msg)
    print *,'program exits with exit code ', code
    
!#if (defined (RS6000))
!    call exit_ (code)
!#else
!    call exit (code)
!#endif
  end subroutine exit_with_status


END MODULE cat_debias
