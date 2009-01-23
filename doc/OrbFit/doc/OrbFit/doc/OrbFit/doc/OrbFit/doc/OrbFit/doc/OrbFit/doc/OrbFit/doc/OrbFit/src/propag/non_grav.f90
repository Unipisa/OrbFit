! ==========MODULE non_grav============
! CONTAINS
! comet_non_grav

MODULE non_grav
USE fund_const
IMPLICIT NONE
PRIVATE

PUBLIC comet_non_grav_symm, comet_non_grav_asymm
! private common data
      INTEGER inongrav ! flag to select non grav, 0=no, 1=sym, 2=asym

! non gravitational model terms A1, A2,A3, and Delta T delay from perihelion passage
      DOUBLE PRECISION A1ng,A2ng,A3ng,DTdelay ! Units of A1,A2 and A3 are in AU/d**2
      PUBLIC inongrav, A1ng, A2ng,A3ng,DTdelay
      
CONTAINS 

! ===========================================================
! COMET NON GRAV. TERMS added 29/10/2008 - Symmetric Model
! Marsden and Sekanina 1973 - "Comets and nongravitational
! forces. V.

      SUBROUTINE comet_non_grav_symm(x,v,nongrav)
!  interface: INPUT
        DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: x,v ! position and velocity, heliocentric
!  interface: OUTPUT
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(3) :: nongrav
! end interface
!        DOUBLE PRECISION :: vvec(3),vv2, rr, vsize, prscal, factor, conv
        DOUBLE PRECISION :: rr,rr2,rv,trans_size,trans(3),xu(3),tu(3),hu(3)
        DOUBLE PRECISION :: vsize, prscal
        DOUBLE PRECISION :: al1,rr0,mesp,nesp,kesp,grfac
! =========================================================
         rr=vsize(x)
         rr0=2.808         ! in AU - scale distance for water ice(frost line)
         al1=0.111262      ! normalizing constant
         mesp=2.15         ! esponent m of the formula
         nesp=5.093        ! esponent n of the formula
         kesp=4.6142       ! esponent k of the formula
         grfac=al1*((rr/rr0)**(-mesp))*(1+(rr/rr0)**nesp)**(-kesp)
         xu=x/rr
         rr2=prscal(x,x)
         rv =prscal(x,v)
         CALL prvec(x,v,hu)
         hu=hu/vsize(hu)
         trans=v-(x*rv/rr2)
         trans_size=vsize(trans)
         tu=trans/trans_size
         nongrav=grfac*(a1ng*xu+a2ng*tu+a3ng*hu)
      END SUBROUTINE comet_non_grav_symm

! ===========================================================
! COMET NON GRAV. TERMS added 07/11/2008 - Asymmetric Model
! Yeomans and Chodas 1989 - "An Asymmetric Outgassing Model
! for Cometary Nongravitational Accelerations"

      SUBROUTINE comet_non_grav_asymm(x,v,nongrav)
        USE ever_pitkin
!  interface: INPUT
        DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: x,v ! position and velocity, heliocentric
!  interface: OUTPUT
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(3) :: nongrav
! end interface
!        DOUBLE PRECISION :: vvec(3),vv2, rr, vsize, prscal, factor, conv
        DOUBLE PRECISION :: rr,rr2,rv,trans_size,trans(3),xu(3),tu(3),hu(3)
        DOUBLE PRECISION :: vsize, prscal
        DOUBLE PRECISION :: al1,rr0,mesp,nesp,kesp,grfac
        DOUBLE PRECISION :: xdelay(6), t0
! =========================================================
        t0=0
        CALL fser_propag(x,v,t0,-DTdelay,gms,xdelay(1:3),xdelay(4:6))
        rr=vsize(xdelay(1:3))
        rr0=2.808         ! in AU - scale distance for water ice(frost line)
        al1=0.111262      ! normalizing constant
        mesp=2.15         ! esponent m of the formula
        nesp=5.093        ! esponent n of the formula
        kesp=4.6142       ! esponent k of the formula
        grfac=al1*((rr/rr0)**(-mesp))*(1+(rr/rr0)**nesp)**(-kesp)
        xu=x/rr
        rr2=prscal(x,x)
        rv =prscal(x,v)
        CALL prvec(x,v,hu)
        hu=hu/vsize(hu)
        trans=v-(x*rv/rr2)
        trans_size=vsize(trans)
        tu=trans/trans_size
        nongrav=grfac*(a1ng*xu+a2ng*tu+a3ng*hu)
      END SUBROUTINE comet_non_grav_asymm
END MODULE non_grav
