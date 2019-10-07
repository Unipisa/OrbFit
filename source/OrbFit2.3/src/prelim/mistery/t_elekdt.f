*
* TEST of subroutine elekdt
*
      PROGRAM prova
      IMPLICIT NONE

      INTEGER ll,i
      DOUBLE PRECISION gk,g,gm,rm,xl,dt,dt1
      DOUBLE PRECISION elem0(6),elem(6),elemi(6),xv(6),xva(6)
      CHARACTER*4 eltype,eltypi

* Gauss' gravitational constant
      gk=0.01720209895d0
      g=gk**2

      OPEN(1,FILE='t_elekdt.in',STATUS='OLD')
      READ(1,100) eltype
      READ(1,*) elem0
      READ(1,*) rm
      READ(1,*) dt
      READ(1,*) xl,ll
 100  FORMAT(A)

      gm=g*(1.D0+1/rm)
      DO 1 i=1,6
      elem(i)=elem0(i)
 1    CONTINUE

      CALL forcin(1,0,1.D0,rm)
      CALL ekcc1(elem0,eltype,xv,gm,0.d0)
      dt1=dt
      CALL ra15(xv(1),xv(4),dt1,xl,ll,3,-2,15)
      CALL ccek1(elemi,eltypi,xv,gm)

      CALL elekdt(elem,eltype,gm,dt)
      CALL ekcc1(elem,eltype,xva,gm,0.d0)

      WRITE(*,201) eltype,elem0
      WRITE(*,201) eltype,elem
      WRITE(*,201) eltypi,elemi

      WRITE(*,*)

      WRITE(*,200) xv
      WRITE(*,200) xva

 200  FORMAT(1P,6E20.12)
 201  FORMAT(A,1P,6E20.12)

      END
