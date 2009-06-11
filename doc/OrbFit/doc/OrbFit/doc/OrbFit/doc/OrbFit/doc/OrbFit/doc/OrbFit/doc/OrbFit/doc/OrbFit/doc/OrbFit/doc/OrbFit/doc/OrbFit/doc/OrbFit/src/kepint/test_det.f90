program test_det
  implicit none
  INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
  integer:: h,k
  COMPLEX(KIND=qkind) :: ev_Sylv_j(22,22),det_Sylv_j,testc
  REAL(KIND=qkind) :: theta

!  write(*,*) qcos(3.14q0)
!  read(*,*)

  ev_Sylv_j=0.q0
  do h=1,22
     do k=1,22
        if(h.eq.k)then
           ev_Sylv_j(h,k)=(1.q0,0.q0)
        endif
!        write(*,*)ev_Sylv_j(h,k)
     enddo
  enddo
!     CALL cdetcomp22_QP(ev_Sylv_j,det_Sylv_j)
!     write(*,*) 'deteminant=',det_Sylv_j

     testc=QCMPLX(1.01234567890123456789q0,2.01234567890123456789q0)
!     write(*,*)'cnum=',testc
!     write(*,100)'real=',qreal(testc)
!     write(*,100)'imag=',aimag(testc)
100  FORMAT(a5,1x,f30.20)

     DO h=1,100
        theta = 0.01q0*h
        write(*,*)'trig formula:',qsqrt(qsin(theta)**2+qcos(theta)**2)-1.q0
     ENDDO

end program test_det
