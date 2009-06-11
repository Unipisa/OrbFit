program conta
  implicit none
  INTEGER:: m,n,count,nmax
  count=0
  DO m=0,20
     IF (m.le.4) THEN
        nmax=20
     ELSEIF(m.gt.4) THEN
        nmax = 2*((24-m)/2)
     ENDIF
     DO n = 0,nmax
        count=count+1     
     ENDDO
  ENDDO
  write(*,*)'count=',count  

end program conta
