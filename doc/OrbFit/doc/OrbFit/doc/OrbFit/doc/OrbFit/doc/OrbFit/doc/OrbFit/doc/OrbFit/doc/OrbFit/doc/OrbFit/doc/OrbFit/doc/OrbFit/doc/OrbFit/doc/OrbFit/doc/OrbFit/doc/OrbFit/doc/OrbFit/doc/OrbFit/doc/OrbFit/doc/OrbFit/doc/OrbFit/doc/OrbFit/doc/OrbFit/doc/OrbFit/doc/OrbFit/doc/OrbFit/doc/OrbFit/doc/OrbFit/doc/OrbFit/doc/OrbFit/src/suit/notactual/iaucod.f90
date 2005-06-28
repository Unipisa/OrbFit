      SUBROUTINE iaucod(mpccod,iaudes,error) 
      IMPLICIT NONE 
                                                                        
      CHARACTER*(*) iaudes,mpccod 
      LOGICAL error 
                                                                        
      INTEGER ln,i,temp 
      CHARACTER*2 head 
      CHARACTER*3 tail 
                                                                        
      INTEGER lench 
      LOGICAL isnum 
      EXTERNAL lench,isnum 
                                                                        
      CHARACTER*12 numfield, desfield 
                                                                        
      error=.false. 
      iaudes=' ' 
                                                                        
      ln=lench(mpccod) 
      IF(ln.LE.0) GOTO 10 
                                                                        
      numfield=mpccod(1:5) 
      call rmsp(numfield,ln) 
! Numbered asteroids                                                    
      if(ln.ne.0) then 
         iaudes=numfield 
         do i=1,ln-1 
            if (iaudes(i:i).eq.'0') then 
               iaudes(i:i)=' ' 
            else 
               goto 123 
            endif 
         enddo 
  123    call rmsp(iaudes,ln) 
         return 
      endif 
                                                                        
! Unnumbered asteroids                                                  
      desfield=mpccod(6:12) 
                                                                        
      if(desfield(3:3).eq.'S')then 
! Survey Asteroid                                                       
         iaudes=desfield(4:7)//desfield(1:1)//'-'//desfield(2:2) 
                                                                        
      elseif (desfield(1:1).eq.'I' .or.                                 &
     &        desfield(1:1).eq.'J' .or.                                 &
     &        desfield(1:1).eq.'K')then                                 
! Temporary designation                                                 
         if(desfield(5:6).eq.'00')then 
!           1999AA = J99A00A                                            
            tail='' 
         elseif(desfield(5:5).eq.'0')then 
!           1999AA1 = J99A01A                                           
            tail=desfield(6:6) 
         elseif(isnum(desfield(5:5)))then 
!           1999AA12 = J99A12A                                          
            tail=desfield(5:6) 
         else 
!           1999AA103 = J99AA3A                                         
            temp=ichar(desfield(5:5))-55 
!           1999AA363 = J99Aa3A                                         
            if (temp.gt.35) temp = temp - 6  
            write(head,103) temp 
            tail=head//desfield(6:6) 
         endif 
         temp=ichar(desfield(1:1))-55 
         write(head,103) temp 
  103    format(I2) 
         iaudes=head//desfield(2:3)//desfield(4:4)//desfield(7:7)//tail 
      else 
! Unknown type                                                          
         write(*,*)'cannot understand MPC designation: ',mpccod 
         goto 10 
      endif 
      return 
                                                                        
   10 CONTINUE 
      iaudes=mpccod 
      error=.true. 
                                                                        
      END SUBROUTINE iaucod 
