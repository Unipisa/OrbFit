      double precision function dcabs1(z)
       complex*16 z
c      complex*16 z,zz
c      double precision t(2)
c      equivalence (zz,t(1))
c      zz = z
c      dcabs1 = dabs(t(1)) + dabs(t(2))
c       dcabs1=abs(real(z))+abs(real(z*(0,1)))
       dcabs1=abs(DBLE(z))+abs(DBLE(z*(0,1)))

       return
      end
