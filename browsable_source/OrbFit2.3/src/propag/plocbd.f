c =====================================================
c PLOAE graphic (GNUPLOT) output of multiple solutions
c =====================================================
      SUBROUTINE ploae(t0,a,e,aa,ee,sigma,npo,titnam)
      implicit none
c =========INPUT=================
      double precision aa,ee,t0,sigma
      integer npo
      double precision a(npo),e(npo)
      character*80 titnam
c =====================================================
c captions
      character*60 ylab
      character*60 xlab,title
      integer idev,istyle
      integer le
c =====================================================
      call rmsp(titnam,le)
      write(title,100)titnam(1:le),t0,sigma
 100  format('Asteroid ',a13,', epoch ',f13.6,
     + ' MJD, level',f6.2)
      xlab='Semimajor axis AU'
      ylab='Eccentricity'
 2    call getdev(idev)
      if(idev.eq.0)return
      istyle=1
      call plotob(a,e,aa,ee,npo,xlab,ylab,title,idev,istyle)
      goto 2
      end
c =====================================================
c GETDEV get plotting device code for passing to plotob  
c =====================================================
      SUBROUTINE getdev(idev)
      implicit none
      integer idev
      CHARACTER*20 menunam
      CHARACTER*70 s8,s9,s10
      menunam='plotdevice'
      CALL menu(idev,menunam,7,
     +     'Select gnuplot output device:=',
     +     'X Windows (X11)=',
     +     'MS Windows=',
     +     'VGA=',
     +     'Tektronix=',
     +     'PostScript file (black & white)=',
     +     'PostScript file (color)=',
     +     'Pipe to PostScript printer (''|lpr'')=',
     +     s8,s9,s10)
      return
      end
c =====================================================
c PLOCBD graphic (GNUPLOT) output of confidence boundary
c =====================================================
      SUBROUTINE plocbd(titnam,alpha,delta,sigma,tobs,da,dd,npo,iob1)
      implicit none
c =========INPUT=================
      double precision alpha,delta,tobs,sigma
      integer npo,iob1
      double precision da(npo),dd(npo)
      character*80 titnam
c =====================================================
c captions
      character*60 ylab,xlab,ylabt,xlabt,title
      integer idev,istyle
      integer le
c trigonometric constants
      include 'trig.h'
      include 'jplhdr.h'
c rescaled plotting points
      integer npox,i,xlen,ylen,lench
      parameter (npox=4000)
      double precision xd(npox),yd(npox),scale
c =====================================================
      call rmsp(titnam,le)
      write(title,100)titnam(1:le),tobs,sigma
 100  format('Asteroid ',a13,', time',f13.6,
     + ' MJD, level',f6.2)
      if(iob1/1000.eq.1)then
         xlabt=' Differences in right ascension (DEG)  from '
         ylabt=' Differences in declination (DEG)  from '
         scale=degrad
         xlen=lench(xlabt)
         ylen=lench(ylabt)
         write(xlab,101)xlabt(1:xlen),alpha*scale
         write(ylab,101)ylabt(1:ylen),delta*scale
 101     format(a,' ',f8.4)
      elseif(iob1/1000.eq.2)then
         xlabt=' Differences in range (km)  from '
         ylabt=' Differences in range rate (km/day)  from '
         scale=au
         xlen=lench(xlabt)
         ylen=lench(ylabt)
         write(xlab,102)xlabt(1:xlen),alpha*scale
         write(ylab,102)ylabt(1:ylen),delta*scale
 102     format(a,' ',f15.4)
      else
         stop 'plocbd: internal error'
      endif
      do  i=1,npo
         xd(i)=da(i)*scale
         yd(i)=dd(i)*scale
      enddo
 2    call getdev(idev)
      if(idev.eq.0)return
      istyle=1
      call plotob(xd,yd,0.d0,0.d0,npo,xlab,ylab,title,idev,istyle)
      goto 2
      end
c =====================================================
c PLOOBS graphic (GNUPLOT) output of confidence boundary
c compared with the real observation at the same time
c =====================================================
      SUBROUTINE ploobs(titnam,alpha,delta,sigma,tobs,da,dd,npo,
     +     aobs,dobs)
      implicit none
c =========INPUT=================
      double precision alpha,delta,sigma,tobs,aobs,dobs
      integer npo
      double precision da(npo),dd(npo)
      character*80 titnam
c =====================================================
c captions
      character*60 ylab
      character*60 xlab,title
      integer idev,istyle
      integer le
c trigonometric constants
      include 'trig.h'
c rescaled plotting points
      integer npox,i
      parameter (npox=4000)
      double precision xd(npox),yd(npox),xx,yy
c =====================================================
      call rmsp(titnam,le)
      write(title,100)titnam(1:le),tobs,sigma
 100  format('Asteroid ',a13,', time',f13.6,
     + ' MJD, level',f6.2)
      write(xlab,101)degrad*alpha
 101  format(' Differences in right ascension (DEG)  from ',f8.4)
      write(ylab,102)degrad*delta
 102  format(' Differences in declination (DEG)  from ',f8.4)
      do  i=1,npo
        xd(i)=da(i)*degrad
        yd(i)=dd(i)*degrad
      enddo
c residual observed-computed
      xx=(aobs-alpha)*degrad
      yy=(dobs-delta)*degrad
c residual in alpha reduced to interval -180 180
      if(xx.gt.180)xx=xx-360
      if(xx.lt.-180)xx=xx+360
 2    call getdev(idev)
      if(idev.eq.0)return
      istyle=1
      call plotob(xd,yd,xx,yy,npo,xlab,ylab,title,idev,istyle)
      goto 2
      end
c **************************************************************
c PLOTOB
c plotting routine using GNUPLOT; line plus point
c
      SUBROUTINE plotob(xd,yd,xx,yy,n,xlab,ylab,title,idev,istyle)
c maximum number of data points
      parameter (nx=36000)
      real x(nx),y(nx),x1,y1
      double precision xd(n),yd(n),xx,yy
      character*60 ylab
      character*60 xlab,title
      character*62 label
      character*10 tempfi
      character*6 style
c write temporary file with data: boundary
      tempfi='confi.tmp'
      open(22,file=tempfi,status='unknown')
      do 1 i=1,n
        x(i)=xd(i)
        y(i)=yd(i)
 1      write(22,*)x(i),y(i)
      close(22)
c write temporary file with data: real observation
      tempfi='obs.tmp'
      open(22,file=tempfi,status='unknown')
      x1=xx
      y1=yy
      write(22,*)x1,y1
      close(22)
c plot style
      if(istyle.eq.1)then
        style='lines'
      elseif(istyle.eq.0)then
        style='dots'
      elseif(istyle.lt.0)then
        style='points'
      else
        write(*,*)' style code ',istyle,' not known'
        return
      endif
c  preparation of command file
      open(22,file='giffv.gnu')
c  xlable,ylabel, title
      write(22,*)'set nokey'
      nc=lench(xlab)
      label='"'//xlab(1:nc)//'"'
      write(22,*)'set xlabel ',label(1:nc+2)
      nc=lench(ylab)
      label='"'//ylab(1:nc)//'"'
      write(22,*)'set ylabel ',label(1:nc+2)
      nc=lench(title)
      label='"'//title(1:nc)//'"'
      write(22,*)'set title ',label(1:nc+2)
c  set graphics terminal type and output device
      if(idev.eq.1)then
         write(22,*)'set terminal X11'
         write(22,*)'set output'
      elseif(idev.eq.2)then
         write(22,*)'set terminal windows'
         write(22,*)'set output'
      elseif(idev.eq.3)then
         write(22,*)'set terminal vgalib'
         write(22,*)'set output'
      elseif(idev.eq.4)then
         write(22,*)'set terminal tek40xx'
         write(22,*)'set output'
      elseif(idev.eq.5)then
         write(22,*)'set terminal postscript monochrome'
         write(22,*)'set output "giffv.ps"'
         write(*,*)'Generating Postscript file ''giffv.ps''.'
      elseif(idev.eq.6)then
         write(22,*)'set terminal postscript color'
         write(22,*)'set output "giffv.ps"'
         write(*,*)'Generating Postscript file ''giffv.ps''.'
      elseif(idev.eq.7)then
         write(22,*)'set terminal postscript monochrome'
         write(22,*)'set output "|lpr"'
         write(*,*)'Sending plot to lpr.'
      else
         write(*,*)' this device flag ',idev,' not known'
         return
      endif
      write(22,123)style
 123  format('plot ''confi.tmp''  with ',a6,',''obs.tmp'' with points') 
c  pause for screen images
      if(idev.le.4)then
         write(22,*)'pause -1'
      endif
      close(22)
c  on SONY NEWS system is a function
c     ii=system('gnuplot giffv.gnu')
c     write(*,*)ii
c  on IBM RISC system is a subroutine
      if(idev.eq.4)then
         ii=system('xterm -t -e gnuplot giffv.gnu')
      else
         ii=system('gnuplot giffv.gnu')
      endif
c  hard copy (if not done already)
c     if(idev.eq.-7)then
c        ii=system('lpr giffv.ps')
c        write(*,*)ii
c        call system('lpr giffv.ps')
c     elseif(idev.eq.-5)then
c        call system('lpr giffv.ps -h')
c     endif
      return
      end







