c =======================================
c  FMUOBS
c =======================================
c output of multiple observations 
c ===============INTERFACE========================
      SUBROUTINE fmuobs(tc,gmag,iobs,ids,t1,tut1,sigma,eqm,hmu,
     +     aobs,dobs,iff,imim,imip,imi0,titnam,filnam,iun20)
      IMPLICIT NONE
c reference orbit: epoch, G magnitude 
      DOUBLE PRECISION tc,gmag
c first last, and reference index of multiple orbits
      INTEGER imim,imip,imi0
c multiple solutions elements, magnitudes,sigma value
      DOUBLE PRECISION eqm(6,imip),hmu(imip),sigma
c actual observations
      DOUBLE PRECISION aobs,dobs
c strings with asteroid names, output unit
      CHARACTER*80 titnam
      CHARACTER*60 filnam
      INTEGER iun20
c station code, flag for use of act.obs., observation type
      INTEGER ids,iff,iobs
c target time (TDT, UTC)
      DOUBLE PRECISION t1,tut1
c =================END INTERFACE====================
      INCLUDE 'parmul.h'
      DOUBLE PRECISION alm(mulx),dem(mulx)
      INTEGER ng,i,npop
      DOUBLE PRECISION eqm1(6,mulx),amagn(mulx),eq1(6),alpha,delta
c trig constants
      INCLUDE 'trig.h'
c observation auxiliary data (e.g. distance)
      INCLUDE 'phase.h'
c multiple data for confidence boundary
      INCLUDE 'npoint.h'
c covariance (here it is dummy)
      integer icov
      double precision gamad(2,2),axes(2,2),sig(2)
c ====================================================
c this routine does not handle confidece elllipses 
      icov=1
c first compute nominal prediction
      CALL proele('EQU',tc,eqm(1,imi0),t1,eq1)
      CALL preobs('EQU',t1,ids,t1,eq1,iobs,alpha,delta,hmu(imi0),gmag
     + ,amagn(imi0))
      WRITE(*,*)' nominal solution '
      CALL outobc(iun20,iobs,ids,tut1,alpha,delta,amagn(imi0),adot,ddot,
     +     elo,dis,icov,gamad,sig,axes)
      alm(imi0)=0.d0
      dem(imi0)=0.d0
      adotv(imi0)=adot
      ddotv(imi0)=ddot
      disv(imi0)=dis
c loop on existing multiple solutions: first forward, then backward
      DO 147 i=imi0+1,imip
         CALL proele('EQU',tc,eqm(1,i),t1,eqm1(1,i))
         CALL preobs('EQU',t1,ids,t1,eqm1(1,i),iobs,alm(i),dem(i),
     +        hmu(i),gmag,amagn(i))
         IF(alm(i).gt.pig)alm(i)=alm(i)-dpig
         disv(i)=dis
         adotv(i)=adot
         ddotv(i)=ddot
         WRITE(*,*)' alternate obs.no. ',i
         CALL outobc(iun20,iobs,ids,tut1,alm(i),dem(i),amagn(i),
     +     adot,ddot,elo,dis,icov,gamad,sig,axes)
         alm(i)=alm(i)-alpha
         dem(i)=dem(i)-delta
c keep count of lost revolutions in alpha
         IF(i.eq.imi0)THEN
c initialize revolution counter
            ng=0
         ELSE
c update revolution counter
            CALL angupd(alm(i),alm(i-1),ng)
         ENDIF
 147  CONTINUE
      DO 148 i=imi0-1,imim,-1
         CALL proele('EQU',tc,eqm(1,i),t1,eqm1(1,i))
         CALL preobs('EQU',t1,ids,t1,eqm1(1,i),iobs,alm(i),dem(i),
     +        hmu(i),gmag,amagn(i))
         disv(i)=dis
         adotv(i)=adot
         ddotv(i)=ddot
         WRITE(*,*)' alternate obs.no. ',i
         CALL outobc(iun20,iobs,ids,tut1,alm(i),dem(i),amagn(i),
     +     adot,ddot,elo,dis,icov,gamad,sig,axes)
         alm(i)=alm(i)-alpha
         dem(i)=dem(i)-delta
c keep count of lost revolutions in alpha
         IF(i.eq.imi0)THEN
c initialize revolution counter
            ng=0
         ELSE
            CALL angupd(alm(i),alm(i+1),ng)
         ENDIF
 148  CONTINUE
c ===============================================
c output multiple prediction of observations
      CALL outmul(titnam,filnam,tut1,sigma,alpha,delta,
     +        alm,dem,hmu,imim,imip,imi0,iff,aobs,dobs,iobs)
      RETURN
      END
c
c ====================================================
c FMUPRO multiple state propagation for FITOBS
c ====================================================
      SUBROUTINE fmupro(iun20,imim,imip,t0,eqm,hmu,gm,cm,tr,eq1,g1,c1)
      IMPLICIT NONE
c maximum number of alternate solutions
      INCLUDE 'parmul.h'
c =================INPUT=========================================
c min and max index of alternate orbits
      INTEGER imim,imip
c orbits, covariance, normal matrices, abs.magnitude
      DOUBLE PRECISION eqm(6,imip-imim+1),gm(6,6,imip-imim+1),
     +     cm(6,6,imip-imim+1),hmu(imip-imim+1)
c requirements on covariance: assumed available in gm,cm
c output units
      INTEGER iun20
c epoch time, target time
      DOUBLE PRECISION t0,tr
c ================OUTPUT=================================
c elements, covariance and normal matrix at epoch tr
      DOUBLE PRECISION eq1(6,imip),g1(6,6,imip),
     +        c1(6,6,imip)
c ================END INTERFACE==========================
c keplerian elements, mean motion, opposition effect
c     DOUBLE PRECISION ekr(6),enne,gmag
c ======== output moid =====================
      DOUBLE PRECISION moid(mulx), dnp(mulx), dnm(mulx)
      INTEGER iconv(mulx)
c loop indexes
      INTEGER j,i
c =====================================================================
c main loop
c propagation to time tr
      DO j=imim,imip
        WRITE(*,*)' orbit ',j
c       CALL proelc('EQU',t0,eqm(1,j),gm(1,1,j),cm(1,1,j),
c    +          tr,eq1(1,j),g1(1,1,j),c1(1,1,j))
        CALL proele('EQU',t0,eqm(1,j),tr,eq1(1,j))
c orbital distance
        CALL nomoid(tr,eq1(1,j),moid(j),
     +              iconv(j),dnp(j),dnm(j))
      ENDDO
c =====================================================================
c summary table
c =====================================================================
      CALL tee(iun20,'SUMMARY OF MULTIPLE SOLUTIONS=')
      WRITE(iun20,223) tr
      WRITE(*,223) tr
 223  FORMAT(' elements at time ',f8.1,' (MJD):')
      CALL tee(iun20,
     +  'no.,     a      h      k      p      q      lambda=') 
      DO i=imim,imip
        WRITE(*,144)i,(eqm(j,i),j=1,6)
 144    FORMAT(i3,6f12.8)
        WRITE(iun20,144)i,(eqm(j,i),j=1,6)
      ENDDO
      CALL tee(iun20,'no.,  magn,  MOID ,  nod+  ,  nod-=')
      DO i=imim,imip
        WRITE(*,145)i,hmu(i),moid(i),dnp(i),dnm(i),iconv(i)
        WRITE(iun20,145)i,hmu(i),moid(i),dnp(i),dnm(i),iconv(i)
 145    FORMAT(i3,2x,f5.2,1x,f8.5,1x,f8.5,1x,f8.5,1x,i2)
      ENDDO
      RETURN
      END
c
c =======================================
c  FMUPLO
c =======================================
c graphic output of multiple orbits
c ===============INTERFACE========================
      SUBROUTINE fmuplo(eqm,tc,numb,eqc,titnam,sigma)
      IMPLICIT NONE
      INTEGER numb
      CHARACTER*80 titnam
      DOUBLE PRECISION tc,sigma,eqc(6),eqm(6,numb)
c ============END INTERFACE==================
c number of alternate solutions, maximum
      INCLUDE 'parmul.h'
      DOUBLE PRECISION a(mulx),e(mulx),aa,ee
      INTEGER i
c ===========================================
c a-e plot of multiple solutions
      DO i=1,numb
        a(i)=eqm(1,i)
        e(i)=sqrt(eqm(2,i)**2+eqm(3,i)**2)
      ENDDO
      aa=eqc(1)
      ee=sqrt(eqc(2)**2+eqc(3)**2)
      CALL ploae(tc,a,e,aa,ee,sigma,numb,titnam)
      RETURN
      END
c
c =======================================
c mult_input initializes storage of multiple solution data
c =======================================
      SUBROUTINE mult_input(catname,eqa,cc,gg,hmu,tcat,
     +  m1,m2,m0,vel_inf,ok)
      IMPLICIT NONE
c ------------INPUT------------------
c file with catalog of multiple solutions
      CHARACTER*160 catname
c ------------OUTPUT-----------------
c multiple output arrays
      INCLUDE 'parmul.h'
      DOUBLE PRECISION eqa(6,mulx)
      DOUBLE PRECISION gg(6,6,mulx),cc(6,6,mulx),tcat,hmu(mulx)
c index range of multiple solution 
      INTEGER m1,m2,m0
c velocity at infinite with respect to Earth (circular approx)
      DOUBLE PRECISION vel_inf
c succesful input
      LOGICAL ok
c ------------END INTERFACE---------
c asteroid name
c     CHARACTER*9 astname
c -----for call to rdorb----
c names
      CHARACTER*19 name0
      CHARACTER*9 name1
      INTEGER le
c magnitude slope parameter, mass, epoch
      double precision sl,mass,t
c type of orbital elements (KEP/EQU/CAR),
      character*3 eltype 
c reference system type (EQUM/EQUT/ECLM),
      character*4 rsys
c epoch specification (J2000/OFDATE)
      character*6 epoch
c record number
      integer no
c avalaibility of covariance and normal matrices
      logical defcov,defnor
c end of file
      LOGICAL eof
c indexes of multiple solutions
      INTEGER imul(mulx),norb,j
c velocity w.r.to Earth for each orbit
      DOUBLE PRECISION v_inf(mulx),v_infty,v_max,v_min
c temporary input arrays
      DOUBLE PRECISION eq(6)
      DOUBLE PRECISION g(6,6),c(6,6),h
c loop indexes
      INTEGER i
c system dependencies
      INCLUDE 'sysdep.h'
c -------------------------------------------------
      v_min=100.d0
      v_max=0.d0
c opening and reading multiple solution catalog
      CALL rmsp(catname,le)
      INQUIRE(file=catname(1:le),exist=ok)
      IF(.not.ok)THEN
          WRITE(*,*)' file ',catname(1:le),' not found'
          RETURN
      ENDIF
      CALL oporbf(catname(1:le))
      DO i=1,mulx
         CALL rdorb(name0,eq,eltype,t,g,defcov,
     +        c,defnor,h,sl,mass,rsys,epoch,no,eof) 
         IF(eof)THEN
            norb=i-1
            GOTO 2
         ENDIF
c control on time
         IF(i.eq.1)THEN
            tcat=t
         ELSE
            IF(t.ne.tcat)THEN
               WRITE(*,*)'mult_input: time discrepancy from tcat=',tcat
               WRITE(*,*)'mult_input: at record ',i,' t=',t
               ok=.false.
               RETURN
            ENDIF
         ENDIF
c control on elements type
         IF(eltype.ne.'EQU')THEN
            WRITE(*,*)'mult_input: non equinoctal, but of type ',eltype
            ok=.false.
            RETURN
         ENDIF
c availability of covariance
         IF (.NOT.(defnor.AND.defcov)) THEN
            WRITE(*,*)'mult_input',  name0, ' matrices not avalaible'
            ok=.false.
            RETURN
         ENDIF
c handling of name; but note that the multiple solutions are assumed to be orderd
c and with consecutive indexes!!! that is imul(i)=imul(1)+i-1
         CALL splinam(name0,name1,imul(i))
         IF(imul(i).ne.imul(1)+i-1)THEN
            WRITE(*,*)'mult_input: indexes not in order: ',
     +           imul(i),' at record ',i,' should be ',imul(1)+i-1
c            RETURN
         ENDIF
         j=imul(i)
c copy in output arrays
         CALL vcopy(6,eq,eqa(1,j))
         hmu(j)=h
         CALL mcopy(6,6,c,cc(1,1,J))
         CALL mcopy(6,6,g,gg(1,1,J))
c v_infinity computation
         v_inf(i)=v_infty(eqa(1,i))
         v_min=min(v_min,v_inf(i))
         v_max=max(v_max,v_inf(i))
      ENDDO
c increase nmax
      WRITE(*,*)'mult_input: increase mulx, file is longer than',mulx
 2    CONTINUE
      norb=i-1
      ok=norb.gt.0
      m1=imul(1)
      m2=imul(norb)
      WRITE(*,*)' mult_input: input of ',norb,' multiple solutions'
      WRITE(*,*)' indexes between ', m1, ' and ', m2
      WRITE(*,*)' max of v_infty ',v_max,' min ',v_min
      vel_inf=v_min
 3    WRITE(*,*)' which one is the nominal?'
      READ(*,*)m0
      IF(m0.lt.m1.or.m0.gt.m2)THEN
          WRITE(*,*)' OUT OF RANGE ',m1,m2
          GOTO 3
      ELSE
         WRITE(*,*)' nominal is ',m0
      ENDIF
      RETURN
      END
c
c =====================================================================
c  OUTMUL
c =====================================================================
c output multiple observations
c =============INTERFACE===============================================
      SUBROUTINE outmul(titnam,filnam,t1,sigma,alpha,delta,
     +     alm,dem,hmagn,imim,imip,imi0,iff,aobs,dobs,iobs)
      IMPLICIT NONE
c =============INPUT===================================================
c file name
      CHARACTER*80 titnam
      CHARACTER*60 filnam
c first and last index of existing multiple solutions, 
c index of reference one, control for closed curve,obs type
      INTEGER imim,imip,imi0,iff,iobs
c observation time MJD, sigma value, nominal prediction
      DOUBLE PRECISION t1,sigma,alpha,delta
c observation: predicted  value alpha, delta, magnitude, actual
      DOUBLE PRECISION alm(imip),dem(imip),hmagn(imip),aobs,dobs
c =============END INTERFACE===========================================
c multiple data for confidence boundary
      INCLUDE 'npoint.h'
c trig constants
      INCLUDE 'trig.h'
c max no. orbits
      INCLUDE 'parmul.h'
c conversion to sessagesimal
      DOUBLE PRECISION seca,secd
      INTEGER inta,mina,intd,mind
      CHARACTER*1 siga,sigd
c conversion of time
      INTEGER iy,imo,iday
      DOUBLE PRECISION hour
c scalar temporaries, differences
      DOUBLE PRECISION dee,daa,ado,ddo
c file name
      INTEGER le
      CHARACTER*80 file
c loop indexes, units
      INTEGER n,iun7
c ======================================================================
c open output file
      CALL rmsp(filnam,le)
      file=filnam(1:le)//'.cbd'
      CALL filopn(iun7,file,'unknown')
c date and sigma value
      CALL mjddat(t1,iday,imo,iy,hour)
      WRITE(iun7,297)iday,imo,iy,hour,sigma
 297  FORMAT(i3,i3,i5,f8.4,f5.2)
c line of variations
      DO n=imim,imip
        daa=alpha+alm(n)
        daa=mod(daa,dpig)
        IF(daa.lt.0.d0)daa=daa+dpig
        IF(daa.gt.dpig)daa=daa-dpig
        daa=daa*degrad/15
        IF(daa.lt.0.d0.or.daa.gt.24.d0)THEN
           WRITE(*,*)' outmul: daa out of range ', daa
        ENDIF
        CALL sessag(daa,siga,inta,mina,seca)
        dee=(delta+dem(n))*degrad
        CALL sessag(dee,sigd,intd,mind,secd)
c proper motion in arcsec/hour
        ado=adotv(n)*secrad/24.d0
        ddo=ddotv(n)*secrad/24.d0
c output
        IF(siga.eq.'+')siga=' '
        WRITE(iun7,396)n,siga,inta,mina,seca,sigd,intd,mind,secd,
     +       disv(n),ado,ddo,hmagn(n)
 396    FORMAT(i3,1x,a1,i2,1x,i2,1x,f4.1,2x,a1,i2,1x,i2,1x,f4.1,
     +       1x,f8.5,1x,f8.2,1x,f8.2,1x,f5.2)
      ENDDO
      CALL filclo(iun7,' ')
c graphics output
      IF(iff.eq.1)THEN
         CALL plocbd(titnam,alpha,delta,sigma,t1,
     +         alm(imim),dem(imim),imip-imim+1,iobs)
      ELSEIF(iff.eq.2)THEN
         CALL ploobs(titnam,alpha,delta,sigma,t1,
     +         alm(imim),dem(imim),imip-imim+1,
     +                 aobs,dobs)
      ENDIF
      RETURN
      END










