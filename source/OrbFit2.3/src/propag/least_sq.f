c =============MODULE least_sq==================
c PUBLIC ROUTINES 
c                fdifco
c INTERNAL ROUTINES
c                outcov
c                radres
c MODULE CONTAINS
c ROUTINES
c                difcor 
c                sincor
c                minsol
c                  invmat
c                srtoss
c                  unsort
c                rescov
c                difini
c
c                reject
c                  entopp
c                rejini
c
c                magest
c
c                meanti
c                bizarre
c
c  HEADERS
c least_sq.o: comdif.h jplhdr.h mag.h phase.h proout.h
c               
c           codode.h private
c           comrej.h private, used by orbfit.x 
c                       and others in idconf2 (only for setting autrej)
c                               
c           iunrad.h  unit for radar residuals
c            parobx.h only to dimension weights and residuals
c           trig.h public
c
* Copyright (C) 1998 by OrbFit Consortium
* Version: December 15, 1997 Steven Chesley
c =========================================
c    FDIFCO
c =========================================
c differential corrections interface routine for FITOBS
c =========INTERFACE=============================
      SUBROUTINE fdifco(iarc,obs0,ini0,ok,cov0,
     +     t0,eq0,m,objid,iobs,tau,tutm,idsta,aln,den,rmsa,rmsd,
     +     rmsmag,smag,h0,sel,nobs,
     +     rwofi0,iun20,iun8,eqc,g0,c0,csino0,delno0,succ)
      IMPLICIT NONE
c =========INPUT==============================
c arc: 1,2 arcs, 3 together
      INTEGER iarc
c data availability
      LOGICAL obs0,ini0
c initial conditions (first guess)
      DOUBLE PRECISION t0,eq0(6)
c ======observations====
c number of observations, station codes, obs. type, object id
      INTEGER m,idsta(m),iobs(m)
      CHARACTER*(*) objid(m)
c observation times (ET, UT), alpha, delta, a priori rms
      DOUBLE PRECISION tau(m),tutm(m),aln(m),den(m),rmsa(m),rmsd(m)
c magnitudes (string), a priori rms'
      CHARACTER*6 smag(m)
      DOUBLE PRECISION rmsmag(m)
c output units 
      CHARACTER*60 rwofi0
      INTEGER iun20,iun8
c ==========INPUT/OUTPUT=============================
c selection flags, no. observations used and weights
      INTEGER sel(m),nobs
      INCLUDE 'parobx.h'
c magnitude (a priori, estimate)
      DOUBLE PRECISION h0
c ========= OUTLIER REJECTION ===========
      INCLUDE 'comrej.h'
      INCLUDE 'comdif.h'
c =========OUTPUT=============================
      INCLUDE 'verbosity.h'
c data available, differential correction convergent,matrices available
      LOGICAL ok,succ,cov0
c corrected orbital elements
c WARNING: if accepted, they are also copied in eq0
      DOUBLE PRECISION eqc(6)
c normal and covariance matrix (warning: not full if ncor.lt.6)
      DOUBLE PRECISION c0(6,6),g0(6,6)
c norm of residuals, of last correction
      DOUBLE PRECISION csino0,delno0
c =========END INTERFACE======================
c no. solve for variables, list, pseudoNewton/Newton
      INTEGER ncor,icor(6),inew,inter
c rms of residuals
      DOUBLE PRECISION rms,sigma
c controls for iterations to be passed to difcor
      DOUBLE PRECISION delcr
c residuals
      DOUBLE PRECISION csir(nob2x),resa(nobx),resd(nobx),x2(nobx)
      DOUBLE PRECISION w(nob2x),ratio
      PARAMETER (ratio=0.6666666d0)
      INTEGER itsav
c fit residuals, rms
      DOUBLE PRECISION resmag(nobx),rmsh
      LOGICAL radar
c data to be used to output radar reports
      INCLUDE 'iunrad.h'
c discarded observations?
      LOGICAL disc
c loop indexes, counters, modes, units
      INTEGER i,nobnew,idif,iunf
c scalar temporaries
      DOUBLE PRECISION ra,rd,recov,disca
c characters for menu 
      CHARACTER*20 menunam
      CHARACTER*70 s5,s6,s7,s8,s9,s10
c trigonometric constatnts
      INCLUDE 'trig.h'
c =====================================================================
c check availability of observations and initial condition
      CALL cheobs(obs0,ini0,ok)
      IF(.not.ok) RETURN
c =====================================================================
c check availability of JPL ephemerides and ET-UT table
      CALL chetim(tau(1),tau(m),ok)
      IF(.not.ok) RETURN
c ok, go on with arc
      IF(.not.batch.and.verb_dif.gt.9)THEN
         IF(iarc.eq.1)THEN
            CALL tee(iun20,' FIRST ARC=')
         ELSEIF(iarc.eq.2)THEN
            CALL tee(iun20,' SECOND ARC=')
         ELSEIF(iarc.eq.4)THEN
            CALL tee(iun20,' BOTH ARCS TOGETHER=')
         ELSEIF(iarc.eq.0)THEN
            CONTINUE
         ELSE
            WRITE(*,*)'FDIFCO: this should not happen, iarc=',iarc
         ENDIF
      ENDIF
c =====================================================================
c choose method, how many elements to correct (interactively???)
      IF(batch)THEN
c leave autoreject control autrej as it is set 
c in options file or difini.def
c
c solve for all elements anyway
         inter=0
         itsav=itmax
         CALL whicor(inter,icor,ncor,inew)
      ELSE
         menunam='difcomod'
         CALL menu(idif,menunam,4,' select correction and reject mode=',
     +      'correct all, autoreject=',
     +      'correct all, manual reject=',
     +      'correct only some elements (reject is manual)=',
     +      'compute residuals and covariance w/o correcting=',
     +      s5,s6,s7,s8,s9,s10)
         itsav=itmax
         if(idif.eq.0)then
            ok=.false.
            return
         endif
         if(idif.eq.4)itmax=0
c auto/manual reject
         IF(idif.eq.1)THEN
            autrej=.true.
         ELSE
            autrej=.false.
         ENDIF
c which elements to correct
         IF(idif.ne.3)THEN
            inter=0
         ELSE
            inter=1
         ENDIF 
         CALL whicor(inter,icor,ncor,inew)
      ENDIF
c default iteration control parameters
      delcr=1.d-3
c weights (non-zero!)
      call fitwgt(rmsa,rmsd,den,sel,iobs,w,m,.false.)
c differential corrections
      IF(batch)THEN
         iunf=-iun20
      ELSE
         iunf=iun20
      ENDIF
      CALL difcor(m,w,sel,t0,iobs,tau,idsta,eq0,aln,den,icor,inew,
     +         iunf,delcr,eqc,g0,c0,csino0,delno0,csir,x2,succ)
      itmax=itsav
c estimate magnitude here
      IF(succ)THEN
         CALL magest(smag,rmsmag,sel,m,h0,resmag,rmsh)
      ELSE
         IF(verb_dif.gt.9)
     +        WRITE(*,*)' magnitude cannot be estimated w/o new orbit'
         rmsh=-1.d0
         DO i=1,m
            resmag(i)=1.d9
         ENDDO
      ENDIF
c residuals and weights written on .rwo file
      DO i=1,m
         resa(i)=csir(2*i-1)
         resd(i)=csir(2*i)
      ENDDO
c output of residuals, weights and observatiosn file
c warning: when the orbit is hyperbolic, these are the residuals of
c iteration one.
c when the divergence is mild (e.g. target function continues to increase,
c too many iterations) should be (TO BE FIXED)
      CALL wrirwo(rwofi0,objid,iobs,tutm,idsta,
     +     aln,rmsa,resa,den,rmsd,resd,smag,rmsmag,resmag,rmsh,
     +     sel,x2,m,csino0)
c Output radar residuals
      IF(.not.batch)THEN
         iunrad=iun20
         iunrare=iun20
c ELSE iunrad,iunrare is set in the main
      ENDIF
      CALL radres(iobs,tutm,resa,rmsa,resd,rmsd,sel,m,radar)
c is covariance available?
      IF(.not.succ)THEN
         cov0=.false.
      ELSE
c covariance matrix
         IF(ncor.eq.6)THEN
            cov0=.true.
         ELSE
            cov0=.false.
         ENDIF
c output covariance in .fga file, not in batch
         IF(.not.batch.and.verb_dif.gt.9)THEN 
            IF(iarc.eq.1)then
               WRITE(iun8,*) 'COVARIANCE MATRIX FOR FIRST ARC'
            ELSEIF(iarc.eq.2)then
               WRITE(iun8,*) 'COVARIANCE MATRIX FOR SECOND ARC'
            ELSEIF(iarc.eq.4)THEN
               WRITE(iun8,*)'COVARIANCE MATRIX FOR BOTH ARCS'
            ENDIF
            CALL outcov(iun8,icor,g0,c0)
         ENDIF
c improved orbital elements are accepted
         CALL vcopy(6,eqc,eq0)
c observations used
         nobs=0
         DO i=1,m
           IF(sel(i).gt.0)nobs=nobs+1
         ENDDO
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Manually discard big residuals
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         IF(.not.autrej .and. .not.batch)THEN
         rms=csino0
         WRITE(*,*)' RMS of weighed residuals is ',rms
         WRITE(*,*)' Discard residuals bigger than sigma*rms'
         WRITE(*,*)' give sigma; sigma.le.0 not to discard'
         WRITE(*,*)
         READ(*,*)sigma
         IF(sigma.gt.0.d0)THEN
            recov=sigma*rms*ratio
            disca=sigma*rms
            disc=.false.
            DO i=1,m
               ra=abs(resa(i))*sqrt(w(2*i-1))
               rd=abs(resd(i))*sqrt(w(2*i))
               IF(ra.gt.disca.or.rd.gt.disca)THEN
c discard observation
                  IF(sel(i).ne.0)THEN
                     WRITE(*,330)tau(i),idsta(i),resa(i)*secrad,
     +                    resd(i)*secrad,sel(i)
                     disc=.true.
                  ENDIF
                  WRITE(iun20,330)tau(i),idsta(i),resa(i)*secrad,
     +                 resd(i)*secrad,sel(i)
 330              FORMAT('OBS. at',f11.4,' from ',i3,
     +                 ' residuals ',2f10.3,' TO BE DISCARDED',i2) 
                  sel(i)=0
               ELSEIF(ra.lt.recov.and.rd.lt.recov)THEN
c recover previously discarded observation
                  IF(sel(i).eq.0)THEN
                     WRITE(*,340)tau(i),idsta(i),resa(i)*secrad,
     +                    resd(i)*secrad,sel(i)
                     WRITE(iun20,340)tau(i),idsta(i),resa(i)*secrad,
     +                    resd(i)*secrad,sel(i)
 340                 FORMAT('OBSERVATION at',f11.4,' from ',i3,
     +                    ' residuals ',2f13.4,' RECOVERED',i2)
                     sel(i)=1
                     disc=.true.
                  ENDIF
               ENDIF
            ENDDO
            IF(.not.disc)THEN
               WRITE(*,*)'NOTHING TO BE DISCARDED/RECOVERED'
               WRITE(*,*)' at sigma level=',sigma
            ENDIF
            nobnew=0
            DO i=1,m
               IF(sel(i).gt.0)nobnew=nobnew+1
            ENDDO
            WRITE(*,*)'OBSERVATIONS USED=',nobs,' TO BE USED=',nobnew
         ENDIF
         ENDIF
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c End Manual discard
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ENDIF
      RETURN
      END
c ==========================================================
c Output radar residuals
      SUBROUTINE radres(iobs,tutm,resa,rmsa,resd,rmsd,
     +     sel,m,radar)
      IMPLICIT NONE
c number of observations, obs. type,
      INTEGER m,iobs(m),sel(m)
c observation times, residuals, a priori rms
      DOUBLE PRECISION tutm(m),resa(m),resd(m),rmsa(m),rmsd(m)
c radar data flag
      LOGICAL radar
      INCLUDE 'verbosity.h'
c scalar temporaries, indexes
c radar residuals
      DOUBLE PRECISION rsum,rrsum,rwsum,rrwsum,hour
      INTEGER ii,rcnt,rrcnt,iday,imonth,iyear
c unit for output; if positive, include standard output
      INTEGER iunp
      LOGICAL radrange,radrate
c needed to use au
      INCLUDE 'jplhdr.h'
c data to be used to output radar reports
      INCLUDE 'iunrad.h'
c ===========================================================
      rsum=0d0
      rrsum=0d0
      rwsum=0d0
      rrwsum=0d0
      rcnt=0
      rrcnt=0
      radar=.false.
      radrange=.false.
      radrate=.false.
      do ii=1,m
c     if radar...
         if(iobs(ii)/1000.eq.2)then
            radar=.true.
c convert time
               CALL mjddat(tutm(ii),iday,imonth,iyear,hour)
c     range:               
            if(mod(iobs(ii),100).ne.3)then
c write residual in range
               radrange=.true.
               IF(verb_dif.gt.5)WRITE(iunrare,100)namast,iyear,
     +              imonth,iday,hour,resa(ii)*au,rmsa(ii)*au,sel(ii)
 100  FORMAT(a9,1x,i4,1x,i2,1x,i2,1x,'1',1x,f5.2,1x,f11.5,1x,f9.5,1x,i1)
               rsum=rsum+(resa(ii)/rmsa(ii))**2
               rwsum=rwsum+(1d0/rmsa(ii))**2
               rcnt=rcnt+1
            endif
c     range rate:              
            if(mod(iobs(ii),100).ne.2)then
               radrate=.true.
c write residual in range rate
               IF(verb_dif.gt.5)WRITE(iunrare,101)namast,iyear,
     +              imonth,iday,hour,resd(ii)*au,rmsd(ii)*au,sel(ii)
 101  FORMAT(a9,1x,i4,1x,i2,1x,i2,1x,'2',1x,f5.2,1x,f11.5,1x,f9.5,1x,i1)
               rrsum=rsum+(resd(ii)/rmsd(ii))**2
               rrwsum=rrwsum+(1d0/rmsd(ii))**2
               rrcnt=rrcnt+1
            endif
         endif
      enddo
      iunp=abs(iunrad)
      IF(radrange)THEN 
         IF(iunrad.gt.0)THEN
            WRITE (*,456)
     +           'range  RMS      =',sqrt(rsum/rwsum)*au,' km',
     +           'range weight =',sqrt(1d0/rwsum*rcnt)*au,' km'
         ENDIF
         WRITE (iunp,456)
     +        'range        =',sqrt(rsum/rwsum)*au,' km',
     +        'range weight =',sqrt(1d0/rwsum*rcnt)*au,' km'
 456     FORMAT(a,f10.5,a/a,f10.5,a)
      ENDIF
      IF(radrate)THEN
         IF(iunrad.gt.0)THEN
         WRITE(*,457)
     +           'range rate RMS   =',sqrt(rrsum/rrwsum)*au,' km/day',
     +           'range rate wt=',sqrt(1d0/rrwsum*rrcnt)*au,' km/day'
         ENDIF
         WRITE (iunp,456) 
     +        'range rate   =',sqrt(rrsum/rrwsum)*au,' km/day',
     +        'range rate wt=',sqrt(1d0/rrwsum*rrcnt)*au,' km/day'
 457     FORMAT(a,f10.5,a/a,f10.5,a/)
      ENDIF
      RETURN
      END
c =====================================================================
c OUTCOV
c =====================================================================
c  output of covariance and normal matrix, computation of eigenvalues
c   input: iun   = output unit
c          icor  = flags >0 solved for varaiable
c          gamma = covariance matrix
c          c     = normal matrix
c WARNING: if not all the variables are solved for, the corresponding
c          rows and columns of both gamma and c are set to zero.
c          These null lines are removed in the output.
      SUBROUTINE outcov(iun,icor,gamma,c)
      implicit none
c output unit, error flag
      integer iun,ierr
c which variables have been solved for
      integer icor(6),ncor
c loop indexes
      integer i,j,ii,jj
c covariance and normal matrices
      double precision gamma(6,6),c(6,6)
      double precision gam(36),cc(36),eigv(36)
c eigenvalues, eigenvectors
      double precision eigval(6),fv1(6),fv2(6)
c packing of matrices: no. lines and columns
      ncor=0
      do 1 i=1,6
        if(icor(i).gt.0)ncor=ncor+1
 1    continue
      if(ncor.eq.0)then
         write(*,*)' no correction, icor=',icor
         return
      endif
      write(iun,*)
      write(iun,*)'no. solve for=',ncor
      write(iun,*)icor
c packing of matrices:
      ii=0
      do 10 i=1,6
        jj=0
        if(icor(i).gt.0)then
          ii=ii+1
          do 11 j=1,6
            if(icor(j).gt.0)then
               jj=jj+1
               gam(ii+(jj-1)*ncor)=gamma(i,j)
               cc(ii+(jj-1)*ncor)=c(i,j)
            endif
 11       continue
        endif
 10   continue
c output covariance
      write(iun,*)
      write(iun,*) 'COVARIANCE MATRIX'
      do 2 j=1,ncor
        write(iun,109) (gam(i+(j-1)*ncor),i=1,ncor)
 109    format(6e24.16)
 2    continue
c eigenvalues
      call rs(ncor,ncor,gam,eigval,1,eigv,fv1,fv2,ierr)
      write(iun,*)
      write(iun,*) 'EIGENVALUES '
      write(iun,109) (eigval(i),i=1,ncor)
      write(iun,*)
      write(iun,*) 'EIGENVECTORS'
      do 3 j=1,ncor
        write(iun,109) (eigv(i+(j-1)*ncor),i=1,ncor)
 3    continue
c normal matrix
      write(iun,*)
      write(iun,*) 'NORMAL MATRIX'
      do 4 j=1,ncor
        write(iun,109) (cc(i+(j-1)*ncor),i=1,ncor)
 4    continue
      write(iun,*)
      return 
      end
c
* Copyright (C) 1998 by OrbFit Consortium
c ===================================================================
c DIFCOR  differential corrector
c ===================================================================
c version 1.8 Steven Chesley, Dec. 15, 1998
c        works on equinoctal elements: 
c        a,h=e sin(long.peri),k=e cos(long.peri),
c        p=tg I/2 sin(long.node),q=tg I/2 cos(long.node), lambda=mean long.
c Input: m observations number
c        w weights (only diagonal matrix)
c        sel selection flags
c        t0 epoch time for asteroid elements
c        iobs type of observations 1=alpha,delta 2=r,rdot
c        tau observations time vector
c        ioco station codes vector
c        eq0 asteroid equinoctal elements at time t0 (first guess)
c        al,de real observations vectors
c        icor(6) flags .ne.0 to correct this element, 0 to leave it as in eq0 
c        inew =2 for pseudo-Newton (use this), =1 for Newton (not well tested)
c        iunf = unit file for output; if <0, no residuals output
c        itmax=max. no iterations, if =0 then only calculate covariance...
c        itgmax= max no iterations with target function paralysed/increasing
c        delcr=control for stop due to small corrections
c        divrat=control for paralyzed target function
c Output eq corrected orbital elements at time t0
c        gamma covariance matrix
c        gtwg inverse of covariance matrix
c           warning: if only some elements are corrected, only the 
c              corresponding entries in gamma and gtwg are nonzero 
c        csinor residuals norm
c        delnor differential corrections norm
c        csir residuals (in radians) 
c                r.a. first, then declination, for each obs
c        x2 chi**2 value for each observation
c        succ logical success flag
c =============INTERFACE===== =========================================
      SUBROUTINE difcor(m,w,sel,t0,iobs,tau,ioco,eq0,al,de,icor,inew,
     +     iunf,delcr,eq,gamma,gtwg,csinor,delnor,csir,x2,succ)
c =====================================================================
      IMPLICIT NONE
c include files
      INCLUDE 'parobx.h'
      INCLUDE 'codode.h'
      INCLUDE 'trig.h'
      INCLUDE 'comdif.h'
      INCLUDE 'comrej.h'
      INCLUDE 'verbosity.h'
c ================input data==========================
c no. observations, observatory codes, selection flag, obs. type 
      INTEGER m, ioco(m), sel(m), iobs(m)
c  times, alpha, delta, 
      DOUBLE PRECISION tau(m),al(m),de(m)
c weights
      DOUBLE PRECISION w(2*m)
c unit file to write output
      INTEGER iunf
c controls 
      INTEGER inew,icor(6)
c  corr. and residuals  controls
      DOUBLE PRECISION delcr
c epoch time, initial equinoctal elements 
      DOUBLE PRECISION t0, eq0(6)
c ================output ==========================
c corrected equinoctal elements 
      DOUBLE PRECISION eq(6)
c normal and covar. matrix
      DOUBLE PRECISION gtwg(6,6),gamma(6,6)
c success flag
      LOGICAL succ
c  corr. and residuals norm
      DOUBLE PRECISION delnor,csinor
c residuals
      DOUBLE PRECISION csir(2*m),x2(m)
c =============END INTERFACE============================================
c input data sorted, in scalar form 
      DOUBLE PRECISION tsort(nobx),als(nobx),des(nobx)
      INTEGER iocos(nobx),iposs(nobx),sels(nobx),iobsrt(nobx)
c weights,no observations used (w.ne.0)
      DOUBLE PRECISION ws(nob2x),wsrej(nob2x),wsec(nob2x)
c chi**2 rms for each obs
      DOUBLE PRECISION x2s(nobx)
c residuals, condition number,  scalar residuals
      DOUBLE PRECISION csi(nob2x),ra,rd
c residuals norm: previous and first (zeroth) iterations
      DOUBLE PRECISION csino1,csinor0
c differential correction, eccentricity,scaling
      INTEGER nused,nsolv
      DOUBLE PRECISION ecc,mu,rescov
c loop variables: it iteration number, itg=iter. with increase of q,
c ittodo= itmax, but forced to 1 for itmax=0 (no correction)
      INTEGER it,itg,itrej
c loop indexes: j=1,m; i,k=1,6
      INTEGER j,k,i
c Flag for only computing covariance, not changing orbit.
      LOGICAL matonly
c unit for output
      INTEGER iun
c outlier rejection
      INTEGER nmod
c magnitude data
      include 'mag.h'
c function for control of bizarre orbit, deserving to give up
      LOGICAL bizarre
c ====================================================================
c ================== INITIALIZATION ==================================
c ====================================================================
c Are commons loaded?
      IF(iicdif.ne. 36) STOP 'difcor: internal error(1)'
      IF(iicrej.ne. 36) STOP 'difcor: internal error(2)'
c sort of times and reordering of obs. and weights
      CALL srtoss(t0,iobs,tau,al,de,ioco,sel,w,m,iobsrt,tsort,iposs,
     +   als,des,iocos,sels,ws)
c count number of solve-for variables
      nsolv=0
      DO  j=1,6
        IF(icor(j).ne.0)nsolv=nsolv+1
      ENDDO
c create weight vector with zeros for rejected obs
      DO i=1,m
         IF(sels(i).eq.0)THEN
             wsrej(2*i-1)=0
             wsrej(2*i)=0
         ELSE
             wsrej(2*i-1)=ws(2*i-1)
             wsrej(2*i)=ws(2*i)
         ENDIF
      ENDDO
c Initialisation with starting value for elements
      DO  k=1,6
        eq(k)=eq0(k)
      ENDDO
      iun=abs(iunf)
      IF(iunf.gt.0.and.verb_dif.gt.9)THEN
         write(*,220) eq
         write(iun,220) eq
 220     format(' starting values'/6f13.7)
      ENDIF
c norm of corrections for "zeroth" iteration
      csino1=1.0d25
c control for covariance only pass
      matonly=.false.
      IF(itmax.eq.0)THEN
         itmax=1
         matonly=.true.
      ENDIF
c ====================================================================
c ================== BEGIN OUTLIER LOOP ==============================
c ====================================================================
      DO itrej=1,itmaxr
c counter for increasing/paralyzed iterations
         itg=0
c +++++++++++++++++++++ Begin Inner Loop +++++++++++++++++++++++++++++
         DO it=1,itmax
c ================== single iteration ==============================
            CALL sincor(m,iobsrt,tsort,als,des,iocos,t0,icor,inew,
     +           wsrej,iun,matonly,eq,delnor,csi,csinor,gamma,gtwg)
c save first pass (sorted) data
c These will be passed out in case of hyperbolicity
            IF(it.eq.1 .and. itrej.eq.1)THEN
               CALL unsort(iposs,m,2*m,csi,csir,sel,sels,tsort,tau)
               csinor0=csinor
c cleanup the chi2
               DO j=1,m
                  x2(j)=0.d0
                  x2s(j)=0.d0
               ENDDO
            ENDIF
c for computation of covariance without correction, there is no rejection
c and no failure possible
            IF(matonly)THEN
               succ=.true.
               IF(iunf.gt.0)write(*,*)'One pass only. No corrections.'
               IF(verb_dif.ge.9)write(iun,*)
     +                 'One pass only. No corrections applied.'
               GOTO 70
            ENDIF
c control against hyperbolic and bizarre orbits
            IF(bizarre(eq))THEN
               ecc=sqrt(eq(2)**2+eq(3)**2)
               IF(verb_dif.ge.9)
     +              write(*,*)' bizarre; e=',ecc,' a=',eq(1)
               write(iun,*)' bizarre; e=',ecc,' a=',eq(1)
               succ=.false.
c     return saved rms
               csinor=csinor0
               return
            endif
c norm of the correction: the zeros DO not matter!
            IF(iunf.gt.0)THEN
               IF(verb_dif.ge.9)write(*,200)it,csinor,delnor,eq
               write(iun,200)it,csinor,delnor,eq
 200          format(' *** iteration ',i3/,'  RMS residuals =',1p,d12.4,
     +           '   norm corr =',d12.4/,'  new elem values'/0p,6f13.7/)
            ENDIF
c ====== DECISION #1 ===================================================
c Check if we need another iteration
c Small corrections to elements?
            IF(delnor.lt.delrej)THEN
               succ=.true. 
               GOTO 77
            ENDIF
c Target function increasing/paralyzed?
            IF(csinor.gt.csino1*1.1d0)THEN
               itg=itg+1
               IF(iunf.gt.0)THEN
                  IF(verb_dif.ge.9)THEN
                     write(*,*)' target function increasing '
                     write(*,*)
                  ENDIF
                  write(iun,*)' target function increasing '
               ENDIF
               succ=.false.
            ELSEIF(csinor.gt.csino1*divrat)THEN
               itg=itg+1
               IF(iunf.gt.0)THEN
                  IF(verb_dif.ge.9)THEN
                     write(*,*)' target function paralyzed '
                     write(*,*)
                  ENDIF
                  write(iun,*)' target function paralyzed '
               ENDIF
               succ=.true.
            ELSE
               itg=0
            ENDIF
            IF(itg.gt.itgmax)THEN
c              go to outlier rejection if paralyzed
               IF(succ) GOTO 77
c              otherwise (increasing) exit with failure
               GOTO 70
            ENDIF
            csino1=csinor
c +++++++++++++++++++++ End Inner Loop +++++++++++++++++++++++++++++
         ENDDO
         succ=.false.
         IF(iunf.gt.0)THEN
            IF(verb_dif.ge.9)write(*,*)' too many iterations'
            write(iun,*)' too many iterations'
         ENDIF
c           Must get computed orbit and selection flags to agree:
         if(autrej) call sincor(m,iobsrt,tsort,als,des,iocos,t0
     +        ,icor,inew,
     +        wsrej,iun,matonly,eq,delnor,csi,csinor,gamma,gtwg)
         GOTO 70
c ================== AUTOMATIC OUTLIER REJECTION =====================
 77      CONTINUE
         IF(autrej)THEN
         IF(iunf.gt.0.and.verb_dif.ge.9)write(*,*)'REJECTING...'
c call with -iun to avoid logging rejection stats
            CALL reject(iunf,csinor,sels,csi,tsort,x2s,ws,m,
     +        gamma,icor,nmod)
c update weight vector with zeros for rejected obs
            DO i=1,m
                IF(sels(i).eq.0)THEN
                    wsrej(2*i-1)=0
                    wsrej(2*i)=0
                ELSE
                    wsrej(2*i-1)=ws(2*i-1)
                    wsrej(2*i)=ws(2*i)
                ENDIF
            ENDDO
            IF(iunf.gt.0)THEN
               IF(verb_dif.ge.9)write(*,201)itrej,it,nmod
               write(iun,201)itrej,it,nmod
 201           format('Outlier rejection pass ',i3,
     +         ' after ',i3,' single differential correction passes.'/
     +              'There were ',i3,' changes.')
            ENDIF
         ELSE
            DO i=1,m
               x2s(i)=0.d0
            ENDDO
            nmod=0
         ENDIF
c go to next loop if no modifications or if not rejecting anyway
         IF(nmod.eq.0 .or. .not.autrej)THEN
            GOTO 75
         ENDIF
c ================== END OF OUTLIER LOOP ==============================
      ENDDO
c ====================================================================
c ================== BEGIN FINAL LOOP ================================
c ====================================================================
 75   CONTINUE
c      WRITE(*,*)delnor,delcr,succ
      IF(delnor.gt.delcr)THEN
         IF(iunf.gt.0)THEN
            IF(verb_dif.ge.9)write(*,*)'Final convergence loop after ',
     +           itrej,' outlier rejection loops.'
            write(iun,*)'Entering final convergence loop after ',
     +           itrej,' outlier rejection loops.'
         ENDIF
c counter for increasing/paralyzed iterations
         itg=0
         DO it=1,itmax
c ================== single iteration ==================================
            CALL sincor(m,iobsrt,tsort,als,des,iocos,t0,icor,inew,wsrej,
     +           iun,matonly,eq,delnor,csi,csinor,gamma,gtwg)
c control against hyperbolic and bizarre orbits
            IF(bizarre(eq))THEN
               ecc=sqrt(eq(2)**2+eq(3)**2)
               IF(verb_dif.ge.9)write(*,*)' bizarre; e=',ecc,' a=',eq(1)
               write(iun,*)' bizarre; e=',ecc,' a=',eq(1)
               succ=.false.
c cleanup the chi2
               DO j=1,m
                  x2(j)=0.d0
               ENDDO
c the residuals of the first iteration will be passed up
               RETURN
            ENDIF
c norm of the correction: the zeros DO not matter!
            IF(verb_dif.ge.5)write(*,200)it,csinor,delnor,eq
            write(iun,200)it,csinor,delnor,eq
c ====== DECISION # 2 ===============================================
c Check if we need another iteration
c Small corrections to elements?
            IF(delnor.lt.delcr)THEN
               IF(iunf.gt.0)THEN
                  IF(verb_dif.ge.9)THEN
                     write(*,*) 'Done. Corrections small after ',
     +                 it,' passes.'
                     write(*,*)
                  ENDIF
                  write(iun,*)'Done. Corrections small after ',
     +                 it,' passes.'
               ENDIF
               succ=.true.
               goto 70
            ENDIF
c     Target function increasing/paralyzed?
            IF(csinor.gt.csino1*1.1d0)THEN
               itg=itg+1
               IF(iunf.gt.0)THEN
                  IF(verb_dif.ge.9)THEN
                     write(*,*)' target function increasing '
                     write(*,*)
                  ENDIF
                  write(iun,*)' target function increasing '
               ENDIF
               succ=.false.
            ELSEIF(csinor.gt.csino1*divrat)THEN
               itg=itg+1
               IF(iunf.gt.0)THEN
                  IF(verb_dif.ge.9)THEN
                  write(*,*)' target function paralyzed '
                  write(*,*)
                  ENDIF
                  write(iun,*)' target function paralyzed '
               ENDIF
               succ=.true.
            ELSE
               itg=0
            ENDIF
            IF(itg.gt.itgmax)THEN
                  IF(verb_dif.ge.5)WRITE(*,*)
     +              'Done. Target funct. not decreasing after ',
     +              it,' passes.'
               IF(iunf.gt.0)write(iunf,*)
     +                 'Done. Target funct. not decreasing after ',
     +              it,' passes.'
               goto 70
            ENDIF
            csino1=csinor
c ====================== END FINAL LOOP =============================
         ENDDO
         succ=.false.
         IF(verb_dif.ge.9)write(*,*)' too many iterations'
         write(iun,*)' too many iterations'
      ENDIF
c ====================== CLEAN UP AND RETURN ========================
 70   CONTINUE
c covariance (and normal matrix) rescaling
      nused=0
      DO i=1,2*m
         IF(ws(i).GT.0.d0) nused=nused+1
      ENDDO
      IF(nsolv .gt. nused)  STOP 'difcor: internal error, nsolv>nused'
c To be done better .... maybe done
      mu=rescov(nsolv,nused,csinor)
      DO  i=1,6
         DO  j=1,6
            gamma(i,j)=gamma(i,j)*mu**2
            gtwg(i,j)=gtwg(i,j)/mu**2
         ENDDO
      ENDDO
c reordering the residuals for output
      CALL unsort(iposs,m,2*m,csi,csir,sel,sels,tsort,tau)
c unsort x2 and dmagn, dsunv,disv,phav,adotv,ddotv
c (data for magnitude computation)
      do i=1,m
         j=iposs(i)
         IF(autrej.and..not.matonly)THEN
            x2(j)=x2s(i)
         ENDIF
         dsuna(j)=dsunas(i)
         disa(j)=disas(i)
         phaa(j)=phaas(i)
         dmagn(j)=dmagns(i)
         adotmv(j)=adots(i)
         ddotmv(j)=ddots(i)
      enddo
c=====================
      return
      end
* Copyright (C) 1998 by OrbFit Consortium
* Version: December 11, 1997 Steven Chesley
c Single iteration of differential correction of orbital elements
c
c Input: m      observations number
c        iobsrt sorted obs. type 1000's=astrometry, 2000's=radar
c        tsort  sorted observation times
c        als    sorted right ascensions
c        des    sorted declination
c        iocos  sorted observatory codes
c        t0     epoch time for asteroid elements
c        icor(6) flags .ne.0 to correct this element, 0 to leave it as in eq0 
c        inew =2 for pseudo-Newton (use this), =1 for Newton (not well tested)
c        ws     sorted weight vector
c        iun    unit for output
c        matonly = .true. if only need covariance; do not update orbit
c
c Output eq     corrected orbital elements at time t0
c        delnor norm of the corrections
c        csi    residuals
c        csinor norm of residuals
c        gamma covariance matrix
c        gtwg inverse of covariance matrix
c           warning: if only some elements are corrected, only the 
c              corresponding entries in gamma and gtwg are nonzero
c
      SUBROUTINE sincor(m,iobsrt,tsort,als,des,iocos,t0,icor,inew,ws,
     +              iun,matonly,eq,delnor,csi,csinor,gamma,gtwg)
      IMPLICIT NONE

      INCLUDE 'parobx.h'
c elongation,distance to Earth, distance to Sun (to compute magnitude)
      include 'phase.h'
c magnitude data
      include 'mag.h'

      INTEGER m,iocos(m),inew,iun,iobsrt(m)
      DOUBLE PRECISION tsort(m),eq(6),t0,csi(2*m),csinor,ws(2*m)
      DOUBLE PRECISION gtwg(6,6),gamma(6,6),delnor,als(m),des(m)
      LOGICAL matonly

      INTEGER j,iocj,ider,i,icor(6),k,no
      INTEGER nused
      DOUBLE PRECISION tauj
      DOUBLE PRECISION alj,dej,cond,deq(6)
c first and second derivatives of alpha, delta w.r. to elements
      DOUBLE PRECISION dade(6),ddde(6),ddade(6,6),dddde(6,6)
c control of two body approximation (must be false)
      LOGICAL twobo
c second derivatives of observations
      DOUBLE PRECISION h(nob2x,6,6)
c DOUBLE PRECISION functions
      DOUBLE PRECISION snormd,snorm,appmag,pridif

      INCLUDE 'codode.h'

c assign ider depending from inew
      IF(inew.eq.2)THEN
        ider=1
      ELSEIF(inew.eq.1)THEN
        ider=2
      ELSE
        WRITE(*,*)'sincor: inew not correct',inew,icor
        STOP ' sincor: inew '
      ENDIF
      twobo=.false.
c Compute observations and derivatives
      CALL set_restart(.true.)
      DO 61 j=1,m
          tauj=tsort(j)
          iocj=iocos(j)
          IF(iobsrt(j)/1000.eq.1)THEN
             CALL alfdel(eq,t0,tauj,iocj,alj,dej,dade,ddde,ider,twobo,
     +                ddade,dddde)
c  compute magnitude difference (apparent minus absolute)
             dmagns(j)=appmag(0.d0,gmagc,dsun,dis,pha)
             phaas(j)=pha
             dsunas(j)=dsun
             disas(j)=dis
             adots(j)=adot
             ddots(j)=ddot
          ELSEIF(iobsrt(j)/1000.eq.2)THEN
             CALL rrdot(eq,iobsrt(j),t0,tauj,iocj,alj,dej,dade,ddde,
     +             ider,twobo)
          ELSE
             WRITE(*,*)'sincor: iobs= ',iobsrt(j), ' not known'
             STOP 'sincor: iobs'
          ENDIF
          CALL set_restart(.false.)
c Compute residuals, form matrix g
          IF(iobsrt(j)/1000.eq.1)THEN
             csi(2*j-1)=pridif(als(j),alj)
             csi(2*j)=pridif(des(j),dej)
          ELSE
             csi(2*j-1)=als(j)-alj
             csi(2*j)=des(j)-dej
          ENDIF
          DO 62 i=1,6
              g(2*j-1,i)=dade(i)
              g(2*j,i)=ddde(i)
              IF(inew.eq.1)THEN
                  DO k=1,6
                      h(2*j-1,i,k)=ddade(i,k)
                      h(2*j,i,k)=dddde(i,k)
                  ENDDO
              END IF
 62       continue
 61   continue
      CALL set_restart(.true.)
c Compute solution of linear least squares
      no=2*m
      CALL minsol(csi,no,ws,g,h,inew,icor,iun,gtwg,deq,gamma,cond)
c Norm of the residuals
      csinor=snormd(csi,ws,no,nused)
c Norm of the corrections
      delnor=snorm(deq,gtwg,6,6)
c Update solution (only solve-for variables) 
      IF(.not.matonly)THEN
          DO k=1,6
              IF(icor(k).ne.0) eq(k)=eq(k)+deq(k)
          ENDDO
      END IF

      END
c ===========================================================
c
c                     M I N S O L
c version 1.3, A. Milani, 28/6/1997
c
c least squares solver, including optional handling of second
c derivatives and normal matrix normalisation
c      INPUT :   csi   =  vector of residuals (of alpha and delta)
c                no    =  double of observations number
c                w     =  weight matrix 
c                g     =  matrix of first derivatives of residuals
c                          multiplied by -1
c                h     =  matrix of second derivatives of residuals
c                          multiplied by -1
c                inew  = if inew=2, h is not used; if inew=1, h is used
c                icor(k)=  .ne.0 to correct this element, 0 to let it as it was
c                iunf   = unit for output file
c      OUTPUT:   gtwg  =  normal matrix GtWG 
c                dx0   =  corrections vector
c                gamma =  inverse matrix
c                cond  =  conditioning number of gamma
c           warning: if only some elements are corrected, only the 
c              corresponding entries in gamma and gtwg are nonzero in output
c ===========================================================
      SUBROUTINE minsol(csi,no,w,g,h,inew,icor,iunf,
     +                  gtwg,dx0,gamma,cond)
c ===========================================================
      implicit none
c ===========================================================
      include 'parobx.h'
c dimension of elements vector
      INTEGER nd
      parameter (nd=6)
c unit for output file
      INTEGER iunf
c controls for corrections, no. of solve for variables
      INTEGER icor(nd),ndc
c number of (scalar) observations
      INTEGER no
c residuals, weights,derivatives
      DOUBLE PRECISION csi(no),w(no),g(nob2x,nd),h(nob2x,nd,nd)
c output: normal matrix, differential correction, covariance matrix
      DOUBLE PRECISION gtwg(nd,nd),dx0(nd),gamma(nd,nd)
c reduced normal matrix and covariance matrix
      DOUBLE PRECISION cr(nd,nd),gr(nd,nd)
c right hand side of the normal equation, normalisation factors
      DOUBLE PRECISION gtwcsi(nd)
c condition number, control for positive-definite
      DOUBLE PRECISION cond
c loop indexes: j,k=1,nd; i=1,no; jj,kk=1,ndc
      INTEGER j,k,i,jj,kk
c controls: newton/pseudo newton 
      INTEGER inew
c temporary variables and workspace 
      DOUBLE PRECISION tmp
      INTEGER iun
      INCLUDE 'verbosity.h'
c ===========================================================
c control on dimensioning
      if (no.gt.nob2x) THEN
         write(*,*) 'no>nob2x in minsol',no,nob2x
         stop 'minsol: no > nob2x'
      ENDIF
c ===========================================================
c normal matrix GtWG of the pseudo-Newton method
      DO 1 j=1,nd
        DO 2 k=1,nd
          gtwg(j,k)=0.d0
          IF(icor(k).ne.0.and.icor(j).ne.0)THEN
             DO  i=1,no
               gtwg(j,k)=gtwg(j,k)+g(i,j)*w(i)*g(i,k)
             ENDDO
          ENDIF
 2      CONTINUE
 1    CONTINUE
c ===========================================================
c Adding second derivatives for full Newton's method
      IF(inew.eq.1)THEN
         DO 24 j=1,nd
           DO 26 k=1,nd
             tmp=0.d0
             IF(icor(k).ne.0.and.icor(j).ne.0)THEN
                DO  i=1,no
                  tmp=tmp-csi(i)*w(i)*h(i,j,k)
                ENDDO
             ENDIF
             gtwg(j,k)=gtwg(j,k)+tmp
 26        CONTINUE
 24      CONTINUE
      ENDIF
c ===========================================================
c simplyfied version if 6 elements to be determined
      ndc=0
      DO  j=1,nd
        IF(icor(j).gt.0) ndc=ndc+1
      ENDDO
      IF(ndc.eq.6)THEN
         CALL invmat(gamma,nd,nd,gtwg,cond,iunf)
         goto 78
      ENDIF
c ===========================================================
c squeeze normal matrix into a smaller one
      jj=0
      DO 7 j=1,nd
        IF(icor(j).ne.0)THEN
           jj=jj+1
           kk=0
           DO 8 k=1,nd
             IF(icor(k).ne.0)THEN
                kk=kk+1
                cr(jj,kk)=gtwg(j,k)
             ENDIF
 8         CONTINUE
        ENDIF
 7    CONTINUE
      IF(jj.eq.kk.and.jj.ne.0)THEN
         ndc=jj
      ELSE
         write(*,*)' minsol, this should not happen ',jj,kk
         stop 'minsol: jj,kk '
      ENDIF        
c ===========================================================
c Cholewski inversion
      CALL invmat(gr,nd,ndc,cr,cond,iunf)
c ===========================================================
c  covariance matrix
c  warning: rows and columns not involved in the correction
c  are set to zero (should be infinite!)
      kk=0
      DO 14 k=1,nd
        IF(icor(k).ne.0)THEN
           kk=kk+1
           jj=0
           DO 16 j=1,nd
             IF(icor(j).ne.0)THEN
                jj=jj+1
                gamma(k,j)=gr(kk,jj)
             ELSE
                gamma(k,j)=0.d0
             ENDIF
 16        CONTINUE
        ELSE
           DO  j=1,nd
             gamma(k,j)=0.d0
           ENDDO
        ENDIF
 14   CONTINUE
c =======================
 78   continue
      IF(cond.gt.1.d12)THEN
         iun=abs(iunf)
         WRITE(iun,100)cond
 100  FORMAT('Conditioning number: ',1p,d12.5)
         IF(verb_dif.ge.9)WRITE(*,100)cond
      ENDIF
c ===========================================================
c Computation of vectors GtWcsi and DX0
c   warning: correction is zero automatically for non solved-for variables 
       DO 17 j=1,nd
         gtwcsi(j)=0.d0
         DO  i=1,no
           gtwcsi(j)=gtwcsi(j)+g(i,j)*w(i)*csi(i)
         ENDDO
 17    CONTINUE
       DO 18 k=1,nd
         dx0(k)=0.d0
         DO j=1,nd
           dx0(k)=dx0(k)+gamma(k,j)*gtwcsi(j)
         ENDDO
 18    CONTINUE
       RETURN
       END
c ================================================
c  INVMAT  ndim x ndim  matrix inversion 
c in this version it is assumed that the input matrix a is symmetric,
c definite positive;
c and the output matrix is symmetric, definite positive
c If this is not the case, a warning is issued on the standard ouput
c ===========INTERFACE=============================
      SUBROUTINE invmat(c,nx,ndim,a,cond,iunf)
      implicit none
      INTEGER nx,ndim,iunf
      DOUBLE PRECISION a(nx,ndim),c(nx,ndim),cond
c ========END INTERFACE=============================
      INTEGER nxx
c warning: here it is assumed that ndim never exceeds 6; 
c to be changed for bigger matrices
      parameter(nxx=6)
      DOUBLE PRECISION an(nxx),v(nxx)
      INTEGER i,j,indp,nb,iun
      DOUBLE PRECISION err,omax,omin,da,eps
      DOUBLE PRECISION roff
      logical sym
      INCLUDE 'verbosity.h'
c ==================================================
c check that the matrix is indeed symmetric
      sym=.true.
      iun=abs(iunf)
      eps=1.d2*roff(nb)
      DO 1 i=1,ndim
        DO 2 j=1,i-1
          da=abs((a(i,j)-a(j,i))/(a(i,j)+a(j,i)))
          IF(da.gt.100*eps)THEN
             write(iun,*)'invmat: ',i,j,a(i,j),a(j,i),da,100*eps
             IF(verb_dif.ge.5)
     +            write(*,*)'invmat: ',i,j,a(i,j),a(j,i),da,100*eps
             sym=.false.
          ENDIF
 2      CONTINUE
 1    CONTINUE
      IF(.not.sym)THEN
         write(iun,*)'invmat: input matrix not symmetric'
         IF(verb_dif.ge.5)write(*,*)'invmat: input matrix not symmetric'
      ENDIF
c ==========================================================
c Tcholewski
c ==========================================================
c normalisation of columns of matrix to be inverted
      DO  i=1,ndim
        an(i)=sqrt(abs(a(i,i)))
      ENDDO
      DO 88 i=1,ndim
        DO  j=1,ndim
          c(i,j)=a(i,j)/(an(i)*an(j))
        ENDDO
 88   CONTINUE
      DO 86 i=ndim+1,nx
        DO  j=1,ndim
          c(i,j)=0.d0
        ENDDO
 86   CONTINUE
c first Cholewsky factorisation
      err=eps
      CALL tchol(c,nx,ndim,indp,err)
      IF(indp.eq.0)THEN
c Control of conditioning number of the inverted matrix
        omax=c(1,1)
        omin=c(1,1)
        DO 5 i=2,ndim
          if (c(i,i).gt.omax) THEN
             omax=c(i,i)
          ENDIF
          if (c(i,i).lt.omin) THEN
              omin=c(i,i)
          ENDIF
 5      CONTINUE
        cond=(omax/omin)**2
        CALL inver(c,v,nx,ndim)
c unnormalize the matrix by norm of columns
        DO  i=1,ndim
          DO  j=1,ndim
            c(i,j)=c(i,j)/(an(i)*an(j))
          ENDDO
        ENDDO
      ELSE
        IF(iunf.gt.0)THEN
           write(iun,*)' matrix is not positive definite'
           write(iun,*)' pivot number ',indp,' is ',c(indp,indp)
           IF(verb_dif.gt.5)THEN
              write(*,*)' matrix is not positive definite'
              write(*,*)' pivot number ',indp,' is ',c(indp,indp)
           ENDIF
        ENDIF
      ENDIF
      RETURN
      END
c Copyright 1998, 2001 Orbfit Consortium
c Modified Dec. 2 1998 by Steve Chesley to include sorting of selection flags.
c Modified 2001 by A. Milani to muse quick sort
c SRTOSS: remake of srtarc
c Ordinamento dei tempi di osservazione di un arco rispetto a t0
c INPUT:    t0       -  Tempo centrale dell'arco
c           iobs     -  Observation type
c           toss(j)  -  Tempi di osservazione (j=1,noss)
c           al(j)    -  Ascensioni rette (j=1,noss)
c           de(j)    -  Declinazioni (j=1,noss)
c           sel(i)   -  Selection Flags
c           w(j)     -  Pesi delle osservazioni (j=1,2*noss)
c           noss     -  Numero di osservazioni 
c OUTPUT: 
c         iobsrt   -  Obs. types sorted
c         tsort(k) -  Tempi ordinati (prima i tempi >t0, in ordine crescente,
c                       poi quelli <t0, in ordine decrescente) (k=1,ntot)
c         iposs(k) -  Mappa: I_{noss} --> I_{noss}
c           als(j) -  Ascensioni rette riordinate
c           des(j) -  Declinazioni riordinate
c           sels(i)-  Sorted Selection Flags
c           ws(j)  -  Pesi delle osservazioni riordinati
      SUBROUTINE srtoss(t0,iobs,toss,al,de,ioco,sel,w,noss,iobsrt,tsort,
     +     iposs2,als,des,iocos,sels,ws)
      implicit none
      include 'parobx.h'
c observations number, counters
      integer noss,k,j,itt,ii,kk
      double precision toss(noss),al(noss),de(noss),w(2*noss)
      double precision tsort(noss),tt,t0
      integer iobs(noss),iobsrt(noss)
      double precision als(noss),des(noss),ws(2*noss)
      integer sel(noss),sels(noss)
      integer  iposs2(noss),ioco(noss),iocos(noss)
      logical change,order
      integer iposs(nobx)
c****************
c   static memory not required
c****************
c Sort: generating the indexes, toss is not changed
      CALL heapsort(toss,noss,iposs)
c Reordering of right ascensions declinations and weigths
c observatory codes, observation types and Selection Flags
c first the ones with time of observation > t0
      kk=0
      DO  ii=1,noss
         IF(toss(iposs(ii)).ge.t0)THEN
            kk=kk+1
            als(kk)=al(iposs(ii))
            des(kk)=de(iposs(ii))
            iocos(kk)=ioco(iposs(ii))
            sels(kk)=sel(iposs(ii))
            ws(2*kk-1)=w(2*iposs(ii)-1)
            ws(2*kk)=w(2*iposs(ii))
            iobsrt(kk)=iobs(iposs(ii))
            tsort(kk)=toss(iposs(ii))
            iposs2(kk)=iposs(ii)
         ENDIF
      ENDDO
c now those with time of observation < t0
      DO  ii=noss,1,-1
         IF(toss(iposs(ii)).lt.t0)THEN
            kk=kk+1
            als(kk)=al(iposs(ii))
            des(kk)=de(iposs(ii))
            iocos(kk)=ioco(iposs(ii))
            sels(kk)=sel(iposs(ii))
            ws(2*kk-1)=w(2*iposs(ii)-1)
            ws(2*kk)=w(2*iposs(ii))
            iobsrt(kk)=iobs(iposs(ii))
            tsort(kk)=toss(iposs(ii))
            iposs2(kk)=iposs(ii)
         ENDIF
      ENDDO
      return
      end
c==================================================
c================     UNSORT     ==================
c==================================================
c UNSORT: reordering of residuals in original order
c
c Input:   iposs
c          noss
c          nos2
c          csi
c          tsort
c          toss
c          sels
c
c Output:  csir
c          sel
c
c==================================================
      SUBROUTINE unsort(iposs,noss,nos2,csi,csir,sel,sels,tsort,toss)
      implicit none
      integer noss,nos2,iposs(noss),n,k
      double precision csi(nos2),csir(nos2)
      double precision tsort(noss),toss(noss),tmp
      integer sel(noss),sels(noss)
      do  n=1,noss
        k=iposs(n)
        tmp=toss(k)-tsort(n)
        if(tmp.ne.0.d0)then
           write(*,*)'sorting problem ',tmp,n,tsort(n),k,toss(k)
        endif
        sel(k)=sels(n)
        csir(2*k)=csi(2*n)
        csir(2*k-1)=csi(2*n-1)
      enddo
      return
      end
c =========================================================
c RESCOV: rescaling factor of covariance matrix
c to be used as follows:
c     gamma ---> gamma *rescov**2
c     normatr--> normatr/rescov**2
      DOUBLE PRECISION FUNCTION rescov(nsolv,nused,csinor)
      IMPLICIT NONE
c input: number of parameters to be solved, nomber of scaalr obs,
      INTEGER nsolv,nused
c norm of residuals (RMS relative to the given observatory weights)
      DOUBLE PRECISION csinor
      IF(nsolv.lt.nused)THEN
         IF(csinor.gt.1.d0)THEN
            rescov=csinor*sqrt(float(nused)/float(nused-nsolv))
         ELSE
            rescov=sqrt(float(nused)/float(nused-nsolv))
         ENDIF
      ELSE
         rescov=1.d0
      ENDIF
      RETURN
      END
* Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 3, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         D I F I N I                           *
*  *                                                               *
*  *        Initialization of options for iteration control        *
*  *                     in ROUTINE DIFCOR                         *
*  *                                                               *
*  *****************************************************************
*
      SUBROUTINE difini
      IMPLICIT NONE

* Common blocks to be initialized:
      INCLUDE 'comdif.h'

      LOGICAL found,fail1,fail

      fail=.false.

      batch=.false.
      CALL rdnlog('difcor.','batch',batch,.false.,found,
     +            fail1,fail)

      itmax=20
      CALL rdnint('difcor.','nit_max',itmax,.false.,found,
     +            fail1,fail)

      itgmax=5
      CALL rdnint('difcor.','nitg_max',itgmax,.false.,found,
     +            fail1,fail)

      divrat=0.999d0
      CALL rdnrea('difcor.','div_cntr',divrat,.false.,found,
     +            fail1,fail)

      IF(fail) STOP '**** difini: abnormal end ****'

      iicdif=36

      END
c
* Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)
*                    and Steve Chesley (chesley@@dm.unipi.it)
* Version: December 15, 1998 Steven Chesley
* hacked for radar: A. Milani, January 16, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R E J E C T                           *
*  *                                                               *
*  *            Rejection and recovery of outliers                 *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNILOG    -  File unit for logging rejections (<0 for no log)
*           CSINOR    -  Norm of residuals
*           SEL       -  Selection flag (0 = not selected)
*           CSI       -  O-C residuals
*           W         -  Observation weights
*           NOBS      -  Number of observations
*           GAMMA     -  Covariance matrix of orbital elements
*           ICOR      -  List of solve-for parameters
*
* OUTPUT:   SEL       -  Updated selection flags
*           X2        -  CHI**2 residual
*           NMOD      -  Number of modifications to the outlier list
*                        (additions or deletions)
*
* WARNING: in some cases the weights are zero, because a scalar observation
*          has been squeezed in a two-dimensional observation
*          This is the case for obs. types 3,4 (radar)

      SUBROUTINE reject(unilog,csinor,sel,csi,mjd,x2,w,
     +     nobs,gamma,icor,nmod)
      IMPLICIT NONE

      INCLUDE 'parobx.h'
c ===========units for err,pro,clo files ========
      INCLUDE 'proout.h'
      INCLUDE 'verbosity.h'
* Output of debugging information
      LOGICAL debug
      PARAMETER (debug=.false.)

      INTEGER nobs,icor(6),nmod,unilog
      INTEGER sel(nobs)
      DOUBLE PRECISION csinor,csi(2*nobs),gamma(6,6),w(2*nobs),x2(nobs)
      DOUBLE PRECISION mjd(nobs)

* NEEDED common blocks:
      INCLUDE 'comrej.h'
      INCLUDE 'codode.h'

      INTEGER iob,ipa,ipd,i,k,nsel,minobs,selp,nrej,nrec
      DOUBLE PRECISION sa,sd,sad,gamga(6),gamgd(6),covres(2,2),wres(2,2)
      DOUBLE PRECISION det,x2max
      DOUBLE PRECISION covr,wrr
c calendar date variables
      INTEGER iyear,imonth,iday
      DOUBLE PRECISION hour
      CHARACTER*16 date
      LOGICAL entopp
c fudge factor, stuck loops
      DOUBLE PRECISION fudge
      INTEGER nit,oldmod,oldnob

      SAVE nit,oldmod,oldnob
      DATA nit/0/

      IF(iicrej.NE.36) STOP '**** reject: internal error (01) ****'
c     log=.true.
c     IF(unilog.le.0) log=.false.

      nrej=0
      nrec=0
*****************************************************
* First compute chi**2 for each obs
*****************************************************
* Minimum number of observations
*     First count solved for params
      minobs=6
      DO 10 i=1,6
         IF(icor(i).NE.1) minobs=minobs-1
 10   CONTINUE
*     Next set minobs
      minobs=nint(0.5d0*minobs)
      IF(nobs.le.minobs)THEN
         IF(unilog.gt.0)WRITE(unilog,*)
     +        'No autorejection with only ',nobs,' observations.'
         IF(verb_rej.ge.5)WRITE(*,*)
     +        'No autorejection with only ',nobs,' observations.'
         nmod=0
         RETURN
      ENDIF

      x2max=0
      nsel=0
      IF(unilog.gt.0)WRITE(unilog,300) csinor
      IF(verb_rej.ge.9) WRITE(*,300) csinor
 300  FORMAT('==== OUTLIER REJECTION ======='/
     +        'Residual norm =',1P,E13.5,0P)
      DO 1 iob=1,nobs
* Pointers to RA and DEC observations
         ipa=2*iob-1
         ipd=2*iob
* Expected variance of fit residuals
         DO 2 i=1,6
            sa=0
            sd=0
            DO 3 k=1,6
               sa=sa+gamma(i,k)*g(ipa,k)
               sd=sd+gamma(i,k)*g(ipd,k)
 3          CONTINUE
            gamga(i)=sa
            gamgd(i)=sd
 2       CONTINUE
c separate case with zero weight in one of the two components
         IF(w(ipa).gt.0.d0.and.w(ipd).gt.0.d0)THEN
c both components to be used
            sa=0
            sd=0
            sad=0
            DO  i=1,6
               sa =sa +g(ipa,i)*gamga(i)
               sd =sd +g(ipd,i)*gamgd(i)
               sad=sad+g(ipa,i)*gamgd(i)
            ENDDO
            IF(sel(iob).EQ.0) THEN
               covres(1,1)=1/w(ipa)+sa
               covres(2,2)=1/w(ipd)+sd
               covres(1,2)=+sad
            ELSE
               covres(1,1)=1/w(ipa)-sa
               covres(2,2)=1/w(ipd)-sd
               covres(1,2)=-sad
            END IF
            covres(2,1)=covres(1,2)
* Chi-square value
            CALL inv22(covres,wres,det)
            IF (det.eq.0.d0)THEN
               WRITE(ierrou,*) 'WARNING: reject.f, det=0.'
               numerr=numerr+1
            ENDIF
            x2(iob)=csi(ipa)*wres(1,1)*csi(ipa)
     +           +csi(ipd)*wres(2,2)*csi(ipd)
     +           +2*csi(ipa)*wres(1,2)*csi(ipd)
            
         ELSEIF(w(ipa).eq.0.d0)THEN
c only second component to be used
            sd=0
            DO  i=1,6
               sa =sa +g(ipa,i)*gamga(i)
               sd =sd +g(ipd,i)*gamgd(i)
               sad=sad+g(ipa,i)*gamgd(i)
            ENDDO
            IF(sel(iob).EQ.0) THEN
               covr=1/w(ipd)+sd
            ELSE
               covr=1/w(ipd)-sd
            END IF
* Chi-square value
            IF (covr.eq.0.d0)THEN
               WRITE(ierrou,*) 'WARNING: reject.f, covr=0. (1)'
               numerr=numerr+1
            ENDIF
            wrr=1.d0/covr
            x2(iob)=csi(ipd)*wrr*csi(ipd)
         ELSEIF(w(ipd).eq.0.d0)THEN
c only first component to be used
            sa=0
            DO  i=1,6
               sa =sa +g(ipa,i)*gamga(i)
            ENDDO
            IF(sel(iob).EQ.0) THEN
               covr=1/w(ipa)+sa
            ELSE
               covr=1/w(ipa)-sa
            END IF
* Chi-square value
            IF (covr.eq.0.d0)THEN
               WRITE(ierrou,*) 'WARNING: reject.f, covr=0. (2)'
               numerr=numerr+1
            ENDIF
            wrr=1.d0/covr
            x2(iob)=csi(ipa)*wrr*csi(ipa)
         ELSE
            WRITE(*,*)'reject: this should not happen ',
     +           w(ipa),' ',w(ipd)
            STOP 'reject: w '
         ENDIF
         IF (x2(iob).lt.0.d0)THEN
            WRITE(ierrou,*) 'WARNING: reject.f, x2=',x2(iob)
            If(verb_rej.ge.5)WRITE(*,*) 'WARNING: reject.f, x2=',x2(iob)
            numerr=numerr+1
            x2(iob)=0.d0
         ENDIF
         IF(sel(iob).NE.0) THEN
            x2max=MAX(x2(iob),x2max)
            nsel=nsel+1
         END IF
 1    CONTINUE
*****************************************************
* Second perform tests
*****************************************************

      IF(nsel.LE.0) STOP '**** reject: internal error (02) ****'

*     For <=18 obs. or for inf. loop -> reject one at a time
      IF(nsel .le. minobs*6) THEN
         IF(unilog.gt.0) WRITE(unilog,*) 'Rejecting no more than one.'
      ELSEIF(nit.ge.4) THEN
         IF(unilog.gt.0) WRITE(unilog,*) 'Trying to get unstuck.'
      ELSE
         x2max=x2max*x2frac
      ENDIF
      IF(unilog.gt.0)WRITE(unilog,302) x2max
      IF(verb_rej.ge.9) WRITE(*,302) x2max
 302  FORMAT('Chi2 threshold =',1P,E13.5)
c This is a "fudge" factor to make it more difficult to reject when 
c there are very few observations. If there are more than ~50 then
c this has little effect.
      fudge=400.d0 * (1.2)**(-nsel)

      DO 5 iob=1,nobs
         selp=sel(iob)
         IF(selp.EQ.0) THEN
            IF(x2(iob).LE.x2rec + 0.75*fudge) THEN
               call mjddat(mjd(iob),iday,imonth,iyear,hour)
               write(date,'(i4,a1,i2.2,a1,f8.5)')
     +              iyear,'/',imonth,'/',iday+hour/24d0
               IF(verb_rej.ge.9)write(*,*) 'Recover: ',date,x2(iob)
               sel(iob)=1
               nrec=nrec+1
            END IF
         ELSE
            IF(x2(iob).GE.x2max .AND. x2(iob).GT.x2rej + fudge) THEN
c check to see if this is the only remaining obs in an opposition
               call mjddat(mjd(iob),iday,imonth,iyear,hour)
               write(date,'(i4,a1,i2.2,a1,f8.5)')
     +              iyear,'/',imonth,'/',iday+hour/24d0
               IF(verb_rej.ge.9)write(*,*) 'Reject : ',date,x2(iob)
               if(entopp(iob,sel,mjd,nobs))then
                  if(rejopp)then
                     write(*,*)
     +                '*** WARNING: Rejecting last obs. in opp. Date:',
     +                  date
                     write(ierrou,*)
     +                    'Rejecting last obs. in opp. Date:',date
                     numerr=numerr+1
                     sel(iob)=0
                     nrej=nrej+1                     
                  else
                     write(*,*)
     +            '*** WARNING: Did not reject last obs. in opp. Date:',
     +                    date
                     write(ierrou,*)
     +                    'Did not reject last obs. in opp. Date:',date
                     numerr=numerr+1
                  endif
               else
                  sel(iob)=0
                  nrej=nrej+1
               endif
            END IF
         END IF
         IF(debug.and.unilog.gt.0) 
     +        WRITE(unilog,301) iob,x2(iob),selp,sel(iob)
 5    CONTINUE
 301  FORMAT(I5,1P,E13.5,0P,2I3)
      IF(unilog.gt.0)WRITE(unilog,303) nobs,nsel,nrej,nrec
      IF(verb_rej.ge.9) WRITE(*,303) nobs,nsel,nrej,nrec
 303  FORMAT('No. of observations:'/
     +       10X,'total     =',I5/
     +       10X,'selected  =',I5/
     +       10X,'rejected  =',I5/
     +       10X,'recovered =',I5)

      nmod=nrej+nrec
c     Is it the same no. of modifications with the same number of obs?
      IF(oldmod.EQ.nmod .AND. nobs.EQ.oldnob)THEN
         nit=nit+1
      ELSE
         nit=0
      ENDIF
      oldmod=nmod
      oldnob=nobs

      END
c
c decide if index is the last selected obs in an opposition
      logical function entopp(index,sel,mjd,nobs)      
      integer index,nobs,sel(nobs),i
      double precision mjd(nobs)
      entopp=.false.
      do i=1,nobs
c return if we find a different selected obs that is close enough
         if(i.ne.index.and.
     +        abs(mjd(i)-mjd(index)).le.180d0.and.
     +        sel(i).eq.1) return
      enddo
      entopp=.true.
      return
      end
* Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 14, 1998 Steve Chesley
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R E J I N I                           *
*  *                                                               *
*  *        Initialization of options for outlier rejection        *
*  *                                                               *
*  *****************************************************************
*
      SUBROUTINE rejini
      IMPLICIT NONE

* Common blocks to be initialized:
      INCLUDE 'comrej.h'

      LOGICAL found,fail1,fail

      fail=.false.

      autrej=.true.
      CALL rdnlog('reject.','auto',autrej,.false.,found,
     +            fail1,fail)

      rejopp=.false.
      CALL rdnlog('reject.','rejopp',rejopp,.false.,found,
     +            fail1,fail)

      x2rej=8.0D0
      CALL rdnrea('reject.','chi2_reject',x2rej,.false.,found,
     +            fail1,fail)

      x2rec=7.0D0
      CALL rdnrea('reject.','chi2_recover',x2rec,.false.,found,
     +            fail1,fail)

      x2frac=0.25D0
      CALL rdnrea('reject.','chi2_frac',x2frac,.false.,found,
     +            fail1,fail)

      delrej=1.D-1
      CALL rdnrea('reject.','conv_cntr',delrej,.false.,found,
     +            fail1,fail)

      itmaxr=15
      CALL rdnint('reject.','nit_max',itmaxr,.false.,found,
     +            fail1,fail)

      fomax=50.D0
      CALL rdnrea('reject.','max_perc',fomax,.false.,found,
     +            fail1,fail)
      fomax=fomax/100

c Magnitude chi^2 values
      x2mrej=8.0D0
      CALL rdnrea('reject.','chi2_mag_rej',x2mrej,.false.,found,
     +            fail1,fail)

      x2mrec=7.0D0
      CALL rdnrea('reject.','chi2_mag_rec',x2mrec,.false.,found,
     +            fail1,fail)

      IF(fail) STOP '**** rejini: abnormal end ****'

      iicrej=36
      
      RETURN
      END
* Copyright 1999 Orbfit Consortium
* Written by Milani and Chesley
* Last Modified 1/20/99 by S. Chesley
************************************************************************
* ROUTINE MAGEST
************************************************************************
* INPUTS:  smag -     magnitude strings e.g., '15.55V'
*          rmsmag -   a priori rms of the observation (NOTE: rmsmag will 
*                        be changed to 99.99 for outlier rejection)
*          sel -      selection flag for the astrometric observation
*          m -        total number of observations
* OUTPUTS: h0 -       absolute magnitude
*          resmag -   residual of the mag observation (NOTE: If resmag =  
*                        1.d9 then the mag obs was not used in the fit.)
*          rmsh -     weighted RMS of the residuals
************************************************************************
      SUBROUTINE magest(smag,rmsmag,sel,m,h0,resmag,rmsh)
      IMPLICIT NONE
c number of observations
      INTEGER m
c magnitudes: estimated values, observed (string), a priori rms
      DOUBLE PRECISION h0
      CHARACTER*6 smag(m)
      DOUBLE PRECISION rmsmag(m)
c fit residuals, rms
      DOUBLE PRECISION resmag(m),rmsh
c selection flags, no. observations used and weights
      INTEGER sel(m)
c ===========END INTERFACE=================
c get nobx
      INCLUDE 'parobx.h'
c get x2mrej,x2mrec
      INCLUDE 'comrej.h'
c interactive?
c     INCLUDE 'comdif.h'
c observed magnitudes, mean
      DOUBLE PRECISION obsm(nobx),hmean,hsum,wsum,wrsum
      DOUBLE PRECISION tmprms(nobx),hnew,chi2
      LOGICAL avail(nobx)
      CHARACTER*1 col
      CHARACTER*5 smag5
      INTEGER numrej,numchg
c debugging:
      LOGICAL verbose
c common magnitude data
c dmagn is difference apparent minus absolute magnitude
      include 'mag.h'
c error file and number
      INCLUDE 'proout.h'
      INCLUDE 'verbosity.h'
c function to compute default mag RMS
      DOUBLE PRECISION magrms
c loop indexes
      INTEGER j,icount,nrej
c ======================================================================
      verbose = .false.
c First create a vector of appmags and rms's
      DO j=1,m
         READ(smag(j),101)col
 101     FORMAT(5x,a1)
         avail(j)=smag(j).ne.'      '.and.
     +        (col.eq.'B'.or.col.eq.'V'.or.col.eq.'R'.or.
     +        col.eq.'I'.or.col.eq.' ')
         IF(.not.(col.eq.'B'.or.col.eq.'V'.or.col.eq.'R'.or.
     +        col.eq.'I'.or.col.eq.' '))THEN
            IF(ierrou.gt.0)THEN
               WRITE(ierrou,*) 'Unknown Color:',col,' Obs #',j
               numerr=numerr+1
            ELSE
               WRITE(*,*) 'Unknown Color:',col,' Obs #',j
            ENDIF
         ENDIF
         IF(avail(j))THEN
c Get the default rms (to be used for recovery)
c The last three arguments are dummies.
            tmprms(j)=magrms(smag(j),0d0,0,'0')
c Read the magnitude and make color corrections.
            smag5=smag(j)
            READ(smag5,*,ERR=3)obsm(j)
            IF(col.eq.'V')THEN
               CONTINUE
            ELSEIF(col.eq.' ')THEN
               obsm(j)=obsm(j)-0.8d0
            ELSEIF(col.eq.'B')THEN
               obsm(j)=obsm(j)-0.8d0
            ELSEIF(col.eq.'R')THEN
               obsm(j)=obsm(j)+0.4d0
            ELSEIF(col.eq.'I')THEN
               obsm(j)=obsm(j)+0.8d0
            ENDIF
         ENDIF
      ENDDO
c Start Outlier loop
      DO nrej=1,10
         hsum=0.d0
         icount=0
         wsum=0.d0
c main loop
         DO j=1,m
            IF(avail(j))THEN
               icount=icount+1
               hsum=hsum+(obsm(j)-dmagn(j))/rmsmag(j)
               wsum=wsum+1.d0/rmsmag(j)
            ENDIF
         ENDDO
c jump out if no good observations
         IF(icount.eq.0)THEN
            If(verb_dif.ge.5)
     +           WRITE(*,*)' magnitude cannot be estimated w/o obs'
            rmsh=-1.d0
            RETURN
         ENDIF             
c Compute H, etc.
         hmean=hsum/wsum
         rmsh=0.d0
         wrsum=0.d0
         DO j=1,m
            IF(avail(j))THEN
               resmag(j)=obsm(j)-dmagn(j)-hmean
               rmsh=rmsh+(resmag(j)/rmsmag(j))**2
               wrsum=wrsum+(1.d0/rmsmag(j))**2
            ELSE
c              if data not used
               resmag(j)=1.d9
            ENDIF
         ENDDO
C Divide by icount or wsum here?
         rmsh=sqrt(rmsh/wrsum)
         IF(verbose) write(*,*)'H,RMS',hmean,rmsh
c Handle outliers in a separate loop
         numrej=0
         numchg=0
c Never reject if there are less than 5 obs.
         IF(icount.le.4) GOTO 10
         DO j=1,m
            IF(avail(j))THEN
               IF(rmsmag(j).gt.99.98d0)THEN
C              test for recovery
                  hnew=(hsum+(obsm(j)-dmagn(j))/tmprms(j))/
     +                 (wsum+1.d0/tmprms(j))
                  chi2=(((obsm(j)-dmagn(j))-hnew)/tmprms(j))**2
                  IF(chi2.lt.x2mrec)THEN
                     numchg=numchg+1
                     rmsmag(j)=tmprms(j)
                  ENDIF
                  IF(verbose) write(*,*)'REC:Hnew,chi2,rms,res',
     +                 hnew,chi2,rmsmag(j),resmag(j)
               ELSE
C              test for rejection
                  hnew=(hsum-(obsm(j)-dmagn(j))/rmsmag(j))/
     +                 (wsum-1.d0/rmsmag(j))
                  chi2=(((obsm(j)-dmagn(j))-hnew)/rmsmag(j))**2
                  IF(chi2.gt.x2mrej)THEN
                     numchg=numchg+1
                     rmsmag(j)=99.99d0
                  ENDIF
                  IF(verbose) write(*,*)'REJ:Hnew,chi2,rms,res',
     +                 hnew,chi2,rmsmag(j),resmag(j)
               ENDIF
C              count rejection
               IF(rmsmag(j).gt.99.98d0) numrej=numrej+1
            ENDIF
         ENDDO
         IF(numchg.eq.0) GOTO 10
      ENDDO
      If(verb_dif.ge.5)
     +     WRITE(*,*)'WARNING: Did not finish rejecting photometry...'
************************************************************************
 10   IF(verb_dif.ge.9)THEN
         WRITE(*,190)h0,hmean,icount,rmsh,numrej,nrej
 190     FORMAT('Absol Magnitude: Old=',f5.2,' New=',f5.2,/
     +        'no. magnitude obs.=',i4, ' with RMS ',f5.2,/
     +        'no. rejected obs.=',i4,' in ',i4,' passes.')
         WRITE(*,*)
      ENDIF
      h0=hmean
      RETURN

c parsing error
 3    If(verb_dif.ge.5)
     +     WRITE(*,*)'magest: error in parsing magn. no ',j,' ',smag(j)
      STOP
      END
c =========================================================
c MEANTI time in the middle of the observed arc
c =====================================================
      DOUBLE PRECISION FUNCTION meanti(tau,rmsa,rmsd,m)
      IMPLICIT NONE
      INTEGER m
      DOUBLE PRECISION tau(m),rmsa(m),rmsd(m)
      INTEGER i
      DOUBLE PRECISION tw,w
      meanti=0.d0
      IF(m.le.0)THEN
         WRITE(*,*)'meanti: no data '
         RETURN
      ENDIF
      tw=0.d0
      DO i=1,m
        w=1.d0/(rmsa(i)**2+rmsd(i)**2)
        meanti=meanti+tau(i)*w
        tw=tw+w
      ENDDO
      meanti=meanti/tw
      RETURN
      END
c =============================
c BIZARRE
c logical function deciding when it is a good idea to give up 
c differential corrections, both for difcor and difvin.
c uses parameters set in bizset
      LOGICAL FUNCTION bizarre(eq)
      IMPLICIT NONE
c input equinoctal elements
      DOUBLE PRECISION eq(6)
c controls
      DOUBLE PRECISION ecclim0, samin0,samax0,phmin0,ahmax0
      COMMON / bizcon /ecclim0, samin0,samax0,phmin0,ahmax0 
c local variables
      DOUBLE PRECISION ecc,a
c ==========================
      ecc=sqrt(eq(2)**2+eq(3)**2)
      a=eq(1)
      bizarre=.false.
      IF(ecc.ge.ecclim0)bizarre=.true.
      IF(a.ge.samax0.or.a.le.samin0)bizarre=.true.
      IF(a*(1.d0-ecc).le.phmin0.or.a*(1.d0+ecc).ge.ahmax0)bizarre=.true.
      RETURN
      END
c ===============================
c BIZSET
c assign/default values of bizarre orbit control parameters
      SUBROUTINE bizset(ecclim,samin,samax,phmin,ahmax)
      IMPLICIT NONE
c controls
      DOUBLE PRECISION ecclim, samin,samax,phmin,ahmax
      DOUBLE PRECISION ecclim0, samin0,samax0,phmin0,ahmax0
      COMMON / bizcon /ecclim0, samin0,samax0,phmin0,ahmax0 
c check if zero, then select default
      IF(ecclim.eq.0.d0)THEN
c exclude only almost hyperbolic
         ecclim0=0.99d0
      ELSE
         ecclim0=ecclim
      ENDIF
      IF(samin.eq.0.d0)THEN
c exclude inner Vulcanoids
         samin=0.1d0
      ELSE
         samin0=samin
      ENDIF
      IF(samax.eq.0.d0)THEN
c exclude long period comets
         samax0=200.d0
      ELSE
         samax0=samax
      ENDIF
      IF(phmin.eq.0.d0)THEN
c exclude sun dipping comets
         phmin0=0.001d0
      ELSE
         phmin0=phmin
      ENDIF
      IF(ahmax.eq.0.d0)THEN
c exclude long periodic comets
         ahmax0=400.d0
      ELSE
         ahmax0=ahmax
      ENDIF
      RETURN
      END
