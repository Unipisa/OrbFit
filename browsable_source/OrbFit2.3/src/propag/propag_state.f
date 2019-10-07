c==========MODULE PROPAGATOR==============================
c PUBLIC ROUTINES
c              set_restart   to operate on restart control
c              propag        general purpose propagator
c MODULE CONTAINS:
c ROUTINES
c                inivar
c                varwra
c              propin
c              rkimp
c                fct
c              rkstep
c              bdnstev
c              catst
c                kintrp
c              legnum
c              inipro
c              selste
c                compco
c                zed
c                bessel
c              ra15v
c                rasust
c                rapred
c                rabeta
c                bintrp
c                radcon
c                invaxv
c                vawrxv
c                varunw
c              propa2
c
c  HEADERS
c propag_state.o: model.h proout.h 
c
c               restart.h private
c               dx0de.h  with close_app.f
c               nvarx.h        "
c               parint.h       "
c               comint.h       "
c               rkcoef.h       "
c               closta.h       "
c               sunmass.h, trig.h public
c               iclap.h only for iclap
c               parbep.h masses.h only for names
c
c               mscoef.h  internal
c 
c
c
c ==========================================================
c SET_RESTART
c with argument logical, true or false
c if true forces restart anyway
c if false, restart is avoided it it is a continuation orbit
c as controlled in propag
c ========================================================
      SUBROUTINE set_restart(res_log)
      IMPLICIT NONE
      LOGICAL res_log
      INCLUDE 'restart.h'
      restar=res_log
      RETURN
      END
c last update  March 1999 by A. Milani
c =====================================================================
c PROPAG vers. 1.9
c =====================================================================
c N+1-body problem propagator
c
c This version uses JPL ephemerides as source for the planetary positions
c
c  WARNING: the input elements east and the output cartesian
c           coordinates and derivatives are in ecliptic (mean of J2000.0)
c           coordinates.
c =====================================================================
c
c  input: t0 epoch time
c         t2 prediction time
c         east equinoctal orbital elements vector a,h,k,p,q,lambda 
c               at time t0
c         ider 0=no partials 1= partials dx/deast (3x6 matrix)
c              2=second derivatives (3x6x6 tensor) 2-body approximation
c  output:
c        xast position and velocity vector in heliocentric cartesian 
c               coordinates for the asteroid at time t1
c        xea  position and velocity vector in heliocentric cartesian 
c               coordinates for the Earth at time t1
c        dxde first derivatives of position vector with respect to elements 
c
c        ddxde second derivatives of position vector with respect to elem.
c
c This routine cannot compute n-body second derivatives: if ider.ge.2,
c the matrix ddxde is computed by prop2b. 
c
c ================INTERFACE===========================================
      SUBROUTINE propag(t0,east,t2,xast,xea,ider,dxde,ddxde)
      implicit none
c ===============INPUT==============================
c elements (equinoctal), epoch of the elements (MJD), target epoch
      double precision east(6),t0,t2
c order of derivatives 0,1,2
      integer ider
c ===============OUTPUT=============================
c cartesia coord., at time t2, of the asteroid, of the Earth, derivatives 
      double precision xast(6),xea(6),dxde(6,6),ddxde(3,6,6)
c ===========END INTERFACE===========================================
c ========INCLUDE HEADERS===================
c options common
      include 'model.h'
      include 'iclap.h'
      include 'comint.h'
c ======== constant of gravitation ==============
      include 'sunmass.h'
c =============STATE VECTOR AND DERIVATIVES======================
c main state vector for only one asteroid + variational eq.
      INCLUDE 'nvarx.h'
      double precision y1(nvarx),y2(nvarx)
c cartesian pos. of asteroid, 2-body approximation
      double precision xast2b(6)
c matrices of partial derivatives
c derivatives with respect to elements stored for strclan
      INCLUDE 'dx0de.h'
      double precision dxdx0(6,6),ddummy(6,6)
c =============RESTART CONTROL=========================
c times: current, epoch at previous call
      double precision t1,t0s
c asteroid elements at the previous call
      double precision easts(6)
c restart control
      integer nfl
c dirty trick to force restart anyway
      include 'restart.h'
c startup from a non-close approaching state at each integration restart
      INCLUDE 'closta.h'
c close approach at initial time monitoring
      INTEGER iplam
      DOUBLE PRECISION texit
c ===================================================
c stepsize: as given by selste, current, previous
      double precision hgiv,h,hs
c flag for interpolation of position only or pos and vel.
      integer istate
      common/cstate/istate
c =============
c integers for loop control
      integer i
c integers for dimensions
      integer nv,nv1,nvar,nvar2
c other reals: eccentricity, mean motion 
      double precision ecc,enne
c dirty trick to pass ider to right hand side
      integer ide
      common/deriv/ide
c initialisation control (to have dummy initial conditions)
      integer lflag
* **************************************
c static memory allocation
      save
c
      data lflag/0/
c **************************************
c  dummy initial conditions to force restart the first time
      if(lflag.eq.0)then
         lflag=1      
c Store fictitious epoch time and elements for the asteroid
c (to be sure they are different the first time)
         t0s=-1d+55
         do 72 i=1,6
 72        easts(i)=0.d0
c masses, distance control and so on needed the first time
         call masjpl
      endif
c =====================================================================
c JPL Earth vector at observation time
      CALL earcar(t2,xea,1)
c =====================================================================
c Check if time and/or elements changed: we could not need to compute 
c initial conditions for the asteroid
c     WRITE(*,*)t0,t1,t2,restar
      if(t0.ne.t0s.or.easts(1).ne.east(1).or.easts(2).ne.east(2).
     +   or.easts(3).ne.east(3).or.easts(4).ne.east(4).
     +   or.easts(5).ne.east(5).or.easts(6).ne.east(6)
     +   .or.abs(t2-t1).gt.abs(t2-t0).or.restar
     +   )then
c     WRITE(*,*)east

*  new arc; propin needs to be informed
         nfl=0
c also cloapp needs to be informed
         clost=.true.
         CALL clotest(t0,east,iplam,texit)  
c         IF(iplam.ne.0)WRITE(*,*)'initial close app. with plan. ',iplam 
c =====================================================================
c Compute asteroid cartesian elements at time $t0$
         call prop2b(t0,east,t0,xast,gms,1,dx0de,ddxde)
         DO  i=1,6
           easts(i)=east(i)
         ENDDO
* =====================================================================
c choices about integration method and stepsize
         if(imet.eq.0) then
c automatical choise of numerical integration method
           call selmet(east)
         else
            icmet=imet
         endif
c read masses from JPL header (stored in a common...)
* ****** moved March 12, 1999 ***********
c alignement of masses might have changed, because of different list of 
c asteroids (to avoid self perturbation) 
         call masjpl
* ****** added on Sat Jun 14 1997 ******
c control if both position and vel. are needed in dpleph
         if(icmet.eq.3.and.iclap.eq.1) then
           istate=2
         elseif(irel.eq.1)then
           istate=2
         else
           istate=1
         endif
* ***************************************
* selection of an appropriate stepsize 
         if(icmet.ne.3)then
c for the multistep:
c the starter uses the same stepsize
            ecc=sqrt(east(2)**2+east(3)**2) 
            enne=sqrt(gms/east(1)**3)
            call selste(ecc,enne,error,mms,hms,hgiv)
            h=hgiv
         elseif(icmet.eq.3)then
c for Everhart: step selection is automatic, hev is the maximum
            h=hev
         endif
c sign control:
         if(t2-t0.lt.0.d0)then
            h=-h
         endif
c =====================================================================
c position vector dimension
         nv=3
c derivatives required?
         ide=ider
         if(ider.ge.1)then
c Variational equations: 3 X 6 components
            nv1=18
         else
            nv1=0
         endif
c total number of variables (and total number of positions)
         nvar2=nv+nv1
         nvar=nvar2*2
c =====================================================================
c Vector y1, length 6, contains at least
c the asteroid position and velocity 
         DO  i=1,3
            y1(i)=xast(i)
            y1(i+nvar2)=xast(i+3)
         ENDDO
         if(ider.ge.1)then
c in this case the vector y1 contain also the matrices 
c of partial derivatives; then we need to
c initialise the  variation matrix as the 6 X 6 identity
            CALL inivar(y1,nvar2,nvar)
         endif
c =====================================================================
c Restart from $t0$
         t0s=t0
         t1=t0
      else
c old arc:
c stepsize was already decided
         if(icmet.eq.3)then
            h=hev
         else
            h=hgiv
         endif
c sign control:
         if(t2-t1.lt.0.d0)then
            h=-h
         endif
c need to check that direction is not changed
         if(h*hs.lt.0.d0)then
c restart is necessary because of U-turn
            nfl=0
         else
c restart not necessary; propin can go on 
            nfl=1
         endif
      endif
c =====================================================================
c propagator to compute asteroid and planets orbits
      call propin(nfl,y1,t1,t2,y2,h,nvar2,nvar)
      hs=h
c =====================================================================
c Asteroid coordinates
      do  i=1,3
          xast(i)=y2(i)
          xast(i+3)=y2(i+nvar2)
      enddo
      if(ider.ge.1)then
c if partial derivative are required
c First parzial derivatives: rewrap vector into 6x6 matrix
         CALL varwra(y2,dxdx0,nvar,nvar2)
c =====================================================================
c Chain rule: we compute dxde=dxdx0*dx0de
         call mulmat(dxdx0,6,6,dx0de,6,6,dxde)
c second derivatives can be computed only in 2 body approximation
         if(ider.ge.2)then
            call prop2b(t0,east,t1,xast2b,gms,2,ddummy,ddxde)
         endif
      endif
      return
      end
c =====================================================================
c inivar
c
c initialise the  variation matrix as the 6 X 6 identity
c ================================
      SUBROUTINE inivar(y1,nvar2,nvar)
      IMPLICIT NONE
      INTEGER nvar,nvar2
      DOUBLE PRECISION y1(nvar)
c end interfface
      INTEGER i,j,iii,ij        
      iii=3
      do 7 j=1,6
         do  i=1,3
            ij=i+3*(j-1)
            y1(iii+ij)=0.d0
            y1(iii+ij+nvar2)=0.d0
            if(i.eq.j)then
               y1(iii+ij)=1.d0
            elseif(j.eq.i+3)then
               y1(iii+ij+nvar2)=1.d0
            endif
         enddo
 7    continue
      RETURN
      END   
c ====================================================
c varwar
c
c First parzial derivatives: rewrap vector into 6x6 matrix
c ====================================================
      SUBROUTINE varwra(y2,dxdx0,nvar,nvar2)
      IMPLICIT NONE
      INTEGER nvar,nvar2
      DOUBLE PRECISION dxdx0(6,6),y2(nvar)
      INTEGER i,j,ij
c ====================================================
      DO j=1,3
         DO  i=1,3
            ij=i+3*(j-1)+3
            dxdx0(i,j)=y2(ij)
            dxdx0(i,j+3)=y2(ij+9)
            dxdx0(i+3,j)=y2(ij+nvar2)
            dxdx0(i+3,j+3)=y2(ij+9+nvar2)
         ENDDO 
      ENDDO
      RETURN
      END
c **********************************************************
c  PROPIN   ORBFIT version, Jun 13 1997
c  purpose: propagator/interpolator; the state vector y1
c           at time t1 is propagated to y2 at time t2
c  input:
c      nfl: if 0, integration needs to be restarted
c      t1: initial time
c      y1: state vector at t1
c      nfl: flag $>0$ if continuing, $=0$ if it has to restart
c      t2: final time
c      h: stepsize (only if it is fixed; ra15 uses hev as initial step)
c      nvar2: nvar/2
c      nvar: length(y1)
c  output:
c      t2: time to which y2 refers (different by .lt.deltos from input)
c      y2: state vector at time t2
c      t1,y1: time and state vector after last integration step
c     
c  osservazione: se nelle chiamate successive alla prima t1 ed y1 non
c                sono stati modificati dal programma chiamante,l'integra
c                zione prosegue a partire da t1,y1 usando i dati dei pas
c                si precedenti accumulati nelle cataste ck (per rkimp)
c                e dd, delta (per il multistep)
c                se invece t1,y1 sono stati modificati, il che risulta
c                da nfl=0,l'integrazione riprende con un passo iniziale.
c                se pero' t1,y1 non sono modificati,i tempi t2 devono es
c                sere in successione monotona consistente con h
c  subroutine chiamate:
c       ra15: metodo di everhart
c       force: external per chiamata secondo membro 
c       legnum :lettura coefficienti r-k
c       rkimp :propagatore runge-kutta
c       kintrp :interpolatore di ck
c       bdnste :multistep BDS che fa nstep passi
c       catst :caricatore catasta differenze e somme
c       coeff :calcolo coefficienti multistep
c       rkstep :variatore automatico del passo
c       fct :secondo membro equazioni differenziali (fisso; riduce
c            al primo ordine quanto fornito da force)
c   codice indici : i=1,nvar, j=1,isrk
c **********************************************************
      SUBROUTINE propin(nfl,y1,t1,t2,y2,h,nvar2,nvar)
      implicit none
c headers
      INCLUDE 'model.h'
      INCLUDE 'nvarx.h'
      INCLUDE 'parint.h'
      INCLUDE 'comint.h'
      INCLUDE 'proout.h'
      INCLUDE 'mscoef.h'
c state vectors,times, stepsize
      double precision y1(nvarx),y2(nvarx),t1,t2,h,h2
c  workspace per rkimp
      double precision ck(ismax,nvarx),ck1(ismax,nvarx),ck2(ismax,nvarx)
      double precision ep(itmax),dery(nvarx)
c  workspace for bdstep
      double precision dd(nvar2x,4),delta(nvar2x,m2max)
c  close approaching planet
      double precision xxpla(6)
c internal times
      double precision t0,tint,din
c controls
      double precision sdelto,epin
      integer nfl
c dimensions
      integer nvar,nvar2,ndim
c counters
      integer nstep,nrk,npas
c integers
      integer lflag, idc,j,j1,n,lf,it,lit,i,m
c type of equation 
      integer nclass
c arrays to store state transition matrix
      DOUBLE PRECISION stm(6,6),stmout(6,6),stm0(6,6),tcur
c name of right hand side routine
      external force
* **************************************
c static memory allocation
      save
c **********************************************************
c  controllo metodi propagazione:
c           imet=1 (multistep) =2 (runge-kutta) =3 (everhart)
c    runge-kutta:
c           isrk=ordine del runge-kutta (diviso 2)
c           h=passo di integrazione
c           eprk=controllo di convergenza nell'equaz. implicita del rk
c           lit1,lit2=iterazioni di gauss-seidel al primo passo e dopo
c           iusci=flag uscita controlli numerici
c    multistep:
c           mms=ordine del multistep ('m' nell'art. Cel.Mech)
c nnnn      ipc=iteraz correttore; epms=controllo
c nnnn      iusms=flag uscita controlli numerici
c    everhart:
c nnnn      isrk=ordine (per ora solo 15)
c           h=passo (fisso se ll.le.0)
c           ll=controllo 10**(-ll); se ll.gt.0, 
c              scelta automatica del passo
c           iusci=flag uscita contr. num.; 
c                  se iusci.ge.0, uscita passo cambiato
      data lflag/0 /
c *******************************************************
      IF(abs(t1-t2).lt.deltos)THEN
c safety against rounding off accumulation (if times are equal)
c        t2=t1
         DO  i=1,nvar
              y2(i)=y1(i)
         ENDDO
         RETURN
      ENDIF         
      if(icmet.eq.3)then
c  everhart method (propagates to exact time):
         if(icrel.gt.0)then
            nclass=2
         else
            nclass=-2
         endif
c accumulate state transition matrix
         IF(nvar.gt.6)THEN 
            CALL vawrxv(y1,y1(nvar2+1),stm0,nvar2)
            CALL invaxv(y1,y1(nvar2+1),nvar2)
         ENDIF
 666     CONTINUE
         call ra15(y1,y1(nvar2+1),t1,t2,tcur,nvar2,nclass,idc)
         IF(tcur.eq.t2)THEN
c  current time and state is the final one;
            IF(nvar.gt.6)THEN 
c  but the accumulated state
c  transition matrix stm0 has to be used: stmout=stm*stm0
               CALL vawrxv(y1,y1(nvar2+1),stm,nvar2)
               CALL mulmat(stm,6,6,stm0,6,6,stmout)
               CALL varunw(stmout,y1,y1(nvar2+1),nvar2)
            ENDIF
            do i=1,nvar
               y2(i)=y1(i)
            enddo
            t1=t2
            return
         ELSE
c write message
c               WRITE(*,*)'end close approach to planet',idc,tcur
c propagation has been interrupted because of a close approach
            IF(nvar.gt.6)THEN
c setup the close approach record with derivatives
               CALL strclan3(stm0)
c multiply the accumulated state transition matrix
               CALL vawrxv(y1,y1(nvar2+1),stm,nvar2)
               CALL mulmat(stm,6,6,stm0,6,6,stmout)
               CALL mcopy(6,6,stmout,stm0)
               CALL invaxv(y1,y1(nvar2+1),nvar2)
            ENDIF
            t1=tcur
            GOTO 666
         ENDIF
      endif
c **********************************************************
c set m (order of multistep method)
      m=mms
      if(lflag.eq.0)then
c  inizializzazione
c  contapassi:npas caricamento delta
c  nrk caricamento ck
         npas=0
         t0=t1
         nrk=0
c  dimensione del sec membro eq ordine 2
         nvar2=nvar/2
c  dimensione usata per controllo di convergenza (no variational eq.)
         ndim=3
c  inizializzazioni per il multistep
         j1=0
         h2=h*h
c  fine inizializzazioni
         lflag=1
      else
c  controllo se nuovo arco
         if(nfl.eq.0)then
            npas=0
            t0=t1
            nrk=0
c  inizializzazioni per il multistep
            j1=0
            h2=h*h
         endif
      endif
*******************************************************
c  controllo direzione del tempo
      if((t2-t1)*h.lt.0.d0) then
          write(*,999)t1,t2,h
 999      format('propin: from t1=',f12.4,' to t2=',f12.4,
     $         ' with step h=',d12.4)
          write(*,*)'propin: this should not happen'
          STOP
c          h=-h
      endif
c **********************************************************
c  main loop:
c  controllo di funzione: propagatore o interpolatore
  5   if(dabs(t2-t1).lt.dabs(h)-deltos)goto 90
c
c  scelta propagatore
      if(npas.lt.m.or.icmet.ge.2)then
c **********************************************************
c  runge kutta implicito fa un solo passo
c
 24      continue
c  passo iniziale?
         if(nrk.le.0)then
c
c  passo iniziale: store secondo membro
          call force(y1,y1(nvar2+1),t1,delta(1,m+1),nvar2,idc,xxpla,0,1)
c close approach control
            CALL clocms(idc,t1,xxpla)
c  passo iniziale: inizializzazione di ck a zero
            do 11 j=1,isrk
            do 11 i=1,nvar
 11           ck(j,i)=0.d0
            lit=lit1
         else
c  passo non iniziale : interpolazione dei ck
            do 21 j=1,isrk
            do 21 n=1,nvar
 21           ck1(j,n)=ck(j,n)
            call kintrp(ck1,ck,isrk,nvar)
            lit=lit2
         endif
c  un passo del rk
 22      call rkimp(t1,h,y1,dery,ck,isrk,y2,lit,force,
     +            nvar,eprk,ep,lf,ndim)
c  controllo di avvenuta convergenza
         if(lf.le.0)then
c  caso di non convergenza
c           call camrk(ep,npas,nrk,lf)
            CALL rkstep(ep,npas,nrk,lf,h)
            h2=h*h
            if(lf.eq.2)goto 24
            goto 22
         endif
c  passo del rk adottato
         npas=npas+1
         nrk=nrk+1
         t1=t0+h*npas
         do 25 i=1,nvar
 25        y1(i)=y2(i)
         if(iusci.gt.0)then
            it=iabs(lf)
            if(npas.eq.iusci*(npas/iusci))write(ipirip,1000)npas,
     $                (ep(i),i=1,it)
 1000       format(' npas',i6,' ep ',5d12.3/(5d12.3))
         endif
         if(icmet.lt.2)then
c  preparazione per il multistep
            call force(y1,y1(nvar2+1),t1,delta(1,m+1-npas),nvar2,
     +            idc,xxpla,0,1)
c close approach control
            CALL clocms(idc,t1,xxpla)
c store in differences array
            if(npas.eq.m)call catst(m,m+1,b,a,nvar2,nvar2x,nvar,
     $                              delta,dd,y1,h,h2)
         endif
         goto 5
c ***************************************************************
c  propagatore multistep fa nstep passi, l'ultimo con la velocita'
      else
c  occorrono altri passi; quanti?
         sdelto=deltos*dabs(h)/h
         din=(t2-t1+sdelto)/h
         nstep=din
c
c  propagazione per nstep passi
c
c  versione predittore soltanto
c
 32      call bdnste(t1,y1,h,h2,nstep,m,j1,dd,delta,nvar2,nvar2x,
     $                 nvar,force)
c  passo del ms adottato
         npas=npas+nstep
         t1=t0+h*npas
      endif
c  nuovo passo
      goto 5
c ***************************************************************
c
c  interpolazione
 90   continue
c  controllo se occorre interpolare
      if(dabs(t2-t1).le.deltos)then
c  t2 e' circa t1 : non occorre interpolare
c        t2=t1
         DO n=1,nvar
           y2(n)=y1(n)
         ENDDO
      else
c  occorre interpolare
c  con il propagatore runge-kutta usato con passo variato
         do 92 j=1,isrk
         do 92 n=1,nvar
 92         ck2(j,n)=0.d0
         lit=lit1
         epin=eprk
         tint=t2-t1
 94      call rkimp(t1,tint,y1,dery,ck2,isrk,y2,lit,force,
     +   nvar,epin,ep,lf,ndim)
         if(lf.le.0)then
c  interpolatore impazzito
            it=iabs(lf)
            write(ipirip,1002)tint,epin,(ep(j),j=1,it)
 1002       format( ' interpolatore impazzito,passo h=',
     $      f10.7,' controllo =',d12.3/' ep= ',5 d12.3/(5d12.3/))
            epin=epin*10.d0
            goto 94
         endif
      endif
*************************************************
      return
      end
c ***************************************************************
c  RKIMP ORBFIT
c
c  scopo : propagatore che compie un passo con
c       il metodo di runge kutta implicito (formula di gauss)
c       ad is passi. il numero di iterazioni non deve superare
c       lit .  fct e' il secondo membro, che riduce al primo
c       ordine il vero primo membro fct2
c  input:
c      t1 : tempo iniziale
c      h : lunghezza del passo da eseguire
c      y1(nvar) : vettore di stato al tempo t1
c      ck(ismax,nvar) : valori iniziali per l'equazione implicita
c      epsi : controllo di convergenza per l'equazione implicita
c      ndim : no. variabili da usare per la norma da cfr. con epsi
c  output :
c      t1: invariato
c      y3(nvar) : stato al tempo t1+h (non puo' avere lo stesso
c               indirizzo di y1)
c      ep(i) : controlli numerici;per i=1,lit converg.in
c              gauss-seidel
c      dery(nvar) : spazio di lavoro per il calcolo del secondo membro
c      lf : flag che e' >0 se c'e' stata soddisfacente convergenza
c              altrimenti segnala problemi
c  codice indici: i=1,nvar; j=1,is; id=1,ndim; it=indice di iterazione
      SUBROUTINE rkimp(t1,h,y1,dery,ck,is,y3,lit,fct2,
     +     nvar,epsi,ep,lfleps,ndim)
      IMPLICIT NONE
c coefficients
      INCLUDE 'parint.h'
      INCLUDE 'rkcoef.h'
c  dimensioni variabili
      integer nvar
      double precision y1(nvar),dery(nvar),y3(nvar)
      double precision ep(itmax),ck(ismax,nvar),t(ismax)
c name of right hand side
      external fct2
c  tempo, passo
      double precision t1,h
c  indici di iterazione, di sottopasso, di dimensione, flag
      integer it,lit,j,i,is,jj,ndim,id,lfleps
c  controlli
      double precision epsi
c controllo memoria
      INTEGER ips,imem      
c  temporanei
      double precision de
c****************
c   static memory not required
c****************
c ===============================================
      do 8 j=1,is
 8      t(j)=t1+h*c(j)
      do 7 it=1,itmax
 7      ep(it)=0.d0
c
c  gauss-seidel per i ck
c  inizio iterazioni per i ck
      ips=0
      it=1
c  main loop
 1    do 11 j=1,is
        do 12 i=1,nvar
          de=0.d0
          do 13 jj=1,is
 13         de=de+a(j,jj)*ck(jj,i)
 12       y3(i)=de*h+y1(i)
        imem=j
        call fct(t(j),y3,dery,nvar,fct2,ips,imem)
        do 14 i=1,ndim
 14        ep(it)=ep(it)+dabs(dery(i)-ck(j,i))
        do 15 id=1,nvar
 15        ck(j,id)=dery(id)
 11   continue
c  controllo se le iterazioni g-s sono finite
      ep(it)=ep(it)/is
      lfleps=it
      if(ep(it).gt.epsi)then
         if(it.ge.lit)then
c  troppe iterazioni in gauss-seidel
c  il nuovo valore non viene calcolato
            lfleps=-it
            return
          else
            it=it+1
            ips=-1
            goto 1
          endif
      endif
c
c  calcolo nuovo punto y3
      do 41 i=1,nvar
        de=0.d0
        do 41 j=1,is
          de=de+b(j)*ck(j,i)
 41       y3(i)=y1(i)+h*de
      return
      end
c ===========================================================
c   FCT
c   right hand side routine
c   reduces equations to order 1
c   starting from accelerations computed by force
      SUBROUTINE fct(t,y,dery,nvar,fct2,ips,imem)
      implicit none
      integer nvar,nvar2,idc,i
      double precision y(nvar),dery(nvar)
      double precision xxpla(6)
      double precision t
      INTEGER ips,imem
c name of right hand side
      external fct2
      INCLUDE 'proout.h'
      INCLUDE 'parbep.h'
      INCLUDE 'masses.h'
      INCLUDE 'iclap.h'
c****************
c   static memory not required
c****************
      nvar2=nvar/2
      call fct2(y,y(nvar2+1),t,dery(nvar2+1),nvar2,idc,xxpla,ips,imem)
      IF(iorb.eq.11)THEN
         if(iclap.ne.0.and.idc.ne.0)then
* to be improved with a real close approach subroutine like closapp/falsi
            write(*,*)'t =',t,' close approach to planet=',
     +           ordnam(idc)
            write(iuncla,*)'t =',t,' close approach to planet=',
     +           ordnam(idc)
         endif
      ELSEIF(iorb.eq.9)THEN 
         if(idc.ne.0)then
            write(*,*)'t =',t,' close approach code=',idc
            write(iuncla,*)'t =',t,' close approach code =',idc
         endif
      ENDIF
      do i=1,nvar2
        dery(i)=y(nvar2+i)
      enddo
      return
      end
c ==========================================================
c RKSTEP
c automated RKG stepsize change
c ==========================================================
      SUBROUTINE rkstep(ep,npas,nrk,lf,h)
      IMPLICIT NONE
c     INCLUDE 'model.h'
      INCLUDE 'comint.h'
      INCLUDE 'parint.h'
      INCLUDE 'proout.h'
      DOUBLE PRECISION ep(itmax),h
      INTEGER npas,nrk,lf,l,i
c  printout data on aborted step
      l=iabs(lf)
      write (*,*) 'Non-convergence in rk-gauss. See .pro file.'
      write (ipirip,100)npas,nrk,lf,ep(l),eprk,h,isrk,(ep(i),i=1,l)
 100  format(' non convergence in rk at step ',i4,' nrk ',i3,' lf= ',i5
     */' last control ',d14.5,' convergence required ',d12.3/
     *' stepsize ',d14.6,'  order 2*',i2/
     *' controls ',5d12.3/(5d12.3/))
c change stepsize and retry
      h=0.8d0*h
c with h changed, kintrp can not be used;
c also the catatst prepared for the mulktistep needs to be redone
      nrk=0
      lf=2
      npas=0
      RETURN
      END
c ******************************************************************
c  b  d  n  s  t  e
c
c  propagatore multistep
c  bds ordine m+2
c  esegue nstep passi, all'ultimo calcola anche le velocita'
c  con aggiornamento catasta a flip-flop
c
c  input:
c           h=passo; h2=h*h
c           m=max ord diff da tenere (predittore ha ord.m+3)
c           c=Cow pred b=Cow corr f=Ad pred a=Ad corr
c           nvar=dim di y1; nvar2=nvar/2
c           fct2=nome secondo membro (forma con solo der. 2e)
c  input e output:
c           t1=tempo (esce invariato)
c           y1=vett stato pos+vel (esce aggiornato se converge)
c           j1=0,1 indice flip=flop per l'indirizzamento
c           delta=differenze; dd=somme
c
c  versione predittore soltanto
      SUBROUTINE bdnste(t1,y1,h,h2,nstep,m,j1,dd,delta,nvar2,nvar2x,
     $                  nvar,fct2)
      implicit none
c header files
      include 'parint.h'
      include 'mscoef.h'
      include 'model.h'
*   number of steps, order-2,flipflop control,
      integer nstep,m,j1
*   time, stepsize
      double precision t1,tt,h,h2
c name of right hand side
      external fct2
*   workspace
      integer nvar,nvar2,nvar2x
      double precision delta(nvar2x,m2max),dd(nvar2x,4),y1(nvar)
*   close approaching planet  
      integer idc
      double precision xxpla(6)
*   indexes for address computations
      integer id1,id2,in1,in2,kmax,kmin,kj1,kj2,kj2p
*   loop indexes
      integer i,k,n     
*   scalar temporary
      double precision dtemp
* **************************************
c static memory not required
* **************************************
c
c  loop sul numero di passi
      tt=t1
      do 1 n=1,nstep
c
c  indirizzamento controllato dal flip-flop j1
      if(j1.eq.0)then
         id2=2
         id1=1
         in2=4
         in1=3
         kmax=m+1
         kmin=1
         kj1=0
         kj2=m+1
      else
         id2=4
         id1=3
         in2=2
         in1=1
         kmax=2*m+2
         kmin=m+2
         kj1=m+1
         kj2=0
      endif
c
c  predittore
c
      if(m.eq.10)then
         do  i=1,nvar2
           y1(i)=h2*(delta(i,kmax)*c(kmax)+delta(i,kmax-1)*c(kmax-1)+
     +         delta(i,kmax-2)*c(kmax-2)+delta(i,kmax-3)*c(kmax-3)+
     +         delta(i,kmax-4)*c(kmax-4)+delta(i,kmax-5)*c(kmax-5)+
     +         delta(i,kmax-6)*c(kmax-6)+delta(i,kmax-7)*c(kmax-7)+
     +         delta(i,kmax-8)*c(kmax-8)+delta(i,kmax-9)*c(kmax-9)+
     +         delta(i,kmax-10)*c(kmax-10)+dd(i,id2))
         enddo
      elseif(m.eq.6)then
         do  i=1,nvar2
           y1(i)=h2*(delta(i,kmax)*c(kmax)+delta(i,kmax-1)*c(kmax-1)+
     +         delta(i,kmax-2)*c(kmax-2)+delta(i,kmax-3)*c(kmax-3)+
     +         delta(i,kmax-4)*c(kmax-4)+delta(i,kmax-5)*c(kmax-5)+
     +         delta(i,kmax-6)*c(kmax-6)+dd(i,id2))
         enddo
      else
         do 11 i=1,nvar2
          dtemp=0.d0
          do 12 k=kmax,kmin,-1
 12         dtemp=dtemp+delta(i,k)*c(k)
 11         y1(i)=(dtemp+dd(i,id2))*h2
      endif
      tt=tt+h
      IF(icrel.gt.0)THEN
         do  i=1,nvar2
            dtemp=0.d0
            do k=kmax,kmin,-1
               dtemp=dtemp+delta(i,k)*f(k)
            enddo
            y1(i+nvar2)=h*(dtemp+dd(i,id1))
         enddo
      ENDIF
c
c   accelerazione nel nuovo punto
      kj2p=kj2+1
      call fct2(y1,y1(nvar2+1),tt,delta(1,kj2p),nvar2,idc,xxpla,0,1)
c ==============================
c close approach control
      CALL clocms(idc,tt,xxpla)
c ===============================
c   aggiornamento catasta
      do 15 i=1,nvar2
       dd(i,in1)=dd(i,id1)+delta(i,kj2p)
       dd(i,in2)=dd(i,id2)+dd(i,in1)
 15   continue
      do 16 k=1,m
        do 16 i=1,nvar2
 16       delta(i,k+1+kj2)=delta(i,k+kj2)-delta(i,k+kj1)
c ============================================
c   passo concluso
      j1=1-j1
c   aggiungere un correttore per controllo?
c
c   fine loop sul numero di passi
 1    continue
c   calcolo velocita', se non gia' fatto
      IF(icrel.eq.0)THEN
         do  i=1,nvar2
            dtemp=0.d0
            do k=kmax,kmin,-1
               dtemp=dtemp+delta(i,k)*f(k)
            enddo
            y1(i+nvar2)=h*(dtemp+dd(i,id1))
         enddo
      ENDIF
      return
      end
c *******************************************************************
c  {\bf catst} ORB8V
c   riempimento tavola differenze
      SUBROUTINE catst(m,m1,b,a,nvar2,nvar2x,nvar,delta,dd,y1,h,h2)
      implicit none
c order, number of variables (positions only)
      integer m,m1,nvar2,nvar2x,nvar
c multistep coefficients, table of differences, sums, state vector
      double precision b(m1),a(m1),delta(nvar2x,m1),
     +        dd(nvar2x,2),y1(nvar)
c loop indexes
      integer i,k,l
      double precision h,h2,dtemp
c****************
c   static memory not required
c****************
      do 2 i=1,nvar2
        do 2 k=1,m
          do 2 l=1,m-k+1
 2        delta(i,m+2-l)=delta(i,m+1-l)-delta(i,m+2-l)
c   prima e seconda somma determinate dalle condizioni iniziali
      do 3 i=1,nvar2
        dtemp=y1(i+nvar2)/h
        do 4 k=0,m
 4        dtemp=dtemp-a(k+1)*delta(i,k+1)
        dd(i,1)=dtemp
        dtemp=y1(i)/h2+dtemp
        do 5 k=0,m
 5        dtemp=dtemp-b(k+1)*delta(i,k+1)
        dd(i,2)=dtemp
 3    continue
      return
      end
c ***************************************************************
c  {\bf  kintrp} ORB8V
c
c  scopo : predizione dei valori di ck per interpolazione da quelli
c          del passo precedente (con polinomio di grado is-1)
c          serve per avere condizioni iniziali vicine al punto unito
c          nel procedimento iterativo per risolvere le equazioni
c          implicite la cui soluzione e' la matrice ck.
c  input :
c       ck1(ismax,nvar) : ck al passo precedente
c       is : numero di stadi del metodo rk
c       nvar : numero di variabili nell'equazione di moto
c  output:
c       ck(ismax,nvar) : valori interpolati
c  osservazione : non va usato nel passo iniziale e se rispetto al
c                 passo precedente e' cambiato il passo, oppure is.
c****************
c   static memory not required
c****************
      SUBROUTINE kintrp(ck1,ck,is,nvar)
      implicit none
      include 'parint.h'
c coefficients
      INCLUDE 'rkcoef.h'
c  dimensioni dipendenti da ismax (qui=ismax)
      integer nvar
      double precision ck1(ismax,nvar),ck(ismax,nvar)
c  numero di step intermedi
      integer is
c  indici di loop
      integer j,n,jj
c  temporanei
      double precision de
      do 1 j=1,is
        do 1 n=1,nvar
        de=0.d0
        do 2 jj=1,is
 2        de=de+a1(j,jj)*ck1(jj,n)
 1      ck(j,n)=de
      return
      end
c ***************************************************************
c  LEGNUM
c ***************************************************************
c reads  Runge--Kutta--gauss coefficients 
c to be read in file ./lib/rk.coef
c  is=required number of substeps; isfl=0 if found, otherwise
c  uses closest available
c ==============INTERFACE===========================
      SUBROUTINE legnum(is,isfl)
      implicit none
c =========INPUT============
      integer is
c ========OUTPUT============
      integer isfl
c ========HEADERS===========
      include 'parint.h'
c ======================================
c  coefficients RKG for rkimp
      INCLUDE 'rkcoef.h'
c input unit, current is
      integer iun,is1
c loop indexes
      integer i,j,jj
c skip trick
      character*1 cc(ismax)
c****************
c   static memory not required
c****************
c reads RKG coefficients rk
      isfl=-1
      call filopl(iun,'rk.coe')
 198  read(iun,100,end=199)is1
 100  format(6x,i4)
c  control on is compatible wtih parameter ismax
      if(is1.gt.ismax.or.is1.le.0)goto 199
      if(is1.eq.is)then
         read(iun,101)(c(j),j=1,is1)
 101     format(7x,5d24.16)
         read(iun,102)(b(j),j=1,is1)
 102     format(7x,5d24.16)
         do 103 j=1,is1
           read(iun,104)(a(i,j),i=1,is1)
 103     continue
 104     format(3x,5d24.16)
         do 105 j=1,is1
           read(iun,106)(a1(i,j),i=1,is1)
 106       format(4x,5d24.16)
 105     continue
         isfl=0
         goto 199
       else
         read(iun,201)(cc(j),j=1,is1)
 201     format(7x,5(23x,a1))
         read(iun,201)(cc(j),j=1,is1)
         do 203 j=1,is1
           read(iun,204)(cc(jj),jj=1,is1)
 204     format(3x,5(23x,a1))
 203     continue
         do 205 j=1,is1
           read(iun,206)(cc(jj),jj=1,is1)
 206     format(4x,5(23x,a1))
 205     continue
         isfl=is1
         goto 198
      endif
c end read
 199  call filclo(iun,' ')  
      return
      end
c ================================================================
c INIPRO
c ================================================================
c this subroutine reads from file 'propag.def' all the propagator options;
c  WARNING: we advise the user against changing the file propag.def;
c                do it at your risk !....
c  options:
c  control of propagation methods:
c           imet =1 (multistep) =2 (runge-kutta) =3 (everhart)
c             =0 automatic (multistep for main belt, Everhart for high
c                 eccentricity and/or planet crossing)
c    Runge-Kutta-Gauss: (imet=2; also used as starter for imet=1)
c           isrk = order of the method is isrk*2
c           h = integration step
c           eprk = convergence control in the solution of the implicit eq.
c           lit1,lit2 = no. of Gauss-Seidel iterations (first step, afterwards)
c 
c    multistep:
c           mms = multistep number of previous steps
c         WARNING: order is mms+2; e.g in Milani and Nobili, 1988, m=mms+2
c 
c    everhart:
c           h = stepsize (fixed if  llev.le.0; initial if llev.gt.0))
c           llev = control is 10**(-llev)
c 
c           iusci = control of output numerical parameters
c ================================================================
      SUBROUTINE inipro
* ************************************
      implicit none
* model parameters
      include 'model.h'
      include 'parint.h'
      include 'comint.h'
      include 'proout.h'
      include 'mscoef.h'
      logical fail,fail1,found
      double precision h
      integer iork,iork_c,iord,isfl
      integer j
c****************
c   static memory not required (used only once)
c****************
      fail=.false.
*   read options from already assigned input channel
      call rdnint('propag.','imet',imet,.true.,found,fail1,fail)
      call rdnint('propag.','llev',llev,.true.,found,fail1,fail)
      call rdnrea('propag.','hev',hev,.true.,found,fail1,fail)
      call rdnrea('propag.','h',h,.true.,found,fail1,fail)
      hms=h
      call rdnrea('propag.','deltos',deltos,.true.,found,fail1,fail)
      call rdnrea('propag.','error',error,.true.,found,fail1,fail)
      call rdnint('propag.','iord',iord,.true.,found,fail1,fail)
* iord= multistep order; mms=number of back steps required to start
      mms=iord-2
      call rdnrea('propag.','epms',epms,.true.,found,fail1,fail)
c options for RKGauss
      call rdnint('propag.','iork',iork,.true.,found,fail1,fail)
      isrk=iork/2
      call rdnrea('propag.','eprk',eprk,.true.,found,fail1,fail)
      call rdnint('propag.','lit1',lit1,.true.,found,fail1,fail)
      call rdnint('propag.','lit2',lit2,.true.,found,fail1,fail)
c options for Radau
      call rdnint('propag.','lit1_r',lit1_r,.true.,found,fail1,fail)
      call rdnint('propag.','lit2_r',lit2_r,.true.,found,fail1,fail)
c options for both Radau and RKGauss during close approach
      call rdnrea('propag.','eprk_c',eprk_c,.true.,found,fail1,fail)
      call rdnint('propag.','iork_c',iork_c,.true.,found,fail1,fail)
      isrk_c=iork_c/2
      call rdnint('propag.','lit1_c',lit1_c,.true.,found,fail1,fail)
      call rdnint('propag.','lit2_c',lit2_c,.true.,found,fail1,fail)
c interactive/verbosity options, obsolescent
      call rdnint('propag.','iusci',iusci,.true.,found,fail1,fail)
      call rdnint('propag.','icha',icha,.true.,found,fail1,fail)
      if(fail)stop '**** inipro: abnormal end ****'
*   initialize multistep and implicit Runge-Kutta-Gauss
c ********************************************************************
c  input coefficients for Runge-Kutta
      if(imet.ne.3)then
         call legnum(isrk,isfl)
         if(isfl.ne.0)then
            write(*,997)isrk,isfl
 997        format(' required isrk=',i4,'  found only up to ',i4)
            stop
         endif
      endif
c ********************************************************************
c   calcolo coefficienti multistep
c   cambio notazione-nell'input iord=m nell'art. cel.mech.
c   d'ora in poi mms come in revtst, orbit8a
      if(mms.gt.mmax)then
         write(ipirip,998)mms,mmax
 998     format(' chiesto mms=',i4,'  spazio solo per ',i4)
         stop 998
      endif
c   calcolo coefficienti predittore
c   c=Cow pred f=Cow corr b=Ad pred a=Ad corr
c   warning: one order more than used, because c(m+2) is required
c   by the error formula
      call compco(mms+2,c,f,b,a)
      cerr=c(mms+2)
c  duplicazione per ridurre i calcoli di indici nel multistep
      do  j=1,mms+1
           c(j+mms+1)=c(j)
           f(j+mms+1)=f(j)
c  warning--i coeff per il correttore non sono pronti(dupl+primo)
      enddo
      return
      end
* last update Fri Feb 21 1997
* =====================================================================
* SELSTE: selection of the stepsize for multistep
* =====================================================================
*  input: 
*         ecc eccentricity
*         enne mean motion
*  output: 
*         h recommended stepsize for multistep
* =====================================================================
      SUBROUTINE selste(ecc,enne,error,mms,hmax,h)
      implicit none
* input (hms max stepsize from common)(mms multistep order from common)
      INCLUDE 'parint.h'
      INCLUDE 'trig.h'
      INCLUDE 'mscoef.h'
      INCLUDE 'proout.h'
c input/output
      double precision ecc,enne,error,hmax,h
      integer mms
c scalar temporaries
      double precision eps,econv,hh,z,emul,step,err
      integer ila,igr,nb
c functions
      double precision roff
c***************
c   static memory not required
c****************
* =====================================================================
* the truncation error will be of the order of error
* per revolution squared; however, error cannot be less than
* machine accuracy, otherwise the rounding off would be the dominant
* source of integration error.
      eps=roff(nb)
      eps=max(eps,error)
      hh=hmax
* =====================================================================
* for each orbit, compute the truncation error in longitude
* according to Milani and Nobili, 1988, Celest. Mech. 43, 1--34
* and select a stepsize giving a truncation  error equivalent to
* one roundoff per revolution squared 
      econv=1.d0
c      write(ipirip,*)' Selection of the best stepsize'
      call  zed(ecc,mms+2,z,econv,ila,igr)
* if the series with Bessel functions is divergent, e.g. for e=0.2
* with iord=12, some wild guess is used because the error
* estimation formula is no good anyway; hope in your luck to have
* a not too bad orbit
*     if(ila.ge.20)z=1.d6
      if(mod(mms,2).eq.0)then
              emul=dpig**2*3/2*cerr*z
              step=(eps/emul)**(1.d0/(mms+3))/enne
              err=emul*(enne*hh)**(mms+3)
      else
              emul=dpig**2*3/4*mms*cerr*z
              step=(eps/emul)**(1.d0/(mms+4))/enne
              err=emul*(enne*hh)**(mms+4)
      endif
c     write(ipirip,*)enne,ecc,z
* =====================================================================
* now compare with h given in input
      h=min(step,hmax)
      write(ipirip,110)hmax,h
 110  format(' max. step required ',f9.6,' selected ',f9.6)
      return
      end
c **********************************************************
c   {\bf compco} ORB8V
c
c   calcolo coefficienti multistep
c   cs=Cow pred bs=Cow corr fs=Ad pred as=Ad corr
c   versione con calcoli real*8 (trasportabile)
c   per calcoli in quadrupla precisione cambiare:
c      - tutti i .d0 in .q0
c      - real*16 a(mmax) etc
c   m1=m+1
c   i coefficienti per le differenze fino
c   all'ordine iord-2=m sono calcolati e restituiti
c   negli array con indice da 1 a m+1
      SUBROUTINE compco(m1,cs,fs,bs,as)
      implicit none
      include 'parint.h'
      integer m,m1,i,j,i1,i2,k
      double precision a(mmax),b(mmax),c(mmax),d(mmax)
      double precision as(m1),bs(m1),cs(m1),fs(m1)
c****************
c   static memory not required (used only once)
c****************
      m=m1-1
      a(1)=1.d0
      do 10 i=2,m+3
        k=i-1
        a(i)=0.d0
        do 10 j=1,k
 10       a(i)=a(i)-a(j)/(k-j+2.d0)
      do 11 j=1,m+3
        b(j)=0.d0
        do 11 k=1,j
 11       b(j)=b(j)+a(j+1-k)*a(k)
      do 12 i=1,m+3
        c(i)=0.d0
        d(i)=0.d0
        do 12 j=1,i
          c(i)=c(i)+a(j)
 12       d(i)=d(i)+b(j)
      do 17 i=2,m+2
        i1=i-1
        as(i1)=a(i)
 17     fs(i1)=c(i)
      do 18 i=3,m+3
        i2=i-2
        bs(i2)=b(i)
 18     cs(i2)=d(i)
c     write(6,200)(cs(i),fs(i),bs(i),as(i),i=1,m+1)
c200  format(4f19.16)
      return
      end
c**********************************************************************
c  {\bf zed}  ORB8V
c Calcolo di Z(m,e) -- correzione all'errore di troncamento
c dovuta alla eccentricita'. m e' l'ordine del multistep nella
c forma Stormer predictor; in LONGSTOP m=12.
c eps controllo di convergenza.
c
c Routine che calcola Z(m,e);imax numero massimo iterazioni (dato 20)
c**********************************************************************
      SUBROUTINE zed(e,m,f,eps,i,igr)
      implicit double precision (a-h,o-z)
      INCLUDE 'proout.h'
      data imax/20/
      f=0.d0
      adf=0.d0
      if(e.gt.0.2d0)then
         f=1.d6
         write(ipirip,*)' with e>0.2 select the stepsize by hand'
         return
      elseif(e.gt.0.3d0)then
         f=1.d6
         write(ipirip,*)' with e>0.3 you should not use the multistep'
         return
      endif
      do 1 i=1,imax
        x=i*e
        ri=i
        cie=(bessel(i-1,x)-bessel(i+1,x))/i
        sie=(bessel(i-1,x)+bessel(i+1,x))/i
        c2ie=cie*cie
        s2ie=sie*sie
        sum=(c2ie+s2ie)/2.d0
        df=sum*ri**(m+4)
        f=f+df
cc      write(*,100)i,df
cc100   format(i5,d18.6)
        vadf=adf
        adf=dabs(df)
        if(adf.gt.vadf)igr=i
        if(adf.lt.eps.and.i.gt.1)return
 1    continue
      write(ipirip,101)i-1,df,igr
 101  format(' non convergence in zed; last term=',i5,d18.6/
     +       ' max was for ',i4)
      f=1.d6
      return
      end
c**********************************************************************
c  {\bf bessel}  ORB8V
c Function che calcola la funzione di Bessel J(i,x).
c imax numero massimo di iterazioni (dato 20)
c**********************************************************************
      double precision function bessel(i,x)
      implicit double precision (a-h,o-z)
      INCLUDE 'proout.h'
      data imax/20/
      data epbs/1.d-10/
      x2=x/2.d0
      x2p=x2
      ifat=1
      do 1 j=2,i
        x2p=x2p*x2
 1      ifat=ifat*j
      if(i.eq.0)x2p=1.d0
      jfat=1
      jpifat=1
      ifl=1
      bessel=0.d0
      do 2 j=1,imax
        dbess=x2p*ifl/ifat
        dbess=dbess/jfat
        dbess=dbess/jpifat
cc      write(*,*)dbess
        bessel=bessel+dbess
cc      write(*,*)i,j,ifat,jfat,jpifat
        if(dabs(dbess).lt.epbs)return
        jpifat=jpifat*(i+j)
        jfat=jfat*j
        ifl=-ifl
        x2p=x2p*x2*x2
 2    continue
      write(ipirip,101)j-1,dbess
 101  format('bessel: non convergence; last term=',i5,d18.6)
      return
      end
c ========================================================
c RA15V
c Vers 2.29 last update 9 January 2002 A. Milani
c with partial recomputation of right hand side.
c
c structured version; subroutines are in rasubs.f
c
c  integrator radau by e. everhart, physics department, university of denver
c  this 15th-order version, called ra15, is written out for faster execution.
c  y'=f(y,t) is  nclass=1,  y"=f(y,t) is nclass= -2,  y"=f(y',y,t) is nclass=2
c  tfin is t(final); tini is t(initial).
c  nv = the number of simultaneous differential equations.
c  the dimensioning below assumes nv will not be larger than 60.
c  ll controls sequence size. thus ss=10**(-ll) controls the size of a term.
c  a typical ll-value is in the range 6 to 12 for this order 15 program.
c  however, if ll.lt.0 then xl is the constant sequence size used.
c  x and v enter as the starting position-velocity vector, at time tini,
c  and are output as the final position-velocity vector at time tfin.
c =========================INTERFACE===========================
      SUBROUTINE ra15(x,v,tini,tfin,tcur,nv,nclass,idcend)
      IMPLICIT NONE
c =================INPUT============================
c initial and final time, current time at output 
c (in case propagation is interrupted)
      DOUBLE PRECISION tini,tfin,tcur
c dimension of state variables
      INTEGER nv
c error control
      INTEGER ll
c initial stepsize
      DOUBLE PRECISION xl
c equation type
      INTEGER nclass
c =================OUPUT============================
c state variables
      INCLUDE 'nvarx.h' 
      DOUBLE PRECISION x(nvar2x),v(nvar2x)
c close approach end flag
      INTEGER idcend
c =====================END INTERFACE===========================
c total number of function calls
      INTEGER nf
c state variables (WARNING: there is a fixed limit here!)
      INTEGER ndx,nvx
      PARAMETER (ndx=60,nvx=2*ndx)
c right hand side at begin/end of step
      DOUBLE PRECISION f1(nvx)
c memory control flag
      INTEGER ips,imem
c storage arrays
      DOUBLE PRECISION b(7,nvx),g(7,nvx),e(7,nvx),bd(7,nvx)
c ===============CLOSE APPROACH CONTROL====================
c headers
      INCLUDE 'iclap.h'
      INCLUDE 'proout.h'
      INTEGER idc
c positon and velocities of planet which has close-encounter
c whith the asteroid
      DOUBLE PRECISION xpla(6)
      LOGICAL cloend
c logical flag for variational eq., state transition matrix accumulation
      LOGICAL variaz
c ========================================================
c constants of the integration method
      DOUBLE PRECISION h(8),w(7),u(7),c(21),d(21),r(21),w1
      COMMON/raconst/h,w,u,c,d,r,w1
c ========================================================
c logical control
      LOGICAL npq,nsf,nper,ncl,nes
c ========================================================
c time and stepsize
c initial stepsize with sign,time direction, final and current time from tini 
      DOUBLE PRECISION xldir,dir,tf,tm
c stepsize control: first guess, current stepsize, t or t**2
      DOUBLE PRECISION tp,t,t2
c error control, abs. val. stepsize, truncation error parameter 
      DOUBLE PRECISION ss,tval,hv
c error control: convergence
      INCLUDE 'parint.h'
      INCLUDE 'comint.h'
      DOUBLE PRECISION ep(itmax)
c intermediate times
      DOUBLE PRECISION s,q
c scalar temporaries
      DOUBLE PRECISION gk,temp
c loop indexes k=1,nv; l=1,7; j=2,8
      INTEGER k,l,j
c iteration count for implicit equations m=1,ni 
      INTEGER m
c iteration number control
      INTEGER ni
c count of steps done
      INTEGER ns
c count of step reductions
      INTEGER ncount
c ========================================================
c scalar variables
c constants
      DOUBLE PRECISION pw
c static memory allocation
      SAVE
c ===============to remove ORBFIT commons========================
c     INTEGER lit1,lit2,itmax,iusci,ipirip
c     DOUBLE PRECISION eprk
c     parameter (itmax=20)
c     lit1=10
c     lit2=4
c     ipirip=10
c     iusci=0
c     eprk=1.d-12
c     iclap=0
c ========================================================
c begin execution
c ========================================================
c setup of logical controls
c  y'=f(y,t)  ncl=.true.    y"=f(y,t)  ncl=.false.   y"=f(y',y,t) ncl=.false.
c  nclass=1   npq=.true.    nclass= -2 npq=.true.    nclass= 2    npq=.false.
      npq=nclass.lt.2
      ncl=nclass.eq.1
c set to zero storage arrays (for first order equations only)
      DO  k=1,nv
        IF(ncl) v(k)=0.d0
      ENDDO
c  nper is .true. only on last sequence of the integration.
      nper=.false.
c  nsf is .false. on starting sequence, otherwise .true.
      nsf=.false.
c  nes is .true. only if ll is negative. then the sequence size is xl.
      ll=llev
      xl=hev
      nes=ll.lt.0
c variaz is true if the state transition matrix is being propagated;
      variaz=nv.gt.3
c ===============================================================
c  evaluate the constants in the h-, w-, u-, c-, d-, and r-vectors
c ===============================================================
      CALL radcon(ncl)
      pw=1.d0/9.d0
c ==========================================================
c initialisation
c ==========================================================
c direction of time
      dir=1.d0
      tf=tfin-tini
      if(tf.lt.0.d0) dir=-1.d0
      xldir=dir*abs(xl)
      ss=10.**(-ll)
c ===============================================================
c  the statements above are used only once in an integration to set up the
c  constants.  the  next set in
c  a reasonable estimate to tp based on experience. same sign as dir.
c  an initial first sequence size can be set with xl even with ll positive.
c ===============================================================
      IF(xldir.ne.0.d0)then
        tp=xldir
      ELSE
         IF(nes)then
            write(*,*)' ra15v: fixed stepsize xl=0; ll=',ll
            stop
         ELSE
            tp=0.1d0*dir
         ENDIF
      ENDIF
      IF(tp/tf.gt.1.d0)then
         tp=tf
      ELSEIF(tp/tf.gt.0.5d0)then
         tp=0.5d0*tf
      ENDIF
      ncount=0
c ============================================================
c information on initial state of the integrator
      IF(iusci.gt.100)write(ipirip,999)tini,tfin,tp*dir
 999  format(' ra15: from tini=',f12.4,' to tfin=',f12.4,
     $   ' with max. step h=',1p,d12.4)
c ===============================================================
c  line 4000 is the starting place of the first sequence.
c ===============================================================
 4000 ns=0
      nf=0
c ===============================================================
c force initial value for first step
c ================================================================
      tm=0.d0
      ips=0
      imem=1
      CALL force(x, v, tm+tini, f1,nv,idc,xpla,ips,imem)
      nf=nf+1
c
 3000 CONTINUE
c set to zero storage arrays
c when extrapolation from previous step is not to be used
      do  k=1,nv
        do  l=1,7
          e(l,k)=0.d0
          b(l,k)=0.d0
          bd(l,k)=0.d0
        enddo
      enddo
c ===============================================================
c line 722  begins every step after the first.
 722  CONTINUE
c ===============================================================
c compute array g from b
      CALL rabeta(nv,b,d,g)
c ===============================================================
c set time variables
      t=tp
      t2=t*t
      IF(ncl) t2=t
      tval=abs(t)
c ===============================================================
c  loop 175 is lit1_r iterations on first step and lit2_r iterations therafter.
      IF(nes)THEN
c problem: if the use of fixed step size is forced by setting in input ll<0
c (in the option file, propag.llev ,0), then the numebr of iterations
c used for close approaches extends to the entire integration. We leave the
c option for distinction of the first step from the following one, although this
c implies taht lit1_c has two uses, one for fixed step radau and one
c for single step RKG used by falsi.
         IF(nsf)THEN
            ni=lit1_c
         ELSE
            ni=lit2_c
         ENDIF
      ELSE
         IF(nsf)THEN
            ni=lit2_r
         ELSE
            ni=lit1_r
         ENDIF
      ENDIF
      do 175 m=1,ni
c =========do on substep========================= 
        CALL rasust(m,t,t2,tm,tini,x,v,b,f1,nv,ncl,npq,g,ep(m),nf)
c =====================================================
c  iteration of step is over.
c =====================================================
c iteration control based on norm of differences
        IF(m.eq.1)goto 175
        IF(ep(m)/ep(1).lt.eprk)then
c prepare for stepsize determination of next step
           IF(.not.nes)then
c if variable stepsize, step size control
              hv=0.d0
              do 635 k=1,nv
 635            hv=dmax1(hv,abs(b(7,k)))
              hv=hv*w(7)/tval**7
           ENDIF
           IF(iusci.gt.1000)write(ipirip,996)ni,tini+tm+t,ep(1)
     +            ,(ep(j)/ep(1),j=2,m)
 996       format('ra15v: good convergence iter ',i3,' time=',f10.2,
     +    ' controls:'/(5d12.4/)) 
           goto 176
        ENDIF
c =========end do on iterations========================= 
 175  continue
c ======================================================
c bad convergence of iterations
      IF(iusci.gt.0)write(ipirip,997)ni,eprk,ep(1)
     +            ,(ep(j)/ep(1),j=2,ni)
 997  format('ra15v: bad convergence iter ',i3,' eprk=',d10.2,
     +    ' controls:'/(5d12.4/))   
      IF(nes)then
         write(*,*)tm,ep
         write(*,*)' ra15v: non convergence with fixed step ',t 
c         stop
      ELSE
         tp=.8d0*tp
         goto 3000
      ENDIF
c ======================================================
c satisfactory convergence; continue as in old version
 176  IF (.not.nsf)then
c ======================================================
c block executed only for first step
c ======================================================
         IF(nes)then
c fixed stepsize
            tp=xldir
         ELSE
c variable stepsize: compute initial stepsize
            tp=(ss/hv)**pw*dir
c not too quick increase 
            IF(tp/t.gt.1.4d0) tp=t*1.4d0
c not too quick decrease
            IF(tp/t.lt.1.d0)then
               tp=.8d0*tp
               ncount=ncount+1
               IF(ncount.gt.20) then
c too many shortenings of the stepsize
                  write(ipirip,*)' ra15v: ncount.gt.10',xl,tp
                  stop
               ENDIF
        IF(ncount.gt.1.and.iusci.gt.100)write(ipirip,888)ncount,t,tp
 888           format (2x,2i3,2d18.10)
c  restart with 0.8x sequence size if new size called for is smaller than
c  originally chosen starting sequence size on first sequence.
               go to 4000
            ELSE
c  initial step accepted
               IF(iusci.gt.100)write(ipirip,*)'ra15v: second step=',tp 
            ENDIF
         ENDIF
         nsf=.true.
c ======================================================
         IF(iusci.gt.100)write(ipirip,*)' tp ',tp
      ENDIF
c =====================================================
      CALL rapred(ncl,x,v,t,t2,f1,b,nv)
      DO k=1,nv
        IF(abs(x(k)).gt.1.d15)THEN
            WRITE(*,*)' rapred ',k,x(k),v(k),f1(k)
        ENDIF
        IF(abs(v(k)).gt.1.d15)THEN
            WRITE(*,*)' rapred ',k,x(k),v(k),f1(k)
        ENDIF
      ENDDO
c =====================================================
c current time update
      tm=tm+t
      ns=ns+1
c  return if done.
      IF(nper.or.abs(t-tf).lt.1.d-10)then
         tcur=tfin
         IF(iusci.gt.100)write(ipirip,*) 'nf,ns,tfin ',nf,ns,tfin
         return
      ENDIF
c  control on size of next sequence and adjust last sequence to exactly
c  cover the integration span. nper=.true. set on last sequence.
c ===============================================================
c force initial value for next step
c ================================================================
      idcend=idc
      ips=0
      imem=1
      CALL force(x,v,tm+tini,f1,nv,idc,xpla,ips,imem)
      nf=nf+1
c ================================================================
c close approach control
      IF(iclap.eq.1)THEN
c problem: if the integration begins inside a clsoe approach, then the first step
c is done with variable stepsize; will the stepsize adjust automatically
c to a value clsoe to the true-anomaly controlled one???
         CALL cloapp(tm+tini,t,x,v,nv,idc,xpla,xldir,dir,nes,cloend)
         IF(cloend.and.variaz)THEN
            tcur=tm+tini
            return
         ELSE
c            idcend=0
         ENDIF
      ENDIF
c ================================================================
c check if the integration required is finished
      IF(nes)THEN
         tp=xldir
      ELSE
         tp=dir*(ss/hv)**pw
         IF(tp/t.gt.1.4d0) tp=t*1.4d0
c         WRITE(*,*)' next step ',tp
         IF(iusci.gt.100)write(ipirip,*)'step=',tp
      ENDIF
      IF(dir*(tm+tp).gt.dir*tf-1.d-8)THEN
         tp=tf-tm
         nper=.true.
      ENDIF
c ================================================================
c  new step, no extrapolation, recompute right hand side
      IF(.not.nsf) GOTO 4000
c  if the integration continues without discontinuities, 
c  extrapolate b-values for next step.
      q=tp/t
      CALL bintrp(q,b,e,bd,nv,ns)
c  new step
c     IF(.not.nsf) GOTO 3000
      go to 722
      end
c ================================================
c rasust
c
c iterative convergence of the right hand sides at intermediate points
c ================================================
      SUBROUTINE rasust(m,t,t2,tm,tini,x,v,b,f1,nv,ncl,npq
     +        ,g,epsi,nf)
      IMPLICIT NONE
c iteration number
      INTEGER m
c control to be used for convergence
      DOUBLE PRECISION epsi
c dimension of state vector
      INTEGER nv
c logical flag for first order diff. eq., for second order diff.eq.
      LOGICAL ncl,npq
c
      DOUBLE PRECISION x(nv),v(nv)
c right hand side at begining of step
      DOUBLE PRECISION f1(nv)
c storage arrays
      DOUBLE PRECISION b(7,nv),g(7,nv)
c times
      DOUBLE PRECISION t,t2,tm,tini      
c =====================END INTERFACE===========================
c memory control
      INTEGER ips,imem
c constants of the integration method
      DOUBLE PRECISION h(8),w(7),u(7),c(21),d(21),r(21),w1
      COMMON/raconst/h,w,u,c,d,r,w1
c state variables 
      INTEGER ndx,nvx
      PARAMETER (ndx=60,nvx=2*ndx)
c x,v temporary
      DOUBLE PRECISION y(nvx),z(nvx)
c error control: convergence
      INCLUDE 'parint.h'
      INCLUDE 'comint.h'
      DOUBLE PRECISION ck(nvx,8),fj(nvx)
c total number of function calls
      INTEGER nf
c loop indexes
      INTEGER j,k,jd
c scalar temporaries
      DOUBLE PRECISION q,s,temp,gk
c close approach flags: detected, which planet
      INTEGER idc
c positon and velocities of planet which has close-encounter
c whith the asteroid
      DOUBLE PRECISION xpla(6)
c memory model static
      SAVE
c initialize previous iteration array to zero if first iteration
      IF(m.eq.1)THEN
         DO  k=1,nv
            DO j=2,8
               ck(k,j)=0.d0
            ENDDO
         ENDDO
      ENDIF
      epsi=0.d0 
c =========do on substep========================= 
c  loop 174 is for each substep within a sequence.
      do 174 j=2,8
          jd=j-1
c         jdm=j-2
          s=h(j)
          q=s
          IF(ncl) q=1.d0
c ===============================================================
c  use eqs. (2.9) and (2.10) of text to predict positions at each substep.
c  these collapsed series are broken into two parts because an otherwise
c  excellent  compiler could not handle the complicated expression.
c ===============================================================
          do 130 k=1,nv
            temp=w(3)*b(3,k)+s*(w(4)*b(4,k)+s*(w(5)*b(5,k)
     +        +s*(w(6)*b(6,k)+s*w(7)*b(7,k))))
            y(k)=x(k)+q*(t*v(k)+t2*s*(f1(k)*w1+s*(w(1)*b(1,k)
     +          +s*(w(2)*b(2,k)+s*temp))))
            IF(.not.npq)then
c  next are calculated the velocity predictors need for general class ii.
               temp=u(3)*b(3,k)+s*(u(4)*b(4,k)+s*(u(5)*b(5,k)
     +              +s*(u(6)*b(6,k)+s*u(7)*b(7,k))))
               z(k)=v(k)+s*t*(f1(k)+s*(u(1)*b(1,k)+s*(u(2)*b(2,k)
     +             +s*temp)))
            ENDIF
 130      continue
c ==================================================================
c  find forces at each substep.
c ==================================================================
          IF(m.eq.1)THEN
             ips=0
          ELSE
             ips=-1
          ENDIF
          imem=j
          CALL force(y,z,tm+s*t+tini,fj,nv,idc,xpla,ips,imem)
          do  k=1,nv
            epsi=epsi+dabs(fj(k)-ck(k,j))
          enddo
          do  k=1,nv
            ck(k,j)=fj(k)
          enddo
          nf=nf+1
c ================do on components===================
          do 171 k=1,nv
c ===============================================================
c  find g-value for the force fj found at the current substep. this
c  section, including the many-branched goto, uses eq. (2.4) of text.
c ===============================================================
            temp=g(jd,k)
            gk=(fj(k)-f1(k))/s
            IF(j.le.2)then
               g(1,k)=gk
            ELSEIF(j.eq.3)then
                g(2,k)=(gk-g(1,k))*r(1)
            ELSEIF(j.eq.4)then
                g(3,k)=((gk-g(1,k))*r(2)-g(2,k))*r(3)
            ELSEIF(j.eq.5)then
                g(4,k)=(((gk-g(1,k))*r(4)-g(2,k))*r(5)-g(3,k))*r(6)
            ELSEIF(j.eq.6)then
                g(5,k)=((((gk-g(1,k))*r(7)-g(2,k))*r(8)-g(3,k))*r(9)-
     +                    g(4,k))*r(10)
            ELSEIF(j.eq.7)then
                g(6,k)=(((((gk-g(1,k))*r(11)-g(2,k))*r(12)
     +                -g(3,k))*r(13)-g(4,k))*r(14)-g(5,k))*r(15)
            ELSEIF(j.eq.8)then
                g(7,k)=((((((gk-g(1,k))*r(16)-g(2,k))*r(17)
     +                -g(3,k))*r(18)-g(4,k))*r(19)
     +                -g(5,k))*r(20)-g(6,k))*r(21)
            ENDIF
c ===============================================================
c  upgrade all b-values
c ===============================================================
            temp=g(jd,k)-temp
            b(jd,k)=b(jd,k)+temp
c ===============================================================
c  temp is now the improvement on g(jd,k) over its former value.
c  now we upgrade the b-value using this dfference in the one term.
c  this section is based on eq. (2.5).
c ===============================================================
            IF(j.eq.3)then
                b(1,k)=b(1,k)+c(1)*temp
            ELSEIF(j.eq.4)then
                b(1,k)=b(1,k)+c(2)*temp
                b(2,k)=b(2,k)+c(3)*temp
            ELSEIF(j.eq.5)then
                b(1,k)=b(1,k)+c(4)*temp
                b(2,k)=b(2,k)+c(5)*temp
                b(3,k)=b(3,k)+c(6)*temp
            ELSEIF(j.eq.6)then
                b(1,k)=b(1,k)+c(7)*temp
                b(2,k)=b(2,k)+c(8)*temp
                b(3,k)=b(3,k)+c(9)*temp
                b(4,k)=b(4,k)+c(10)*temp
            ELSEIF(j.eq.7)then
                b(1,k)=b(1,k)+c(11)*temp
                b(2,k)=b(2,k)+c(12)*temp
                b(3,k)=b(3,k)+c(13)*temp
                b(4,k)=b(4,k)+c(14)*temp
                b(5,k)=b(5,k)+c(15)*temp
            ELSEIF(j.eq.8)then
                b(1,k)=b(1,k)+c(16)*temp
                b(2,k)=b(2,k)+c(17)*temp
                b(3,k)=b(3,k)+c(18)*temp
                b(4,k)=b(4,k)+c(19)*temp
                b(5,k)=b(5,k)+c(20)*temp
                b(6,k)=b(6,k)+c(21)*temp
            ENDIF
c =========end do on components========================= 
 171      continue
c =========end do on substeps========================= 
 174    continue
        RETURN
        END

c =====================================================
c rapred
c
c finds new x and v values at end of sequence using eqs. (2.11),(2.12)
c =====================================================
      SUBROUTINE rapred(ncl,x,v,t,t2,f1,b,nv)
      IMPLICIT NONE
c flag for first order equations
      LOGICAL  ncl
c actual dimension
      INTEGER nv 
c stepsize, squared for second order eq.
      DOUBLE PRECISION t,t2
c state variables 
      DOUBLE PRECISION x(nv),v(nv)
c right hand side at beginning of step
      DOUBLE PRECISION f1(nv)
c storage arrays
      DOUBLE PRECISION b(7,nv)
c end interface
c loop indexes
      INTEGER k
c ========================================================
c constants of the integration method
      DOUBLE PRECISION h(8),w(7),u(7),c(21),d(21),r(21),w1
      COMMON/raconst/h,w,u,c,d,r,w1
c ============================================
      DO k=1,nv
        x(k)=x(k)+v(k)*t+t2*(f1(k)*w1+b(1,k)*w(1)+b(2,k)*w(2)
     +   +b(3,k)*w(3)+b(4,k)*w(4)+b(5,k)*w(5)+b(6,k)*w(6)+b(7,k)*w(7))
        IF(.not.ncl)THEN
           v(k)=v(k)+t*(f1(k)+b(1,k)*u(1)+b(2,k)*u(2)+b(3,k)*u(3)
     +       +b(4,k)*u(4)+b(5,k)*u(5)+b(6,k)*u(6)+b(7,k)*u(7))
        ENDIF
      ENDDO
      RETURN
      END
c ===============================================================
c rabeta
c
c  find new beta-values from the predicted b-values, 
c  following eq. (2.7) in text.
c ===============================================================
      SUBROUTINE rabeta(nv,b,d,g)
      IMPLICIT NONE
c number of equations
      INTEGER nv
c storage arrays
      DOUBLE PRECISION b(7,nv),g(7,nv),d(21)
c loop indexes
      INTEGER k
c =========================
      DO k=1,nv
        g(1,k)=b(1,k)+d(1)*b(2,k)+d(2)*b(3,k)+
     +             d(4)*b(4,k)+d( 7)*b(5,k)+d(11)*b(6,k)+d(16)*b(7,k)
        g(2,k)=            b(2,k)+d(3)*b(3,k)+
     +             d(5)*b(4,k)+d( 8)*b(5,k)+d(12)*b(6,k)+d(17)*b(7,k)
        g(3,k)=            b(3,k)+
     +             d(6)*b(4,k)+d( 9)*b(5,k)+d(13)*b(6,k)+d(18)*b(7,k)
        g(4,k)=         b(4,k)+d(10)*b(5,k)+d(14)*b(6,k)+d(19)*b(7,k)
        g(5,k)=                      b(5,k)+d(15)*b(6,k)+d(20)*b(7,k)
        g(6,k)=                                   b(6,k)+d(21)*b(7,k)
        g(7,k)=                                                b(7,k)
      ENDDO
      RETURN
      END  
c ===============================================================
c  bintrp
c
c  now predict b-values for next step. the predicted values from the preceding
c  sequence were saved in the e-matrix. te correction bd between the actual
c  b-values found and these predicted values is applied in advance to the
c  next sequence. the gain in accuracy is significant. using eqs. (2.13):
c ===============================================================
      SUBROUTINE bintrp(q,b,e,bd,nv,ns)
      IMPLICIT NONE
c =============INPUT===============
c rescaled time
      DOUBLE PRECISION q
c dimension of state vector, number of steps completed
      INTEGER nv,ns
c =============INPUT AND OUTPUT=============
c storage arrays
      DOUBLE PRECISION e(7,nv),b(7,nv),bd(7,nv)
c =============END INTERFACE================
      INTEGER k,j,l
c ===============================================================
      DO 39 k=1,nv
        IF(ns.ne.1)then
           DO  j=1,7
             bd(j,k)=b(j,k)-e(j,k)
           ENDDO

        ENDIF
        e(1,k)=      q*(b(1,k)+ 2.d0*b(2,k)+ 3.d0*b(3,k)+
     x           4.d0*b(4,k)+ 5.d0*b(5,k)+ 6.d0*b(6,k)+ 7.d0*b(7,k))
        e(2,k)=                q**2*(b(2,k)+ 3.d0*b(3,k)+
     y           6.d0*b(4,k)+10.d0*b(5,k)+15.d0*b(6,k)+21.d0*b(7,k))
        e(3,k)=                             q**3*(b(3,k)+
     z           4.d0*b(4,k)+10.d0*b(5,k)+20.d0*b(6,k)+35.d0*b(7,k))
        e(4,k)= q**4*(b(4,k)+ 5.d0*b(5,k)+15.d0*b(6,k)+35.d0*b(7,k))
        e(5,k)=              q**5*(b(5,k)+ 6.d0*b(6,k)+21.d0*b(7,k))
        e(6,k)=                           q**6*(b(6,k)+ 7.d0*b(7,k))
        e(7,k)=                                         q**7*b(7,k)
        DO  l=1,7
          b(l,k)=e(l,k)+bd(l,k)
        ENDDO
 39   ENDDO
      RETURN
      END
c ===============================================================
c  RADCON assigns h spacings and
c  evaluates the constants in the w-, u-, c-, d-, and r-vectors
c ===========begin interface=========================================
      SUBROUTINE radcon(ncl)
      IMPLICIT NONE
      LOGICAL ncl
      DOUBLE PRECISION h(8),w(7),u(7),c(21),d(21),r(21),w1
      COMMON/raconst/h,w,u,c,d,r,w1
c ===========end interface===========================================
c loop indexes
      INTEGER n,k,la,lc,lb,ld,le,l
c scalar real
      DOUBLE PRECISION ww
      DOUBLE PRECISION half,one
c integers
      INTEGER nw(8)
      DATA half, one/0.5d0, 1.0d0/
      DATA nw/0,0,1,3,6,10,15,21/
c ===================================================================
c  these h values are the gauss-radau spacings, scaled to the range 0 to 1,
c  for integrating to order 15.
      h(1)=0.d0
      h(2)= .05626256053692215d0
      h(3)= .18024069173689236d0
      h(4)=.35262471711316964d0 
      h(5)=.54715362633055538d0
      h(6)= .73421017721541053d0
      h(7)=.88532094683909577d0
      h(8)= .97752061356128750d0
c  the sum of the h-values should be 3.73333333333333333
c ====================================================================
c  evaluate the constants in the w-, u-, c-, d-, and r-vectors
      do 14 n=2,8
        ww=n+n*n
        IF(ncl) ww=n
        w(n-1)=one/ww
        ww=n
  14    u(n-1)=one/ww
      w1=half
      IF(ncl) w1=one
      c(1)=-h(2)
      d(1)=h(2)
      r(1)=one/(h(3)-h(2))
      la=1
      lc=1
      do 73 k=3,7
        lb=la
        la=lc+1
        lc=nw(k+1)
        c(la)=-h(k)*c(lb)
        c(lc)=c(la-1)-h(k)
        d(la)=h(2)*d(lb)
        d(lc)=-c(lc)
        r(la)=one/(h(k+1)-h(2))
        r(lc)=one/(h(k+1)-h(k))
        IF(k.eq.3) go to 73
        do 72 l=4,k
          ld=la+l-3
          le=lb+l-4
          c(ld)=c(le)-h(k)*c(le+1)
          d(ld)=d(le)+h(l-1)*d(le+1)
  72      r(ld)=one/(h(k+1)-h(l-1))
  73  continue
      return 
      end
c ================================
c invaxv
c
c initialise the  variation matrix as the 6 X 6 identity
c ================================
      SUBROUTINE invaxv(x,v,nvar2)
      IMPLICIT NONE
      INTEGER nvar2
      DOUBLE PRECISION x(nvar2),v(nvar2)
c end interfface
      INTEGER i,j,iii,ij        
      iii=3
      do 7 j=1,6
         do  i=1,3
            ij=i+3*(j-1)
            x(iii+ij)=0.d0
            v(iii+ij)=0.d0
            if(i.eq.j)then
               x(iii+ij)=1.d0
            elseif(j.eq.i+3)then
               v(iii+ij)=1.d0
            endif
         enddo
 7    continue
      RETURN
      END   
c ====================================================
c vawrxv
c
c First parzial derivatives: rewrap vector into 6x6 matrix
c ====================================================
      SUBROUTINE vawrxv(x,v,dxdx0,nvar2)
      IMPLICIT NONE
      INTEGER nvar2
      DOUBLE PRECISION dxdx0(6,6),x(nvar2),v(nvar2)
      INTEGER i,j,ij
c ====================================================
      DO j=1,3
         DO  i=1,3
            ij=i+3*(j-1)+3
            dxdx0(i,j)=x(ij)
            dxdx0(i,j+3)=x(ij+9)
            dxdx0(i+3,j)=v(ij)
            dxdx0(i+3,j+3)=v(ij+9)
         ENDDO 
      ENDDO
      RETURN
      END
c ====================================================
c varunw
c
c First parzial derivatives: unwrap 6x6 matrix into vector
c ====================================================
      SUBROUTINE varunw(dxdx0,x,v,nvar2)
      IMPLICIT NONE
      INTEGER nvar2
      DOUBLE PRECISION dxdx0(6,6),x(nvar2),v(nvar2)
      INTEGER i,j,ij
c ====================================================
      DO j=1,3
         DO  i=1,3
            ij=i+3*(j-1)+3
            x(ij)=dxdx0(i,j)
            x(ij+9)=dxdx0(i,j+3)
            v(ij)=dxdx0(i+3,j)
            v(ij+9)=dxdx0(i+3,j+3)
         ENDDO 
      ENDDO
      RETURN
      END
c last update  Jun 14 1997 by A. Milani
c =====================================================================
c PROPA2 
c =====================================================================
c 2-body problem propagator
c
c This version uses JPL ephemerides as source for the Earth position
c
c  WARNING: the input elements east and the output cartesian
c           coordinates and derivatives are in ecliptic (mean of J2000.0)
c           coordinates.
c =====================================================================
c
c  input: t0 epoch time
c         t2 prediction time
c         east equinoctal orbital elements vector a,h,k,p,q,lambda 
c               at time t0
c         ider 0=no partials 1= partials dx/deast (3x6 matrix)
c              2=second derivatives (3x6x6 tensor) all in 2-body approximation
c  output:
c        xast position and velocity vector in heliocentric cartesian 
c               coordinates for the asteroid at time t1
c        xea  position and velocity vector in heliocentric cartesian 
c               coordinates for the Earth at time t1
c        dxde first derivatives of position vector with respect to elements 
c
c        ddxde second derivatives of position vector with respect to elem.
c
c ======================INTERFACE=======================================
      SUBROUTINE propa2(t0,east,t2,xast,xea,ider,dxde,ddxde)
      implicit none
c =========INPUT========================================================
c times
      double precision t0,t2,east(6)
      integer ider
c =========OUTPUT=======================================================
      double precision xast(6),xea(6),dxde(6,6),ddxde(3,6,6)
c =====================================================================
c options common
      include 'model.h'
c ======== constant of gravitation ==============
      include 'sunmass.h'
c workspace for rotations, strings to define ref. systems
      double precision rrd(6),et(2),rot(3,3)
c integers for call to JPl routines
      integer ntarg,ncent,istate
c initialisation: needed only for rotation matrix rot
* **************************************
c static memory allocation only for:
      save rot,lflag
* **************************************
      integer lflag
      data lflag /0/
c =====================================================================
c JPL Earth vector at observation time
      et(1)=2400000.5d0
      et(2)=t2
      ntarg=3
      ncent=11
* ****** added on Sat Jun 14 1997 ******
* first istate need to be=2  (also vel. required for aberration)
      istate=2
* **************************************
      call dpleph(et,ntarg,ncent,rrd,istate)
* Change of reference system EQUM00 ---> ECLM00
      if(lflag.eq.0)then
         call rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0)
         lflag=1
      endif
      call prodmv(xea,rot,rrd)
      call prodmv(xea(4),rot,rrd(4))
*******************
c two body propagation
      call prop2b(t0,east,t2,xast,gms,ider,dxde,ddxde)
c
      return
      end
c =====================================================================




