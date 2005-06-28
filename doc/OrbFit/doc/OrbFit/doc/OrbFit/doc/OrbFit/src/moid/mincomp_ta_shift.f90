! ******************************************************************
! **************  PROGRAM   M I N _ C O M P _ T A ******************
! ******** test for the subroutine compute_minima_ta.f90 ***********
! ************* written by GIOVANNI F. GRONCHI (2001) **************
! ********** Department of Mathematics, UNIVERSITY of PISA *********
! ******************************************************************
! last modified 24/06/2005 GFG
! ==================================================================
  PROGRAM min_comp_ta 
    USE fund_const
    USE orbit_elements
    USE output_control
    USE critical_points_shift, ONLY : poldeg 
    IMPLICIT NONE 
! ==================================================================
! asteroid counter                                                  
    INTEGER :: astnum,loop 
    INTEGER, PARAMETER :: nastmax = 2000
    INTEGER, PARAMETER :: num_obj = 1
    CHARACTER*9 :: name, namev(1) ! asteroid name 
! orbital elements                                    
    TYPE(orbit_elem) :: elcarpl ! cartesian elements of the planet
    TYPE(orbit_elem) :: elcar   ! cartesian elements of the asteroid/comet
! for compute_minima_ta
    REAL(KIND=8) :: carpl(6) ! cartesian coords of the planet
    REAL(KIND=8) :: car(6)   ! cartesian coords of the asteroid/comet
!    INTEGER,PARAMETER :: poldeg=16
    DOUBLE PRECISION :: cmin(3,poldeg),cplmin(3,poldeg) 
    DOUBLE PRECISION :: D2(poldeg) ! squared distance function 
                                                  ! at minimum points 
    INTEGER :: nummin ! number of relative minima found 
! auxiliary
    TYPE(orbit_elem) :: eqcom   ! equinoctal elements of the comet
    TYPE(orbit_elem) :: eqp     ! equinoctal elements of the Earth
    DOUBLE PRECISION,DIMENSION(6) :: eqtmp,eqptmp ! auxiliary
    DOUBLE PRECISION :: elem1(6)
    DOUBLE PRECISION :: t1,t1v(1)! current time (JD) 
    INTEGER :: fail_flag ! for coo_cha
! for subroutine rdelem                                             
    DOUBLE PRECISION :: cov0(6,6),norm0(6,6),mass,massv(1) 
    DOUBLE PRECISION :: hmag,gmag,enne 
    DOUBLE PRECISION :: hmagv(1),gmagv(1)
    CHARACTER*3 eltype,eltypev(1) 
    CHARACTER*80 comele,comelev(1) 
! successful input flags                                            
    LOGICAL defcov,defelem,neodys_elem 
    LOGICAL defcovv(1),defelemv(1)
! available data                                                    
    LOGICAL ok 
! file names                                                        
    CHARACTER*60 file,elefil,eledir,elemfiles(4) 
    INTEGER le,lnam,num_files 
    CHARACTER*6 progna 
    character*80 :: run 
! file identifiers                                                  
    INTEGER :: iunlog,iunit,iunhand 
! loop index                                                        
    INTEGER :: i,j 
! Earth orbital elements
!    DOUBLE PRECISION :: eqp(6),ekpl(6) 
    DOUBLE PRECISION :: apl,epl,Ipl,omegapl,Ompl,lpl 
!                                                                       
    CHARACTER*60 filcat 
    INTEGER :: iplam 
! =================================================================

! input options                                                     
!    progna='statco' 
    progna='fitobs'
    run='moid_comets'
! read options
    CALL trivopt(progna,run,iunlog)
    CALL rmodel
    CALL trange
! verbosity control
!    verb_moid = 20
! =============================                                     
! check for non-existing options                                    
    call rmsp(progna,le) 
    file=progna(1:le)//'.key' 
    CALL rdklst(file) 
    CALL chkkey 
! ==============================                                    
! read option for physical model and integration method             
!    CALL rmodel
    eledir='./' 
! CALL optpro(progna,eledir)                                        
! =======================================================           
! check files                                                       
    CALL filopn(iunlog,'mincomp.log','unknown') 
    CALL filopn(iunhand,'mincomp.names.fail','unknown') 
    astnum = 0 
! =========== ENTER MAIN LOOP ===========================           
    DO 1 loop=1,nastmax 
       astnum = astnum + 1 
! read next name                                                    
       write(*,*) 'Name (Ctrl-D to quit)?' 
5      READ(*,100,end=999) name 
!     5       READ(iunlst,100,end=999) name                             
100    FORMAT(a9) 
       CALL rmsp(name,lnam) 
       IF(name(1:1).eq.'!')THEN 
          WRITE(iunhand,*) name,'Commented near line ',astnum 
          GOTO 5 
       ENDIF
       namev(1) = name
       defelem=.false. 
       defelemv(1)=.false.
       defcov=.false. 
       defcovv(1)=.false. 
! put flag in each log file                                         
       WRITE(iunlog,*)'***************************************' 
       WRITE(iunlog,*)' Processing ',name,' number ',astnum 
       WRITE(iunlog,*)'***************************************' 
! input file                                                        
       CALL filnam(eledir,name,'eq1',elefil,le) 
       INQUIRE(file=elefil(1:le),exist=neodys_elem) 
       IF(neodys_elem)THEN 
          elemfiles(1) = elefil(1:le) 
          num_files = 1 
       ELSE 
          WRITE(*,*)' elements not found for ',name 
          GOTO 2 
       ENDIF
       CALL rdelem(iunlog,namev,num_obj,elemfiles,num_files,defelemv, &
            & defcovv,eltypev,t1v,eqtmp,cov0,norm0,massv,hmagv,gmagv, &
            & comelev) 
       IF(.not.defelemv(1))THEN 
          WRITE(iunlog,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
          WRITE(iunlog,*)'!!WARNING! WARNING! WARNING! WARNING!!' 
          WRITE(iunlog,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
          WRITE(iunlog,*)'Asteroid ',name,' number ',astnum,   &
               & ' not found in asteroid catalog.'             
          WRITE(iunhand,*)name,'not in catalogs.' 
          GOTO 2 
       ENDIF
       eltype = eltypev(1)
       t1=t1v(1)
       mass=massv(1)

! comet elements
       eqcom=undefined_orbit_elem
       eqcom%coord(1:6)=eqtmp(1:6)
       eqcom%coo='EQU'
       eqcom%t=t1
!       write(*,*)'asteroid/comet elems (equinoctal):',eqcom%coord(1:6)
! conversion into cartesian elems
       CALL coo_cha(eqcom,'CAR',elcar,fail_flag) 
       elcar%coo='COM'
       elcar%t=elcar%coord(6)
!       write(*,*)'asteroid/comet elems (cartesian):',elcar%coord(1:6)
       
! Earth elements
       CALL earth(t1,eqptmp(1:6))
       eqp=undefined_orbit_elem
       eqp%coord(1:6)=eqptmp(1:6)
       eqp%coo='EQU'
       eqp%t=t1
       CALL coo_cha(eqp,'CAR',elcarpl,fail_flag) !conversion into cartesian elems
!       write(*,*)'time (MJD):',elcarpl%t
!       write(*,*)'Earth elems (cartesian):',elcarpl%coord(1:6)

       car(1:6) = elcar%coord(1:6)
       carpl(1:6) = elcarpl%coord(1:6)

!     select Earth as planet                                            
       iplam = 3 
                                                                        
!     =======================================================           
       CALL compute_minima_ta(car,carpl,iplam,cmin,cplmin,D2,nummin) 
!     =======================================================           
!       write(*,*)'number of minimum points for ',name(1:lnam), ' =',nummin 
                                                                        
! output file                                                       
       CALL filnam('.',name,'minpts',elefil,le) 
       CALL filopn(iunit,elefil(1:le),'unknown') 

       WRITE(*,*)'##################################################' 
       WRITE(*,*)'######### M I N I M U M    P O I N T S ###########' 
       WRITE(*,*)'##################################################' 
       WRITE(*,*)'=======================================================&
            &======================='                                          
!     loop on number of stationary points                               
       DO j = 1,nummin 
          write(*,*) 
!     writing on screen                                                 
          WRITE(*,*)'minimum point label =',j,'  local MOID =',dsqrt(D2(j&
               &        ))
          WRITE(*,*)' '
          WRITE(*,*)'Asteroid cartesian coordinates at this minimum point'   
          WRITE(*,*)'      x             y              z'
          WRITE(*,107)cmin(1,j),cmin(2,j),cmin(3,j)
          WRITE(*,*)'Earth cartesian coordinates at this minimum point ' 
          WRITE(*,*)'      xpl           ypl            zpl'
          WRITE(*,107)cplmin(1,j),cplmin(2,j),cplmin(3,j)
          WRITE(*,*)'====================================================&
               &==========================' 
!     writing on file iunit.minpts                                      
          WRITE(iunit,106)j,dsqrt(D2(j)) 
          WRITE(iunit,107)cmin(1,j),cmin(2,j),cmin(3,j)
          WRITE(iunit,107)cplmin(1,j),cplmin(2,j),cplmin(3,j)
106       FORMAT(2x,i2,15x,f13.8) 
107       FORMAT(2x,f11.6,2x,f11.6,2x,f11.6,2x,f11.6,2x,f11.6,2x,f11.6) 
                                                                        
!     check on the distances                                            
!         write(*,*)'-----------------------------------------'         
!         WRITE(*,*)'controllo di D2(',j,')=',dsqrt((cmin(1,j)- &
!              & cplmin(1,j))**2 + (cmin(2,j)-cplmin(2,j))**2 + &
!              & (cmin(3,j)-cplmin(3,j))**2)
!         write(*,*)'-----------------------------------------'         
                                                                        
!     END MAIN DO LOOP                                                  
       ENDDO

! close output file                                                 
       CALL filclo(iunit,' ') 
! clean up err and clo files                                        
2      CONTINUE 
! end of main loop                                                  
1   ENDDO
    STOP ' Warning, too many objects. Reached end of loop.' 
! EOF exit point                                                    
999 CONTINUE 
    CALL filclo(iunlog,' ') 
    CALL filclo(iunhand,' ') 
    
  END PROGRAM min_comp_ta
