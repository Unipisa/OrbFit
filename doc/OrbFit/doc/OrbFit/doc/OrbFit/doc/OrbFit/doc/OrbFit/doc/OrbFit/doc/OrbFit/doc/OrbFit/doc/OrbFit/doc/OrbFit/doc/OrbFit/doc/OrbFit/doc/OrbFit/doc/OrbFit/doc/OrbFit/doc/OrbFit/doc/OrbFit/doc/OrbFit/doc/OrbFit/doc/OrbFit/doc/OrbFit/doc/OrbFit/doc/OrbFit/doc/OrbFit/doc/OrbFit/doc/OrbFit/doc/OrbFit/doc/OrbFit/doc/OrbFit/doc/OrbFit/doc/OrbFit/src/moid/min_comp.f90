! ******************************************************************
! *****************  PROGRAM   M I N _ C O M P  ********************
! ******************************************************************
! ************* written by GIOVANNI F. GRONCHI (2001) **************
! ******************************************************************
! *************** E-MAIL gronchi@mail.dm.unipi.it ******************
! ********** Department of Mathematics, UNIVERSITY of PISA *********
! ******************************************************************
! last modified July 2003 (GFG)
! ==================================================================
  PROGRAM min_comp 
    USE fund_const
    USE output_control
    IMPLICIT NONE 
! ==================================================================
! asteroid counter                                                  
    INTEGER :: astnum,loop 
    INTEGER, PARAMETER :: nastmax = 2000
    INTEGER, PARAMETER :: num_obj = 1
    CHARACTER*9 :: name, namev(1) ! asteroid name 
! orbital elements                                                  
    DOUBLE PRECISION :: elkep(6),elem0(6) 
! for subroutine rdelem                                             
    DOUBLE PRECISION :: cov0(6,6),norm0(6,6),mass,massv(1) 
    DOUBLE PRECISION :: hmag,gmag,enne 
    DOUBLE PRECISION :: hmagv(1),gmagv(1)
    CHARACTER*3 eltype,eltypev(1) 
    CHARACTER*80 comele,comelev(1) 
! successful input flags                                            
    LOGICAL defcov,defelem,neodys_elem 
    LOGICAL defcovv(1),defelemv(1)
! current time (JD)                                                 
    DOUBLE PRECISION :: t1,t1v(1) 
! elements                                                          
    DOUBLE PRECISION :: elem1(6) 
! available data                                                    
    LOGICAL ok 
! file names                                                        
    CHARACTER*60 file,elefil,eledir,elemfiles(4) 
    INTEGER le,lnam,num_files 
    CHARACTER*6 progna 
! file identifiers                                                  
    INTEGER :: iunlog,iunit,iunhand,iun1,iun2,iun3,iun4 
! loop index                                                        
    INTEGER :: i 
! Earth orbital elements                                            
    DOUBLE PRECISION :: eqp(6),ekpl(6) 
    DOUBLE PRECISION :: apl,epl,Ipl,omegapl,Ompl,lpl 
!                                                                       
    CHARACTER*60 filcat 
! =================================================================


! input options                                                     
    progna='statco' 
! =============================                                     
! read option for propagator                                        
    CALL libini 
    CALL namini 
! default options are only for propag                               
    CALL filopl(iunit,'propag.def') 
    CALL rdnam(iunit) 
    CALL filclo(iunit,' ') 
! =============================                                     
! check for non-existing options                                    
    call rmsp(progna,le) 
    file=progna(1:le)//'.key' 
    CALL rdklst(file) 
    CALL chkkey 
! ==============================                                    
! read option for physical model and integration method             
    CALL rmodel
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
       WRITE(*,*)'***************************************' 
       WRITE(*,*)'  Processing ',name,' number ',astnum 
       WRITE(*,*)'***************************************' 
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
            & defcovv,eltypev,t1v,elem1,cov0,norm0,massv,hmagv,gmagv, &
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

       CALL coocha(elem1,eltype,gms,elkep,'KEP',enne) 
! output file                                                       
       CALL filnam('.',name,'minpts',elefil,le) 
       CALL filopn(iunit,elefil(1:le),'unknown') 
! elements of the Earth                                             
       CALL earth(t1,eqp) 
       CALL coocha(eqp,'EQU',gms,ekpl,'KEP',enne) 
                                                                        
! ******************************************************            
       CALL writemin(name,lnam,iunit,ekpl,elkep) 
! ******************************************************            
                                                                        
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
    
  END PROGRAM min_comp
