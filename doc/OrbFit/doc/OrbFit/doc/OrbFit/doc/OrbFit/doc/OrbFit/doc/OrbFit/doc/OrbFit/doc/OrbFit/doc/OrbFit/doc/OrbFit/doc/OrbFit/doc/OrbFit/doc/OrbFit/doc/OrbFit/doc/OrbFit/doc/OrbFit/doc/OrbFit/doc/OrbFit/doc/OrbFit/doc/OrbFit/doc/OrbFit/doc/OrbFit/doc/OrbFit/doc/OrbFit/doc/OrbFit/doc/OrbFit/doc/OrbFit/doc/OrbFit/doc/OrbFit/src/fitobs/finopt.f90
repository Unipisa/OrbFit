! ===================================================================   
! FINOPT                                                                
! ===================================================================   
! input options, for the propagator and the specific main program       
! input: progna = program name (6 characters)                           
!        run    = run identifier (80 characters, up to 76 non-blank)    
SUBROUTINE finopt(progna,run,astna0,astnap,error_model) 
  USE output_control
  implicit none 
  character*6,INTENT(IN) :: progna 
  character*80,INTENT(IN) :: run 
  CHARACTER*(*),INTENT(OUT) :: astna0,astnap 
  CHARACTER*(20),INTENT(OUT) :: error_model ! error model file name
! hidden output through bizset: controls for bizarre orbits
  DOUBLE PRECISION ecclim,samin,samax,phmin,ahmax,qmax 
! ==========END INTERFACE============================================   
  integer le
  character*90 file 
  LOGICAL ireq,found
  CHARACTER*60 comment 
  CHARACTER*100 filnam
  INTEGER iunout, iuncovar ! now only local names
! =============================
  CALL initopt(progna,run,'fop')                                         
! read option for physical model and integration method                 
  CALL rmodel 
! ==============================                                        
! initialisations for Gauss method                                      
  CALL iodini 
! ===================================================================== 
! Output files: for control and results, for covariance  
  filnam=run//'.fou' 
  CALL rmsp(filnam,le) 
  call filopn(iunout,filnam,'UNKNOWN')
  iun_log=iunout ! later iunout to be abolished
  filnam=run//'.fga' 
  CALL rmsp(filnam,le) 
  call filopn(iuncovar,filnam,'UNKNOWN')
  iun_covar=iuncovar
! ===================================================================== 
! Output files: for errors, close approaches, for propagator parameters 
  filnam=run//'.err' 
  CALL rmsp(filnam,le)
  call filopn(ierrou,filnam,'UNKNOWN') 
  numerr=0 
  filnam=run//'.clo' 
  CALL rmsp(filnam,le)
  call filopn(iuncla,filnam,'UNKNOWN') 
  numcla=0
  filnam=run//'.pro' 
  CALL rmsp(filnam,le)
  call filopn(ipirip,filnam,'UNKNOWN')
! =============================                                         
! asteroid name(s)                                                      
! =============SOME ASTEROID NAME NEEDED===================
  ireq=.false.
  astna0=run
  comment='first arc asteroid name'
  CALL input_cha_opt(progna,'astna0',astna0,ireq,found,comment,iunout)
  CALL rmsp(astna0,le)
  ireq=.false. 
  astnap=' '
  comment='second arc asteroid name'
  CALL input_cha_opt(progna,'astnap',astnap,ireq,found,comment,iunout)
  CALL rmsp(astnap,le)
! =============SELECTED ERROR MODEL===================   
  error_model=' '
  ireq=.false.
  comment='error model file'
  CALL input_cha_opt(progna,'error_model',error_model,ireq,found,comment,iunout)
! =============bizarre control=======================   
! ================limit eccentricity============================            
  ireq=.false. 
  ecclim=0.d0  ! leave at the bizset default
  comment='limit eccentricity'
  CALL input_rea_opt(progna,'ecclim',ecclim,ireq,found,comment,iunout)
! ================limit semimajor axis============================            
  ireq=.false. 
  samax=0.d0 !leave at the bizset default
  comment='max semimajor axis'
  CALL input_rea_opt(progna,'samax',samax,ireq,found,comment,iunout)
  ireq=.false. 
  samin=0.d0 !leave at the bizset default
  comment='min semimajor axis'
  CALL input_rea_opt(progna,'samin',samin,ireq,found,comment,iunout)
! ================limit perihelion, aphelion============================
  ireq=.false. 
  ahmax=0.d0 !leave at the bizset default
  comment='max aphelion'
  CALL input_rea_opt(progna,'ahmax',ahmax,ireq,found,comment,iunout)
  ireq=.false. 
  phmin=0.d0 !leave at the bizset default
  comment='min perihelion'
  CALL input_rea_opt(progna,'phmin',phmin,ireq,found,comment,iunout)
  qmax=0.d0 !leave at the bizset default
  comment='max perihelion'
  CALL input_rea_opt(progna,'qmax',qmax,ireq,found,comment,iunout)
! ================control on bizarre orbits================ 
  CALL bizset(ecclim,samin,samax,phmin,ahmax,qmax) 
END SUBROUTINE finopt
