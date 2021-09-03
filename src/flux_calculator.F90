!------------------------------------------------------------------------
! Copyright 2010, CERFACS, Toulouse, France.
! All rights reserved. Use is subject to OASIS3 license terms.
!=============================================================================
!
!
PROGRAM flux_calculator
  !
  ! Use our own modules
  USE flux_calculator_basic     ! type definitions and simple functions e.g. for allocating arrays
  USE flux_calculator_io        ! NetCDF IO functions for flux_calculator
  USE flux_calculator_prepare   ! Functions to check if we have all we need to do the calculations
  USE flux_calculator_calculate ! Functions to actually do the calculations

  USE flux_calculator_parse_arg
  USE flux_calculator_create_namcouple

  ! Use OASIS communication library
  USE mod_oasis

  IMPLICIT NONE
  INCLUDE 'mpif.h'
  !
  ! By default OASIS3 exchanges data in double precision.
  ! To exchange data in single precision with OASIS3, 
  ! the coupler has to be compiled with CPP key "use_realtype_single" 
  ! and the model with CPP key "NO_USE_DOUBLE_PRECISION"

  CHARACTER(len=6),  PARAMETER   :: comp_name = 'flxcal' ! Component name (6 characters) same as in the namcouple    
  CHARACTER(len=50), PARAMETER   :: t_grid_filename='mappings/t_grid_exchangegrid.nc' 
  CHARACTER(len=50), PARAMETER   :: u_grid_filename='mappings/u_grid_exchangegrid.nc'
  CHARACTER(len=50), PARAMETER   :: v_grid_filename='mappings/v_grid_exchangegrid.nc'
  CHARACTER(len=50), PARAMETER   :: regrid_u_to_t_filename='mappings/regrid_u_grid_to_t_grid.nc'
  CHARACTER(len=50), PARAMETER   :: regrid_t_to_u_filename='mappings/regrid_t_grid_to_u_grid.nc'
  CHARACTER(len=50), PARAMETER   :: regrid_v_to_t_filename='mappings/regrid_v_grid_to_t_grid.nc'
  CHARACTER(len=50), PARAMETER   :: regrid_t_to_v_filename='mappings/regrid_t_grid_to_v_grid.nc'

  ! General MPI and output variables
  INTEGER :: mype, npes ! rank and  number of pe
  INTEGER :: localComm  ! local MPI communicator and Initialized
  CHARACTER(len=128) :: comp_out             ! name of the output log file 
  CHARACTER(len=3)   :: chout                ! 3-digit rank number (mype)
  INTEGER :: ierror, rank

  ! OASIS related variables
  INTEGER               :: comp_id             ! component identification
  INTEGER, DIMENSION(3) :: part_id             ! id for my grid partition (t_grid, u_grid, v_grid)
  INTEGER               :: var_nodims(2)       ! number of grid dimensions, number of bundles
  INTEGER               :: var_type            ! data type of variable to define
  INTEGER               :: var_actual_shape(2) ! local dimensions of the arrays to the pe,  2 x field rank (= 2 because fields are of rank = 1)

  ! Namelist variables
  INTEGER                                                    :: timestep            = 0  ! coupling timestep in seconds
  INTEGER                                                    :: num_timesteps       = 0  ! number of time steps in this run
  CHARACTER(50)                                              :: name_atmos_model    = ''
  CHARACTER(50), DIMENSION(MAX_BOTTOM_MODELS)                :: name_bottom_model   = ''
  CHARACTER(1), DIMENSION(MAX_BOTTOM_MODELS)                 :: letter_bottom_model = ''
  INTEGER, DIMENSION(MAX_BOTTOM_MODELS)                      :: num_tasks_per_model = 0
  INTEGER, DIMENSION(MAX_BOTTOM_MODELS, MAX_TASKS_PER_MODEL) :: num_t_grid_cells    = 0  ! if num_tasks_per_model=1 this is automatically read from the grid 
  INTEGER, DIMENSION(MAX_BOTTOM_MODELS, MAX_TASKS_PER_MODEL) :: num_u_grid_cells    = 0
  INTEGER, DIMENSION(MAX_BOTTOM_MODELS, MAX_TASKS_PER_MODEL) :: num_v_grid_cells    = 0
  
  CHARACTER(len=4), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES, MAX_VARS) :: name_bottom_var_t = 'none'  ! define which input and output variables are sent
  CHARACTER(len=4), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES, MAX_VARS) :: name_bottom_var_u = 'none'  ! possible choices:
  CHARACTER(len=4), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES, MAX_VARS) :: name_bottom_var_v = 'none'  ! ALBE AMOI AMOM FARE FICE PATM PSUR QATM TATM TSUR UATM VATM
  CHARACTER(len=4), DIMENSION(MAX_VARS)                                       :: name_atmos_var_t  = 'none'  ! also fluxes are possible to pass them on, 
  CHARACTER(len=4), DIMENSION(MAX_VARS)                                       :: name_atmos_var_u  = 'none'  ! or auxiliary variables typically calculated can be given,
  CHARACTER(len=4), DIMENSION(MAX_VARS)                                       :: name_atmos_var_v  = 'none'  ! these are: QSUR
  CHARACTER(len=4), DIMENSION(MAX_VARS)                                       :: name_send_t       = 'none'  ! possible choices:
  CHARACTER(len=4), DIMENSION(MAX_VARS)                                       :: name_send_u       = 'none'  ! HLAT HSEN MEVA MPRE MRAI MSNO RBBR RLWD RLWU RSWD RSWU UMOM VMOM
  CHARACTER(len=4), DIMENSION(MAX_VARS)                                       :: name_send_v       = 'none'
  LOGICAL,          DIMENSION(MAX_VARS)                                       :: send_to_atmos_t = .TRUE.    ! indicates that this flux will only be sent once to the coupler, not for every surface_type
  LOGICAL,          DIMENSION(MAX_VARS)                                       :: send_to_atmos_u = .TRUE.
  LOGICAL,          DIMENSION(MAX_VARS)                                       :: send_to_atmos_v = .TRUE.
  LOGICAL,          DIMENSION(MAX_BOTTOM_MODELS, MAX_VARS)                    :: send_to_bottom_t = .TRUE.    ! indicates that this flux will only be sent once to the coupler, not for every surface_type
  LOGICAL,          DIMENSION(MAX_BOTTOM_MODELS, MAX_VARS)                    :: send_to_bottom_u = .TRUE.
  LOGICAL,          DIMENSION(MAX_BOTTOM_MODELS, MAX_VARS)                    :: send_to_bottom_v = .TRUE.
  LOGICAL,          DIMENSION(MAX_BOTTOM_MODELS, MAX_VARS)                    :: send_uniform_t = .FALSE.    ! indicates that this flux will only be sent once to the coupler, not for every surface_type
  LOGICAL,          DIMENSION(MAX_BOTTOM_MODELS, MAX_VARS)                    :: send_uniform_u = .FALSE.
  LOGICAL,          DIMENSION(MAX_BOTTOM_MODELS, MAX_VARS)                    :: send_uniform_v = .FALSE.
  
  CHARACTER(len=4), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES, MAX_VARS) :: regrid_u_to_t      = 'none'  ! remap these variables, fluxes etc. from one grid to another
  CHARACTER(len=4), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES, MAX_VARS) :: regrid_v_to_t      = 'none'  
  CHARACTER(len=4), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES, MAX_VARS) :: regrid_t_to_u      = 'none'  
  CHARACTER(len=4), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES, MAX_VARS) :: regrid_t_to_v      = 'none'  

  REAL (kind=wp), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES, MAX_VARS) :: val_bottom_var_t = -1.0e20    ! give a value in the namelist to avoid reading
  REAL (kind=wp), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES, MAX_VARS) :: val_bottom_var_u = -1.0e20    
  REAL (kind=wp), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES, MAX_VARS) :: val_bottom_var_v = -1.0e20  
  REAL (kind=wp), DIMENSION(MAX_VARS)                                       :: val_atmos_var_t  = -1.0e20  
  REAL (kind=wp), DIMENSION(MAX_VARS)                                       :: val_atmos_var_u  = -1.0e20  
  REAL (kind=wp), DIMENSION(MAX_VARS)                                       :: val_atmos_var_v  = -1.0e20
  REAL (kind=wp), DIMENSION(MAX_VARS)                                       :: val_flux_t       = 0.0        ! give a value in the namelist that can be locally 
  REAL (kind=wp), DIMENSION(MAX_VARS)                                       :: val_flux_u       = 0.0        ! overridden by calculation
  REAL (kind=wp), DIMENSION(MAX_VARS)                                       :: val_flux_v       = 0.0

  CHARACTER(len=20), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES) :: which_spec_vapor_surface_t = 'none'    ! 'none', 'copy', 'CCLM'
  CHARACTER(len=20), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES) :: which_spec_vapor_surface_u = 'none'    ! 'none', 'copy', 'CCLM'
  CHARACTER(len=20), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES) :: which_spec_vapor_surface_v = 'none'    ! 'none', 'copy', 'CCLM'

  CHARACTER(len=20), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES) :: which_flux_mass_evap           = 'none'    ! 'none' 'zero', 'copy', 'CCLM'           -> MEVA
  CHARACTER(len=20), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES) :: which_flux_heat_latent         = 'none'    ! 'none' 'zero', 'copy', 'ice', 'water'   -> HLAT
  CHARACTER(len=20), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES) :: which_flux_heat_sensible       = 'none'    ! 'none' 'zero', 'copy', 'CCLM'           -> HSEN
  CHARACTER(len=20), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES) :: which_flux_momentum            = 'none'    ! 'none' 'zero', 'copy', 'CCLM'           -> UMOM, VMOM
  CHARACTER(len=20), DIMENSION(MAX_BOTTOM_MODELS, MAX_SURFACE_TYPES) :: which_flux_radiation_blackbody = 'none'    ! 'none' 'StBo'                           -> RBBR

  NAMELIST /input/ timestep, num_timesteps,                                  &
                   verbosity_level,                                          &
                   name_atmos_model,                                         &
                   name_bottom_model, letter_bottom_model,                   &
                   num_tasks_per_model,                                      &
                   num_t_grid_cells, num_u_grid_cells, num_v_grid_cells,     &
                   name_bottom_var_t, name_bottom_var_u, name_bottom_var_v,  &
                   name_atmos_var_t, name_atmos_var_u, name_atmos_var_v,     & 
                   name_send_t, name_send_u, name_send_v,                    &
                   regrid_t_to_u, regrid_u_to_t,                             &
                   regrid_t_to_v, regrid_v_to_t,                             &
                   send_to_atmos_t, send_to_atmos_u, send_to_atmos_v,        &
                   send_to_bottom_t, send_to_bottom_u, send_to_bottom_v,     &
                   send_uniform_t, send_uniform_u, send_uniform_v,           &
                   val_bottom_var_t, val_bottom_var_u, val_bottom_var_v,     &
                   val_atmos_var_t, val_atmos_var_u, val_atmos_var_v,        & 
                   val_flux_t, val_flux_u, val_flux_v,                       &
                   which_spec_vapor_surface_t, which_spec_vapor_surface_u,   &
                   which_spec_vapor_surface_v,                               &
                   which_flux_mass_evap,                                     &  
                   which_flux_heat_latent, which_flux_heat_sensible,         &
                   which_flux_momentum, which_flux_radiation_blackbody

  ! local variables derived from namelist variables in init procedure
  INTEGER          :: i,j,k,n_timestep
  INTEGER          :: current_time
  CHARACTER(len=4) :: myname
  REAL(kind=wp)    :: myval
  INTEGER                                 :: num_bottom_models   = 1
  INTEGER                                 :: my_bottom_model     = 0
  CHARACTER(len=1)                        :: my_bottom_letter    = ''
  INTEGER, DIMENSION(3)                   :: grid_size           = [0, 0, 0]  ! how many t_grid, u_grid, v_grid exchangegrid cells are calculated locally
  INTEGER, DIMENSION(3)                   :: grid_size_global    = [0, 0, 0]  ! how large is the entire t_grid_exchangegrid
  INTEGER, DIMENSION(3)                   :: grid_offset         = [0, 0, 0]  ! how many grid cells do we skip because lower MPI ranks calculate them
  TYPE(realarray2), DIMENSION(3)          :: grid_longitude_global
  TYPE(realarray2), DIMENSION(3)          :: grid_latitude_global
  TYPE(realarray2), DIMENSION(3)          :: grid_area_global
  INTEGER                                 :: num_surface_types   = 0  ! how many surface types are actually given in my local bottom model
  INTEGER, DIMENSION(MAX_SURFACE_TYPES)   :: num_bottom_vars_t   = 0  ! how many variables from the bottom models are read from each grid (t,u,v)
  INTEGER, DIMENSION(MAX_SURFACE_TYPES)   :: num_bottom_vars_u   = 0  
  INTEGER, DIMENSION(MAX_SURFACE_TYPES)   :: num_bottom_vars_v   = 0
  INTEGER                                 :: num_atmos_vars_t    = 0  ! how many variables from the atmospheric model are read from each grid (t,u,v)
  INTEGER                                 :: num_atmos_vars_u    = 0  
  INTEGER                                 :: num_atmos_vars_v    = 0 
  INTEGER                                 :: num_fluxes_t        = 0  ! how many calculated fluxes will be sent on each grid (t,u,v)
  INTEGER                                 :: num_fluxes_u        = 0 
  INTEGER                                 :: num_fluxes_v        = 0 
  INTEGER                                 :: num_input_fields    = 0
  INTEGER                                 :: num_output_fields   = 0

  TYPE(local_fields_type), DIMENSION(0:MAX_SURFACE_TYPES, 3), TARGET :: local_field   ! first index zero is for entire cell, second index (1..3) is for grid (t,u,v)
  TYPE(io_fields_type),    DIMENSION(:), ALLOCATABLE                 :: input_field
  TYPE(io_fields_type),    DIMENSION(:), ALLOCATABLE                 :: output_field

  TYPE(sparse_regridding_matrix) :: regrid_u_to_t_matrix, regrid_t_to_u_matrix
  TYPE(sparse_regridding_matrix) :: regrid_v_to_t_matrix, regrid_t_to_v_matrix

  LOGICAL :: generate_namcouple = .FALSE.
  
  !###############################################################################
  !# STEP 1:  INITIALIZATION                                                     #
  !###############################################################################

  !###############################################################################
  !# STEP 1.1:  READING NAMELIST                                                 #
  !###############################################################################
  OPEN(10,file='flux_calculator.nml')
  READ(10,nml=input)
  WRITE(*,nml=input)
  CLOSE(10)

  ! Initialize the idx_???? variables which store the index of a variable name
  CALL init_varname_idx

  ! Check if we should generate the namcouple file from here,
  ! if yes, nothing else will be done
  IF (find_argument("--generate_namcouple")) THEN
    generate_namcouple = .TRUE.
  ENDIF


  !###############################################################################
  !# STEP 1.2:  OASIS INITIALIZATION                                             #
  !###############################################################################



  CALL MPI_Init(ierror)

! if we generate the namcouple from here, we must not use OASIS coupler
IF (.NOT. generate_namcouple) THEN
  !!!!!!!!!!!!!!!!! OASIS_INIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL oasis_init_comp (comp_id, comp_name, ierror )   ! get component id
  IF (ierror /= 0) THEN
      WRITE(0,*) 'oasis_init_comp abort by flux_calculator compid ',comp_id
      CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_init_comp')
  ENDIF
ENDIF
  !
  ! Unit for output messages : one file for each process
  CALL MPI_Comm_Rank ( MPI_COMM_WORLD, rank, ierror )   ! get my own rank (globally) - we will not use it
  IF (ierror /= 0) THEN
      WRITE(0,*) 'MPI_Comm_Rank abort by flux_calculator compid ',comp_id
      CALL oasis_abort(comp_id,comp_name,'Failed to call MPI_Comm_Rank')
  ENDIF
  !
  w_unit = 100 + rank
  WRITE(chout,'(I3)') w_unit
  comp_out=comp_name//'.out_'//chout
  !
  OPEN(w_unit,file=TRIM(comp_out),form='formatted')
  WRITE (w_unit,*) '-----------------------------------------------------------'
  WRITE (w_unit,*) TRIM(comp_name), ' Running with reals compiled as kind =',wp
  WRITE (w_unit,*) 'I am component ', TRIM(comp_name), ' rank :',rank
  WRITE (w_unit,*) '----------------------------------------------------------'
  CALL flush(w_unit)

! if we generate the namcouple from here, we must not use OASIS coupler
IF (.NOT. generate_namcouple) THEN
  !
  !!!!!!!!!!!!!!!!! OASIS_GET_LOCALCOMM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  localComm = MPI_COMM_WORLD
  CALL oasis_get_localcomm ( localComm, ierror )   ! get local communicator - we will only use that for finding out "mype"
  IF (ierror /= 0) THEN
      WRITE (w_unit,*) 'oasis_get_localcomm abort by flux_calculator compid ',comp_id
      CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_get_localcomm')
  ENDIF
ELSE
  localComm = MPI_COMM_WORLD
ENDIF
  !
  ! Get MPI size and rank
  CALL MPI_Comm_Size ( localComm, npes, ierror )   ! get number of PEs running flux_calculator
  IF (ierror /= 0) THEN
      WRITE(w_unit,*) 'MPI_comm_size abort by flux_calculator compid ',comp_id
      CALL oasis_abort(comp_id,comp_name,'Failed to call MPI_Comm_Size')
  ENDIF
  !
  CALL MPI_Comm_Rank ( localComm, mype, ierror )   ! get mype (my rank within the flux_calculator instances)
  IF (ierror /= 0) THEN
      WRITE (w_unit,*) 'MPI_Comm_Rank abort by flux_calculator compid ',comp_id
      CALL oasis_abort(comp_id,comp_name,'Failed to call MPI_Comm_Rank')
  ENDIF
  !
  WRITE(w_unit,*) 'I am the ', TRIM(comp_name), ' with local rank ', mype
  WRITE (w_unit,*) 'Number of processors :',npes
 
  !###############################################################################
  !# STEP 1.3:  FIND OUT WHO I AM                                                #
  !###############################################################################

  ! Find out how many bottom models we have
  num_bottom_models=0
  DO i=1,MAX_BOTTOM_MODELS
    IF (TRIM(name_bottom_model(i)) /= '') THEN 
      num_bottom_models=num_bottom_models+1
      IF (num_bottom_models < i) THEN
        WRITE(0,*) 'ERROR: Vector name_bottom_model in flux_calculator.nml must not contain gaps.'
        CALL oasis_abort(comp_id,comp_name,'namelist error')
      ENDIF
    ENDIF
  ENDDO
  
  ! Find out which bottom model this instance of flux_calculator will run
  i=0   ! present pe
  j=1   ! present model
  k=1   ! present task per model
  DO WHILE (i < mype)
    IF (k < num_tasks_per_model(j)) THEN
      k=k+1; j=j+1; i=i+1
    ELSE IF (k == num_tasks_per_model(j)) THEN
      k=1; j=j+1; i=i+1
    ENDIF
    IF (num_tasks_per_model(j)==0) THEN
      j=-1; i=mype
    ENDIF
    IF (j > num_bottom_models) THEN
      j=-1; i=mype
    ENDIF
  ENDDO
  IF (j < 0) THEN
    WRITE (w_unit,*) 'Too many MPI instances of flux_calculator were started. Check your num_tasks_per_model settings in flux_calculator.nml.'
    CALL oasis_abort(comp_id,comp_name,'Failed because of too many MPI instances being called.')
  ELSE
    my_bottom_model = j
  ENDIF
  WRITE(w_unit,*) 'I am calculating fluxes for bottom model ', my_bottom_model
  WRITE(w_unit,*) '                                which is ', trim(name_bottom_model(my_bottom_model))
  ! Find out bottom model letter
  my_bottom_letter = letter_bottom_model(my_bottom_model)
  IF (trim(my_bottom_letter)=='') THEN 
    WRITE (w_unit,*) 'ERROR: letter_bottom_model(',my_bottom_model,') must not be empty.'
    CALL oasis_abort(comp_id,comp_name,'Failed because letter_bottom_model contained empty values.')
  ENDIF
  IF ((my_bottom_letter=='A') .OR. (my_bottom_letter=='R') .OR. (my_bottom_letter=='S')) THEN 
    WRITE (w_unit,*) 'ERROR: letter_bottom_model(',my_bottom_model,') has value ',my_bottom_letter,' which is reserved (A=atmosphere, R=receive, S=send).'
    CALL oasis_abort(comp_id,comp_name,'Failed because letter_bottom_model contained reserved values.')
  ENDIF

  ! Read the exchange grid files and find out how many points this MPI instance has to allocate
  CALL read_scrip_grid_dimensions(t_grid_filename, mype, npes, num_bottom_models, num_tasks_per_model, & ! input
                                  num_t_grid_cells,                                                    & ! output if it is zero, input otherwise
                                  grid_size_global(1), grid_size(1), grid_offset(1),                   & ! output
                                  grid_longitude_global(1), grid_latitude_global(1), grid_area_global(1)) 
  CALL read_scrip_grid_dimensions(u_grid_filename, mype, npes, num_bottom_models, num_tasks_per_model, & ! input
                                  num_u_grid_cells,                                                    & ! output if it is zero, input otherwise
                                  grid_size_global(2), grid_size(2), grid_offset(2),                   & ! output
                                  grid_longitude_global(2), grid_latitude_global(2), grid_area_global(2)) 
  CALL read_scrip_grid_dimensions(v_grid_filename, mype, npes, num_bottom_models, num_tasks_per_model, & ! input
                                  num_v_grid_cells,                                                    & ! output if it is zero, input otherwise
                                  grid_size_global(3), grid_size(3), grid_offset(3),                   & ! output
                                  grid_longitude_global(3), grid_latitude_global(3), grid_area_global(3)) 
  WRITE(w_unit,*) 'My portion of the exchange grid has these sizes:'
  WRITE(w_unit,*) '    t_grid=', grid_size(1)
  WRITE(w_unit,*) '    u_grid=', grid_size(2)
  WRITE(w_unit,*) '    v_grid=', grid_size(3)
  ! read in the regridding matrices
  CALL read_regridding_matrix(regrid_u_to_t_filename,grid_size(2),grid_offset(2),grid_size(1),grid_offset(1),&
                              regrid_u_to_t_matrix)
  CALL read_regridding_matrix(regrid_t_to_u_filename,grid_size(1),grid_offset(1),grid_size(2),grid_offset(2),&
                              regrid_t_to_u_matrix)
  CALL read_regridding_matrix(regrid_v_to_t_filename,grid_size(3),grid_offset(3),grid_size(1),grid_offset(1),&
                              regrid_v_to_t_matrix)
  CALL read_regridding_matrix(regrid_t_to_v_filename,grid_size(1),grid_offset(1),grid_size(3),grid_offset(3),&
                              regrid_t_to_v_matrix)

  !###############################################################################
  !# STEP 1.4:  FIND OUT WHAT I WILL RECEIVE FROM THE COUPLER                    #
  !###############################################################################
  ! The data will be stored in local_field(surface_type,grid_number)%var(var_number)%field(grid_cell_number)
  ! To speed up the input, we will, however, create a list of pointers:
  !   input_field(num_input_field)%field => local_field(surface_type,grid_number)%var(var_number)%field
  ! We will later walk through the array input_field(:) to get the data and they will end up in local_field(:)

  DO i=0,MAX_SURFACE_TYPES
    WRITE(numtype(i),'(I0.2)') i   ! Create a list of strings '00', '01', '02', ...
  ENDDO

  DO i=1,MAX_SURFACE_TYPES   ! Walk through all names of input variables given in the "name_bottom_var" parameter in the namelist
    DO j=1,3
      call nullify_localvars(local_field(i,j))   ! set some POINTER arrays to NULL for initialization
    ENDDO
  ENDDO

  ! Count how many input fields we have in total, to find out how large we need to allocate the input array
  num_input_fields = 0
  num_surface_types = 0
  DO i=1,MAX_SURFACE_TYPES   ! walk through all names of input variables given in the "name_bottom_var" parameter in the namelist
    DO j=1,MAX_VARS
      myname = name_bottom_var_t(my_bottom_model,i,j)
      IF (trim(myname) /= 'none') THEN  
        num_surface_types = max(num_surface_types, i)
        ! check if a constant value was given
        myval = val_bottom_var_t(my_bottom_model,i,j)
        IF (myval < -0.99e20) THEN ! no constant value given
          num_input_fields = num_input_fields + 1
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DO j=1,MAX_VARS
    myname = name_atmos_var_t(j)
    IF (trim(myname) /= 'none') THEN  
      ! check if a constant value was given
      myval = val_atmos_var_t(j)
      IF (myval < -0.99e20) THEN ! no constant value given
        num_input_fields = num_input_fields + 1
      ENDIF
    ENDIF
  ENDDO
  DO i=1,MAX_SURFACE_TYPES   ! the same for u grid
    DO j=1,MAX_VARS
      myname = name_bottom_var_u(my_bottom_model,i,j)
      IF (trim(myname) /= 'none') THEN  
        num_surface_types = max(num_surface_types, i)
        ! check if a constant value was given
        myval = val_bottom_var_u(my_bottom_model,i,j)
        IF (myval < -0.99e20) THEN ! no constant value given
          num_input_fields = num_input_fields + 1
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DO j=1,MAX_VARS
    myname = name_atmos_var_u(j)
    IF (trim(myname) /= 'none') THEN  
      ! check if a constant value was given
      myval = val_atmos_var_u(j)
      IF (myval < -0.99e20) THEN ! no constant value given
        num_input_fields = num_input_fields + 1
      ENDIF
    ENDIF
  ENDDO
  DO i=1,MAX_SURFACE_TYPES   ! the same for v grid
    DO j=1,MAX_VARS
      myname = name_bottom_var_v(my_bottom_model,i,j)
      IF (trim(myname) /= 'none') THEN  
        num_surface_types = max(num_surface_types, i)
        ! check if a constant value was given
        myval = val_bottom_var_v(my_bottom_model,i,j)
        IF (myval < -0.99e20) THEN ! no constant value given
          num_input_fields = num_input_fields + 1
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DO j=1,MAX_VARS
    myname = name_atmos_var_v(j)
    IF (trim(myname) /= 'none') THEN  
      ! check if a constant value was given
      myval = val_atmos_var_v(j)
      IF (myval < -0.99e20) THEN ! no constant value given
        num_input_fields = num_input_fields + 1
      ENDIF
    ENDIF
  ENDDO
  ALLOCATE(input_field(num_input_fields))

  ! Now we will:
  ! (a) allocate the required arrays in local_field(:)
  ! (b) set the pointers in input_field(:)
  ! (c) determine the variable names for OASIS
  num_input_fields = 0
  DO i=1,MAX_SURFACE_TYPES   ! walk through all names of input variables given in the "name_bottom_var" parameter in the namelist
    DO j=1,MAX_VARS
      myname = name_bottom_var_t(my_bottom_model,i,j)
      IF (trim(myname) /= 'none') THEN  
        ! some variable is given here - allocate it
        CALL allocate_localvar(myname, grid_size(1), local_field(i,1))
        ! check if a constant value was given
        myval = val_bottom_var_t(my_bottom_model,i,j)
        IF (myval > -0.99e20) THEN 
          ! a constant value is given in the namelist - set the whole array to this value
          CALL init_localvar(myname, myval, local_field(i,1))
        ELSE IF (myval < -1.99e20) THEN
          ! a value of -2.0e20 means that we will use the value of the first surface_type -> set a pointer to that field
          CALL distribute_input_field(myname, 1, 1, 0, num_surface_types, local_field)   
        ELSE
          ! no constant value is given, so we will need to receive this field from the coupler
          ! we call a subroutine that sets a pointer to the correct local_field in the input_field list
          CALL add_input_field(myname, my_bottom_letter, i, 1, local_field(i,1), num_input_fields, input_field)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DO j=1,MAX_VARS  ! walk through all names of input variables given in the "name_atmos_var" parameter in the namelist
    myname = name_atmos_var_t(j)
    IF (trim(myname) /= 'none') THEN  
      ! some variable is given here - allocate it
      CALL allocate_localvar(myname, grid_size(1), local_field(0,1))   ! index 0 means "all surface types"
      ! check if a constant value was given
      myval = val_atmos_var_t(j)
      IF (myval > -0.99e20) THEN 
        ! a constant value is given in the namelist - set the whole array to this value
        CALL init_localvar(myname, myval, local_field(0,1)) 
      ELSE
        ! no constant value is given, so we will need to receive this field from the coupler
        ! we call a subroutine that sets a pointer to the correct local_field in the input_field list
        CALL add_input_field(myname, 'A', 0, 1, local_field(0,1), num_input_fields, input_field)
      ENDIF
      ! since we need this variable for calculations in each surface_type, set pointers
      CALL distribute_input_field(myname, 1, 0, 0, num_surface_types, local_field)
    ENDIF
  ENDDO
  DO i=1,MAX_SURFACE_TYPES   ! the same for u grid: bottom ...
    DO j=1,MAX_VARS
      myname = name_bottom_var_u(my_bottom_model,i,j)
      IF (trim(myname) /= 'none') THEN  
        ! some variable is given here - allocate it
        CALL allocate_localvar(myname, grid_size(2), local_field(i,2))
        ! check if a constant value was given
        myval = val_bottom_var_u(my_bottom_model,i,j)
        IF (myval > -0.99e20) THEN 
          ! a constant value is given in the namelist - set the whole array to this value
          CALL init_localvar(myname, myval, local_field(i,2))
        ELSE IF (myval < -1.99e20) THEN
          ! a value of -2.0e20 means that we will use the value of the first surface_type -> set a pointer to that field
          CALL distribute_input_field(myname, 2, 1, 0, num_surface_types, local_field)   
        ELSE
          ! no constant value is given, so we will need to receive this field from the coupler
          ! we call a subroutine that sets a pointer to the correct local_field in the input_field list
          CALL add_input_field(myname, my_bottom_letter, i, 2, local_field(i,2), num_input_fields, input_field)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DO j=1,MAX_VARS  ! ... and atmos
    myname = name_atmos_var_u(j)
    IF (trim(myname) /= 'none') THEN  
      ! some variable is given here - allocate it
      CALL allocate_localvar(myname, grid_size(2), local_field(0,2))   ! index 0 means "all surface types"
      ! check if a constant value was given
      myval = val_atmos_var_u(j)
      IF (myval > -0.99e20) THEN 
        ! a constant value is given in the namelist - set the whole array to this value
        CALL init_localvar(myname, myval, local_field(0,2)) 
      ELSE
        ! no constant value is given, so we will need to receive this field from the coupler
        ! we call a subroutine that sets a pointer to the correct local_field in the input_field list
        CALL add_input_field(myname, 'A', 0, 2, local_field(0,2), num_input_fields, input_field)
      ENDIF
      ! since we need this variable for calculations in each surface_type, set pointers
      CALL distribute_input_field(myname, 2, 0, 0, num_surface_types, local_field)
    ENDIF
  ENDDO
  DO i=1,MAX_SURFACE_TYPES   ! the same for v grid: bottom ...
    DO j=1,MAX_VARS
      myname = name_bottom_var_v(my_bottom_model,i,j)
      IF (trim(myname) /= 'none') THEN  
        ! some variable is given here - allocate it
        CALL allocate_localvar(myname, grid_size(3), local_field(i,3))
        ! check if a constant value was given
        myval = val_bottom_var_v(my_bottom_model,i,j)
        IF (myval > -0.99e20) THEN 
          ! a constant value is given in the namelist - set the whole array to this value
          CALL init_localvar(myname, myval, local_field(i,3)) 
        ELSE IF (myval < -1.99e20) THEN
          ! a value of -2.0e20 means that we will use the value of the first surface_type -> set a pointer to that field
          CALL distribute_input_field(myname, 3, 1, 0, num_surface_types, local_field)   
        ELSE
          ! no constant value is given, so we will need to receive this field from the coupler
          ! we call a subroutine that sets a pointer to the correct local_field in the input_field list
          CALL add_input_field(myname, my_bottom_letter, i, 3, local_field(i,3), num_input_fields, input_field)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DO j=1,MAX_VARS  ! ... and atmos
    myname = name_atmos_var_v(j)
    IF (trim(myname) /= 'none') THEN  
      ! some variable is given here - allocate it
      CALL allocate_localvar(myname, grid_size(3), local_field(0,3))   ! index 0 means "all surface types"
      ! check if a constant value was given
      myval = val_atmos_var_v(j)
      IF (myval > -0.99e20) THEN 
        ! a constant value is given in the namelist - set the whole array to this value
        CALL init_localvar(myname, myval, local_field(0,3)) 
      ELSE
        ! no constant value is given, so we will need to receive this field from the coupler
        ! we call a subroutine that sets a pointer to the correct local_field in the input_field list
        CALL add_input_field(myname, 'A', 0, 3, local_field(0,3), num_input_fields, input_field)
      ENDIF
      ! since we need this variable for calculations in each surface_type, set pointers
      CALL distribute_input_field(myname, 3, 0, 0, num_surface_types, local_field)
    ENDIF
  ENDDO

  WRITE(w_unit,*) 'Number of surface types = ',num_surface_types
  WRITE(w_unit,*) '    '
  WRITE(w_unit,*) 'I will read these variables from the coupler:'
  DO i=1,num_input_fields
    WRITE(w_unit,*) '    ',input_field(i)%name,'   on grid ',input_field(i)%which_grid
  ENDDO

  !###############################################################################
  !# STEP 1.5:  FIND OUT WHICH INPUT VARIABLES I SHALL REGRID TO OTHER GRIDS,    #
  !#            ALLOCATE THESE ON THE DESTINATION GRID, AND SET put_to_t_grid /  #
  !#                                                            get_from_t_grid  #
  !###############################################################################
  WRITE(w_unit,*) 'I will regrid these variables after reading from the coupler:'
  ! this output will come from subroutine prepare_regridding
  
  DO i=1,num_input_fields
    CALL prepare_regridding(input_field(i)%idx, input_field(i)%surface_type, local_field, my_bottom_model, &
                            regrid_u_to_t, regrid_v_to_t, regrid_t_to_u, regrid_t_to_v, grid_size)
  ENDDO

  !###############################################################################
  !# STEP 1.6:  FIND OUT WHAT I SHALL CALCULATE,                                 #
  !#            AND IF I WILL HAVE RECEIVED ALL VARIABLES I NEED FOR THAT        #
  !###############################################################################

  ! Call prepare_??? functions. These will do two things:
  ! (a) They make sure that all input variables for the calculations are present and exit with a helpful error message otherwise.
  ! (b) They allocate the arrays for their output variables.

  WRITE(w_unit,*) 'I will regrid these variables after calculation:'
  ! this output will come from subroutine prepare_regridding

  ! AUXILIARY VARIABLES
  !   QSUR:
  DO i=1,num_surface_types
    CALL prepare_spec_vapor_surface(i, 1, which_spec_vapor_surface_t(my_bottom_model,i), grid_size(1), local_field)  ! QSUR on t_grid 
    CALL prepare_spec_vapor_surface(i, 2, which_spec_vapor_surface_u(my_bottom_model,i), grid_size(2), local_field)  ! QSUR on u_grid
    CALL prepare_spec_vapor_surface(i, 3, which_spec_vapor_surface_v(my_bottom_model,i), grid_size(3), local_field)  ! QSUR on v_grid
  ENDDO
  CALL prepare_regridding(idx_QSUR, 0, local_field, my_bottom_model, regrid_u_to_t, regrid_v_to_t, regrid_t_to_u, regrid_t_to_v, grid_size)

  ! MASS FLUXES
  !   MEVA:
  DO i=1,num_surface_types
    CALL prepare_flux_mass_evap(i, 1, which_flux_mass_evap(my_bottom_model,i), grid_size(1), local_field) ! MEVA on t_grid
  ENDDO
  CALL prepare_regridding(idx_MEVA, 0, local_field, my_bottom_model, regrid_u_to_t, regrid_v_to_t, regrid_t_to_u, regrid_t_to_v, grid_size)

  ! HEAT FLUXES
  !   HLAT:
  DO i=1,num_surface_types
    CALL prepare_flux_heat_latent(i, 1, which_flux_heat_latent(my_bottom_model,i), grid_size(1), local_field) ! HLAT on t_grid
  ENDDO
  CALL prepare_regridding(idx_HLAT, 0, local_field, my_bottom_model, regrid_u_to_t, regrid_v_to_t, regrid_t_to_u, regrid_t_to_v, grid_size)
  !   HSEN:
  DO i=1,num_surface_types
    CALL prepare_flux_heat_sensible(i, 1, which_flux_heat_sensible(my_bottom_model,i), grid_size(1), local_field) ! HSEN on t_grid
  ENDDO
  CALL prepare_regridding(idx_HSEN, 0, local_field, my_bottom_model, regrid_u_to_t, regrid_v_to_t, regrid_t_to_u, regrid_t_to_v, grid_size)

  ! RADIATION FLUXES
  !   RBBR:
  DO i=1,num_surface_types
    CALL prepare_flux_radiation_blackbody(i, 1, which_flux_radiation_blackbody(my_bottom_model,i), grid_size(1), local_field) ! RBBR on t_grid
  ENDDO
  CALL prepare_regridding(idx_RBBR, 0, local_field, my_bottom_model, regrid_u_to_t, regrid_v_to_t, regrid_t_to_u, regrid_t_to_v, grid_size)

  ! ! RSDD: redistribution on different surface types
  ! DO i=1,num_surface_types
  !   CALL prepare_distribute_radiation_flux(i, 1, 'test', grid_size(1), local_field) ! RSDD on t_grid
  ! ENDDO
  ! CALL prepare_regridding(idx_RSDD, 0, local_field, my_bottom_model, regrid_u_to_t, regrid_v_to_t, regrid_t_to_u, regrid_t_to_v, grid_size)
  
  ! MOMENTUM FLUXES
  !   UMOM:
  DO i=1,num_surface_types
    CALL prepare_flux_momentum_east(i, 2, which_flux_momentum(my_bottom_model,i), grid_size(2), local_field) ! UMOM on u_grid
  ENDDO
  CALL prepare_regridding(idx_UMOM, 0, local_field, my_bottom_model, regrid_u_to_t, regrid_v_to_t, regrid_t_to_u, regrid_t_to_v, grid_size)
  !   VMOM:
  DO i=1,num_surface_types
    CALL prepare_flux_momentum_north(i, 3, which_flux_momentum(my_bottom_model,i), grid_size(3), local_field) ! VMOM on v_grid
  ENDDO
  CALL prepare_regridding(idx_VMOM, 0, local_field, my_bottom_model, regrid_u_to_t, regrid_v_to_t, regrid_t_to_u, regrid_t_to_v, grid_size)

  !###############################################################################
  !# STEP 1.7:  FIND OUT WHAT I SHALL SEND,                                      #
  !#            AND IF I WILL HAVE CALCULATED EVERYTHING I NEED FOR THAT         #
  !###############################################################################
  ! Count how many output fields we have in total, to find out how large we need to allocate the output array
  ! We may double-count some because they are input fields already, but that does not matter too much, 
  ! since we will later count them exactly. 
  num_output_fields = 0
  DO j=1,MAX_VARS
    myname = name_send_t(j)
    IF (trim(myname) /= 'none') THEN  
      IF (send_uniform_t(my_bottom_model,j)) THEN
        num_output_fields = num_output_fields + 1
      ELSE
        num_output_fields = num_output_fields + (1+num_surface_types)
      ENDIF
    ENDIF
  ENDDO
  DO j=1,MAX_VARS
    myname = name_send_u(j)
    IF (trim(myname) /= 'none') THEN  
      IF (send_uniform_u(my_bottom_model,j)) THEN
        num_output_fields = num_output_fields + 1
      ELSE
        num_output_fields = num_output_fields + (1+num_surface_types)
      ENDIF
    ENDIF
  ENDDO
  DO j=1,MAX_VARS
    myname = name_send_v(j)
    IF (trim(myname) /= 'none') THEN  
      IF (send_uniform_v(my_bottom_model,j)) THEN
        num_output_fields = num_output_fields + 1
      ELSE
        num_output_fields = num_output_fields + (1+num_surface_types)
      ENDIF
    ENDIF
  ENDDO
  ALLOCATE(output_field(num_output_fields))

  ! Now we will call the function add_output_field that will:
  ! (a) check if the required arrays in local_field(1:num_surface_types,:) are already allocated
  ! (b) allocate the arrays in local_field(0,:) (for aggregation over surface types)
  ! (c) set the pointers in output_field(:)
  ! (d) determine the variable names for OASIS
  num_output_fields = 0
  DO j=1,MAX_VARS
    myname = name_send_t(j)
    IF (trim(myname) /= 'none') THEN  
      IF (send_to_atmos_t(j)) THEN
        IF (send_to_bottom_t(my_bottom_model,j)) THEN
          CALL add_output_field(myname, 'A', 0, 1, num_surface_types, grid_size(1), send_uniform_t(my_bottom_model,j), val_flux_t(j), &
          local_field, num_input_fields, input_field, num_output_fields, output_field)
        ELSE
          CALL add_output_field(myname, 'A', 0, 1, num_surface_types, grid_size(1), .FALSE., val_flux_t(j), &
                              local_field, num_input_fields, input_field, num_output_fields, output_field)
        ENDIF
      ENDIF 
      IF (send_to_bottom_t(my_bottom_model,j)) THEN
        IF (send_uniform_t(my_bottom_model,j)) THEN
          CALL add_output_field(myname, my_bottom_letter, 1, 1, num_surface_types, grid_size(1), .TRUE., val_flux_t(j), &
                                local_field, num_input_fields, input_field, num_output_fields, output_field)
        ELSE
          DO i=1,num_surface_types
            CALL add_output_field(myname, my_bottom_letter, i, 1, num_surface_types, grid_size(1), .FALSE., val_flux_t(j), &
                                  local_field, num_input_fields, input_field, num_output_fields, output_field)
          ENDDO
        ENDIF
      ENDIF
    ENDIF
    myname = name_send_u(j)
    IF (trim(myname) /= 'none') THEN  
      IF (send_to_atmos_u(j)) THEN
        IF (send_to_bottom_u(my_bottom_model,j)) THEN
          CALL add_output_field(myname, 'A', 0, 2, num_surface_types, grid_size(2), send_uniform_u(my_bottom_model,j), val_flux_u(j), &
                              local_field, num_input_fields, input_field, num_output_fields, output_field)
        ELSE
          CALL add_output_field(myname, 'A', 0, 2, num_surface_types, grid_size(2), .TRUE., val_flux_u(j), &
                              local_field, num_input_fields, input_field, num_output_fields, output_field)
        ENDIF
      ENDIF 
      IF (send_to_bottom_u(my_bottom_model,j)) THEN
        IF (send_uniform_u(my_bottom_model,j)) THEN
          CALL add_output_field(myname, my_bottom_letter, 1, 2, num_surface_types, grid_size(2), .TRUE., val_flux_u(j), &
                                local_field, num_input_fields, input_field, num_output_fields, output_field)
        ELSE
          DO i=1,num_surface_types
            CALL add_output_field(myname, my_bottom_letter, i, 2, num_surface_types, grid_size(2), .FALSE., val_flux_u(j), &
                                  local_field, num_input_fields, input_field, num_output_fields, output_field)
          ENDDO
        ENDIF
      ENDIF
    ENDIF
    myname = name_send_v(j)
    IF (trim(myname) /= 'none') THEN  
      IF (send_to_atmos_v(j)) THEN
        IF (send_to_bottom_v(my_bottom_model,j)) THEN
          CALL add_output_field(myname, 'A', 0, 3, num_surface_types, grid_size(3), send_uniform_v(my_bottom_model,j), val_flux_v(j), &
          local_field, num_input_fields, input_field, num_output_fields, output_field)
        ELSE
          CALL add_output_field(myname, 'A', 0, 3, num_surface_types, grid_size(3), .TRUE., val_flux_v(j), &
                              local_field, num_input_fields, input_field, num_output_fields, output_field)
        ENDIF
      ENDIF 
      IF (send_to_bottom_v(my_bottom_model,j)) THEN
        IF (send_uniform_v(my_bottom_model,j)) THEN
          CALL add_output_field(myname, my_bottom_letter, 1, 3, num_surface_types, grid_size(3), .TRUE., val_flux_v(j), &
                                local_field, num_input_fields, input_field, num_output_fields, output_field)
        ELSE
          DO i=1,num_surface_types
            CALL add_output_field(myname, my_bottom_letter, i, 3, num_surface_types, grid_size(3), .FALSE., val_flux_v(j), &
                                  local_field, num_input_fields, input_field, num_output_fields, output_field)
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  WRITE(w_unit,*) 'I will send these variables to the coupler:'
  DO i=1,num_output_fields
    WRITE(w_unit,*) '    ',output_field(i)%name,'   on grid ',output_field(i)%which_grid
  ENDDO

  !###############################################################################
  !# STEP 1.7.1:  generate namcouple file                                        #
  !###############################################################################
IF (generate_namcouple) THEN
  CALL create_namcouple(input_field, num_input_fields, output_field, num_output_fields, &
                            name_atmos_model, name_bottom_model, letter_bottom_model,       &
                            timestep, num_timesteps)

  CALL mpi_finalize(ierror)
  STOP ! if we generate the namcouple from here, we are done
ENDIF

  !###############################################################################
  !# STEP 1.8:  OASIS GRID INITIALIZATION                                        #
  !###############################################################################

  IF (mype == 0) THEN
    CALL oasis_start_grids_writing(i)                                           ! i is never used again
    CALL oasis_write_grid('flxt', grid_size_global(1), 1, grid_longitude_global(1)%field, grid_latitude_global(1)%field)
    CALL oasis_write_grid('flxu', grid_size_global(2), 1, grid_longitude_global(2)%field, grid_latitude_global(2)%field)
    CALL oasis_write_grid('flxv', grid_size_global(3), 1, grid_longitude_global(3)%field, grid_latitude_global(3)%field)
    !CALL oasis_write_corner(cgrid, nlon, nlat,nc,globalgrid_clo,globalgrid_cla)
    !CALL oasis_write_mask(cgrid,nlon, nlat,indice_mask(:,:))     
    CALL oasis_write_area('flxt', grid_size_global(1), 1, grid_area_global(1)%field)
    CALL oasis_write_area('flxu', grid_size_global(2), 1, grid_area_global(2)%field)
    CALL oasis_write_area('flxv', grid_size_global(3), 1, grid_area_global(3)%field)
    CALL oasis_terminate_grids_writing()
  ENDIF
  
  !###############################################################################
  !# STEP 1.9:  OASIS VARIABLE INITIALIZATION                                    #
  !###############################################################################
  ! define grid partitions
  DO i=1,3  ! loop over grids
    CALL oasis_def_partition (part_id(i), [1, grid_offset(i), grid_size(i)], ierror)   ! get part_id by putting my part of the parallelization layout
  ENDDO

  ! define input variables
  DO i=1,3  ! loop over grids
    var_nodims(1) = 1    ! Rank of the field array is 1
    var_nodims(2) = 1    ! Bundles always 1 for OASIS3
    var_type = OASIS_Real
    var_actual_shape(1) = 1
    var_actual_shape(2) = grid_size(i)
    DO j=1,num_input_fields
      IF (input_field(j)%which_grid==i) THEN
        CALL oasis_def_var (input_field(j)%id, input_field(j)%name, part_id(i),         &
                            var_nodims, OASIS_In, var_actual_shape, var_type, ierror)
        IF (ierror /= 0) THEN
          WRITE (w_unit,*) 'oasis_def_var abort by flux_calculator compid ',comp_id
          CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_def_var with varname='//TRIM(input_field(j)%name))
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  ! define output variables
  DO i=1,3  ! loop over grids
    var_nodims(1) = 1    ! Rank of the field array is 1
    var_nodims(2) = 1    ! Bundles always 1 for OASIS3
    var_type = OASIS_Real
    var_actual_shape(1) = 1
    var_actual_shape(2) = grid_size(i)
    DO j=1,num_output_fields
      IF (output_field(j)%which_grid==i) THEN
        CALL oasis_def_var (output_field(j)%id, output_field(j)%name, part_id(i),         &
                            var_nodims, OASIS_Out, var_actual_shape, var_type, ierror)
        IF (ierror /= 0) THEN
          WRITE (w_unit,*) 'oasis_def_var abort by flux_calculator compid ',comp_id
          CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_def_var with varname='//TRIM(output_field(j)%name))
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  WRITE (w_unit,*) 'Calling oasis_enddef ...'
  CALL oasis_enddef ( ierror )
  IF (ierror /= 0) THEN
    WRITE (w_unit,*) 'oasis_enddef abort by flux_calculator compid ',comp_id
    CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_enddef')
  ENDIF
  WRITE (w_unit,*) 'Finished initialization.'

  !###############################################################################
  !# STEP 2:  TIME LOOP                                                          #
  !###############################################################################
  
  DO n_timestep = 1,num_timesteps
    current_time = (n_timestep - 1) * timestep
    IF (verbosity_level >= VERBOSITY_LEVEL_STANDARD) THEN
      WRITE (w_unit,*) 'Time since start = ',current_time,' seconds.'
      CALL FLUSH(w_unit)
    ENDIF

    CALL MPI_BARRIER(MPI_COMM_WORLD, i) 
    !#############################################################################
    !# STEP 2.1:  GET EARLY INPUT FROM COUPLER                                   #
    !#############################################################################

    DO i=1,3  ! loop over grids
      DO j=1,num_input_fields
        IF ((input_field(j)%which_grid==i) .AND. (input_field(j)%early==.TRUE.)) THEN
          IF (verbosity_level >= VERBOSITY_LEVEL_DEBUG) THEN
            WRITE (w_unit,*) '  try to get ',input_field(j)%name,' at runtime=',current_time,' seconds.'
          ENDIF
          CALL oasis_get(input_field(j)%id, current_time, input_field(j)%field, ierror)
          IF (verbosity_level >= VERBOSITY_LEVEL_DEBUG) THEN
            WRITE(w_unit,*) '  received  ',input_field(j)%name,' at runtime=',current_time,' seconds:'
            WRITE(w_unit,*) '      range = ',MINVAL(input_field(j)%field),MAXVAL(input_field(j)%field)
          ENDIF
          IF ( ierror .NE. OASIS_Ok .AND. ierror .LT. OASIS_Recvd) THEN
            WRITE (w_unit,*) 'oasis_get abort by flux_calculator compid ',comp_id
            WRITE (w_unit,*) 'could not receive  ',input_field(j)%name,' at runtime=',current_time,' seconds.'
            CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_get')
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    DO j=1,num_input_fields
      IF (input_field(j)%early==.TRUE.) THEN
        CALL do_regridding(input_field(j)%idx, input_field(j)%surface_type, local_field, &
                           regrid_u_to_t_matrix, regrid_v_to_t_matrix, regrid_t_to_u_matrix, regrid_t_to_v_matrix)
      ENDIF
    ENDDO

    !#############################################################################
    !# STEP 2.2:  DO THE EARLY CALCULATIONS                                      #
    !#############################################################################

    CALL calc_flux_radiation_blackbody(my_bottom_model, num_surface_types, which_flux_radiation_blackbody, grid_size, local_field) 
    CALL do_regridding(idx_RBBR, 0, local_field, regrid_u_to_t_matrix, regrid_v_to_t_matrix, regrid_t_to_u_matrix, regrid_t_to_v_matrix)

    !#############################################################################
    !# STEP 2.3:  SEND EARLY RESULTS TO COUPLER                                  #
    !#############################################################################

    DO i=1,3  ! loop over grids
      DO j=1,num_output_fields
        IF ((output_field(j)%which_grid==i) .AND. (output_field(j)%early==.TRUE.)) THEN
          IF (output_field(j)%surface_type == 0) THEN
            IF (ASSOCIATED(local_field(0,i)%var(output_field(j)%idx)%field) .AND. &
                ASSOCIATED(local_field(2,i)%var(output_field(j)%idx)%field) ) THEN
              IF (verbosity_level >= VERBOSITY_LEVEL_DEBUG) THEN
                WRITE (w_unit,*) ' Averaging ',output_field(j)%name,' at runtime=',current_time,' seconds.'
              ENDIF
              CALL average_across_surface_types(i,output_field(j)%idx,num_surface_types,grid_size,local_field)
            ENDIF
          ENDIF
          IF (verbosity_level >= VERBOSITY_LEVEL_DEBUG) THEN
            WRITE (w_unit,*) '  try to put ',output_field(j)%name,' at runtime=',current_time,' seconds.'
            WRITE(w_unit,*) '      range = ',MINVAL(output_field(j)%field),MAXVAL(output_field(j)%field)
          ENDIF
          CALL oasis_put(output_field(j)%id, current_time, output_field(j)%field, ierror)
          IF (verbosity_level >= VERBOSITY_LEVEL_DEBUG) THEN
            WRITE(w_unit,*) '  sent        ',output_field(j)%name,' at runtime=',current_time,' seconds:'
          ENDIF
          IF ( ierror .NE. OASIS_Ok .AND. ierror .LT. OASIS_Recvd) THEN
            WRITE (w_unit,*) 'oasis_get abort by flux_calculator compid ',comp_id
            WRITE (w_unit,*) 'could not send ',output_field(j)%name,' at runtime=',current_time,' seconds.'
            CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_put')
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    !#############################################################################
    !# STEP 2.4:  GET NORMAL INPUT FROM COUPLER                                  #
    !#############################################################################

    DO i=1,3  ! loop over grids
      DO j=1,num_input_fields
        IF ((input_field(j)%which_grid==i) .AND. (input_field(j)%early==.FALSE.)) THEN
          IF (verbosity_level >= VERBOSITY_LEVEL_DEBUG) THEN
            WRITE (w_unit,*) '  try to get ',input_field(j)%name,' at runtime=',current_time,' seconds.'
          ENDIF
          CALL oasis_get(input_field(j)%id, current_time, input_field(j)%field, ierror)
          IF (verbosity_level >= VERBOSITY_LEVEL_DEBUG) THEN
            WRITE(w_unit,*) '  received  ',input_field(j)%name,' at runtime=',current_time,' seconds:'
            WRITE(w_unit,*) '      range = ',MINVAL(input_field(j)%field),MAXVAL(input_field(j)%field)
          ENDIF
          IF ( ierror .NE. OASIS_Ok .AND. ierror .LT. OASIS_Recvd) THEN
            WRITE (w_unit,*) 'oasis_get abort by flux_calculator compid ',comp_id
            WRITE (w_unit,*) 'could not receive  ',input_field(j)%name,' at runtime=',current_time,' seconds.'
            CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_get')
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    DO j=1,num_input_fields
      IF (input_field(j)%early==.FALSE.) THEN
        CALL do_regridding(input_field(j)%idx, input_field(j)%surface_type, local_field, &
                           regrid_u_to_t_matrix, regrid_v_to_t_matrix, regrid_t_to_u_matrix, regrid_t_to_v_matrix)
      ENDIF
    ENDDO

    !#############################################################################
    !# STEP 2.5:  DO THE NORMAL CALCULATIONS                                     #
    !#############################################################################

    CALL calc_spec_vapor_surface(my_bottom_model, num_surface_types, 1, which_spec_vapor_surface_t, grid_size, local_field)
    CALL calc_spec_vapor_surface(my_bottom_model, num_surface_types, 2, which_spec_vapor_surface_u, grid_size, local_field)
    CALL calc_spec_vapor_surface(my_bottom_model, num_surface_types, 3, which_spec_vapor_surface_v, grid_size, local_field)
    CALL do_regridding(idx_QSUR, 0, local_field, regrid_u_to_t_matrix, regrid_v_to_t_matrix, regrid_t_to_u_matrix, regrid_t_to_v_matrix)

    CALL calc_flux_mass_evap(my_bottom_model, num_surface_types, which_flux_mass_evap, grid_size, local_field)
    CALL do_regridding(idx_MEVA, 0, local_field, regrid_u_to_t_matrix, regrid_v_to_t_matrix, regrid_t_to_u_matrix, regrid_t_to_v_matrix)

    CALL calc_flux_heat_latent(my_bottom_model, num_surface_types, which_flux_heat_latent, grid_size, local_field)
    CALL do_regridding(idx_HLAT, 0, local_field, regrid_u_to_t_matrix, regrid_v_to_t_matrix, regrid_t_to_u_matrix, regrid_t_to_v_matrix)

    CALL calc_flux_heat_sensible(my_bottom_model, num_surface_types, which_flux_heat_sensible, grid_size, local_field)
    CALL do_regridding(idx_HSEN, 0, local_field, regrid_u_to_t_matrix, regrid_v_to_t_matrix, regrid_t_to_u_matrix, regrid_t_to_v_matrix)

    CALL calc_flux_momentum_east(my_bottom_model, num_surface_types, 2, which_flux_momentum, grid_size, local_field)
    CALL do_regridding(idx_UMOM, 0, local_field, regrid_u_to_t_matrix, regrid_v_to_t_matrix, regrid_t_to_u_matrix, regrid_t_to_v_matrix)
    CALL calc_flux_momentum_north(my_bottom_model, num_surface_types, 3, which_flux_momentum, grid_size, local_field)
    CALL do_regridding(idx_VMOM, 0, local_field, regrid_u_to_t_matrix, regrid_v_to_t_matrix, regrid_t_to_u_matrix, regrid_t_to_v_matrix)

    CALL distribute_shortwave_radiation_flux(my_bottom_model, num_surface_types, grid_size, local_field)
    !CALL do_regridding(idx_RSDD, 0, local_field, regrid_u_to_t_matrix, regrid_v_to_t_matrix, regrid_t_to_u_matrix, regrid_t_to_v_matrix)


    !#############################################################################
    !# STEP 2.6:  SEND NORMAL RESULTS TO COUPLER                                 #
    !#############################################################################

    DO i=1,3  ! loop over grids
      DO j=1,num_output_fields
        IF ((output_field(j)%which_grid==i) .AND. (output_field(j)%early==.FALSE.)) THEN
          IF (output_field(j)%surface_type == 0) THEN
            IF (ASSOCIATED(local_field(0,i)%var(output_field(j)%idx)%field) .AND. &
                ASSOCIATED(local_field(2,i)%var(output_field(j)%idx)%field) ) THEN
              IF (verbosity_level >= VERBOSITY_LEVEL_DEBUG) THEN
                WRITE (w_unit,*) ' Averaging ',output_field(j)%name,' at runtime=',current_time,' seconds.'
              ENDIF
              CALL average_across_surface_types(i,output_field(j)%idx,num_surface_types,grid_size,local_field)
            ENDIF
          ENDIF
          IF (verbosity_level >= VERBOSITY_LEVEL_DEBUG) THEN
            WRITE (w_unit,*) '  try to put ',output_field(j)%name,' at runtime=',current_time,' seconds.'
            WRITE(w_unit,*) '      range = ',MINVAL(output_field(j)%field),MAXVAL(output_field(j)%field)
          ENDIF
          CALL oasis_put(output_field(j)%id, current_time, output_field(j)%field, ierror)
          IF (verbosity_level >= VERBOSITY_LEVEL_DEBUG) THEN
            WRITE(w_unit,*) '  sent        ',output_field(j)%name,' at runtime=',current_time,' seconds:'
          ENDIF
          IF ( ierror .NE. OASIS_Ok .AND. ierror .LT. OASIS_Recvd) THEN
            WRITE (w_unit,*) 'oasis_get abort by flux_calculator compid ',comp_id
            WRITE (w_unit,*) 'could not send ',output_field(j)%name,' at runtime=',current_time,' seconds.'
            CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_put')
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  
  ENDDO  ! end of time loop
  WRITE (w_unit,*) 'Finished time loop.'

  !###############################################################################
  !# STEP 3:  FINALIZATION                                                       #
  !###############################################################################
  
!   CALL flux_calculator_init()

!   !
!   ! Global grid parameters : 
!   INTEGER :: nlon, nlat     ! dimensions in the 2 directions of space
!   INTEGER :: ntot           ! total dimension
!   INTEGER :: il_paral_size
!   INTEGER :: nc             ! number of corners
!   INTEGER :: indi_beg, indi_end, indj_beg, indj_end
!   CHARACTER(len=4)   :: cgrid='torc'
!   !
!   DOUBLE PRECISION, DIMENSION(:,:), POINTER   :: globalgrid_lon,globalgrid_lat ! lon, lat of the points
!   DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: globalgrid_clo,globalgrid_cla ! lon, lat of the corners
!   DOUBLE PRECISION, DIMENSION(:,:), POINTER   :: globalgrid_srf ! surface of the grid meshes
!   INTEGER, DIMENSION(:,:), POINTER            :: indice_mask ! mask, 0 == valid point, 1 == masked point  
!   !
!   !
!   INTEGER, DIMENSION(:), ALLOCATABLE :: il_paral ! Decomposition for each proc
!   !
!   INTEGER :: ierror, rank, w_unit
!   INTEGER :: i, j
!   INTEGER :: FILE_Debug=2
!   !
!   ! Names of exchanged Fields
!   CHARACTER(len=8), PARAMETER :: var_name1 = 'FSENDOCN' ! 8 characters field sent by model1 to model2
!   CHARACTER(len=8), PARAMETER :: var_name2 = 'FRECVOCN' ! 8 characters field received by model1 from model2
!   !
!   ! Used in oasis_def_var and oasis_def_var
!   INTEGER                   :: var_id(2) 
!   INTEGER                   :: var_nodims(2) 
!   INTEGER                   :: var_type
!   !
!   REAL (kind=wp), PARAMETER :: field_ini = -1. ! initialisation of received fields
!   !
!   INTEGER               ::  ib
!   INTEGER, PARAMETER    ::  il_nb_time_steps = 6 ! number of time steps
!   INTEGER, PARAMETER    ::  delta_t = 3600       ! time step
!   !
!   !
!   INTEGER                 :: il_flag  ! Flag for grid writing by proc 0
!   !
!   INTEGER                 :: itap_sec ! Time used in oasis_put/get
!   !
!   ! Grid parameters definition
!   INTEGER                 :: part_id  ! use to connect the partition to the variables
!                                       ! in oasis_def_var
!   INTEGER                 :: var_actual_shape(4) ! local dimensions of the arrays to the pe
!                                                  ! 2 x field rank (= 4 because fields are of rank = 2)
!   !
!   ! Exchanged local fields arrays
!   ! used in routines oasis_put and oasis_get
!   REAL (kind=wp), POINTER :: field1_send(:,:)
!   REAL (kind=wp), POINTER :: field2_recv(:,:)
!   !
!   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   !   INITIALISATION 
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   !
  
!   !
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   !  GRID DEFINITION 
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   !
!   ! Reading global grid netcdf file
!   !
!   ! Reading dimensions of the global grid
!   CALL read_dimgrid(nlon,nlat,t_grid_filename,w_unit)
!   nc=4
!   !
!   ! Allocation
!   ALLOCATE(globalgrid_lon(nlon,nlat), STAT=ierror )
!   IF ( ierror /= 0 ) WRITE(w_unit,*) 'Error allocating globalgrid_lon'
!   ALLOCATE(globalgrid_lat(nlon,nlat), STAT=ierror )
!   IF ( ierror /= 0 ) WRITE(w_unit,*) 'Error allocating globalgrid_lat'
!   ALLOCATE(globalgrid_clo(nlon,nlat,nc), STAT=ierror )
!   IF ( ierror /= 0 ) WRITE(w_unit,*) 'Error allocating globalgrid_clo'
!   ALLOCATE(globalgrid_cla(nlon,nlat,nc), STAT=ierror )
!   IF ( ierror /= 0 ) WRITE(w_unit,*) 'Error allocating globalgrid_cla'
!   ALLOCATE(globalgrid_srf(nlon,nlat), STAT=ierror )
!   IF ( ierror /= 0 ) WRITE(w_unit,*) 'Error allocating globalgrid_srf'
!   ALLOCATE(indice_mask(nlon,nlat), STAT=ierror )
!   IF ( ierror /= 0 ) WRITE(w_unit,*) 'Error allocating indice_mask'
!   !
!   ! Reading of the longitudes, latitudes, longitude and latitudes of the corners, mask of the global grid
!   CALL read_grid(nlon,nlat,nc, t_grid_filename, w_unit, &
!                  globalgrid_lon,globalgrid_lat, &
!                  globalgrid_clo,globalgrid_cla, &
!                  globalgrid_srf, &
!                  indice_mask)
!   !
!   ! (Global) grid definition for OASIS
!   ! Writing of the file grids.nc and masks.nc by the processor 0 from the grid read in
!   !
!   IF (mype == 0) THEN
!       !
!       ! Mask inversion to follow (historical) OASIS convention (0=not masked;1=masked)
!       WHERE(indice_mask == 1) 
!           indice_mask=0
!       ELSEWHERE
!           indice_mask=1
!       END WHERE
!       !
!       ! TOCOMPLETE - Put here OASIS grid, corner, areas and mask writing calls !
!      CALL oasis_start_grids_writing(il_flag)
!      CALL oasis_write_grid(cgrid,nlon, nlat, globalgrid_lon, globalgrid_lat)
!      CALL oasis_write_corner(cgrid, nlon, nlat,nc,globalgrid_clo,globalgrid_cla)
!      CALL oasis_write_mask(cgrid,nlon, nlat,indice_mask(:,:))     
!      CALL oasis_write_area(cgrid,nlon, nlat,globalgrid_srf)
!      CALL oasis_terminate_grids_writing()
!      !
!      call flush(w_unit)

!   ENDIF
!   !
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   !  PARTITION DEFINITION 
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!   !
!   ! Definition of the partition of the grid (calling oasis_def_partition)
!   ntot=nlon*nlat
! #ifdef DECOMP_APPLE
!   il_paral_size = 3
! #elif defined DECOMP_BOX
!   il_paral_size = 5
! #endif
!   ALLOCATE(il_paral(il_paral_size))
!   WRITE(w_unit,*) 'After allocate il_paral, il_paral_size', il_paral_size
!   call flush(w_unit)
!   !
!   CALL decomp_def (il_paral,il_paral_size,nlon,nlat,mype,npes,w_unit)
!   WRITE(w_unit,*) 'After decomp_def, il_paral = ', il_paral(:)
!   call flush(w_unit)

!   CALL oasis_def_partition (part_id, il_paral, ierror)   ! get part_id by putting my part of the parallelization layout

!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! DEFINITION OF THE LOCAL FIELDS  
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   !
!   !!!!!!!!!!!!!!! !!!!!!!!! OASIS_DEF_VAR !!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !
!   !  Define transient variables
!   !
!   var_nodims(1) = 2    ! Rank of the field array is 2
!   var_nodims(2) = 1    ! Bundles always 1 for OASIS3
!   var_type = OASIS_Real
!   !
!   var_actual_shape(1) = 1
!   var_actual_shape(2) = il_paral(3)
!   var_actual_shape(3) = 1 
! #ifdef DECOMP_APPLE
!   var_actual_shape(4) = 1
! #elif defined DECOMP_BOX
!   var_actual_shape(4) = il_paral(4)
! #endif
!   !
!   ! Declaration of the field associated with the partition
!   !
!   ! TOCOMPLETE - Put here OASIS call to declare the coupling fields
!   !              FRECVOCN, FSENDOCN
!   ! var_name1 = 'FSENDOCN'
!   ! var_name2 = 'FRECVOCN'
!   CALL oasis_def_var (var_id(1),var_name1, part_id, &
!      var_nodims, OASIS_Out, var_actual_shape, var_type, ierror)
!   IF (ierror /= 0) THEN
!       WRITE (w_unit,*) 'oasis_def_var abort by flux_calculator compid ',comp_id
!       CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_def_var with varname='//TRIM(var_name1))
!   ENDIF
!   CALL oasis_def_var (var_id(2),var_name2, part_id, &
!      var_nodims, OASIS_In, var_actual_shape, var_type, ierror)
!   IF (ierror /= 0) THEN
!       WRITE (w_unit,*) 'oasis_def_var abort by flux_calculator compid ',comp_id
!       CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_def_var with varname='//TRIM(var_name2))
!   ENDIF

!   !
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   !         TERMINATION OF DEFINITION PHASE 
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   !  All processes involved in the coupling must call oasis_enddef; 
!   !  here all processes are involved in coupling
!   !
!   !!!!!!!!!!!!!!!!!! OASIS_ENDDEF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!   CALL oasis_enddef ( ierror )
!   IF (ierror /= 0) THEN
!       WRITE (w_unit,*) 'oasis_enddef abort by flux_calculator compid ',comp_id
!       CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_enddef')
!   ENDIF
!   !
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! SEND AND RECEIVE ARRAYS 
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   !
!   ! Allocate the fields send and received by the model
!   !
!   !
!   ALLOCATE(field1_send(var_actual_shape(2), var_actual_shape(4)), STAT=ierror )
!   IF ( ierror /= 0 ) WRITE(w_unit,*) 'Error allocating field1_send'
!   !
!   ALLOCATE(field2_recv(var_actual_shape(2), var_actual_shape(4)), STAT=ierror )
!   IF ( ierror /= 0 ) WRITE(w_unit,*) 'Error allocating field2_recv'
!   !
!   DEALLOCATE(il_paral)
!   !
!   !!!!!!!!!!!!!!!!!!!!!!!!OASIS_PUT/OASIS_GET !!!!!!!!!!!!!!!!!!!!!! 
!   !
!   indi_beg=1 ; indi_end=nlon
!   indj_beg=((nlat/npes)*mype)+1 
!   !
!   IF (mype .LT. npes - 1) THEN
!       indj_end = (nlat/npes)*(mype+1)
!   ELSE
!       indj_end = nlat 
!   ENDIF
!   !
!   ! Data exchange 
!   ! 
!   ! Time loop
!   DO ib=1, il_nb_time_steps
!     itap_sec = delta_t * (ib-1) ! Time
!     !
!     ! Get FRECVOCN
!     ! TOCOMPLETE - Put here the OASIS call to receive FRECVOCN (field2_recv)
!     ! Let's suppose here that FRECVOCN contains BC needed for the timestep

!     field2_recv=field_ini
!     CALL oasis_get(var_id(2),itap_sec, field2_recv, ierror)
!     IF (FILE_Debug >= 2) THEN
!         WRITE(w_unit,*) 'tcx recvf1 ',itap_sec,MINVAL(field2_recv),MAXVAL(field2_recv)
!     ENDIF
!     IF ( ierror .NE. OASIS_Ok .AND. ierror .LT. OASIS_Recvd) THEN
!         WRITE (w_unit,*) 'oasis_get abort by flux_calculator compid ',comp_id
!         CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_get')
!     ENDIF
!     !
!     !

!     !
!     ! Here the model computes its timestep
!     !
!     CALL function_sent(var_actual_shape(2), var_actual_shape(4), &
!                        RESHAPE(globalgrid_lon(indi_beg:indi_end,indj_beg:indj_end),&
!                                (/ var_actual_shape(2), var_actual_shape(4) /)), &
!                        RESHAPE(globalgrid_lat(indi_beg:indi_end,indj_beg:indj_end),&
!                                (/ var_actual_shape(2), var_actual_shape(4) /)), &
!                                 field1_send,ib)
!     !
!     ! Send FSENDOCN

!     IF (FILE_Debug >= 2) THEN
!         WRITE(w_unit,*) 'tcx sendf ',itap_sec,MINVAL(field1_send),MAXVAL(field1_send)
!     ENDIF
!     CALL oasis_put(var_id(1),itap_sec, field1_send, ierror)
!     IF ( ierror .NE. OASIS_Ok .AND. ierror .LT. OASIS_Sent) THEN
!         WRITE (w_unit,*) 'oasis_put abort by model1 compid ',comp_id
!         CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_put')
!     ENDIF
!     !
!     !
!   ENDDO
!   !
!   WRITE (w_unit,*) 'End of the program'
!   CALL flush(w_unit)
!   CLOSE (w_unit)
!   !
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   !         TERMINATION 
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   !
!   !!!!!!!!!!!!!!!!!! OASIS_TERMINATE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !
!   ! Collective call to terminate the coupling exchanges
!   !
!   IF (FILE_Debug >= 2) THEN
!       WRITE (w_unit,*) 'End of the program, before oasis_terminate'
!       CALL FLUSH(w_unit)
!   ENDIF
  !
  !!!!!!!!!!!!!!!!!! OASIS_ENDDEF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Collective call to terminate the coupling exchanges
  !
  CALL oasis_terminate (ierror)
  IF (ierror /= 0) THEN
      WRITE (w_unit,*) 'oasis_terminate abort by flux_calculator compid ',comp_id
      CALL oasis_abort(comp_id,comp_name,'Failed to call oasis_terminate')
  ENDIF
  !
  !
  ! TOCOMPLETE - Put here the OASIS call to terminate the coupling
  !
  CALL mpi_finalize(ierror)

END PROGRAM flux_calculator
!
