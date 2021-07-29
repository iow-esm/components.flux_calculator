MODULE flux_calculator_basic
! This module contains:
!   - type definitions
!   - size limits for arrays
!   - the list of variable names and their indexes
!   - basic auxiliary functions

    IMPLICIT NONE

!!!!!!!!!! FUNCTIONS DEFINED IN THIS MODULE
    PUBLIC add_input_field
    PUBLIC allocate_localvar
    PUBLIC distribute_input_field
    PUBLIC do_regridding
    PUBLIC init_localvar
    PUBLIC init_varname_idx
    PUBLIC nullify_localvars
    PUBLIC prepare_regridding

!!!!!!!!!! NOW EVERYTHING ELSE
#ifdef NO_USE_DOUBLE_PRECISION
    INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(6,37)   ! real
#elif defined USE_DOUBLE_PRECISION
    INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307) ! double
#endif

    INTEGER, PARAMETER :: MAX_BOTTOM_MODELS   = 10       ! how many models we can couple to the atmospheric model
    INTEGER, PARAMETER :: MAX_SURFACE_TYPES   = 10       ! how many surface types (land or ice classes) are allowed
    INTEGER, PARAMETER :: MAX_TASKS_PER_MODEL = 1000     ! on how many PEs per bottom model can the flux calculator run
    INTEGER, PARAMETER :: MAX_VARS            = 100      ! how many variables can a model send/receive on each grid
    INTEGER, PARAMETER :: MAX_INPUT_FIELDS    = 1000     ! how many variables can flux_calculator receive in total
    INTEGER, PARAMETER :: MAX_OUTPUT_FIELDS   = 1000     ! how many variables can flux_calculator send in total


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! If you add a varname ????, you need to do 4 changes:         !
    ! (1) increase MAX_VARNAMES                                    !
    ! (2) add '????' to varnames                                   !
    ! (3) declare idx_varname                                      !
    ! (4) add a line in init_varname_idx                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER, PARAMETER                                   :: MAX_VARNAMES = 32
    CHARACTER(len=4), PARAMETER, DIMENSION(MAX_VARNAMES) :: varnames = [ &
        'ALBE', 'ALBA', 'AMOI', 'AMOM', 'FARE', 'FICE', 'PATM', 'PSUR', &
        'QATM', 'TATM', 'TSUR', 'UATM', 'VATM', 'U10M', 'V10M', &  ! variables read in
        'QSUR',                                                 &  ! auxiliary variables calculated
        'HLAT', 'HSEN',                                         &  ! heat fluxes
        'MEVA', 'MPRE', 'MRAI', 'MSNO',                         &  ! mass fluxes
        'RBBR', 'RLWD', 'RLWU', 'RSID', 'RSIU', 'RSIN', 'RSDD', 'RSDR', &  ! radiation fluxes, the last are: Shortwave_Indirect_Down, Shortwave_Indirect_Up, Shortwave_Indirect_Net, Shortwave_Direct_Down
        'UMOM', 'VMOM']                                            ! momentum fluxes

    INTEGER :: idx_ALBE, idx_ALBA, idx_AMOI, idx_AMOM, idx_FARE, idx_FICE, idx_PATM, idx_PSUR
    INTEGER :: idx_QATM, idx_TATM, idx_TSUR, idx_UATM, idx_VATM, idx_U10M, idx_V10M
    INTEGER :: idx_QSUR
    INTEGER :: idx_HLAT, idx_HSEN
    INTEGER :: idx_MEVA, idx_MPRE, idx_MRAI, idx_MSNO
    INTEGER :: idx_RBBR, idx_RLWD, idx_RLWU, idx_RSID, idx_RSIU, idx_RSIN, idx_RSDD, idx_RSDR
    INTEGER :: idx_UMOM, idx_VMOM

    CHARACTER(len=2), DIMENSION(0:MAX_SURFACE_TYPES)  :: numtype

    CHARACTER(len=6), DIMENSION(3) :: grid_name = ['t_grid', 'u_grid', 'v_grid']

    INTEGER :: w_unit  ! a logfile to write the progress and error messages

      ! control amount of output
    ENUM, BIND(c)
        ENUMERATOR :: VERBOSITY_LEVEL_ERROR = 0
        ENUMERATOR :: VERBOSITY_LEVEL_STANDARD
        ENUMERATOR :: VERBOSITY_LEVEL_DEBUG
    ENDENUM

#ifdef IOW_ESM_DEBUG
    INTEGER :: verbosity_level = VERBOSITY_LEVEL_DEBUG   ! 1 = standard, 2 = debug output
#else
    INTEGER :: verbosity_level = VERBOSITY_LEVEL_STANDARD   ! 1 = standard, 2 = debug output
#endif
  
    TYPE integerarray
        INTEGER (kind=4), DIMENSION(:), POINTER :: field     
        LOGICAL                                 :: allocated
    END TYPE integerarray

    TYPE realarray
        REAL (kind=wp), DIMENSION(:), POINTER :: field      
        LOGICAL                               :: allocated
        LOGICAL                               :: put_to_t_grid   ! whether this shall be regridded to the t_grid after reading/calculation here
        LOGICAL                               :: put_to_u_grid   ! whether this shall be regridded to the u_grid after reading/calculation here
        LOGICAL                               :: put_to_v_grid   ! whether this shall be regridded to the v_grid after reading/calculation here
    END TYPE realarray

    TYPE realarray2
        REAL (kind=wp), DIMENSION(:,:), POINTER :: field    
        LOGICAL                                 :: allocated
    END TYPE realarray2
    
    ! define a type to store the variables
    TYPE local_fields_type
        ! will be allocated with number of gridcells
        TYPE(realarray), DIMENSION(MAX_VARNAMES) :: var             ! these arrays will be allocated if the variable is present  
    END TYPE local_fields_type
    
    ! define a type to store the input/output fields
    TYPE io_fields_type
        REAL (kind=wp), DIMENSION(:), POINTER   :: field      ! arrays will just be pointers to local_fields arrays 
        CHARACTER(len=8)                        :: name       ! will store the OASIS variable name
        INTEGER                                 :: which_grid ! 1=t_grid, 2=u_grid, 3=v_grid
        INTEGER                                 :: id         ! the id returned by oasis_def_var
        LOGICAL                                 :: early      ! whether or not we communicate these variable before all others
        INTEGER                                 :: surface_type ! surface type, used for regridding of input data
        INTEGER                                 :: idx          ! idx_???? number of variable, used for regridding of input data
    END TYPE io_fields_type

    ! define a type to store a regridding matrix (e.g. t_grid to u_grid)
    TYPE sparse_regridding_matrix
        INTEGER            :: num_elements
        TYPE(integerarray) :: src_index
        TYPE(integerarray) :: dst_index
        TYPE(realarray)    :: weight
    END TYPE sparse_regridding_matrix

    CONTAINS
    
    SUBROUTINE add_input_field(myname, my_letter, surface_type, which_grid, local_field, num_input_fields, input_field)
        CHARACTER(len=4),                      INTENT(IN)    :: myname               ! e.g. "PATM"
        CHARACTER(len=1),                      INTENT(IN)    :: my_letter     ! 'A' for atmosphere, my_bottom_letter for bottom model
        INTEGER,                               INTENT(IN)    :: surface_type
        INTEGER,                               INTENT(IN)    :: which_grid           ! 1=t_grid, 2=u_grid, 3=v_grid
        TYPE(local_fields_type),               INTENT(INOUT) :: local_field          ! here we want to set the pointer to
        INTEGER,                               INTENT(INOUT) :: num_input_fields     
        type(io_fields_type),    DIMENSION(:), INTENT(INOUT) :: input_field
        INTEGER                                :: i
        LOGICAL                                :: found
        CHARACTER(len=2)                       :: appendstring         ! two-digit number of surface type
        appendstring = numtype(surface_type)
        found=.FALSE.
        DO i=1,MAX_VARNAMES
            IF (varnames(i) == myname) THEN
                found=.TRUE.
                num_input_fields = num_input_fields + 1
                IF (num_input_fields > MAX_INPUT_FIELDS) THEN
                    WRITE (w_unit,*) "Maximum number of input fields exceeded. Consider increasing MAX_INPUT_FIELDS and then recompile flux_calculator."
                    CALL mpi_finalize(1)
                ELSE
                    input_field(num_input_fields)%name="R"//my_letter//myname//appendstring
                    input_field(num_input_fields)%which_grid = which_grid
                    input_field(num_input_fields)%field => local_field%var(i)%field
                    input_field(num_input_fields)%early=.FALSE.
                    IF ((myname=="FARE") .OR. (myname=="TSUR") .OR. (myname=="ALBE")) THEN
                        input_field(num_input_fields)%early=.TRUE.
                    ENDIF
                    input_field(num_input_fields)%surface_type = surface_type
                    input_field(num_input_fields)%idx = i
                ENDIF
            ENDIF
        ENDDO
        IF (.NOT. found) THEN
            WRITE (w_unit,*) "Could not add input field for variable ",myname," because flux_calculator does not know this variable."
            CALL mpi_finalize(1)
        ENDIF
    END SUBROUTINE add_input_field

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    SUBROUTINE add_output_field(myname, my_letter, surface_type, which_grid, num_surface_types, grid_size, &
                                uniform, default_value, local_field, num_input_fields, input_field, num_output_fields, output_field)
        CHARACTER(len=4),                         INTENT(IN)    :: myname               ! e.g. "PATM"
        CHARACTER(len=1),                         INTENT(IN)    :: my_letter            ! "A" for atmospheric model, my_bottom_letter for bottom model
        INTEGER,                                  INTENT(IN)    :: surface_type         ! 0=entire cell, 1,2,3,... = specific surface type
        INTEGER,                                  INTENT(IN)    :: which_grid           ! 1=t_grid, 2=u_grid, 3=v_grid
        INTEGER,                                  INTENT(IN)    :: num_surface_types    ! how many does my bottom model have
        INTEGER,                                  INTENT(IN)    :: grid_size           
        LOGICAL,                                  INTENT(IN)    :: uniform              ! whether this flux is the same for each surface_type
        REAL(kind=wp),                            INTENT(IN)    :: default_value        ! if no method calculated the flux
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field          ! here we want to set the pointer to
        INTEGER,                                  INTENT(IN)    :: num_input_fields     
        type(io_fields_type),    DIMENSION(:),    INTENT(IN)    :: input_field
        INTEGER,                                  INTENT(INOUT) :: num_output_fields     
        type(io_fields_type),    DIMENSION(:),    INTENT(INOUT) :: output_field
        INTEGER                                :: i,j
        LOGICAL                                :: found, found_all_fluxes, found_all_areas
        CHARACTER(len=2)                       :: appendstring         ! two-digit number of surface type
        appendstring=numtype(surface_type)
        ! First check if this already exists as input field. In this case we do not need the output field.
        found=.FALSE.
        DO i=1,num_input_fields
          IF (trim(input_field(i)%name)=="R"//my_letter//myname//appendstring) found=.TRUE.
        ENDDO
        IF (.NOT. found) THEN ! Okay field does not come as input
            DO i=1,MAX_VARNAMES
                IF (varnames(i) == myname) THEN
                    found=.TRUE.
                    ! Do we want output for surface_type=0 (the entire exchange grid cell)?
                    IF (surface_type==0) THEN
                        ! Yes, output for the entire cell
                        IF (uniform) THEN
                            ! Flux is the same everywhere so we only need to find one surface_type for which it is calculated
                            DO j=1,num_surface_types
                                IF (ASSOCIATED(local_field(j,which_grid)%var(i)%field) .AND. (.NOT. ASSOCIATED(local_field(0,which_grid)%var(i)%field))) THEN
                                    local_field(0,which_grid)%var(i)%field => local_field(j,which_grid)%var(i)%field
                                ENDIF
                            ENDDO
                        ELSE
                            ! Flux will differ, so we will need a value for each surface_type and an area fraction of that surface_type
                            found_all_fluxes = .TRUE.
                            found_all_areas  = .TRUE.
                            DO j=1,num_surface_types
                                IF (.NOT. ASSOCIATED(local_field(j,which_grid)%var(i)%field))        found_all_fluxes=.FALSE.
                                IF (.NOT. ASSOCIATED(local_field(j,which_grid)%var(idx_FARE)%field)) found_all_areas=.FALSE.
                            ENDDO
                            IF (.NOT. found_all_areas) THEN
                                ! fractional area is not provided - this is fatal
                                WRITE (w_unit,*) "ERROR: Output field ",myname," has not been defined as uniform (flux_?_uniform=.FALSE.)."
                                WRITE (w_unit,*) "To calculate its average value across different surface_types, their fractional area (FARE) must be given but is missing."
                                CALL mpi_finalize(1)
                            ELSE 
                                IF (found_all_fluxes) THEN
                                    ! this is nice, we found everything we need, so just allocate
                                    IF (.NOT. ASSOCIATED(local_field(0,which_grid)%var(i)%field)) THEN
                                        ALLOCATE(local_field(0,which_grid)%var(i)%field(grid_size))
                                        local_field(0,which_grid)%var(i)%allocated=.TRUE.
                                    ENDIF
                                ENDIF
                            ENDIF
                        ENDIF
                        IF (.NOT. ASSOCIATED(local_field(0,which_grid)%var(i)%field)) THEN
                            ! Okay we have not found a pointer yet or allocated the field we want to send.
                            ! But we also have not encountered a serious error. 
                            ! So we may just be missing the calculation of this flux.
                            ! This will just give a warning message because the fluxes have standard values (typically zero)
                            WRITE (w_unit,*) "WARNING: Flux ",myname," cannot be calculated for surface_type=0  =>  set to ",default_value
                            ALLOCATE(local_field(0,which_grid)%var(i)%field(grid_size))
                            local_field(0,which_grid)%var(i)%field = default_value
                            local_field(0,which_grid)%var(i)%allocated=.TRUE.
                        ENDIF
                    ELSE
                        IF (uniform) THEN
                            IF (.NOT. ASSOCIATED(local_field(surface_type,which_grid)%var(i)%field)) THEN
                                !okay it is not calculated for this surface_type, but probably for another surface_type, so we can use its values
                                DO j=1,num_surface_types
                                    IF (ASSOCIATED(local_field(j,which_grid)%var(i)%field) .AND. (.NOT. ASSOCIATED(local_field(0,which_grid)%var(i)%field))) THEN
                                        local_field(0,which_grid)%var(i)%field => local_field(j,which_grid)%var(i)%field
                                    ENDIF
                                ENDDO
                            ENDIF
                        ENDIF
                        ! We want output for one specific surface type only
                        IF (.NOT. ASSOCIATED(local_field(surface_type,which_grid)%var(i)%field)) THEN
                            ! We are missing the calculation of this flux.
                            ! This will just give a warning message because the fluxes have standard values (typically zero)
                            WRITE (w_unit,*) "WARNING: Flux ",myname," cannot be calculated for surface_type=",surface_type,"  =>  set to ",default_value
                            ALLOCATE(local_field(surface_type,which_grid)%var(i)%field(grid_size))
                            local_field(surface_type,which_grid)%var(i)%field = default_value
                            local_field(surface_type,which_grid)%var(i)%allocated=.TRUE.
                        ENDIF
                    ENDIF
                    num_output_fields = num_output_fields + 1
                    IF (num_output_fields > MAX_OUTPUT_FIELDS) THEN
                        WRITE (w_unit,*) "Maximum number of output fields exceeded. Consider increasing MAX_OUTPUT_FIELDS and then recompile flux_calculator."
                        CALL mpi_finalize(1)
                    ELSE
                        output_field(num_output_fields)%name="S"//my_letter//myname//appendstring
                        output_field(num_output_fields)%which_grid = which_grid
                        output_field(num_output_fields)%field => local_field(surface_type,which_grid)%var(i)%field
                        output_field(num_output_fields)%early=.FALSE.
                        IF ((myname=="RBBR") .OR. (myname=="TSUR") .OR. (myname=="FICE") .OR. (myname=="ALBE")) THEN
                            output_field(num_output_fields)%early=.TRUE.
                        ENDIF
                        output_field(num_output_fields)%surface_type = surface_type
                        output_field(num_output_fields)%idx = i
                    ENDIF
                ENDIF
            ENDDO
            IF (.NOT. found) THEN
                WRITE (w_unit,*) "Could not add output field for variable ",myname," because flux_calculator does not know this variable."
                CALL mpi_finalize(1)
            ENDIF
        ENDIF
    END SUBROUTINE add_output_field

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    SUBROUTINE allocate_localvar(myname, length, local_field)
        CHARACTER(len=4),        INTENT(IN)    :: myname  ! e.g. "PATM"
        INTEGER,                 INTENT(IN)    :: length
        TYPE(local_fields_type), INTENT(INOUT) :: local_field
        INTEGER                                :: i
        LOGICAL                                :: found
        found=.FALSE.
        DO i=1,MAX_VARNAMES
            IF (varnames(i) == myname) THEN
                found=.TRUE.
                ALLOCATE(local_field%var(i)%field(length))
                local_field%var(i)%allocated = .TRUE.
                local_field%var(i)%put_to_t_grid    = .FALSE.
                local_field%var(i)%put_to_u_grid    = .FALSE.
                local_field%var(i)%put_to_v_grid    = .FALSE.
            ENDIF
        ENDDO
        IF (.NOT. found) THEN
            WRITE (w_unit,*) "Could not allocate local variable ",myname," because flux_calculator does not know this variable."
            CALL mpi_finalize(1)
        ENDIF
    END SUBROUTINE allocate_localvar

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    SUBROUTINE init_localvar(myname, value, local_field)
        CHARACTER(len=4),        INTENT(IN)    :: myname  ! e.g. "PATM"
        REAL(kind=wp),           INTENT(IN)    :: value
        TYPE(local_fields_type), INTENT(INOUT) :: local_field
        INTEGER                                :: i
        LOGICAL                                :: found
        found=.FALSE.
        DO i=1,MAX_VARNAMES
            IF (varnames(i) == myname) THEN
                found=.TRUE.
                local_field%var(i)%field = value
            ENDIF
        ENDDO
        IF (.NOT. found) THEN
            WRITE (w_unit,*) "Could not set local variable ",myname," to constant value because flux_calculator does not know this variable."
            CALL mpi_finalize(1)
        ENDIF
    END SUBROUTINE init_localvar

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    SUBROUTINE distribute_input_field(myname, which_grid, from_surface_type, to_surface_type, num_surface_types, local_field)
        CHARACTER(len=4),        INTENT(IN)                     :: myname  ! e.g. "PATM"
        INTEGER,                 INTENT(IN)                     :: which_grid
        INTEGER,                 INTENT(IN)                     :: from_surface_type
        INTEGER,                 INTENT(IN)                     :: to_surface_type
        INTEGER,                 INTENT(IN)                     :: num_surface_types
        TYPE(local_fields_type), INTENT(INOUT), DIMENSION(0:,:) :: local_field
        INTEGER :: i, j
        LOGICAL :: found
        found=.FALSE.
        DO i=1,MAX_VARNAMES
            IF (varnames(i) == myname) THEN
                found=.TRUE.
                DO j=1,num_surface_types
                    IF ((j == to_surface_type) .OR. ((j /= from_surface_type) .AND. (to_surface_type == 0))) THEN
                        local_field(j,which_grid)%var(i)%field => local_field(from_surface_type,which_grid)%var(i)%field
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        IF (.NOT. found) THEN
            WRITE (w_unit,*) "Could not distribute local variable ",myname," to other surface_types because flux_calculator does not know this variable."
            CALL mpi_finalize(1)
        ENDIF
    END SUBROUTINE distribute_input_field

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    SUBROUTINE prepare_regridding(varidx, surface_type, local_field, my_bottom_model,                 &
                                  regrid_u_to_t, regrid_v_to_t, regrid_t_to_u, regrid_t_to_v,         &
                                  grid_size)
        INTEGER,                                   INTENT(IN)    :: varidx
        INTEGER,                                   INTENT(IN)    :: surface_type ! set to 0 for all surface types
        TYPE(local_fields_type), DIMENSION(0:,:),  INTENT(INOUT) :: local_field
        INTEGER,                                   INTENT(IN)    :: my_bottom_model
        CHARACTER(len=4),        DIMENSION(:,:,:), INTENT(IN)    :: regrid_u_to_t
        CHARACTER(len=4),        DIMENSION(:,:,:), INTENT(IN)    :: regrid_v_to_t
        CHARACTER(len=4),        DIMENSION(:,:,:), INTENT(IN)    :: regrid_t_to_u
        CHARACTER(len=4),        DIMENSION(:,:,:), INTENT(IN)    :: regrid_t_to_v
        INTEGER, DIMENSION(:),                     INTENT(IN)    :: grid_size

        INTEGER :: i, j, k
        INTEGER :: from_grid, to_grid
        
        i=varidx
        from_grid=2; to_grid=1
        DO j=1,MAX_SURFACE_TYPES   
            IF ((j==surface_type) .OR. (surface_type==0)) THEN
                DO k=1,MAX_VARS
                    IF (varnames(i) == trim(regrid_u_to_t(my_bottom_model,j,k))) THEN
                        IF (local_field(j,to_grid)%var(i)%allocated) THEN
                            ! This should not be the case - the variable already exists on the destination grid
                            WRITE (w_unit,*) "Could not regrid local variable ",varnames(i)," from ",grid_name(from_grid)," to ",grid_name(to_grid),&
                                                " as requested in the namelist, because it already exists on that grid."
                            CALL mpi_finalize(1)
                        ELSE
                            ALLOCATE(local_field(j,to_grid)%var(i)%field(grid_size(to_grid)))
                            local_field(j,to_grid)%var(i)%allocated = .TRUE.
                            local_field(j,from_grid)%var(i)%put_to_t_grid  = .TRUE.
                            WRITE(w_unit,*) '    u_grid -> t_grid: ',varnames(i),' for surface type ',j
                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        from_grid=3; to_grid=1
        DO j=1,MAX_SURFACE_TYPES 
            IF ((j==surface_type) .OR. (surface_type==0)) THEN
                DO k=1,MAX_VARS
                    IF (varnames(i) == trim(regrid_v_to_t(my_bottom_model,j,k))) THEN
                        IF (local_field(j,to_grid)%var(i)%allocated) THEN
                            ! This should not be the case - the variable already exists on the destination grid
                            WRITE (w_unit,*) "Could not regrid local variable ",varnames(i)," from ",grid_name(from_grid)," to ",grid_name(to_grid),&
                                                " as requested in the namelist, because it already exists on that grid."
                            CALL mpi_finalize(1)
                        ELSE
                            ALLOCATE(local_field(j,to_grid)%var(i)%field(grid_size(to_grid)))
                            local_field(j,to_grid)%var(i)%allocated = .TRUE.
                            local_field(j,from_grid)%var(i)%put_to_t_grid  = .TRUE.
                            WRITE(w_unit,*) '    v_grid -> t_grid: ',varnames(i),' for surface type ',j
                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        from_grid=1; to_grid=2
        DO j=1,MAX_SURFACE_TYPES  
            IF ((j==surface_type) .OR. (surface_type==0)) THEN
                DO k=1,MAX_VARS
                    IF (varnames(i) == trim(regrid_t_to_u(my_bottom_model,j,k))) THEN
                        IF (local_field(j,to_grid)%var(i)%allocated) THEN
                            ! This should not be the case - the variable already exists on the destination grid
                            WRITE (w_unit,*) "Could not regrid local variable ",varnames(i)," from ",grid_name(from_grid)," to ",grid_name(to_grid),&
                                                " as requested in the namelist, because it already exists on that grid."
                            CALL mpi_finalize(1)
                        ELSE
                            ALLOCATE(local_field(j,to_grid)%var(i)%field(grid_size(to_grid)))
                            local_field(j,to_grid)%var(i)%allocated = .TRUE.
                            local_field(j,from_grid)%var(i)%put_to_u_grid  = .TRUE.
                            WRITE(w_unit,*) '    t_grid -> u_grid: ',varnames(i),' for surface type ',j
                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        from_grid=1; to_grid=3
        DO j=1,MAX_SURFACE_TYPES 
            IF ((j==surface_type) .OR. (surface_type==0)) THEN
                DO k=1,MAX_VARS
                    IF (varnames(i) == trim(regrid_t_to_v(my_bottom_model,j,k))) THEN
                        IF (local_field(j,to_grid)%var(i)%allocated) THEN
                            ! This should not be the case - the variable already exists on the destination grid
                            WRITE (w_unit,*) "Could not regrid local variable ",varnames(i)," from ",grid_name(from_grid)," to ",grid_name(to_grid),&
                                                " as requested in the namelist, because it already exists on that grid."
                            CALL mpi_finalize(1)
                        ELSE
                            ALLOCATE(local_field(j,to_grid)%var(i)%field(grid_size(to_grid)))
                            local_field(j,to_grid)%var(i)%allocated = .TRUE.
                            local_field(j,from_grid)%var(i)%put_to_v_grid  = .TRUE.
                            WRITE(w_unit,*) '    t_grid -> v_grid: ',varnames(i),' for surface type ',j
                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
        ENDDO 
    END SUBROUTINE prepare_regridding

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    SUBROUTINE do_regridding(varidx, surface_type, local_field, &
                             regrid_u_to_t_matrix, regrid_v_to_t_matrix, regrid_t_to_u_matrix, regrid_t_to_v_matrix)
        INTEGER,                                   INTENT(IN)    :: varidx
        INTEGER,                                   INTENT(IN)    :: surface_type  ! 0 stands for all surface types
        TYPE(local_fields_type), DIMENSION(0:,:),   INTENT(INOUT) :: local_field
        TYPE(sparse_regridding_matrix),            INTENT(IN)    :: regrid_u_to_t_matrix
        TYPE(sparse_regridding_matrix),            INTENT(IN)    :: regrid_v_to_t_matrix
        TYPE(sparse_regridding_matrix),            INTENT(IN)    :: regrid_t_to_u_matrix
        TYPE(sparse_regridding_matrix),            INTENT(IN)    :: regrid_t_to_v_matrix

        INTEGER :: i, j, k, from, to
        INTEGER :: from_grid, to_grid
        
        i=varidx
        DO j=1,MAX_SURFACE_TYPES
            IF ((j==surface_type) .OR. (surface_type==0)) THEN
                IF (local_field(j,2)%var(i)%put_to_t_grid) THEN
                    from=2; to=1
                    local_field(j,to)%var(i)%field = 0.0
                    DO k=1,regrid_u_to_t_matrix%num_elements
                        local_field(j,to)%var(i)%field(regrid_u_to_t_matrix%dst_index%field(k)) =       & ! sparse matrix multiplication
                            local_field(j,to)%var(i)%field(regrid_u_to_t_matrix%dst_index%field(k)) +   & 
                            local_field(j,from)%var(i)%field(regrid_u_to_t_matrix%src_index%field(k)) * & 
                            regrid_u_to_t_matrix%weight%field(k)
                    ENDDO
                ENDIF
                IF (local_field(j,3)%var(i)%put_to_t_grid) THEN
                    from=3; to=1
                    local_field(j,to)%var(i)%field = 0.0
                    DO k=1,regrid_v_to_t_matrix%num_elements
                        local_field(j,to)%var(i)%field(regrid_v_to_t_matrix%dst_index%field(k)) =       & ! sparse matrix multiplication
                            local_field(j,to)%var(i)%field(regrid_v_to_t_matrix%dst_index%field(k)) +   & 
                            local_field(j,from)%var(i)%field(regrid_v_to_t_matrix%src_index%field(k)) * & 
                            regrid_v_to_t_matrix%weight%field(k)
                    ENDDO
                ENDIF
                IF (local_field(j,1)%var(i)%put_to_u_grid) THEN
                    from=1; to=2
                    local_field(j,to)%var(i)%field = 0.0
                    DO k=1,regrid_t_to_u_matrix%num_elements
                        local_field(j,to)%var(i)%field(regrid_t_to_u_matrix%dst_index%field(k)) =       & ! sparse matrix multiplication
                            local_field(j,to)%var(i)%field(regrid_t_to_u_matrix%dst_index%field(k)) +   & 
                            local_field(j,from)%var(i)%field(regrid_t_to_u_matrix%src_index%field(k)) * & 
                            regrid_t_to_u_matrix%weight%field(k)
                    ENDDO
                ENDIF
                IF (local_field(j,1)%var(i)%put_to_v_grid) THEN
                    from=1; to=3
                    local_field(j,to)%var(i)%field = 0.0
                    DO k=1,regrid_t_to_v_matrix%num_elements
                        local_field(j,to)%var(i)%field(regrid_t_to_v_matrix%dst_index%field(k)) =       & ! sparse matrix multiplication
                            local_field(j,to)%var(i)%field(regrid_t_to_v_matrix%dst_index%field(k)) +   & 
                            local_field(j,from)%var(i)%field(regrid_t_to_v_matrix%src_index%field(k)) * & 
                            regrid_t_to_v_matrix%weight%field(k)
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
            
    END SUBROUTINE do_regridding

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    SUBROUTINE init_varname_idx
        INTEGER :: i
        DO i=1,MAX_VARNAMES
            ! INPUT VARIABLES:
            IF (varnames(i) == 'ALBE')    idx_ALBE = i    ! Albedo of the surface                                    (1)
            IF (varnames(i) == 'ALBA')    idx_ALBA = i    ! averaged Albedo from atmosphere                          (1)
            IF (varnames(i) == 'AMOI')    idx_AMOI = i    ! Diffusion coefficient for moisture                       (1)
            IF (varnames(i) == 'AMOM')    idx_AMOM = i    ! Diffusion coefficient for momentum                       (1)
            IF (varnames(i) == 'FARE')    idx_FARE = i    ! Fraction of grid cell area covered by this surface_type  (1)
            IF (varnames(i) == 'FICE')    idx_FICE = i    ! Fraction of area covered by ice                          (1)
            IF (varnames(i) == 'PATM')    idx_PATM = i    ! Pressure in lowest atmospheric grid cell                 (Pa)
            IF (varnames(i) == 'PSUR')    idx_PSUR = i    ! Pressure at surface                                      (Pa)
            IF (varnames(i) == 'QATM')    idx_QATM = i    ! Moisture content of lowest atmospheric grid cell         (kg/kg)
            IF (varnames(i) == 'TATM')    idx_TATM = i    ! Temperature in lowest atmospheric grid cell              (K)
            IF (varnames(i) == 'TSUR')    idx_TSUR = i    ! Temperature at the surface                               (K)
            IF (varnames(i) == 'UATM')    idx_UATM = i    ! Eastward velocity in lowest atmospheric grid cell        (m/s)
            IF (varnames(i) == 'VATM')    idx_VATM = i    ! Northward velocity in lowest atmospheric grid cell       (m/s)
            IF (varnames(i) == 'U10M')    idx_U10M = i    ! Eastward velocity 10m above the surface                  (m/s)
            IF (varnames(i) == 'V10M')    idx_V10M = i    ! Northward velocity 10m above the surface                 (m/s)
            ! AUXILIARY VARIABLES:
            IF (varnames(i) == 'QSUR')    idx_QSUR = i    ! Moisture content directly above the surface              (kg/kg)
            ! FLUXES: (all are positive upward, such that e.g. precipitation is always negative)
            IF (varnames(i) == 'HLAT')    idx_HLAT = i    ! Heat flux (latent)                                       (W/m2)
            IF (varnames(i) == 'HSEN')    idx_HSEN = i    ! Heat flux (sensible)                                     (W/m2)
            IF (varnames(i) == 'MEVA')    idx_MEVA = i    ! Mass flux (evaporation minus condensation)               (kg/m2/s)
            IF (varnames(i) == 'MPRE')    idx_MPRE = i    ! Mass flux (total precipitation, negative values)         (kg/m2/s)
            IF (varnames(i) == 'MRAI')    idx_MRAI = i    ! Mass flux (liquid precipitation "rain", neg. val.)       (kg/m2/s)
            IF (varnames(i) == 'MSNO')    idx_MSNO = i    ! Mass flux (solid precipitation "snow", neg. val.)        (kg/m2/s)
            IF (varnames(i) == 'RBBR')    idx_RBBR = i    ! Radiation flux (blackbody radiation of surface)          (W/m2)
            IF (varnames(i) == 'RLWD')    idx_RLWD = i    ! Radiation flux (longwave downward, negative values)      (W/m2)
            IF (varnames(i) == 'RLWU')    idx_RLWU = i    ! Radiation flux (longwave upward)                         (W/m2)
            IF (varnames(i) == 'RSID')    idx_RSID = i    ! Radiation flux (shortwave indirect downward, neg. val.)  (W/m2)
            IF (varnames(i) == 'RSIU')    idx_RSIU = i    ! Radiation flux (shortwave indirect upward)               (W/m2)
            IF (varnames(i) == 'RSIN')    idx_RSIN = i    ! Radiation flux (shortwave indirect net upward)           (W/m2)
            IF (varnames(i) == 'RSDD')    idx_RSDD = i    ! Radiation flux (shortwave directed downward, neg. val.)  (W/m2)
            IF (varnames(i) == 'RSDR')    idx_RSDR = i    ! Radiation flux redistributed on surface types (shortwave directed downward, neg. val.)  (W/m2)
            IF (varnames(i) == 'UMOM')    idx_UMOM = i    ! Upward flux of eastward momentum                         (N/m2)
            IF (varnames(i) == 'VMOM')    idx_VMOM = i    ! Upward flux of northward momentum                        (N/m2)
        ENDDO
    END SUBROUTINE init_varname_idx

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    SUBROUTINE nullify_localvars(local_field)
        TYPE(local_fields_type), INTENT(INOUT) :: local_field
        INTEGER                                :: i
        DO i=1,MAX_VARNAMES
            NULLIFY(local_field%var(i)%field) ! this ensures that ASSOCIATED(local_field%var(i)%field)==.FALSE.
            local_field%var(i)%allocated     = .FALSE.
            local_field%var(i)%put_to_t_grid = .FALSE.
            local_field%var(i)%put_to_u_grid = .FALSE.
            local_field%var(i)%put_to_v_grid = .FALSE.
        ENDDO
    END SUBROUTINE nullify_localvars
    
    
END MODULE flux_calculator_basic