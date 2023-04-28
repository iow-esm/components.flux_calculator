MODULE bias_corrections

USE netcdf, ONLY : &
  nf90_open,               &
  nf90_close,              &
  nf90_get_var, &
  nf90_inq_varid, &
  nf90_get_att, &
  NF90_NOWRITE, &
  NF90_NOERR
  
IMPLICIT NONE

PUBLIC initialize_bias_corrections

ENUM, BIND(C)
    ENUMERATOR :: E_MASS_EVAP_CORRECTION = 1
    ENUMERATOR :: E_N_CORRECTIONS = 1
ENDENUM

! initializes names for corrextion -> corresponds to out variable that is corrected
CHARACTER(len=10), PARAMETER, DIMENSION(E_N_CORRECTIONS) :: corrections_names = [ &
    'mass_evap' &
    ] 

INTEGER :: &
    init_date

REAL, ALLOCATABLE :: &
    corrections(:,:,:) ! (variable, month, space)

LOGICAL  ::  &
    lcorrections(E_N_CORRECTIONS) = .FALSE.

CONTAINS

SUBROUTINE process_input_corrections (n, errstat)
    
! Parameter list:
  INTEGER, INTENT (IN)      ::        &
    n           ! Unit number for Namelist INPUT file

  INTEGER, INTENT (OUT)   ::        &
    errstat        ! error status variable

    
! Local variables: 
  INTEGER   ::  &
    i

! Define the namelist group
  NAMELIST /correctionsctl/ init_date, lcorrections
                      
!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_oasisctl
!------------------------------------------------------------------------------

  errstat = 0

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------
  READ (n, correctionsctl, IOSTAT=errstat)

  IF (errstat /= 0) THEN

    errstat = 1

    DO i = 1, E_N_CORRECTIONS
      lcorrections(i) = .FALSE.
    ENDDO   

    WRITE(*,*) "Could not read INPUT_BIAS set lcorrections to default: ", lcorrections 
    RETURN
  ENDIF
    
    
    !------------------------------------------------------------------------------
    !- Section 4: Check values for errors and consistency
    !------------------------------------------------------------------------------
    
  DO i = 1, E_N_CORRECTIONS
      IF ( lcorrections(i) /= .FALSE. .and.  lcorrections(i) /= .TRUE.) THEN
          WRITE (*,*) ' ERROR  *** Wrong value ', lcorrections(i) , ' for correction ', i
          errstat = 2
          WRITE (*,*) ' ERROR  *** Error while checking values of the namelist oasisctl *** '
          RETURN
      ENDIF
  ENDDO
    
    
  !------------------------------------------------------------------------------
  !- Section 6: Output of the namelist variables and their default values
  !------------------------------------------------------------------------------
    
  WRITE (*,*) 'Apply bias corrections for:'

  DO i = 1, E_N_CORRECTIONS
      WRITE (*,*)  corrections_names(i), lcorrections(i)
  ENDDO
  !------------------------------------------------------------------------------
  !- End of the Subroutine
  !------------------------------------------------------------------------------
    
END SUBROUTINE process_input_corrections 

SUBROUTINE read_namelist_corrections(filename, errstat)
    
  CHARACTER (LEN=  *),      INTENT(IN)           ::                      &
    filename      ! error message

  INTEGER, INTENT(OUT) :: errstat
  
  INTEGER :: n
  
  n = 10
  ! -----------------------------------------------------------------
  ! 1 Open NAMELIST-INPUT file
  ! ----------------------------------------------------------------

  errstat = 0

  WRITE  (*,*) '    INPUT OF THE NAMELIST FOR BIAS CORRECTIONS'

  OPEN (n, FILE=filename, FORM='FORMATTED', STATUS='UNKNOWN', IOSTAT=errstat)
  IF (errstat /= 0) THEN
    errstat = 1
    WRITE(*,*) ' ERROR    *** Error while opening file '//filename//' *** '
    RETURN
  ENDIF
        
    ! -----------------------------------------------------------------
    ! 2 read the NAMELIST group oasisctl
    ! ----------------------------------------------------------------
    
  CALL process_input_corrections (n, errstat)
    
  IF (errstat /= 0) THEN
    WRITE (*,*) ' ERROR *** Wrong values occured in NAMELIST group /correctionsctl/ *** '
    errstat = 2
    RETURN
  ENDIF
    
    ! -----------------------------------------------------------------
    ! 3 Close NAMELIST-INPUT file
    ! ----------------------------------------------------------------
  
    
  CLOSE (n, STATUS='KEEP', IOSTAT=errstat)
  IF (errstat /= 0) THEN
    WRITE(*,*) ' ERROR *** while closing file '//filename//'*** '
    errstat  = 4
  ENDIF

    !------------------------------------------------------------------------------
    !- End of the Subroutine
    !------------------------------------------------------------------------------
    
END SUBROUTINE read_namelist_corrections

SUBROUTINE initialize_bias_corrections(filename, grid_offset, grid_size)
  
  CHARACTER (LEN=  *),      INTENT(IN)           ::                      &
    filename      ! name of input file

  INTEGER, INTENT(IN) :: &
    grid_offset, &
    grid_size

  INTEGER :: &
    istatus,                 & ! NetCDF status
    ncfileid, ncvarid,       & ! NetCDF IDs
    nerror, &
    i, j

  REAL:: &
    fillvalue    
    
  CHARACTER(LEN=128) :: correction_filename

  CHARACTER(LEN=2) :: yj

  ! Read namelist for corrections
  CALL read_namelist_corrections(filename, nerror)

  ! allocate resources for each process
  ALLOCATE(corrections(E_N_CORRECTIONS, 12, grid_size))

  ! intialize with zero
  corrections(:,:,:) = 0.0

  ! Read in correction for each month
  DO i = 1, E_N_CORRECTIONS
      IF (.NOT. lcorrections(i)) THEN
          CYCLE
      ENDIF

      DO j = 1, 12
          WRITE (yj,'(I2.2)') j 
          yj = ADJUSTL(yj)

          correction_filename = 'corrections/'//TRIM(corrections_names(i))//'-'//TRIM(yj)//'.nc'

          istatus = nf90_open(TRIM(correction_filename), NF90_NOWRITE, ncfileid)
          IF (istatus /= NF90_NOERR) THEN
              WRITE(*,*) 'Could not open ', TRIM(correction_filename), ' for bias correction. Unset correction.'
              CYCLE
          ENDIF

          istatus = nf90_inq_varid(ncfileid, TRIM(corrections_names(i)) , ncvarid)
          IF (istatus /= NF90_NOERR) THEN
              WRITE(*,*) 'Could not get varid for variable '//TRIM(corrections_names(i))//'. Unset correction.'
              CYCLE
          ENDIF

          istatus = nf90_get_var(ncfileid, ncvarid, corrections(i,j,:), &
                              (/ grid_offset/),      &
                              (/ grid_size/))
          IF (istatus /= NF90_NOERR) THEN
              WRITE(*,*) 'Could not get variable '//TRIM(corrections_names(i))//'. Unset correction.'
              corrections(i,j,:) = 0.0
              CYCLE
          ENDIF    

          istatus = nf90_get_att(ncfileid, ncvarid, "_FillValue", fillvalue)
          IF (istatus /= NF90_NOERR) THEN
              WRITE(*,*) 'Could not get fill value. Unset correction.'
              corrections(i,j,:) = 0.0
              CYCLE
          ENDIF             

          istatus = nf90_close(ncfileid)
          IF (istatus /= NF90_NOERR) THEN
              WRITE(*,*) 'Could not close ', TRIM(correction_filename), 'for bias correction.'
              CYCLE
          ENDIF

          WHERE (corrections(i,j,:) == fillvalue) corrections(i,j,:) = 0.0
      ENDDO

  ENDDO

  WRITE (*,*) "Read in bias correction fields."

END SUBROUTINE initialize_bias_corrections

SUBROUTINE finalize_bias_corrections

    DEALLOCATE(corrections)

END SUBROUTINE finalize_bias_corrections

END MODULE bias_corrections