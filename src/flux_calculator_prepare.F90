MODULE flux_calculator_prepare
! This module contains functions that 
!    (1) check whether the calculations are possible, i.e. all required variables will be known,
!    (2) allocate the arrays to store the output results

    use flux_calculator_basic

    IMPLICIT NONE

!!!!!!!!!! FUNCTIONS DEFINED IN THIS MODULE
    PUBLIC prepare_spec_vapor_surface
    PUBLIC prepare_flux_mass_evap
    PUBLIC prepare_flux_radiation_blackbody

!!!!!!!!!! NOW EVERYTHING ELSE

    CONTAINS

    SUBROUTINE do_prepare_calculation(missing_field, idx, myvarname, surface_type, which_grid, method, grid_size, local_field)
        ! This subroutine will be called after the check which input fields are present
        CHARACTER(len=200),                       INTENT(IN)    :: missing_field
        INTEGER,                                  INTENT(IN)    :: idx
        CHARACTER(len=4),                         INTENT(IN)    :: myvarname
        INTEGER,                                  INTENT(IN)    :: surface_type         
        INTEGER,                                  INTENT(IN)    :: which_grid           ! 1=t_grid, 2=u_grid, 3=v_grid
        CHARACTER(len=20),                        INTENT(IN)    :: method              
        INTEGER,                                  INTENT(IN)    :: grid_size            
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field          ! pass the entire field because method can be "copy"
        IF (TRIM(missing_field) /= '') THEN 
            ! something went wrong - print the error message and exit
            WRITE (w_unit,*) "Error calculating ",myvarname," for surface_type ",surface_type," on the grid ",grid_name(which_grid),":"
            WRITE (w_unit,*) "    For method ",method," we are lacking the following variables: ",trim(missing_field)
            CALL mpi_finalize(1)
        ELSE
            ! success - all variables present
            IF (trim(method)=='copy') THEN
                ! no need to allocate - we can just set a pointer
                local_field(surface_type,which_grid)%var(idx)%field => local_field(1,which_grid)%var(idx)%field
            ELSE
                ALLOCATE(local_field(surface_type,which_grid)%var(idx)%field(grid_size))
                local_field(surface_type,which_grid)%var(idx)%allocated = .TRUE.
            ENDIF
        ENDIF
    END SUBROUTINE do_prepare_calculation

    !!!!!!!!!! AUXILIARY VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE prepare_spec_vapor_surface(surface_type, which_grid, method, grid_size, local_field) 
        INTEGER,                                  INTENT(IN)    :: surface_type         
        INTEGER,                                  INTENT(IN)    :: which_grid           ! 1=t_grid, 2=u_grid, 3=v_grid
        CHARACTER(len=20),                        INTENT(IN)    :: method              
        INTEGER,                                  INTENT(IN)    :: grid_size            ! to allocate the array to the correct size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field          ! pass the entire field because method can be "copy"
        CHARACTER(len=4), PARAMETER :: myvarname = 'QSUR'
        CHARACTER(len=200)          :: missing_field
        missing_field=''
        IF (trim(method) /= 'none') THEN 
            IF (trim(method)=='copy') THEN ! copy from surface_type=1
                IF (.NOT. ASSOCIATED( local_field(1,which_grid)%var(idx_QSUR)%field )) missing_field=myvarname//' for surface_type=1 '
            ELSEIF (trim(method)=='CCLM') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_FICE)%field )) missing_field=trim(missing_field)//' FICE'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_PSUR)%field )) missing_field=trim(missing_field)//' PSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TSUR)%field )) missing_field=trim(missing_field)//' TSUR'
            ELSE
                WRITE (w_unit,*) "Error calculating ",myvarname," for surface_type ",surface_type," on the grid ",grid_name(which_grid),":"
                WRITE (w_unit,*) "    Method ",method," is not known. "
                CALL mpi_finalize(1)
            ENDIF
            CALL do_prepare_calculation(missing_field, idx_QSUR, myvarname, surface_type, which_grid, method, grid_size, local_field)
        ENDIF
    END SUBROUTINE prepare_spec_vapor_surface

    !!!!!!!!!! MASS FLUXES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE prepare_flux_mass_evap(surface_type, which_grid, method, grid_size, local_field) 
        INTEGER,                                  INTENT(IN)    :: surface_type         
        INTEGER,                                  INTENT(IN)    :: which_grid           ! 1=t_grid, 2=u_grid, 3=v_grid
        CHARACTER(len=20),                        INTENT(IN)    :: method              
        INTEGER,                                  INTENT(IN)    :: grid_size            ! to allocate the array to the correct size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field          ! pass the entire field because method can be "copy"
        CHARACTER(len=4), PARAMETER :: myvarname = 'MEVA'
        CHARACTER(len=200)          :: missing_field
        missing_field=''
        IF (trim(method) /= 'none') THEN 
            IF (trim(method)=='copy') THEN ! copy from surface_type=1
                IF (.NOT. ASSOCIATED( local_field(1,which_grid)%var(idx_MEVA)%field )) missing_field=myvarname//' for surface_type=1 '
            ELSEIF (trim(method)=='zero') THEN
                ! nothing to be done
            ELSEIF (trim(method)=='CCLM') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_AMOI)%field )) missing_field=trim(missing_field)//' AMOI'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_PSUR)%field )) missing_field=trim(missing_field)//' PSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_QATM)%field )) missing_field=trim(missing_field)//' QATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_QSUR)%field )) missing_field=trim(missing_field)//' QSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_UATM)%field )) missing_field=trim(missing_field)//' UATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_VATM)%field )) missing_field=trim(missing_field)//' VATM'
            ELSEIF (trim(method)=='MOM5') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_CMOI)%field )) missing_field=trim(missing_field)//' CMOI'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_PSUR)%field )) missing_field=trim(missing_field)//' PSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_QATM)%field )) missing_field=trim(missing_field)//' QATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_QSUR)%field )) missing_field=trim(missing_field)//' QSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_UATM)%field )) missing_field=trim(missing_field)//' UATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_VATM)%field )) missing_field=trim(missing_field)//' VATM'
            ELSEIF (trim(method)=='RCO') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_QATM)%field )) missing_field=trim(missing_field)//' QATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_QSUR)%field )) missing_field=trim(missing_field)//' TSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_UATM)%field )) missing_field=trim(missing_field)//' UATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_VATM)%field )) missing_field=trim(missing_field)//' VATM'       
            ELSE
                WRITE (w_unit,*) "Error calculating ",myvarname," for surface_type ",surface_type," on the grid ",grid_name(which_grid),":"
                WRITE (w_unit,*) "    Method ",method," is not known. "
                CALL mpi_finalize(1)
            ENDIF
            CALL do_prepare_calculation(missing_field, idx_MEVA, myvarname, surface_type, which_grid, method, grid_size, local_field)
        ENDIF
    END SUBROUTINE prepare_flux_mass_evap

    !!!!!!!!!! HEAT FLUXES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE prepare_flux_heat_latent(surface_type, which_grid, method, grid_size, local_field) 
        INTEGER,                                  INTENT(IN)    :: surface_type         
        INTEGER,                                  INTENT(IN)    :: which_grid           ! 1=t_grid, 2=u_grid, 3=v_grid
        CHARACTER(len=20),                        INTENT(IN)    :: method              
        INTEGER,                                  INTENT(IN)    :: grid_size            ! to allocate the array to the correct size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field          ! pass the entire field because method can be "copy"
        CHARACTER(len=4), PARAMETER :: myvarname = 'HLAT'
        CHARACTER(len=200)          :: missing_field
        missing_field=''
        IF (trim(method) /= 'none') THEN 
            IF (trim(method)=='copy') THEN ! copy from surface_type=1
                IF (.NOT. ASSOCIATED( local_field(1,which_grid)%var(idx_HSEN)%field )) missing_field=myvarname//' for surface_type=1 '
            ELSEIF (trim(method)=='zero') THEN
                ! nothing to be done
            ELSEIF (trim(method)=='water') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_MEVA)%field )) missing_field=trim(missing_field)//' MEVA'
            ELSEIF (trim(method)=='ice') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_MEVA)%field )) missing_field=trim(missing_field)//' MEVA'
            ELSE
                WRITE (w_unit,*) "Error calculating ",myvarname," for surface_type ",surface_type," on the grid ",grid_name(which_grid),":"
                WRITE (w_unit,*) "    Method ",method," is not known. "
                CALL mpi_finalize(1)
            ENDIF
            CALL do_prepare_calculation(missing_field, idx_HLAT, myvarname, surface_type, which_grid, method, grid_size, local_field)
        ENDIF
    END SUBROUTINE prepare_flux_heat_latent

    SUBROUTINE prepare_flux_heat_sensible(surface_type, which_grid, method, grid_size, local_field) 
        INTEGER,                                  INTENT(IN)    :: surface_type         
        INTEGER,                                  INTENT(IN)    :: which_grid           ! 1=t_grid, 2=u_grid, 3=v_grid
        CHARACTER(len=20),                        INTENT(IN)    :: method              
        INTEGER,                                  INTENT(IN)    :: grid_size            ! to allocate the array to the correct size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field          ! pass the entire field because method can be "copy"
        CHARACTER(len=4), PARAMETER :: myvarname = 'HSEN'
        CHARACTER(len=200)          :: missing_field
        missing_field=''
        IF (trim(method) /= 'none') THEN 
            IF (trim(method)=='copy') THEN ! copy from surface_type=1
                IF (.NOT. ASSOCIATED( local_field(1,which_grid)%var(idx_HSEN)%field )) missing_field=myvarname//' for surface_type=1 '
            ELSEIF (trim(method)=='zero') THEN
                ! nothing to be done
            ELSEIF (trim(method)=='CCLM') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_AMOI)%field )) missing_field=trim(missing_field)//' AMOI'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_PATM)%field )) missing_field=trim(missing_field)//' PATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_PSUR)%field )) missing_field=trim(missing_field)//' PSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_QSUR)%field )) missing_field=trim(missing_field)//' QSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_UATM)%field )) missing_field=trim(missing_field)//' UATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_VATM)%field )) missing_field=trim(missing_field)//' VATM'
            ELSEIF (trim(method)=='MOM5') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_CHEA)%field )) missing_field=trim(missing_field)//' CHEA'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_PATM)%field )) missing_field=trim(missing_field)//' PATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_PSUR)%field )) missing_field=trim(missing_field)//' PSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_QSUR)%field )) missing_field=trim(missing_field)//' QSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_UATM)%field )) missing_field=trim(missing_field)//' UATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_VATM)%field )) missing_field=trim(missing_field)//' VATM'
            ELSEIF (trim(method)=='RCO') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_UATM)%field )) missing_field=trim(missing_field)//' UATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_VATM)%field )) missing_field=trim(missing_field)//' VATM'
            ELSE
                WRITE (w_unit,*) "Error calculating ",myvarname," for surface_type ",surface_type," on the grid ",grid_name(which_grid),":"
                WRITE (w_unit,*) "    Method ",method," is not known. "
                CALL mpi_finalize(1)
            ENDIF
            CALL do_prepare_calculation(missing_field, idx_HSEN, myvarname, surface_type, which_grid, method, grid_size, local_field)
        ENDIF
    END SUBROUTINE prepare_flux_heat_sensible

    !!!!!!!!!! MOMENTUM FLUXES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE prepare_flux_momentum_east(surface_type, which_grid, method, grid_size, local_field) 
        INTEGER,                                  INTENT(IN)    :: surface_type         
        INTEGER,                                  INTENT(IN)    :: which_grid           ! 1=t_grid, 2=u_grid, 3=v_grid
        CHARACTER(len=20),                        INTENT(IN)    :: method              
        INTEGER,                                  INTENT(IN)    :: grid_size            ! to allocate the array to the correct size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field          ! pass the entire field because method can be "copy"
        CHARACTER(len=4), PARAMETER :: myvarname = 'UMOM'
        CHARACTER(len=200)          :: missing_field
        missing_field=''
        IF (trim(method) /= 'none') THEN 
            IF (trim(method)=='copy') THEN ! copy from surface_type=1
                IF (.NOT. ASSOCIATED( local_field(1,which_grid)%var(idx_UMOM)%field )) missing_field=myvarname//' for surface_type=1 '
            ELSEIF (trim(method)=='zero') THEN
                ! nothing to be done
            ELSEIF (trim(method)=='CCLM') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_AMOM)%field )) missing_field=trim(missing_field)//' AMOM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_PSUR)%field )) missing_field=trim(missing_field)//' PSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_QSUR)%field )) missing_field=trim(missing_field)//' QSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_UATM)%field )) missing_field=trim(missing_field)//' UATM'
            ELSEIF (trim(method)=='MOM5') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_CMOM)%field )) missing_field=trim(missing_field)//' CMOM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_PSUR)%field )) missing_field=trim(missing_field)//' PSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_QSUR)%field )) missing_field=trim(missing_field)//' QSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_UATM)%field )) missing_field=trim(missing_field)//' UATM'
            ELSEIF (trim(method)=='RCO') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_UATM)%field )) missing_field=trim(missing_field)//' UATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_VATM)%field )) missing_field=trim(missing_field)//' VATM'                
            ELSE
                WRITE (w_unit,*) "Error calculating ",myvarname," for surface_type ",surface_type," on the grid ",grid_name(which_grid),":"
                WRITE (w_unit,*) "    Method ",method," is not known. "
                CALL mpi_finalize(1)
            ENDIF
            CALL do_prepare_calculation(missing_field, idx_UMOM, myvarname, surface_type, which_grid, method, grid_size, local_field)
        ENDIF
    END SUBROUTINE prepare_flux_momentum_east

    SUBROUTINE prepare_flux_momentum_north(surface_type, which_grid, method, grid_size, local_field) 
        INTEGER,                                  INTENT(IN)    :: surface_type         
        INTEGER,                                  INTENT(IN)    :: which_grid           ! 1=t_grid, 2=u_grid, 3=v_grid
        CHARACTER(len=20),                        INTENT(IN)    :: method              
        INTEGER,                                  INTENT(IN)    :: grid_size            ! to allocate the array to the correct size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field          ! pass the entire field because method can be "copy"
        CHARACTER(len=4), PARAMETER :: myvarname = 'VMOM'
        CHARACTER(len=200)          :: missing_field
        missing_field=''
        IF (trim(method) /= 'none') THEN 
            IF (trim(method)=='copy') THEN ! copy from surface_type=1
                IF (.NOT. ASSOCIATED( local_field(1,which_grid)%var(idx_VMOM)%field )) missing_field=myvarname//' for surface_type=1 '
            ELSEIF (trim(method)=='zero') THEN
                ! nothing to be done
            ELSEIF (trim(method)=='CCLM') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_AMOM)%field )) missing_field=trim(missing_field)//' AMOM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_PSUR)%field )) missing_field=trim(missing_field)//' PSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_QSUR)%field )) missing_field=trim(missing_field)//' QSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_UATM)%field )) missing_field=trim(missing_field)//' VATM'
            ELSEIF (trim(method)=='MOM5') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_CMOM)%field )) missing_field=trim(missing_field)//' CMOM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_PSUR)%field )) missing_field=trim(missing_field)//' PSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_QSUR)%field )) missing_field=trim(missing_field)//' QSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TATM)%field )) missing_field=trim(missing_field)//' TSUR'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_UATM)%field )) missing_field=trim(missing_field)//' VATM'
            ELSEIF (trim(method)=='RCO') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_UATM)%field )) missing_field=trim(missing_field)//' UATM'
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_VATM)%field )) missing_field=trim(missing_field)//' VATM'                
            ELSE                
                WRITE (w_unit,*) "Error calculating ",myvarname," for surface_type ",surface_type," on the grid ",grid_name(which_grid),":"
                WRITE (w_unit,*) "    Method ",method," is not known. "
                CALL mpi_finalize(1)
            ENDIF
            CALL do_prepare_calculation(missing_field, idx_VMOM, myvarname, surface_type, which_grid, method, grid_size, local_field)
        ENDIF
    END SUBROUTINE prepare_flux_momentum_north

    !!!!!!!!!! RADIATION FLUXES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE prepare_flux_radiation_blackbody(surface_type, which_grid, method, grid_size, local_field) 
        INTEGER,                                  INTENT(IN)    :: surface_type         
        INTEGER,                                  INTENT(IN)    :: which_grid           ! 1=t_grid, 2=u_grid, 3=v_grid
        CHARACTER(len=20),                        INTENT(IN)    :: method              
        INTEGER,                                  INTENT(IN)    :: grid_size            ! to allocate the array to the correct size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field          ! pass the entire field because method can be "copy"
        CHARACTER(len=4), PARAMETER :: myvarname = 'RBBR'
        CHARACTER(len=200)          :: missing_field
        missing_field=''
        IF (trim(method) /= 'none') THEN 
            IF (trim(method)=='copy') THEN ! copy from surface_type=1
                IF (.NOT. ASSOCIATED( local_field(1,which_grid)%var(idx_RBBR)%field )) missing_field=myvarname//' for surface_type=1 '
            ELSEIF (trim(method)=='zero') THEN
                ! nothing to be done
            ELSEIF (trim(method)=='StBo') THEN 
                IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_TSUR)%field )) missing_field=trim(missing_field)//' TSUR'
            ELSE
                WRITE (w_unit,*) "Error calculating ",myvarname," for surface_type ",surface_type," on the grid ",grid_name(which_grid),":"
                WRITE (w_unit,*) "    Method ",method," is not known. "
                CALL mpi_finalize(1)
            ENDIF
            CALL do_prepare_calculation(missing_field, idx_RBBR, myvarname, surface_type, which_grid, method, grid_size, local_field)
        ENDIF
    END SUBROUTINE prepare_flux_radiation_blackbody

    ! SUBROUTINE prepare_distribute_radiation_flux(surface_type, which_grid, method, grid_size, local_field) 
    !     INTEGER,                                  INTENT(IN)    :: surface_type         
    !     INTEGER,                                  INTENT(IN)    :: which_grid           ! 1=t_grid, 2=u_grid, 3=v_grid
    !     CHARACTER(len=20),                        INTENT(IN)    :: method              
    !     INTEGER,                                  INTENT(IN)    :: grid_size            ! to allocate the array to the correct size
    !     TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field          ! pass the entire field because method can be "copy"
    !     CHARACTER(len=4), PARAMETER :: myvarname = 'RBBR'
    !     CHARACTER(len=200)          :: missing_field
    !     missing_field=''
    !     IF (trim(method) /= 'none') THEN 
    !         IF (trim(method)=='copy') THEN ! copy from surface_type=1
    !             IF (.NOT. ASSOCIATED( local_field(1,which_grid)%var(idx_RSDD)%field )) missing_field=myvarname//' for surface_type=1 '
    !         ELSEIF (trim(method)=='zero') THEN
    !             ! nothing to be done
    !         ELSEIF (trim(method)=='test') THEN 
    !             IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_RSDD)%field )) missing_field=trim(missing_field)//' RSDD'
    !             IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_ALBE)%field )) missing_field=trim(missing_field)//' ALBE'
    !             IF (.NOT. ASSOCIATED( local_field(surface_type,which_grid)%var(idx_FARE)%field )) missing_field=trim(missing_field)//' FARE'
    !         ELSE
    !             WRITE (w_unit,*) "Error calculating ",myvarname," for surface_type ",surface_type," on the grid ",grid_name(which_grid),":"
    !             WRITE (w_unit,*) "    Method ",method," is not known. "
    !             CALL mpi_finalize(1)
    !         ENDIF
    !         CALL do_prepare_calculation(missing_field, idx_RSDD, myvarname, surface_type, which_grid, method, grid_size, local_field)
    !     ENDIF
    ! END SUBROUTINE prepare_distribute_radiation_flux
END MODULE flux_calculator_prepare