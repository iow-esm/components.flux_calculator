MODULE flux_calculator_calculate
! This module contains functions that do the calculations using the flux library functions

    use flux_calculator_basic
    use flux_library

    IMPLICIT NONE

!!!!!!!!!! FUNCTIONS DEFINED IN THIS MODULE
    PUBLIC calc_spec_vapor_surface
    PUBLIC calc_flux_mass_evap
    PUBLIC calc_flux_radiation_blackbody
    PUBLIC average_across_surface_types

!!!!!!!!!! NOW EVERYTHING ELSE

    CONTAINS

    !!!!!!!!!! AUXILIARY VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE calc_spec_vapor_surface(my_bottom_model, num_surface_types, which_grid, methods, grid_size, local_field) 
        ! calculates specific water vapor content on t_grid
        INTEGER,                                  INTENT(IN)    :: my_bottom_model
        INTEGER,                                  INTENT(IN)    :: num_surface_types
        INTEGER,                                  INTENT(IN)    :: which_grid
        CHARACTER(len=20),       DIMENSION(:,:),  INTENT(IN)    :: methods
        INTEGER,                 DIMENSION(:),    INTENT(IN)    :: grid_size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field
        CHARACTER(len=4), PARAMETER :: myvarname = 'QSUR'
        CHARACTER(len=20)           :: method
        INTEGER                     :: i, j

        DO i=1,num_surface_types
            method = methods(my_bottom_model,i)
            IF (trim(method) /= 'none') THEN 
                IF (trim(method)=='CCLM') THEN
                    DO j=1,grid_size(which_grid) 
                        CALL spec_vapor_surface_cclm(local_field(i,which_grid)%var(idx_QSUR)%field(j), &
                                                     local_field(i,which_grid)%var(idx_FICE)%field(j), &
                                                     local_field(i,which_grid)%var(idx_PSUR)%field(j), &
                                                     local_field(i,which_grid)%var(idx_TSUR)%field(j))
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
    END SUBROUTINE calc_spec_vapor_surface

    !!!!!!!!!! MASS FLUXES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE calc_flux_mass_evap(my_bottom_model, num_surface_types, methods, grid_size, local_field) 
        ! calculates specific water vapor content on t_grid
        INTEGER,                                  INTENT(IN)    :: my_bottom_model
        INTEGER,                                  INTENT(IN)    :: num_surface_types
        CHARACTER(len=20),       DIMENSION(:,:),  INTENT(IN)    :: methods
        INTEGER,                 DIMENSION(:),    INTENT(IN)    :: grid_size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field
        CHARACTER(len=4), PARAMETER :: myvarname = 'MEVA'
        CHARACTER(len=20)           :: method
        INTEGER                     :: i, j

        DO i=1,num_surface_types
            method = methods(my_bottom_model,i)
            IF (trim(method) /= 'none') THEN 
                IF (trim(method)=='zero') THEN
                    local_field(i,1)%var(idx_MEVA)%field(:) = 0.0
                ELSEIF (trim(method)=='CCLM') THEN
                    DO j=1,grid_size(1) 
                        CALL flux_mass_evap_cclm(local_field(i,1)%var(idx_MEVA)%field(j), &
                                                 local_field(i,1)%var(idx_AMOI)%field(j), &
                                                 local_field(i,1)%var(idx_PSUR)%field(j), &
                                                 local_field(i,1)%var(idx_QATM)%field(j), &
                                                 local_field(i,1)%var(idx_QSUR)%field(j), &
                                                 local_field(i,1)%var(idx_TATM)%field(j), &
                                                 local_field(i,1)%var(idx_UATM)%field(j), &
                                                 local_field(i,1)%var(idx_VATM)%field(j))
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
        CALL average_across_surface_types(1,idx_MEVA,num_surface_types,grid_size,local_field)
    END SUBROUTINE calc_flux_mass_evap

    !!!!!!!!!! HEAT FLUXES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE calc_flux_heat_latent(my_bottom_model, num_surface_types, methods, grid_size, local_field) 
        ! calculates latent heat flux on t_grid
        INTEGER,                                  INTENT(IN)    :: my_bottom_model
        INTEGER,                                  INTENT(IN)    :: num_surface_types
        CHARACTER(len=20),       DIMENSION(:,:),  INTENT(IN)    :: methods
        INTEGER,                 DIMENSION(:),    INTENT(IN)    :: grid_size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field
        CHARACTER(len=4), PARAMETER :: myvarname = 'HLAT'
        CHARACTER(len=20)           :: method
        INTEGER                     :: i, j

        DO i=1,num_surface_types
            method = methods(my_bottom_model,i)
            IF (trim(method) /= 'none') THEN 
                IF (trim(method)=='zero') THEN
                    local_field(i,1)%var(idx_MEVA)%field(:) = 0.0
                ELSEIF (trim(method)=='water') THEN
                    DO j=1,grid_size(1) 
                        CALL flux_heat_latent_water(local_field(i,1)%var(idx_HLAT)%field(j), &
                                                    local_field(i,1)%var(idx_MEVA)%field(j))
                    ENDDO
                ELSEIF (trim(method)=='ice') THEN
                    DO j=1,grid_size(1) 
                        CALL flux_heat_latent_ice(local_field(i,1)%var(idx_HLAT)%field(j), &
                                                  local_field(i,1)%var(idx_MEVA)%field(j))
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
        CALL average_across_surface_types(1,idx_HLAT,num_surface_types,grid_size,local_field)
    END SUBROUTINE calc_flux_heat_latent

    SUBROUTINE calc_flux_heat_sensible(my_bottom_model, num_surface_types, methods, grid_size, local_field) 
        ! calculates sensible heat flux on t_grid
        INTEGER,                                  INTENT(IN)    :: my_bottom_model
        INTEGER,                                  INTENT(IN)    :: num_surface_types
        CHARACTER(len=20),       DIMENSION(:,:),  INTENT(IN)    :: methods
        INTEGER,                 DIMENSION(:),    INTENT(IN)    :: grid_size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field
        CHARACTER(len=4), PARAMETER :: myvarname = 'HSEN'
        CHARACTER(len=20)           :: method
        INTEGER                     :: i, j

        DO i=1,num_surface_types
            method = methods(my_bottom_model,i)
            IF (trim(method) /= 'none') THEN 
                IF (trim(method)=='zero') THEN
                    local_field(i,1)%var(idx_MEVA)%field(:) = 0.0
                ELSEIF (trim(method)=='CCLM') THEN
                    DO j=1,grid_size(1) 
                        CALL flux_heat_sensible_cclm(local_field(i,1)%var(idx_HSEN)%field(j), &
                                                     local_field(i,1)%var(idx_AMOI)%field(j), &
                                                     local_field(i,1)%var(idx_PATM)%field(j), &
                                                     local_field(i,1)%var(idx_PSUR)%field(j), &
                                                     local_field(i,1)%var(idx_QATM)%field(j), &
                                                     local_field(i,1)%var(idx_TATM)%field(j), &
                                                     local_field(i,1)%var(idx_TSUR)%field(j), &
                                                     local_field(i,1)%var(idx_UATM)%field(j), &
                                                     local_field(i,1)%var(idx_VATM)%field(j))
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
        CALL average_across_surface_types(1,idx_HSEN,num_surface_types,grid_size,local_field)
    END SUBROUTINE calc_flux_heat_sensible

    !!!!!!!!!! MOMENTUM FLUXES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE calc_flux_momentum_east(my_bottom_model, num_surface_types, which_grid, methods, grid_size, local_field) 
        ! calculates sensible heat flux on t_grid
        INTEGER,                                  INTENT(IN)    :: my_bottom_model
        INTEGER,                                  INTENT(IN)    :: num_surface_types
        INTEGER,                                  INTENT(IN)    :: which_grid
        CHARACTER(len=20),       DIMENSION(:,:),  INTENT(IN)    :: methods
        INTEGER,                 DIMENSION(:),    INTENT(IN)    :: grid_size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field
        CHARACTER(len=4), PARAMETER :: myvarname = 'UMOM'
        CHARACTER(len=20)           :: method
        INTEGER                     :: i, j
        REAL                        :: dummy ! since northward momentum may be required on another grid

        DO i=1,num_surface_types
            method = methods(my_bottom_model,i)
            IF (trim(method) /= 'none') THEN 
                IF (trim(method)=='zero') THEN
                    local_field(i,which_grid)%var(idx_UMOM)%field(:) = 0.0
                ELSEIF (trim(method)=='CCLM') THEN
                    DO j=1,grid_size(which_grid) 
                        CALL flux_momentum_cclm(local_field(i,which_grid)%var(idx_UMOM)%field(j), &
                                                dummy,                                            &
                                                local_field(i,which_grid)%var(idx_AMOM)%field(j), &
                                                local_field(i,which_grid)%var(idx_PSUR)%field(j), &
                                                local_field(i,which_grid)%var(idx_QSUR)%field(j), &
                                                local_field(i,which_grid)%var(idx_TSUR)%field(j), &
                                                local_field(i,which_grid)%var(idx_UATM)%field(j), &
                                                local_field(i,which_grid)%var(idx_VATM)%field(j))
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
        CALL average_across_surface_types(1,idx_UMOM,num_surface_types,grid_size,local_field)
    END SUBROUTINE calc_flux_momentum_east

    SUBROUTINE calc_flux_momentum_north(my_bottom_model, num_surface_types, which_grid, methods, grid_size, local_field) 
        ! calculates sensible heat flux on t_grid
        INTEGER,                                  INTENT(IN)    :: my_bottom_model
        INTEGER,                                  INTENT(IN)    :: num_surface_types
        INTEGER,                                  INTENT(IN)    :: which_grid
        CHARACTER(len=20),       DIMENSION(:,:),  INTENT(IN)    :: methods
        INTEGER,                 DIMENSION(:),    INTENT(IN)    :: grid_size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field
        CHARACTER(len=4), PARAMETER :: myvarname = 'VMOM'
        CHARACTER(len=20)           :: method
        INTEGER                     :: i, j
        REAL                        :: dummy ! since northward momentum may be required on another grid

        DO i=1,num_surface_types
            method = methods(my_bottom_model,i)
            IF (trim(method) /= 'none') THEN 
                IF (trim(method)=='zero') THEN
                    local_field(i,which_grid)%var(idx_VMOM)%field(:) = 0.0
                ELSEIF (trim(method)=='CCLM') THEN
                    DO j=1,grid_size(which_grid) 
                        CALL flux_momentum_cclm(dummy,                                            &
                                                local_field(i,which_grid)%var(idx_VMOM)%field(j), &
                                                local_field(i,which_grid)%var(idx_AMOM)%field(j), &
                                                local_field(i,which_grid)%var(idx_PSUR)%field(j), &
                                                local_field(i,which_grid)%var(idx_QSUR)%field(j), &
                                                local_field(i,which_grid)%var(idx_TSUR)%field(j), &
                                                local_field(i,which_grid)%var(idx_UATM)%field(j), &
                                                local_field(i,which_grid)%var(idx_VATM)%field(j))
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
        CALL average_across_surface_types(1,idx_VMOM,num_surface_types,grid_size,local_field)
    END SUBROUTINE calc_flux_momentum_north

    !!!!!!!!!! RADIATION FLUXES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE calc_flux_radiation_blackbody(my_bottom_model, num_surface_types, methods, grid_size, local_field) 
        ! calculates blackbody radiation on t_grid
        INTEGER,                                  INTENT(IN)    :: my_bottom_model
        INTEGER,                                  INTENT(IN)    :: num_surface_types
        CHARACTER(len=20),       DIMENSION(:,:),  INTENT(IN)    :: methods
        INTEGER,                 DIMENSION(:),    INTENT(IN)    :: grid_size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field
        CHARACTER(len=4), PARAMETER :: myvarname = 'RBBR'
        CHARACTER(len=20)           :: method
        INTEGER                     :: i, j

        DO i=1,num_surface_types
            method = methods(my_bottom_model,i)
            IF (trim(method) /= 'none') THEN 
                IF (trim(method)=='zero') THEN
                    local_field(i,1)%var(idx_RBBR)%field(:) = 0.0
                ELSEIF (trim(method)=='StBo') THEN
                    DO j=1,grid_size(1) 
                        CALL flux_radiation_blackbody_StBo(local_field(i,1)%var(idx_RBBR)%field(j), &
                                                           local_field(i,1)%var(idx_TSUR)%field(j) )
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
        CALL average_across_surface_types(1,idx_RBBR,num_surface_types,grid_size,local_field)
    END SUBROUTINE calc_flux_radiation_blackbody

    SUBROUTINE distribute_shortwave_radiation_flux(my_bottom_model, num_surface_types, grid_size, local_field) 
        ! calculates blackbody radiation on t_grid
        INTEGER,                                  INTENT(IN)    :: my_bottom_model
        INTEGER,                                  INTENT(IN)    :: num_surface_types
        INTEGER,                 DIMENSION(:),    INTENT(IN)    :: grid_size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field
        INTEGER                     :: i, j

        DO i=1,num_surface_types
            local_field(i,1)%var(idx_RSDD)%field(:) = 0.0
            DO j=1,grid_size(1) 
                CALL distribute_radiation_flux(local_field(i,1)%var(idx_RSDD)%field(j), &
                                                    local_field(0,1)%var(idx_RSDD)%field(j), &
                                                    local_field(0,1)%var(idx_ALBE)%field(j), &
                                                    local_field(i,1)%var(idx_ALBE)%field(j), &
                                                    local_field(i,1)%var(idx_FARE)%field(j))
            ENDDO
        ENDDO
        !CALL average_across_surface_types(1,idx_RSDD,num_surface_types,grid_size,local_field)
    END SUBROUTINE distribute_shortwave_radiation_flux

    !!!!!!!!!! AVERAGING ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE average_across_surface_types(which_grid, my_idx, num_surface_types, grid_size, local_field)
        INTEGER,                                  INTENT(IN)    :: which_grid
        INTEGER,                                  INTENT(IN)    :: my_idx
        INTEGER,                                  INTENT(IN)    :: num_surface_types
        INTEGER,                 DIMENSION(:),    INTENT(IN)    :: grid_size
        TYPE(local_fields_type), DIMENSION(0:,:), INTENT(INOUT) :: local_field
        INTEGER                                                 :: i,j
        
        IF (local_field(0,1)%var(my_idx)%allocated) THEN
            local_field(0,1)%var(my_idx)%field=0.0
            DO i=1,num_surface_types
                DO j=1,grid_size(which_grid)
                    local_field(0,1)%var(my_idx)%field(j) = local_field(0,1)%var(my_idx)%field(j) + &
                                                            local_field(i,1)%var(my_idx)%field(j)*local_field(i,1)%var(idx_FARE)%field(j)
                ENDDO
            ENDDO
        ENDIF
    END SUBROUTINE average_across_surface_types
END MODULE flux_calculator_calculate