!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculating auxiliary variable: specific vapor content !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reproduces formulae used in CCLM                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hagen Radtke, 2020                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module aux_specific_vapor_content

  use flux_constants, only: default_values_type, default_values, prec

  implicit none ; private

  ! expose all functions
  public spec_vapor_surface_cclm

contains

  subroutine spec_vapor_surface_cclm( &
        specific_vapor_content_surface,         & ! RESULT       (kg/kg)
        fraction_ice,                           & ! f_{ice}      (must be 0.0 or 1.0)
        pressure_surface,                       & ! p_s          (Pa)
        temperature_surface,                    & ! T_s          (K)
        gas_constant_air_new,                   & ! R_d          (J/kg/K)
        gas_constant_vapor_new                  & ! R_v          (J/kg/K)
    )

    real(prec), intent(out) :: specific_vapor_content_surface   ! RESULT       (kg/kg)
    real(prec), intent(in)  :: fraction_ice                     ! f_{ice}      (must be 0.0 or 1.0) 
    real(prec), intent(in)  :: pressure_surface                 ! p_s          (Pa)
    real(prec), intent(in)  :: temperature_surface              ! T_s          (K)
    real(prec), intent(in), optional :: gas_constant_air_new    ! R_d          (J/kg/K)
    real(prec), intent(in), optional :: gas_constant_vapor_new  ! R_v          (J/kg/K)

    real(prec)            :: gas_constant_air          ! (J/kg/K)
    real(prec)            :: gas_constant_vapor        ! (J/kg/K)
    real(prec)            :: alpha                     ! (1)
    real(prec), parameter :: alpha_water = 17.2693882  ! (1)
    real(prec), parameter :: alpha_ice   = 21.8745584  ! (1)
    real(prec), parameter :: T_1 = 273.16              ! (K)
    real(prec)            :: T_2                       ! (K)
    real(prec), parameter :: T_2_water = 35.86         ! (K)
    real(prec), parameter :: T_2_ice   =  7.66         ! (K)
    real(prec), parameter :: p_0 = 610.78              ! (Pa)
    real(prec)            :: saturation_vapor_pressure ! (Pa)

    if (PRESENT(gas_constant_air_new)) then 
      gas_constant_air = gas_constant_air_new 
    else 
      gas_constant_air = default_values%gas_constant_air
    end if
    if (PRESENT(gas_constant_vapor_new)) then 
      gas_constant_vapor = gas_constant_vapor_new 
    else 
      gas_constant_vapor = default_values%gas_constant_vapor
    end if

    ! determine alpha and T_2 by multiplication to avoid if-clause for more efficient vectorization
    alpha = alpha_water + (alpha_ice-alpha_water) * fraction_ice
    T_2   = T_2_water   + (T_2_ice  -T_2_water  ) * fraction_ice
    
    saturation_vapor_pressure =                        & ! saturation vapor pressure (Pa)
        p_0 * exp(alpha * (temperature_surface - T_1) / (temperature_surface - T_2) )    
    
    specific_vapor_content_surface =                   & ! RESULT (kg/kg)
        (gas_constant_air/gas_constant_vapor) *saturation_vapor_pressure / &
        (pressure_surface - (1.0 - gas_constant_air/gas_constant_vapor)*saturation_vapor_pressure)

  end subroutine spec_vapor_surface_cclm

end module aux_specific_vapor_content