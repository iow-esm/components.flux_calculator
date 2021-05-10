!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mass exchange fluxes by evaporation/condensation!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reproduces formulae used in CCLM                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hagen Radtke, 2020                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module flux_mass_evap

  use flux_constants, only: default_values_type, default_values, prec

  implicit none ; private

  ! expose all functions
  public flux_mass_evap_cclm

contains

  subroutine flux_mass_evap_cclm( &
        flux_mass_evap,                         & ! RESULT       (kg/m2/s)
        diffusion_coefficient_moisture,         & ! a_{moisture} (1)
        pressure_surface,                       & ! p_s          (Pa)
        specific_vapor_content_atmos,           & ! q_{v,a}      (kg/kg)
        specific_vapor_content_surface,         & ! q_{v,s}      (kg/kg)
        temperature_surface,                    & ! T_s          (K)
        u_atmos,                                & ! u_a          (m/s)
        v_atmos,                                & ! v_a          (m/s)
        u_min_evap_new,                         & ! u_{min,evap} (m/s)
        gas_constant_air_new,                   & ! R_d          (J/kg/K)
        gas_constant_vapor_new                  & ! R_v          (J/kg/K)
    )

    real(prec), intent(out) :: flux_mass_evap                   ! RESULT       (kg/m2/s)
    real(prec), intent(in)  :: diffusion_coefficient_moisture   ! a_{moisture} (1)
    real(prec), intent(in)  :: pressure_surface                 ! p_s          (Pa)
    real(prec), intent(in)  :: specific_vapor_content_atmos     ! q_{v,a}      (kg/kg)
    real(prec), intent(in)  :: specific_vapor_content_surface   ! q_{v,s}      (kg/kg)
    real(prec), intent(in)  :: temperature_surface              ! T_s          (K)
    real(prec), intent(in)  :: u_atmos                          ! u_a          (m/s)
    real(prec), intent(in)  :: v_atmos                          ! v_a          (m/s)
    real(prec), intent(in), optional :: u_min_evap_new          ! u_{min,evap} (m/s)     
    real(prec), intent(in), optional :: gas_constant_air_new    ! R_d          (J/kg/K)
    real(prec), intent(in), optional :: gas_constant_vapor_new  ! R_v          (J/kg/K)

    real(prec) :: u_min_evap
    real(prec) :: gas_constant_air
    real(prec) :: gas_constant_vapor
    real(prec) :: T_tilde
    real(prec) :: vel
    real(prec) :: flux_air

    if (PRESENT(u_min_evap_new)) then 
      u_min_evap = u_min_evap_new 
    else 
      u_min_evap = default_values%u_min_evap
    end if
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
    

    T_tilde = temperature_surface * ( 1.0 +                    &  ! Temperature that dry air would need to have
                (gas_constant_vapor/gas_constant_air - 1.0) *  &  ! to show the same p*V as the moist air at the surface (K)
                specific_vapor_content_surface ) 

    vel = sqrt(u_atmos*u_atmos + v_atmos*v_atmos)                  ! atmospheric velocity (m/s)

    flux_air = diffusion_coefficient_moisture *                   &  ! mass exchange rate of air (kg/m2/s)
                 max(vel, u_min_evap) * pressure_surface /        &
                 (gas_constant_air * T_tilde)
        
    flux_mass_evap = flux_air * (specific_vapor_content_surface - &  ! mass flux of water (kg/m2/s)
                       specific_vapor_content_atmos)   

  end subroutine flux_mass_evap_cclm

end module flux_mass_evap