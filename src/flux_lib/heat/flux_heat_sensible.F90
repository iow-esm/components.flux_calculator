!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flux library for different options of           !
! calculating exchange fluxes between atmosphere  !
! and ocean/ice/land.                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reproduces formulae used in CCLM                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hagen Radtke, 2020                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module flux_heat_sensible

  use flux_constants, only: default_values_type, default_values, prec

  implicit none ; private

  ! expose all functions
  public flux_heat_sensible_cclm
  public flux_heat_sensible_mom5
  public flux_heat_sensible_rco

contains

  subroutine flux_heat_sensible_cclm( &
        flux_heat_sensible,                     & ! RESULT       (W/m2)
        diffusion_coefficient_moisture,         & ! a_{moisture} (1)
        pressure_atmos,                         & ! p_a          (Pa)
        pressure_surface,                       & ! p_s          (Pa)
        specific_vapor_content_surface,         & ! q_{v,s}      (kg/kg)
        temperature_atmos,                      & ! T_a          (K)
        temperature_surface,                    & ! T_s          (K)
        u_atmos,                                & ! u_a          (m/s)
        v_atmos,                                & ! v_a          (m/s)
        heat_capacity_air_new,                  & ! c_d          (J/kg/K)
        u_min_evap_new,                         & ! u_{min,evap} (m/s)
        gas_constant_air_new,                   & ! R_d          (J/kg/K)
        gas_constant_vapor_new                  & ! R_v          (J/kg/K)
    )

    real(prec), intent(out) :: flux_heat_sensible               ! RESULT       (W/m2)
    real(prec), intent(in)  :: diffusion_coefficient_moisture   ! a_{moisture} (1)
    real(prec), intent(in)  :: pressure_atmos                   ! p_a          (Pa)
    real(prec), intent(in)  :: pressure_surface                 ! p_s          (Pa)
    real(prec), intent(in)  :: specific_vapor_content_surface   ! q_{v,s}      (kg/kg)
    real(prec), intent(in)  :: temperature_atmos                ! T_a          (K)
    real(prec), intent(in)  :: temperature_surface              ! T_s          (K)
    real(prec), intent(in)  :: u_atmos                          ! u_a          (m/s)
    real(prec), intent(in)  :: v_atmos                          ! v_a          (m/s)
    real(prec), intent(in), optional :: heat_capacity_air_new   ! c_d          (J/kg/K)     
    real(prec), intent(in), optional :: u_min_evap_new          ! u_{min,evap} (m/s)     
    real(prec), intent(in), optional :: gas_constant_air_new    ! R_d          (J/kg/K)
    real(prec), intent(in), optional :: gas_constant_vapor_new  ! R_v          (J/kg/K)

    real(prec) :: heat_capacity_air
    real(prec) :: u_min_evap
    real(prec) :: gas_constant_air
    real(prec) :: gas_constant_vapor
    real(prec) :: T_tilde
    real(prec) :: vel
    real(prec) :: flux_air
    real(prec) :: EF

    if (PRESENT(heat_capacity_air_new)) then 
        heat_capacity_air = heat_capacity_air_new 
    else 
        heat_capacity_air = default_values%heat_capacity_air
    end if
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

    vel = sqrt(u_atmos*u_atmos + v_atmos*v_atmos)                 ! atmospheric velocity (m/s)

    flux_air = diffusion_coefficient_moisture *                &  ! mass exchange rate of air (kg/m2/s)
                 max(vel, u_min_evap) * pressure_surface /     &
                 (gas_constant_air * T_tilde)

    EF = (pressure_surface/pressure_atmos) **                  &  ! Exner function (1)
           (gas_constant_air / heat_capacity_air)   
        
    flux_heat_sensible = flux_air * heat_capacity_air *        &  ! sensible heat flux (W/m2)
                           (temperature_surface - temperature_atmos * EF)   
  end subroutine flux_heat_sensible_cclm

  subroutine flux_heat_sensible_mom5( &
    flux_heat_sensible,                     & ! RESULT       (W/m2)
    diffusion_coefficient_moisture,         & ! a_{moisture} (1)
    pressure_atmos,                         & ! p_a          (Pa)
    pressure_surface,                       & ! p_s          (Pa)
    specific_vapor_content_surface,         & ! q_{v,s}      (kg/kg)
    temperature_atmos,                      & ! T_a          (K)
    temperature_surface,                    & ! T_s          (K)
    u_atmos,                                & ! u_a          (m/s)
    v_atmos                                 & ! v_a          (m/s)
)

    real(prec), intent(out) :: flux_heat_sensible               ! RESULT       (W/m2)
    real(prec), intent(in)  :: diffusion_coefficient_moisture   ! a_{moisture} (1)
    real(prec), intent(in)  :: pressure_atmos                   ! p_a          (Pa)
    real(prec), intent(in)  :: pressure_surface                 ! p_s          (Pa)
    real(prec), intent(in)  :: specific_vapor_content_surface   ! q_{v,s}      (kg/kg)
    real(prec), intent(in)  :: temperature_atmos                ! T_a          (K)
    real(prec), intent(in)  :: temperature_surface              ! T_s          (K)
    real(prec), intent(in)  :: u_atmos                          ! u_a          (m/s)
    real(prec), intent(in)  :: v_atmos                          ! v_a          (m/s)


    call flux_heat_sensible_cclm(flux_heat_sensible,                     & ! RESULT       (W/m2)
                                 diffusion_coefficient_moisture,         & ! a_{moisture} (1)
                                 pressure_atmos,                         & ! p_a          (Pa)
                                 pressure_surface,                       & ! p_s          (Pa)
                                 specific_vapor_content_surface,         & ! q_{v,s}      (kg/kg)
                                 temperature_atmos,                      & ! T_a          (K)
                                 temperature_surface,                    & ! T_s          (K)
                                 u_atmos,                                & ! u_a          (m/s)
                                 v_atmos                                 & ! v_a          (m/s)
    )
  
end subroutine flux_heat_sensible_mom5  

subroutine flux_heat_sensible_rco( &
  flux_heat_sensible,                     & ! RESULT       (W/m2)
  temperature_atmos,                      & ! T_a          (K)
  temperature_surface,                    & ! T_s          (K)
  u_atmos,                                & ! u_a          (m/s)
  v_atmos                                 & ! v_a          (m/s)
)

  real(prec), intent(out) :: flux_heat_sensible               ! RESULT       (W/m2)
  real(prec), intent(in)  :: temperature_atmos                ! T_a          (K)
  real(prec), intent(in)  :: temperature_surface              ! T_s          (K)
  real(prec), intent(in)  :: u_atmos                          ! u_a          (m/s)
  real(prec), intent(in)  :: v_atmos                          ! v_a          (m/s)

  real(prec) :: rho_a = 1.225 ! air density [kg / m^3] 
  real(prec) :: c_pa = 1.008E+03 ! specific heat capacity of air [J / (kg K)]
  real(prec) :: c_aw ! transfer coefficient for sensible heat (Stanton number) [1]
  real(prec) :: vel ! absolute value of wind speed

  ! get Stanton number according to temperature difference
  IF (temperature_atmos .lt. temperature_surface) THEN
    c_aw = 1.13E-03 ! unstable
  ELSE
    c_aw = 0.66E-03 ! stable
  ENDIF

  vel = sqrt(u_atmos*u_atmos + v_atmos*v_atmos)                  ! atmospheric velocity (m/s)

  flux_heat_sensible = rho_a * c_pa * c_aw * vel * (temperature_surface - temperature_atmos)

end subroutine flux_heat_sensible_rco
  
end module flux_heat_sensible


