!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Momentum exchange fluxes                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reproduces formulae used in CCLM                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hagen Radtke, 2020                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module flux_momentum

  use flux_constants, only: default_values_type, default_values, prec

  implicit none ; private

  ! expose all functions
  public flux_momentum_cclm

contains

  subroutine flux_momentum_cclm( &
        flux_momentum_east,                     & ! RESULT       (N/m2)
        flux_momentum_north,                    & ! RESULT       (N/m2)
        diffusion_coefficient_momentum,         & ! a_{momentum} (1)
        pressure_surface,                       & ! p_s          (Pa)
        specific_vapor_content_surface,         & ! q_{v,s}      (kg/kg)
        temperature_surface,                    & ! T_s          (K)
        u_atmos,                                & ! u_a          (m/s)
        v_atmos,                                & ! v_a          (m/s)
        gas_constant_air_new,                   & ! R_d          (J/kg/K)
        gas_constant_vapor_new                  & ! R_v          (J/kg/K)
    )

    real(prec), intent(out) :: flux_momentum_east               ! RESULT       (N/m2)
    real(prec), intent(out) :: flux_momentum_north              ! RESULT       (N/m2)
    real(prec), intent(in)  :: diffusion_coefficient_momentum   ! a_{momentum} (1)
    real(prec), intent(in)  :: pressure_surface                 ! p_s          (Pa)
    real(prec), intent(in)  :: specific_vapor_content_surface   ! q_{v,s}      (kg/kg)
    real(prec), intent(in)  :: temperature_surface              ! T_s          (K)
    real(prec), intent(in)  :: u_atmos                          ! u_a          (m/s)
    real(prec), intent(in)  :: v_atmos                          ! v_a          (m/s)
    real(prec), intent(in), optional :: gas_constant_air_new    ! R_d          (J/kg/K)
    real(prec), intent(in), optional :: gas_constant_vapor_new  ! R_v          (J/kg/K)

    real(prec) :: gas_constant_air
    real(prec) :: gas_constant_vapor
    real(prec) :: T_tilde
    real(prec) :: vel
    real(prec) :: flux_air

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

    flux_air = diffusion_coefficient_momentum * vel *           &  ! mass exchange rate of air (kg/m2/s)
                 pressure_surface / (gas_constant_air * T_tilde)
                     
    flux_momentum_east  = - flux_air * u_atmos                      ! upward flux of eastward momentum (N/m2)
    flux_momentum_north = - flux_air * v_atmos                      ! upward flux of northward momentum (N/m2)

  end subroutine flux_momentum_cclm

end module flux_momentum