!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Default values for constants of flux library    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hagen Radtke, 2020                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module flux_constants

  implicit none

  integer, parameter :: prec=selected_real_kind(15, 307)  ! double precision

  type default_values_type

    real(prec) :: heat_capacity_air         = 1005.0         ! c_d          heat capacity of dry air at constant pressure     (J/kg/K)

    real(prec) :: latent_heat_freezing      = 0.334e6        ! \Delta H_f   latent heat of freezing                           (J/kg)
    real(prec) :: latent_heat_vaporization  = 2.501e6        ! \Delta H_v   latent heat of vaporization                       (J/kg)
    real(prec) :: latent_heat_sublimation   = 2.835e6        ! \Delta H_s   latent heat of sumblimation                       (J/kg) 

    real(prec) :: gravity_acc               = 9.81           ! g            gravitational acceleration                        (m/s2)

    real(prec) :: gas_constant_air          = 287.05         ! R_d          dry air gas constant                              (J/kg/K)
    real(prec) :: gas_constant_vapor        = 461.51         ! R_v          water vapor gas constant                          (J/kg/K)

    real(prec) :: stefan_boltzmann_constant = 5.67e-8        ! \sigma       Stefan-Boltzmann constant                         (W/m2/K4)

    real(prec) :: u_min_evap                = 0.01           ! u_{min,evap} minimum velocity for CCLM evaporation calculation (m/s)

  end type default_values_type

  type(default_values_type), parameter :: default_values = default_values_type()

end module flux_constants


