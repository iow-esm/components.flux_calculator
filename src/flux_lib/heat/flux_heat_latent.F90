!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flux library for different options of           !
! calculating exchange fluxes between atmosphere  !
! and ocean/ice/land.                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reproduces formulae used in CCLM                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hagen Radtke, 2020                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module flux_heat_latent

  use flux_constants, only: default_values_type, default_values, prec

  implicit none ; private

  ! expose all functions
  public flux_heat_latent_ice
  public flux_heat_latent_water

contains

  subroutine flux_heat_latent_ice(&
    flux_heat_latent,                       & ! RESULT       (W/m2)
    flux_mass_evap,                         & !              (kg/m2/s)
    latent_heat_sublimation_new             & ! \Delta H_s   (J/kg)
  )

    real(prec), intent(out) :: flux_heat_latent                       ! RESULT       (W/m2)
    real(prec), intent(in)  :: flux_mass_evap                         !              (kg/m2/s)
    real(prec), intent(in), optional :: latent_heat_sublimation_new   ! \Delta H_s   (J/kg)

    real(prec) :: latent_heat_sublimation

    if (PRESENT(latent_heat_sublimation_new)) then 
      latent_heat_sublimation = latent_heat_sublimation_new
    else 
      latent_heat_sublimation = default_values%latent_heat_sublimation
    end if

    flux_heat_latent = flux_mass_evap * latent_heat_sublimation
    
  end subroutine flux_heat_latent_ice

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine flux_heat_latent_water(&
    flux_heat_latent,                       & ! RESULT       (W/m2)
    flux_mass_evap,                         & !              (kg/m2/s)
    latent_heat_vaporization_new            & ! \Delta H_v   (J/kg)
  )

    real(prec), intent(out) :: flux_heat_latent                       ! RESULT       (W/m2)
    real(prec), intent(in)  :: flux_mass_evap                         !              (kg/m2/s)
    real(prec), intent(in), optional :: latent_heat_vaporization_new  ! \Delta H_v   (J/kg)

    real(prec) :: latent_heat_vaporization

    if (PRESENT(latent_heat_vaporization_new)) then 
      latent_heat_vaporization = latent_heat_vaporization_new
    else 
      latent_heat_vaporization = default_values%latent_heat_vaporization
    end if

    flux_heat_latent = flux_mass_evap * latent_heat_vaporization
    
  end subroutine flux_heat_latent_water

end module flux_heat_latent


