!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flux library for different options of           !
! calculating exchange fluxes between atmosphere  !
! and ocean/ice/land.                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reproduces formulae used in CCLM                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hagen Radtke, 2020                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module flux_radiation_blackbody

  use flux_constants, only: default_values_type, default_values, prec

  implicit none ; private

  ! expose all functions
  public flux_radiation_blackbody_StBo    ! Stefan-Boltzmann's law (assume ideal black body)

contains

  subroutine flux_radiation_blackbody_StBo(&
    flux_radiation_blackbody,               & ! RESULT       (W/m2)
    temperature_surface,                    & ! T_s          (K)
    stefan_boltzmann_constant_new           & ! \sigma       (W/m2/K4)
  )

    real(prec), intent(out) :: flux_radiation_blackbody               ! RESULT       (W/m2)
    real(prec), intent(in)  :: temperature_surface                    ! T_s          (K)
    real(prec), intent(in), optional :: stefan_boltzmann_constant_new ! \sigma       (W/m2/K4)

    real(prec) :: stefan_boltzmann_constant

    if (PRESENT(stefan_boltzmann_constant_new)) then 
        stefan_boltzmann_constant = stefan_boltzmann_constant_new
    else 
        stefan_boltzmann_constant = default_values%stefan_boltzmann_constant
    end if

    flux_radiation_blackbody = stefan_boltzmann_constant * temperature_surface**4
    
  end subroutine flux_radiation_blackbody_StBo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module flux_radiation_blackbody