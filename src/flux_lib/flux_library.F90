!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flux library for different options of           !
! calculating exchange fluxes between atmosphere  !
! and ocean/ice/land.                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reproduces formulae used in CCLM                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hagen Radtke, 2020                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module flux_library

  ! import auxiliary variable calculation from different modules
  ! calculation of specific vapor content
  use aux_specific_vapor_content, only: spec_vapor_surface_cclm

  ! import flux calculation functions from different modules:
  !   mass fluxes
  use flux_mass_evap,             only: flux_mass_evap_cclm
  !   heat fluxes
  use flux_heat_latent,           only: flux_heat_latent_ice, flux_heat_latent_water
  use flux_heat_sensible,         only: flux_heat_sensible_cclm
  !   radiation fluxes
  use flux_radiation_blackbody,   only: flux_radiation_blackbody_StBo
  !   momentum fluxes
  use flux_momentum,              only: flux_momentum_cclm

  implicit none ; private

  ! expose all functions
  public flux_heat_latent_ice
  public flux_heat_latent_water
  public flux_heat_sensible_cclm
  public flux_mass_evap_cclm
  public flux_momentum_cclm
  public flux_radiation_blackbody_StBo
  public spec_vapor_surface_cclm

!contains

end module flux_library


