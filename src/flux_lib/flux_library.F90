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
  use flux_mass_evap,                 only: flux_mass_evap_cclm, flux_mass_evap_mom5, flux_mass_evap_meier
  !   heat fluxes
  use flux_heat_latent,               only: flux_heat_latent_ice, flux_heat_latent_water
  use flux_heat_sensible,             only: flux_heat_sensible_cclm, flux_heat_sensible_mom5, flux_heat_sensible_meier
  !   radiation fluxes
  use flux_radiation_blackbody,       only: flux_radiation_blackbody_StBo
  use distribute_radiation_flux_mod,  only: distribute_radiation_flux
  !   momentum fluxes
  use flux_momentum,                  only: flux_momentum_cclm, flux_momentum_mom5, flux_momentum_meier

  implicit none ; private

  ! expose all functions
  public flux_heat_latent_ice
  public flux_heat_latent_water
  public flux_heat_sensible_cclm
  public flux_heat_sensible_mom5
  public flux_heat_sensible_meier
  public flux_mass_evap_cclm
  public flux_mass_evap_mom5
  public flux_mass_evap_meier
  public flux_momentum_cclm
  public flux_momentum_mom5
  public flux_momentum_meier
  public flux_radiation_blackbody_StBo
  public distribute_radiation_flux
  public spec_vapor_surface_cclm

!contains

end module flux_library


