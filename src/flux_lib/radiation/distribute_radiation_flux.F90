module distribute_radiation_flux_mod

    use flux_constants, only: prec
  
    implicit none ; private
  
    ! expose all functions
    public distribute_radiation_flux    ! 
  
  contains
  
    subroutine distribute_radiation_flux(   &
      flux_radiation_surface_type,          & ! RESULT, surface-dependent radiation flux sent to bottom 
      flux_radiation_averaged,              & ! surface-independent radiation flux from atmosphere    
      albedo_averaged,                      & ! surface-independent albedo from atmosphere  
      albedo_surface_type                   & ! surface-dependent albedo received from ocean
    )
  
      real(prec), intent(out) :: flux_radiation_surface_type               ! RESULT     
      real(prec), intent(in)  :: flux_radiation_averaged, albedo_averaged, albedo_surface_type

  
      ! apply surface-type-dependent albedo and get rid of averaged albedo
      flux_radiation_surface_type = flux_radiation_averaged * (1.0 - albedo_surface_type) / (1.0 - albedo_averaged)

    end subroutine distribute_radiation_flux
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  end module distribute_radiation_flux_mod