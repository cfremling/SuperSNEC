module mixing_utils
  implicit none
contains

  subroutine build_effective_delta_mass(zones, mass_coords, delta_mass_out)
    implicit none
    integer, intent(in) :: zones
    real*8, intent(in) :: mass_coords(zones)
    real*8, intent(out) :: delta_mass_out(zones)
    integer :: i
    real*8, parameter :: min_spacing = 1.0d-20

    if (zones <= 0) return

    if (zones == 1) then
       delta_mass_out(1) = max(min_spacing, mass_coords(1))
       return
    endif

    delta_mass_out(1) = mass_coords(2) - mass_coords(1)
    if (delta_mass_out(1) <= 0.0d0) delta_mass_out(1) = min_spacing

    do i = 2, zones-1
       delta_mass_out(i) = 0.5d0*(mass_coords(i+1) - mass_coords(i-1))
       if (delta_mass_out(i) <= 0.0d0) delta_mass_out(i) = min_spacing
    enddo

    delta_mass_out(zones) = mass_coords(zones) - mass_coords(zones-1)
    if (delta_mass_out(zones) <= 0.0d0) delta_mass_out(zones) = min_spacing

  end subroutine build_effective_delta_mass

end module mixing_utils
