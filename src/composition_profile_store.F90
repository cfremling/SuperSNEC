module composition_profile_store
  use mixing_utils, only: build_effective_delta_mass
  implicit none
  real*8, parameter :: tiny_mass_interval = 1.0d-30
  integer :: stored_profile_zones = 0
  integer :: stored_profile_species = 0
  real*8, allocatable :: stored_profile_mass(:)
  real*8, allocatable :: stored_profile_comp(:,:)
  real*8, allocatable :: stored_profile_delta(:)
  real*8, allocatable :: stored_profile_edges(:)
contains

  subroutine store_processed_profile(pmass, pcomp, zones, species)
    implicit none
    integer, intent(in) :: zones, species
    real*8, intent(in) :: pmass(zones)
    real*8, intent(in) :: pcomp(zones,species)

    call release_stored_profile()

    if (zones <= 0 .or. species <= 0) then
       stored_profile_zones = 0
       stored_profile_species = 0
       return
    endif

    allocate(stored_profile_mass(zones))
    allocate(stored_profile_comp(zones,species))
    allocate(stored_profile_delta(zones))
    allocate(stored_profile_edges(zones+1))
    stored_profile_mass(1:zones) = pmass(1:zones)
    stored_profile_comp(1:zones,1:species) = pcomp(1:zones,1:species)
    call build_effective_delta_mass(zones, pmass, stored_profile_delta)
    call build_profile_edges(zones)
    stored_profile_zones = zones
    stored_profile_species = species
  end subroutine store_processed_profile

  subroutine release_stored_profile()
    if (allocated(stored_profile_mass)) deallocate(stored_profile_mass)
    if (allocated(stored_profile_comp)) deallocate(stored_profile_comp)
    if (allocated(stored_profile_delta)) deallocate(stored_profile_delta)
    if (allocated(stored_profile_edges)) deallocate(stored_profile_edges)
    stored_profile_zones = 0
    stored_profile_species = 0
  end subroutine release_stored_profile

  logical function has_stored_profile()
    has_stored_profile = stored_profile_zones > 0 .and. stored_profile_species > 0 &
         & .and. allocated(stored_profile_mass) .and. allocated(stored_profile_comp) &
         & .and. allocated(stored_profile_edges)
  end function has_stored_profile

  subroutine map_stored_profile_to_grid(target_edges, target_points, nspecies, comp_out)
    implicit none
    integer, intent(in) :: target_points, nspecies
    real*8, intent(in) :: target_edges(target_points)
    real*8, intent(out) :: comp_out(target_points,nspecies)
    integer :: species, new_idx
    integer :: old_idx
    real*8 :: left_edge, right_edge, denom, accumulator
    real*8 :: segment_end, overlap
    real*8 :: clip_left, clip_right
    real*8 :: outer_edge

    if (.not.has_stored_profile()) return
    if (nspecies > stored_profile_species) then
       write(*,*) 'ERROR: requested species exceed stored profile content'
       stop
    endif

    if (target_points <= 0) return

    outer_edge = stored_profile_edges(stored_profile_zones+1)

    do species=1,nspecies
       old_idx = 1
       do new_idx=1,target_points-1
          left_edge = target_edges(new_idx)
          right_edge = target_edges(new_idx+1)

          if (right_edge <= stored_profile_edges(1)) then
             comp_out(new_idx,species) = stored_profile_comp(1,species)
             cycle
          endif
          if (left_edge >= outer_edge) then
             comp_out(new_idx,species) = stored_profile_comp(stored_profile_zones,species)
             cycle
          endif

          clip_left = max(left_edge, stored_profile_edges(1))
          clip_right = min(right_edge, outer_edge)
          denom = max(tiny_mass_interval, clip_right - clip_left)
          accumulator = 0.0d0
          left_edge = clip_left
          right_edge = clip_right

          do while (old_idx < stored_profile_zones .and. &
               stored_profile_edges(old_idx+1) <= left_edge + tiny_mass_interval)
             old_idx = old_idx + 1
          enddo

          do
             if (left_edge >= right_edge - tiny_mass_interval) exit
             segment_end = min(right_edge, stored_profile_edges(old_idx+1))
             overlap = max(0.0d0, segment_end - left_edge)
             accumulator = accumulator + stored_profile_comp(old_idx,species) * overlap
             left_edge = segment_end
             if (left_edge >= right_edge - tiny_mass_interval) exit
             if (old_idx < stored_profile_zones) old_idx = old_idx + 1
          enddo

          comp_out(new_idx,species) = accumulator/denom
       enddo

       if (target_points .ge. 2) then
          comp_out(target_points,species) = comp_out(target_points-1,species)
       else
          comp_out(target_points,species) = stored_profile_comp(1,species)
       endif
    enddo
  end subroutine map_stored_profile_to_grid

  subroutine build_profile_edges(zones)
    implicit none
    integer, intent(in) :: zones
    integer :: i

    if (zones <= 0) return

    stored_profile_edges(1) = stored_profile_mass(1)
    if (zones >= 1) then
       stored_profile_edges(1) = stored_profile_mass(1) - 0.5d0*stored_profile_delta(1)
    endif

    do i=1,zones
       stored_profile_edges(i+1) = stored_profile_edges(i) + stored_profile_delta(i)
    enddo

    stored_profile_edges(1) = min(stored_profile_edges(1), stored_profile_mass(1))
    stored_profile_edges(zones+1) = max(stored_profile_edges(zones+1), stored_profile_mass(zones))

    do i=2,zones+1
       if (stored_profile_edges(i) <= stored_profile_edges(i-1)) then
          stored_profile_edges(i) = stored_profile_edges(i-1) + tiny_mass_interval
       endif
    enddo
  end subroutine build_profile_edges

end module composition_profile_store
