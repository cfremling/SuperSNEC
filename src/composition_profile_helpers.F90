module composition_profile_helpers
  implicit none
contains

  subroutine load_raw_composition_profile(prof_name, profile_zones, pmass, pradius, pcomp)
    use blmod, only: comp_details, ncomps
    implicit none
    character(*), intent(in) :: prof_name
    integer, intent(out) :: profile_zones
    real*8, allocatable, intent(out) :: pmass(:)
    real*8, allocatable, intent(out) :: pradius(:)
    real*8, allocatable, intent(out) :: pcomp(:,:)
    integer :: ibuffer, i
    integer, parameter :: unit = 667

    open(unit=unit,file=trim(prof_name),status='unknown',form='formatted',action='read')
    read(unit,*) profile_zones, ibuffer
    allocate(pmass(profile_zones))
    allocate(pradius(profile_zones))
    allocate(pcomp(profile_zones,ncomps))

    read(unit,*) comp_details(1:ncomps,1)
    read(unit,*) comp_details(1:ncomps,2)

    do i=1,profile_zones
       read(unit,*) pmass(i),pradius(i),pcomp(i,1:ncomps)
    enddo

    close(unit)
  end subroutine load_raw_composition_profile

  subroutine identify_major_species(H_number, He_number, C_number, O_number, Ni_number)
    use blmod, only: comp_details, ncomps
    implicit none
    integer, intent(out) :: H_number, He_number, C_number, O_number, Ni_number
    integer :: l

    H_number = 0
    He_number = 0
    C_number = 0
    O_number = 0
    Ni_number = 0

    do l=1,ncomps
      if(comp_details(l,2).eq.1.0d0 .and. comp_details(l,1).eq.1.0d0) then
          H_number = l
      else if(comp_details(l,2).eq.2.0d0 .and. comp_details(l,1).eq.4.0d0) then
          He_number = l
      else if(comp_details(l,2).eq.6.0d0 .and. comp_details(l,1).eq.12.0d0) then
          C_number = l
      else if(comp_details(l,2).eq.8.0d0 .and. comp_details(l,1).eq.16.0d0) then
          O_number = l
      else if(comp_details(l,2).eq.28.0d0 .and. comp_details(l,1).eq.56.0d0) then
          Ni_number = l
      end if
    enddo
  end subroutine identify_major_species

  subroutine normalize_raw_composition(pcomp, profile_zones, ncomps_in, Ni_number)
    implicit none
    integer, intent(in) :: profile_zones, ncomps_in, Ni_number
    real*8, intent(inout) :: pcomp(profile_zones,ncomps_in)
    integer :: i,l
    real*8 :: sum_initial, denom

    do i=1,profile_zones
        sum_initial = sum(pcomp(i,1:ncomps_in))
        if(sum_initial.le.0.d0) then
            write(*,*) 'ERROR: raw composition sum <= 0 at zone ', i
            stop
        endif
        if(Ni_number.eq.0) then
            pcomp(i,1:ncomps_in) = pcomp(i,1:ncomps_in)/sum_initial
        else
            denom = sum_initial - pcomp(i,Ni_number)
            if(denom.le.0.d0) then
                write(*,*) 'ERROR: invalid Ni normalization at zone ', i
                stop
            endif
            do l=1,ncomps_in
                if(l.ne.Ni_number) then
                    pcomp(i,l) = pcomp(i,l)*(1.d0-pcomp(i,Ni_number))/denom
                endif
            enddo
        endif
    enddo
  end subroutine normalize_raw_composition

  subroutine verify_mass_fraction_closure(abundances, zones, nspecies, tag)
    implicit none
    integer, intent(in) :: zones, nspecies
    real*8, intent(in) :: abundances(zones,nspecies)
    character(*), intent(in), optional :: tag
    real*8, parameter :: tolerance = 1.0d-8
    real*8 :: sum_val
    integer :: i

    do i=1,zones
        sum_val = sum(abundances(i,1:nspecies))
        if(abs(sum_val - 1.0d0).gt.tolerance) then
            if(present(tag)) then
                write(*,*) 'ERROR: mass fractions lost normalization after ', trim(tag), &
                            ' at zone ', i, '. Sum = ', sum_val
            else
                write(*,*) 'ERROR: mass fractions lost normalization at zone ', i, &
                            '. Sum = ', sum_val
            endif
            stop
        endif
    enddo
  end subroutine verify_mass_fraction_closure

end module composition_profile_helpers
