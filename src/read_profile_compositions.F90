subroutine read_profile_compositions(prof_name)

  use blmod, only: comp, cmass, H_number, He_number, C_number, &
                    O_number, Ni_number, mass, ye, abar, metallicity, &
                    comp_details, ncomps, delta_mass
  use composition_profile_helpers, only: load_raw_composition_profile, &
       identify_major_species, normalize_raw_composition, &
       verify_mass_fraction_closure
  use composition_profile_store, only: store_processed_profile, map_stored_profile_to_grid
  use parameters
  use physical_constants
  use mixing_utils, only: build_effective_delta_mass
  implicit none

!input:
  character(*) :: prof_name

!local:
  real*8,allocatable :: pmass(:),pradius(:),pcomp(:,:)
  real*8 :: Ni_mass_fraction
  real*8 :: mix_outer_mass
  real*8 :: pmass_inner, pmass_outer
  real*8 :: ejecta_inner_mass, ejecta_outer_mass
  real*8 :: m_ej, q_mix
  real*8 :: diag_ni_total, diag_ni_above, diag_ni_below
  real*8 :: smooth_nominal_mass_msun
  real*8 :: smooth_ni_nominal_mass_msun
  integer :: smooth_number_iterations
  integer :: smooth_ni_number_iterations
  logical :: ni_override_active, preserve_ni_during_general
  real*8, allocatable :: profile_delta_mass(:)
  real*8, allocatable :: saved_ni(:)
  real*8 :: sum_other, target_frac, scale_frac
  integer :: i_excise, smooth_zones

  integer :: profile_zones
  integer :: i,l,j
  character(len=1024) :: filename


!------------------------------------------------------------------------------

  call load_raw_composition_profile(prof_name, profile_zones, pmass, pradius, pcomp)
  write(*,*) "We have ",profile_zones, "composition profile zones."

  call identify_major_species(H_number, He_number, C_number, O_number, Ni_number)

  if(Ni_switch.eq.1 .and. Ni_number.eq.0) then
      write(*,*) 'radioactive Ni is absent in the composition profile'
      write(*,*) 'please, add a column for it (see documentation of the code)'
      write(*,*) 'or put Ni_switch = 0'
      stop
  endif

  !seed Ni by hand as a step function before renormalization and boxcar
  pmass_inner = minval(pmass)
  pmass_outer = maxval(pmass)
  ejecta_inner_mass = max(pmass_inner, mass(1))
  ejecta_outer_mass = max(pmass_outer, ejecta_inner_mass)
  m_ej = ejecta_outer_mass - ejecta_inner_mass

  if (Ni_by_hand.ne.0 .and. Ni_number.ne.0) then
    q_mix = Ni_mix_fraction
    if (q_mix.le.0.d0 .or. m_ej.le.0.d0) then
        write(*,*) 'WARNING: Ni_mix_fraction <= 0 or no ejecta; skipping Ni_by_hand.'
    else
        if (q_mix.gt.1.d0) q_mix = 1.0d0
        mix_outer_mass = ejecta_inner_mass + q_mix * m_ej
        Ni_mass_fraction = Ni_mass*msun / (q_mix * m_ej)
        write(*,*) 'Ni mixing (step): q_mix =', q_mix, &
                   ' mix_outer_mass =', mix_outer_mass/msun, ' Msun'
        do i=1, profile_zones
            if (pmass(i) .lt. ejecta_inner_mass) then
                pcomp(i,Ni_number) = 0.0d0
            else if (pmass(i) .le. mix_outer_mass) then
                pcomp(i,Ni_number) = Ni_mass_fraction
            else
                pcomp(i,Ni_number) = 0.0d0
            endif
        enddo
    endif
  end if

  call normalize_raw_composition(pcomp, profile_zones, ncomps, Ni_number)

  ! --- Ni mass diagnostics: after kernel + normalize ---
  if (Ni_number.ne.0) then
     allocate(profile_delta_mass(profile_zones))
     call build_effective_delta_mass(profile_zones, pmass, profile_delta_mass)
     diag_ni_total = 0.d0
     diag_ni_above = 0.d0
     diag_ni_below = 0.d0
     do i = 1, profile_zones
        diag_ni_total = diag_ni_total + pcomp(i,Ni_number)*profile_delta_mass(i)
        if (pmass(i) .ge. mass(1)) then
           diag_ni_above = diag_ni_above + pcomp(i,Ni_number)*profile_delta_mass(i)
        else
           diag_ni_below = diag_ni_below + pcomp(i,Ni_number)*profile_delta_mass(i)
        endif
     enddo
     write(*,'(A)')         '  --- Ni mass diagnostics (after kernel + normalize) ---'
     write(*,'(A,ES12.5,A)') '    Total in profile:    ', diag_ni_total/msun, ' Msun'
     write(*,'(A,ES12.5,A)') '    Above excision:      ', diag_ni_above/msun, ' Msun'
     write(*,'(A,ES12.5,A)') '    Below excision:      ', diag_ni_below/msun, ' Msun'
  endif

  if(boxcar_smoothing) then
      ! Only smooth profile zones above the excision boundary so that
      ! the boxcar cannot diffuse species into the discarded core region.
      i_excise = 1
      do i = 1, profile_zones
         if (pmass(i) .ge. mass(1)) then
            i_excise = i
            exit
         endif
      enddo
      smooth_zones = profile_zones - i_excise + 1
      smooth_nominal_mass_msun = boxcar_nominal_mass_msun
      smooth_number_iterations = boxcar_number_iterations
      smooth_ni_nominal_mass_msun = smooth_nominal_mass_msun
      smooth_ni_number_iterations = smooth_number_iterations
      ni_override_active = .false.
      if (boxcar_smooth_ni .ne. 0 .and. Ni_number .ne. 0) then
         if (boxcar_ni_nominal_mass_msun.gt.0.d0) then
            smooth_ni_nominal_mass_msun = boxcar_ni_nominal_mass_msun
            ni_override_active = .true.
         endif
         if (boxcar_ni_number_iterations.gt.0) then
            smooth_ni_number_iterations = boxcar_ni_number_iterations
            ni_override_active = .true.
         endif
      endif
      preserve_ni_during_general = (boxcar_smooth_ni.eq.0 .or. ni_override_active) .and. &
           Ni_number.ne.0

      if (preserve_ni_during_general) then
         allocate(saved_ni(profile_zones))
         saved_ni(:) = pcomp(:, Ni_number)
      endif
      call boxcar_profile(pmass(i_excise:), pcomp(i_excise:,:), &
                          smooth_zones, ncomps, &
                          smooth_nominal_mass_msun*msun, &
                          smooth_number_iterations)

      if (preserve_ni_during_general) then
         pcomp(:, Ni_number) = min(saved_ni(:), 0.95d0)
         do i = i_excise, profile_zones
            sum_other = 0.0d0
            do j = 1, ncomps
               if (j /= Ni_number) sum_other = sum_other + pcomp(i,j)
            enddo
            target_frac = max(0.0d0, 1.0d0 - pcomp(i, Ni_number))
            if (sum_other > 0.0d0) then
               scale_frac = target_frac / sum_other
               do j = 1, ncomps
                  if (j /= Ni_number) pcomp(i,j) = pcomp(i,j) * scale_frac
               enddo
            endif
         enddo
      endif

      if (ni_override_active) then
         call boxcar_profile(pmass(i_excise:), pcomp(i_excise:,Ni_number:Ni_number), &
                             smooth_zones, 1, &
                             smooth_ni_nominal_mass_msun*msun, &
                             smooth_ni_number_iterations)
         pcomp(:, Ni_number) = min(pcomp(:, Ni_number), 0.95d0)
         do i = i_excise, profile_zones
            sum_other = 0.0d0
            do j = 1, ncomps
               if (j /= Ni_number) sum_other = sum_other + pcomp(i,j)
            enddo
            target_frac = max(0.0d0, 1.0d0 - pcomp(i, Ni_number))
            if (sum_other > 0.0d0) then
               scale_frac = target_frac / sum_other
               do j = 1, ncomps
                  if (j /= Ni_number) pcomp(i,j) = pcomp(i,j) * scale_frac
               enddo
            endif
         enddo
      endif

      if (preserve_ni_during_general) then
         deallocate(saved_ni)
      endif
  endif

  ! --- Ni mass diagnostics: after boxcar smoothing ---
  if (Ni_number.ne.0) then
     diag_ni_total = 0.d0
     diag_ni_above = 0.d0
     diag_ni_below = 0.d0
     do i = 1, profile_zones
        diag_ni_total = diag_ni_total + pcomp(i,Ni_number)*profile_delta_mass(i)
        if (pmass(i) .ge. mass(1)) then
           diag_ni_above = diag_ni_above + pcomp(i,Ni_number)*profile_delta_mass(i)
        else
           diag_ni_below = diag_ni_below + pcomp(i,Ni_number)*profile_delta_mass(i)
        endif
     enddo
     write(*,'(A)')         '  --- Ni mass diagnostics (after boxcar smoothing) ---'
     write(*,'(A,ES12.5,A)') '    Total in profile:    ', diag_ni_total/msun, ' Msun'
     write(*,'(A,ES12.5,A)') '    Above excision:      ', diag_ni_above/msun, ' Msun'
     write(*,'(A,ES12.5,A)') '    Below excision (LOST):', diag_ni_below/msun, ' Msun'
     deallocate(profile_delta_mass)
  endif

  call store_processed_profile(pmass, pcomp, profile_zones, ncomps)
  call map_stored_profile_to_grid(mass, imax, ncomps, comp)

  deallocate(pmass)
  deallocate(pradius)
  deallocate(pcomp)

!------------------------------------------------------------------------------

  call verify_mass_fraction_closure(comp, imax, ncomps, 'composition mapping')

  ! --- Ni mass diagnostics: on grid ---
  if (Ni_number.ne.0) then
     diag_ni_total = 0.d0
     do i = 1, imax
        diag_ni_total = diag_ni_total + comp(i,Ni_number)*delta_mass(i)
     enddo
     write(*,'(A)')         '  --- Ni mass diagnostics (on grid) ---'
     write(*,'(A,ES12.5,A)') '    Ni on grid:          ', diag_ni_total/msun, ' Msun'
     write(*,'(A,ES12.5,A)') '    Expected (parameter):', Ni_mass, ' Msun'
     write(*,'(A,F8.2,A)')   '    Ratio grid/expected: ', &
          100.d0*diag_ni_total/(Ni_mass*msun), ' %'
     diag_ni_total = 0.d0
     do i = 1, imax
        diag_ni_total = diag_ni_total + delta_mass(i)
     enddo
     write(*,'(A,ES12.5,A)') '    Ejecta mass on grid: ', diag_ni_total/msun, ' Msun'
  endif

  !set Y_e and abar
  do i=1,imax
      ye(i) = sum(comp(i,1:ncomps)* &
          (comp_details(1:ncomps,2)/comp_details(1:ncomps,1)))
      abar(i) = 1.0d0/sum(comp(i,1:ncomps)/comp_details(1:ncomps,1))
  enddo


  !output the initial fractions of some elements
  do l=1,ncomps
  if(l.eq.H_number) then !hydrogen
      filename = trim(adjustl(outdir))//"/H_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  else if(l.eq.He_number) then !helium
      filename = trim(adjustl(outdir))//"/He_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  else if(l.eq.C_number) then !carbon
      filename = trim(adjustl(outdir))//"/C_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  else if(l.eq.O_number) then !oxygen
      filename = trim(adjustl(outdir))//"/O_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  else if(l.eq.Ni_number) then !radioactive Ni
      filename = trim(adjustl(outdir))//"/Ni_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  end if
  enddo


  !calculate metallicity as 1 - He_fraction - H_fraction, output initial value
  do i=1,imax
      if(H_number.ne.0 .and. He_number.ne.0) then
          metallicity(i) = 1 - comp(i,H_number) - comp(i,He_number)
      else if(H_number.ne.0 .and. He_number.eq.0) then
          metallicity(i) = 1 - comp(i,H_number)
      else if(H_number.eq.0 .and. He_number.ne.0) then
          metallicity(i) = 1 - comp(i,He_number)
      else
          metallicity(i) = 1
      end if
  enddo

  filename = trim(adjustl(outdir))//"/metallicity_init.dat"
  call output_screenshot(metallicity,filename,imax)


end subroutine read_profile_compositions
