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

  real*8 :: Ni_mass_g, m_cut, m_ej, q_mix, delta_q, q, q0
  real*8 :: sum_fdm, scaling, Ni_uniform_frac
  real*8, allocatable :: fNi(:)
  real*8, allocatable :: profile_delta_mass(:)
  real*8 :: pmass_inner, pmass_outer
  real*8 :: ejecta_inner_mass, ejecta_outer_mass
  real*8 :: mix_outer_mass, mix_outer_mass2
  real*8 :: q_mix2, f_outer
  real*8 :: Ni_core_mass, Ni_outer_mass
  real*8 :: support_core_mass, support_outer_mass
  real*8 :: core_fraction, outer_fraction
  real*8 :: diag_ni_total, diag_ni_above, diag_ni_below
  real*8 :: smooth_nominal_mass_msun
  real*8 :: smooth_ni_nominal_mass_msun
  integer :: smooth_number_iterations
  integer :: smooth_ni_number_iterations
  logical :: ni_override_active, preserve_ni_during_general
  integer :: i_excise, smooth_zones

  real*8, allocatable :: saved_ni(:)
  real*8 :: sum_other, target_frac, scale_frac
  integer :: j

  integer :: profile_zones
  integer :: i,l

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

  !---------------------------------------------------------------------------
  ! Place 56Ni using the selected kernel (Ni_mix_kernel):
  !   1 = step function (uniform box)
  !   2 = flat-core + exponential tail
  !   3 = two-component boxcar (inner core + outer redistributed shell)
  ! Ni_mix_fraction (0-1) sets q_mix: the fraction of ejecta mass over
  ! which Ni is distributed.  Boxcar smoothing (if enabled) is applied
  ! afterward to smooth the profile.
  !---------------------------------------------------------------------------
  pmass_inner = minval(pmass)
  pmass_outer = maxval(pmass)
  ejecta_inner_mass = max(pmass_inner, mass(1))
  ejecta_outer_mass = max(pmass_outer, ejecta_inner_mass)

  ! When the user controls Ni manually (Ni_by_hand=1), always zero out the
  ! progenitor's native Ni-56 first.  The progenitor profile from MESA can
  ! have up to ~80% Ni in the core, which causes EOS solver crashes if left
  ! in place.  The requested Ni_mass is then seeded by the kernel below.
  ! With Ni_mass=0 this gives a clean zero-Ni composition.
  if (Ni_by_hand.ne.0 .and. Ni_number.ne.0) then
     do i = 1, profile_zones
        pcomp(i, Ni_number) = 0.d0
     enddo
  endif

  if (Ni_by_hand.ne.0 .and. Ni_number.ne.0 .and. Ni_mass.gt.0.d0) then

     Ni_mass_g = Ni_mass*msun
     m_cut = ejecta_inner_mass
     m_ej  = ejecta_outer_mass - m_cut
     if (m_ej.le.0.d0) then
        write(*,*) 'ERROR in read_profile_compositions: m_ej <= 0'
        stop
     endif

     q_mix = Ni_mix_fraction
     if (q_mix.le.0.d0) then
        write(*,*) 'WARNING: Ni_mix_fraction <= 0; skipping Ni_by_hand remapping.'
        goto 100
     endif
     if (q_mix.gt.1.d0) q_mix = 1.0d0

     mix_outer_mass = m_cut + q_mix * m_ej

     if (Ni_mix_kernel.ne.1 .and. Ni_mix_kernel.ne.2 .and. Ni_mix_kernel.ne.3) then
        write(*,*) 'WARNING: unsupported Ni_mix_kernel =', Ni_mix_kernel, &
                   '; using kernel=2 (exponential).'
        Ni_mix_kernel = 2
     endif

     if (Ni_mix_kernel .eq. 1) then
        !---------------------------------------------------------------
        ! Kernel 1: Step function (uniform box)
        ! If the requested Ni mass cannot fit in the requested q_mix
        ! (i.e. X_Ni would exceed 1.0), expand q_mix until it fits.
        !---------------------------------------------------------------
        Ni_uniform_frac = Ni_mass_g / (q_mix * m_ej)
        if (Ni_uniform_frac .gt. 1.0d0) then
           q0 = Ni_mass_g / m_ej   ! minimum q_mix for X_Ni = 1.0
           q0 = min(q0, 1.0d0)
           write(*,'(A,F8.4,A,F8.4,A,F8.4)') &
                '  WARNING: Ni_mass cannot fit in q_mix=', q_mix, &
                ' (would need X_Ni=', Ni_uniform_frac, &
                '); expanding to q_mix=', q0
           q_mix = q0
           mix_outer_mass = m_cut + q_mix * m_ej
           Ni_uniform_frac = Ni_mass_g / (q_mix * m_ej)
        endif
        write(*,*) 'Ni mixing (step): kernel=1, q_mix =', q_mix, &
                   ' mix_outer_mass =', mix_outer_mass/msun, ' Msun'
        write(*,'(A,F8.5)') '    X_Ni = ', Ni_uniform_frac
        do i = 1, profile_zones
           if (pmass(i) .ge. m_cut .and. pmass(i) .le. mix_outer_mass) then
              pcomp(i,Ni_number) = Ni_uniform_frac
           else
              pcomp(i,Ni_number) = 0.d0
           endif
        enddo

     else if (Ni_mix_kernel .eq. 2) then
        !---------------------------------------------------------------
        ! Kernel 2: Flat core + exponential tail
        ! Flat at X_Ni=const for q in [0, q0] where q0 = q_mix/2,
        ! then exponential decay exp(-(q-q0)/delta_q) for q in (q0, q_mix].
        ! delta_q = (q_mix - q0)/3, so the profile drops to ~5% at q=q_mix.
        ! Numerically normalized to conserve Ni_mass exactly.
        !---------------------------------------------------------------
        allocate(profile_delta_mass(profile_zones))
        call build_effective_delta_mass(profile_zones, pmass, profile_delta_mass)
        allocate(fNi(profile_zones))

        write(*,*) 'Ni mixing (exponential): kernel=2, q_mix =', q_mix, &
                   ' mix_outer_mass =', mix_outer_mass/msun, ' Msun'

        q0 = 0.5d0*q_mix
        delta_q = (q_mix - q0)/3.d0
        if (delta_q.le.1.d-6) delta_q = 1.d-6

        sum_fdm = 0.d0
        do i = 1, profile_zones
           if (pmass(i) .lt. m_cut .or. pmass(i) .gt. mix_outer_mass) then
               fNi(i) = 0.d0
               cycle
           endif
           q = (pmass(i) - m_cut)/m_ej
           if (q.lt.0.d0) q = 0.d0
           if (q.gt.1.d0) q = 1.d0

           if (q.le.q0) then
              fNi(i) = 1.d0
           else if (q.le.q_mix) then
              fNi(i) = exp( -(q - q0)/delta_q )
           else
              fNi(i) = 0.d0
           endif

           sum_fdm = sum_fdm + fNi(i)*profile_delta_mass(i)
        enddo

        if (sum_fdm.le.0.d0) then
           write(*,*) 'WARNING: sum_fdm <= 0, skipping Ni_by_hand remapping.'
        else
           scaling = Ni_mass_g/sum_fdm
           do i = 1, profile_zones
              if (pmass(i) .lt. m_cut .or. pmass(i) .gt. mix_outer_mass) cycle
              pcomp(i,Ni_number) = scaling*fNi(i)
           enddo
        endif

        deallocate(fNi)
        deallocate(profile_delta_mass)

     else if (Ni_mix_kernel .eq. 3) then
        !---------------------------------------------------------------
        ! Kernel 3: Two-component boxcar with outer redistribution.
        ! Component 1 fills q in [0, q_mix] with (1-f_outer)*M_Ni.
        ! Component 2 fills q in (q_mix, q_mix2] with f_outer*M_Ni.
        !---------------------------------------------------------------
        allocate(profile_delta_mass(profile_zones))
        call build_effective_delta_mass(profile_zones, pmass, profile_delta_mass)

        q_mix2 = Ni_mix_component2_extent
        if (q_mix2.lt.0.d0) q_mix2 = q_mix
        if (q_mix2.gt.1.d0) q_mix2 = 1.0d0
        f_outer = Ni_mix_component2_fraction
        if (f_outer.lt.0.d0) f_outer = 0.d0
        if (f_outer.gt.1.d0) f_outer = 1.d0

        if (q_mix2.lt.q_mix) then
           write(*,*) 'WARNING: kernel=3 requires component2 extent >= Ni_mix_fraction.'
           write(*,*) '         Clamping Ni_mix_component2_extent to', q_mix
           q_mix2 = q_mix
        endif

        ! Expand regions if Ni mass cannot fit (X_Ni would exceed 1.0).
        ! Core needs (1-f_outer)*Ni_mass in q_mix*m_ej.
        ! Outer needs f_outer*Ni_mass in (q_mix2-q_mix)*m_ej.
        q0 = (1.d0 - f_outer) * Ni_mass_g / m_ej
        if (q0 .gt. q_mix) then
           write(*,'(A,F8.4,A,F8.4)') &
                '  WARNING: core Ni cannot fit in q_mix=', q_mix, &
                '; expanding to q_mix=', min(q0, 1.0d0)
           q_mix = min(q0, 1.0d0)
           mix_outer_mass = m_cut + q_mix * m_ej
        endif
        if (f_outer .gt. 0.d0) then
           q0 = q_mix + f_outer * Ni_mass_g / m_ej
           if (q0 .gt. q_mix2) then
              write(*,'(A,F8.4,A,F8.4)') &
                   '  WARNING: outer Ni cannot fit in q_mix2=', q_mix2, &
                   '; expanding to q_mix2=', min(q0, 1.0d0)
              q_mix2 = min(q0, 1.0d0)
           endif
        endif

        mix_outer_mass2 = m_cut + q_mix2 * m_ej
        Ni_outer_mass = f_outer * Ni_mass_g
        Ni_core_mass = Ni_mass_g - Ni_outer_mass
        if (Ni_core_mass.lt.0.d0) Ni_core_mass = 0.d0

        support_core_mass = 0.d0
        support_outer_mass = 0.d0
        do i = 1, profile_zones
           if (pmass(i) .lt. m_cut .or. pmass(i) .gt. mix_outer_mass2) cycle
           q = (pmass(i) - m_cut)/m_ej
           if (q.lt.0.d0) q = 0.d0
           if (q.gt.1.d0) q = 1.d0
           if (q.le.q_mix) then
              support_core_mass = support_core_mass + profile_delta_mass(i)
           else if (q.le.q_mix2) then
              support_outer_mass = support_outer_mass + profile_delta_mass(i)
           endif
        enddo

        if (Ni_core_mass.gt.0.d0 .and. support_core_mass.le.0.d0) then
           write(*,*) 'ERROR: kernel=3 core support mass is zero.'
           write(*,*) '       Increase Ni_mix_fraction.'
           stop
        endif

        if (Ni_outer_mass.gt.0.d0 .and. support_outer_mass.le.0.d0) then
           write(*,*) 'WARNING: kernel=3 outer shell has zero support mass.'
           write(*,*) '         Moving component-2 Ni mass back into core component.'
           Ni_core_mass = Ni_core_mass + Ni_outer_mass
           Ni_outer_mass = 0.d0
        endif

        core_fraction = 0.d0
        outer_fraction = 0.d0
        if (Ni_core_mass.gt.0.d0) core_fraction = Ni_core_mass / support_core_mass
        if (Ni_outer_mass.gt.0.d0) outer_fraction = Ni_outer_mass / support_outer_mass

        write(*,*) 'Ni mixing (two-component): kernel=3, q_core =', q_mix, &
                   ' q_outer =', q_mix2, ' f_outer =', f_outer
        write(*,*) '    mix_outer_mass(core)  =', mix_outer_mass/msun, ' Msun'
        write(*,*) '    mix_outer_mass(shell) =', mix_outer_mass2/msun, ' Msun'
        write(*,*) '    M_core =', Ni_core_mass/msun, ' Msun, M_outer =', &
                   Ni_outer_mass/msun, ' Msun'

        do i = 1, profile_zones
           if (pmass(i) .lt. m_cut .or. pmass(i) .gt. mix_outer_mass2) then
              pcomp(i,Ni_number) = 0.d0
              cycle
           endif
           q = (pmass(i) - m_cut)/m_ej
           if (q.lt.0.d0) q = 0.d0
           if (q.gt.1.d0) q = 1.d0

           if (q.le.q_mix) then
              pcomp(i,Ni_number) = core_fraction
           else if (q.le.q_mix2) then
              pcomp(i,Ni_number) = outer_fraction
           else
              pcomp(i,Ni_number) = 0.d0
           endif
        enddo

        if (maxval(pcomp(:,Ni_number)) .gt. 0.95d0) then
           write(*,'(A,F8.4)') '  WARNING: kernel=3 peak Ni fraction exceeds 0.95 (peak=', &
                maxval(pcomp(:,Ni_number))
        endif

        deallocate(profile_delta_mass)

     endif

  endif
100 continue

  call normalize_raw_composition(pcomp, profile_zones, ncomps, Ni_number)

  ! --- Ni mass diagnostics: after kernel + normalization ---
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

      write(*,'(A,F7.4,A,I3)') '  Boxcar smoothing (all species): nominal_mass=', &
           smooth_nominal_mass_msun, ' Msun, iterations=', smooth_number_iterations
      if (ni_override_active) then
         write(*,'(A,F7.4,A,I3)') '  Boxcar smoothing (Ni only):    nominal_mass=', &
              smooth_ni_nominal_mass_msun, ' Msun, iterations=', smooth_ni_number_iterations
      endif

      ! Save Ni profile before general boxcar if we either preserve Ni
      ! or re-smooth Ni in a dedicated second pass.
      if (preserve_ni_during_general) then
         allocate(saved_ni(profile_zones))
         saved_ni(:) = pcomp(:, Ni_number)
      endif

      call boxcar_profile(pmass(i_excise:), pcomp(i_excise:,:), &
                          smooth_zones, ncomps, &
                          smooth_nominal_mass_msun*msun, &
                          smooth_number_iterations)

      ! Restore pre-boxcar Ni and renormalize others if requested.
      if (preserve_ni_during_general) then
         if (maxval(saved_ni) .gt. 0.95d0) then
            write(*,'(A,F8.4)') '  WARNING: Ni fraction capped at 0.95 ' // &
                 '(peak was', maxval(saved_ni)
            write(*,'(A)')      '           Ni_mix_fraction is too small ' // &
                 'for given Ni_mass. Increase Ni_mix_fraction.'
         endif
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

      ! Optional second-pass Ni-only smoothing with dedicated parameters.
      if (ni_override_active) then
         call boxcar_profile(pmass(i_excise:), pcomp(i_excise:,Ni_number:Ni_number), &
                             smooth_zones, 1, &
                             smooth_ni_nominal_mass_msun*msun, &
                             smooth_ni_number_iterations)
         if (maxval(pcomp(:,Ni_number)) .gt. 0.95d0) then
            write(*,'(A,F8.4)') '  WARNING: Ni fraction capped at 0.95 after Ni-only smoothing ' // &
                 '(peak was', maxval(pcomp(:,Ni_number))
         endif
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

  call verify_mass_fraction_closure(comp, imax, ncomps, 'composition mapping (custom)')

  ! --- Ni mass diagnostics: after grid mapping ---
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
  if (output_mode .eq. 0) then
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


  !output all species fractions for validation (one line per zone, all ncomps columns)
  filename = trim(adjustl(outdir))//"/comp_init_all.dat"
  open(unit=668,file=trim(adjustl(filename)),status='unknown',form='formatted')
  write(668,'(A,I4)') '# ncomps = ', ncomps
  do i=1,imax
     write(668,'(I5.4)',advance='no') i
     do l=1,ncomps
        write(668,'(E25.16E3)',advance='no') comp(i,l)
     enddo
     write(668,*)
  enddo
  close(668)
  endif ! output_mode .eq. 0

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

  if (output_mode .eq. 0) then
     filename = trim(adjustl(outdir))//"/metallicity_init.dat"
     call output_screenshot(metallicity,filename,imax)
  endif


end subroutine read_profile_compositions
