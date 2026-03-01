module adaptive_grid_runtime

  use blmod, only: mass, cmass, delta_mass, delta_cmass, r, cr, vel, rho, temp, &
        p, eps, ye, abar, metallicity, &
        tau, comp, ncomps, Ni_heating, Ni_deposit_function, Ni_total_luminosity, Ni_energy_rate, &
        Ni_number, H_number, He_number, comp_details, time, time_Ni, &
        time_Ni_last_fired, shockpos, shockpos_prev, diag_ni_eval_remap
  use parameters, only: imax, grid_update_interval_days, Ni_switch, Ni_period, &
       grid_mode, grid_surface_alpha, grid_relax_days, grid_min_cell_frac, grid_debug, &
       grid_adaptive_interval
  use physical_constants, only: pi
  use composition_profile_store, only: has_stored_profile, map_stored_profile_to_grid
  use composition_profile_helpers, only: verify_mass_fraction_closure
  implicit none

  real*8, parameter :: seconds_per_day = 86400.0d0
  real*8, parameter :: tiny_mass_interval = 1.0d-30
  real*8 :: next_grid_update_time = huge(1.0d0)
  real*8 :: last_remap_time = 0.0d0
  logical :: runtime_enabled = .false.
  logical :: is_adaptive_mode = .false.
  integer :: diag_remap_count = 0
  ! Photosphere tracking for physics-driven interval
  real*8 :: prev_q_photo = -1.0d0
  real*8 :: prev_q_photo_time = 0.0d0
  real*8 :: smoothed_dqdt = 0.0d0

  ! Persistent work arrays — allocated once at init, reused every remap
  real*8, allocatable :: remap_r_old(:), remap_vel_old(:)
  real*8, allocatable :: remap_temp_old(:), remap_eps_old(:)
  real*8, allocatable :: remap_ye_old(:), remap_abar_old(:), remap_met_old(:)
  real*8, allocatable :: remap_Ni_heat_old(:), remap_Ni_dep_old(:)
  real*8, allocatable :: remap_comp_old(:,:)
  ! Shared remap geometry — computed once per remap, reused by all scalar remaps
  real*8, allocatable :: remap_centers(:), remap_inv_dc(:)

  contains

  subroutine initialize_adaptive_grid(current_time)
    real*8, intent(in) :: current_time
    character(len=len(grid_mode)) :: mode

    mode = normalize_mode(grid_mode)
    if (len_trim(mode).eq.0) mode = 'adaptive_runtime'
    is_adaptive_mode = (mode .eq. 'adaptive_runtime')

    if (imax .le. 2 .or. .not.is_adaptive_mode) then
       runtime_enabled = .false.
       next_grid_update_time = huge(1.0d0)
       return
    endif

    runtime_enabled = .true.
    last_remap_time = current_time
    next_grid_update_time = current_time + &
         max(grid_update_interval_days, 0.01d0) * seconds_per_day

    ! Allocate persistent work arrays once
    if (.not.allocated(remap_r_old)) then
       allocate(remap_r_old(imax))
       allocate(remap_vel_old(imax))
       allocate(remap_temp_old(imax))
       allocate(remap_eps_old(imax))
       allocate(remap_ye_old(imax))
       allocate(remap_abar_old(imax))
       allocate(remap_met_old(imax))
       allocate(remap_Ni_heat_old(imax))
       allocate(remap_Ni_dep_old(imax))
       allocate(remap_centers(imax-1))
       allocate(remap_inv_dc(imax-1))
       if (ncomps .gt. 0) allocate(remap_comp_old(imax,ncomps))
    endif
  end subroutine initialize_adaptive_grid

  subroutine maybe_update_adaptive_grid(current_time)
    ! grid_adaptive_interval selects the remap scheduling mode:
    !   0     = fixed interval (grid_update_interval_days, default)
    !   1-999 = linear: interval = max(base, time_days / N)  [e.g. 15 = t/15]
    !  -1     = photosphere-driven: interval ~ q_thresh / |dq_photo/dt|
    !            with q_thresh=0.01 and t/20 safety cap
    !  -2     = photosphere-driven (uncapped): q_thresh=0.01, max 5d
    !  -3     = early-dense: 0.25d for t<3d, 0.5d for 3-10d, then t/20
    real*8, intent(in) :: current_time
    real*8 :: m_photo, total_mass, q_photo
    real*8 :: dt_since_remap
    real*8 :: adaptive_interval, time_days
    real*8 :: dt_photo, abs_dqdt, q_thresh
    integer :: i

    if (.not.runtime_enabled) return
    if (current_time + 1.0d-9 .lt. next_grid_update_time) return

    ! Elapsed time since last remap (used for relaxation alpha)
    dt_since_remap = current_time - last_remap_time

    ! Find photosphere mass coordinate (tau = 2/3, scanning inward)
    m_photo = mass(imax)
    do i=imax-1,1,-1
       if (tau(i) .ge. 0.6667d0) then
          m_photo = 0.5d0*(mass(i) + mass(i+1))
          exit
       endif
    enddo

    total_mass = mass(imax) - mass(1)
    if (total_mass .le. tiny_mass_interval) return
    q_photo = (m_photo - mass(1)) / total_mass

    call rebuild_lagrangian_grid(q_photo, dt_since_remap)

    last_remap_time = current_time
    time_days = current_time / seconds_per_day

    ! --- Compute next remap interval based on mode ---
    if (grid_adaptive_interval .gt. 0 .and. &
        grid_adaptive_interval .le. 999) then
       ! Linear: interval = max(base, time_days / N)
       adaptive_interval = max(grid_update_interval_days, &
            time_days / dble(grid_adaptive_interval))

    else if (grid_adaptive_interval .eq. -1 .or. &
             grid_adaptive_interval .eq. -2) then
       ! Photosphere-driven: interval = q_thresh / |dq_photo/dt|
       ! Tracks photosphere recession rate and remaps only when it
       ! would shift by more than q_thresh between updates.
       q_thresh = 0.01d0
       if (prev_q_photo .ge. 0.0d0 .and. &
            current_time .gt. prev_q_photo_time + 1.0d0) then
          dt_photo = (current_time - prev_q_photo_time) / seconds_per_day
          abs_dqdt = abs(q_photo - prev_q_photo) / dt_photo
          ! Exponential smoothing (alpha=0.5)
          if (smoothed_dqdt .le. 0.0d0) then
             smoothed_dqdt = abs_dqdt
          else
             smoothed_dqdt = 0.5d0*abs_dqdt + 0.5d0*smoothed_dqdt
          endif
       endif
       prev_q_photo = q_photo
       prev_q_photo_time = current_time
       if (smoothed_dqdt .gt. 1.0d-10) then
          adaptive_interval = max(grid_update_interval_days, &
               min(q_thresh / smoothed_dqdt, 5.0d0))
       else
          ! No rate estimate yet or photosphere static: use base interval
          adaptive_interval = grid_update_interval_days
       endif
       ! Mode -1: also cap at t/20 for safety
       if (grid_adaptive_interval .eq. -1) then
          adaptive_interval = min(adaptive_interval, &
               max(grid_update_interval_days, time_days / 20.0d0))
       endif

    else if (grid_adaptive_interval .eq. -3) then
       ! Early-dense: very frequent at early times, then t/20 at late times
       if (time_days .lt. 3.0d0) then
          adaptive_interval = 0.25d0
       else if (time_days .lt. 10.0d0) then
          adaptive_interval = 0.5d0
       else
          adaptive_interval = max(grid_update_interval_days, time_days / 20.0d0)
       endif

    else
       ! Default: fixed interval
       adaptive_interval = max(grid_update_interval_days, 0.01d0)
    endif

    next_grid_update_time = current_time + adaptive_interval * seconds_per_day
  end subroutine maybe_update_adaptive_grid

  subroutine rebuild_lagrangian_grid(q_photo, dt_elapsed)
    ! Piecewise uniform+geometric grid with:
    !   - inner_frac derived from photosphere position (q_photo)
    !   - alpha derived from grid_surface_alpha, scaled by q_photo
    !   - Exponential relaxation toward target (grid_relax_days)
    !   - Global shift limiter (max_allowed_rel_shift = 0.50)
    real*8, intent(in) :: q_photo
    real*8, intent(in) :: dt_elapsed   ! seconds since last remap
    real*8 :: mass_target(imax), mass_old(imax), mass_new(imax)
    real*8 :: shock_mass, shock_prev_mass
    real*8 :: inner_frac, r_eff, alpha
    real*8 :: relax_alpha, eff_relax_days, relax_scale
    real*8 :: span, mstart, mend
    real*8 :: max_shift, shift_tol, dt_days
    real*8 :: abs_shift, local_width, max_rel_shift
    integer :: i, n_transition, n_inner, n_outer
    real*8 :: dm_inner, outer_sum, rpower
    ! Maximum relative shift any boundary may move per remap (fraction of
    ! local cell width).  Limits remap-induced interpolation diffusion.
    real*8, parameter :: max_allowed_rel_shift = 0.50d0

    if (imax .le. 2) return

    ! Scale grid_relax_days with imax:
    !   imax <= 80:  eff = 0.4 * grid_relax_days  (=2 at baseline 5)
    !   imax = 100:  eff = 1.0 * grid_relax_days  (=5 at baseline 5)
    !   imax >= 200: eff = 2.0 * grid_relax_days  (=10 at baseline 5)
    !   Piecewise linear in between.
    if (imax .le. 80) then
       relax_scale = 0.4d0
    else if (imax .le. 100) then
       relax_scale = 0.4d0 + 0.6d0 * dble(imax - 80) / 20.0d0
    else if (imax .le. 200) then
       relax_scale = 1.0d0 + 1.0d0 * dble(imax - 100) / 100.0d0
    else
       relax_scale = 2.0d0
    endif
    eff_relax_days = grid_relax_days * relax_scale

    ! Compute relaxation alpha from actual elapsed time and e-folding timescale
    dt_days = dt_elapsed / seconds_per_day
    if (eff_relax_days .gt. 0.0d0) then
       relax_alpha = 1.0d0 - exp(-dt_days / eff_relax_days)
    else
       relax_alpha = 1.0d0
    endif
    relax_alpha = min(max(relax_alpha, 0.0d0), 1.0d0)
    if (relax_alpha .le. 0.0d0) return

    mstart = mass(1)
    mend = mass(imax)
    span = mend - mstart
    if (span .le. tiny_mass_interval) return

    ! Derive inner_frac from photosphere position
    inner_frac = 0.15d0 + 0.30d0 * min(max(q_photo, 0.0d0), 1.0d0)
    inner_frac = min(max(inner_frac, 0.10d0), 0.60d0)

    n_transition = nint(inner_frac * dble(imax))
    n_transition = max(3, min(imax-3, n_transition))
    n_inner = n_transition - 1
    n_outer = imax - n_transition

    ! Derive per-step ratio from alpha (total log concentration).
    ! Scale alpha by q_photo: full concentration when photosphere is at
    ! surface (q~1), reduced when deep (q~0.35).
    alpha = max(grid_surface_alpha, 0.1d0) &
         * max(0.35d0, min(1.0d0, q_photo))
    if (n_outer .gt. 0) then
       r_eff = exp(-alpha / dble(n_outer))
       if (grid_min_cell_frac .gt. 0.0d0 .and. n_outer .gt. 1) then
          r_eff = max(r_eff, &
               grid_min_cell_frac ** (1.0d0 / dble(n_outer - 1)))
       endif
       r_eff = max(0.01d0, min(0.9999d0, r_eff))
    else
       r_eff = 1.0d0
    endif

    ! Build piecewise target: uniform inner + geometric outer
    outer_sum = (1.0d0 - r_eff**n_outer) / (1.0d0 - r_eff)
    dm_inner = span / (dble(n_inner) + outer_sum)

    mass_target(1) = mstart
    do i=2,n_transition
       mass_target(i) = mass_target(i-1) + dm_inner
    enddo
    rpower = 1.0d0
    do i=n_transition+1,imax
       mass_target(i) = mass_target(i-1) + dm_inner * rpower
       rpower = rpower * r_eff
    enddo
    mass_target(imax) = mend

    ! Always count remaps (used by end-of-run log)
    diag_remap_count = diag_remap_count + 1

    ! Diagnostic dump (before blending, shows raw target)
    if (grid_debug .ge. 1) then
       block
         character(len=80) :: diag_fname
         write(diag_fname,'(A,I4.4,A)') 'grid_diag_', diag_remap_count, '.dat'
         open(unit=199, file=trim(diag_fname), status='replace')
         write(199,'(A,ES14.6)') '# time_s = ', time
         write(199,'(A,ES14.6)') '# time_d = ', time/seconds_per_day
         write(199,'(A,F10.6)') '# q_photo = ', q_photo
         write(199,'(A,F10.6)') '# inner_frac = ', inner_frac
         write(199,'(A,F10.6)') '# alpha = ', alpha
         write(199,'(A,F10.6)') '# r_eff = ', r_eff
         write(199,'(A,I6)') '# n_inner = ', n_inner
         write(199,'(A,I6)') '# n_outer = ', n_outer
         write(199,'(A,F10.6)') '# relax_alpha = ', relax_alpha
         write(199,'(A,I6)') '# remap_number = ', diag_remap_count
         write(199,'(A)') '# i  mass_cur  mass_target  dm_cur  dm_target  tau'
         do i = 1, imax-1
            write(199,'(I5,5ES16.7)') i, mass(i), mass_target(i), &
                 mass(i+1)-mass(i), mass_target(i+1)-mass_target(i), &
                 tau(i)
         enddo
         write(199,'(I5,5ES16.7)') imax, mass(imax), mass_target(imax), &
              0.0d0, 0.0d0, tau(imax)
         close(199)
       end block
    endif

    ! Limit relax_alpha globally so no boundary moves more than
    ! max_allowed_rel_shift of its local cell width per remap.
    ! Global (not per-cell) scaling preserves grid smoothness.
    max_rel_shift = 0.0d0
    do i=2,imax-1
       abs_shift = abs(relax_alpha * (mass_target(i) - mass(i)))
       local_width = mass(i+1) - mass(i-1)
       if (local_width .gt. tiny_mass_interval) then
          max_rel_shift = max(max_rel_shift, abs_shift / local_width)
       endif
    enddo
    if (max_rel_shift .gt. max_allowed_rel_shift) then
       relax_alpha = relax_alpha * max_allowed_rel_shift / max_rel_shift
    endif

    ! Blend current grid toward target
    mass_old = mass
    mass_new(1) = mass_old(1)
    mass_new(imax) = mass_old(imax)
    do i=2,imax-1
       mass_new(i) = (1.0d0 - relax_alpha) * mass_old(i) &
                    + relax_alpha * mass_target(i)
    enddo

    ! Avoid unnecessary remaps
    max_shift = 0.0d0
    do i=2,imax-1
       max_shift = max(max_shift, abs(mass_new(i) - mass_old(i)))
    enddo
    shift_tol = 1.0d-4 * span
    if (max_shift .le. shift_tol) return

    ! Enforce strict monotonicity (forward pass)
    do i=2,imax
       if (mass_new(i) .le. mass_new(i-1)) then
          mass_new(i) = mass_new(i-1) + tiny_mass_interval
       endif
    enddo
    ! Backward pass
    do i=imax-1,1,-1
       if (mass_new(i) .ge. mass_new(i+1)) then
          mass_new(i) = mass_new(i+1) - tiny_mass_interval
       endif
    enddo

    shock_mass = mass_old(max(1, min(imax, shockpos)))
    shock_prev_mass = mass_old(max(1, min(imax, shockpos_prev)))

    call remap_state(mass_old, mass_new, shock_mass, shock_prev_mass)
  end subroutine rebuild_lagrangian_grid

  subroutine remap_state(mass_old, mass_new, shock_mass, shock_prev_mass)
    real*8, intent(in) :: mass_old(imax)
    real*8, intent(in) :: mass_new(imax)
    real*8, intent(in) :: shock_mass, shock_prev_mass
    real*8 :: dm_l, dm_r
    integer :: i

    ! Save old state into persistent arrays
    remap_r_old = r
    remap_vel_old = vel
    remap_temp_old = temp
    remap_eps_old = eps
    remap_ye_old = ye
    remap_abar_old = abar
    remap_met_old = metallicity
    remap_Ni_heat_old = Ni_heating
    remap_Ni_dep_old = Ni_deposit_function
    if (ncomps .gt. 0) remap_comp_old = comp

    ! Interpolate boundary-defined quantities (r, vel)
    call interpolate_boundaries(remap_r_old, mass_old, mass_new, r)
    call interpolate_boundaries(remap_vel_old, mass_old, mass_new, vel)

    ! Update mass grid and derived differences
    mass = mass_new
    do i=1,imax-1
       delta_mass(i) = mass(i+1) - mass(i)
       cmass(i) = mass(i) + 0.5d0*delta_mass(i)
    enddo
    cmass(imax) = mass(imax) + 0.5d0*(mass(imax) - mass(imax-1))
    delta_mass(imax) = delta_mass(imax-1)
    do i=1,imax-1
       delta_cmass(i) = cmass(i+1) - cmass(i)
    enddo
    delta_cmass(imax) = delta_cmass(imax-1)

    ! Precompute old cell centers and validity flags for conservative remap
    do i=1,imax-1
       remap_centers(i) = 0.5d0*(mass_old(i) + mass_old(i+1))
    enddo
    remap_inv_dc(1) = 0.0d0
    remap_inv_dc(imax-1) = 0.0d0
    do i=2,imax-2
       dm_l = remap_centers(i) - remap_centers(i-1)
       dm_r = remap_centers(i+1) - remap_centers(i)
       if (dm_l .gt. tiny_mass_interval .and. dm_r .gt. tiny_mass_interval) then
          remap_inv_dc(i) = 1.0d0
       else
          remap_inv_dc(i) = 0.0d0
       endif
    enddo

    ! Conservatively remap temp and eps
    call remap_cell_scalar_fast(remap_temp_old, mass_old, mass_new, temp)
    call remap_cell_scalar_fast(remap_eps_old, mass_old, mass_new, eps)

    ! Remap composition
    if (has_stored_profile() .and. ncomps .gt. 0) then
       call map_stored_profile_to_grid(mass, imax, ncomps, comp)
       call verify_mass_fraction_closure(comp, imax, ncomps, 'adaptive remap')
       call rebuild_composition_scalars()
    else
       call remap_cell_scalar_fast(remap_ye_old, mass_old, mass_new, ye)
       call remap_cell_scalar_fast(remap_abar_old, mass_old, mass_new, abar)
       call remap_cell_scalar_fast(remap_met_old, mass_old, mass_new, metallicity)
       if (ncomps .gt. 0) then
          call remap_composition(remap_comp_old, mass_old, mass_new)
       endif
    endif

    if (Ni_switch .ne. 0) then
       call remap_cell_scalar_fast(remap_Ni_heat_old, mass_old, mass_new, Ni_heating)
       call remap_cell_scalar_fast(remap_Ni_dep_old, mass_old, mass_new, Ni_deposit_function)
    else
       Ni_heating(:) = 0.0d0
       Ni_deposit_function(:) = 0.0d0
    endif

    ! Recompute geometry (rho, cr) from new mass grid + interpolated r
    call rebuild_geometry()

    ! Boundary values
    temp(imax) = 0.0d0
    eps(imax) = 0.0d0
    p(imax) = 0.0d0
    rho(imax) = 0.0d0

    if (Ni_switch .ne. 0) then
       Ni_total_luminosity = Ni_energy_rate * sum(Ni_deposit_function(1:imax)*delta_mass(1:imax))
       time_Ni = min(time_Ni, time)
       diag_ni_eval_remap = diag_ni_eval_remap + 1
    else
       Ni_total_luminosity = 0.0d0
    endif

    shockpos = locate_mass_index(shock_mass, mass)
    shockpos_prev = locate_mass_index(shock_prev_mass, mass)
  end subroutine remap_state

  subroutine interpolate_boundaries(values_old, mass_old, mass_new, values_new)
    real*8, intent(in) :: values_old(imax)
    real*8, intent(in) :: mass_old(imax)
    real*8, intent(in) :: mass_new(imax)
    real*8, intent(out) :: values_new(imax)
    integer :: idx_old, j
    real*8 :: frac, denom

    idx_old = 1
    do j=1,imax
       if (mass_new(j) <= mass_old(1)) then
          values_new(j) = values_old(1)
       else if (mass_new(j) >= mass_old(imax)) then
          values_new(j) = values_old(imax)
       else
          do while (idx_old < imax-1 .and. mass_old(idx_old+1) < mass_new(j))
             idx_old = idx_old + 1
          enddo
          denom = mass_old(idx_old+1) - mass_old(idx_old)
          if (denom <= tiny_mass_interval) then
             frac = 0.0d0
          else
             frac = (mass_new(j) - mass_old(idx_old))/denom
          endif
          frac = max(0.0d0, min(1.0d0, frac))
          values_new(j) = values_old(idx_old) + frac*(values_old(idx_old+1) - values_old(idx_old))
       endif
    enddo
  end subroutine interpolate_boundaries

  subroutine remap_cell_scalar_fast(values_old, mass_old, mass_new, values_new)
    real*8, intent(in) :: values_old(imax)
    real*8, intent(in) :: mass_old(imax)
    real*8, intent(in) :: mass_new(imax)
    real*8, intent(out) :: values_new(imax)
    real*8 :: slopes(imax-1)
    real*8 :: grad_left, grad_right
    integer :: new_idx, old_idx, i
    real*8 :: left_edge, right_edge, denom, accumulator
    real*8 :: segment_end, center_i, slope_i, a, b

    slopes(1) = 0.0d0
    slopes(imax-1) = 0.0d0
    do i=2,imax-2
       if (remap_inv_dc(i) .gt. 0.5d0) then
          grad_left  = (values_old(i) - values_old(i-1)) / &
               (remap_centers(i) - remap_centers(i-1))
          grad_right = (values_old(i+1) - values_old(i)) / &
               (remap_centers(i+1) - remap_centers(i))
          if (grad_left * grad_right .le. 0.0d0) then
             slopes(i) = 0.0d0
          else if (abs(grad_left) .lt. abs(grad_right)) then
             slopes(i) = grad_left
          else
             slopes(i) = grad_right
          endif
       else
          slopes(i) = 0.0d0
       endif
    enddo

    old_idx = 1
    do new_idx=1,imax-1
       left_edge = mass_new(new_idx)
       right_edge = mass_new(new_idx+1)
       denom = max(tiny_mass_interval, right_edge - left_edge)
       accumulator = 0.0d0

       do
          if (left_edge >= right_edge - tiny_mass_interval) exit
          do while (old_idx < imax-1 .and. mass_old(old_idx+1) <= left_edge + tiny_mass_interval)
             old_idx = old_idx + 1
          enddo
          segment_end = min(right_edge, mass_old(old_idx+1))

          center_i = remap_centers(old_idx)
          slope_i = slopes(old_idx)

          a = left_edge
          b = segment_end
          accumulator = accumulator + values_old(old_idx) * (b - a) &
               + slope_i * 0.5d0 * ((b - center_i)**2 - (a - center_i)**2)

          left_edge = segment_end
          if (left_edge >= right_edge - tiny_mass_interval) exit
          if (old_idx < imax-1) old_idx = old_idx + 1
       enddo

       values_new(new_idx) = accumulator/denom
    enddo
    values_new(imax) = values_new(imax-1)
  end subroutine remap_cell_scalar_fast

  subroutine remap_composition(comp_old, mass_old, mass_new)
    real*8, intent(in) :: comp_old(imax,ncomps)
    real*8, intent(in) :: mass_old(imax)
    real*8, intent(in) :: mass_new(imax)
    real*8 :: work_old(imax)
    real*8 :: work_new(imax)
    integer :: species, i
    real*8 :: norm

    do species=1,ncomps
       work_old(:) = comp_old(:,species)
       call remap_cell_scalar_fast(work_old, mass_old, mass_new, work_new)
       comp(:,species) = work_new
    enddo

    do i=1,imax-1
       norm = sum(comp(i,1:ncomps))
       if (norm .gt. tiny_mass_interval) then
          comp(i,1:ncomps) = comp(i,1:ncomps)/norm
       endif
    enddo
    comp(imax,1:ncomps) = comp(imax-1,1:ncomps)
  end subroutine remap_composition

  subroutine rebuild_composition_scalars()
    integer :: i

    if (ncomps .le. 0) then
       ye(1:imax) = 0.0d0
       abar(1:imax) = 1.0d0
       metallicity(1:imax) = 0.0d0
       return
    endif

    do i=1,imax
       ye(i) = sum(comp(i,1:ncomps) * &
            (comp_details(1:ncomps,2)/comp_details(1:ncomps,1)))
       abar(i) = 1.0d0/sum(comp(i,1:ncomps)/comp_details(1:ncomps,1))
    enddo

    do i=1,imax
       if (H_number.ne.0 .and. He_number.ne.0) then
          metallicity(i) = 1.0d0 - comp(i,H_number) - comp(i,He_number)
       else if (H_number.ne.0 .and. He_number.eq.0) then
          metallicity(i) = 1.0d0 - comp(i,H_number)
       else if (H_number.eq.0 .and. He_number.ne.0) then
          metallicity(i) = 1.0d0 - comp(i,He_number)
       else
          metallicity(i) = 1.0d0
       endif
    enddo
  end subroutine rebuild_composition_scalars

  subroutine rebuild_geometry()
    integer :: i
    real*8 :: shell_volume

    do i=1,imax-1
       shell_volume = (4.0d0*pi/3.0d0)*(r(i+1)**3 - r(i)**3)
       if (shell_volume .gt. tiny_mass_interval) then
          rho(i) = delta_mass(i)/shell_volume
       else if (i .gt. 1) then
          rho(i) = rho(i-1)
       else
          rho(i) = rho(i+1)
       endif
    enddo
    rho(imax) = 0.0d0

    do i=1,imax-1
       cr(i) = ((r(i)**3 + r(i+1)**3)*0.5d0)**(1.0d0/3.0d0)
    enddo
    cr(imax) = r(imax) + (r(imax) - cr(imax-1))
  end subroutine rebuild_geometry

  integer function locate_mass_index(target, mass_array)
    real*8, intent(in) :: target
    real*8, intent(in) :: mass_array(imax)
    integer :: lo, hi, mid

    if (target <= mass_array(1)) then
       locate_mass_index = 1
       return
    endif
    if (target >= mass_array(imax)) then
       locate_mass_index = imax
       return
    endif

    lo = 1
    hi = imax
    do while (hi - lo .gt. 1)
       mid = (lo + hi) / 2
       if (mass_array(mid) .le. target) then
          lo = mid
       else
          hi = mid
       endif
    enddo
    locate_mass_index = lo
  end function locate_mass_index

  logical function runtime_grid_selected()
    runtime_grid_selected = is_adaptive_mode
  end function runtime_grid_selected

  pure function normalize_mode(raw) result(clean)
    character(len=*), intent(in) :: raw
    character(len=len(raw)) :: clean
    integer :: i, code

    clean = trim(adjustl(raw))
    do i=1,len(clean)
       code = iachar(clean(i:i))
       if (code.ge.iachar('A') .and. code.le.iachar('Z')) then
          clean(i:i) = achar(code + 32)
       endif
    enddo
  end function normalize_mode

end module adaptive_grid_runtime
