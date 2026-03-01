subroutine blstep

  use blmod, only: rho, temp, ye, abar, eps, p, cs2, vel, r, cr, scratch_step, &
                    time, dtime, shockpos_stop, nt, tdump_scalar, &
                    diag_scratch_count, diag_nickel_called, diag_analysis_called
  use parameters
  use physical_constants
  implicit none

  real*8 :: rho_save(imax),temp_save(imax),vel_save(imax)
  real*8 :: p_save(imax),ye_save(imax),cs2_save(imax), abar_save(imax)
  real*8 :: dedt_save(imax), r_save(imax), cr_save(imax)
  real*8 :: eps_save(imax)

  integer :: iterations

!------------------------------------------------------------------------------

  ! save old values in case we have to redo the step
  rho_save   = rho
  temp_save  = temp
  ye_save    = ye
  abar_save    = abar
  eps_save    = eps
  p_save     = p
  cs2_save   = cs2
  vel_save  = vel
  r_save     = r
  cr_save    = cr

  iterations = 0
  diag_scratch_count = 0
  diag_nickel_called = .false.
  diag_analysis_called = .false.

  do while(scratch_step.or.(iterations.eq.0))

     iterations = iterations + 1

     if(scratch_step) then
        ! redo step with half the timestep; this is 
        ! sometimes necessary if the solver does not converge
        if (iterations.gt. 10) then
           stop "Stopping evolution. Time step repeated 10 times without luck"
        endif
        write(6,*) "Scratching entire step...", nt
        write(6,"(A25,1P10E15.6)") "time,dt,dt_new: ",time,dtime,dtime/2.0d0
        dtime = dtime / 2.0d0
        rho = rho_save
        temp = temp_save
        ye = ye_save
        abar = abar_save
        eps = eps_save
        p = p_save
        cs2 = cs2_save
        vel = vel_save
        r = r_save
        cr = cr_save
        diag_scratch_count = diag_scratch_count + 1
        scratch_step = .false.
     endif

     ! select between pure hydro or
     ! radiation+hydro solver
     if(radiation) then
        call hydro_rad
     else
        call hydro
     end if

          if (time.ge.bomb_tend .and. shockpos_stop.eq.0) then
        call shock_capture
     endif

     ! Only do the heavy diagnostics when we actually need scalar output.
     ! time is the "old" time here (same one used in snec.F90 for tdump_scalar).
     if (.not.sedov) then
        if (time.eq.0.0d0 .or. time.ge.tdump_scalar) then
           call analysis
           diag_analysis_called = .true.
        endif
     endif

     ! Energy conservation: we already did an initial call in problem.F90.
     ! Here we just sample it on the scalar dump cadence.
     if (time.ge.tdump_scalar) then
        call conservation_compute_energies
     endif

  enddo

end subroutine blstep
