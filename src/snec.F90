program snec

  use blmod, only: dtime, dtime_p, time, nt, ntstart, tstart,   &
     tdump, tdump_scalar, rho, tdump_check, time_Ni, &
     diag_hydro_rad_iters, diag_scratch_count, diag_nickel_called, &
     diag_analysis_called, diag_ni_eval_total, diag_ni_eval_remap
  use parameters
  use outinfomod, only: outinfo_count
  use adaptive_timing, only: compute_adaptive_periods
  use adaptive_grid_runtime, only: maybe_update_adaptive_grid, diag_remap_count
  implicit none

  logical :: OutputFlag = .false.
  logical :: OutputFlagScalar = .false.
  logical :: OutputFlagCheck = .false.
  integer :: timing_unit = 198
  integer(8) :: t_count_start, t_count_end, t_count_rate
  real*8 :: blstep_walltime

!------------------------------------------------------------------------------

  write(*,*)

  write(*,*) "***********************************"
  write(*,*) "* Supernova Explosion Code (SNEC) *"
  write(*,*) "***********************************"

  write(*,*)
! *****************************************************
! INITIALIZATION
! *****************************************************

  call input_parser

  call problem

  call artificial_viscosity

! output before first timestep
  call output_all(0)
  call output_all(1)
  call output_all(2)

  call timestep

  time = tstart
  nt = ntstart

  call compute_adaptive_periods(time, dtime, dtout, dtout_scalar, Ni_period)

  tdump_check = tstart+dtout_check
  tdump_scalar = time+dtout_scalar
  tdump = time+dtout

  call system_clock(count_rate=t_count_rate)
  if (output_mode .eq. 0) then
     open(unit=timing_unit, file=trim(adjustl(outdir))//'/blstep_timing.dat', &
          status='replace')
     write(timing_unit,'(A)') '# nt  time_s  time_d  dtime  blstep_ms  iters  scratches  Ni  analysis'
  endif

! *****************************************************
! MAIN LOOP
! *****************************************************

  IntegrationLoop: do

     dtime_p = dtime
     ! determine dt
     call timestep

     call compute_adaptive_periods(time, dtime, dtout, dtout_scalar, Ni_period)

     if (time.lt.tdump) then
        tdump = min(tdump, time + dtout)
     endif
     if (time.lt.tdump_scalar) then
        tdump_scalar = min(tdump_scalar, time + dtout_scalar)
     endif
     if (time.lt.time_Ni) then
        time_Ni = min(time_Ni, time + Ni_period)
     endif


     if(ntinfo.gt.0) then
        if(mod(nt,ntinfo).eq.0) then
           ! print useful info to stdout
           call outinfo
        endif
     endif

     if((time+dtime).gt.tend) dtime = tend-time

     ! actual integration step
     call system_clock(t_count_start)
     call blstep
     call system_clock(t_count_end)
     blstep_walltime = dble(t_count_end - t_count_start) / dble(t_count_rate) * 1.0d3
     if (output_mode .eq. 0) then
        write(timing_unit,'(I8,3ES16.7,F12.4,2I5,2L3)') nt, time, time/86400.0d0, dtime, &
             blstep_walltime, diag_hydro_rad_iters, diag_scratch_count, &
             diag_nickel_called, diag_analysis_called
     endif

     ! increment timestep
     nt = nt + 1

     ! various output related things
     if (ntout.gt.0) then
        if ( mod(nt,ntout) .eq. 0) OutputFlag = .true.
     endif

     if (ntout_scalar.gt.0) then
        if ( mod(nt,ntout_scalar) .eq. 0 ) OutputFlagScalar = .true.
     endif

     if (ntout_check.gt.0) then
        if ( mod(nt,ntout_check) .eq. 0 ) OutputFlagCheck = .true.
     endif

     if ( time.ge.tdump) then
        tdump=time+dtout
        OutputFlag = .true.
     endif

     if ( time.ge.tdump_scalar) then
        tdump_scalar=time+dtout_scalar
        OutputFlagScalar = .true.
     endif

     if ( time.ge.tdump_check) then
        tdump_check=tdump_check+dtout_check
        OutputFlagCheck = .true.
     endif

     ! increment time
     time = time+dtime

     call maybe_update_adaptive_grid(time)

     if (OutputFlag) then
        call output_all(0)
        call output_all(1)
        OutputFlag = .false.
     endif

     if (OutputFlagScalar) then
        call output_all(2)
        OutputFlagScalar = .false.
     endif

     if (OutputFlagCheck) then
        OutputFlagCheck = .false.
     endif

     if (time.eq.tend) then
        write(*,'(A,I8,A)') " Done! :-) tend reached after ", nt, " timesteps"
        write(*,'(A,I6,A,I6,A)') " Ni evals: ", diag_ni_eval_total, &
             "  (remap-forced resets: ", diag_ni_eval_remap, ")"
        write(*,'(A,I6)') " Grid remaps: ", diag_remap_count
        call output_all(0)
        call output_all(1)
        call output_all(2)
        call write_run_summary(nt)
        exit
     endif

     if (nt.ge.ntmax) then
        write(*,*) "Done! :-) ntmax reached"
        write(*,'(A,I6,A,I6,A)') " Ni evals: ", diag_ni_eval_total, &
             "  (remap-forced resets: ", diag_ni_eval_remap, ")"
        write(*,'(A,I6)') " Grid remaps: ", diag_remap_count
        call output_all(0)
        call output_all(1)
        call output_all(2)
        call write_run_summary(nt)
        exit
     endif

  enddo IntegrationLoop

  if (output_mode .eq. 0) close(timing_unit)

end program snec
