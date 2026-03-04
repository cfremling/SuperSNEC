subroutine input_parser
!
! This routine parses the "parameters" file and sets various flags
! that are kept in the modules blmod and parameters.
!
  use blmod, only: gravity_switch, wipe_outdir
  use parameters
  implicit none

  character(len=128) :: cpstring
  character(len=500) :: rmstring
  logical :: opt
  logical :: outdirthere
  integer :: i_mode, c_mode
  integer :: legacy_boxcar_mode3_number_iterations
  real*8 :: legacy_boxcar_mode3_nominal_mass_msun

!------------------------------------------------------------------------------

  opt = .false.

!****************************** LAUNCH ****************************************

  call get_string_parameter('outdir',outdir,opt)

!************************* STELLAR PROFILE ************************************


  call get_string_parameter('profile_name',profile_name,opt)
  call get_string_parameter('comp_profile_name',composition_profile_name,opt)


!***************************** EXPLOSION **************************************

  call get_string_parameter('initial_data',initial_data,opt)

  if(initial_data.eq."Piston_Explosion") then
     call get_double_parameter('piston_vel',piston_vel,opt)
     call get_double_parameter('piston_tstart',piston_tstart,opt)
     call get_double_parameter('piston_tend',piston_tend,opt)
  endif

  if(initial_data.eq."Thermal_Bomb") then
     call get_double_parameter('final_energy',final_energy,opt)
     call get_double_parameter('bomb_tstart',bomb_tstart,opt)
     call get_double_parameter('bomb_tend',bomb_tend,opt)
     call get_double_parameter('bomb_mass_spread',bomb_mass_spread,opt)
     call get_integer_parameter('bomb_start_point',bomb_start_point,opt)
     call get_integer_parameter('bomb_mode',bomb_mode,.true.)
     if(bomb_mode.eq.-666) then
        bomb_mode = 1 ! default
     endif
  endif

!******************************** GRID ****************************************

   call get_integer_parameter('imax',imax,opt)
   call get_double_parameter('grid_update_interval_days',grid_update_interval_days,opt)
   call get_string_parameter('grid_mode',grid_mode,.true.)
   do i_mode=1,len_trim(grid_mode)
      c_mode = iachar(grid_mode(i_mode:i_mode))
      if (c_mode.ge.iachar('A') .and. c_mode.le.iachar('Z')) then
         grid_mode(i_mode:i_mode) = achar(c_mode + 32)
      endif
   enddo
   if (len_trim(grid_mode).eq.0) grid_mode = 'adaptive_runtime'
   if (trim(adjustl(grid_mode)).ne.'adaptive_runtime' .and. &
        trim(adjustl(grid_mode)).ne.'legacy_pattern') then
      write(*,*) '******* Unsupported grid_mode: ', trim(grid_mode)
      write(*,*) '******* Supported modes: adaptive_runtime, legacy_pattern'
      stop 1
   endif
   call get_string_parameter('grid_pattern_file',grid_pattern_file,.true.)
   if (len_trim(grid_pattern_file).eq.0) grid_pattern_file = 'tables/GridPattern.dat'
   call get_double_parameter('grid_surface_alpha',grid_surface_alpha,.true.)
   call get_double_parameter('grid_relax_days',grid_relax_days,.true.)
   call get_double_parameter('grid_min_cell_frac',grid_min_cell_frac,.true.)
   call get_integer_parameter('grid_debug',grid_debug,.true.)
   call get_double_parameter('grid_remap_qphoto_stop',grid_remap_qphoto_stop,.true.)
   call get_logical_parameter('mass_excision',mass_excision,opt)
  if(mass_excision) then
     call get_double_parameter('mass_excised',mass_excised,opt)
  end if

!****************************** EVOLUTION *************************************

  call get_logical_parameter('radiation',radiation,opt)
  call get_integer_parameter('eoskey',eoskey,opt)
  call get_integer_parameter('Ni_switch',Ni_switch,opt)
  call get_double_parameter('Ni_mass',Ni_mass,opt)
  call get_double_parameter('Ni_mix_fraction',Ni_mix_fraction,opt)
  call get_integer_parameter('Ni_mix_kernel',Ni_mix_kernel,.true.)
  if (Ni_mix_kernel.eq.-666) Ni_mix_kernel = 2 ! default: exponential
  call get_double_parameter('Ni_mix_component2_fraction',Ni_mix_component2_fraction,.true.)
  call get_double_parameter('Ni_mix_component2_extent',Ni_mix_component2_extent,.true.)
  if (Ni_mix_component2_extent.lt.0.0d0) Ni_mix_component2_extent = Ni_mix_fraction
  if (Ni_mix_component2_fraction.lt.0.0d0) then
     write(*,*) 'WARNING: Ni_mix_component2_fraction < 0; clamping to 0.'
     Ni_mix_component2_fraction = 0.0d0
  else if (Ni_mix_component2_fraction.gt.1.0d0) then
     write(*,*) 'WARNING: Ni_mix_component2_fraction > 1; clamping to 1.'
     Ni_mix_component2_fraction = 1.0d0
  endif
  if (Ni_mix_component2_extent.lt.0.0d0) then
     write(*,*) 'WARNING: Ni_mix_component2_extent < 0; clamping to 0.'
     Ni_mix_component2_extent = 0.0d0
  else if (Ni_mix_component2_extent.gt.1.0d0) then
     write(*,*) 'WARNING: Ni_mix_component2_extent > 1; clamping to 1.'
     Ni_mix_component2_extent = 1.0d0
  endif
  call get_double_parameter('Ni_period',Ni_period,opt)
  Ni_period_min_target = Ni_period
  call get_double_parameter('Ni_period_max',Ni_period_max_target,opt)
  if (Ni_period_max_target.le.0.0d0) then
     Ni_period_max_target = 10.0d0 * Ni_period_min_target
  endif
  call get_double_parameter('Ni_fractional_change',Ni_fractional_change_target,opt)
  if (Ni_fractional_change_target.le.0.0d0) then
     Ni_fractional_change_target = 0.25d0
  endif
  call get_integer_parameter('Ni_quad_npoints',Ni_quad_npoints,opt)
  if (Ni_quad_npoints.le.0) Ni_quad_npoints = 150
  call get_integer_parameter('Ni_by_hand',Ni_by_hand,.true.)
  if (Ni_by_hand.eq.-666) Ni_by_hand = 1 ! set to default value

  call get_integer_parameter('saha_ncomps',saha_ncomps,opt)
  
  call get_logical_parameter('boxcar_smoothing',boxcar_smoothing,opt)
  call get_integer_parameter('boxcar_smooth_ni',boxcar_smooth_ni,opt)
  call get_double_parameter('boxcar_nominal_mass_msun',boxcar_nominal_mass_msun,.true.)
  if (boxcar_nominal_mass_msun.le.0.0d0) then
     write(*,*) 'WARNING: boxcar_nominal_mass_msun <= 0; using default 0.4'
     boxcar_nominal_mass_msun = 0.4d0
  endif
  call get_integer_parameter('boxcar_number_iterations',boxcar_number_iterations,.true.)
  if (boxcar_number_iterations.eq.-666) then
     boxcar_number_iterations = 4
  else if (boxcar_number_iterations.le.0) then
     write(*,*) 'WARNING: boxcar_number_iterations <= 0; using default 4'
     boxcar_number_iterations = 4
  endif
  call get_double_parameter('boxcar_ni_nominal_mass_msun',boxcar_ni_nominal_mass_msun,.true.)
  legacy_boxcar_mode3_nominal_mass_msun = -1.23456789d307
  call get_double_parameter('boxcar_mode3_nominal_mass_msun',legacy_boxcar_mode3_nominal_mass_msun,.true.)
  if (boxcar_ni_nominal_mass_msun.eq.-1.0d0 .and. &
       legacy_boxcar_mode3_nominal_mass_msun.ne.-1.23456789d307) then
     boxcar_ni_nominal_mass_msun = legacy_boxcar_mode3_nominal_mass_msun
     write(*,*) 'WARNING: boxcar_mode3_nominal_mass_msun is deprecated; use boxcar_ni_nominal_mass_msun.'
  endif
  if (boxcar_ni_nominal_mass_msun.eq.0.0d0) then
     write(*,*) 'WARNING: boxcar_ni_nominal_mass_msun = 0; disabling secondary Ni-smoothing override.'
     boxcar_ni_nominal_mass_msun = -1.0d0
  else if (boxcar_ni_nominal_mass_msun.lt.0.0d0 .and. &
       boxcar_ni_nominal_mass_msun.ne.-1.0d0) then
     write(*,*) 'WARNING: boxcar_ni_nominal_mass_msun < 0; disabling secondary Ni-smoothing override.'
     boxcar_ni_nominal_mass_msun = -1.0d0
  endif

  call get_integer_parameter('boxcar_ni_number_iterations',boxcar_ni_number_iterations,.true.)
  if (boxcar_ni_number_iterations.eq.-666) then
     boxcar_ni_number_iterations = -1
  endif
  legacy_boxcar_mode3_number_iterations = -666
  call get_integer_parameter('boxcar_mode3_number_iterations',legacy_boxcar_mode3_number_iterations,.true.)
  if (boxcar_ni_number_iterations.eq.-1 .and. legacy_boxcar_mode3_number_iterations.ne.-666) then
     boxcar_ni_number_iterations = legacy_boxcar_mode3_number_iterations
     write(*,*) 'WARNING: boxcar_mode3_number_iterations is deprecated; use boxcar_ni_number_iterations.'
  endif
  if (boxcar_ni_number_iterations.eq.0) then
     write(*,*) 'WARNING: boxcar_ni_number_iterations = 0; disabling secondary Ni-smoothing override.'
     boxcar_ni_number_iterations = -1
  else if (boxcar_ni_number_iterations.lt.-1) then
     write(*,*) 'WARNING: boxcar_ni_number_iterations < -1; disabling secondary Ni-smoothing override.'
     boxcar_ni_number_iterations = -1
  endif

  call get_double_parameter('opacity_floor_envelope',of_env,opt)
  call get_double_parameter('opacity_floor_core',of_core,opt)

!********************** WHEN TO DO THINGS *************************************

  call get_integer_parameter('ntmax',ntmax,opt)
  call get_double_parameter('tend',tend,opt)
  call get_double_parameter('dtout',dtout,opt)
  dtout_late_target = dtout
  call get_double_parameter('dtout_mid',dtout_mid_target,opt)
  call get_double_parameter('dtout_fast',dtout_fast_target,opt)
  call get_double_parameter('dtout_scalar',dtout_scalar,opt)
  dtout_scalar_late_target = dtout_scalar
  call get_double_parameter('dtout_scalar_mid',dtout_scalar_mid_target,opt)
  call get_double_parameter('dtout_scalar_fast',dtout_scalar_fast_target,opt)
  call get_double_parameter('output_mid_transition_days',output_mid_transition_days,opt)
  call get_double_parameter('output_late_transition_days',output_late_transition_days,opt)
  call get_double_parameter('dtout_check',dtout_check,opt)
  call get_integer_parameter('ntout',ntout,opt)
  call get_integer_parameter('ntout_scalar',ntout_scalar,opt)
  call get_integer_parameter('ntout_check',ntout_check,opt)
  call get_integer_parameter('ntinfo',ntinfo,opt)
  call get_double_parameter('dtmin',dtmin,opt)
  call get_double_parameter('dtmax',dtmax,opt)
  
!*************************** SMOOTH NI LUMINOSITY *****************************

  call get_integer_parameter('smooth_ni_luminosity',smooth_ni_luminosity,.true.)
  if (smooth_ni_luminosity.eq.-666) smooth_ni_luminosity = 1
  write(*,'(A,I2)') ' Smooth Ni luminosity: ', smooth_ni_luminosity

!*************************** SOLVER TOLERANCES ********************************

  call get_double_parameter('epstol_hydro',epstol_hydro,.true.)
  if (epstol_hydro .le. 0.0d0) epstol_hydro = 1.0d-4
  call get_double_parameter('epstol_rad',epstol_rad,.true.)
  if (epstol_rad .le. 0.0d0) epstol_rad = 1.0d-4
  call get_integer_parameter('itmax_hydro',itmax_hydro,.true.)
  if (itmax_hydro .le. 0) itmax_hydro = 100
  call get_integer_parameter('itmax_rad',itmax_rad,.true.)
  if (itmax_rad .le. 0) itmax_rad = 300
  write(*,'(A,ES9.2,A,ES9.2,A,I4,A,I4)') &
       ' Solver: epstol_hydro=', epstol_hydro, &
       '  epstol_rad=', epstol_rad, &
       '  itmax_hydro=', itmax_hydro, &
       '  itmax_rad=', itmax_rad

!************************ NI RAY-TRACING OPTIMIZATION *************************

  call get_integer_parameter('ni_raytrace_opt',ni_raytrace_opt,.true.)
  if (ni_raytrace_opt.eq.-666) ni_raytrace_opt = 1

  call get_integer_parameter('ni_ray_interp',ni_ray_interp,.true.)
  if (ni_ray_interp.eq.-666) ni_ray_interp = 1

!****************************** OUTPUT MODE ************************************

  call get_integer_parameter('output_mode',output_mode,.true.)
  if (output_mode.eq.-666) output_mode = 0
  write(*,'(A,I2)') ' Output mode: ', output_mode

!********************************** TEST **************************************

  call get_logical_parameter('sedov',sedov,opt)
  if(sedov) then
     gravity_switch = 0
  else
     gravity_switch = 1
  end if

!******************************************************************************

! check if output directory exists
#if __INTEL_COMPILER
  inquire(directory=trim(adjustl(outdir)),exist=outdirthere)
#else
  inquire(file=trim(adjustl(outdir)),exist=outdirthere)
#endif
  if(.not.outdirthere) then
     write(6,*) "*** Output directory does not exist."
     write(6,*) "Please create the output directory: ", trim(adjustl(outdir))
     stop
  endif

! wipe output dir if requested:
  if(wipe_outdir) then
     write(*,*) "Removing output directory contents: ", trim(outdir)
     write(rmstring,*) "rm -rf ", trim(outdir), '/*'
     call system(rmstring)
  endif
  
! copy parameter file
  cpstring="cp parameters "//trim(adjustl(outdir))
  call system(cpstring)


!######## subroutines to get parse the parameters file ########################

contains
    subroutine get_string_parameter(parname,par,opt)

        implicit none
        logical opt
        character*(*) parname
        character*(*) par
        character*(200) line_string
        integer i,j,l,ll
        character*(200) temp_string
	integer :: par_len, ncopy

        open(unit=27,file='parameters',status='unknown')

        10 continue
        read(27,'(a)',end=19,err=10) line_string
        ! separator is an equal sign '=', # is comment
        i = index(line_string,'=')
        j = index(line_string,'#')

        if (i.eq.0.or.j.eq.1) goto 10
        !   if the whole line is a comment or there is no
        !   equal sign, then go on to the next line    

        if(j.gt.0.and.j.lt.i) goto 10
        !   if there is an equal sign, but it is in a comment
        !   then go on to the next line

        ! is this the right parameter? If not, cycle
        temp_string=trim(adjustl(line_string(1:i-1)))
        l=len(parname)
        if(parname.ne.temp_string(1:l)) goto 10

        !  If there is a comment in the line, exclude it!
        l = len(line_string)
        if (j.gt.0) l = j - 1
	
	! Clear the output string first
	par = ' '

	! Only copy if there *is* a substring
	if (l > i) then
 	  par_len = len(par)
	   ncopy   = min(par_len, l - i)
	   par(1:ncopy) = line_string(i+1:i+ncopy)
	endif

        ! now remove potential crap!
        do ll=1,len(par)
            if(par(ll:ll).eq.'\t') par(ll:ll) = ' '
            if(par(ll:ll).eq.'"') par(ll:ll) = ' '
            if(par(ll:ll).eq."'") par(ll:ll) = ' '
        enddo
        ! adjust left...
        par = adjustl(par)
        ! get rid of trailing blanks
        j = index(par," ")
        par = par(1:j-1)


        ! now look for " or ' and remove them
        j=index(par,'"')
        if(j.ne.0) stop "No quotes in my strings, please!"

        j=index(par,"'")
        if(j.ne.0) stop "No quotes in my strings, please!"

        close(27)
        return

        19 continue
        if(.not.opt) then
            write(6,*) "Fatal problem in input parser:"
            write(6,*) "Parameter ",parname
            write(6,*) "could not be read!"
            write(6,*) 
            call flush(6)
            stop
        else
            par = "NOTTHERE"
            close(27)
        endif

    end subroutine get_string_parameter

    subroutine get_double_parameter(parname,par,opt)

        implicit none
        logical opt
        character(*) parname
        character*256 line_string
        real*8 par

        call get_string_parameter(parname,line_string,opt)

        if(opt) then
            if(trim(adjustl(line_string)) .eq. 'NOTTHERE') return
        endif

        if(index(line_string,'.').eq.0) then
            write(6,*) "Uh. Bad double parameter ",trim(parname)
            write(6,*) "Please check input file!"
            call flush(6)
            stop
        endif

        read(line_string,"(e20.15)") par

    end subroutine get_double_parameter

    subroutine get_integer_parameter(parname,par,opt)

        implicit none
        logical opt
        character(*) parname
        character*256 line_string
        integer par

        call get_string_parameter(parname,line_string,opt)

        if((opt).and.trim(adjustl(line_string)) &
             .eq."NOTTHERE") then
           par = -666
        else
           read(line_string,"(i10)") par
        endif

    end subroutine get_integer_parameter

    subroutine get_logical_parameter(parname,par,opt)

        implicit none
        logical opt
        character*(*) parname
        character*(50) value_string
        integer temp_par
        logical par


        call get_string_parameter(parname,value_string,opt)

        if(opt) then
            ! don't try to set the parameter if it is not
            ! in the input file
            if(value_string .eq. "NOTTHERE") then
                write(6,*) "*** Parameter ",trim(adjustl(parname)), &
                    "not found in input file. Using default value."
                return
            endif
        endif
        read(value_string,"(i10)") temp_par

        if(temp_par.ne.0) then
            par = .true.
        else
            par = .false.
        endif

    end subroutine get_logical_parameter

end subroutine input_parser
