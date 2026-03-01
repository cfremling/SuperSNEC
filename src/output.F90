subroutine output_all(modeflag)

  use blmod
  use parameters
  use physical_constants
  implicit none

  character(len=1024) :: filename
  character(len=256) :: basename
  integer :: modeflag

!------------------------------------------------------------------------------

  if(modeflag.eq.0) then

  ! meant for checkpoints; not used at the moment

  else if(modeflag.eq.1) then

    if (output_mode .eq. 0) then

    filename = trim(adjustl(outdir))//"/vel.xg"
    call output_single_mass(vel,filename)

    filename = trim(adjustl(outdir))//"/rho.xg"
    call output_single_mass(rho,filename)

    filename = trim(adjustl(outdir))//"/ye.xg"
    call output_single_mass(ye,filename)

    filename = trim(adjustl(outdir))//"/press.xg"
    call output_single_mass(p,filename)

    filename = trim(adjustl(outdir))//"/cs2.xg"
    call output_single_mass(cs2,filename)

    filename = trim(adjustl(outdir))//"/Q.xg"
    call output_single_mass(Q,filename)

    filename = trim(adjustl(outdir))//"/eps.xg"
    call output_single_mass(eps,filename)

    filename = trim(adjustl(outdir))//"/mass.xg"
    call output_single_radius(mass,filename)

    filename = trim(adjustl(outdir))//"/temp.xg"
    call output_single_mass(temp,filename)

    filename = trim(adjustl(outdir))//"/lum.xg"
    call output_single_mass(lum,filename)

    filename = trim(adjustl(outdir))//"/tau.xg"
    call output_single_mass(tau,filename)

    filename = trim(adjustl(outdir))//"/delta_time.xg"
    call output_single_mass(delta_time,filename)

    filename = trim(adjustl(outdir))//"/radius.xg"
    call output_single_mass(r,filename)

    filename = trim(adjustl(outdir))//"/kappa.xg"
    call output_single_mass(kappa,filename)

    filename = trim(adjustl(outdir))//"/kappa_table.xg"
    call output_single_mass(kappa_table,filename)

    filename = trim(adjustl(outdir))//"/logR_op.xg"
    call output_single_mass(logR_op,filename)

    filename = trim(adjustl(outdir))//"/logT.xg"
    call output_single_mass(logT,filename)

    filename = trim(adjustl(outdir))//"/p_rad.xg"
    call output_single_mass(p_rad,filename)

    filename = trim(adjustl(outdir))//"/Ni_deposit_function.xg"
    call output_single_mass(Ni_deposit_function,filename)

    filename = trim(adjustl(outdir))//"/He_1.xg"
    call output_single_mass(ion_fractions(He_number,1,:),filename)

    filename = trim(adjustl(outdir))//"/He_2.xg"
    call output_single_mass(ion_fractions(He_number,2,:),filename)

    filename = trim(adjustl(outdir))//"/He_3.xg"
    call output_single_mass(ion_fractions(He_number,3,:),filename)

    filename = trim(adjustl(outdir))//"/H_1.xg"
    call output_single_mass(ion_fractions(H_number,1,:),filename)

    filename = trim(adjustl(outdir))//"/H_2.xg"
    call output_single_mass(ion_fractions(H_number,2,:),filename)

    filename = trim(adjustl(outdir))//"/free_electron_frac.xg"
    call output_single_mass(free_electron_frac,filename)

    filename = trim(adjustl(outdir))//"/E_shell.xg"
    call output_single_mass(E_shell,filename)

    filename = trim(adjustl(outdir))//"/time_diff.xg"
    call output_single_mass(time_diff,filename)

    filename = trim(adjustl(outdir))//"/time_exp.xg"
    call output_single_mass(time_exp,filename)

    filename = trim(adjustl(outdir))//"/photosphere_tracer.xg"
    call output_single_mass(photosphere_tracer,filename)

    endif ! output_mode .eq. 0

  else if(modeflag.eq.2) then

     ! These files are always written (needed by fitting codes)
     filename = trim(adjustl(outdir))//"/lum_observed.dat"
     call output_scalar(lum_observed,filename)

     filename = trim(adjustl(outdir))//"/vel_photo.dat"
     call output_scalar(vel_photo,filename)

     ! These files are only written in full output mode
     if (output_mode .eq. 0) then

     filename = trim(adjustl(outdir))//"/T_eff.dat"
     call output_scalar(T_eff,filename)

     filename = trim(adjustl(outdir))//"/Ni_total_luminosity.dat"
     call output_scalar(Ni_total_luminosity,filename)

     filename = trim(adjustl(outdir))//"/index_photo.dat"
     call output_integer(index_photo,filename)

     filename = trim(adjustl(outdir))//"/lum_photo.dat"
     call output_scalar(lum_photo,filename)

     filename = trim(adjustl(outdir))//"/mass_photo.dat"
     call output_scalar(mass_photo,filename)

     filename = trim(adjustl(outdir))//"/rad_photo.dat"
     call output_scalar(rad_photo,filename)

     filename = trim(adjustl(outdir))//"/opacity_corrupted.dat"
     call output_integer(opacity_corrupted,filename)

     filename = trim(adjustl(outdir))//"/index_lumshell.dat"
     call output_integer(index_lumshell,filename)

     filename = trim(adjustl(outdir))//"/mass_lumshell.dat"
     call output_scalar(mass_lumshell,filename)

     endif ! output_mode .eq. 0

  endif


end subroutine output_all

! *******************************************************************

subroutine output_single_mass(var,filename)
  
  use blmod, only: mass,time
  use parameters

  implicit none
  real*8 var(*)
  character(len=100) filename
  integer nt
  integer i



  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")
  write(666,*) '"Time = ',time

  do i=1,imax
     write(666,"(1P20E29.20E3)") mass(i),var(i)
  enddo
  write(666,*) " "
  write(666,*) " "
  close(666)



end subroutine output_single_mass

! *******************************************************************

subroutine output_single_radius(var,filename)
  
  use blmod, only: r,time
  use parameters

  implicit none
  real*8 var(*)
  character(*) filename
  integer nt
  integer i


  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")
  write(666,*) '"Time = ',time

  do i=1,imax
     write(666,"(1P20E19.10E3)") r(i),var(i)
  enddo
  write(666,*) " "
  write(666,*) " "
  close(666)



end subroutine output_single_radius
    
! *******************************************************************

subroutine output_single_mass_integer(var,filename)

  use blmod, only: mass,time
  use parameters

  implicit none
  integer var(*)
  character(*) filename
  integer nt
  integer i


  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")
  write(666,*) '"Time = ',time

  do i=1,imax
     write(666,"(E19.10E3, I5.4)") mass(i), var(i)
  enddo
  write(666,*) " "
  write(666,*) " "
  close(666)


end subroutine output_single_mass_integer
! ******************************************************************

subroutine output_central(var,filename)
  
  use blmod, only: r,mass,time
  use parameters

  implicit none
  real*8 var(*)
  character(*) filename
  integer nt
  integer i

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")

  write(666,"(1P20E19.10E3)") time,var(1)
  
  close(666)


end subroutine output_central
    
! ******************************************************************

subroutine output_outer(var,filename)
  
  use blmod, only: r,mass,time
  use parameters

  implicit none
  real*8 var(*)
  character(*) filename
  integer nt
  integer i

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                        form='formatted',position="append")

  write(666,"(1P20E19.10E3)") time,var(imax)

  close(666)

end subroutine output_outer


! *******************************************************************
subroutine output_scalar(var,filename)

  use blmod, only: time

  implicit none
  real*8 var
  character(len=100) filename
  integer nt
  integer i

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")

  write(666,"(1P20E19.10E3)") time,var
  
  close(666)

end subroutine output_scalar
    
! *******************************************************************
subroutine output_integer(var,filename)

  use blmod, only: time

  implicit none
  integer var
  character(len=100) filename
  integer nt
  integer i

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")
  
  write(666,"(E19.10E3, I5.4)") time,var

  close(666)

end subroutine output_integer


!******************************************************************************
!output of the variable versus grid point number for a given moment of time
!used to output the initial values of some variables
subroutine output_screenshot(var,filename,imaximum)

  implicit none
  real*8 var(*)
  character(len=100) filename
  integer nt
  integer i
  integer imaximum

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
       form='formatted',position="append")
  
  do i=1, imaximum
      write(666,"(I5.4, E25.16E3)") i,var(i)
  enddo

  close(666)

end subroutine output_screenshot


! *******************************************************************
    subroutine generate_filename(varname,outdir,time,nt,suffix,fname)

      implicit none


      real*8 time
      integer nt
      character(*) varname
      character(len=256) outdir
      character*(*) suffix
      character*(*) fname
      character*(400) aa
      character(len=100) outtime
      character(len=20) outnt
      integer i,ii

      aa=" "
      fname=" "
!      write(aa,"(a32,'_',f10.8,'_nt',i6,'.dat')") varname,time,nt
      write(outnt,"(i10.10)") nt
!      write(aa,"(a32,'_nt',i6)") varname,nt

      fname = trim(adjustl(outdir))//"/"//trim(adjustl(varname))//"_nt_"//outnt
      write(outtime,"(f14.7)") time
!      write(*,*) aa

!      ii=0
!      do i=1,80
!         if(aa(i:i).ne.' ') then
!            ii=ii+1
!            fname(ii:ii)=aa(i:i)
!         endif
!      enddo

      fname = trim(adjustl(fname))//"_time_"//trim(adjustl(outtime))

    end subroutine generate_filename


!******************************************************************************
! write_run_summary: Write a self-contained run metadata file at end of run.
! Called once at program exit. Key-value format for easy machine parsing.
!******************************************************************************
subroutine write_run_summary(nsteps)

  use blmod, only: time, photosphere_fell_on_the_center, &
       diag_ni_eval_total, diag_ni_eval_remap
  use parameters
  use physical_constants, only: msun
  use adaptive_grid_runtime, only: diag_remap_count
  implicit none

  integer, intent(in) :: nsteps
  integer :: u
  character(len=1024) :: filepath

  u = 199
  filepath = trim(adjustl(outdir))//"/run_summary.dat"

  open(unit=u, file=trim(adjustl(filepath)), status='replace', form='formatted')

  write(u,'(A)') '# SuperSNEC run summary'
  write(u,'(A)') '# key = value'
  write(u,'(A)') '#'

  ! Model identity
  write(u,'(A)') '# --- Model identity ---'
  write(u,'(A,A)')     'profile_name = ', trim(adjustl(profile_name))
  write(u,'(A,A)')     'initial_data = ', trim(adjustl(initial_data))

  ! Physical parameters
  write(u,'(A)') '# --- Physical parameters ---'
  write(u,'(A,ES14.6)') 'Ni_mass = ', Ni_mass
  write(u,'(A,F10.6)')  'Ni_mix_fraction = ', Ni_mix_fraction
  write(u,'(A,I4)')     'Ni_mix_kernel = ', Ni_mix_kernel
  if (initial_data .eq. 'Thermal_Bomb') then
     write(u,'(A,ES14.6)') 'final_energy = ', final_energy
     write(u,'(A,ES14.6)') 'bomb_mass_spread = ', bomb_mass_spread
     write(u,'(A,I6)')     'bomb_start_point = ', bomb_start_point
  endif

  ! Grid configuration
  write(u,'(A)') '# --- Grid configuration ---'
  write(u,'(A,I6)')    'imax = ', imax
  write(u,'(A,A)')     'grid_mode = ', trim(adjustl(grid_mode))
  write(u,'(A,F10.4)') 'grid_surface_alpha = ', grid_surface_alpha
  write(u,'(A,I6)')    'grid_adaptive_interval = ', grid_adaptive_interval

  ! Speed parameters
  write(u,'(A)') '# --- Speed parameters ---'
  write(u,'(A,ES14.6)') 'Ni_period = ', Ni_period
  write(u,'(A,ES14.6)') 'Ni_period_max = ', Ni_period_max_target
  write(u,'(A,F10.6)')  'Ni_fractional_change = ', Ni_fractional_change_target
  write(u,'(A,I6)')     'Ni_quad_npoints = ', Ni_quad_npoints
  write(u,'(A,ES10.2)') 'epstol_hydro = ', epstol_hydro
  write(u,'(A,ES10.2)') 'epstol_rad = ', epstol_rad
  write(u,'(A,I6)')     'itmax_hydro = ', itmax_hydro
  write(u,'(A,I6)')     'itmax_rad = ', itmax_rad
  write(u,'(A,ES14.6)') 'dtmax = ', dtmax

  ! Run statistics
  write(u,'(A)') '# --- Run statistics ---'
  write(u,'(A,I10)')    'timesteps = ', nsteps
  write(u,'(A,I10)')    'Ni_evals = ', diag_ni_eval_total
  write(u,'(A,I10)')    'Ni_evals_remap = ', diag_ni_eval_remap
  write(u,'(A,I10)')    'grid_remaps = ', diag_remap_count
  write(u,'(A,ES14.6)') 'tend_actual = ', time
  write(u,'(A,I2)')     'photosphere_fell = ', photosphere_fell_on_the_center

  ! Output configuration
  write(u,'(A)') '# --- Output configuration ---'
  write(u,'(A,I2)')  'output_mode = ', output_mode
  write(u,'(A,I2)')  'eoskey = ', eoskey

  close(u)

end subroutine write_run_summary
