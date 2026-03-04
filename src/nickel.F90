subroutine nickel_heating

  use blmod, only: ye, comp, Ni_heating, Ni_total_luminosity, Ni_number, &
                   r, cr, rho, time, Ni_deposit_function, Ni_energy_rate, delta_mass
  use parameters
  use physical_constants
  implicit none

!local:

  real*8 :: r_Ni !boundary of the region with non-negligible mass fraction of Ni
  real*8 :: delta_r
  real*8 :: r_x, r_max
  real*8 :: r_j, comp_Ni_j
  real*8 :: delta_tau_j

  real*8 :: th, delta_th
  real*8 :: th_min(imax)

  real*8 :: I_prime, I_prime_av

  integer :: i, i_Ni

  integer :: index, index_hint

  ! precomputed gamma-ray opacity per zone (used when ni_raytrace_opt >= 1)
  real*8 :: kappa_gamma(imax)

  ! evaluation point per zone (cell center when ni_ray_interp, boundary otherwise)
  real*8 :: r_eval(imax)

  ! linear interpolation for ni_ray_interp (using cell centers)
  real*8 :: interp_frac
  integer :: idx_lo, idx_hi

  ! parameter constants
  real*8, parameter :: nimin = 1.0d-5
  real*8, parameter :: th_max = pi

!------------------------------------------------------------------------------
!based on the work of Swartz et al., ApJ 446:766 (1995)


  if(Ni_switch.eq.0) then
     !****** heating by Ni is not taken into account *******

     Ni_heating(:) = 0.0d0
     Ni_deposit_function(:) = 0.0d0
     Ni_total_luminosity = 0.0d0

  else
     !*** solve for the local heating due to the radioactive decay of Ni ****

     i_Ni = 0
     do i=imax, 1, -1
        if(comp(i,Ni_number).gt.nimin) then
           r_Ni = r(i)
           i_Ni = i
           exit
        endif
     enddo

     ! Evaluation points default to zone boundaries r(i) (original SNEC)
     r_eval(1:imax) = r(1:imax)

     ! find the limits of integration with respect to the polar angle
     if(i_Ni.eq.0) then
        write(*,*) 'mass fraction of Ni is lower than NIMIN in every grid point'
        write(*,*) 'try reducing NIMIN'
        stop
     else if(i_Ni.eq.1) then
        th_min(i_Ni) = pi*0.5d0
        do i = i_Ni+1, imax
           th_min(i) = acos(-sqrt(r(i)*r(i)-r_Ni*r_Ni)/r(i))
        end do
     else if(i_Ni.eq.imax) then
        th_min(1:i_Ni-1) = 0.0d0
        th_min(i_Ni) = pi*0.5d0
     else
        th_min(1:i_Ni-1) = 0.0d0
        th_min(i_Ni) = pi*0.5d0
        do i = i_Ni+1, imax
           th_min(i) = acos(-sqrt(r(i)*r(i)-r_Ni*r_Ni)/r(i))
        end do
     end if


     ! precompute gamma-ray opacity for each zone (optimized path)
     if (ni_raytrace_opt .ge. 1) then
        do i=1, imax
           kappa_gamma(i) = ye(i) * 0.06d0 * rho(i)
        end do
     endif

     ! find the deposition function at each grid point
     do i=1, imax
        th = th_max
        I_prime_av = 0
        delta_th = (th_max-th_min(i))/Ni_quad_npoints
        index_hint = 1

        do while(th.gt.(th_min(i)+1.d-14))
           r_max = -r(i)*cos(th) + sqrt((r(i)*cos(th))**2-(r(i)**2-r_Ni**2))
           delta_r = r_max/Ni_quad_npoints
           I_prime = 0
           r_x = r_max

           do while(r_x.gt.0)
              r_j = sqrt(r(i)*r(i) + r_x*r_x + 2.0d0*r(i)*r_x*cos(th))

              if(r_j.le.r(1)) then !inside the excised region
                 delta_tau_j = 0.0d0
                 comp_Ni_j = 0.0d0
              else if(r_j.ge.r(imax-1)) then
                 if (ni_raytrace_opt .ge. 1) then
                    delta_tau_j = delta_r * kappa_gamma(imax-1)
                 else
                    delta_tau_j = delta_r * ye(imax-1) * 0.06d0 * rho(imax-1)
                 endif
                 comp_Ni_j = comp(imax-1,Ni_number)
              else
                 if (ni_raytrace_opt .ge. 1) then
                    ! sequential hunt from previous index (replaces bisection)
                    do while (index_hint > 1 .and. r_j .le. r(index_hint))
                       index_hint = index_hint - 1
                    end do
                    do while (index_hint < imax-1 .and. r_j .gt. r(index_hint+1))
                       index_hint = index_hint + 1
                    end do
                    index = index_hint
                    if (ni_ray_interp .ge. 1) then
                       ! Interpolate cell-center quantities between their
                       ! actual positions cr(i), not zone boundaries r(i)
                       if (r_j .lt. cr(index) .and. index .gt. 1) then
                          idx_lo = index - 1
                          idx_hi = index
                       else if (index .lt. imax-1) then
                          idx_lo = index
                          idx_hi = index + 1
                       else
                          idx_lo = index
                          idx_hi = index
                       endif
                       if (idx_lo .eq. idx_hi) then
                          delta_tau_j = delta_r * kappa_gamma(index)
                          comp_Ni_j = comp(index,Ni_number)
                       else
                          interp_frac = (r_j - cr(idx_lo)) / (cr(idx_hi) - cr(idx_lo))
                          interp_frac = max(0.0d0, min(1.0d0, interp_frac))
                          delta_tau_j = delta_r * ((1.0d0-interp_frac)*kappa_gamma(idx_lo) &
                               + interp_frac*kappa_gamma(idx_hi))
                          comp_Ni_j = (1.0d0-interp_frac)*comp(idx_lo,Ni_number) &
                               + interp_frac*comp(idx_hi,Ni_number)
                       endif
                    else
                       delta_tau_j = delta_r * kappa_gamma(index)
                       comp_Ni_j = comp(index,Ni_number)
                    endif
                 else
                    ! original SNEC: linear search + inline kappa computation
                    index = 1
                    do while (index < imax-1 .and. r_j .gt. r(index+1))
                       index = index + 1
                    end do
                    if (ni_ray_interp .ge. 1) then
                       if (r_j .lt. cr(index) .and. index .gt. 1) then
                          idx_lo = index - 1
                          idx_hi = index
                       else if (index .lt. imax-1) then
                          idx_lo = index
                          idx_hi = index + 1
                       else
                          idx_lo = index
                          idx_hi = index
                       endif
                       if (idx_lo .eq. idx_hi) then
                          delta_tau_j = delta_r * ye(index)*0.06d0*rho(index)
                          comp_Ni_j = comp(index,Ni_number)
                       else
                          interp_frac = (r_j - cr(idx_lo)) / (cr(idx_hi) - cr(idx_lo))
                          interp_frac = max(0.0d0, min(1.0d0, interp_frac))
                          delta_tau_j = delta_r * ( &
                               (1.0d0-interp_frac)*ye(idx_lo)*0.06d0*rho(idx_lo) &
                               + interp_frac*ye(idx_hi)*0.06d0*rho(idx_hi))
                          comp_Ni_j = (1.0d0-interp_frac)*comp(idx_lo,Ni_number) &
                               + interp_frac*comp(idx_hi,Ni_number)
                       endif
                    else
                       delta_tau_j = delta_r * ye(index) * 0.06d0 * rho(index)
                       comp_Ni_j = comp(index,Ni_number)
                    endif
                 endif
              end if

              I_prime = (I_prime-comp_Ni_j)*exp(-delta_tau_j) + comp_Ni_j
              r_x = r_x - delta_r
           end do

           I_prime_av = I_prime_av + I_prime*sin(th)*delta_th*0.5d0
           th = th - delta_th
        end do

        Ni_deposit_function(i) = I_prime_av
     end do

     ! rate of energy release per gram of radioactive material
     Ni_energy_rate = &
          3.24d10*exp(-time*overtau_Ni) + 7.29d9*exp(-time*overtau_Co)

     ! local rate of gamma-ray energy deposition at a given grid point
     Ni_heating(1:imax) = Ni_energy_rate*Ni_deposit_function(1:imax)

     ! total energy per second, deposited to the model by gamma-rays
     Ni_total_luminosity = &
          Ni_energy_rate*sum(Ni_deposit_function(1:imax)*delta_mass(1:imax))

  end if

end subroutine nickel_heating
