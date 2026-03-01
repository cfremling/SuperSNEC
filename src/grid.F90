subroutine grid

  use blmod, only: mass, cmass, delta_mass, delta_cmass
  use parameters, only: imax, grid_mode, grid_pattern_file
  implicit none

  character(len=32) :: mode
  character(len=256) :: pattern_path
  real*8 :: mass_start, mass_end

  if (imax.le.1) return

  mode = normalize_mode(grid_mode)
  pattern_path = trim(adjustl(grid_pattern_file))
  if (len_trim(pattern_path).eq.0) pattern_path = 'tables/GridPattern.dat'

  mass_start = mass(1)
  mass_end = mass(imax)

  if (mode .eq. 'legacy_pattern') then
     if (.not.build_legacy_pattern_grid(pattern_path, mass_start, mass_end)) then
        call legacy_pattern_failure('Unable to build legacy pattern grid.')
     endif
  else if (mode .eq. 'adaptive_runtime') then
     call build_surface_concentrated_grid(mass_start, mass_end)
  else
     write(*,*) '******* Unsupported grid_mode: ', trim(grid_mode)
     write(*,*) '******* Supported modes: adaptive_runtime, legacy_pattern'
     stop 1
  endif

  call finalize_differences()

contains

  function normalize_mode(raw) result(clean)
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

  subroutine build_surface_concentrated_grid(mstart, mend)
    ! Piecewise grid: uniform inner half + geometric outer half.
    ! The per-step ratio is derived from grid_surface_alpha:
    !   r_eff = exp(-alpha / n_outer)
    ! This keeps the total surface concentration bounded (~exp(alpha))
    ! regardless of imax, ensuring proper convergence.
    ! Minimum cell width is enforced by clamping r_eff from below.
    use parameters, only: grid_surface_alpha, grid_min_cell_frac
    real*8, intent(in) :: mstart, mend
    real*8 :: r_eff, r_min, alpha, span, dm_inner, outer_sum, rpower
    integer :: i, n_transition, n_inner, n_outer

    span = mend - mstart
    n_transition = imax / 2
    n_inner = n_transition - 1
    n_outer = imax - n_transition

    alpha = max(grid_surface_alpha, 0.1d0)
    if (n_outer .gt. 0) then
       r_eff = exp(-alpha / dble(n_outer))
       ! Enforce minimum cell width: thinnest cell = dm_inner * r^(n-1).
       ! We need r^(n-1) >= grid_min_cell_frac, so r >= frac^(1/(n-1)).
       if (grid_min_cell_frac .gt. 0.0d0 .and. n_outer .gt. 1) then
          r_min = grid_min_cell_frac ** (1.0d0 / dble(n_outer - 1))
          r_eff = max(r_eff, r_min)
       endif
       r_eff = max(0.01d0, min(0.9999d0, r_eff))
    else
       r_eff = 1.0d0
    endif

    outer_sum = (1.0d0 - r_eff**n_outer) / (1.0d0 - r_eff)
    dm_inner = span / (dble(n_inner) + outer_sum)

    mass(1) = mstart
    do i=2,n_transition
       mass(i) = mass(i-1) + dm_inner
    enddo
    rpower = 1.0d0
    do i=n_transition+1,imax
       mass(i) = mass(i-1) + dm_inner * rpower
       rpower = rpower * r_eff
    enddo
    mass(imax) = mend
  end subroutine build_surface_concentrated_grid

  logical function build_legacy_pattern_grid(path, mstart, mend)
    character(len=*), intent(in) :: path
    real*8, intent(in) :: mstart, mend
    real*8, allocatable :: pattern(:)
    character(len=512) :: line
    integer :: i, ios, lines
    integer, parameter :: unit = 666
    logical :: exists

    build_legacy_pattern_grid = .false.

    if (imax .le. 1) then
       if (imax .eq. 1) mass(1) = mstart
       return
    endif

    inquire(file=path, exist=exists)
    if (.not.exists) then
       call legacy_pattern_failure('Unable to read pattern file: '//trim(path))
    endif

    open(unit=unit, file=path, status='unknown', action='read', iostat=ios)
    if (ios .ne. 0) then
       call legacy_pattern_failure('Unable to read pattern file: '//trim(path))
    endif

    lines = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios .ne. 0) exit
       line = adjustl(line)
       if (len_trim(line).eq.0) cycle
       if (line(1:1) .eq. '#') cycle
       lines = lines + 1
    enddo
    close(unit)

    if (lines .ne. imax) then
       call legacy_pattern_count_error(trim(path), lines)
    endif

    allocate(pattern(imax))
    open(unit=unit, file=path, status='unknown', action='read', iostat=ios)
    if (ios .ne. 0) then
        deallocate(pattern)
        call legacy_pattern_failure('Unable to read pattern file: '//trim(path))
    endif

    i = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios .ne. 0) exit
       line = adjustl(line)
       if (len_trim(line).eq.0) cycle
       if (line(1:1) .eq. '#') cycle
       i = i + 1
       if (i .gt. imax) then
          close(unit)
          deallocate(pattern)
          call legacy_pattern_count_error(trim(path), i)
       endif
       read(line, *) pattern(i)
    enddo
    close(unit)

    if (i .ne. imax) then
       deallocate(pattern)
       call legacy_pattern_count_error(trim(path), i)
    endif

    mass(1) = mstart
    do i=2,imax
       mass(i) = mass(i-1) + (pattern(i) - pattern(i-1))*(mend - mstart)
    enddo

    deallocate(pattern)
    build_legacy_pattern_grid = .true.
  end function build_legacy_pattern_grid

  subroutine legacy_pattern_failure(message)
    character(len=*), intent(in) :: message
    write(*,*) '******* '//trim(message)
    stop 1
  end subroutine legacy_pattern_failure

  subroutine legacy_pattern_count_error(path, count)
    character(len=*), intent(in) :: path
    integer, intent(in) :: count

    write(*,*) '******* Number of lines in the file '//trim(path)
    write(*,*) '******* does not coincide with the number of grid points.'
    write(*,*) '******* Lines in file: ', count
    write(*,*) '******* Requested grid points: ', imax
    write(*,*) '******* Please, adjust one of the two.'
    stop 1
  end subroutine legacy_pattern_count_error

  subroutine finalize_differences()
    integer :: i

    do i=1,imax-1
       cmass(i) = mass(i) + 0.5d0*(mass(i+1) - mass(i))
    enddo
    cmass(imax) = mass(imax) + 0.5d0*(mass(imax) - mass(imax-1))

    do i=1,imax-1
       delta_mass(i) = mass(i+1) - mass(i)
       delta_cmass(i) = cmass(i+1) - cmass(i)
    enddo
    delta_mass(imax) = delta_mass(imax-1)
    delta_cmass(imax) = delta_cmass(imax-1)
  end subroutine finalize_differences

end subroutine grid
