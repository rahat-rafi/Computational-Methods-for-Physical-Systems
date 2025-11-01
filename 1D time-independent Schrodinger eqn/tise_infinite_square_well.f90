program schrodinger_infinite_well
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 300)

    ! Physical constants and parameters
    real(dp), parameter :: hbar = 1.0_dp, mass = 1.0_dp, L = 1.0_dp
    integer, parameter  :: n_points = 1001   ! any number >= 2 works; uniform grid assumed
    integer, parameter  :: max_states = 5, n_vars = 2
    real(dp), parameter :: tol_mismatch = 1.0e-10_dp, tol_kspan = 1.0e-12_dp

    ! Variables
    real(dp) :: k, k_low, k_high, mismatch, mismatch_low, mismatch_high, energy
    real(dp) :: y(n_vars), x(n_points)
    real(dp) :: analytic_k, analytic_energy, abs_error_k, rel_error_k, abs_error_E, rel_error_E
    integer :: n_state, i, iter
    logical :: converged

    ! Initialize grid
    do i = 1, n_points
        x(i) = (i-1) * L / (n_points-1)
    end do

    print *, '=================================================='
    print *, 'Infinite Square Well - Shooting Method'
    print *, 'Using your RK4N and DERIVS functions'
    print *, '=================================================='
    print *, 'Well width L =', L
    print *, 'Grid points  =', n_points, ' (uniform spacing required)'
    print *, 'Searching for first', max_states, 'eigenstates'
    print *, ''

    do n_state = 1, max_states
        ! Analytic-centered bracket
        analytic_k = n_state * 3.141592653589793_dp / L
        k_low  = max(0.1_dp, analytic_k - 1.0_dp)
        k_high = analytic_k + 1.0_dp

        call shooting_method(k_low,  L, n_points, y, mismatch_low)
        call shooting_method(k_high, L, n_points, y, mismatch_high)

        ! Expand bracket if needed
        do i = 1, 40
            if (mismatch_low * mismatch_high < 0.0_dp) exit
            k_low  = max(0.1_dp, k_low  - 0.5_dp)
            k_high = k_high + 0.5_dp
            call shooting_method(k_low,  L, n_points, y, mismatch_low)
            call shooting_method(k_high, L, n_points, y, mismatch_high)
        end do
        if (mismatch_low * mismatch_high > 0.0_dp) then
            print *, '  Failed to bracket a root for state n =', n_state
            cycle
        end if

        print '(A,I0,A)', 'Searching for state n = ', n_state, ' ...'

        ! Bisection
        converged = .false.
        do iter = 1, 200
            k = 0.5_dp * (k_low + k_high)
            call shooting_method(k, L, n_points, y, mismatch)

            if (abs(mismatch) < tol_mismatch) then
                converged = .true.
                exit
            end if

            if (mismatch_low * mismatch < 0.0_dp) then
                k_high = k
                mismatch_high = mismatch
            else
                k_low = k
                mismatch_low = mismatch
            end if

            if (k_high - k_low < tol_kspan) then
                converged = .true.
                k = 0.5_dp*(k_low+k_high)
                exit
            end if
        end do

        if (.not. converged) then
            print *, '  Failed to converge for state n =', n_state
            cycle
        end if

        ! Energies and error report
        energy = (hbar**2 * k**2) / (2.0_dp * mass)
        analytic_k = n_state * 3.141592653589793_dp / L
        analytic_energy = (hbar**2 * analytic_k**2) / (2.0_dp * mass)
        abs_error_k = abs(k - analytic_k)
        rel_error_k = abs_error_k / analytic_k
        abs_error_E = abs(energy - analytic_energy)
        rel_error_E = abs_error_E / analytic_energy

        print '(A,F16.12)', '  Found k = ', k
        print '(A,F16.12)', '  Energy E = ', energy
        print '(A,F16.12,A)', '  Analytic k = ', analytic_k, ' (exact)'
        print '(A,ES10.2)', '  |k - k_exact|  = ', abs_error_k
        print '(A,ES10.2)', '  rel. error k   = ', rel_error_k
        print '(A,ES10.2)', '  |E - E_exact|  = ', abs_error_E
        print '(A,ES10.2)', '  rel. error E   = ', rel_error_E
        print '(A,I3,A)',   '  Converged in ', iter, ' iterations'
        print *, ''

        ! Save wavefunction (normalized) with columns & per-point errors
        call save_wavefunction_with_error(k, energy, analytic_k, analytic_energy, &
             abs_error_k, rel_error_k, abs_error_E, rel_error_E, x, n_points, n_state)
    end do

    print *, 'Potential saved to: isw_potential.dat'
    print *, 'Wavefunctions saved to: isw_n*.dat'

contains

! ------------------------------------------------------------
! Potential (0 inside; very large at walls)
! ------------------------------------------------------------
real(dp) function infinite_well_potential(xpos, L)
    implicit none
    real(dp), intent(in) :: xpos, L
    if (xpos <= 0.0_dp .or. xpos >= L) then
        infinite_well_potential = 1.0e20_dp
    else
        infinite_well_potential = 0.0_dp
    end if
end function infinite_well_potential

! ------------------------------------------------------------
! Shooting method: integrate y'' = -k^2 y from 0 to L
! Mismatch is psi(L) (psi(0)=0, psi'(0)=1)
! ------------------------------------------------------------
subroutine shooting_method(k, L, n_steps, y, mismatch)
    implicit none
    real(dp), intent(in)  :: k, L
    integer, intent(in)   :: n_steps
    real(dp), intent(out) :: y(2), mismatch
    real(dp) :: dx
    real(dp) :: err(2)
    integer :: i

    dx = L / (n_steps - 1)

    y(1) = 0.0_dp   ! psi(0)
    y(2) = 1.0_dp   ! psi'(0)

    do i = 1, n_steps-1
        call rk4n(dx, y, 2, err, k)
    end do

    mismatch = y(1)   ! want psi(L)=0
end subroutine shooting_method

subroutine rk4n(dx, y, n, err, k)
    implicit none
    real(dp), intent(in) :: dx
    integer, intent(in) :: n
    real(dp), intent(inout) :: y(n)
    real(dp), intent(out) :: err(n)
    real(dp), intent(in) :: k
    
    real(dp) :: y1(n), y2(n), y3(n)
    real(dp) :: k1(n), k2(n), k3(n), k4(n)
    real(dp) :: rksum
    integer :: ivar
    
    call derivs(y, k1, n, k)
    k1 = k1 * dx
    y1 = y + 0.5_dp * k1
    
    call derivs(y1, k2, n, k)
    k2 = k2 * dx
    y2 = y + 0.5_dp * k2
    
    call derivs(y2, k3, n, k)
    k3 = k3 * dx
    y3 = y + k3
    
    call derivs(y3, k4, n, k)
    k4 = k4 * dx
    
    do ivar = 1, n
        rksum = (1.0_dp/6.0_dp) * (k1(ivar) + 2.0_dp*k2(ivar) + 2.0_dp*k3(ivar) + k4(ivar))
        err(ivar) = abs(rksum - k2(ivar))      ! embedded midpoint vs RK4 step error
        y(ivar)   = y(ivar) + rksum
    end do
end subroutine rk4n

subroutine derivs(y, dydt, n, k)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: y(n)
    real(dp), intent(out) :: dydt(n)
    real(dp), intent(in) :: k
    
    dydt(1) = y(2)          ! psi'
    dydt(2) = -k**2 * y(1)  ! psi''
end subroutine derivs

real(dp) function integrate_simpson_any(x, y) result(intval)
    implicit none
    real(dp), intent(in) :: x(:), y(:)
    integer :: n, i, m
    real(dp) :: h, s_odd, s_even, trap

    n = size(x)
    if (size(y) /= n) error stop 'integrate_simpson_any: x and y must be same length'
    if (n < 2)        error stop 'integrate_simpson_any: need at least 2 points'

    ! uniform spacing check
    h = x(2) - x(1)
    do i = 3, n
        if (abs((x(i)-x(i-1)) - h) > 1.0e-12_dp) error stop 'integrate_simpson_any: x not uniform'
    end do

    m = n - 1   ! number of intervals

    if (m == 1) then
        ! Just trapezoid for 2 points
        intval = 0.5_dp*h*(y(1) + y(2))
        return
    end if

    if (mod(m,2) == 0) then
        ! Even number of intervals: pure composite Simpson on all
        s_odd  = 0.0_dp
        s_even = 0.0_dp
        do i = 2, n-1, 2
            s_odd = s_odd + y(i)
        end do
        do i = 3, n-2, 2
            s_even = s_even + y(i)
        end do
        intval = (h/3.0_dp)*( y(1) + 4.0_dp*s_odd + 2.0_dp*s_even + y(n) )
    else
        ! Odd number of intervals: Simpson on first (m-1) intervals, trap on last
        s_odd  = 0.0_dp
        s_even = 0.0_dp
        do i = 2, n-2, 2
            s_odd = s_odd + y(i)
        end do
        do i = 3, n-3, 2
            s_even = s_even + y(i)
        end do
        intval = (h/3.0_dp)*( y(1) + 4.0_dp*s_odd + 2.0_dp*s_even + y(n-1) )
        trap   = 0.5_dp*h*( y(n-1) + y(n) )
        intval = intval + trap
    end if
end function integrate_simpson_any

subroutine normalize_wavefunction(x, psi)
    implicit none
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: psi(:)
    real(dp), allocatable :: f(:)
    real(dp) :: normint, norm
    integer :: n

    n = size(x)
    allocate(f(n))
    f = psi*psi

    normint = integrate_simpson_any(x, f)
    if (normint <= 0.0_dp) then
        deallocate(f)
        error stop 'normalize: nonpositive integral'
    end if

    norm = 1.0_dp / sqrt(normint)
    psi  = psi * norm
    deallocate(f)
end subroutine normalize_wavefunction

subroutine save_wavefunction_with_error(k, energy, analytic_k, analytic_energy, &
                           abs_error_k, rel_error_k, abs_error_E, rel_error_E, &
                           x, n_points, n_state)
    implicit none
    real(dp), intent(in) :: k, energy, analytic_k, analytic_energy
    real(dp), intent(in) :: abs_error_k, rel_error_k, abs_error_E, rel_error_E
    integer, intent(in)  :: n_points, n_state
    real(dp), intent(in) :: x(n_points)

    character(len=64) :: filename
    real(dp) :: y(2), step_e(2), dx
    real(dp), allocatable :: psi(:), V(:), step_err(:), psi_exact(:), err_point(:)
    integer :: i, file_unit
    integer :: nfac
    real(dp), parameter :: pi = 3.141592653589793_dp

    dx = x(2) - x(1)
    allocate(psi(n_points), V(n_points), step_err(n_points), psi_exact(n_points), err_point(n_points))

    ! Integrate IVP to get psi (unnormalized) and per-step embedded error
    y(1) = 0.0_dp  ! psi(0)
    y(2) = 1.0_dp  ! psi'(0)
    psi(1) = y(1)
    V(1)   = infinite_well_potential(x(1), L)
    step_err(1) = 0.0_dp

    do i = 1, n_points-1
        call rk4n(dx, y, 2, step_e, k)
        psi(i+1)      = y(1)
        V(i+1)        = infinite_well_potential(x(i+1), L)
        step_err(i+1) = step_e(1)   ! embedded step error on psi
    end do

    ! Normalize numerical psi
    call normalize_wavefunction(x, psi)

    ! Exact normalized eigenfunction for this state: sqrt(2/L) * sin(n*pi*x/L)
    ! Determine n from analytic_k = n*pi/L  => n = round(analytic_k * L / pi)
    nfac = nint(analytic_k * L / pi)
    do i = 1, n_points
        psi_exact(i) = sqrt(2.0_dp/L) * sin(nfac * pi * x(i) / L)
        err_point(i) = psi(i) - psi_exact(i)
    end do

    write(filename, '(A,I0,A)') 'isw_n', n_state, '.dat'
    open(newunit=file_unit, file=trim(filename), status='replace')
    write(file_unit, '(A,I0)')         '# Infinite Square Well - State n = ', n_state
    write(file_unit, '(A,F20.12)')     '# Numerical k        = ', k
    write(file_unit, '(A,F20.12)')     '# Analytic  k        = ', analytic_k
    write(file_unit, '(A,ES15.6)')     '# |k - k_exact|      = ', abs_error_k
    write(file_unit, '(A,ES15.6)')     '# rel. error k       = ', rel_error_k
    write(file_unit, '(A,F20.12)')     '# Numerical Energy   = ', energy
    write(file_unit, '(A,F20.12)')     '# Analytic  Energy   = ', analytic_energy
    write(file_unit, '(A,ES15.6)')     '# |E - E_exact|      = ', abs_error_E
    write(file_unit, '(A,ES15.6)')     '# rel. error E       = ', rel_error_E
    write(file_unit, '(A,F10.6)')      '# Well width L       = ', L
    write(file_unit, '(A)')            '# columns: x   V(x)   psi_num   psi_exact   psi_num^2   (psi_num-psi_exact)   step_err_psi'

    do i = 1, n_points
        write(file_unit, '(7ES20.10)') x(i), V(i), psi(i), psi_exact(i), psi(i)*psi(i), err_point(i), step_err(i)
    end do

    close(file_unit)
    deallocate(psi, V, step_err, psi_exact, err_point)
end subroutine save_wavefunction_with_error

end program schrodinger_infinite_well
