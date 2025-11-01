program sho_write_pw_n0_7
  implicit none
  integer,  parameter :: dp = selected_real_kind(15,300)

  ! ===== Physics & nondimensionalization =====
  ! ξ = x/a,  v(ξ)=ξ^2
  ! γ = (m ω a^2)/ħ;  V0 = ½ m ω^2 a^2;  ε = E/V0;  ε_n^(exact) = (2n+1)/γ
  real(dp), parameter :: hbar  = 1.0_dp
  real(dp), parameter :: mass  = 1.0_dp
  real(dp), parameter :: omega = 1.0_dp
  real(dp), parameter :: a     = 1.0_dp
  real(dp), parameter :: V0    = 0.5_dp*mass*omega*omega*a*a
  real(dp), parameter :: gamma = (mass*omega*a*a)/hbar

  ! ===== Grid / numerics =====
  integer,  parameter :: n_points = 2001          ! full grid size on [-Xmax, Xmax]
  real(dp), parameter :: Xmax     = 8.0_dp
  real(dp), parameter :: TolE     = 1.0e-10_dp    ! bisection tol for epsilon
  real(dp), parameter :: PAD0     = 0.6_dp        ! initial +/- bracket around guess
  real(dp), parameter :: PADG     = 1.6_dp        ! geometric expand factor
  real(dp), parameter :: EMAX     = 200.0_dp      ! safety cap

  real(dp) :: dx, xi(n_points), vgrid(n_points)
  integer  :: i, n, i0

  ! visible to derivs()
  real(dp) :: eps_current

  ! ----- grid & potential -----
  dx = (2.0_dp*Xmax)/real(n_points-1,dp)
  do i=1,n_points
    xi(i)    = -Xmax + dx*real(i-1,dp)
    vgrid(i) = xi(i)*xi(i)
  end do
  i0 = nearest_index(xi, 0.0_dp)

  print '(A)',       'Quantum HO (dimensionless): writing pw_n0..pw_n7'
  print '(A,F10.6)', 'gamma = ', gamma
  print '(A,F10.6)', 'Xmax  = ', Xmax
  print '(A,I0)',    'grid  = ', n_points

  do n = 0, 7
    call build_and_write_one(n, xi, vgrid, dx)
  end do

  print *, 'Done: pw_n0.dat ... pw_n7.dat'

contains
  subroutine build_and_write_one(nsho, xi, vgrid, dx)
    implicit none
    integer,  intent(in) :: nsho
    real(dp), intent(in) :: xi(:), vgrid(:), dx

    integer :: N, i0, j, u
    logical :: evenp
    real(dp) :: eps, eps_exact, En_phys
    real(dp), allocatable :: psi(:), psi_ex(:), errp(:), steperr(:)
    character(len=64) :: fname
    real(dp) :: pi, gqrt, Hn, normfac_ex

    N     = size(xi)
    i0    = nearest_index(xi, 0.0_dp)
    evenp = (mod(nsho,2) == 0)

    ! ---- Find epsilon numerically (parity mismatch at ξ=0) ----
    call bracket_and_bisect_level(nsho, evenp, eps)

    ! ---- Build ψ on full grid: -Xmax -> 0, then 0 -> +Xmax with parity ICs ----
    allocate(psi(N), psi_ex(N), errp(N), steperr(N))
    call build_full_mode(eps, evenp, xi, dx, psi, steperr)

    ! ---- Normalize numerically (Simpson) ----
    call normalize_simpson(xi, psi)

    ! ---- Analytic ψ_n(ξ) for comparison: ψ_n = (γ/π)^{1/4}/√(2^n n!) e^{-γ ξ^2/2} H_n(√γ ξ)
    pi   = 4.0_dp*atan(1.0_dp)
    gqrt = sqrt(gamma)
    normfac_ex = (gamma**0.25_dp) / (pi**0.25_dp * sqrt( two_pow(nsho) * real(fact_int(nsho),dp) ))
    do j=1,N
      call hermiteH(nsho, gqrt*xi(j), Hn)
      psi_ex(j) = normfac_ex * exp(-0.5_dp*(gqrt*xi(j))**2) * Hn
    end do
    errp = psi - psi_ex

    eps_exact = (2.0_dp*real(nsho,dp) + 1.0_dp)/gamma
    En_phys   = hbar*omega*(real(nsho,dp)+0.5_dp)

    ! ---- write output: pw_n<n>.dat (10-line header) ----
    write(fname,'(A,I0,A)') 'pw_n', nsho, '.dat'
    open(newunit=u, file=trim(fname), status='replace')
    write(u,'(A,I0)')         '# Harmonic Oscillator - state n = ', nsho
    write(u,'(A,I0)')         '# Parity (even=+1, odd=-1) = ', merge(1,-1,evenp)
    write(u,'(A,ES15.7)')     '# gamma = ', gamma
    write(u,'(A,ES15.7)')     '# epsilon (numeric) = ', eps
    write(u,'(A,ES15.7)')     '# epsilon (analytic)= ', eps_exact
    write(u,'(A,ES15.7)')     '# E_n (physical)    = ', En_phys
    write(u,'(A,F10.6)')      '# Xmax = ', Xmax
    write(u,'(A,I0)')         '# grid = ', N
    write(u,'(A)')            '# columns: xi   v(xi)   psi_num   psi_exact   psi_num^2   (psi_num-psi_exact)   step_err_psi'
    do j=1,N
      write(u,'(7ES20.10)') xi(j), vgrid(j), psi(j), psi_ex(j), psi(j)*psi(j), errp(j), steperr(j)
    end do
    close(u)

    deallocate(psi, psi_ex, errp, steperr)
  end subroutine build_and_write_one

  subroutine bracket_and_bisect_level(nsho, evenp, eps_root)
    integer, intent(in) :: nsho
    logical, intent(in) :: evenp
    real(dp), intent(out) :: eps_root
    real(dp) :: guess, pad, a, b, Fa, Fb, mid, Fmid
    integer :: it

    guess = (2.0_dp*real(nsho,dp) + 1.0_dp)/gamma
    pad   = max(PAD0, 0.25_dp*guess)
    a     = max(1.0e-12_dp, guess - pad)
    b     = guess + pad

    call parity_mismatch(a, evenp, Fa)
    call parity_mismatch(b, evenp, Fb)

    it=0
    do while (Fa*Fb > 0.0_dp .and. b < EMAX)
      pad = pad*PADG
      a = max(1.0e-12_dp, guess - pad)
      b = guess + pad
      call parity_mismatch(a, evenp, Fa)
      call parity_mismatch(b, evenp, Fb)
      it = it+1; if (it>80) exit
    end do
    if (Fa*Fb > 0.0_dp) then
      eps_root = guess      ! fallback; should rarely happen
      return
    end if

    do
      if (abs(b-a) <= TolE) exit
      mid = 0.5_dp*(a+b)
      call parity_mismatch(mid, evenp, Fmid)
      if (Fa*Fmid <= 0.0_dp) then
        b  = mid; Fb = Fmid
      else
        a  = mid; Fa = Fmid
      end if
    end do
    eps_root = 0.5_dp*(a+b)
  end subroutine bracket_and_bisect_level

  ! march from ξ=-Xmax to 0 with arbitrary slope; mismatch = ψ(0) (odd) or ψ'(0) (even)
  subroutine parity_mismatch(eps, evenp, M)
    real(dp), intent(in)  :: eps
    logical, intent(in)   :: evenp
    real(dp), intent(out) :: M
    real(dp) :: y(3), e(3), x
    integer  :: i, i0loc

    eps_current = eps
    y = [0.0_dp, 1.0_dp, -Xmax]     ! ψ(-Xmax)=0, ψ' arbitrary, ξ=-Xmax
    i0loc = nearest_index(xi, 0.0_dp)
    x = y(3)
    do i = 2, i0loc
      call rk4n(dx, y, 3, e, gamma) ! step to the right by +dx
    end do
    if (evenp) then
      M = y(2)                      ! want ψ'(0)=0
    else
      M = y(1)                      ! want ψ(0)=0
    end if
  end subroutine parity_mismatch

  subroutine build_full_mode(eps, evenp, x, dx, psi, step_err)
    real(dp), intent(in)    :: eps, x(:), dx
    logical,  intent(in)    :: evenp
    real(dp), intent(out)   :: psi(:), step_err(:)

    integer :: N, i, i0loc
    real(dp) :: y(3), e(3), psi0, dpsi0

    eps_current = eps
    N = size(x)
    i0loc = nearest_index(x, 0.0_dp)
    psi = 0.0_dp; step_err = 0.0_dp

    ! left: -Xmax -> 0
    y = [0.0_dp, 1.0_dp, -Xmax]
    psi(1) = y(1)
    do i = 2, i0loc
      call rk4n(dx, y, 3, e, gamma)
      psi(i)      = y(1)
      step_err(i) = e(1)
    end do
    psi0  = y(1)
    dpsi0 = y(2)

    ! right: 0 -> +Xmax with parity ICs
    if (evenp) then
      y = [psi0, 0.0_dp, 0.0_dp]
    else
      y = [0.0_dp, dpsi0, 0.0_dp]
    end if
    psi(i0loc) = y(1)
    do i = i0loc+1, N
      call rk4n(dx, y, 3, e, gamma)
      psi(i)      = y(1)
      step_err(i) = e(1)
    end do
  end subroutine build_full_mode

  subroutine rk4n(dx, y, n, err, k)
    real(dp), intent(in)    :: dx
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: y(n)
    real(dp), intent(out)   :: err(n)
    real(dp), intent(in)    :: k       ! here: k = gamma
    real(dp) :: y1(n), y2(n), y3(n)
    real(dp) :: k1(n), k2(n), k3(n), k4(n)
    real(dp) :: rksum
    integer  :: iv

    call derivs(y, k1, n, k);  k1 = k1*dx;  y1 = y + 0.5_dp*k1
    call derivs(y1, k2, n, k); k2 = k2*dx;  y2 = y + 0.5_dp*k2
    call derivs(y2, k3, n, k); k3 = k3*dx;  y3 = y + k3
    call derivs(y3, k4, n, k); k4 = k4*dx

    do iv=1,n
      rksum     = (k1(iv) + 2.0_dp*k2(iv) + 2.0_dp*k3(iv) + k4(iv))/6.0_dp
      err(iv)   = abs(rksum - k2(iv))
      y(iv)     = y(iv) + rksum
    end do
  end subroutine rk4n

  ! State vector: y(1)=ψ, y(2)=ψ', y(3)=ξ ; ξ' = 1 ; ψ'' = γ^2(ξ^2 − ε) ψ
  subroutine derivs(y, dydt, n, k)
    integer,  intent(in) :: n
    real(dp), intent(in) :: y(n)
    real(dp), intent(out):: dydt(n)
    real(dp), intent(in) :: k     ! k = gamma
    dydt(1) = y(2)
    dydt(2) = (k*k) * ( y(3)*y(3) - eps_current ) * y(1)
    dydt(3) = 1.0_dp
  end subroutine derivs

  ! Composite Simpson over a uniform grid
  real(dp) function simpson(x, f) result(S)
    real(dp), intent(in) :: x(:), f(:)
    integer :: n, i, m
    real(dp) :: h, so, se, trap
    n = size(x); if (size(f)/=n) error stop 'simpson: size mismatch'
    if (n<2) error stop 'simpson: need >=2 points'
    h = x(2)-x(1)
    do i=3,n
      if (abs((x(i)-x(i-1))-h) > 1.0e-12_dp) error stop 'simpson: nonuniform grid'
    end do
    m = n-1
    if (m==1) then
      S = 0.5_dp*h*(f(1)+f(2)); return
    end if
    if (mod(m,2)==0) then
      so=0; se=0
      do i=2,n-1,2; so=so+f(i); end do
      do i=3,n-2,2; se=se+f(i); end do
      S = (h/3.0_dp)*(f(1)+4.0_dp*so+2.0_dp*se+f(n))
    else
      so=0; se=0
      do i=2,n-2,2; so=so+f(i); end do
      do i=3,n-3,2; se=se+f(i); end do
      S = (h/3.0_dp)*(f(1)+4.0_dp*so+2.0_dp*se+f(n-1))
      trap = 0.5_dp*h*(f(n-1)+f(n))
      S = S + trap
    end if
  end function simpson

  subroutine normalize_simpson(x, psi)
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: psi(:)
    real(dp), allocatable   :: f(:)
    real(dp) :: I
    allocate(f(size(x))); f = psi*psi
    I = simpson(x, f)
    if (I <= 0.0_dp) then
      deallocate(f); error stop 'normalize: nonpositive integral'
    end if
    psi = psi / sqrt(I)
    deallocate(f)
  end subroutine normalize_simpson

  ! Physicists' Hermite H_n(z) by three-term recurrence
  subroutine hermiteH(n, z, Hn)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: z
    real(dp), intent(out) :: Hn
    integer :: k
    real(dp) :: Hkm1, Hk, Hkp1
    if (n==0) then
      Hn = 1.0_dp
    else if (n==1) then
      Hn = 2.0_dp*z
    else
      Hkm1 = 1.0_dp
      Hk   = 2.0_dp*z
      do k=1,n-1
        Hkp1 = 2.0_dp*z*Hk - 2.0_dp*real(k,dp)*Hkm1
        Hkm1 = Hk
        Hk   = Hkp1
      end do
      Hn = Hk
    end if
  end subroutine hermiteH

  pure real(dp) function two_pow(n) result(v)
    integer, intent(in) :: n
    integer :: i
    v = 1.0_dp
    do i=1,n
      v = 2.0_dp*v
    end do
  end function two_pow

  pure integer function fact_int(n) result(fv)
    integer, intent(in) :: n
    integer :: k
    fv = 1
    if (n>1) then
      do k=2,n
        fv = fv*k
      end do
    end if
  end function fact_int

  integer function nearest_index(x, x0) result(idx)
    real(dp), intent(in) :: x(:), x0
    integer :: j, N
    real(dp) :: best, d
    N = size(x); idx = 1; best = abs(x(1)-x0)
    do j=2,N
      d = abs(x(j)-x0)
      if (d < best) then
        best = d; idx = j
      end if
    end do
  end function nearest_index

end program sho_write_pw_n0_7

