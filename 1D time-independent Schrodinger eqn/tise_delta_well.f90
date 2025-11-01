program tise_delta_well
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)

  ! ===== Dimensionless setup =====
  ! d^2 psi / d xi^2 = gamma^2 * ( v(xi) - eps ) * psi
  ! delta well: v(xi) = -lambda * delta(xi)
  ! Jump: psi'(0^+) - psi'(0^-) = - gamma^2 * lambda * psi(0)
  !
  ! Bound state:
  !   eps_b   = - (gamma^2 * lambda^2) / 4
  !   kappa   = gamma * sqrt(-eps_b) = (gamma^2 * lambda) / 2

  ! Physical scale choices (arbitrary here; only gamma, lambda matter)
  real(dp), parameter :: hbar  = 1.0_dp
  real(dp), parameter :: mass  = 1.0_dp
  real(dp), parameter :: a     = 1.0_dp
  real(dp), parameter :: V0    = 1.0_dp
  real(dp), parameter :: gamma = sqrt(2.0_dp*mass*a*a*V0)/hbar   ! = sqrt(2)

  ! Delta strength (dimensionless)
  real(dp), parameter :: lambda = 1.0_dp  ! >0 gives one bound state

  ! ===== Grid / numerics =====
  integer,  parameter :: n_points   = 2001
  real(dp), parameter :: Xmax       = 8.0_dp
  real(dp), parameter :: pole_guard = 1.0e-14_dp

  real(dp) :: dx, xi(n_points), vgrid(n_points)
  integer  :: i

  ! Visible to derivs()
  real(dp) :: eps_current

  ! Energetics
  real(dp) :: eps_b, kappa

  ! ---- grid & (display) potential ----
  dx = (2.0_dp*Xmax)/real(n_points-1,dp)
  do i=1,n_points
    xi(i)    = -Xmax + dx*real(i-1,dp)
    vgrid(i) = 0.0_dp                ! delta cannot be represented on a grid; keep 0
  end do

  ! ---- bound-state parameters ----
  eps_b = - (gamma*gamma * lambda*lambda) / 4.0_dp
  kappa =   (gamma*gamma * lambda)       / 2.0_dp
  eps_current = eps_b

  print *, 'Delta well v(x) = - lambda delta(x): single bound state'
  print '(A,ES12.5)', 'gamma  = ', gamma
  print '(A,ES12.5)', 'lambda = ', lambda
  print '(A,ES12.5)', 'eps_b  = ', eps_b
  print '(A,ES12.5)', 'kappa  = ', kappa
  print '(A,F10.6)', 'Xmax   = ', Xmax
  print '(A,I6)',    'grid   = ', n_points

  call build_and_write_delta_bound(xi, vgrid, gamma, dx, eps_b, kappa)

  print *, 'Done: dw_n0.dat (single bound state)'

contains
  ! ===================== BUILD ψ AND WRITE FILE =====================
  subroutine build_and_write_delta_bound(xi, vgrid, gamma, dx, eps_b, kappa)
    real(dp), intent(in) :: xi(:), vgrid(:), gamma, dx, eps_b, kappa
    integer :: N, i0, j, u
    real(dp) :: yR(3), eR(3), yL(3), eL(3)
    real(dp), allocatable :: psi(:), psi_ex(:), errp(:), steperr(:)
    character(len=64) :: fname

    N = size(xi)
    allocate(psi(N), psi_ex(N), errp(N), steperr(N))
    psi = 0.0_dp; psi_ex = 0.0_dp; errp = 0.0_dp; steperr = 0.0_dp

    ! index nearest to xi=0
    i0 = nearest_index(xi, 0.0_dp)

    ! ---- Right side: start at 0+ with psi(0)=1, psi'(0+) = -kappa
    yR = [ 1.0_dp, -kappa, 0.0_dp ]   ! y(1)=psi, y(2)=psi', y(3)=xi
    psi(i0)      = yR(1)
    steperr(i0)  = 0.0_dp
    do j = i0+1, N
      call rk4n(dx, yR, 3, eR, gamma)   ! free-space ODE with eps_current=eps_b
      psi(j)     = yR(1)
      steperr(j) = eR(1)
    end do

    ! ---- Left side: start at 0- with psi(0)=1, psi'(0-) = +kappa
    yL = [ 1.0_dp, +kappa, 0.0_dp ]
    do j = i0-1, 1, -1
      call rk4n(-dx, yL, 3, eL, gamma)
      psi(j)     = yL(1)
      steperr(j) = eL(1)
    end do

    ! ---- normalize numerically (Simpson)
    call normalize_simpson(xi, psi)

    ! ---- exact psi for delta: sqrt(kappa) * exp(-kappa*|xi|)
    do j=1,N
      psi_ex(j) = sqrt(kappa) * exp( -kappa * abs(xi(j)) )
    end do
    ! (Optional) normalize psi_ex too (it already is)
    ! error
    errp = psi - psi_ex

    ! ---- write file (same columns)
    write(fname,'(A)') 'dw_n0.dat'
    open(newunit=u, file=trim(fname), status='replace')
    write(u,'(A)')        '# Delta Well (attractive) - single bound state n = 0'
    write(u,'(A,ES15.7)') '# gamma = ', gamma
    write(u,'(A,ES15.7)') '# lambda = ', lambda
    write(u,'(A,ES15.7)') '# epsilon (bound) = ', eps_b
    write(u,'(A,ES15.7)') '# kappa (decay)   = ', kappa
    write(u,'(A,F10.6)')  '# Xmax = ', Xmax
    write(u,'(A,I0)')     '# grid = ', N
    write(u,'(A)')        '# columns: xi   v(xi)   psi_num   psi_exact   psi_num^2   (psi_num-psi_exact)   step_err_psi'
    do j=1,N
      write(u,'(7ES20.10)') xi(j), vgrid(j), psi(j), psi_ex(j), psi(j)*psi(j), errp(j), steperr(j)
    end do
    close(u)

    deallocate(psi, psi_ex, errp, steperr)
  end subroutine build_and_write_delta_bound

  subroutine rk4n(dx, y, n, err, k)
    real(dp), intent(in)    :: dx
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: y(n)
    real(dp), intent(out)   :: err(n)
    real(dp), intent(in)    :: k
    real(dp) :: y1(n), y2(n), y3(n)
    real(dp) :: k1(n), k2(n), k3(n), k4(n)
    real(dp) :: rksum
    integer  :: iv

    call derivs(y, k1, n, k);  k1 = k1*dx;  y1 = y + 0.5_dp*k1
    call derivs(y1, k2, n, k); k2 = k2*dx;  y2 = y + 0.5_dp*k2
    call derivs(y2, k3, n, k); k3 = k3*dx;  y3 = y + k3
    call derivs(y3, k4, n, k); k4 = k4*dx

    do iv=1,n
      rksum   = (k1(iv) + 2.0_dp*k2(iv) + 2.0_dp*k3(iv) + k4(iv))/6.0_dp
      err(iv) = abs(rksum - k2(iv))
      y(iv)   = y(iv) + rksum
    end do
  end subroutine rk4n

  ! Free-space derivs: v(xi)=0 everywhere (delta handled via jump in ICs)
  ! y(1)=psi, y(2)=psi', y(3)=xi;  xi' = 1
  subroutine derivs(y, dydt, n, k)
    integer,  intent(in) :: n
    real(dp), intent(in) :: y(n)
    real(dp), intent(out):: dydt(n)
    real(dp), intent(in) :: k
    dydt(1) = y(2)
    dydt(2) = (k*k) * ( - eps_current ) * y(1)   ! v=0 => (v - eps) = -eps
    dydt(3) = 1.0_dp
  end subroutine derivs

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
    real(dp) :: I, Nrm
    allocate(f(size(x))); f = psi*psi
    I = simpson(x, f)
    if (I <= 0.0_dp) then
      deallocate(f); error stop 'normalize: nonpositive integral'
    end if
    Nrm = 1.0_dp/sqrt(I)
    psi = psi * Nrm
    deallocate(f)
  end subroutine normalize_simpson

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

end program tise_delta_well

