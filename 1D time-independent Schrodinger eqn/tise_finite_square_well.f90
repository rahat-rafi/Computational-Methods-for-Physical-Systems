program schrodinger_finite_well_midpoint
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)

  ! ===== Dimensionless setup (PHY-711) =====
  ! xi = x/a, well: |xi|<1 (inside), v(xi)=0 inside, v(xi)=1 outside
  ! gamma = sqrt(2 m a^2 V0 / hbar^2), epsilon = E/V0 in (0,1) for bound states
  real(dp), parameter :: hbar = 1.0_dp, mass = 1.0_dp, a = 1.0_dp, V0 = 20.0_dp
  real(dp), parameter :: gamma = sqrt(2.0_dp*mass*a*a*V0)/hbar

  integer,  parameter :: n_points = 1001
  real(dp), parameter :: Xmax = 2.0_dp          ! domain in xi: [-Xmax, +Xmax]
  real(dp), parameter :: tol_mismatch = 1.0e-10_dp, tol_eps_span = 1.0e-12_dp
  real(dp), parameter :: eps_floor = 1.0e-12_dp, pole_guard = 1.0e-14_dp

  real(dp) :: dx, xi(n_points), vgrid(n_points)
  integer  :: i, found, s

  ! eigen storage
  integer,  parameter :: MAXF = 64
  real(dp) :: eps_found(MAXF)     ! epsilon = E/V0
  integer  :: par_found(MAXF)     ! +1 even, -1 odd

  ! ----- grid & v(xi) -----
  dx = (2.0_dp*Xmax)/real(n_points-1,dp)
  do i=1,n_points
    xi(i) = -Xmax + dx*real(i-1,dp)
    if (abs(xi(i)) < 1.0_dp) then
      vgrid(i) = 0.0_dp
    else
      vgrid(i) = 1.0_dp
    end if
  end do

  print *, 'Finite Square Well (dimensionless): midpoint shooting + exact ψ'
  print '(A,F8.4)', 'gamma = ', gamma
  print '(A,F8.4)', 'Xmax  = ', Xmax
  print '(A,I6)',   'grid  = ', n_points

  call scan_and_collect_roots( eps_floor, 1.0_dp-eps_floor, 8000, &
                               found, eps_found, par_found, gamma, dx )

  if (found==0) then
    print *, 'No bound states; increase V0 or a.'
    stop
  end if
  if (found>MAXF) found = MAXF

  do s=1,found
    call build_and_save_mode_with_exact( s, eps_found(s), par_found(s), xi, vgrid, gamma, dx )
  end do

  print *, 'Potential saved to: fsw_potential.dat'
  print *, 'Wavefunctions saved to: fsw_n*.dat'

contains
  ! ------------------------------------------------------------
  ! ------------------------------------------------------------
  ! Mismatch for parity at wall xi=+1, using exact step-to-wall.
  ! F(eps) = (psi'(1)/psi(1)) + kappa, with
  !     k    = gamma*sqrt(eps)         (inside)
  !     kappa = gamma*sqrt(1-eps)      (outside)
  ! Inside ODE (v=0): psi'' + k^2 psi = 0 ; ICs at xi=0:
  ! even: psi=1, psi'=0 ;  odd: psi=0, psi'=1
  ! ------------------------------------------------------------
  real(dp) function mismatch_parity(eps, parity, gamma, dx) result(F)
    implicit none
    real(dp), intent(in) :: eps, gamma, dx
    integer,  intent(in) :: parity
    real(dp) :: k, kap, y(2), err(2), xi_cur, dx_hit

    if (eps<=0.0_dp .or. eps>=1.0_dp) then
      F = huge(1.0_dp) * sign(1.0_dp, 0.5_dp - eps)
      return
    end if

    k   = gamma*sqrt(eps)
    kap = gamma*sqrt(1.0_dp-eps)

    if (parity==1) then
      y = [1.0_dp, 0.0_dp]
    else
      y = [0.0_dp, 1.0_dp]
    end if

    ! march with full steps until next step would pass the wall
    xi_cur = 0.0_dp
    do while (xi_cur + dx < 1.0_dp - 1.0e-14_dp)
      call rk4n(dx, y, 2, err, k)
      xi_cur = xi_cur + dx
    end do
    ! short step to land exactly at xi=1
    dx_hit = 1.0_dp - xi_cur
    if (dx_hit > 0.0_dp) call rk4n(dx_hit, y, 2, err, k)

    if (abs(y(1))<pole_guard) then
      F = sign(1.0_dp,y(2)) * huge(1.0_dp)
    else
      F = (y(2)/y(1)) + kap
    end if
  end function mismatch_parity

  subroutine scan_and_collect_roots(eps_lo, eps_hi, Nscan, found, eps_out, par_out, gamma, dx)
    implicit none
    real(dp), intent(in) :: eps_lo, eps_hi, gamma, dx
    integer,  intent(in) :: Nscan
    integer,  intent(out):: found
    real(dp), intent(out):: eps_out(:)
    integer,  intent(out):: par_out(:)

    integer :: p, i
    real(dp) :: e0, e1, de, m0, m1, er

    found = 0
    de = (eps_hi - eps_lo)/real(Nscan,dp)

    do p=1,2  ! 1->even(+1) 2->odd(-1)
      e0 = eps_lo
      m0 = mismatch_parity(e0, merge(1,-1,p==1), gamma, dx)
      do i=1,Nscan
        e1 = min(eps_lo + de*real(i,dp), eps_hi)
        m1 = mismatch_parity(e1, merge(1,-1,p==1), gamma, dx)

        if (is_finite(m0) .and. is_finite(m1)) then
          if (m0*m1 < 0.0_dp) then
            call bisect_eps(e0, e1, er, merge(1,-1,p==1), gamma, dx)
            if (er>0.0_dp) then
              found = found + 1
              eps_out(found) = er
              par_out(found) = merge(1,-1,p==1)
              if (found >= size(eps_out)) exit
            end if
            m0 = m1; e0 = e1
          else
            m0 = m1; e0 = e1
          end if
        else
          m0 = m1; e0 = e1
        end if
      end do
    end do

    call sort_levels(found, eps_out, par_out)
  end subroutine scan_and_collect_roots

  logical function is_finite(x)
    implicit none
    real(dp), intent(in) :: x
    is_finite = (x < huge(x) .and. x > -huge(x))
  end function is_finite

  subroutine bisect_eps(elo, ehi, eroot, parity, gamma, dx)
    implicit none
    real(dp), intent(in) :: elo, ehi, gamma, dx
    integer,  intent(in) :: parity
    real(dp), intent(out):: eroot
    real(dp) :: lo, hi, mid, flo, fhi, fmid

    lo = elo; hi = ehi
    flo = mismatch_parity(lo, parity, gamma, dx)
    fhi = mismatch_parity(hi, parity, gamma, dx)
    if (.not.(is_finite(flo) .and. is_finite(fhi))) then; eroot=-1.0_dp; return; end if
    if (flo*fhi > 0.0_dp) then; eroot=-1.0_dp; return; end if

    do
      mid = 0.5_dp*(lo+hi)
      fmid = mismatch_parity(mid, parity, gamma, dx)
      if (abs(fmid) < tol_mismatch) exit
      if (.not.is_finite(fmid)) then
        lo = mid; flo = mismatch_parity(lo, parity, gamma, dx)
      else if (flo*fmid < 0.0_dp) then
        hi = mid; fhi = fmid
      else
        lo = mid; flo = fmid
      end if
      if (hi-lo < tol_eps_span) exit
    end do
    eroot = 0.5_dp*(lo+hi)
  end subroutine bisect_eps

  subroutine sort_levels(n, earr, parr)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(inout) :: earr(:)
    integer,  intent(inout) :: parr(:)
    integer :: i,j,ptmp
    real(dp) :: etmp
    do i=2,n
      etmp = earr(i); ptmp = parr(i)
      j = i-1
      do while (j>=1 .and. earr(j) > etmp)
        earr(j+1)=earr(j); parr(j+1)=parr(j)
        j = j - 1
      end do
      earr(j+1)=etmp; parr(j+1)=ptmp
    end do
  end subroutine sort_levels

  subroutine build_and_save_mode_with_exact(n_state, eps, parity, xi, vgrid, gamma, dx)
    implicit none
    integer,  intent(in) :: n_state, parity
    real(dp), intent(in) :: eps, gamma, dx
    real(dp), intent(in) :: xi(:), vgrid(:)

    integer :: N, i, i0, j, j1, u
    real(dp) :: k, kap, y(2), step_e(2), xi_cur, dx_hit, psi_wall, sign_par
    real(dp) :: A_exact
    real(dp), allocatable :: psi(:), psi_exact(:), err_point(:), step_err(:)
    character(len=64) :: fname

    N = size(xi)
    allocate(psi(N), psi_exact(N), err_point(N), step_err(N))

    k   = gamma*sqrt(eps)
    kap = gamma*sqrt(1.0_dp-eps)
    sign_par = merge(1.0_dp, -1.0_dp, parity==1)

    ! center index (nearest xi=0)
    i0 = 1
    do i=2,N
      if (abs(xi(i)) < abs(xi(i0))) i0=i
    end do

    ! ICs at xi=0 by parity
    if (parity==1) then
      y = [1.0_dp, 0.0_dp]    ! even
    else
      y = [0.0_dp, 1.0_dp]    ! odd
    end if

    psi(i0)      = y(1)
    step_err(i0) = 0.0_dp

    ! ---- integrate to just before wall using full dx ----
    j = i0
    xi_cur = xi(j)
    do while (xi_cur + dx < 1.0_dp - 1.0e-14_dp .and. j+1 <= N)
      call rk4n(dx, y, 2, step_e, k)
      j = j + 1
      xi_cur = xi_cur + dx
      psi(j)      = y(1)
      step_err(j) = step_e(1)
    end do

    ! ---- short step to land exactly at xi=1 ----
    dx_hit = 1.0_dp - xi_cur
    if (dx_hit > 0.0_dp) call rk4n(dx_hit, y, 2, step_e, k)
    psi_wall = y(1)

    ! First grid node at/after the wall, then fill tail for xi > 1
    j1 = j + 1
    if (j1 <= N) then
      psi(j1)      = psi_wall * exp( -kap * ( xi(j1) - 1.0_dp ) )
      step_err(j1) = 0.0_dp
    end if
    do i = j1+1, N
      psi(i)      = psi_wall * exp( -kap * ( xi(i) - 1.0_dp ) )
      step_err(i) = 0.0_dp
    end do

    ! Mirror to xi<0 with parity
    do i = 1, i0-1
      psi(i)      = sign_par * psi(2*i0 - i)
      step_err(i) = step_err(2*i0 - i)
    end do

    ! Normalize numeric ψ
    call normalize_wavefunction(xi, psi)

    ! ---- exact, normalized ψ (a=1 in ξ units) ----
    call exact_amplitude(parity, k, kap, 1.0_dp, A_exact)
    do i=1,N
      if (abs(xi(i)) <= 1.0_dp) then
        if (parity==1) then
          psi_exact(i) = A_exact * cos(k*xi(i))
        else
          psi_exact(i) = A_exact * sin(k*xi(i))
        end if
      else
        if (parity==1) then
          psi_exact(i) = A_exact * cos(k) * exp( -kap * (abs(xi(i)) - 1.0_dp) )
        else
          psi_exact(i) = sign(1.0_dp, xi(i)) * A_exact * sin(k) * exp( -kap * (abs(xi(i)) - 1.0_dp) )
        end if
      end if
    end do
    call normalize_wavefunction(xi, psi_exact)

    ! error = psi_num - psi_exact
    err_point = psi - psi_exact

    write(fname,'(A,I0,A)') 'fsw_n', n_state, '.dat'
    open(newunit=u, file=trim(fname), status='replace')
    write(u,'(A,I0)') '# Finite Square Well - State n = ', n_state
    write(u,'(A,I0)') '# Parity (even=+1, odd=-1) = ', parity
    write(u,'(A,ES15.6)') '# epsilon = E/V0 = ', eps
    write(u,'(A,ES15.6)') '# k (inside)     = ', k
    write(u,'(A,ES15.6)') '# kappa (outside)= ', kap
    if (parity==1) then
      write(u,'(A,ES15.6)') '# residual k*tan(k)-kappa ≈ ', k*tan(k) - kap
    else
      write(u,'(A,ES15.6)') '# residual k*cot(k)+kappa ≈ ', k*(1.0_dp/tan(k)) + kap
    end if
    write(u,'(A)') '# columns: xi   v(xi)   psi_num   psi_exact   psi_num^2   (psi_num-psi_exact)   step_err_psi'
    do i=1,N
      write(u,'(7ES20.10)') xi(i), vgrid(i), psi(i), psi_exact(i), psi(i)*psi(i), err_point(i), step_err(i)
    end do
    close(u)

    deallocate(psi, psi_exact, err_point, step_err)
  end subroutine build_and_save_mode_with_exact

  ! Exact normalization amplitude (aw is half-width in ξ; here aw=1)
  ! Even: A = 1/sqrt( aw + sin(2k aw)/(2k) + cos^2(k aw)/kappa )
  ! Odd : A = 1/sqrt( aw - sin(2k aw)/(2k) + sin^2(k aw)/kappa )
  subroutine exact_amplitude(parity, k, kappa, aw, A)
    implicit none
    integer,  intent(in) :: parity
    real(dp), intent(in) :: k, kappa, aw
    real(dp), intent(out):: A
    real(dp) :: denom
    if (parity==1) then
      denom = aw + sin(2.0_dp*k*aw)/(2.0_dp*k) + (cos(k*aw)**2)/kappa
    else
      denom = aw - sin(2.0_dp*k*aw)/(2.0_dp*k) + (sin(k*aw)**2)/kappa
    end if
    A = 1.0_dp / sqrt(denom)
  end subroutine exact_amplitude

  subroutine rk4n(dx, y, n, err, k)
    implicit none
    real(dp), intent(in) :: dx
    integer, intent(in)  :: n
    real(dp), intent(inout) :: y(n)
    real(dp), intent(out)   :: err(n)
    real(dp), intent(in)    :: k        ! k = gamma*sqrt(epsilon)

    real(dp) :: y1(n), y2(n), y3(n)
    real(dp) :: k1(n), k2(n), k3(n), k4(n)
    real(dp) :: rksum
    integer :: ivar

    call derivs(y, k1, n, k)
    k1 = k1*dx
    y1 = y + 0.5_dp*k1

    call derivs(y1, k2, n, k)
    k2 = k2*dx
    y2 = y + 0.5_dp*k2

    call derivs(y2, k3, n, k)
    k3 = k3*dx
    y3 = y + k3

    call derivs(y3, k4, n, k)
    k4 = k4*dx

    do ivar=1,n
      rksum     = (k1(ivar) + 2.0_dp*k2(ivar) + 2.0_dp*k3(ivar) + k4(ivar))/6.0_dp
      err(ivar) = abs(rksum - k2(ivar))   ! diagnostic only
      y(ivar)   = y(ivar) + rksum
    end do
  end subroutine rk4n

  subroutine derivs(y, dydt, n, k)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: y(n)
    real(dp), intent(out):: dydt(n)
    real(dp), intent(in) :: k  ! IMPORTANT: k = gamma*sqrt(epsilon)
    ! Dimensionless inside-well ODE (v=0):  psi'' + k^2 psi = 0
    dydt(1) = y(2)          ! psi'
    dydt(2) = -(k*k)*y(1)   ! psi''
  end subroutine derivs

  real(dp) function integrate_simpson_any(x, y) result(intval)
    implicit none
    real(dp), intent(in) :: x(:), y(:)
    integer :: n, i, m
    real(dp) :: h, s_odd, s_even, trap

    n = size(x)
    if (size(y)/=n) error stop 'integrate_simpson_any: size mismatch'
    if (n<2)        error stop 'integrate_simpson_any: need >=2 points'

    h = x(2) - x(1)
    do i=3,n
      if (abs((x(i)-x(i-1)) - h) > 1.0e-12_dp) error stop 'integrate_simpson_any: nonuniform grid'
    end do
    m = n-1

    if (m==1) then
      intval = 0.5_dp*h*(y(1)+y(2)); return
    end if

    if (mod(m,2)==0) then
      s_odd=0; s_even=0
      do i=2,n-1,2; s_odd = s_odd + y(i);   end do
      do i=3,n-2,2; s_even= s_even + y(i); end do
      intval = (h/3.0_dp)*( y(1) + 4.0_dp*s_odd + 2.0_dp*s_even + y(n) )
    else
      s_odd=0; s_even=0
      do i=2,n-2,2; s_odd = s_odd + y(i);   end do
      do i=3,n-3,2; s_even= s_even + y(i); end do
      intval = (h/3.0_dp)*( y(1) + 4.0_dp*s_odd + 2.0_dp*s_even + y(n-1) )
      trap   = 0.5_dp*h*( y(n-1) + y(n) )
      intval = intval + trap
    end if
  end function integrate_simpson_any

  subroutine normalize_wavefunction(x, psi)
    implicit none
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: psi(:)
    real(dp), allocatable   :: f(:)
    real(dp) :: normint, norm
    integer :: n
    n = size(x)
    allocate(f(n)); f = psi*psi
    normint = integrate_simpson_any(x, f)
    if (normint<=0.0_dp) then
      deallocate(f); error stop 'normalize: nonpositive integral'
    end if
    norm = 1.0_dp/sqrt(normint)
    psi  = psi*norm
    deallocate(f)
  end subroutine normalize_wavefunction

end program schrodinger_finite_well_midpoint

