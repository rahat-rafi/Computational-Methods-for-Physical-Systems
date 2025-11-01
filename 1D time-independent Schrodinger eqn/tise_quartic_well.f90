program tise_quartic_well
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)

  ! ===== Dimensionless setup =====
  ! d^2 psi / d xi^2 = gamma^2 * ( v(xi) - eps ) * psi
  ! quartic: v(xi) = xi^4
  real(dp), parameter :: hbar  = 1.0_dp
  real(dp), parameter :: mass  = 1.0_dp
  real(dp), parameter :: a     = 1.0_dp
  real(dp), parameter :: V0    = 1.0_dp
  real(dp), parameter :: gamma = sqrt(2.0_dp*mass*a*a*V0)/hbar   ! = sqrt(2)

  ! ===== Grid / numerics =====
  integer,  parameter :: n_points     = 2001
  real(dp), parameter :: Xmax         = 8.0_dp
  real(dp), parameter :: delta_mid    = 1.0_dp   ! buffer past turning point
  real(dp), parameter :: edge_buf     = 1.0_dp   ! xi_R = Xmax - edge_buf
  real(dp), parameter :: pole_guard   = 1.0e-14_dp

  ! Energy search / tolerance
  integer,  parameter :: max_states   = 8        ! n = 0..7
  real(dp), parameter :: eps_floor    = 1.0e-8_dp
  real(dp), parameter :: eps_max      = 200.0_dp
  real(dp), parameter :: dE_init      = 0.25_dp
  real(dp), parameter :: tol_mismatch = 1.0e-10_dp
  real(dp), parameter :: tol_eps_span = 1.0e-12_dp

  ! Work arrays
  real(dp) :: dx, xi(n_points), vgrid(n_points)
  integer  :: i, s, found
  real(dp) :: eps_found(max_states), eps_even( (max_states+1)/2 ), eps_odd( max_states/2 )
  integer  :: par_found(max_states), n_quantum(max_states)

  ! Visible to derivs()
  real(dp) :: eps_current

  ! ---- grid & potential ----
  dx = (2.0_dp*Xmax)/real(n_points-1,dp)
  do i=1,n_points
    xi(i)    = -Xmax + dx*real(i-1,dp)
    vgrid(i) = v_of(xi(i))
  end do

  print *, 'Quartic Well v(x)=x^4: midpoint shooting + Wronskian mismatch'
  print '(A,F10.6)', 'gamma = ', gamma
  print '(A,F10.6)', 'Xmax  = ', Xmax
  print '(A,I6)',    'grid  = ', n_points

  ! ---- find even & odd, then merge to n=0..7 ----
  call scan_parity_levels(+1, size(eps_even), eps_even)   ! even states
  call scan_parity_levels(-1, size(eps_odd),  eps_odd)    ! odd states
  call merge_levels(eps_even, eps_odd, eps_found, par_found, n_quantum, found)

  if (found < max_states) then
    print *, 'WARNING: only found ', found, ' states'
  end if

  do s=1,found
    write(*,'(A,I0,A,ES14.7,A)') ' n=', n_quantum(s), '  eps=', eps_found(s), &
         merge(' (even)',' (odd) ', par_found(s)==+1)
  end do

  ! ---- build and write files ----
  do s=1,found
    call build_and_write_one( n_quantum(s), eps_found(s), par_found(s), xi, vgrid, gamma, dx )
  end do

  print *, 'Done: qw_n0.dat ... qw_n7.dat'

contains
  pure real(dp) function v_of(x) result(v)
    real(dp), intent(in) :: x
    v = x*x * x*x     ! x^4
  end function v_of

  real(dp) function mismatch_quartic(eps, parity, gamma, dx) result(F)
    real(dp), intent(in) :: eps, gamma, dx
    integer,  intent(in) :: parity
    real(dp) :: xi_R, xi_m, xi_cur, dx_hit
    real(dp) :: yL(3), eL(3), yR(3), eR(3)

    if (eps <= 0.0_dp) then
      F = -huge(1.0_dp);  return
    end if

    eps_current = eps
    xi_R = Xmax - edge_buf
    xi_m = eps**0.25_dp + delta_mid
    xi_m = max( 1.0_dp, xi_m )
    xi_m = min( xi_m, xi_R - 4.0_dp*dx )

    ! ---- Left: from 0 -> xi_m with parity ICs ----
    if (parity==+1) then
      yL = [1.0_dp, 0.0_dp, 0.0_dp]     ! even
    else
      yL = [0.0_dp, 1.0_dp, 0.0_dp]     ! odd
    end if
    xi_cur = 0.0_dp
    do while (xi_cur + dx < xi_m - 1.0e-14_dp)
      call rk4n(dx, yL, 3, eL, gamma)
      xi_cur = xi_cur + dx
    end do
    dx_hit = xi_m - xi_cur
    if (dx_hit > 0.0_dp) call rk4n(dx_hit, yL, 3, eL, gamma)
    if (abs(yL(1)) < pole_guard) then
      F = sign(1.0_dp, yL(2)) * huge(1.0_dp); return
    end if

    ! ---- Right: start at xi_R with WKB log-derivative, go to xi_m ----
    yR(3) = xi_R
    yR(1) = 1.0_dp
    yR(2) = -gamma * sqrt( max( v_of(xi_R) - eps, 0.0_dp ) ) * yR(1)
    xi_cur = xi_R
    do while (xi_cur - dx > xi_m + 1.0e-14_dp)
      call rk4n(-dx, yR, 3, eR, gamma)
      xi_cur = xi_cur - dx
    end do
    dx_hit = xi_m - xi_cur
    if (dx_hit /= 0.0_dp) call rk4n(dx_hit, yR, 3, eR, gamma)
    if (abs(yR(1)) < pole_guard) then
      F = -sign(1.0_dp, yR(2)) * huge(1.0_dp); return
    end if

    F = yL(1)*yR(2) - yL(2)*yR(1)    ! Wronskian
  end function mismatch_quartic

  subroutine scan_parity_levels(parity, need, eps_out)
    integer,  intent(in)  :: parity, need
    real(dp), intent(out) :: eps_out(need)
    real(dp) :: e_prev, e_cur, F_prev, F_cur, dE
    integer  :: foundp

    foundp = 0
    dE     = dE_init
    e_prev = eps_floor
    F_prev = mismatch_quartic(e_prev, parity, gamma, dx)

    do while (foundp < need .and. e_prev < eps_max)
      e_cur = min(e_prev + dE, eps_max)
      F_cur = mismatch_quartic(e_cur, parity, gamma, dx)

      if (is_finite(F_prev) .and. is_finite(F_cur)) then
        if (F_prev*F_cur <= 0.0_dp) then
          call bisect_eps(e_prev, e_cur, eps_out(foundp+1), parity, gamma, dx)
          if (eps_out(foundp+1) > 0.0_dp) then
            foundp = foundp + 1
            e_prev = eps_out(foundp) + 0.02_dp   ! move away from the root
            F_prev = mismatch_quartic(e_prev, parity, gamma, dx)
            cycle
          end if
        end if
      end if

      ! march on; gently increase step
      e_prev = e_cur
      F_prev = F_cur
      dE     = min( dE*1.15_dp, 5.0_dp )
    end do

    if (foundp < need) then
      eps_out(foundp+1:need) = -1.0_dp
    end if
  end subroutine scan_parity_levels

  subroutine merge_levels(e_even, e_odd, e_all, p_all, n_all, found)
    real(dp), intent(in)  :: e_even(:), e_odd(:)
    real(dp), intent(out) :: e_all(:)
    integer,  intent(out) :: p_all(:), n_all(:), found
    integer :: ne, no, i, j, k
    real(dp), allocatable :: tmpE(:)
    integer,  allocatable :: tmpP(:), tmpQN(:)   ! quantum-number labels

    ne = count(e_even > 0.0_dp)
    no = count(e_odd  > 0.0_dp)

    allocate(tmpE(ne+no), tmpP(ne+no), tmpQN(ne+no))
    k = 0
    do i=1,ne
      k=k+1; tmpE(k)=e_even(i); tmpP(k)=+1; tmpQN(k)=2*(i-1)    ! n = 0,2,4,6...
    end do
    do j=1,no
      k=k+1; tmpE(k)=e_odd(j);  tmpP(k)=-1; tmpQN(k)=2*(j-1)+1  ! n = 1,3,5,7...
    end do

    call sort_by_energy(k, tmpE, tmpP, tmpQN)

    found = min(k, size(e_all))
    e_all(1:found) = tmpE(1:found)
    p_all(1:found) = tmpP(1:found)
    n_all(1:found) = tmpQN(1:found)

    deallocate(tmpE, tmpP, tmpQN)
  end subroutine merge_levels

  subroutine sort_by_energy(ncnt, E, P, QN)
    integer, intent(in)     :: ncnt
    real(dp), intent(inout) :: E(:)
    integer,  intent(inout) :: P(:), QN(:)
    integer :: i, j, pt, qnt
    real(dp) :: et
    do i = 2, ncnt
      et  = E(i);  pt  = P(i);  qnt = QN(i)
      j = i - 1
      do while (j >= 1 .and. E(j) > et)
        E(j+1)  = E(j)
        P(j+1)  = P(j)
        QN(j+1) = QN(j)
        j = j - 1
      end do
      E(j+1)  = et
      P(j+1)  = pt
      QN(j+1) = qnt
    end do
  end subroutine sort_by_energy

  logical function is_finite(x)
    real(dp), intent(in) :: x
    is_finite = (x < huge(x) .and. x > -huge(x))
  end function is_finite

  subroutine bisect_eps(elo, ehi, eroot, parity, gamma, dx)
    real(dp), intent(in) :: elo, ehi, gamma, dx
    integer,  intent(in) :: parity
    real(dp), intent(out):: eroot
    real(dp) :: lo, hi, mid, flo, fhi, fmid

    lo = elo; hi = ehi
    flo = mismatch_quartic(lo, parity, gamma, dx)
    fhi = mismatch_quartic(hi, parity, gamma, dx)
    if (.not.(is_finite(flo) .and. is_finite(fhi))) then; eroot=-1.0_dp; return; end if
    if (flo*fhi > 0.0_dp) then; eroot=-1.0_dp; return; end if

    do
      mid  = 0.5_dp*(lo+hi)
      fmid = mismatch_quartic(mid, parity, gamma, dx)
      if (abs(fmid) < tol_mismatch) exit
      if (.not.is_finite(fmid)) then
        lo = mid; flo = mismatch_quartic(lo, parity, gamma, dx)
      else if (flo*fmid < 0.0_dp) then
        hi = mid; fhi = fmid
      else
        lo = mid; flo = fmid
      end if
      if (hi-lo < tol_eps_span) exit
    end do
    eroot = 0.5_dp*(lo+hi)
  end subroutine bisect_eps

  subroutine build_and_write_one(nsho, eps, parity, xi, vgrid, gamma, dx)
    integer,  intent(in) :: nsho, parity
    real(dp), intent(in) :: eps, gamma, dx
    real(dp), intent(in) :: xi(:), vgrid(:)

    integer :: N, i0, im, iR, j, u
    real(dp) :: y(3), e(3), yR(3), eR(3), scaleA
    real(dp), allocatable :: psi(:), psi_ex(:), errp(:), steperr(:), pr(:), pr_err(:)
    character(len=64) :: fname
    real(dp) :: xi_m

    N = size(xi)
    allocate(psi(N), psi_ex(N), errp(N), steperr(N), pr(N), pr_err(N))
    psi=0.0_dp; psi_ex=0.0_dp; errp=0.0_dp; steperr=0.0_dp; pr=0.0_dp; pr_err=0.0_dp

    eps_current = eps

    ! indices
    i0 = nearest_index(xi, 0.0_dp)
    xi_m = min( Xmax - edge_buf - 4.0_dp*dx, max(1.0_dp, eps**0.25_dp + delta_mid) )
    im = nearest_index(xi, xi_m)
    iR = N

    ! ---- center -> midpoint (parity ICs)
    if (parity==+1) then
      y = [1.0_dp, 0.0_dp, 0.0_dp]
    else
      y = [0.0_dp, 1.0_dp, 0.0_dp]
    end if
    psi(i0)      = y(1); steperr(i0)=0.0_dp
    do j=i0+1, im
      call rk4n(dx, y, 3, e, gamma)
      psi(j)      = y(1)
      steperr(j)  = e(1)
    end do

    ! ---- right boundary -> midpoint with WKB log-derivative
    yR(3) = xi(iR)
    yR(1) = 1.0_dp
    yR(2) = -gamma * sqrt( max( v_of(yR(3)) - eps, 0.0_dp ) ) * yR(1)
    pr(iR)     = yR(1); pr_err(iR) = 0.0_dp
    do j=iR-1, im, -1
      call rk4n(-dx, yR, 3, eR, gamma)
      pr(j)     = yR(1)
      pr_err(j) = eR(1)
    end do

    ! ---- scale & stitch at midpoint
    if (abs(pr(im)) < pole_guard) then
      scaleA = 1.0_dp
    else if (abs(psi(im)) <= 1.0e-30_dp .or. im<=2 .or. im>=N-1) then
      scaleA = psi(im) / max(pr(im), 1.0e-30_dp)
    else
      scaleA = ( (psi(im+1)-psi(im-1))/(2*dx) ) / max( (pr(im+1)-pr(im-1))/(2*dx), 1.0e-30_dp )
    end if

    do j=im, N
      psi(j)     = scaleA * pr(j)
      steperr(j) = abs(scaleA) * pr_err(j)
    end do

    ! ---- reflect to the left by parity: ψ(-x)=±ψ(x)
    do j=1, i0-1
      psi(j)     = merge(+1.0_dp, -1.0_dp, parity==+1) * psi(2*i0 - j)
      steperr(j) = steperr(2*i0 - j)
    end do

    ! ---- normalize (Simpson)
    call normalize_simpson(xi, psi)

    ! ---- "psi_exact" not available: keep zeros; errp = psi - 0 = psi
    errp   = psi
    psi_ex = 0.0_dp

    ! ---- write file (same column order as your HO files)
    write(fname,'(A,I0,A)') 'qw_n', nsho, '.dat'
    open(newunit=u, file=trim(fname), status='replace')
    write(u,'(A,I0)')     '# Quartic Well - state n = ', nsho
    write(u,'(A,I0)')     '# Parity (even=+1, odd=-1) = ', parity
    write(u,'(A,ES15.7)') '# gamma = ', gamma
    write(u,'(A,ES15.7)') '# epsilon (numeric) = ', eps
    write(u,'(A)')        '# Note: no analytic psi for quartic; psi_exact column is 0.0'
    write(u,'(A,F10.6)')  '# Xmax = ', Xmax
    write(u,'(A,I0)')     '# grid = ', N
    write(u,'(A)')        '# columns: xi   v(xi)   psi_num   psi_exact   psi_num^2   (psi_num-psi_exact)   step_err_psi'
    do j=1,N
      write(u,'(7ES20.10)') xi(j), vgrid(j), psi(j), psi_ex(j), psi(j)*psi(j), errp(j), steperr(j)
    end do
    close(u)

    deallocate(psi, psi_ex, errp, steperr, pr, pr_err)
  end subroutine build_and_write_one

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

  ! y(1)=psi, y(2)=psi', y(3)=xi;  xi' = 1
  subroutine derivs(y, dydt, n, k)
    integer,  intent(in) :: n
    real(dp), intent(in) :: y(n)
    real(dp), intent(out):: dydt(n)
    real(dp), intent(in) :: k
    dydt(1) = y(2)
    dydt(2) = (k*k) * ( v_of(y(3)) - eps_current ) * y(1)
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

end program tise_quartic_well

