program schrodinger_parabolic_well_midpoint
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)

  ! ===== Dimensionless HO setup (LibreTexts/OpenStax) =====
  ! xi = x/a,  v(xi)=xi^2,  gamma = sqrt(2 m a^2 V0)/hbar.
  ! For HO: V0 = (1/2) m ω^2 a^2  =>  gamma = (m ω a^2)/ħ
  ! epsilon = E/V0, and exact  epsilon_n = (2n+1)/gamma.
  real(dp), parameter :: hbar  = 1.0_dp, mass = 1.0_dp, a = 1.0_dp
  real(dp), parameter :: omega = 1.0_dp
  real(dp), parameter :: V0    = 0.5_dp*mass*omega*omega*a*a
  real(dp), parameter :: gamma = sqrt(2.0_dp*mass*a*a*V0)/hbar   ! = m ω a^2 / ħ

  integer,  parameter :: n_points   = 2001
  real(dp), parameter :: Xmax       = 8.0_dp
  real(dp), parameter :: tol_mismatch = 1.0e-10_dp, tol_eps_span = 1.0e-12_dp
  real(dp), parameter :: eps_floor  = 1.0e-12_dp
  real(dp), parameter :: pole_guard = 1.0e-14_dp
  real(dp), parameter :: delta_mid  = 1.2_dp      ! xi_m = sqrt(eps) + delta_mid
  real(dp), parameter :: edge_buf   = 1.0_dp      ! xi_R = Xmax - edge_buf
  integer,  parameter :: max_states = 8

  real(dp) :: dx, xi(n_points), vgrid(n_points)
  integer  :: i, s, found
  real(dp) :: eps_found(2*max_states)
  integer  :: par_found(2*max_states)     ! +1 even, -1 odd
  integer  :: n_quantum(2*max_states)     ! SHO n

  ! Visible to derivs() via host association:
  real(dp) :: eps_current

  ! ----- grid & v(xi)=xi^2 -----
  dx = (2.0_dp*Xmax)/real(n_points-1,dp)
  do i=1,n_points
    xi(i)    = -Xmax + dx*real(i-1,dp)
    vgrid(i) = xi(i)*xi(i)
  end do

  print *, 'Parabolic (HO) Well: midpoint shooting (Wronskian) + exact Hermite ψ'
  print '(A,F10.6)', 'gamma = ', gamma
  print '(A,F10.6)', 'Xmax  = ', Xmax
  print '(A,I6)',    'grid  = ', n_points

  call find_parab_levels(max_states, eps_found, par_found, n_quantum, found, gamma, dx)

  if (found==0) then
    print *, 'No states found (check params)'; stop
  end if

  do s=1,found
    call build_and_save_mode_with_exact( s, eps_found(s), par_found(s), n_quantum(s), &
                                         xi, vgrid, gamma, dx )
  end do

  print *, 'Potential saved to: pw_potential.dat'
  print *, 'Wavefunctions saved to: pw_n*.dat'

contains
  ! ---------------- I/O ----------------
  ! -------- Mismatch: Wronskian at midpoint xi_m (> turning point) --------
  ! Left  (0 -> xi_m): parity ICs at xi=0
  ! Right (xi_R -> xi_m): decaying tail log-derivative  ψ'/ψ = -γ ξ  at xi_R
  ! F(ε) = ψ_L ψ'_R − ψ'_L ψ_R  (zero when left/right solutions are proportional)
  real(dp) function mismatch_parab(eps, parity, gamma, dx) result(F)
    implicit none
    real(dp), intent(in) :: eps, gamma, dx
    integer,  intent(in) :: parity
    real(dp) :: xi_m, xi_R, xi_cur, dx_hit
    real(dp) :: yL(3), eL(3), yR(3), eR(3)

    if (eps <= 0.0_dp) then
      F = -huge(1.0_dp);  return
    end if

    eps_current = eps
    xi_R = Xmax - edge_buf
    xi_m = max( sqrt(max(eps,0.0_dp)) + delta_mid, 1.0_dp )
    xi_m = min( xi_m, xi_R - 4.0_dp*dx )

    ! ---- Left (0 -> xi_m) ----
    if (parity==1) then
      yL = [1.0_dp, 0.0_dp, 0.0_dp]     ! even: ψ(0)=1, ψ'(0)=0
    else
      yL = [0.0_dp, 1.0_dp, 0.0_dp]     ! odd : ψ(0)=0, ψ'(0)=1
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

    ! ---- Right (xi_R -> xi_m) ----
    yR(3) = xi_R
    yR(1) = 1.0_dp
    yR(2) = -gamma*xi_R * yR(1)          ! asymptotic log-derivative
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

    F = yL(1)*yR(2) - yL(2)*yR(1)
  end function mismatch_parab

  ! -------- Root search (even/odd) around exact guesses ε_n=(2n+1)/γ --------
  subroutine find_parab_levels(maxst, eps_out, par_out, n_out, found, gamma, dx)
    implicit none
    integer,  intent(in)  :: maxst
    real(dp), intent(out) :: eps_out(:)
    integer,  intent(out) :: par_out(:), n_out(:), found
    real(dp), intent(in)  :: gamma, dx

    integer :: je, jo, n
    real(dp) :: e_guess, e_lo, e_hi, er
    real(dp), parameter :: pad = 1.0_dp   ! bracket half-width ≈ pad/γ

    found = 0

    ! even: n = 0,2,4,...
    je = 0
    do while (found < size(eps_out) .and. je < maxst)
      n       = 2*je
      e_guess = (2.0_dp*real(n,dp) + 1.0_dp)/gamma
      e_lo    = max(eps_floor, e_guess - pad/gamma)
      e_hi    = e_guess + pad/gamma
      call safe_expand_and_bisect(e_lo, e_hi, +1, er, gamma, dx)
      if (er > 0.0_dp) then
        found = found + 1
        eps_out(found) = er;  par_out(found) = +1;  n_out(found) = n
        write(*,'(A,I0,A,ES14.7,A,ES14.7)') 'even n=', n, '  eps_num=', er, '  eps_exact=', (2.0_dp*real(n,dp)+1.0_dp)/gamma
      end if
      je = je + 1
    end do

    ! odd: n = 1,3,5,...
    jo = 0
    do while (found < size(eps_out) .and. jo < maxst)
      n       = 2*jo + 1
      e_guess = (2.0_dp*real(n,dp) + 1.0_dp)/gamma
      e_lo    = max(eps_floor, e_guess - pad/gamma)
      e_hi    = e_guess + pad/gamma
      call safe_expand_and_bisect(e_lo, e_hi, -1, er, gamma, dx)
      if (er > 0.0_dp) then
        found = found + 1
        eps_out(found) = er;  par_out(found) = -1;  n_out(found) = n
        write(*,'(A,I0,A,ES14.7,A,ES14.7)') 'odd  n=', n, '  eps_num=', er, '  eps_exact=', (2.0_dp*real(n,dp)+1.0_dp)/gamma
      end if
      jo = jo + 1
    end do

    call sort_levels(found, eps_out, par_out, n_out)
  end subroutine find_parab_levels

  subroutine safe_expand_and_bisect(e_lo, e_hi, parity, eroot, gamma, dx)
    implicit none
    real(dp), intent(inout) :: e_lo, e_hi
    integer,  intent(in)    :: parity
    real(dp), intent(out)   :: eroot
    real(dp), intent(in)    :: gamma, dx
    real(dp) :: F_lo, F_hi
    integer :: tries

    F_lo = mismatch_parab(e_lo, parity, gamma, dx)
    F_hi = mismatch_parab(e_hi, parity, gamma, dx)

    tries = 0
    do while ( ( .not. is_finite(F_lo) .or. .not. is_finite(F_hi) .or. F_lo*F_hi > 0.0_dp ) .and. tries < 80 )
      e_lo = max(eps_floor, e_lo - 0.4_dp/gamma)
      e_hi = e_hi + 0.4_dp/gamma
      F_lo = mismatch_parab(e_lo, parity, gamma, dx)
      F_hi = mismatch_parab(e_hi, parity, gamma, dx)
      tries = tries + 1
    end do

    if ( .not.(is_finite(F_lo) .and. is_finite(F_hi)) ) then
      eroot = -1.0_dp; return
    end if
    if (F_lo*F_hi > 0.0_dp) then
      eroot = -1.0_dp; return
    end if

    call bisect_eps(e_lo, e_hi, eroot, parity, gamma, dx)
  end subroutine safe_expand_and_bisect

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
    flo = mismatch_parab(lo, parity, gamma, dx)
    fhi = mismatch_parab(hi, parity, gamma, dx)
    if (.not.(is_finite(flo) .and. is_finite(fhi))) then; eroot=-1.0_dp; return; end if
    if (flo*fhi > 0.0_dp) then; eroot=-1.0_dp; return; end if

    do
      mid  = 0.5_dp*(lo+hi)
      fmid = mismatch_parab(mid, parity, gamma, dx)
      if (abs(fmid) < tol_mismatch) exit
      if (.not.is_finite(fmid)) then
        lo = mid; flo = mismatch_parab(lo, parity, gamma, dx)
      else if (flo*fmid < 0.0_dp) then
        hi = mid; fhi = fmid
      else
        lo = mid; flo = fmid
      end if
      if (hi-lo < tol_eps_span) exit
    end do
    eroot = 0.5_dp*(lo+hi)
  end subroutine bisect_eps

  subroutine sort_levels(n, earr, parr, narr)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(inout) :: earr(:)
    integer,  intent(inout) :: parr(:), narr(:)
    integer :: i,j,ptmp,ntmp
    real(dp) :: etmp
    do i=2,n
      etmp = earr(i); ptmp = parr(i); ntmp = narr(i)
      j = i-1
      do while (j>=1 .and. earr(j) > etmp)
        earr(j+1)=earr(j); parr(j+1)=parr(j); narr(j+1)=narr(j)
        j = j - 1
      end do
      earr(j+1)=etmp; parr(j+1)=ptmp; narr(j+1)=ntmp
    end do
  end subroutine sort_levels

  ! -------------- Build ψ_num (integrate both sides) + exact ψ_n, write file --------------
  subroutine build_and_save_mode_with_exact(idx, eps, parity, nsho, xi, vgrid, gamma, dx)
    implicit none
    integer,  intent(in) :: idx, parity, nsho
    real(dp), intent(in) :: eps, gamma, dx
    real(dp), intent(in) :: xi(:), vgrid(:)

    integer :: N, i, i0, j, u
    real(dp) :: y(3), step_e(3)
    real(dp) :: z, Hn, pi, gqrt, normfac_exact
    real(dp), allocatable :: psi(:), psi_exact(:), err_point(:), step_err(:)
    character(len=64) :: fname
    real(dp) :: eps_exact, abs_err_E, rel_err_E

    N = size(xi)
    allocate(psi(N), psi_exact(N), err_point(N), step_err(N))

    eps_current = eps

    ! center index closest to xi=0
    i0 = 1
    do i=2,N
      if (abs(xi(i)) < abs(xi(i0))) i0=i
    end do

    ! --- integrate to the right from xi=0 ---
    if (parity==1) then
      y = [1.0_dp, 0.0_dp, 0.0_dp]   ! even
    else
      y = [0.0_dp, 1.0_dp, 0.0_dp]   ! odd
    end if
    psi(i0)      = y(1)
    step_err(i0) = 0.0_dp
    do j = i0+1, N
      call rk4n(dx, y, 3, step_e, gamma)
      psi(j)      = y(1)
      step_err(j) = step_e(1)
    end do

    ! --- integrate to the left from xi=0 (NO mirroring) ---
    if (parity==1) then
      y = [1.0_dp, 0.0_dp, 0.0_dp]
    else
      y = [0.0_dp, 1.0_dp, 0.0_dp]
    end if
    do j = i0-1, 1, -1
      call rk4n(-dx, y, 3, step_e, gamma)
      psi(j)      = y(1)
      step_err(j) = step_e(1)
    end do

    ! normalize numerical ψ
    call normalize_wavefunction(xi, psi)

    ! ---- Exact HO ψ_n(ξ) from LibreTexts:
    ! ψ_n(x)=N_n e^{-β^2 x^2/2} H_n(βx) with β=√(mω/ħ).
    ! In ξ: βx = √γ ξ and  N_n = (γ/π)^{1/4} / sqrt(2^n n!)
    pi   = 4.0_dp*atan(1.0_dp)
    gqrt = sqrt(gamma)
    normfac_exact = (gamma**0.25_dp) / (pi**0.25_dp * sqrt( two_pow(nsho) * real(fact_int(nsho),dp) ))
    do i=1,N
      z = gqrt * xi(i)
      call hermiteH_eval(nsho, z, Hn)
      psi_exact(i) = normfac_exact * exp(-0.5_dp*z*z) * Hn
    end do
    ! (Optional) You may skip this; exact N_n already normalizes over (-∞,∞).
    ! call normalize_wavefunction(xi, psi_exact)

    ! energy comparison
    eps_exact = (2.0_dp*real(nsho,dp) + 1.0_dp)/gamma
    abs_err_E = abs( eps - eps_exact )
    rel_err_E = abs_err_E / max(eps_exact, 1.0e-30_dp)

    err_point = psi - psi_exact

    write(fname,'(A,I0,A)') 'pw_n', idx, '.dat'
    open(newunit=u, file=trim(fname), status='replace')
    write(u,'(A,I0)')         '# Parabolic (HO) Well - State index = ', idx
    write(u,'(A,I0)')         '# Parity (even=+1, odd=-1) = ', parity
    write(u,'(A,I0)')         '# SHO quantum number n = ', nsho
    write(u,'(A,ES15.6)')     '# epsilon (numeric)  = ', eps
    write(u,'(A,ES15.6)')     '# epsilon (analytic) = ', eps_exact
    write(u,'(A,ES15.6)')     '# |eps - eps_exact|  = ', abs_err_E
    write(u,'(A,ES15.6)')     '# rel. error (eps)   = ', rel_err_E
    write(u,'(A,F10.6)')      '# Xmax = ', Xmax
    write(u,'(A)')            '# columns: xi   v(xi)   psi_num   psi_exact   psi_num^2   (psi_num-psi_exact)   step_err_psi'
    do i=1,N
      write(u,'(7ES20.10)') xi(i), vgrid(i), psi(i), psi_exact(i), psi(i)*psi(i), err_point(i), step_err(i)
    end do
    close(u)

    deallocate(psi, psi_exact, err_point, step_err)
  end subroutine build_and_save_mode_with_exact

  ! ---- physicists' Hermite H_n(z) via 3-term recurrence ----
  subroutine hermiteH_eval(n, z, Hn)
    implicit none
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
  end subroutine hermiteH_eval

  pure real(dp) function two_pow(n) result(v)
    implicit none
    integer, intent(in) :: n
    integer :: i
    v = 1.0_dp
    do i=1,n
      v = 2.0_dp*v
    end do
  end function two_pow

  pure integer function fact_int(n) result(fv)
    implicit none
    integer, intent(in) :: n
    integer :: k
    fv = 1
    if (n>1) then
      do k=2,n
        fv = fv*k
      end do
    end if
  end function fact_int

  ! ================== REQUIRED NAMES (unchanged) ==================
  subroutine rk4n(dx, y, n, err, k)
    implicit none
    real(dp), intent(in) :: dx
    integer, intent(in)  :: n
    real(dp), intent(inout) :: y(n)
    real(dp), intent(out)   :: err(n)
    real(dp), intent(in)    :: k       ! here: k = gamma

    real(dp) :: y1(n), y2(n), y3(n)
    real(dp) :: k1(n), k2(n), k3(n), k4(n)
    real(dp) :: rksum
    integer :: ivar

    call derivs(y, k1, n, k);  k1 = k1*dx;  y1 = y + 0.5_dp*k1
    call derivs(y1, k2, n, k); k2 = k2*dx;  y2 = y + 0.5_dp*k2
    call derivs(y2, k3, n, k); k3 = k3*dx;  y3 = y + k3
    call derivs(y3, k4, n, k); k4 = k4*dx

    do ivar=1,n
      rksum     = (k1(ivar) + 2.0_dp*k2(ivar) + 2.0_dp*k3(ivar) + k4(ivar))/6.0_dp
      err(ivar) = abs(rksum - k2(ivar))
      y(ivar)   = y(ivar) + rksum
    end do
  end subroutine rk4n

  ! d^2ψ/dξ^2 = γ^2 (ξ^2 − ε) ψ ; carry ξ as y(3) with ξ' = 1
  subroutine derivs(y, dydt, n, k)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: y(n)
    real(dp), intent(out):: dydt(n)
    real(dp), intent(in) :: k     ! k = gamma
    dydt(1) = y(2)
    dydt(2) = (k*k) * ( y(3)*y(3) - eps_current ) * y(1)
    dydt(3) = 1.0_dp
  end subroutine derivs

  ! ----- Composite Simpson over uniform grid -----
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

end program schrodinger_parabolic_well_midpoint

