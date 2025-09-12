program spline_integrals
  implicit none
  
  integer, parameter :: dp = selected_real_kind(15, 300)
  integer, parameter :: n1 = 11            ! integer points: 0..10
  integer, parameter :: n2 = 21            ! half steps: 0 to 10 every 0.5
  real(dp), allocatable :: x1(:), y1(:), y2s(:)
  real(dp), allocatable :: x2(:), y2(:)
  real(dp) :: exact, T1, T2, s, comb
  integer  :: i
  real(dp) :: dum1, dum2
  integer  :: u

  allocate(x1(n1), y1(n1), y2s(n1))
  allocate(x2(n2), y2(n2))

  ! 1) Building two arrays for x, y
  do i = 1, n1
     x1(i) = real(i-1,dp)
     y1(i) = func(x1(i))
  end do

  ! 2) For spline interpolation
  call spline(x1, y1, n1, huge(1.0_dp), huge(1.0_dp), y2s)

  ! Building two additional arrays
  do i = 1, n2
     x2(i) = 0.5_dp * real(i-1,dp)
     call splintd(x1, y1, y2s, n1, x2(i), y2(i), dum1, dum2)
  end do

  ! 5) Finding the exact solution
  exact = F_integral(10.0_dp) - F_integral(0.0_dp)

  ! 5) Trapezoid estimates
  T1 = trapz(x1, y1)
  T2 = trapz(x2, y2)

  ! 4) Simpson (composite) on 11 points (x=0..10, h=1)
  s  = simp11(x1, y1)

  ! 6) Combined estimate (4/3)T2 - (1/3)T1
  comb = (4.0_dp/3.0_dp)*T2 - (1.0_dp/3.0_dp)*T1

  ! Output
  open(newunit=u, file='results_integrals.txt', status='replace', action='write')
  write(u,'(A)') "#  Estimate Value            Error"
  write(u,'(A, F20.8,1X,F20.8)') "Exact integral		", exact, 0.0_dp
  write(u,'(A, F20.8,1X,F20.8)') "Trapeziod (n=11)	", T1,    T1 - exact
  write(u,'(A, F20.8,1X,F20.8)') "Trapezoid (n=21)	", T2,    T2 - exact
  write(u,'(A, F20.8,1X,F20.8)') "Simpson (n=11)		", s,     s  - exact
  write(u,'(A, F20.8,1X,F20.8)') "Combined (4T2/3-T1/3)	", comb,  comb - exact
  close(u)

  print *, 'Wrote results_integrals.txt'
contains

  ! Function f(x)
  pure real(dp) function func(x) result(f)
    real(dp), intent(in) :: x
    f = -3.0e-3_dp*x**4 + 3.0e-3_dp*x**3 + 1.0e-2_dp*x**2 - 1.0e-1_dp*x + 0.2_dp
  end function func

  ! Integral of F(x)
  pure real(dp) function F_integral(x) result(val)
    real(dp), intent(in) :: x
    val = (-3.0e-3_dp/5.0_dp)*x**5 + (3.0e-3_dp/4.0_dp)*x**4 + &
          (1.0e-2_dp/3.0_dp)*x**3 - (1.0e-1_dp/2.0_dp)*x**2 + 0.2_dp*x
  end function F_integral

  ! Trapezoid for arbitrary (x,y) list (nonuniform OK)
  pure real(dp) function trapz(x, y) result(intval)
    real(dp), intent(in) :: x(:), y(:)
    integer :: i, n
    real(dp) :: acc
    n = size(x)
    if (size(y) /= n) error stop 'trapz: size mismatch'
    acc = 0.0_dp
    do i = 1, n-1
       acc = acc + 0.5_dp*(y(i)+y(i+1))*(x(i+1)-x(i))
    end do
    intval = acc
  end function trapz

  ! Composite Simpson specialized for x=0..10 (11 points, h=1)
  pure real(dp) function simp11(x, y) result(intval)
    real(dp), intent(in) :: x(:), y(:)
    integer :: n, i
    real(dp) :: sum4, sum2, h
    n = size(x)
    if (n /= 11 .or. size(y) /= 11) error stop 'simp11: requires 11 points'
    h = 1.0_dp
    sum4 = 0.0_dp
    sum2 = 0.0_dp
    do i = 2, 10, 2
       sum4 = sum4 + y(i)
    end do
    do i = 3, 9, 2
       sum2 = sum2 + y(i)
    end do
    intval = (h/3.0_dp) * ( y(1) + 4.0_dp*sum4 + 2.0_dp*sum2 + y(11) )
  end function simp11

  ! ----------------- Cubic spline precompute (NR-style) -----------------
  subroutine spline(x, y, n, yp1, ypn, y2)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: x(n), y(n), yp1, ypn
    real(dp), intent(out) :: y2(n)
    real(dp), allocatable :: u(:)
    integer :: i, k
    real(dp) :: sig, p, qn, un

    allocate(u(n))

    if (yp1 > 0.99e30_dp) then
       y2(1) = 0.0_dp
       u(1)  = 0.0_dp
    else
       y2(1) = -0.5_dp
       u(1)  = (3.0_dp/(x(2)-x(1))) * ((y(2)-y(1))/(x(2)-x(1)) - yp1)
    end if

    do i = 2, n-1
       sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
       p   = sig*y2(i-1) + 2.0_dp
       y2(i) = (sig - 1.0_dp) / p
       u(i)  = ( 6.0_dp * ( ( (y(i+1)-y(i)) /(x(i+1)-x(i)) ) - &
                             ( (y(i)  -y(i-1))/(x(i)  -x(i-1)) ) ) / (x(i+1)-x(i-1)) - sig*u(i-1) ) / p
    end do

    if (ypn > 0.99e30_dp) then
       qn = 0.0_dp
       un = 0.0_dp
    else
       qn = 0.5_dp
       un = (3.0_dp/(x(n)-x(n-1))) * (ypn - (y(n)-y(n-1))/(x(n)-x(n-1)))
    end if

    y2(n) = (un - qn*u(n-1)) / (qn*y2(n-1) + 1.0_dp)

    do k = n-1, 1, -1
       y2(k) = y2(k)*y2(k+1) + u(k)
    end do

    deallocate(u)
  end subroutine spline

  ! --------------- Spline interpolation with derivatives ----------------
  subroutine splintd(xa, ya, y2a, n, x, y, yprime, ydblprime)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: xa(n), ya(n), y2a(n), x
    real(dp), intent(out) :: y, yprime, ydblprime
    integer :: klo, khi, k
    real(dp) :: h, a, b

    klo = 1
    khi = n
    do
       if (khi - klo <= 1) exit
       k = (khi + klo)/2
       if (xa(k) > x) then
          khi = k
       else
          klo = k
       end if
    end do

    h = xa(khi) - xa(klo)
    if (h == 0.0_dp) error stop 'splintd: bad xa input (duplicate x)'

    a = (xa(khi) - x) / h
    b = (x - xa(klo)) / h

    y = a*ya(klo) + b*ya(khi) + ((a**3 - a)*y2a(klo) + (b**3 - b)*y2a(khi))*(h**2)/6.0_dp
    yprime = (ya(khi)-ya(klo))/h - (3.0_dp*a*a - 1.0_dp)*y2a(klo)*h/6.0_dp + (3.0_dp*b*b - 1.0_dp)*y2a(khi)*h/6.0_dp
    ydblprime = a*y2a(klo) + b*y2a(khi)
  end subroutine splintd

end program spline_integrals

