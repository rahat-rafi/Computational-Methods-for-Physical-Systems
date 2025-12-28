module relax_cfd
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)
contains

  ! psi_ij = (1/5)[ SC + (1/4) SS ] + (1/5) psi_ij + (1/40) SC
  subroutine relax_psi_sweep(psi, vort, nx, ny, wrel, ixst, ixfn, iyst, iyfn, dpmx)
    integer, intent(in) :: nx, ny
    integer, intent(inout) :: ixst, ixfn, iyst, iyfn
    real(dp), intent(inout) :: psi(nx,ny)
    real(dp), intent(in)    :: vort(nx,ny), wrel
    real(dp), intent(out)   :: dpmx
    integer :: i, j, incx, incy
    real(dp) :: oldv, updt, delta, SC, SS, SCO

    dpmx = 0.0_dp

    incx = merge( 1, -1, ixfn>=ixst )
    incy = merge( 1, -1, iyfn>=iyst )

    dpmx = 0.0_dp
    do j = iyst, iyfn, incy
      do i = ixst, ixfn, incx
        oldv = psi(i,j)
        SC  = psi(i+1,j) + psi(i-1,j) + psi(i,j+1) + psi(i,j-1)
        SS  = psi(i+1,j+1) + psi(i-1,j+1) + psi(i-1,j-1) + psi(i+1,j-1)
        SCO = vort(i+1,j) + vort(i-1,j) + vort(i,j+1) + vort(i,j-1)

        updt     = (SC + 0.25_dp*SS)/5.0_dp + vort(i,j)/5.0_dp + SCO/40.0_dp
        psi(i,j) = (1.0_dp - wrel)*oldv + wrel*updt

        delta = abs(psi(i,j) - oldv)
        if (delta > dpmx) dpmx = delta
      end do
    end do
  end subroutine relax_psi_sweep

  ! del^2omega = R * (psi_y omega_x − psi_x omega_y)
  ! omega_ij = (1/5)[ SC + (1/4) SS ] - (3R/10) * F
  ! F = 0.25 * [ (psi_{i,j+1}-psi_{i,j-1})(omega_{i+1,j}-omega_{i-1,j})
  !             - (psi_{i+1,j}-psi_{i-1,j})(omega_{i,j+1}-omega_{i,j-1}) ]
  subroutine relax_vort_sweep(psi, vort, nx, ny, wrel, R, ixst, ixfn, iyst, iyfn, domx)
    integer, intent(in) :: nx, ny
    integer, intent(inout) :: ixst, ixfn, iyst, iyfn
    real(dp), intent(in)    :: psi(nx,ny), wrel, R
    real(dp), intent(inout) :: vort(nx,ny)
    real(dp), intent(out)   :: domx
    integer :: i, j, incx, incy
    real(dp) :: oldv, updt, delta, SCO, SSO, Fterm

    domx = 0.0_dp

    incx = merge( 1, -1, ixfn>=ixst )
    incy = merge( 1, -1, iyfn>=iyst )

    domx = 0.0_dp
    do j = iyst, iyfn, incy
      do i = ixst, ixfn, incx
        oldv = vort(i,j)

        SCO = vort(i+1,j) + vort(i-1,j) + vort(i,j+1) + vort(i,j-1)
        SSO = vort(i+1,j+1) + vort(i-1,j+1) + vort(i-1,j-1) + vort(i+1,j-1)

        Fterm = 0.25_dp * ( (psi(i  ,j+1)-psi(i  ,j-1))*(vort(i+1,j)-vort(i-1,j))  &
                          - (psi(i+1,j)-psi(i-1,j))*(vort(i  ,j+1)-vort(i  ,j-1)) )

        updt      = (SCO + 0.25_dp*SSO)/5.0_dp - 3.0_dp*(R*Fterm)/10.0_dp
        vort(i,j) = (1.0_dp - wrel)*oldv + wrel*updt

        delta = abs(vort(i,j) - oldv)
        if (delta > domx) domx = delta
      end do
    end do
  end subroutine relax_vort_sweep

  ! omega_wall = -2 * (psi_fluid - psi_wall)
  subroutine set_wall_vorticity_from_psi(psi, vort, nx, ny, is_solid, psi_wall)
    integer, intent(in) :: nx, ny
    logical, intent(in) :: is_solid(nx,ny)
    real(dp), intent(in) :: psi(nx,ny), psi_wall
    real(dp), intent(inout) :: vort(nx,ny)
    integer :: i, j
    do j = 1, ny
      do i = 1, nx
        if (is_solid(i,j)) then
          if (i < nx .and. .not.is_solid(i+1,j)) then
            vort(i,j) = -2.0_dp * ( psi(i+1,j) - psi_wall )
          else if (i > 1 .and. .not.is_solid(i-1,j)) then
            vort(i,j) = -2.0_dp * ( psi(i-1,j) - psi_wall )
          else if (j < ny .and. .not.is_solid(i,j+1)) then
            vort(i,j) = -2.0_dp * ( psi(i,j+1) - psi_wall )
          else if (j > 1 .and. .not.is_solid(i,j-1)) then
            vort(i,j) = -2.0_dp * ( psi(i,j-1) - psi_wall )
          else
            vort(i,j) = 0.0_dp
          end if
        end if
      end do
    end do
  end subroutine set_wall_vorticity_from_psi
end module relax_cfd
